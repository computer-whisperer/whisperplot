#include <atomic>
#include <string>
#include <thread>
#include <vector>

#include "thread_mgr.hpp"
#include "entry_sizes.hpp"
#include "buffer.hpp"
#include "encoding.hpp"
#include "plotter.hpp"
#include "chiapos_util.hpp"
#include "status_update.hpp"

using namespace std;

template <PlotConf conf>
template <int8_t table_index>
void Plotter<conf>::phase3ThreadA(
        uint32_t cpu_id,
        atomic<uint64_t>* coordinator,
        vector<Park*>* temporary_parks,
        entries_used_type* entries_used,
        vector<p2_final_positions_type>* final_positions,
        Buffer* output_buffer,
        vector<vector<Park*>>* final_parks,
        uint64_t start_offset)
{
    PinToCpuid(cpu_id);
    vector<uint128_t> temporary_park_line_points(LinePointUIDPackedEntry<table_index>::max_entries_per_row);
    vector<uint128_t> park_line_points(kEntriesPerPark);
    vector<uint8_t> delta_buff(kEntriesPerPark);
    vector<uint64_t> stub_buff(kEntriesPerPark);
	while (true)
	{
		uint64_t park_size_bytes = EntrySizes::CalculateParkSize(conf.K, table_index+1);

		uint32_t park_id = coordinator->fetch_add(1);
		// Output everything for this park id

		uint64_t park_starting_pos = kEntriesPerPark*park_id;
		uint64_t temporary_park_id = 0;

		uint64_t park_entry_i = 0;

		while (temporary_park_id < temporary_parks->size())
		{
		    uint64_t old_end_pos = temporary_parks->at(temporary_park_id)->start_pos + temporary_parks->at(temporary_park_id)->size();
		    uint64_t new_end_pos = old_end_pos;
            if (table_index < 5)
                new_end_pos = (*final_positions)[table_index].get(old_end_pos).getY();
			if (new_end_pos > park_starting_pos)
			{

                temporary_parks->at(temporary_park_id)->readEntries(temporary_park_line_points);
                for (uint64_t i = 0; i < temporary_parks->at(temporary_park_id)->size(); i++)
				{
				    uint64_t old_pos = temporary_parks->at(temporary_park_id)->start_pos+i;
				    uint64_t new_pos = old_pos;
				    if (table_index < 5)
				        new_pos = (*final_positions)[table_index].get(old_pos).getY();
				    if (new_pos < park_starting_pos)
                    {
				        continue;
                    }
				    if ((table_index == 5) || entries_used->get(old_pos).getY())
                    {

				        uint128_t line_point = temporary_park_line_points[i];

				        if (table_index > 0)
                        {
                            // Re-encode line point to compensate for dropped entries in previous table
                            auto res = Encoding::LinePointToSquare(line_point);
                            uint64_t x_new = (*final_positions)[table_index-1].get(res.first).getY();
                            uint64_t y_new = (*final_positions)[table_index-1].get(res.second).getY();
                            line_point = Encoding::SquareToLinePoint(x_new, y_new);
                        }
                        park_line_points[park_entry_i++] = line_point;
                    }
					if (park_entry_i == kEntriesPerPark)
					{
						break;
					}
				}
			}
			if (park_entry_i == kEntriesPerPark)
			{
				break;
			}
            temporary_park_id++;
		}

		if (park_entry_i == 0)
		{
			break;
		}

		// Kick the insertion pointer down the output buffer
		output_buffer->GetInsertionOffset(park_size_bytes);

		// Pad the end with more stuff
		park_line_points.resize(park_entry_i);

		Park* p;
		if (table_index == 0) {
            using cp = CompressedPark<conf.K * 2, conf.K, kStubMinusBits, (uint32_t)(kMaxAverageDeltaTable1*100), (uint32_t)(kRValues[table_index]*100)>;
            cp* p2 = new cp(kEntriesPerPark);
            assert(p2->GetSpaceNeeded() == park_size_bytes);
            p = p2;
            p2->useBuffers(&stub_buff, &delta_buff);
        }
		else {
		    using cp = CompressedPark<conf.K * 2, conf.K, kStubMinusBits, (uint32_t)(kMaxAverageDelta*100), (uint32_t)(kRValues[table_index]*100)>;
            cp* p2 = new cp(kEntriesPerPark);
            assert(p2->GetSpaceNeeded() == park_size_bytes);
            p = p2;
            p2->useBuffers(&stub_buff, &delta_buff);
        }
        p->bind(output_buffer->data + start_offset + park_id * park_size_bytes);
		p->addEntries(park_line_points);



        delete p;
	}
}

template <PlotConf conf>
template <int8_t table_index>
void Plotter<conf>::phase3DoTable()
{
    StatusUpdate::StartSeg("2." + to_string((uint32_t)table_index) + ".0S");

    // Wait for the final positions table to be ready
    final_parks.emplace_back();
    entries_used_type* current_entries_used = nullptr;
    if (table_index < 5)
    {
        phase2b_threads[table_index].join();
        current_entries_used = &(entries_used->at(table_index));
    }
    final_table_begin_pointers[table_index] = bswap_64(*(output_buffer->insert_pos));


    StatusUpdate::StartSeg("2." + to_string((uint32_t)table_index) + ".1M");

    vector<thread> threads;
    std::atomic<uint64_t> coordinator = 0;
    uint64_t start_offset = *output_buffer->insert_pos;
    for (auto & cpu_id : cpu_ids)
    {
        threads.push_back(thread(
                Plotter<conf>::phase3ThreadA<table_index>,
                cpu_id,
                &coordinator,
                &(phase1_graph_parks[table_index]),
                current_entries_used,
                phase2_final_positions,
                output_buffer,
                &final_parks,
                start_offset));
    }

    for (auto &it: threads)
    {
        it.join();
    }

    StatusUpdate::StartSeg("2." + to_string((uint32_t)table_index) + ".2S");

    delete buffers[table_index];
}

template <PlotConf conf>
void Plotter<conf>::phase3()
{
    StatusUpdate::StartSeg("2.-1.1S");


    uint64_t phase_start_seconds = time(nullptr);
	// Try to predict an upper bound to the final file size
	uint64_t predicted_file_size_bytes = 128;// header
	for (uint8_t table_index = 0; table_index < 6; table_index++)
	{
		uint64_t park_size_bytes = EntrySizes::CalculateParkSize(conf.K, table_index+1);
		auto final_park = *(phase1_graph_parks[table_index].end() - 1);
		uint64_t num_parks = (final_park->start_pos + final_park->size() + kEntriesPerPark-1)/kEntriesPerPark;
		predicted_file_size_bytes += park_size_bytes*num_parks;
	}

    uint64_t total_C1_entries = cdiv(((*(phase1_final_parks.end() - 1))->start_pos + (*(phase1_final_parks.end() - 1))->size()), kCheckpoint1Interval);
    uint64_t total_C2_entries = cdiv(total_C1_entries, kCheckpoint2Interval);
    uint32_t size_C3 = EntrySizes::CalculateC3Size(conf.K);
    predicted_file_size_bytes += (total_C1_entries + 1) * (Util::ByteAlign(conf.K) / 8) + (total_C2_entries + 1) * (Util::ByteAlign(conf.K) / 8) + (total_C1_entries)*size_C3;

    output_buffer = new Buffer(predicted_file_size_bytes<<1, filename);
    //output_buffer = new Buffer(predicted_file_size_bytes*1.5);

    // 19 bytes  - "Proof of Space Plot" (utf-8)
    // 32 bytes  - unique plot id
    // 1 byte    - k
    // 2 bytes   - format description length
    // x bytes   - format description
    // 2 bytes   - memo length
    // x bytes   - memo

    output_buffer->InsertString("Proof of Space Plot");
    output_buffer->InsertData((void*)id, kIdLen);
    *(uint8_t*)(output_buffer->data + output_buffer->GetInsertionOffset(1)) = conf.K;

    *(uint16_t*)(output_buffer->data + output_buffer->GetInsertionOffset(2)) = bswap_16(kFormatDescription.size());
    output_buffer->InsertString(kFormatDescription);
    *(uint16_t*)(output_buffer->data + output_buffer->GetInsertionOffset(2)) = bswap_16(memo_size);
    output_buffer->InsertData((void*)memo, memo_size);

    pointer_table_offset = output_buffer->InsertData(final_table_begin_pointers, sizeof(final_table_begin_pointers));


    phase3DoTable<0>();
    phase3DoTable<1>();
    phase3DoTable<2>();
    phase3DoTable<3>();
    phase3DoTable<4>();
    phase3DoTable<5>();
	final_table_begin_pointers[6] = bswap_64(*(output_buffer->insert_pos));

    StatusUpdate::EndSeg();
    cout << "Phase 3 finished in " << time(nullptr) - phase_start_seconds << "s" << endl;
}

#include "explicit_templates.hpp"