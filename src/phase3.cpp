#include <atomic>
#include <string>
#include <thread>
#include <vector>

#include "thread_mgr.hpp"
#include "entry_sizes.hpp"
#include "buffer.hpp"
#include "encoding.hpp"
#include "plotter.hpp"
#include "util.hpp"

using namespace std;

template <uint8_t K>
template <int8_t table_index>
void Plotter<K>::phase3ThreadA(
        uint32_t cpu_id,
        atomic<uint64_t>* coordinator,
        vector<DeltaPark<line_point_delta_len_bits>*>* temporary_parks,
        AtomicPackedArray<BooleanPackedEntry, max_entries_per_graph_table>* entries_used,
        AtomicPackedArray<SimplePackedEntry<uint64_t,K>, max_entries_per_graph_table>* final_positions,
        Buffer* output_buffer)
{
    PinToCpuid(cpu_id);
    vector<uint128_t> temporary_park_line_points(LinePointEntryUIDPackedEntry<K>::max_entries_per_sort_row);
    vector<uint128_t> park_line_points(kEntriesPerPark);
    vector<uint128_t> park_line_points_check(kEntriesPerPark);
	while (true)
	{
		uint64_t park_size_bytes = EntrySizes::CalculateParkSize(K, table_index+1);

		uint32_t park_id = (*coordinator)++;
		// Output everything for this park id

		uint64_t park_starting_pos = kEntriesPerPark*park_id;
		uint64_t temporary_park_id = 0;

		uint64_t park_entry_i = 0;

		while (temporary_park_id < temporary_parks->size())
		{
			if (temporary_parks->at(temporary_park_id)->start_pos + temporary_parks->at(temporary_park_id)->size() > park_starting_pos)
			{

                temporary_parks->at(temporary_park_id)->readEntries(temporary_park_line_points);
                uint128_t prev_new_line_point = 0;
                uint128_t prev_old_line_point = 0;
                for (uint64_t i = 0; i < temporary_parks->at(temporary_park_id)->size(); i++)
				{
				    uint64_t pos = temporary_parks->at(temporary_park_id)->start_pos+i;
				    if (pos < park_starting_pos)
                    {
				        continue;
                    }
				    if (entries_used->read(pos).value)
                    {

				        uint128_t line_point = temporary_park_line_points[i];

				        if (table_index > 0)
                        {
                            // Re-encode line point to compensate for dropped entries in previous table
                            auto res = Encoding::LinePointToSquare(line_point);
                            uint64_t x_new = final_positions->read(res.first).value;
                            uint64_t y_new = final_positions->read(res.second).value;
                            assert(x_new <= res.first);
                            assert(y_new <= res.second);
                            line_point = Encoding::SquareToLinePoint(x_new, y_new);
                            assert(temporary_park_line_points[i]>=prev_old_line_point);
                            assert(line_point>=prev_new_line_point);
                        }
                        park_line_points[park_entry_i++] = line_point;
                        prev_old_line_point = temporary_park_line_points[i];
                        prev_new_line_point = line_point;
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

		// Pad the end with more stuff
		while (park_entry_i < kEntriesPerPark)
        {
            park_line_points[park_entry_i] = park_line_points[park_entry_i-1]+1;
            park_entry_i++;
        }

		Park* p;
		if (table_index == 0) {
            p = new CompressedPark<K * 2, K, kStubMinusBits, kMaxAverageDeltaTable1, kRValues[table_index]>(kEntriesPerPark);
        }
		else {
            p = new CompressedPark<K * 2, K, kStubMinusBits, kMaxAverageDelta, kRValues[table_index]>(kEntriesPerPark);
        }
        p->bind(output_buffer->data + *(output_buffer->insert_pos) + park_id * park_size_bytes);
		p->addEntries(park_line_points);
		p->readEntries(park_line_points_check);

		for (uint64_t i = 0; i < kEntriesPerPark; i++)
        {
		    assert(park_line_points[i] == park_line_points_check[i]);
        }

		delete p;
	}
}

template <uint8_t K>
template <int8_t table_index>
void Plotter<K>::phase3DoTable()
{
    // Wait for the final positions table to be ready
    phase2b_threads[table_index].join();
    final_table_begin_pointers[table_index] = bswap_64(*(output_buffer->insert_pos));
    cout << "Part A"<< (uint32_t)table_index;
    uint64_t start_seconds = time(nullptr);
    vector<thread> threads;
    std::atomic<uint64_t> coordinator = 0;
    for (auto & cpu_id : cpu_ids)
    {
        threads.push_back(thread(
                Plotter<K>::phase3ThreadA<table_index>,
                cpu_id,
                &coordinator,
                &(graph_parks[table_index]),
                &(entries_used[table_index]),
                &(final_positions[table_index-1]),
                output_buffer));
    }

    for (auto &it: threads)
    {
        it.join();
    }
    cout << " (" << time(nullptr) - start_seconds << "s)" << endl;


    // setup output_table_offset for next iteration
    uint64_t park_size_bytes = EntrySizes::CalculateParkSize(K, table_index+1);
    uint64_t num_entries = (*(graph_parks[table_index].end()-1))->start_pos + (*(graph_parks[table_index].end()-1))->size();
    uint64_t parks_needed = (num_entries+kEntriesPerPark-1)/kEntriesPerPark;
    output_buffer->GetInsertionOffset(parks_needed*park_size_bytes);

    delete buffers[table_index];
}

template<uint8_t K>
void Plotter<K>::phase3()
{
    uint64_t phase_start_seconds = time(nullptr);
	// Try to predict an upper bound to the final file size
	uint64_t predicted_file_size_bytes = 128;// header
	for (uint8_t table_index = 0; table_index < 6; table_index++)
	{
		uint64_t park_size_bytes = EntrySizes::CalculateParkSize(K, table_index+1);
		auto final_park = *(graph_parks[table_index].end()-1);
		uint64_t num_parks = (final_park->start_pos + final_park->size() + kEntriesPerPark-1)/kEntriesPerPark;
		predicted_file_size_bytes += park_size_bytes*num_parks;
	}

    uint64_t total_C1_entries = cdiv(((*(final_parks.end()-1))->start_pos + (*(final_parks.end()-1))->size()), kCheckpoint1Interval);
    uint64_t total_C2_entries = cdiv(total_C1_entries, kCheckpoint2Interval);
    uint32_t size_C3 = EntrySizes::CalculateC3Size(K);
    predicted_file_size_bytes += (total_C1_entries + 1) * (Util::ByteAlign(K) / 8) + (total_C2_entries + 1) * (Util::ByteAlign(K) / 8) + (total_C1_entries)*size_C3;

    output_buffer = new Buffer(predicted_file_size_bytes*1.5, filename);
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
    *(uint8_t*)(output_buffer->data + output_buffer->GetInsertionOffset(1)) = K;

    *(uint16_t*)(output_buffer->data + output_buffer->GetInsertionOffset(2)) = bswap_16(kFormatDescription.size());
    output_buffer->InsertString(kFormatDescription);
    *(uint16_t*)(output_buffer->data + output_buffer->GetInsertionOffset(2)) = bswap_16(memo_size);
    output_buffer->InsertData((void*)memo, memo_size);

    uint8_t pointers[12 * 8];
    uint32_t pointers_offset = output_buffer->InsertData(pointers, sizeof(pointers));
    final_table_begin_pointers = (uint64_t*)(output_buffer->data + pointers_offset);

    phase3DoTable<0>();
    phase3DoTable<1>();
    phase3DoTable<2>();
    phase3DoTable<3>();
    phase3DoTable<4>();
    phase3DoTable<5>();
	final_table_begin_pointers[6] = bswap_64(*(output_buffer->insert_pos));

    cout << "Phase 3 finished in " << time(nullptr) - phase_start_seconds << "s" << endl;
}

#include "explicit_templates.hpp"