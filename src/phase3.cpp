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
        vector<TemporaryPark<line_point_delta_len_bits>*>* temporary_parks,
        AtomicPackedArray<BooleanPackedEntry, max_entries_per_graph_table>* entries_used,
        AtomicPackedArray<SimplePackedEntry<uint64_t,K>, max_entries_per_graph_table>* final_positions,
        Buffer* output_buffer)
{
    PinToCpuid(cpu_id);
    vector<uint128_t> temporary_park_line_points(LinePointEntryUIDPackedEntry<K>::max_entries_per_sort_row);
    vector<uint128_t> park_line_points(kEntriesPerPark);
    vector<uint64_t> park_stubs(kEntriesPerPark);
    vector<uint8_t> park_deltas(kEntriesPerPark);
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
                temporary_parks->at(temporary_park_id)->readEntries(temporary_park_line_points.data());
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

        // Since we have approx 2^k line_points between 0 and 2^2k, the average
        // space between them when sorted, is k bits. Much more efficient than storing each
        // line point. This is divided into the stub and delta. The stub is the least
        // significant (k-kMinusStubs) bits, and largely random/incompressible. The small
        // delta is the rest, which can be efficiently encoded since it's usually very
        // small.

		for (uint32_t i = 1; i < park_entry_i; i++)
		{
			uint128_t big_delta = park_line_points[i] - park_line_points[i-1];
            uint64_t stub = big_delta & ((1ULL << (K - kStubMinusBits)) - 1);
            uint64_t small_delta = big_delta >> (K - kStubMinusBits);
            assert(small_delta < 256);
            park_deltas[i-1] = small_delta;
            park_deltas[i-1] = stub;
		}

	    // Parks are fixed size, so we know where to start writing. The deltas will not go over
	    // into the next park.
		uint8_t * park_start = output_buffer->data + *(output_buffer->insert_pos) + park_id * park_size_bytes;
	    uint8_t * dest = park_start;

	    // First entry is static
        uint128_t first_line_point = park_line_points[0] << (128 - 2 * K);
        Util::IntTo16Bytes(dest, first_line_point);
        dest += EntrySizes::CalculateLinePointSize(K);

	    // We use ParkBits instead of Bits since it allows storing more data
	    ParkBits park_stubs_bits;
        for (uint32_t i = 0; i < park_entry_i-1; i++)
        {
	        park_stubs_bits.AppendValue(park_stubs[i], (K - kStubMinusBits));
	    }
	    uint32_t stubs_size = EntrySizes::CalculateStubsSize(K);
	    uint32_t stubs_valid_size = cdiv(park_stubs_bits.GetSize(), 8);
	    park_stubs_bits.ToBytes(dest);
	    memset(dest + stubs_valid_size, 0, stubs_size - stubs_valid_size);
	    dest += stubs_size;

	    // The stubs are random so they don't need encoding. But deltas are more likely to
	    // be small, so we can compress them
	    double R = kRValues[table_index];
	    uint8_t *deltas_start = dest + 2;
	    size_t deltas_size = Encoding::ANSEncodeDeltas(park_deltas, park_entry_i-1, R, deltas_start);

	    if (!deltas_size) {
	        // Uncompressed
	        deltas_size = park_entry_i-1;
	        Util::IntToTwoBytesLE(dest, deltas_size | 0x8000);
	        memcpy(deltas_start, park_deltas.data(), deltas_size);
	    } else {
	        // Compressed
	        Util::IntToTwoBytesLE(dest, deltas_size);
	    }

	    dest += 2 + deltas_size;

        memset(dest, 0x00, park_size_bytes - (dest - park_start));

	    assert(park_size_bytes > (uint64_t)(dest - park_start));
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
                &(final_positions[table_index+1]),
                output_buffer));
    }

    for (auto &it: threads)
    {
        it.join();
    }
    cout << " (" << time(nullptr) - start_seconds << "s)" << endl;


    // setup output_table_offset for next iteration
    uint64_t park_size_bytes = EntrySizes::CalculateParkSize(K, table_index+1);
    uint64_t num_entries = (graph_parks.end()-1)->start_pos + (graph_parks.end()-1)->size();
    uint64_t parks_needed = (num_entries+kEntriesPerPark-1)/kEntriesPerPark;
    output_buffer->GetInsertionOffset(parks_needed*park_size_bytes);

    delete buffers[table_index];
}

template<uint8_t K>
void Plotter<K>::phase3()
{
	// Try to predict an upper bound to the final file size
	uint64_t predicted_file_size_bytes = 128;// header
	for (uint8_t table_index = 0; table_index < 6; table_index++)
	{
		uint64_t park_size_bytes = EntrySizes::CalculateParkSize(K, table_index+1);
		auto final_park = (graph_parks[table_index]->end()-1);
		uint64_t num_parks = (final_park->start_pos + final_park->size())/kEntriesPerPark;
		predicted_file_size_bytes += park_size_bytes*num_parks;
	}

    uint64_t total_C1_entries = cdiv(((final_parks.end()-1)->start_pos + (final_parks.end()-1)->size()), kCheckpoint1Interval);
    uint64_t total_C2_entries = cdiv(total_C1_entries, kCheckpoint2Interval);
    uint32_t size_C3 = EntrySizes::CalculateC3Size(K);
    predicted_file_size_bytes += (total_C1_entries + 1) * (Util::ByteAlign(K) / 8) + (total_C2_entries + 1) * (Util::ByteAlign(K) / 8) + (total_C1_entries)*size_C3;

    output_buffer = new Buffer(predicted_file_size_bytes*1.2, filename);

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
}
