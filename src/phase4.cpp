#include <vector>

#include "encoding.hpp"
#include "entry_sizes.hpp"
#include "plotter.hpp"
#include "util.hpp"

using namespace std;

template <uint8_t K>
void Plotter<K>::phase4()
{
    uint64_t phase_start_seconds = time(nullptr);
    uint32_t P7_park_size = Util::ByteAlign((K + 1) * kEntriesPerPark) / 8;
    uint64_t num_entries = (*(final_parks.end()-1))->start_pos + (*(final_parks.end()-1))->size();
    uint64_t number_of_p7_parks = (num_entries + kEntriesPerPark-1)/kEntriesPerPark;

    uint64_t begin_byte_C1 = bswap_64(final_table_begin_pointers[6]) + number_of_p7_parks * P7_park_size;

    uint64_t total_C1_entries = cdiv(num_entries, kCheckpoint1Interval);
    uint64_t begin_byte_C2 = begin_byte_C1 + (total_C1_entries + 1) * (Util::ByteAlign(K) / 8);
    uint64_t total_C2_entries = cdiv(total_C1_entries, kCheckpoint2Interval);
    uint64_t begin_byte_C3 = begin_byte_C2 + (total_C2_entries + 1) * (Util::ByteAlign(K) / 8);

    uint32_t size_C3 = EntrySizes::CalculateC3Size(K);
    uint64_t end_byte = begin_byte_C3 + (total_C1_entries)*size_C3;

    final_table_begin_pointers[7] = bswap_64(begin_byte_C1);
    final_table_begin_pointers[8] = bswap_64(begin_byte_C2);
    final_table_begin_pointers[9] = bswap_64(begin_byte_C3);
    final_table_begin_pointers[10] = bswap_64(end_byte);

    uint64_t plot_file_reader = 0;
    uint64_t final_file_writer_1 = begin_byte_C1;
    uint64_t final_file_writer_2 = begin_byte_C3;
    uint64_t final_file_writer_3 = bswap_64(final_table_begin_pointers[6]);

    uint64_t prev_y = 0;
    std::vector<Bits> C2;
    uint64_t num_C1_entries = 0;
    std::vector<uint8_t> deltas_to_write;
    //uint32_t right_entry_size_bytes = phase1_table7->entry_len;

    uint8_t *right_entry_buf;
    uint8_t C1_entry_buf[16];
    uint8_t* C3_entry_buf = (uint8_t*)malloc(size_C3);
    uint8_t* P7_entry_buf = (uint8_t*)malloc(P7_park_size);
    assert(C1_entry_buf);
    assert(C3_entry_buf);
    assert(P7_entry_buf);

    std::cout << "\tStarting to write C1 and C3 tables" << std::endl;

    vector<uint128_t> line_points;

    uint64_t f7_position = 0;
    ParkBits to_write_p7;
    // We read each table7 entry, which is sorted by f7, but we don't need f7 anymore. Instead,
    // we will just store pos6, and the deltas in table C3, and checkpoints in tables C1 and C2.
    for (auto & park : final_parks)
    {
        line_points.resize(park->size());
        park->readEntries(line_points);

        for (auto & line_point: line_points)
        {
            uint64_t pos = line_point&((1ULL << K)-1);
            uint64_t y = line_point>>K;

            Bits entry_y_bits = Bits(y, K);

            if (f7_position % kEntriesPerPark == 0 && f7_position > 0) {
                memset(P7_entry_buf, 0, P7_park_size);
                assert(to_write_p7.GetSize()/8 <= P7_park_size);
                to_write_p7.ToBytes(P7_entry_buf);
                memcpy(output_buffer->data + final_file_writer_3, (P7_entry_buf), P7_park_size);
                final_file_writer_3 += P7_park_size;
                to_write_p7 = ParkBits();
            }

            to_write_p7 += ParkBits(pos, K + 1);

            if (f7_position % kCheckpoint1Interval == 0) {
                assert(entry_y_bits.GetSize() <= Util::ByteAlign(K));
                entry_y_bits.ToBytes(C1_entry_buf);
                assert(final_file_writer_1 < output_buffer->data_len);
                memcpy(output_buffer->data + final_file_writer_1, (C1_entry_buf), Util::ByteAlign(K) / 8);
                final_file_writer_1 += Util::ByteAlign(K) / 8;
                if (num_C1_entries > 0) {
                    final_file_writer_2 = begin_byte_C3 + (num_C1_entries - 1) * size_C3;
                    size_t num_bytes =
                            Encoding::ANSEncodeDeltas(deltas_to_write, deltas_to_write.size(), kC3R, C3_entry_buf + 2) + 2;

                    // We need to be careful because deltas are variable sized, and they need to fit
                    assert(size_C3 > num_bytes);

                    // Write the size
                    Util::IntToTwoBytes(C3_entry_buf, num_bytes - 2);

                    memcpy(output_buffer->data + final_file_writer_2, (C3_entry_buf), num_bytes);
                    final_file_writer_2 += num_bytes;
                }
                prev_y = y;
                if (f7_position % (kCheckpoint1Interval * kCheckpoint2Interval) == 0) {
                    C2.emplace_back(std::move(entry_y_bits));
                }
                deltas_to_write.clear();
                ++num_C1_entries;
            } else {
                assert(prev_y <= y);
                if (y == prev_y) {
                    deltas_to_write.push_back(0);
                } else {
                    deltas_to_write.push_back(y - prev_y);
                }
                prev_y = y;
            }

            f7_position++;
        }
    }
    Encoding::ANSFree(kC3R);
    //res.table7_sm.reset();


    // Writes the final park to disk
    assert(to_write_p7.GetSize()/8 < P7_park_size);
    memset(P7_entry_buf, 0, P7_park_size);
    to_write_p7.ToBytes(P7_entry_buf);

    memcpy(output_buffer->data + final_file_writer_3, (P7_entry_buf), P7_park_size);
    final_file_writer_3 += P7_park_size;

    if (!deltas_to_write.empty()) {
        size_t num_bytes = Encoding::ANSEncodeDeltas(deltas_to_write, deltas_to_write.size(), kC3R, C3_entry_buf + 2);
        assert((num_bytes + 2) <= size_C3);
        memset(C3_entry_buf + num_bytes + 2, 0, size_C3 - (num_bytes + 2));
        final_file_writer_2 = begin_byte_C3 + (num_C1_entries - 1) * size_C3;

        // Write the size
        Util::IntToTwoBytes(C3_entry_buf, num_bytes);

        memcpy(output_buffer->data + final_file_writer_2, (C3_entry_buf), size_C3);
        final_file_writer_2 += size_C3;
        Encoding::ANSFree(kC3R);
    }

    Bits(0, Util::ByteAlign(K)).ToBytes(C1_entry_buf);
    memcpy(output_buffer->data + final_file_writer_1, (C1_entry_buf), Util::ByteAlign(K) / 8);
    final_file_writer_1 += Util::ByteAlign(K) / 8;
    std::cout << "\tFinished writing C1 and C3 tables" << std::endl;
    std::cout << "\tWriting C2 table" << std::endl;

    for (Bits &C2_entry : C2) {
        C2_entry.ToBytes(C1_entry_buf);
        memcpy(output_buffer->data + final_file_writer_1, (C1_entry_buf), Util::ByteAlign(K) / 8);
        final_file_writer_1 += Util::ByteAlign(K) / 8;
    }
    Bits(0, Util::ByteAlign(K)).ToBytes(C1_entry_buf);
    memcpy(output_buffer->data + final_file_writer_1, (C1_entry_buf), Util::ByteAlign(K) / 8);
    final_file_writer_1 += Util::ByteAlign(K) / 8;
    std::cout << "\tFinished writing C2 table" << std::endl;

    free(C3_entry_buf);
    free(P7_entry_buf);


    /*
    final_file_writer_1 = res.header_size - 8 * 3;
    uint8_t table_pointer_bytes[8];

    // Writes the pointers to the start of the tables, for proving
    for (int i = 8; i <= 10; i++) {
        Util::IntToEightBytes(table_pointer_bytes, res.final_table_begin_pointers[i]);
        tmp2_disk.Write(final_file_writer_1, table_pointer_bytes, 8);
        final_file_writer_1 += 8;
    }
    */
    delete buffers[6];

    std::cout << "\tFinal table pointers:" << std::endl << std::hex;

    for (int i = 1; i <= 10; i++) {
        std::cout << "\t" << (i < 8 ? "P" : "C") << (i < 8 ? i : i - 7);
        std::cout << ": 0x" << bswap_64(final_table_begin_pointers[i]) << std::endl;
    }
    std::cout << std::dec;

    output_buffer->Truncate(end_byte);
    delete output_buffer;

    cout << "Phase 4 finished in " << time(nullptr) - phase_start_seconds << "s" << endl;
}

#include "explicit_templates.hpp"