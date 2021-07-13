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
#include "bitcopy.hpp"

using namespace std;

// C3 is compressed table of Y deltas, also kicks out GID list and checkpoint list for use later.
template <PlotConf conf>
void Plotter<conf>::Context::writeC3(
        Penguin<RevGIDEntry<conf, 6>, conf.interlace_factor>* gid_penguin_out,
        vector<uint64_t>* y_cpoints_out,
        uint8_t* table_out,
        uint64_t* num_bytes_out,
        uint64_t* num_checkpoints_out)
{
    vector<FwdYCEntry<conf, 5>*> input_penguins;

    for (auto & context : this->plotter->contexts)
    {
        input_penguins.push_back(context->forward_pass_final_yc_penguin);
    }

    vector<uint64_t> bytes_used_by_threads;
    vector<uint64_t> checkpoints_used_by_threads;

    vector<thread> threads;

    for (auto &cpu_id : cpu_ids) {
        threads.push_back(
                bytes_used_by_thread.push_back(0);
                uint64_t* bytes_used = bytes_used_by_thread.end()-1;

                checkpoints_used_by_thread.push_back(0);
                uint64_t* checkpoints_used = checkpoints_used_by_thread.end()-1;
                thread([this, cpu_id, bytes_used, checkpoints_used] {

                    SortedParkedPenguinUnloader<conf, FwdYCEntry<conf, 5>, kCheckpoint1Interval> unloader(input_penguins);
                    const uint64_t park_size_bytes = (kC3BitsPerEntry * kCheckpoint1Interval + 7)/8;

                    vector<uint64_t> deltas_to_write;
                    deltas_to_write.reserve(FwdYCEntry<conf, 5>::max_entries_per_row);

                    uint64_t batch_size = 1024;
                    while (true)
                    {
                        uint32_t park_id = this->plotter->coordinator.fetch_add(batch_size);

                        while (uint64_t i = 0; i < batch_size; i++)
                        {
                            auto park_entries = unloader.getBucketEntries(park_id);

                            y_cpoints_out[park_id] = park_entries[0].getY();
                            (*checkpoints_used)++;

                            deltas_to_write.clear();
                            for (uint32_t j = 1; j < park_entries.size(); j++)
                            {
                                deltas_to_write.push_back(park_entries[i].getY() - park_entries[i-1].getY());
                            }

                            (*y_cpoints)[park_id] = park_entries[0].getY();

                            uint8_t* dest = table_out + park_size_bytes*park_id;

                            size_t num_bytes =
                                    Encoding::ANSEncodeDeltas(deltas_to_write, deltas_to_write.size(), kC3R, dest + 2) + 2;
                            *((uint16_t*)dest) = bswap_16(num_bytes);

                            for (auto & entry : park_entries)
                            {
                                RevGIDEntry<conf, 5> new_entry;
                                new_entry.row = park_id;
                                new_entry.gid = entry.gid;
                                gid_penguin_out.add_entry(new_entry);
                            }

                            park_id++;

                            bytes_used += park_size_bytes;
                        }
                        break;
                    }
                }));
    }
    for (auto &it: threads)
    {
        it.join();
    }

    uint64_t num_bytes_total = 0;
    for (auto & it: bytes_used_by_thread)
    {
        num_bytes_total += it;
    }
    *num_bytes_out = num_bytes_total;

    uint64_t num_checkpoints_total = 0;
    for (auto & it: checkpoints_used_by_thread)
    {
        num_checkpoints_total += it;
    }
    *num_checkpoints_out = num_checkpoints_total;


}

// Checkpoints of Y values used in C3, returns bytes used
template <PlotConf conf>
uint64_t Plotter<conf>::writeC1(vector<uint64_t>* y_cpoints, uint8_t* table_out)
{
    uint64_t bytes_used = 0;
    uint8_t entry_len_bytes = ((conf.K+7)/8)*y_cpoints.size();

    memset(table_out, 0, entry_len_bytes*y_cpoints.size());

    for (auto & y : y_cpoints)
    {
        uint64_t y_big = bswap_64(y);
        bitcopy<0, entry_len_bytes, 64 - conf.K, 8>(table_out + bytes_used, (uint8_t*)&y_big);
        bytes_used += entry_len_bytes;
    }
    return bytes_used;
}

// Checkpoints of Y values used in C1, returns bytes used
template <PlotConf conf>
uint64_t Plotter<conf>::writeC2(vector<uint64_t>* y_cpoints, uint8_t* table_out)
{
    uint64_t bytes_used = 0;
    uint8_t entry_len_bytes = ((conf.K+7)/8)*y_cpoints.size();
    uint32_t num_entries = y_cpoints.size()/kCheckpoint2Interval;

    memset(table_out, 0, entry_len_bytes*num_entries);

    for (uint32_t i = 0; i < num_entries; i++)
    {
        uint64_t y_big = bswap_64(y_cpoints[i*kCheckpoint2Interval]);
        bitcopy<0, entry_len_bytes, 64 - conf.K, 8>(table_out + bytes_used, (uint8_t*)&y_big);
        bytes_used += entry_len_bytes;
    }
    return bytes_used;
}

template <PlotConf conf, uint8_t table_index>
void Plotter<conf>::Context::generate(Penguin<FwdGIDEntry<conf, table_index>, conf.interlace_factor>* gid_penguin, )
{

}

template <PlotConf conf>
void Plotter<conf>::deduplicateAndCompressPlotToFile()
{
    StatusUpdate::StartSeg("2.-1.1S");

    uint64_t phase_start_seconds = time(nullptr);

    output_buffer = new Buffer((1ULL<<conf.K)*28, filename);
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

    auto gid_penguin_6 = new Penguin<RevGIDEntry<conf, 6>, conf.interlace_factor>();

    vector<uint64_t> y_cpoints(((1ULL<<conf.K)*2)/kCheckpoint1Interval);

    // Do C3
    final_table_begin_pointers[9] = bswap_64(*(output_buffer->insert_pos));
    coordinator = 0;
    vector<thread> threads;
    vector<uint64_t> num_bytes_per_context;
    vector<uint64_t> num_checkpoints_per_context;
    for (uint32_t i = 0; i < this->contexts.size(); i++)
    {
        num_bytes_per_context.push_back(0);
        uint64_t* num_bytes = num_bytes_per_context.end()-1;

        num_checkpoints_per_context.push_back(0);
        uint64_t* num_checkpoints = num_checkpoints_per_context.end()-1;

        threads.push_back(thread( [this, i, gid_penguin_6, y_cpoints, num_bytes, num_checkpoints] {
            this->contexts[i]->template writeC3(
                    gid_penguin_6,
                    this->output_buffer.data + (*this->output_buffer->insert_pos),
                    &y_cpoints,
                    &num_bytes,
                    &num_checkpoints);
        }));
    }
    for (auto & thread : threads)
    {
        thread.join();
    }

    // Clear penguins
    for (auto & context : contexts)
    {
        delete context->forward_pass_final_yc_penguin;
    }

    uint64_t num_bytes_used = 0;
    for (auto & it : num_bytes_per_context)
    {
        num_bytes_used += it;
    }
    output_buffer->insert_pos->fetch_add(num_bytes_used);

    uint64_t num_checkpoints = 0;
    for (auto & it : num_checkpoints_per_context)
    {
        num_checkpoints += it;
    }
    y_cpoints.resize(num_checkpoints);

    // C1
    final_table_begin_pointers[7] = bswap_64(*(output_buffer->insert_pos));
    num_bytes_used = writeC1(&y_cpoints, this->output_buffer.data + (*this->output_buffer->insert_pos));
    output_buffer->insert_pos->fetch_add(num_bytes_used);

    // C2
    final_table_begin_pointers[8] = bswap_64(*(output_buffer->insert_pos));
    num_bytes_used = writeC2(&y_cpoints, this->output_buffer.data + (*this->output_buffer->insert_pos));
    output_buffer->insert_pos->fetch_add(num_bytes_used);

}

#include "explicit_templates.hpp"