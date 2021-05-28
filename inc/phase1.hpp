#ifndef PHASE1_HPP
#define PHASE1_HPP

#include <vector>
#include "bucket_page_index.hpp"
#include "buffers.hpp"
#include "park.hpp"
#include "phase1.hpp"

template <uint8_t K>
class Phase1
{

    static void ThreadA(
            std::atomic<uint64_t>& coordinator,
            const uint8_t* id,
            BucketPageIndex<YCBucketEntry<K, -1>>& new_bucket_index
    );
    template <int8_t table_index>
    static void ThreadB(
            std::atomic<uint64_t>& coordinator,
            std::vector<uint32_t>& new_prev_bucket_positions,
            BucketPageIndex<YCBucketEntry<K, table_index-1>>& prev_bucket_index,
            BucketPageIndex<YCBucketEntry<K, table_index>>& new_bucket_index,
            BucketPageIndex<LinepointEntryUIDBucketEntry<K>>& new_linepoint_bucket_index);

    template <int8_t table_index>
    static void ThreadC(
            std::atomic<uint64_t>& coordinator,
            std::vector<uint32_t>& new_entry_positions,
            BucketPageIndex<LinepointEntryUIDBucketEntry<K>>& linepoint_bucket_index,
            Buffer* output_buffer,
            std::vector<TemporaryPark<K+2>*>& temporary_parks_out
    );

    static void ThreadD(
            std::atomic<uint64_t>& coordinator,
            std::vector<uint32_t>& new_entry_positions,
            BucketPageIndex<YCBucketEntry<K, 6>>& linepoint_bucket_index,
            Buffer* output_buffer,
            std::vector<TemporaryPark<finaltable_y_delta_len_bits>*>& temporary_parks_out
    );

    template <int8_t table_index> static void Phase1DoTable(
            uint32_t num_threads,
            void * prev_bucket_index_in,
            void ** next_bucket_index_out,
            std::vector<uint32_t>* prev_remap_table,
            std::vector<uint32_t>** next_remap_table,
            Buffer** new_buffer_out,
            std::vector<TemporaryPark<K+2>*>& temporary_parks_out,
    );

public:

    const uint8_t finaltable_y_delta_len_bits = 4;

    struct Phase1_Results
    {
        std::vector<Buffer*> buffers;
        std::vector<vector<TemporaryPark<K+4>*>> graph_parks;
        std::vector<TemporaryPark<finaltable_y_delta_len_bits>*> final_parks;
    }

    static Phase1_Results Run(const uint8_t* id, uint32_t num_threads);
};



#endif
