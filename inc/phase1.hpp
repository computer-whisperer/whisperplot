#ifndef PHASE1_HPP
#define PHASE1_HPP

#include <vector>
#include "bucket_page_index.hpp"
#include "buffer.hpp"
#include "park.hpp"
#include "phase1.hpp"

template <uint8_t K>
class Phase1
{
public:

    static constexpr uint8_t line_point_delta_len_bits = K + 10;
    static constexpr uint8_t finaltable_y_delta_len_bits = 10;

    std::vector<uint32_t> cpu_ids;

    std::vector<Buffer*> buffers;
    std::vector<std::vector<TemporaryPark<line_point_delta_len_bits>*>> graph_parks;
    std::vector<TemporaryPark<finaltable_y_delta_len_bits>*> final_parks;

    Phase1(const uint8_t* id_in, uint32_t num_threads_in);

private:
    uint32_t num_threads;

    uint32_t d_new_entry_positions[(1ULL<<(K+2))];

    static void ThreadA(
            uint32_t cpu_id,
            std::atomic<uint64_t>* coordinator,
            const uint8_t* id,
            BucketPageIndex<YCBucketEntry<K, -1>>* new_bucket_index
    );

    template <int8_t table_index>
    static void ThreadB(
            uint32_t cpu_id,
            std::atomic<uint64_t>* coordinator,
            std::map<uint32_t, std::vector<uint32_t>>* new_entry_positions,
            std::map<uint32_t, BucketPageIndex<YCBucketEntry<K, table_index-1>>*>* prev_bucket_indexes,
            BucketPageIndex<YCBucketEntry<K, table_index>>* new_bucket_index,
            BucketPageIndex<LinePointEntryUIDBucketEntry<K>>* new_line_point_bucket_index);

    template <int8_t table_index>
    static void ThreadC(
            uint32_t cpu_id,
            std::atomic<uint64_t>* coordinator,
            std::map<uint32_t, std::vector<uint32_t>> * new_entry_positions,
            std::map<uint32_t, BucketPageIndex<LinePointEntryUIDBucketEntry<K>>*> line_point_bucket_index,
            std::vector<TemporaryPark<line_point_delta_len_bits>*>* parks,
            Buffer* buffer
    );

    static void ThreadD(
            uint32_t cpu_id,
            std::atomic<uint64_t>* coordinator,
            std::map<uint32_t, std::vector<uint32_t>> * new_entry_positions,
            std::map<uint32_t, BucketPageIndex<YCBucketEntry<K, 5>>*> line_point_bucket_indexes,
            std::vector<TemporaryPark<finaltable_y_delta_len_bits>*>* parks,
            Buffer* buffer
    );

    template <int8_t table_index>
    vector<BucketPageIndex<YCBucketEntry<K, table_index>>*> DoTable(
            vector<BucketPageIndex<YCBucketEntry<K, table_index-1>>*> prev_bucket_index);

};



#endif
