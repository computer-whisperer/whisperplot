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

    static constexpr uint8_t line_point_delta_len_bits = K + 2;
    static constexpr uint8_t finaltable_y_delta_len_bits = 4;

    std::vector<Buffer*> buffers;
    std::vector<std::vector<TemporaryPark<line_point_delta_len_bits>*>> graph_parks;
    std::vector<TemporaryPark<finaltable_y_delta_len_bits>*> final_parks;

    Phase1(const uint8_t* id_in, uint32_t num_threads_in);

private:
    std::atomic<uint64_t> coordinator;
    uint32_t num_threads;
    const uint8_t* id;
    std::vector<uint32_t> new_entry_positions;

    void ThreadA(
            BucketPageIndex<YCBucketEntry<K, -1>>* new_bucket_index
    );
    template <int8_t table_index>
    void ThreadB(
            BucketPageIndex<YCBucketEntry<K, table_index-1>>* prev_bucket_index,
            BucketPageIndex<YCBucketEntry<K, table_index>>* new_bucket_index,
            BucketPageIndex<LinePointEntryUIDBucketEntry<K>>* new_line_point_bucket_index);

    template <int8_t table_index>
    void ThreadC(
            BucketPageIndex<LinePointEntryUIDBucketEntry<K>>* linepoint_bucket_index
    );

    void ThreadD(
            BucketPageIndex<YCBucketEntry<K, 5>>* line_point_bucket_index
    );

    template <int8_t table_index>
    BucketPageIndex<YCBucketEntry<K, table_index>>* DoTable(
            BucketPageIndex<YCBucketEntry<K, table_index-1>>* prev_bucket_index);

};



#endif
