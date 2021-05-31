#ifndef PHASE1_HPP
#define PHASE1_HPP

#include <vector>
#include "penguin.hpp"
#include "buffer.hpp"
#include "park.hpp"
#include "phase1.hpp"

template <uint8_t K>
class Phase1
{
public:

    static constexpr uint8_t line_point_delta_len_bits = K + 10;
    static constexpr uint8_t finaltable_y_delta_len_bits = K + 10;



    std::vector<Buffer*> buffers;
    std::vector<std::vector<TemporaryPark<line_point_delta_len_bits>*>> graph_parks;
    std::vector<TemporaryPark<finaltable_y_delta_len_bits>*> final_parks;

    Phase1(const uint8_t* id_in, std::vector<uint32_t> num_threads_in);

private:
    std::vector<uint32_t> cpu_ids;
    std::map<uint32_t, std::vector<uint32_t>> d_new_entry_positions;

    static void ThreadA(
            uint32_t cpu_id,
            std::atomic<uint64_t>* coordinator,
            const uint8_t* id,
            Penguin<YCPackedEntry<K, -1>>* new_penguin
    );

    template <int8_t table_index>
    static void ThreadB(
            uint32_t cpu_id,
            std::atomic<uint64_t>* coordinator,
            std::map<uint32_t, std::vector<uint32_t>>* new_entry_positions,
            std::map<uint32_t, Penguin<YCPackedEntry<K, table_index - 1>>*>* prev_penguins,
            Penguin<YCPackedEntry<K, table_index>>* new_yc_penguin,
            Penguin<LinePointEntryUIDPackedEntry<K>>* new_line_point_penguin);

    template <int8_t table_index>
    static void ThreadC(
            uint32_t cpu_id,
            std::atomic<uint64_t>* coordinator,
            std::map<uint32_t, std::vector<uint32_t>> * new_entry_positions,
            std::map<uint32_t, Penguin<LinePointEntryUIDPackedEntry<K>>*> line_point_bucket_index,
            std::vector<TemporaryPark<line_point_delta_len_bits>*>* parks,
            Buffer* buffer
    );

    static void ThreadD(
            uint32_t cpu_id,
            std::atomic<uint64_t>* coordinator,
            std::map<uint32_t, std::vector<uint32_t>> * new_entry_positions,
            std::map<uint32_t, Penguin<YCPackedEntry<K, 5>>*> line_point_bucket_indexes,
            std::vector<TemporaryPark<finaltable_y_delta_len_bits>*>* parks,
            Buffer* buffer
    );

    template <int8_t table_index>
    std::map<uint32_t, Penguin<YCPackedEntry<K, table_index>>*> DoTable(
            std::map<uint32_t, Penguin<YCPackedEntry<K, table_index - 1>>*> prev_bucket_indexes);

};



#endif
