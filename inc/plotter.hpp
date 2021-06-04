#ifndef PHASE1_HPP
#define PHASE1_HPP

#include <vector>
#include "penguin.hpp"
#include "buffer.hpp"
#include "park.hpp"
#include "plotter.hpp"

template <uint8_t K>
class Plotter
{
public:

    static constexpr uint8_t line_point_delta_len_bits = K + 10;
    static constexpr uint8_t finaltable_y_delta_len_bits = K + 10;
    static constexpr uint64_t max_entries_per_graph_table = (1ULL << K)*1.1;



    std::vector<Buffer*> buffers;
    std::vector<std::vector<DeltaPark<line_point_delta_len_bits>*>> phase1_graph_parks;
    std::vector<DeltaPark<finaltable_y_delta_len_bits>*> phase1_final_parks;
    std::string filename;

    inline Plotter(const uint8_t* id_in, const uint8_t* memo_in, uint32_t memo_size_in, std::vector<uint32_t> cpu_ids_in, std::string filename_in)
    {
        id = id_in;
        cpu_ids = cpu_ids_in;
        filename = filename_in;
        memo = memo_in;
        memo_size = memo_size_in;
    }

    void phase1();
    void phase2();
    void phase3();
    void phase4();

    void check();

private:
    Buffer* output_buffer;
    const uint8_t* id;
    const uint8_t* memo;
    uint32_t memo_size;
    std::vector<std::thread> phase2b_threads;
    std::vector<uint32_t> cpu_ids;
    std::map<uint32_t, std::vector<uint64_t>> d_new_entry_positions;
    std::vector<AtomicPackedArray<BooleanPackedEntry, max_entries_per_graph_table>> entries_used;
    std::vector<AtomicPackedArray<SimplePackedEntry<uint64_t,K>, max_entries_per_graph_table>> final_positions;
    uint64_t * final_table_begin_pointers;
    std::vector<std::vector<Park*>> final_parks;

    static void phase1ThreadA(
            uint32_t cpu_id,
            std::atomic<uint64_t>* coordinator,
            const uint8_t* id,
            Penguin<YCPackedEntry<K, -1>>* new_penguin
    );

    template <int8_t table_index>
    static void phase1ThreadB(
            uint32_t cpu_id,
            std::atomic<uint64_t>* coordinator,
            std::map<uint32_t, std::vector<uint64_t>>* new_entry_positions,
            std::map<uint32_t, Penguin<YCPackedEntry<K, table_index - 1>>*>* prev_penguins,
            Penguin<YCPackedEntry<K, table_index>>* new_yc_penguin,
            Penguin<LinePointEntryUIDPackedEntry<K>>* new_line_point_penguin);

    template <int8_t table_index>
    static void phase1ThreadC(
            uint32_t cpu_id,
            std::atomic<uint64_t>* coordinator,
            std::map<uint32_t, std::vector<uint64_t>> * new_entry_positions,
            std::map<uint32_t, Penguin<LinePointEntryUIDPackedEntry<K>>*> line_point_bucket_index,
            std::vector<DeltaPark<line_point_delta_len_bits> *> *parks,
            Buffer* buffer
    );

    static void phase1ThreadD(
            uint32_t cpu_id,
            std::atomic<uint64_t>* coordinator,
            std::map<uint32_t, std::vector<uint64_t>> * new_entry_positions,
            std::map<uint32_t, Penguin<YCPackedEntry<K, 5>>*> line_point_penguins,
            std::vector<DeltaPark<finaltable_y_delta_len_bits> *> *parks,
            Buffer* buffer
    );

    template <int8_t table_index>
    std::map<uint32_t, Penguin<YCPackedEntry<K, table_index>>*> phase1DoTable(
            std::map<uint32_t, Penguin<YCPackedEntry<K, table_index - 1>>*> prev_bucket_indexes);

    template <int8_t table_index>
    static void phase2ThreadA(
            uint32_t cpu_id,
            std::atomic<uint64_t>* coordinator,
            std::vector<DeltaPark<line_point_delta_len_bits>*>* parks,
            AtomicPackedArray<BooleanPackedEntry, max_entries_per_graph_table>* current_entries_used,
            AtomicPackedArray<BooleanPackedEntry, max_entries_per_graph_table>* prev_entries_used
    );

    static void phase2ThreadB(
            AtomicPackedArray<BooleanPackedEntry, max_entries_per_graph_table>* entries_used,
            AtomicPackedArray<SimplePackedEntry<uint64_t,K>, max_entries_per_graph_table>* final_positions);

    template <int8_t table_index>
    void phase2DoTable();

    template <int8_t table_index>
    static void phase3ThreadA(
            uint32_t cpu_id,
            std::atomic<uint64_t>* coordinator,
            std::vector<DeltaPark<line_point_delta_len_bits>*>* temporary_parks,
            AtomicPackedArray<BooleanPackedEntry, max_entries_per_graph_table>* entries_used,
            std::vector<AtomicPackedArray<SimplePackedEntry<uint64_t,K>, max_entries_per_graph_table>>* final_positions,
            Buffer* output_buffer,
            std::vector<std::vector<Park*>>* final_parks,
            uint64_t start_offset);

    template <int8_t table_index>
    void phase3DoTable();
};



#endif
