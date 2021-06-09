

#include <atomic>
#include <thread>
#include <vector>

#include "encoding.hpp"
#include "buffer.hpp"
#include "chiapos_util.hpp"
#include "thread_mgr.hpp"
#include "plotter.hpp"

using namespace std;


template <uint8_t K, uint32_t num_rows>
template <int8_t table_index>
void Plotter<K, num_rows>::phase2ThreadA(
        uint32_t cpu_id,
        std::atomic<uint64_t>* coordinator,
        vector<DeltaPark<line_point_delta_len_bits>*>* parks,
        entries_used_type* current_entries_used,
        entries_used_type* prev_entries_used
        )
{
    PinToCpuid(cpu_id);
	while (true)
	{
		uint32_t park_id = (*coordinator)++;

		if (park_id >= parks->size())
        {
		    break;
        }

		vector<uint128_t> line_points(parks->at(park_id)->size());
		parks->at(park_id)->readEntries(line_points);

		for (uint64_t i = 0; i < line_points.size(); i++)
        {
		    if ((table_index == 5) || (current_entries_used->get(parks->at(park_id)->start_pos + i).getY()))
            {
                auto res = Encoding::LinePointToSquare(line_points[i]);
                prev_entries_used->set(res.first, BooleanPackedEntry(true));
                prev_entries_used->set(res.second, BooleanPackedEntry(true));
            }
        }
	}
}

template <uint8_t K, uint32_t num_rows>
void Plotter<K, num_rows>::phase2ThreadB(
        entries_used_type* entries_used,
        p2_final_positions_type* final_positions
)
{
    uint64_t new_pos = 0;
    for (uint64_t i = 0; i < max_entries_per_graph_table; i++)
    {
        final_positions->set(i, PackedEntry<1, 1ULL << (K+1), 1>(new_pos));
        new_pos += entries_used->get(i).getY();
    }
}


template <uint8_t K, uint32_t num_rows>
template <int8_t table_index>
void Plotter<K, num_rows>::phase2DoTable()
{
    entries_used->at(table_index-1).fill(BooleanPackedEntry(false));
    cout << "Part A"<< (uint32_t)table_index;
    uint64_t start_seconds = time(nullptr);
    vector<thread> threads;
    entries_used_type * current_entries_used = nullptr;
    if(table_index != 5)
    {
        current_entries_used = &(entries_used->at(table_index));
    }
    std::atomic<uint64_t> coordinator = 0;
    for (auto & cpu_id : cpu_ids)
    {
        threads.push_back(thread(
                Plotter<K, num_rows>::phase2ThreadA<table_index>,
                cpu_id,
                &coordinator,
                &(phase1_graph_parks[table_index]),
                current_entries_used,
                &(entries_used->at(table_index-1))));
    }

    for (auto &it: threads)
    {
        it.join();
    }
    cout << " (" << time(nullptr) - start_seconds << "s)" << endl;

    // Start final position sum thread
    phase2b_threads[table_index-1] = thread(
            Plotter<K, num_rows>::phase2ThreadB,
            &(entries_used->at(table_index-1)),
            &(phase2_final_positions->at(table_index - 1)));
}


template <uint8_t K, uint32_t num_rows>
void Plotter<K, num_rows>::phase2()
{
    entries_used = new vector<entries_used_type>(5);
    phase2_final_positions = new vector<p2_final_positions_type>(5);
    phase2b_threads.resize(5);
    // Start final position sum thread
    uint64_t phase_start_seconds = time(nullptr);
    phase2DoTable<5>();
    phase2DoTable<4>();
    phase2DoTable<3>();
    phase2DoTable<2>();
    phase2DoTable<1>();
    cout << "Phase 2 finished in " << time(nullptr) - phase_start_seconds << "s" << endl;
}

#include "explicit_templates.hpp"
#include "packed_array.hpp"
