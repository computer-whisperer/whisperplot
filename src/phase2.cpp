

#include <atomic>
#include <thread>
#include <vector>

#include "encoding.hpp"
#include "buffer.hpp"
#include "util.hpp"
#include "thread_mgr.hpp"
#include "plotter.hpp"

using namespace std;


template <uint8_t K>
template <int8_t table_index>
void Plotter<K>::phase2ThreadA(
        uint32_t cpu_id,
        std::atomic<uint64_t>* coordinator,
        vector<DeltaPark<line_point_delta_len_bits>*>* parks,
        AtomicPackedArray<BooleanPackedEntry, max_entries_per_graph_table>* current_entries_used,
        AtomicPackedArray<BooleanPackedEntry, max_entries_per_graph_table>* prev_entries_used
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
		    if (current_entries_used->read(parks->at(park_id)->start_pos + i).value)
            {
                auto res = Encoding::LinePointToSquare(line_points[i]);
                prev_entries_used->set(res.first, BooleanPackedEntry(true));
                prev_entries_used->set(res.second, BooleanPackedEntry(true));
            }
        }
	}
}

template <uint8_t K>
void Plotter<K>::phase2ThreadB(
        AtomicPackedArray<BooleanPackedEntry, max_entries_per_graph_table>* entries_used,
        AtomicPackedArray<SimplePackedEntry<uint64_t,K>, max_entries_per_graph_table>* final_positions
)
{
    uint64_t new_pos = 0;
    for (uint64_t i = 0; i < entries_used->size(); i++)
    {
        final_positions->append(SimplePackedEntry<uint64_t,K>(new_pos));
        new_pos += entries_used->read(i).value;
    }
}


template <uint8_t K>
template <int8_t table_index>
void Plotter<K>::phase2DoTable()
{
    entries_used[table_index-1].fill(BooleanPackedEntry(false));
    cout << "Part A"<< (uint32_t)table_index;
    uint64_t start_seconds = time(nullptr);
    vector<thread> threads;
    std::atomic<uint64_t> coordinator = 0;
    for (uint32_t i = 0; i < cpu_ids.size(); i++)
    {
        threads.push_back(thread(
                Plotter<K>::phase2ThreadA<table_index>,
                cpu_ids[i],
                &coordinator,
                &(graph_parks[table_index]),
                &(entries_used[table_index]),
                &(entries_used[table_index-1])));
    }

    for (auto &it: threads)
    {
        it.join();
    }
    cout << " (" << time(nullptr) - start_seconds << "s)" << endl;

    // Start final position sum thread
    phase2b_threads[table_index-1] = thread(
            Plotter<K>::phase2ThreadB,
            &(entries_used[table_index-1]),
            &(final_positions[table_index-1]));
}


template <uint8_t K>
void Plotter<K>::phase2()
{
    final_positions.resize(6);
    phase2b_threads.resize(6);
    // Start final position sum thread
    phase2b_threads[5] = thread(
            Plotter<K>::phase2ThreadB,
            &(entries_used[5]),
            &(final_positions[5]));
    uint64_t phase_start_seconds = time(nullptr);
    phase2DoTable<5>();
    phase2DoTable<4>();
    phase2DoTable<3>();
    phase2DoTable<2>();
    phase2DoTable<1>();
    cout << "Phase 2 finished in " << time(nullptr) - phase_start_seconds << "s" << endl;
}

#include "explicit_templates.hpp"