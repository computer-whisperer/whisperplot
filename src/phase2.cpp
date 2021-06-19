

#include <atomic>
#include <thread>
#include <vector>

#include "encoding.hpp"
#include "buffer.hpp"
#include "chiapos_util.hpp"
#include "thread_mgr.hpp"
#include "plotter.hpp"
#include "status_update.hpp"

using namespace std;


template <PlotConf conf>
template <int8_t table_index>
void Plotter<conf>::phase2ThreadA(
        uint32_t cpu_id,
        std::atomic<uint64_t>* coordinator,
        vector<Park*>* parks,
        uint8_t* current_entries_used,
        uint8_t* prev_entries_used
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
		    if ((table_index == 5) || (current_entries_used[parks->at(park_id)->start_pos + i]))
            {
                auto res = Encoding::LinePointToSquare(line_points[i]);
                prev_entries_used[res.first] = 1;
                prev_entries_used[res.second] = 1;
            }
        }
	}
}

template <PlotConf conf>
void Plotter<conf>::phase2ThreadB(
        uint8_t* entries_used,
        p2_final_positions_type* final_positions
)
{
    uint64_t new_pos = 0;
    for (uint64_t i = 0; i < fwd_max_entries_per_graph_table; i++)
    {
        final_positions->set(i, PackedEntry<1, 1ULL << (conf.K+1), 1>(new_pos));
        new_pos += entries_used[i];
    }
}


template <PlotConf conf>
template <int8_t table_index>
void Plotter<conf>::phase2DoTable()
{
    StatusUpdate::StartSeg("1." + to_string((uint32_t)table_index) + ".0S");

    entries_used[table_index-1] = (uint8_t*)malloc(fwd_max_entries_per_graph_table);
    memset(entries_used[table_index-1], 0, fwd_max_entries_per_graph_table);

    uint64_t start_seconds = time(nullptr);

    StatusUpdate::StartSeg("1." + to_string((uint32_t)table_index) + ".1M");

    vector<thread> threads;
    uint8_t * current_entries_used = nullptr;
    if(table_index != 5)
    {
        current_entries_used = entries_used[table_index];
    }
    std::atomic<uint64_t> coordinator = 0;
    for (auto & cpu_id : cpu_ids)
    {
        threads.push_back(thread(
                Plotter<conf>::phase2ThreadA<table_index>,
                cpu_id,
                &coordinator,
                &(phase1_graph_parks[table_index]),
                current_entries_used,
                entries_used[table_index-1]));
    }

    for (auto &it: threads)
    {
        it.join();
    }

    // Start final position sum thread
    phase2b_threads[table_index-1] = thread(
            Plotter<conf>::phase2ThreadB,
            entries_used[table_index-1],
            &(phase2_final_positions->at(table_index - 1)));
}


template <PlotConf conf>
void Plotter<conf>::phase2()
{
    StatusUpdate::StartSeg("1.-1.0S");
    entries_used.resize(5);
    phase2_final_positions = new vector<p2_final_positions_type>(5);
    phase2b_threads.resize(5);
    // Start final position sum thread
    uint64_t phase_start_seconds = time(nullptr);
    phase2DoTable<5>();
    phase2DoTable<4>();
    phase2DoTable<3>();
    phase2DoTable<2>();
    phase2DoTable<1>();
    StatusUpdate::EndSeg();
    cout << "Phase 2 finished in " << time(nullptr) - phase_start_seconds << "s" << endl;
}

#include "explicit_templates.hpp"
#include "packed_array.hpp"
