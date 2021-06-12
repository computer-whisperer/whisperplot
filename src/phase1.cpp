#include <vector>
#include <atomic>
#include <cstring>
#include <ctime>

#include "penguin.hpp"
#include "buffer.hpp"
#include "calculate_bucket.hpp"
#include "encoding.hpp"
#include "park.hpp"
#include "plotter.hpp"
#include "pos_constants.hpp"
#include "thread_mgr.hpp"

using namespace std;

template <uint8_t K, uint32_t num_rows>
void Plotter<K, num_rows>::phase1ThreadA(
        uint32_t cpu_id,
        std::atomic<uint64_t>* coordinator,
        const uint8_t* id,
        Penguin<YCPackedEntry<-1>>* new_penguin
)
{
    PinToCpuid(cpu_id);
    uint64_t batch_size = (1ULL<<kBatchSizes);
    vector<uint64_t> buff(batch_size);
    F1Calculator f1(K, id);
    while (true)
    {
        uint64_t x = coordinator->fetch_add(batch_size);
        if (x >= (1ULL<<K))
            break;
        if ((x+batch_size) > (1ULL<<K))
        {
            batch_size = (1ULL<<K) - x;
        }

        f1.CalculateBuckets(x, batch_size, buff.data());

         for (uint64_t i = 0; i < batch_size; i++)
        {
            YCPackedEntry<-1> entry;
            entry.setY(buff[i]);
            entry.c = x+i;
            new_penguin->addEntry(entry);
        }
    }
}

template <uint8_t K, uint32_t num_rows>
template <int8_t table_index>
void Plotter<K, num_rows>::phase1ThreadB(
        uint32_t cpu_id,
        const uint8_t* id,
        std::atomic<uint64_t> * coordinator,
        p1_buckets_done_type<table_index-1> * bucket_left_done,
        p1_buckets_done_type<table_index-1> * bucket_right_done,
        std::map<uint32_t, phase1_new_positions_type*>* new_entry_positions,
        map<uint32_t, Penguin<YCPackedEntry<table_index - 1>> *> *prev_penguins,
        Penguin<YCPackedEntry<table_index>> * new_yc_penguin,
        Penguin<LinePointUIDPackedEntry<table_index>> * new_line_point_penguin)
{
    using old_yc = YCPackedEntry<table_index - 1>;
    using new_yc = YCPackedEntry<table_index>;
    using temp_yc = YCTempPackedEntry<table_index-1>;
    PinToCpuid(cpu_id);
    F1Calculator f1(K, id);
    vector<vector<uint32_t>> right_map(kBC);
    for (auto& i : right_map) {
        i.clear();
    }
    std::vector<uint32_t> right_map_clean;
    FxCalculator fx(K, table_index+2);

    // These are BC buckets from the chia algo
    constexpr uint32_t bc_bucket_num = (1ULL << (K+kExtraBits))/kBC;
    constexpr uint32_t temp_buckets_needed = old_yc::row_divisor/kBC + 4;
    constexpr uint32_t max_entries_per_bc_bucket = 350;

    vector<PackedArray<temp_yc, max_entries_per_bc_bucket>> temp_buckets(temp_buckets_needed);
    vector<uint32_t> temp_bucket_fill_states(temp_buckets_needed);
    vector<vector<uint64_t>> entry_positions(temp_buckets_needed);

    int64_t latest_bucket_started_load = -1;
    int64_t latest_bucket_finished_load = -1;
    int64_t latest_row_loaded = -1;

    for (auto & i : entry_positions)
    {
        i.reserve(max_entries_per_bc_bucket);
    }

    constexpr uint32_t batchSize = temp_buckets_needed*32;

    while (true)
    {
        uint64_t bucket_id = coordinator->fetch_add(batchSize);
        if (bucket_id > bc_bucket_num-1)
            break;
        for (uint32_t k = 0; k < batchSize; k++)
        {
            if (bucket_id > bc_bucket_num-1)
                break;

            while (latest_bucket_finished_load < (int32_t)(bucket_id+1))
            {
                // Load next row
                uint32_t row_id = latest_row_loaded + 1;

                if (row_id >= num_rows)
                {
                    latest_bucket_finished_load = bucket_id+1;
                    break;
                }

                uint64_t first_value_needed = bucket_id*kBC;
                uint32_t first_row_needed = first_value_needed/old_yc::row_divisor;
                if (row_id < first_row_needed)
                {
                    row_id = first_row_needed;
                }

                uint64_t first_value_contained = row_id*old_yc::row_divisor;
                uint64_t first_bucket_needed = first_value_contained/kBC;

                uint64_t last_value_contained = (row_id+1)*old_yc::row_divisor - 1;
                uint64_t last_bucket_needed = last_value_contained/kBC;

                uint64_t clear_start = latest_bucket_started_load+1;
                if (clear_start < first_bucket_needed)
                {
                    clear_start = first_bucket_needed;
                }

                // Clear temp buckets we will need
                for (uint32_t i = clear_start; i <= last_bucket_needed; i++)
                {
                    temp_bucket_fill_states[i%temp_buckets_needed] = 0;
                    entry_positions[i%temp_buckets_needed].clear();
                }
                latest_bucket_started_load = last_bucket_needed;

                // Load in the row
                for (auto& [numa_node, penguin] : *prev_penguins)
                {
                    for (uint32_t entry_id = 0; entry_id < penguin->getCountInRow(row_id); entry_id++)
                    {
                        auto entry = penguin->readEntry(row_id, entry_id);
                        uint32_t bucket_of_entry = entry.getY()/kBC;
                        uint32_t bucket_entry_id = temp_bucket_fill_states[bucket_of_entry%temp_buckets_needed]++;
                        temp_yc temp_entry;
                        temp_entry.setY(entry.getY());
                        temp_entry.c = entry.c;
                        temp_buckets[bucket_of_entry%temp_buckets_needed].set(bucket_entry_id, temp_entry);
                        auto test_entry = temp_buckets[bucket_of_entry%temp_buckets_needed].get(bucket_entry_id);
                        test_entry.row = bucket_of_entry;
                        if (table_index > 0)
                        {
                            uint64_t uid = penguin->getUniqueIdentifier(row_id, entry_id);
                            uint64_t new_position = (*new_entry_positions)[numa_node]->get(uid).getY();
                            entry_positions[bucket_of_entry%temp_buckets_needed].push_back(new_position);
                        }
                    }
                }

                latest_bucket_finished_load = (last_value_contained+1)/kBC - 1;
                latest_row_loaded = row_id;

                if (GetMaxY(table_index-1) <= last_value_contained)
                {
                    latest_bucket_finished_load = bc_bucket_num;
                }
            }

            bucket_left_done->set(bucket_id, BooleanPackedEntry(true));
            bucket_right_done->set(bucket_id+1, BooleanPackedEntry(true));

            // Check all rows that contribute to either bucket_id and pop them if possible
            uint64_t first_row_possible = (bucket_id*kBC)/old_yc::row_divisor;
            uint64_t last_row_possible = ((bucket_id+2)*kBC - 1)/old_yc::row_divisor;
            for (uint32_t row_id = first_row_possible; row_id <= last_row_possible; row_id++)
            {
                // Check if all buckets in this row are finished
                uint64_t first_bucket_contained = (row_id*old_yc::row_divisor)/kBC;
                uint64_t last_bucket_contained = ((row_id+1)*old_yc::row_divisor - 1)/kBC;
                uint32_t i;
                for (i = first_bucket_contained; i <= last_bucket_contained; i++)
                {
                    if (!bucket_left_done->get(i).getY() || !bucket_right_done->get(i).getY())
                    {
                        break;
                    }
                }
                if (i > last_bucket_contained)
                {
                    // All checks passed, this row can be popped
                    for (auto& [numa_node, penguin] : *prev_penguins) {
                        penguin->popRow(row_id);
                    }
                }
            }


            // Process sort_row and sort_row + 1
            uint16_t parity = bucket_id % 2;
            for (size_t yl : right_map_clean) {
                right_map[yl].clear();
            }
            right_map_clean.clear();
            for (size_t pos_R = 0; pos_R < temp_bucket_fill_states[(bucket_id+1)%temp_buckets_needed]; pos_R++) {
                auto right_entry = temp_buckets[(bucket_id+1)%temp_buckets_needed].get(pos_R);
                right_entry.row = bucket_id;
                uint64_t r_y = right_entry.getY() % kBC;
                if (right_map[r_y].empty())
                {
                    right_map_clean.push_back(r_y);
                }
                right_map[r_y].push_back(pos_R);
            }

            for (size_t pos_L = 0; pos_L < temp_bucket_fill_states[(bucket_id)%temp_buckets_needed]; pos_L++) {
                auto left_entry = temp_buckets[(bucket_id)%temp_buckets_needed].get(pos_L);
                left_entry.row = bucket_id;
                for (uint8_t i = 0; i < kExtraBitsPow; i++) {
                    uint16_t r_target = L_targets[parity][left_entry.getY()%kBC][i];
                    for (auto& pos_R : right_map[r_target]) {
                        auto right_entry = temp_buckets[(bucket_id+1)%temp_buckets_needed].get(pos_R);
                        right_entry.row = bucket_id+1;

                        // Calculate the next F and dump a new entry into the bucket index for the next table
                        uint64_t c_len = kVectorLens[table_index + 2] * K;
                        auto out = fx.CalculateBucket(
                                Bits(left_entry.getY(), K + kExtraBits),
                                Bits(left_entry.c, c_len),
                                Bits(right_entry.c, c_len));

                        new_yc new_entry;

                        if (table_index < 5) {
                            new_entry.setY(out.first.GetValue());
                            uint8_t buff[32];
                            memset(buff, 0, sizeof(buff));
                            out.second.ToBytes(buff);
                            new_entry.c = Util::SliceInt128FromBytes(
                                    buff, 0, kVectorLens[table_index + 3] * K);
                        }
                        else
                        {
                            new_entry.setY(out.first.Slice(0, K).GetValue());
                        }
                        uint64_t new_entry_id = new_yc_penguin->addEntry(new_entry);
                        uint64_t sort_row = new_entry.row;

                        // Emit a new line point for sorting into the phase1 output buffer
                        uint64_t x, y;
                        if (table_index == 0) // First table just gets the two x values
                        {
                            x = left_entry.c;
                            y = right_entry.c;
                        }
                        else
                        {
                            x = entry_positions[bucket_id%temp_buckets_needed][pos_L];
                            y = entry_positions[(bucket_id+1)%temp_buckets_needed][pos_R];
                        }
                        LinePointUIDPackedEntry<table_index> line_point_entry;
                        line_point_entry.setY(Encoding::SquareToLinePoint(x, y));
                        line_point_entry.uid = new_yc_penguin->getUniqueIdentifier(sort_row,
                                                                                         new_entry_id);
                        assert(line_point_entry.uid < (1ULL << (K + 2)));
                        uint32_t test_row_id = line_point_entry.row;
                        uint32_t test_entry_id = new_line_point_penguin->addEntry(line_point_entry);
                    }
                }

            }

            bucket_id++;
        }
    }
}

template <uint8_t K, uint32_t num_rows>
template <int8_t table_index>
void Plotter<K, num_rows>::phase1ThreadC(
        uint32_t cpu_id,
        const uint8_t* id,
        std::atomic<uint64_t>* coordinator,
        std::map<uint32_t, phase1_new_positions_type*>* new_entry_positions,
        map<uint32_t, Penguin<LinePointUIDPackedEntry<table_index>> *> line_point_bucket_indexes,
        vector<Park*> *parks,
        Buffer* buffer)
{
    PinToCpuid(cpu_id);
    F1Calculator f1(K, id);
    vector<pair<uint32_t, LinePointUIDPackedEntry<table_index>>> entries;
    entries.reserve(LinePointUIDPackedEntry<table_index>::max_entries_per_row);
    vector<uint128_t> line_points;
    line_points.reserve(LinePointUIDPackedEntry<table_index>::max_entries_per_row);
    vector<uint128_t> line_points_check(LinePointUIDPackedEntry<table_index>::max_entries_per_row);
    uint64_t last_offset_computed = 0;
    uint64_t total_entries_so_far = 0;
    while (true)
    {
        uint64_t row_id = coordinator->fetch_add(1);
        if (row_id >= num_rows)
            break;

        while (row_id > last_offset_computed)
        {
            for (auto& [numa_node, penguin] : line_point_bucket_indexes) {
                total_entries_so_far += penguin->getCountInRow(last_offset_computed);
            }
            last_offset_computed++;
        }


// Sort all entries in this bucket
        entries.clear();
        uint32_t num_entries = 0;
        for (auto& [numa_node, penguin] : line_point_bucket_indexes)
        {
            uint32_t entries_in_numa = penguin->getCountInRow(row_id);
            for (uint32_t i = 0; i < entries_in_numa; i++)
            {
                auto entry = penguin->readEntry(row_id, i);
                entries.push_back(pair<uint32_t, LinePointUIDPackedEntry<table_index>>(numa_node,entry));
            }
            num_entries += entries_in_numa;
            penguin->popRow(row_id);
        }

// Sort the bucket
        struct {
            bool operator()(pair<uint32_t, LinePointUIDPackedEntry<table_index>> a, pair<uint32_t, LinePointUIDPackedEntry<table_index>> b) const { return a.second.getY() < b.second.getY(); }
        } customLess;
        sort(entries.begin(), entries.begin()+num_entries, customLess);
// We have a sorted entries list! Emit to output buffer and record new positions
        auto new_park = new DeltaPark<line_point_delta_len_bits>(num_entries);
        new_park->start_pos = total_entries_so_far;
        uint64_t num_bytes = new_park->GetSpaceNeeded();
        new_park->bind(buffer->data + buffer->GetInsertionOffset(num_bytes));
        line_points.clear();
        for (uint64_t i = 0; i < num_entries; i++)
        {
            line_points.push_back(entries[i].second.getY());
            assert(entries[i].second.uid < (1ULL<<(K+2)));
            (*new_entry_positions)[entries[i].first]->set(entries[i].second.uid, PackedEntry<1, (1ULL<<(K+1)), 1>(i + total_entries_so_far));
        }
        new_park->addEntries(line_points);
        (*parks)[row_id] = new_park;

        // Check that we can recover the data
        new_park->readEntries(line_points_check);
        for (uint32_t i = 0; i < num_entries; i++)
        {
            assert(line_points[i] == line_points_check[i]);
        }

    }
}

template <uint8_t K, uint32_t num_rows>
void Plotter<K, num_rows>::phase1ThreadD(
        uint32_t cpu_id,
        std::atomic<uint64_t>* coordinator,
        std::map<uint32_t, phase1_new_positions_type*>* new_entry_positions,
        std::map<uint32_t, Penguin<YCPackedEntry<5>>*> line_point_penguins,
        vector<DeltaPark<finaltable_y_delta_len_bits> *> *parks,
        Buffer* buffer)
{
    PinToCpuid(cpu_id);
    vector<YCPackedEntry<5>> entries(YCPackedEntry<5>::max_entries_per_row);
    vector<uint128_t> line_points(YCPackedEntry<5>::max_entries_per_row);
    uint64_t total_entries_so_far = 0;
    uint64_t last_offset_computed = 0;
    while (true)
    {
        uint64_t row_id = coordinator->fetch_add(1);
        if (row_id >= num_rows)
            break;


        while (row_id > last_offset_computed)
        {
            for (auto& [numa_node, penguin] : line_point_penguins) {
                total_entries_so_far += penguin->getCountInRow(last_offset_computed);
            }
            last_offset_computed++;
        }

// Sort all entries in this bucket
        entries.clear();
        uint64_t num_entries = 0;
        for (auto& [numa_node, penguin] : line_point_penguins)
        {
            uint64_t entries_in_numa = penguin->getCountInRow(row_id);
            for (uint64_t i = 0; i < entries_in_numa; i++)
            {
                uint64_t new_pos = (*new_entry_positions)[numa_node]->get(penguin->getUniqueIdentifier(row_id, i)).getY();
                auto e = penguin->readEntry(row_id, i);
                e.setYtemp(e.getY()<<(K+1) | new_pos);
                entries.push_back(e);
            }
            num_entries += entries_in_numa;
            penguin->popRow(row_id);
        }
        //assert(num_entries < YCPackedEntry<K, 5>::max_entries_per_sort_row);

// Sort the bucket
        struct {
            bool operator()(YCPackedEntry<5> a, YCPackedEntry<5> b) const { return a.getY() < b.getY(); }
        } customLess;
        sort(entries.begin(), entries.begin()+num_entries, customLess);

// We have a sorted entries list! Emit to output buffer
        auto new_park = new DeltaPark<finaltable_y_delta_len_bits>(num_entries);
        new_park->start_pos = total_entries_so_far;
        uint64_t num_bytes = new_park->GetSpaceNeeded();
        new_park->bind(buffer->data + buffer->GetInsertionOffset(num_bytes));
        for (uint32_t i = 0; i < num_entries; i++)
        {
            line_points[i] = entries[i].getY();
        }
        new_park->addEntries(line_points);
        (*parks)[row_id] = new_park;
    }
}

template <uint8_t K, uint32_t num_rows, int8_t table_index>
using test_type = typename Plotter<K, num_rows>::template YCPackedEntry<table_index>;

template <uint8_t K, uint32_t num_rows>
template <int8_t table_index>
map<uint32_t, Penguin<test_type<K, num_rows, table_index>>*> Plotter<K, num_rows>::phase1DoTable(
        map<uint32_t, Penguin<YCPackedEntry<table_index-1>>*> prev_bucket_indexes)
{

    map<uint32_t, Penguin<YCPackedEntry<table_index>>*> next_bucket_indexes;
    map<uint32_t, Penguin<LinePointUIDPackedEntry<table_index>>*> new_line_point_bucket_indexes;
    for (auto& numa_node : GetNUMANodesFromCpuIds(cpu_ids))
    {
        next_bucket_indexes[numa_node] = new Penguin<YCPackedEntry<table_index>>();
        new_line_point_bucket_indexes[numa_node] = new Penguin<LinePointUIDPackedEntry<table_index>>();
    }
    auto * bucket_left_done = new p1_buckets_done_type<table_index-1>();
    bucket_left_done->fill(BooleanPackedEntry(false));
    auto * bucket_right_done = new p1_buckets_done_type<table_index-1>();
    bucket_right_done->fill(BooleanPackedEntry(false));
    bucket_left_done->set(0, BooleanPackedEntry(true));
    bucket_right_done->set(GetMaxY(table_index)/kBC - 1, BooleanPackedEntry(true));

    cout << "Part B"<< (uint32_t)table_index;
    uint64_t start_seconds = time(nullptr);
    vector<thread> threads;
    std::atomic<uint64_t> coordinator = 0;
    for (auto & cpu_id : cpu_ids)
    {
        uint32_t numa_node = numa_node_of_cpu(cpu_id);
        threads.push_back(thread(
                Plotter<K, num_rows>::phase1ThreadB<table_index>,
                cpu_id,
                id,
                &coordinator,
                bucket_left_done,
                bucket_right_done,
                &phase1_new_entry_positions,
                &prev_bucket_indexes,
                next_bucket_indexes[numa_node],
                new_line_point_bucket_indexes[numa_node]));
    }

    for (auto &it: threads)
    {
        it.join();
    }

    delete bucket_left_done;
    delete bucket_right_done;

    for (auto& numa_node : GetNUMANodesFromCpuIds(cpu_ids)) {
        delete prev_bucket_indexes[numa_node];
    }
    cout << " (" << time(nullptr) - start_seconds << "s)" << endl;

// Unfortunately single-threaded, but not much can be done to help it
    //new_line_point_bucket_index->SumBucketOffsets();

// Setup the new park list and buffer
    DeltaPark<line_point_delta_len_bits> test_park(LinePointUIDPackedEntry<table_index>::max_entries_per_row);
    buffers.push_back(new Buffer(test_park.GetSpaceNeeded() * LinePointUIDPackedEntry<table_index>::num_rows_v));
    phase1_graph_parks.push_back(vector<Park*>(LinePointUIDPackedEntry<table_index>::num_rows_v));

    coordinator = 0;
    cout << "Part C"<< (uint32_t)table_index;
    start_seconds = time(nullptr);
    threads.clear();
    for (auto & cpu_id : cpu_ids)
    {
        threads.push_back(thread(
                Plotter<K, num_rows>::phase1ThreadC<table_index>,
                cpu_id,
                id,
                &coordinator,
                &phase1_new_entry_positions,
                new_line_point_bucket_indexes,
                &(phase1_graph_parks[table_index]),
                buffers[table_index]
        ));
    }

    for (auto &it: threads)
    {
        it.join();
    }
    cout << " (" << time(nullptr) - start_seconds << "s)" << endl;
    for (auto& numa_node : GetNUMANodesFromCpuIds(cpu_ids)) {
        delete new_line_point_bucket_indexes[numa_node];
    }
    return next_bucket_indexes;
}


template <uint8_t K, uint32_t num_rows>
void Plotter<K, num_rows>::phase1()
{
    load_tables();

    uint64_t phase_start_seconds = time(nullptr);

    map<uint32_t, Penguin<YCPackedEntry<-1>>*> first_bucket_indexes;
    for (auto& numa_node : GetNUMANodesFromCpuIds(cpu_ids))
    {
        first_bucket_indexes[numa_node] = new Penguin<YCPackedEntry<-1>>();
    }

    cout << "Part A";
    uint64_t part_start_seconds = time(nullptr);
    vector<thread> threads;
    std::atomic<uint64_t> coordinator = 0;
    for (auto & cpu_id : cpu_ids)
    {
        uint32_t numa_node = numa_node_of_cpu(cpu_id);
        threads.push_back(thread(
                Plotter<K, num_rows>::phase1ThreadA,
                cpu_id,
                &coordinator,
                id,
                first_bucket_indexes[numa_node]));
    }

    for (auto &it: threads)
    {
        it.join();
    }
    cout << " (" << time(nullptr) - part_start_seconds << "s)" << endl;


    for (auto& numa_node : GetNUMANodesFromCpuIds(cpu_ids))
    {
        phase1_new_entry_positions[numa_node] = (phase1_new_positions_type*) malloc(sizeof(phase1_new_positions_type));
    }

    auto i0 = phase1DoTable<0>(first_bucket_indexes);
    auto i1 = phase1DoTable<1>(i0);
    auto i2 = phase1DoTable<2>(i1);
    auto i3 = phase1DoTable<3>(i2);
    auto i4 = phase1DoTable<4>(i3);
    auto i5 = phase1DoTable<5>(i4);

    DeltaPark<finaltable_y_delta_len_bits> test_park(YCPackedEntry<5>::max_entries_per_row);
    buffers.push_back(new Buffer(test_park.GetSpaceNeeded() * YCPackedEntry<5>::num_rows_v));
    phase1_final_parks.resize(YCPackedEntry<5>::num_rows_v);

    uint64_t num_final_matches = 0;
    for (auto& [numa_node, penguin] : i5)
    {
        for (uint64_t row_id = 0; row_id < num_rows; row_id++)
        {
            num_final_matches += penguin->getCountInRow(row_id);
        }
    }
    cout << "Part D";
    part_start_seconds = time(nullptr);
    threads.clear();
    coordinator = 0;
    for (auto & cpu_id : cpu_ids)
    {
        threads.push_back(thread(
                Plotter<K, num_rows>::phase1ThreadD,
                cpu_id,
                &coordinator,
                &phase1_new_entry_positions,
                i5,
                &phase1_final_parks,
                buffers[6]));
    }
    for (auto &it: threads)
    {
        it.join();
    }
    cout << " (" << time(nullptr) - part_start_seconds << "s)" << endl;
    cout << "Phase 1 finished in " << time(nullptr) - phase_start_seconds << "s" << endl;


    for (auto& numa_node : GetNUMANodesFromCpuIds(cpu_ids))
    {
        free(phase1_new_entry_positions[numa_node]);
    }

    num_final_matches = (*(phase1_final_parks.end() - 1))->start_pos + (*(phase1_final_parks.end() - 1))->size();}

#include "explicit_templates.hpp"