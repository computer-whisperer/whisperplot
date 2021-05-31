#include <vector>
#include <atomic>
#include <cstring>
#include <ctime>

#include "penguin.hpp"
#include "buffer.hpp"
#include "calculate_bucket.hpp"
#include "encoding.hpp"
#include "park.hpp"
#include "phase1.hpp"
#include "pos_constants.hpp"
#include "thread_mgr.hpp"

using namespace std;

template <uint8_t K>
void Phase1<K>::ThreadA(
        uint32_t cpu_id,
        std::atomic<uint64_t>* coordinator,
        const uint8_t* id,
        Penguin<YCPackedEntry<K, -1>>* new_penguin
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
            YCPackedEntry<K, -1> entry;
            entry.setY(buff[i]);
            entry.c = x+i;
            new_penguin->addEntry(entry);
        }
    }
}

template <uint8_t K>
template <int8_t table_index>
void Phase1<K>::ThreadB(
        uint32_t cpu_id,
        std::atomic<uint64_t> * coordinator,
        map<uint32_t, vector<uint32_t>> * new_entry_positions,
        map<uint32_t, Penguin<YCPackedEntry<K, table_index - 1>>*> * prev_penguins,
        Penguin<YCPackedEntry<K, table_index>> * new_yc_penguin,
        Penguin<LinePointEntryUIDPackedEntry<K>> * new_line_point_penguin)
{
    PinToCpuid(cpu_id);
    vector<vector<uint32_t>> right_map(kBC);
    for (auto& i : right_map) {
        i.clear();
    }
    std::vector<uint32_t> right_map_clean;
    FxCalculator fx(K, table_index+2);

    // These are BC buckets from the chia algo
    constexpr uint32_t bc_buckets_per_sort_row = 32;
    constexpr uint32_t bc_bucket_num = (1ULL << (K+kExtraBits))/kBC;
    vector<PackedArray<YCPackedEntry<K, table_index - 1>, 350>> temp_buckets(bc_buckets_per_sort_row*2);
    uint32_t latest_bc_bucket_loaded = 0;
    vector<vector<uint32_t>> entry_positions(bc_buckets_per_sort_row*2);
    for (auto & i : entry_positions)
    {
        i.reserve(350);
    }

    constexpr uint32_t batchSize = bc_buckets_per_sort_row*32;

    while (true)
    {
        uint64_t bucket_id = coordinator->fetch_add(batchSize);
        if (bucket_id > bc_bucket_num-1)
            break;
        for (uint32_t k = 0; k < batchSize; k++)
        {
            if (bucket_id > bc_bucket_num-1)
                break;
            uint16_t parity = bucket_id % 2;

            if (latest_bc_bucket_loaded < bucket_id+1) {
                // Load more buckets, starting at bc bucket sort_row+1
                uint32_t row_id = (bucket_id + 1) / bc_buckets_per_sort_row;
                // Clear the entry position lists for the incoming row
                if (row_id % 2 == 0)
                {
                    for (uint32_t i = 0; i < bc_buckets_per_sort_row; i++)
                    {
                        entry_positions[i].clear();
                        temp_buckets[i].clear();
                    }
                }
                else
                {
                    for (uint32_t i = bc_buckets_per_sort_row; i < bc_buckets_per_sort_row*2; i++)
                    {
                        entry_positions[i].clear();
                        temp_buckets[i].clear();
                    }
                }
                // Load in next set of BC buckets
                for (auto& [numa_node, penguin] : *prev_penguins)
                {
                    for (uint32_t entry_id = 0; entry_id < penguin->getCountInRow(row_id); entry_id++)
                    {
                        auto entry = penguin->readEntry(row_id, entry_id);
                        uint32_t bucket_of_entry = entry.getY()/kBC;
                        temp_buckets[bucket_of_entry%(bc_buckets_per_sort_row*2)].append(entry);
                        if (table_index > 0)
                        {
                            uint64_t uid = penguin->getUniqueIdentifier(row_id, entry_id);
                            uint32_t new_position = (*new_entry_positions)[numa_node][uid];
                            entry_positions[bucket_of_entry%(bc_buckets_per_sort_row*2)].push_back(new_position);
                        }
                    }

                }
                latest_bc_bucket_loaded += bc_buckets_per_sort_row;
            }
            for (auto& [numa_node, penguin] : *prev_penguins) {
                penguin->popRow(bucket_id/ bc_buckets_per_sort_row);
                penguin->popRow((bucket_id + 1) / bc_buckets_per_sort_row);
            }

            // Process sort_row and sort_row + 1
            for (size_t yl : right_map_clean) {
                right_map[yl].clear();
            }
            right_map_clean.clear();
            for (size_t pos_R = 0; pos_R < temp_buckets[(bucket_id+1)%(bc_buckets_per_sort_row*2)].count; pos_R++) {
                auto right_entry = temp_buckets[(bucket_id+1)%(bc_buckets_per_sort_row*2)].read(pos_R);
                right_entry.sort_row = bucket_id / bc_buckets_per_sort_row;
                uint64_t r_y = right_entry.getY() % kBC;
                if (right_map[r_y].empty())
                {
                    right_map_clean.push_back(r_y);
                }
                right_map[r_y].push_back(pos_R);
            }

            for (size_t pos_L = 0; pos_L < temp_buckets[(bucket_id)%(bc_buckets_per_sort_row*2)].count; pos_L++) {
                auto left_entry = temp_buckets[(bucket_id)%(bc_buckets_per_sort_row*2)].read(pos_L);
                left_entry.sort_row = bucket_id / bc_buckets_per_sort_row;
                for (uint8_t i = 0; i < kExtraBitsPow; i++) {
                    uint16_t r_target = L_targets[parity][left_entry.getY()%kBC][i];
                    for (auto& pos_R : right_map[r_target]) {
                        auto right_entry = temp_buckets[(bucket_id+1)%(bc_buckets_per_sort_row*2)].read(pos_R);

                        // Calculate the next F and dump a new entry into the bucket index for the next table
                        uint64_t c_len = kVectorLens[table_index + 2] * K;
                        auto out = fx.CalculateBucket(
                                Bits(left_entry.getY(), K + kExtraBits),
                                Bits(left_entry.c, c_len),
                                Bits(right_entry.c, c_len));

                        YCPackedEntry<K, table_index> new_entry;
                        new_entry.setY(out.first.GetValue());
                        if (table_index < 5) {
                            uint8_t buff[32];
                            memset(buff, 0, sizeof(buff));
                            out.second.ToBytes(buff);
                            new_entry.c = Util::SliceInt128FromBytes(
                                    buff, 0, kVectorLens[table_index + 3] * K);
                        }
                        else
                        {
                            new_entry.setY(out.first.GetValue() % (1ULL << K));
                        }
                        uint64_t new_entry_id = new_yc_penguin->addEntry(new_entry);
                        uint64_t sort_row = new_entry.sort_row;

                        // Emit a new line point for sorting into the phase1 output buffer
                        uint64_t x, y;
                        if (table_index == 0) // First table just gets the two x values
                        {
                            x = left_entry.c;
                            y = right_entry.c;
                        }
                        else
                        {
                            x = entry_positions[bucket_id%(bc_buckets_per_sort_row*2)][pos_L];
                            y = entry_positions[(bucket_id+1)%(bc_buckets_per_sort_row*2)][pos_R];
                        }
                        LinePointEntryUIDPackedEntry<K> line_point_entry;
                        line_point_entry.setLinePoint(Encoding::SquareToLinePoint(x, y));
                        line_point_entry.entry_uid = new_yc_penguin->getUniqueIdentifier(sort_row,
                                                                                         new_entry_id);
                        assert(line_point_entry.entry_uid < (1ULL << (K + 2)));
                        new_line_point_penguin->addEntry(line_point_entry);
                    }
                }

            }

            bucket_id++;
        }
    }
}

template <uint8_t K>
template <int8_t table_index>
void Phase1<K>::ThreadC(
        uint32_t cpu_id,
        std::atomic<uint64_t>* coordinator,
        map<uint32_t, vector<uint32_t>> * new_entry_positions,
        map<uint32_t, Penguin<LinePointEntryUIDPackedEntry<K>>*> line_point_bucket_indexes,
        std::vector<TemporaryPark<line_point_delta_len_bits>*>* parks,
        Buffer* buffer)
{
    PinToCpuid(cpu_id);
    vector<pair<uint32_t, LinePointEntryUIDPackedEntry<K>>> entries;
    entries.reserve(LinePointEntryUIDPackedEntry<K>::max_entries_per_sort_row);
    vector<uint128_t> line_points;
    line_points.reserve(LinePointEntryUIDPackedEntry<K>::max_entries_per_sort_row);
    uint64_t last_offset_computed = 0;
    uint64_t total_entries_so_far = 0;
    while (true)
    {
        uint64_t row_id = coordinator->fetch_add(1);
        if (row_id >= LinePointEntryUIDPackedEntry<K>::num_sort_rows)
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
                entries.push_back(pair<uint32_t, LinePointEntryUIDPackedEntry<K>>(numa_node,
                                                                                  penguin->readEntry(
                                                                                                   row_id, i)));

            }
            num_entries += entries_in_numa;
           // penguin->popRow(sort_row);
        }
        assert(num_entries < LinePointEntryUIDPackedEntry<K>::max_entries_per_sort_row);

// Sort the bucket
        struct {
            bool operator()(pair<uint32_t, LinePointEntryUIDPackedEntry<K>> a, pair<uint32_t, LinePointEntryUIDPackedEntry<K>> b) const { return a.second.getLinePoint() < b.second.getLinePoint(); }
        } customLess;
        sort(entries.begin(), entries.begin()+num_entries, customLess);
// We have a sorted entries list! Emit to output buffer and record new positions
        auto new_park = new TemporaryPark<line_point_delta_len_bits>(num_entries);
        uint64_t num_bytes = new_park->GetSpaceNeeded();
        new_park->Bind(buffer->data + buffer->GetInsertionOffset(num_bytes));
        line_points.clear();
        for (uint32_t i = 0; i < num_entries; i++)
        {
            line_points.push_back(entries[i].second.getLinePoint());
            assert(entries[i].second.entry_uid < (1ULL<<(K+2)));
            (*new_entry_positions)[entries[i].first][entries[i].second.entry_uid] = i + total_entries_so_far;
        }
        new_park->AddEntries(line_points.data());
        (*parks)[row_id] = new_park;
    }
}

template <uint8_t K>
void Phase1<K>::ThreadD(
        uint32_t cpu_id,
        std::atomic<uint64_t>* coordinator,
        map<uint32_t, vector<uint32_t>> * new_entry_positions,
        std::map<uint32_t, Penguin<YCPackedEntry<K, 5>>*> line_point_bucket_indexes,
        std::vector<TemporaryPark<finaltable_y_delta_len_bits>*>* parks,
        Buffer* buffer)
{
    PinToCpuid(cpu_id);
    vector<YCPackedEntry<K, 5>> entries(YCPackedEntry<K, 5>::max_entries_per_sort_row);
    vector<uint128_t> line_points(YCPackedEntry<K, 5>::max_entries_per_sort_row);
    uint64_t total_entries_so_far = 0;
    while (true)
    {
        uint64_t bucket_id = coordinator->fetch_add(1);
        if (bucket_id >= YCPackedEntry<K, 5>::num_sort_rows)
            break;

        for (auto& [numa_node, penguin] : line_point_bucket_indexes) {
            total_entries_so_far += penguin->getCountInRow(bucket_id);
        }
// Sort all entries in this bucket
        entries.clear();
        uint32_t num_entries = 0;
        for (auto& [numa_node, bucket_index] : line_point_bucket_indexes)
        {
            uint32_t entries_in_numa = bucket_index->getCountInRow(bucket_id);
            for (uint32_t i = 0; i < entries_in_numa; i++)
            {
                uint32_t new_pos = (*new_entry_positions)[numa_node][bucket_index->getUniqueIdentifier(bucket_id, i)];
                auto e = bucket_index->readEntry(bucket_id, i);
                e.setY(e.getY()<<K | new_pos);
                entries[num_entries + i] = e;
            }
            num_entries += entries_in_numa;
            bucket_index->popRow(bucket_id);
        }
        //assert(num_entries < YCPackedEntry<K, 5>::max_entries_per_sort_row);

// Sort the bucket
        struct {
            bool operator()(YCPackedEntry<K, 5> a, YCPackedEntry<K, 5> b) const { return a.getY() < b.getY(); }
        } customLess;
        sort(entries.begin(), entries.begin()+num_entries, customLess);

// We have a sorted entries list! Emit to output buffer and record new positions
        auto new_park = new TemporaryPark<finaltable_y_delta_len_bits>(num_entries);
        uint64_t num_bytes = new_park->GetSpaceNeeded();
        new_park->Bind(buffer->data + buffer->GetInsertionOffset(num_bytes));
        for (uint32_t i = 0; i < num_entries; i++)
        {
            line_points[i] = entries[i].getY();
        }
        new_park->AddEntries(line_points.data());
        (*parks)[bucket_id] = new_park;
    }
}

template <uint8_t K>
template <int8_t table_index>
map<uint32_t, Penguin<YCPackedEntry<K, table_index>>*> Phase1<K>::DoTable(
        map<uint32_t, Penguin<YCPackedEntry<K, table_index - 1>>*> prev_bucket_indexes)
{
    map<uint32_t, Penguin<YCPackedEntry<K, table_index>>*> next_bucket_indexes;
    map<uint32_t, Penguin<LinePointEntryUIDPackedEntry<K>>*> new_line_point_bucket_indexes;
    for (auto& numa_node : GetNUMANodesFromCpuIds(cpu_ids))
    {
        d_new_entry_positions[numa_node].resize((1ULL<<(K+2)));
        next_bucket_indexes[numa_node] = new Penguin<YCPackedEntry<K, table_index>>();
        next_bucket_indexes[numa_node]->pops_required_per_row = 2*32;
        new_line_point_bucket_indexes[numa_node] = new Penguin<LinePointEntryUIDPackedEntry<K>>();
    }

    cout << "Part B"<< (uint32_t)table_index;
    uint64_t start_seconds = time(nullptr);
    vector<thread> threads;
    std::atomic<uint64_t> coordinator = 0;
    for (uint32_t i = 0; i < cpu_ids.size(); i++)
    {
        uint32_t numa_node = numa_node_of_cpu(cpu_ids[i]);
        threads.push_back(thread(
                Phase1<K>::ThreadB<table_index>,
                cpu_ids[i],
                &coordinator,
                &d_new_entry_positions,
                &prev_bucket_indexes,
                next_bucket_indexes[numa_node],
                new_line_point_bucket_indexes[numa_node]));
    }

    for (auto &it: threads)
    {
        it.join();
    }

    for (auto& numa_node : GetNUMANodesFromCpuIds(cpu_ids)) {
        delete prev_bucket_indexes[numa_node];
    }
    cout << " (" << time(nullptr) - start_seconds << "s)" << endl;

// Unfortunately single-threaded, but not much can be done to help it
    //new_line_point_bucket_index->SumBucketOffsets();

// Setup the new park list and buffer
    TemporaryPark<line_point_delta_len_bits> test_park(LinePointEntryUIDPackedEntry<K>::max_entries_per_sort_row);
    buffers.push_back(new Buffer(test_park.GetSpaceNeeded() * LinePointEntryUIDPackedEntry<K>::num_sort_rows));
    graph_parks.push_back(vector<TemporaryPark<line_point_delta_len_bits>*>(LinePointEntryUIDPackedEntry<K>::num_sort_rows));

    coordinator = 0;
    cout << "Part C"<< (uint32_t)table_index;
    start_seconds = time(nullptr);
    threads.clear();
    for (uint32_t i = 0; i < cpu_ids.size(); i++)
    {
        threads.push_back(thread(
                Phase1<K>::ThreadC<table_index>,
                cpu_ids[i],
                &coordinator,
                &d_new_entry_positions,
                new_line_point_bucket_indexes,
                &(graph_parks[table_index]),
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


template <uint8_t K>
Phase1<K>::Phase1(const uint8_t* id, std::vector<uint32_t> cpu_ids_in)
{

    load_tables();

    cpu_ids = cpu_ids_in;
    uint64_t phase_start_seconds = time(nullptr);

    map<uint32_t, Penguin<YCPackedEntry<K, -1>>*> first_bucket_indexes;
    for (auto& numa_node : GetNUMANodesFromCpuIds(cpu_ids))
    {
        first_bucket_indexes[numa_node] = new Penguin<YCPackedEntry<K, -1>>();
        first_bucket_indexes[numa_node]->pops_required_per_row = 2*32;
    }

    cout << "Part A";
    uint64_t part_start_seconds = time(nullptr);
    vector<thread> threads;
    std::atomic<uint64_t> coordinator = 0;
    for (uint32_t i = 0; i < cpu_ids.size(); i++)
    {
        uint32_t numa_node = numa_node_of_cpu(cpu_ids[i]);
        threads.push_back(thread(
                Phase1<K>::ThreadA,
                cpu_ids[i],
                &coordinator,
                id,
                first_bucket_indexes[numa_node]));
    }

    for (auto &it: threads)
    {
        it.join();
    }
    cout << " (" << time(nullptr) - part_start_seconds << "s)" << endl;

    auto i0 = DoTable<0>(first_bucket_indexes);
    auto i1 = DoTable<1>(i0);
    auto i2 = DoTable<2>(i1);
    auto i3 = DoTable<3>(i2);
    auto i4 = DoTable<4>(i3);
    auto i5 = DoTable<5>(i4);

    TemporaryPark<finaltable_y_delta_len_bits> test_park(YCPackedEntry<K, 5>::max_entries_per_sort_row);
    buffers.push_back(new Buffer(test_park.GetSpaceNeeded() * YCPackedEntry<K, 5>::num_sort_rows));
    final_parks.resize(YCPackedEntry<K, 5>::num_sort_rows);

    cout << "Part D";
    part_start_seconds = time(nullptr);
    threads.clear();
    coordinator = 0;
    for (uint32_t i = 0; i < cpu_ids.size(); i++)
    {
        threads.push_back(thread(
                Phase1<K>::ThreadD,
                cpu_ids[i],
                &coordinator,
                &d_new_entry_positions,
                i5,
                &final_parks,
                buffers[6]));
    }
    for (auto &it: threads)
    {
        it.join();
    }
    cout << " (" << time(nullptr) - part_start_seconds << "s)" << endl;
    cout << "Phase1 finished in " << time(nullptr) - phase_start_seconds << "s" << endl;
}

template class Phase1<18>;
template class Phase1<22>;
template class Phase1<26>;
template class Phase1<32>;
template class Phase1<33>;