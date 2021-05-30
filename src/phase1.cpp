#include <vector>
#include <atomic>
#include <cstring>
#include <ctime>

#include "bucket_page_index.hpp"
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
        BucketPageIndex<YCBucketEntry<K, -1>>* new_bucket_index
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
            YCBucketEntry<K, -1> entry;
            entry.y = buff[i];
            entry.c = x+i;
            new_bucket_index->InsertEntry(entry);
        }
    }
}

template <uint8_t K>
template <int8_t table_index>
void Phase1<K>::ThreadB(
        uint32_t cpu_id,
        std::atomic<uint64_t> * coordinator,
        map<uint32_t, vector<uint32_t>> * new_entry_positions,
        map<uint32_t, BucketPageIndex<YCBucketEntry<K, table_index-1>>*> * prev_bucket_indexes,
        BucketPageIndex<YCBucketEntry<K, table_index>> * new_bucket_index,
        BucketPageIndex<LinePointEntryUIDBucketEntry<K>> * new_line_point_bucket_index)
{
    PinToCpuid(cpu_id);
    vector<vector<pair<uint32_t, uint32_t>>> right_map(kBC);
    std::vector<uint32_t> right_map_clean;
    FxCalculator fx(K, table_index+2);
    while (true)
    {
        uint64_t bucket_id = coordinator->fetch_add(1);
        if (bucket_id >= YCBucketEntry<K, table_index>::num_buckets)
            break;

        uint16_t parity = bucket_id % 2;
        if (bucket_id < YCBucketEntry<K, table_index>::num_buckets-1)
        {
            for (size_t yl : right_map_clean) {
                right_map[yl].clear();
            }
            right_map_clean.clear();
            for (auto& it : *prev_bucket_indexes)
            {
                for (size_t pos_R = 0; pos_R < it.second->GetCountInBucket(bucket_id + 1); pos_R++) {
                    uint64_t r_y = it.second->ReadEntry(bucket_id + 1, pos_R).y % kBC;
                    if (right_map[r_y].empty())
                    {
                        right_map_clean.push_back(r_y);
                    }
                    right_map[r_y].push_back(pair<uint32_t, uint32_t>(it.first, pos_R));
                }
            }

        }
        for (auto& it : *prev_bucket_indexes)
        {
            for (size_t pos_L = 0; pos_L < it.second->GetCountInBucket(bucket_id); pos_L++) {
                auto left_entry = it.second->ReadEntry(bucket_id, pos_L);
                if (bucket_id < YCBucketEntry<K, table_index>::num_buckets-1) {
                    for (uint8_t i = 0; i < kExtraBitsPow; i++) {
                        uint16_t r_target = L_targets[parity][left_entry.y%kBC][i];
                        for (auto& r_match : right_map[r_target]) {
                            auto right_entry = (*prev_bucket_indexes)[r_match.first]->ReadEntry(bucket_id + 1, r_match.second);

                            // Calculate the next F and dump a new entry into the bucket index for the next table
                            uint64_t c_len = kVectorLens[table_index + 2] * K;
                            auto out = fx.CalculateBucket(
                                    Bits(left_entry.y, K + kExtraBits),
                                    Bits(left_entry.c, c_len),
                                    Bits(right_entry.c, c_len));

                            YCBucketEntry<K, table_index> new_entry;
                            new_entry.y = out.first.GetValue();
                            if (table_index < 5) {
                                uint8_t buff[32];
                                memset(buff, 0, sizeof(buff));
                                out.second.ToBytes(buff);
                                new_entry.c = Util::SliceInt128FromBytes(
                                        buff, 0, kVectorLens[table_index + 3] * K);
                            }
                            else
                            {
                                new_entry.y %= (1ULL << K);
                            }
                            uint64_t new_entry_id = new_bucket_index->InsertEntry(new_entry);
                            uint64_t new_bucket_id = new_entry.get_bucket_id();

                            // Emit a new line point for sorting into the phase1 output buffer
                            uint64_t x, y;
                            if (table_index == 0) // First table just gets the two x values
                            {
                                x = left_entry.c;
                                y = right_entry.c;
                            }
                            else
                            {
                                // use the unique ids for the two values
                                x = new_entry_positions[it.first][(*prev_bucket_indexes)[it.first]->GetUniqueIdentifier(bucket_id, pos_L)];
                                y = new_entry_positions[r_match.first][(*prev_bucket_indexes)[r_match.first]->GetUniqueIdentifier(bucket_id + 1, r_match.second)];
                            }
                            LinePointEntryUIDBucketEntry<K> line_point_entry;
                            line_point_entry.line_point = Encoding::SquareToLinePoint(x, y);
                            line_point_entry.entry_uid = new_bucket_index->GetUniqueIdentifier(new_bucket_id, new_entry_id);
                            assert(line_point_entry.entry_uid < (1ULL << (K + 2)));
                            new_line_point_bucket_index->InsertEntry(line_point_entry);
                        }
                    }
                }
            }
        }
        for (auto& it : *prev_bucket_indexes)
        {
            it.second->PopBucket(bucket_id);
            it.second->PopBucket(bucket_id + 1);
        }
    }
}

template <uint8_t K>
template <int8_t table_index>
void Phase1<K>::ThreadC(
        uint32_t cpu_id,
        std::atomic<uint64_t>* coordinator,
        map<uint32_t, vector<uint32_t>> * new_entry_positions,
        map<uint32_t, BucketPageIndex<LinePointEntryUIDBucketEntry<K>>*> line_point_bucket_indexes,
        std::vector<TemporaryPark<line_point_delta_len_bits>*>* parks,
        Buffer* buffer)
{
    PinToCpuid(cpu_id);
    vector<pair<uint32_t, LinePointEntryUIDBucketEntry<K>>> entries(LinePointEntryUIDBucketEntry<K>::max_entries_per_bucket);
    uint64_t total_entries_so_far = 0;
    while (true)
    {
        uint64_t bucket_id = coordinator->fetch_add(1);
        if (bucket_id >= LinePointEntryUIDBucketEntry<K>::num_buckets)
            break;

        for (auto& [numa_node, bucket_index] : line_point_bucket_indexes) {
            total_entries_so_far += bucket_index->GetCountInBucket(bucket_id);
        }

// Sort all entries in this bucket
        entries.clear();
        uint32_t num_entries = 0;
        for (auto& [numa_node, bucket_index] : line_point_bucket_indexes)
        {
            uint32_t entries_in_numa = bucket_index->GetCountInBucket(bucket_id);
            for (uint32_t i = 0; i < entries_in_numa; i++)
            {
                entries[num_entries + i] = pair<uint32_t, LinePointEntryUIDBucketEntry<K>>(numa_node, bucket_index->ReadEntry(bucket_id, i));

            }
            num_entries += entries_in_numa;
            bucket_index->PopBucket(bucket_id);
        }
        assert(num_entries < LinePointEntryUIDBucketEntry<K>::max_entries_per_bucket);

// Sort the bucket
        struct {
            bool operator()(auto a, auto b) const { return a.second.line_point < b.second.line_point; }
        } customLess;
        sort(entries.begin(), entries.begin()+num_entries, customLess);
// We have a sorted entries list! Emit to output buffer and record new positions
        auto new_park = new TemporaryPark<line_point_delta_len_bits>(num_entries);
        uint64_t num_bytes = new_park->GetSpaceNeeded();
        new_park->Bind(buffer->data + buffer->GetInsertionOffset(num_bytes));
        vector<uint128_t> line_points(LinePointEntryUIDBucketEntry<K>::max_entries_per_bucket);
        for (uint32_t i = 0; i < num_entries; i++)
        {
            line_points[i] = entries[i].second.line_point;
            assert(entries[i].second.entry_uid < (1ULL<<(K+2)));
            new_entry_positions[entries.first][entries[i].second.entry_uid] = i + total_entries_so_far;
        }
        new_park->AddEntries(line_points.data());
        (*parks)[bucket_id] = new_park;
    }
}

template <uint8_t K>
void Phase1<K>::ThreadD(
        uint32_t cpu_id,
        std::atomic<uint64_t>* coordinator,
        map<uint32_t, vector<uint32_t>> * new_entry_positions,
        std::map<uint32_t, BucketPageIndex<YCBucketEntry<K, 5>>*> line_point_bucket_indexes,
        std::vector<TemporaryPark<finaltable_y_delta_len_bits>*>* parks,
        Buffer* buffer)
{
    PinToCpuid(cpu_id);
    vector<pair<uint32_t, YCBucketEntry<K, 5>>> entries(YCBucketEntry<K, 5>::max_entries_per_bucket);
    uint64_t total_entries_so_far = 0;
    while (true)
    {
        uint64_t bucket_id = coordinator->fetch_add(1);
        if (bucket_id >= YCBucketEntry<K, 5>::num_buckets)
            break;

        for (auto& [numa_node, bucket_index] : line_point_bucket_indexes) {
            total_entries_so_far += bucket_index->GetCountInBucket(bucket_id);
        }
// Sort all entries in this bucket
        entries.clear();
        uint32_t num_entries = 0;
        for (auto& [numa_node, bucket_index] : line_point_bucket_indexes)
        {
            uint32_t entries_in_numa = bucket_index->GetCountInBucket(bucket_id);
            for (uint32_t i = 0; i < entries_in_numa; i++)
            {
                entries[num_entries + i] = pair<uint32_t, YCBucketEntry<K, 5>>(numa_node, bucket_index->ReadEntry(bucket_id, i));

            }
            num_entries += entries_in_numa;
            bucket_index->PopBucket(bucket_id);
        }
        assert(num_entries < YCBucketEntry<K, 5>::max_entries_per_bucket);

// Sort the bucket
        struct {
            bool operator()(auto a, auto b) const { return a.second.y < b.second.y; }
        } customLess;
        sort(entries.begin(), entries.begin()+num_entries, customLess);

// We have a sorted entries list! Emit to output buffer and record new positions
        auto new_park = new TemporaryPark<finaltable_y_delta_len_bits>(num_entries);
        uint64_t num_bytes = new_park->GetSpaceNeeded();
        new_park->Bind(buffer->data + buffer->GetInsertionOffset(num_bytes));
        vector<uint128_t> line_points(YCBucketEntry<K, 5>::max_entries_per_bucket);
        for (uint32_t i = 0; i < num_entries; i++)
        {
            line_points[i] = entries[i].y;
        }
        new_park->AddEntries(line_points.data());
        (*parks)[bucket_id] = new_park;
    }
}

template <uint8_t K>
template <int8_t table_index>
map<uint32_t, BucketPageIndex<YCBucketEntry<K, table_index>>*> Phase1<K>::DoTable(
        map<uint32_t, BucketPageIndex<YCBucketEntry<K, table_index-1>>*> prev_bucket_indexes)
{
    map<uint32_t, BucketPageIndex<YCBucketEntry<K, table_index-1>>*> next_bucket_indexes;
    map<uint32_t, BucketPageIndex<LinePointEntryUIDBucketEntry<K>>*> new_line_point_bucket_indexes;
    for (auto& numa_node : GetNUMANodesFromCpuIds(cpu_ids))
    {
        next_bucket_indexes[numa_node] = new BucketPageIndex<YCBucketEntry<K, table_index>>();
        next_bucket_indexes[numa_node]->pops_required_per_bucket = 2;
        new_line_point_bucket_indexes[numa_node] = new BucketPageIndex<LinePointEntryUIDBucketEntry<K>>();
    }

    cout << "Part B"<< (uint32_t)table_index;
    uint64_t start_seconds = time(NULL);
    vector<thread> threads;
    std::atomic<uint64_t> coordinator = 0;
    for (uint32_t i = 0; i < num_threads; i++)
    {
        uint32_t numa_node = numa_node_of_cpu(cpu_ids[i]);
        threads.push_back(thread(
                Phase1<K>::ThreadB<table_index>,
                cpu_ids[i],
                &coordinator,
                static_cast<uint32_t*>(d_new_entry_positions),
                prev_bucket_indexes,
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
    cout << " (" << time(NULL) - start_seconds << "s)" << endl;

// Unfortunately single-threaded, but not much can be done to help it
    //new_line_point_bucket_index->SumBucketOffsets();

// Setup the new park list and buffer
    TemporaryPark<line_point_delta_len_bits> test_park(LinePointEntryUIDBucketEntry<K>::max_entries_per_bucket);
    buffers.push_back(new Buffer(test_park.GetSpaceNeeded() * LinePointEntryUIDBucketEntry<K>::num_buckets));
    graph_parks.push_back(vector<TemporaryPark<line_point_delta_len_bits>*>(LinePointEntryUIDBucketEntry<K>::num_buckets));

    coordinator = 0;
    cout << "Part C"<< (uint32_t)table_index;
    start_seconds = time(NULL);
    threads.clear();
    for (uint32_t i = 0; i < num_threads; i++)
    {
        threads.push_back(thread(
                Phase1<K>::ThreadC<table_index>,
                cpu_ids[i],
                &coordinator,
                static_cast<uint32_t*>(d_new_entry_positions),
                new_line_point_bucket_indexes,
                &(graph_parks[table_index]),
                buffers[table_index]
                ));
    }

    for (auto &it: threads)
    {
        it.join();
    }
    cout << " (" << time(NULL) - start_seconds << "s)" << endl;
    for (auto& numa_node : GetNUMANodesFromCpuIds(cpu_ids)) {
        delete new_line_point_bucket_indexes[numa_node];
    }
    return next_bucket_indexes;
}

template <uint8_t K>
Phase1<K>::Phase1(const uint8_t* id, uint32_t num_threads_in)
{
    for (uint32_t cpuid = 0; cpuid < num_threads_in; cpuid++)
    {
        cpu_ids.push_back(cpuid);
    }


    uint64_t phase_start_seconds = time(nullptr);

    num_threads = num_threads_in;

    map<uint32_t, BucketPageIndex<YCBucketEntry<K, -1>>*> first_bucket_indexes;
    for (auto& numa_node : GetNUMANodesFromCpuIds(cpu_ids))
    {
        first_bucket_indexes[numa_node] = new BucketPageIndex<YCBucketEntry<K, -1>>();
        first_bucket_indexes[numa_node]->pops_required_per_bucket = 2;
    }

    cout << "Part A";
    uint64_t part_start_seconds = time(nullptr);
    vector<thread> threads;
    std::atomic<uint64_t> coordinator = 0;
    for (uint32_t i = 0; i < num_threads; i++)
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

    TemporaryPark<finaltable_y_delta_len_bits> test_park(YCBucketEntry<K, 5>::max_entries_per_bucket);
    buffers.push_back(new Buffer(test_park.GetSpaceNeeded() * YCBucketEntry<K, 5>::num_buckets));
    final_parks.resize(YCBucketEntry<K, 5>::num_buckets);

    cout << "Part D";
    part_start_seconds = time(nullptr);
    threads.clear();
    coordinator = 0;
    for (uint32_t i = 0; i < num_threads; i++)
    {
        threads.push_back(thread(
                Phase1<K>::ThreadD,
                cpu_ids[i],
                &coordinator,
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
template class Phase1<32>;
template class Phase1<33>;