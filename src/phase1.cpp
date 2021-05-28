
#include <vector>
#include <atomic>
#include <cstring>

#include "bucket_page_index.hpp"
#include "buffer.hpp"
#include "calculate_bucket.hpp"
#include "encoding.hpp"
#include "park.hpp"
#include "phase1.hpp"
#include "pos_constants.hpp"

using namespace std;

template <uint8_t K>
void Phase1<K>::ThreadA(
        BucketPageIndex<YCBucketEntry<K, -1>>* new_bucket_index
)
{
    F1Calculator f1(K, id);
    while (true)
    {
        uint64_t batch_size = (1ULL<<kBatchSizes);
        uint64_t x = coordinator.fetch_add(batch_size);
        if (x >= (1ULL<<K))
            break;
        if ((x+batch_size) > (1ULL<<K))
        {
            batch_size = (1ULL<<K) - x;
        }
        uint64_t buff[batch_size];
        f1.CalculateBuckets(x, batch_size, buff);

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
        BucketPageIndex<YCBucketEntry<K, table_index-1>>* prev_bucket_index,
        BucketPageIndex<YCBucketEntry<K, table_index>>* new_bucket_index,
        BucketPageIndex<LinePointEntryUIDBucketEntry<K>>* new_line_point_bucket_index)
{
    vector<vector<uint32_t>> right_map(kBC);
    std::vector<uint32_t> right_map_clean;
    FxCalculator fx(K, table_index+2);
    while (true)
    {
        uint64_t bucket_id = coordinator.fetch_add(1);
        if (bucket_id >= YCBucketEntry<K, table_index>::num_buckets)
            break;

        uint16_t parity = bucket_id % 2;
        if (bucket_id < YCBucketEntry<K, table_index>::num_buckets-1)
        {
            for (size_t yl : right_map_clean) {
                right_map[yl].clear();
            }
            right_map_clean.clear();
            for (size_t pos_R = 0; pos_R < prev_bucket_index->GetCountInBucket(bucket_id+1); pos_R++) {
                uint64_t r_y = prev_bucket_index->ReadEntry(bucket_id+1, pos_R).y%kBC;
                if (right_map[r_y].empty())
                {
                    right_map_clean.push_back(r_y);
                }
                right_map[r_y].push_back(pos_R);
            }
        }
        for (size_t pos_L = 0; pos_L < prev_bucket_index->GetCountInBucket(bucket_id); pos_L++) {
            auto left_entry = prev_bucket_index->ReadEntry(bucket_id, pos_L);
            if (bucket_id < YCBucketEntry<K, table_index>::num_buckets-1) {
                for (uint8_t i = 0; i < kExtraBitsPow; i++) {
                    uint16_t r_target = L_targets[parity][left_entry.y%kBC][i];
                    for (auto& pos_R : right_map[r_target]) {
                        auto right_entry = prev_bucket_index->ReadEntry(bucket_id+1, pos_R);

                        // Calculate the next F and dump a new entry into the bucket index for the next table
                        uint64_t c_len = kVectorLens[table_index + 2] * K;
                        auto out = fx.CalculateBucket(
                                Bits(left_entry.y, K + kExtraBits),
                                Bits(left_entry.c, c_len),
                                Bits(right_entry.c, c_len));

                        YCBucketEntry<K, table_index> new_entry;
                        new_entry.y = out.first.GetValue();
                        if (table_index < 5) {
                            uint8_t buff[16];
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
                            x = new_entry_positions[prev_bucket_index->GetUniqueIdentifier(bucket_id, pos_L)];
                            y = new_entry_positions[prev_bucket_index->GetUniqueIdentifier(bucket_id+1, pos_R)];
                        }
                        LinePointEntryUIDBucketEntry<K> line_point_entry;
                        line_point_entry.line_point = Encoding::SquareToLinePoint(x, y);
                        line_point_entry.entry_uid = new_bucket_index->GetUniqueIdentifier(bucket_id, new_entry_id);
                        assert(line_point_entry.entry_uid < (1ULL << (K + 2)));
                        new_line_point_bucket_index->InsertEntry(line_point_entry);
                    }
                }
            }
        }
        prev_bucket_index->PopBucket(bucket_id);
        prev_bucket_index->PopBucket(bucket_id+1);
    }
}

template <uint8_t K>
template <int8_t table_index>
void Phase1<K>::ThreadC(
        BucketPageIndex<LinePointEntryUIDBucketEntry<K>>* line_point_bucket_index)
{
    while (true)
    {
        uint64_t bucket_id = coordinator.fetch_add(1);
        if (bucket_id >= LinePointEntryUIDBucketEntry<K>::num_buckets)
            break;
// Sort all entries in this bucket
        LinePointEntryUIDBucketEntry<K> entries[LinePointEntryUIDBucketEntry<K>::max_entries_per_bucket];
        uint32_t num_entries = line_point_bucket_index->GetCountInBucket(bucket_id);
        for (uint32_t i = 0; i < num_entries; i++)
        {
            entries[i] = line_point_bucket_index->ReadEntry(bucket_id, i);
        }
// Done with the bucket now
        line_point_bucket_index->PopBucket(bucket_id);

// Sort the bucket
        struct {
            bool operator()(LinePointEntryUIDBucketEntry<K> a, LinePointEntryUIDBucketEntry<K> b) const { return a.line_point < b.line_point; }
        } customLess;
        sort(entries, entries + num_entries, customLess);
// We have a sorted entries list! Emit to output buffer and record new positions
        auto new_park = new TemporaryPark<line_point_delta_len_bits>(num_entries);
        uint64_t num_bytes = new_park->GetSpaceNeeded();
        new_park->Bind(buffers[table_index]->data + buffers[table_index]->GetInsertionOffset(num_bytes));
        uint128_t line_points[LinePointEntryUIDBucketEntry<K>::max_entries_per_bucket];
        for (uint32_t i = 0; i < num_entries; i++)
        {
            line_points[i] = entries[i].line_point;
            new_entry_positions[entries[i].entry_uid] = i + line_point_bucket_index->GetBucketOffset(bucket_id);
        }
        new_park->AddEntries(line_points);
        graph_parks[table_index][bucket_id] = new_park;
    }
}

template <uint8_t K>
void Phase1<K>::ThreadD(
        BucketPageIndex<YCBucketEntry<K, 5>>* line_point_bucket_index)
{
    while (true)
    {
        uint64_t bucket_id = coordinator.fetch_add(1);
        if (bucket_id >= YCBucketEntry<K, 6>::num_buckets)
            break;
// Sort all entries in this bucket
        YCBucketEntry<K, 5> entries[YCBucketEntry<K, 5>::max_entries_per_bucket];
        uint32_t num_entries = line_point_bucket_index->GetCountInBucket(bucket_id);
        for (uint32_t i = 0; i < num_entries; i++)
        {
            entries[i] = line_point_bucket_index->ReadEntry(bucket_id, i);
        }
// Done with the bucket now
        line_point_bucket_index->PopBucket(bucket_id);

// Sort the bucket
        struct {
            bool operator()(YCBucketEntry<K, 5> a, YCBucketEntry<K, 5> b) const { return a.y < b.y; }
        } customLess;
        sort(entries, entries + num_entries, customLess);

// We have a sorted entries list! Emit to output buffer and record new positions
        auto new_park = new TemporaryPark<finaltable_y_delta_len_bits>(num_entries);
        uint64_t num_bytes = new_park->GetSpaceNeeded();
        new_park->Bind(buffers[6]->data + buffers[6]->GetInsertionOffset(num_bytes));
        uint128_t line_points[YCBucketEntry<K, 6>::max_entries_per_bucket];
        for (uint32_t i = 0; i < num_entries; i++)
        {
            line_points[i] = entries[i].y;
        }
        new_park->AddEntries(line_points);
        final_parks[bucket_id] = new_park;
    }
}

template <uint8_t K>
template <int8_t table_index>
BucketPageIndex<YCBucketEntry<K, table_index>>* Phase1<K>::DoTable(
        BucketPageIndex<YCBucketEntry<K, table_index-1>>* prev_bucket_index)
{
    auto next_bucket_index = new BucketPageIndex<YCBucketEntry<K, table_index>>();
    next_bucket_index->pops_required_per_bucket = 2;
    BucketPageIndex<LinePointEntryUIDBucketEntry<K>> new_line_point_bucket_index;

    cout << "Started B "<< table_index << endl;
    vector<thread> threads;
    coordinator = 0;
    for (uint32_t i = 0; i < num_threads; i++)
    {
        threads.push_back(thread([&]{ ThreadB<table_index>(
                prev_bucket_index,
                next_bucket_index,
                &new_line_point_bucket_index);}));
    }

    for (auto &it: threads)
    {
        it.join();
    }
    delete prev_bucket_index;

// Unfortunately single-threaded, but not much can be done to help it
    new_line_point_bucket_index.SumBucketOffsets();

// Setup the new park list and buffer
    TemporaryPark<line_point_delta_len_bits> test_park(LinePointEntryUIDBucketEntry<K>::max_entries_per_bucket);
    buffers.push_back(new Buffer(test_park.GetSpaceNeeded() * LinePointEntryUIDBucketEntry<K>::num_buckets));
    graph_parks.push_back(vector<TemporaryPark<line_point_delta_len_bits>*>(LinePointEntryUIDBucketEntry<K>::num_buckets));

    coordinator = 0;
    cout << "Started C "<< table_index << endl;
    threads.clear();
    for (uint32_t i = 0; i < num_threads; i++)
    {
        threads.push_back(thread([&]{ ThreadC<table_index>(&new_line_point_bucket_index);}));
    }

    for (auto &it: threads)
    {
        it.join();
    }
    return next_bucket_index;
}

template <uint8_t K>
Phase1<K>::Phase1(const uint8_t* id_in, uint32_t num_threads_in)
{
    id = id_in;
    num_threads = num_threads_in;

    auto first_bucket_index = new BucketPageIndex<YCBucketEntry<K, -1>>();

    new_entry_positions.resize(1ULL<<(K+2));

    cout << "Started A"<< endl;
    vector<thread> threads;
    coordinator = 0;
    for (uint32_t i = 0; i < num_threads; i++)
    {
        threads.push_back(thread([&]{ ThreadA(first_bucket_index);}));
    }

    for (auto &it: threads)
    {
        it.join();
    }

    auto i0 = DoTable<0>(first_bucket_index);
    auto i1 = DoTable<1>(i0);
    auto i2 = DoTable<2>(i1);
    auto i3 = DoTable<3>(i2);
    auto i4 = DoTable<4>(i3);
    auto i5 = DoTable<5>(i4);

    cout << "Started D"<< endl;
    threads.clear();
    coordinator = 0;
    for (uint32_t i = 0; i < num_threads; i++)
    {
        threads.push_back(thread([&]{ ThreadD(i5);}));
    }
    for (auto &it: threads)
    {
        it.join();
    }

}

template class Phase1<18>;
template class Phase1<22>;
template class Phase1<32>;
template class Phase1<33>;