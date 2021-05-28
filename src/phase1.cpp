
#include <vector>
#include <atomic>

#include "bucket_page_index.hpp"
#include "buffers.hpp"
#include "calculate_bucket.hpp"
#include "encoding.hpp"
#include "park.hpp"
#include "phase1.hpp"
#include "pos_constants.hpp"

using namespace std;

const uint8_t finaltable_y_delta_len_bits = 4;

template <uint8_t K>
static void Phase1<K>::ThreadA(
        atomic<uint64_t>& coordinator,
        const uint8_t* id,
        BucketPageIndex<YCBucketEntry<K, -1>>& new_bucket_index
)
{
    F1Calculator f1(K, id);
    while (1)
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
            YCBucketEntry<K, table_index> entry;
            entry.y = buff[i];
            entry.c = x+i;
            new_bucket_index.InsertEntry(entry);
        }
    }
}

template <uint8_t K>
template <int8_t table_index>
static void Phase1<K>::ThreadB(
        std::atomic<uint64_t>& coordinator,
        vector<uint32_t>& new_prev_bucket_positions,
        BucketPageIndex<YCBucketEntry<K, table_index-1>>& prev_bucket_index,
        BucketPageIndex<YCBucketEntry<K, table_index>>& new_bucket_index,
        BucketPageIndex<LinepointEntryUIDBucketEntry<K>>& new_linepoint_bucket_index)
{
    vector<vector<uint32_t>> rmap(kBC);
    std::vector<uint32_t> rmap_clean;
    FxCalculator fx(K, table_index+2);
    while (1)
    {
        uint64_t bucket_id = coordinator.fetch_add(1);
        if (bucket_id >= YCBucketEntry<K, table_index>::num_buckets)
            break;

        uint16_t parity = bucket_id % 2;
        if (bucket_id < YCBucketEntry<K, table_index>::num_buckets-1)
        {
            for (size_t yl : rmap_clean) {
                rmap[yl].clear();
            }
            rmap_clean.clear();
            for (size_t pos_R = 0; pos_R < prev_bucket_index.GetCountInBucket(bucket_id+1); pos_R++) {
                uint64_t r_y = prev_bucket_index.ReadEntry(bucket_id+1, pos_R).y;
                if (rmap[r_y].size() == 0)
                {
                    rmap_clean.push_back(r_y);
                }
                rmap[r_y].push_back(pos_R);
            }
        }
        for (size_t pos_L = 0; pos_L < prev_bucket_index.GetCountInBucket(bucket_id); pos_L++) {
            auto left_entry = prev_bucket_index.ReadEntry(bucket_id, pos_L);
            if (bucket_id < YCBucketEntry<K, table_index>::num_buckets-1) {
                for (uint8_t i = 0; i < kExtraBitsPow; i++) {
                    uint16_t r_target = L_targets[parity][left_entry.y_offset][i];
                    for (auto& pos_R : rmap[r_target]) {
                        auto right_entry = prev_bucket_index.ReadEntry(bucket_id+1, pos_R);

                        // Calculate the next F and dump a new entry into the bucket index for the next table
                        uint64_t c_len = kVectorLens[table_index + 2] * K;
                        auto out = fx.CalculateBucket(
                                Bits(left_entry.y, K + kExtraBits),
                                Bits(left_entry.c, c_len),
                                Bits(right_entry.c, c_len));

                        struct YCBucketedEntry new_entry;
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
                        uint64_t new_entry_id = new_bucket_index.InsertEntry(new_entry);


                        // Emit a new linepoint for sorting into the phase1 output buffer
                        uint64_t x, y;
                        if (table_index == 0) // First table just gets the two x values
                        {
                            x = left_entry.c;
                            y = right_entry.c;
                        }
                        else
                        {
                            // use the unique ids for the two values
                            x = new_prev_bucket_positions[prev_bucket_index.GetUniqueIdentifier(bucket_id, pos_L)];
                            y = new_prev_bucket_positions[prev_bucket_index.GetUniqueIdentifier(bucket_id+1, pos_R)];
                        }
                        LinepointEntryUIDBucketEntry<K> linepoint_entry;
                        linepoint_entry.linepoint = Encoding::SquareToLinePoint(x, y);
                        linepoint_entry.entry_uid = new_bucket_index.GetUniqueIdentifier(bucket_id, new_entry_id);
                        assert(linepoint_entry.entry_uid < (1ULL << (K+2)));
                        new_linepoint_bucket_index.InsertEntry(linepoint_entry);
                    }
                }
            }
        }
        prev_bucket_index.PopBucket(bucket_id);
        prev_bucket_index.PopBucket(bucket_id+1);
    }
}

template <uint8_t K>
template <int8_t table_index>
static void Phase1<K>ThreadC(
        std::atomic<uint64_t>& coordinator,
        vector<uint32_t>& new_entry_positions,
        BucketPageIndex<LinepointEntryUIDBucketEntry<K>>& linepoint_bucket_index,
        Buffer* output_buffer,
        vector<TemporaryPark<K+2>*>& temporary_parks_out)
{
while (1)
{
uint64_t bucket_id = coordinator.fetch_add(1);
if (bucket_id >= LinepointEntryUIDBucketEntry<K>::num_buckets)
break;
// Sort all entries in this bucket
LinepointEntryUIDBucketEntry<K> entries[LinepointEntryUIDBucketEntry<K>::max_entries_per_bucket];
uint32_t num_entries = linepoint_bucket_index.GetCountInBucket(bucket_id);
for (uint32_t i = 0; i < num_entries; i++)
{
entries[i] = linepoint_bucket_index.ReadEntry(bucket_id, i);
}
// Done with the bucket now
linepoint_bucket_index.PopBucket(bucket_id);

// Sort the bucket
struct {
    bool operator()(LinepointEntryUIDBucketEntry<K> a, LinepointEntryUIDBucketEntry<K> b) const { return a.y < b.y; }
} customLess;
sort(entries, entries + entry_len, customLess);
// We have a sorted entries list! Emit to output buffer and record new positions
auto new_park = new TemporaryPark<K>(num_entries);
uint64_t num_bytes = new_park.GetSpaceNeeded();
new_park.Bind(output_buffer->data + output_buffer->GetInsertionOffset(num_bytes));
uint128_t linepoints[LinepointEntryUIDBucketEntry<K>::max_entries_per_bucket];
for (uint32_t i = 0; i < num_entries; i++)
{
linepoints[i] = entries[i].linepoint;
new_entry_positions[entries[i].entry_uid] = i + linepoint_bucket_index.GetBucketOffset(bucket_id);
}
new_park.AddEntries(linepoints);
temporary_parks_out[bucket_id] = new_park;
}
};

template <uint8_t K>
static inline void Phase1ThreadD(
        std::atomic<uint64_t>& coordinator,
        vector<uint32_t>& new_entry_positions,
        BucketPageIndex<YCBucketEntry<K, 6>>& linepoint_bucket_index,
        Buffer* output_buffer,
        vector<TemporaryPark<finaltable_y_delta_len_bits>*>& temporary_parks_out
)
{
while (1)
{
uint64_t bucket_id = coordinator.fetch_add(1);
if (bucket_id >= BucketPageIndex<YCBucketEntry<K, 6>>::num_buckets)
break;
// Sort all entries in this bucket
YCBucketEntry<K, 6> entries[YCBucketEntry<K, 6>::max_entries_per_bucket];
uint32_t num_entries = linepoint_bucket_index.GetCountInBucket(bucket_id);
for (uint32_t i = 0; i < num_entries; i++)
{
entries[i] = linepoint_bucket_index.ReadEntry(bucket_id, i);
}
// Done with the bucket now
linepoint_bucket_index.PopBucket(bucket_id);

// Sort the bucket
struct {
    bool operator()(YCBucketEntry<K, 6> a, YCBucketEntry<K, 6> b) const { return a.y < b.y; }
} customLess;
sort(entries, entries + entry_len, customLess);

// We have a sorted entries list! Emit to output buffer and record new positions
auto new_park = new TemporaryPark<K>(num_entries);
uint64_t num_bytes = new_park.GetSpaceNeeded();
new_park.Bind(output_buffer->data + output_buffer->GetInsertionOffset(num_bytes));
uint128_t linepoints[YCBucketEntry<K, 6>::max_entries_per_bucket];
for (uint32_t i = 0; i < num_entries; i++)
{
linepoints[i] = entries[i].y;
}
new_park.AddEntries(linepoints);
}
}

template <uint8_t K, int8_t table_index> static inline void Phase1DoTable(
        uint32_t num_threads,
        void * prev_bucket_index_in,
        void ** next_bucket_index_out,
        vector<uint32_t>* prev_remap_table,
        vector<uint32_t>** next_remap_table,
        Buffer** new_buffer_out,
        vector<TemporaryPark<K+2>*>& temporary_parks_out,
)
{
auto prev_bucket_index = (BucketPageIndex<YCBucketEntry<K, table_index-1>>*)prev_bucket_index_in;
auto next_bucket_index = new BucketPageIndex<YCBucketEntry<K, table_index>>();
next_bucket_index->pops_required_per_bucket = 2;
auto new_linepoint_bucket_index = new BucketPageIndex<LinepointEntryUIDBucketEntry<K>>();

atomic<uint64_t> coordinator = 0;

cout << "Started B "<< table_index << endl;
vector<thread> threads;
for (uint32_t i = 0; i < num_threads; i++)
{
threads.push_back(thread(
        Phase1ThreadB<K, table_index>,
        ref(coordinator),
        prev_remap_table,
        prev_bucket_index,
        next_bucket_index,
        new_linepoint_bucket_index));
}

for (auto &it: threads)
{
it.join();
}
delete prev_bucket_index;
if (prev_remap_table != nullptr)
delete prev_remap_table;

// Unfortunately single-threaded, but not much can be done to help it
new_linepoint_bucket_index.SumBucketOffsets();

// Setup the new park list and buffer
TemporaryPark<K+2> test_park(LinepointEntryUIDBucketEntry<K>::max_entries_per_bucket);
*new_buffer_out = new Buffer(test_park.GetSpaceNeeded()*LinepointEntryUIDBucketEntry<K>::num_buckets);
temporary_parks_out.resize(LinepointEntryUIDBucketEntry<K>::num_buckets);
*next_remap_table = new vector<uint32_t>(1ULL<<(K+2));

coordinator = 0;
cout << "Started C "<< table_index << endl;
threads.clear();
for (uint32_t i = 0; i < num_threads; i++)
{
threads.push_back(thread(
        Phase1ThreadC<K, table_index>,
        ref(coordinator),
        *next_remap_table,
        new_linepoint_bucket_index,
        new_buffer_out,
        ref(temporary_parks_out)
));
}

for (auto &it: threads)
{
it.join();
}


*next_bucket_index_out = (void*)next_bucket_index;
}

template <uint8_t K> struct Phase1_Results
{
    vector<Buffer*> buffers;
    vector<vector<TemporaryPark<K+4>*>> graph_parks;
    vector<TemporaryPark<finaltable_y_delta_len_bits>*> final_parks;
}

template <uint8_t K> inline Phase1_Results<K> Phase1(const uint8_t* id, uint32_t num_threads)
{
    struct Phase1_Results<K> results;

    void * prev_bucket_index = (void *) new BucketPageIndex<YCBucketEntry<K, -1>>();

    atomic<uint64_t> coordinator = 0;
    cout << "Started A"<< endl;
    vector<thread> threads;
    for (uint32_t i = 0; i < num_threads; i++)
    {
        threads.push_back(thread(
                Phase1ThreadA<K>,
                ref(coordinator),
                id,
                ref(prev_bucket_index)));
    }

    for (auto &it: threads)
    {
        it.join();
    }

    vector<uint32_t>* prev_bucket_positions = nullptr;

    for (const int table_index = 0; table_index < 6; table_index++)
    {
        void * next_bucket_index;
        vector<uint32_t>* next_bucket_positions;
        Phase1DoTable<K, table_index>(
                num_threads,
                prev_bucket_index,
                ref(next_bucket_index),
                prev_bucket_positions,
                ref(next_bucket_positions),
                ref(results.buffers[table_index]),
                ref(results.graph_parks[table_index])
        );

        prev_bucket_index = next_bucket_index;
        prev_bucket_positions = next_bucket_positions;
    }

    coordinator = 0;
    cout << "Started D"<< endl;
    threads.clear();
    for (uint32_t i = 0; i < num_threads; i++)
    {
        threads.push_back(thread(
                Phase1ThreadD<K>,
                ref(coordinator),
                ref(prev_bucket_positions),
                ref(prev_bucket_index),
                ref(results.buffers[6]),
                ref(results.final_parks)));
    }
    for (auto &it: threads)
    {
        it.join();
    }

    return results;
}