#include <vector>
#include <atomic>
#include <cstring>
#include <ctime>

#include "penguin.hpp"
#include "buffer.hpp"
#include "calculate_bucket.hpp"
#include "encoding.hpp"
#include "plotter.hpp"
#include "pos_constants.hpp"
#include "thread_mgr.hpp"
#include "bitcopy.hpp"

#include "status_update.hpp"

using namespace std;

template<PlotConf conf>
Penguin<FwdYCEntry<conf, -1>, conf.interlace_factor>* Plotter<conf>::Context::createFirstTable()
{
    auto output_penguin = new Penguin<FwdYCEntry<conf, -1>, conf.interlace_factor>();

    vector<thread> threads;

    for (auto & cpu_id : cpu_ids)
    {
        threads.push_back(thread( [this, cpu_id, output_penguin]
                                  {
            PinToCpuid(cpu_id);
            uint64_t batch_size = (1ULL<<kBatchSizes);
            vector<uint64_t> buff(batch_size);
            F1Calculator f1(conf.K, this->plotter->id.data());
            while (true)
            {
              uint64_t x = this->plotter->coordinator.fetch_add(batch_size);
              if (x >= (1ULL<<conf.K))
                  break;
              if ((x+batch_size) > (1ULL<<conf.K))
              {
                  batch_size = (1ULL<<conf.K) - x;
              }

              f1.CalculateBuckets(x, batch_size, buff.data());

              for (uint64_t i = 0; i < batch_size; i++)
              {
                  FwdYCEntry<conf, -1> entry;
                  entry.setY(buff[i]);
                  uint32_t c = bswap_32(x+i);
                  memcpy(entry.c, &c, 4);
                  entry.gid = x+i;
                  output_penguin->addEntry(entry);
              }
            }}));
    }
    for (auto &it: threads)
    {
        it.join();
    }

    return output_penguin;
}

template<PlotConf conf, uint8_t table_index>
class YCPenguinUnloader
{
    using entry_type = FwdYCEntry<conf, table_index-1>;

    // These are BC buckets from the chia algo


    static constexpr uint32_t max_entries_per_bc_bucket = 350;

    vector<vector<entry_type>> temp_buckets;

    int64_t latest_bucket_started_load = -1;
    int64_t latest_bucket_finished_load = -1;
    int64_t latest_row_loaded = -1;

    const vector<Penguin<entry_type, conf.interlace_factor>*> penguins;

public:
    static constexpr uint32_t bc_bucket_num = (1ULL << (conf.K + kExtraBits)) / kBC;
    static constexpr uint32_t temp_buckets_needed = (conf.GetMaxY(table_index)/conf.num_rows) / kBC + 12;

    explicit YCPenguinUnloader(const vector<Penguin<entry_type, conf.interlace_factor>*> penguins_in)
    : penguins(penguins_in)
    {
        temp_buckets.resize(temp_buckets_needed);
        for (auto &i : temp_buckets)
        {
            i.reserve(max_entries_per_bc_bucket);
        }
    }

    void loadRow(uint32_t row_id)
    {

        uint64_t first_value_contained = entry_type::getFirstYInRow(row_id);
        uint64_t first_bucket_needed = first_value_contained/kBC;

        uint64_t last_value_contained = entry_type::getFirstYInRow(row_id+1) - 1;
        uint64_t last_bucket_needed = last_value_contained/kBC;

        uint64_t clear_start = latest_bucket_started_load+1;
        if (clear_start < first_bucket_needed)
        {
            clear_start = first_bucket_needed;
        }

        // Clear temp buckets we will need
        for (uint32_t i = clear_start; i <= last_bucket_needed; i++)
        {
            temp_buckets[i%temp_buckets_needed].clear();
        }
        latest_bucket_started_load = last_bucket_needed;

        // Load in the row
        for (auto& penguin : penguins)
        {
            for (uint32_t entry_id = 0; entry_id < penguin->getCountInRow(row_id); entry_id++)
            {
                auto entry = penguin->readEntry(row_id, entry_id);
                uint32_t bucket_of_entry = entry.getY()/kBC;
                temp_buckets[bucket_of_entry%temp_buckets_needed].push_back(entry);
            }
        }

        latest_bucket_finished_load = (last_value_contained+1)/kBC - 1;
        latest_row_loaded = row_id;

        if (conf.GetMaxY(table_index-1) <= last_value_contained)
        {
            latest_bucket_finished_load = bc_bucket_num;
        }

        // If this row only contains entries relevant to this batch of buckets, pop it
        /*
        if ((first_bucket_needed > first_bucket_id_in_batch) &&
            (last_bucket_needed < (first_bucket_id_in_batch + batchSize)))
        {
            for (auto& [numa_node, penguin] : *prev_penguins) {
                penguin->popRow(row_id);
            }
        }
         */

        // TODO: pop the borderline entries somehow
    }

    vector<entry_type>& getBucketEntries(uint32_t bucket_id)
    {
        while (latest_bucket_finished_load < (int32_t)(bucket_id))
        {
            // Load next row
            uint32_t row_id = latest_row_loaded + 1;

            if (row_id >= conf.num_rows)
            {
                latest_bucket_finished_load = bucket_id+1;
                break;
            }

            uint64_t first_value_needed = bucket_id*kBC;
            uint32_t first_row_needed = FwdYCEntry<conf, table_index - 1>::getRowFromY(first_value_needed);
            if (row_id < first_row_needed)
            {
                row_id = first_row_needed;
            }

            loadRow(row_id);
        }
        for (auto &it : temp_buckets[bucket_id%temp_buckets_needed])
        {
            assert(it.getY()/kBC == bucket_id);
        }
        return temp_buckets[bucket_id%temp_buckets_needed];
    }
};

template<PlotConf conf, const uint8_t table_index>
class MatchQueue {
    using in_bucket_type = FwdYCEntry<conf, table_index - 1>;
    using out_bucket_type = FwdYCEntry<conf, table_index>;
    using out_gid_type = FwdGIDEntry<conf, table_index>;

    // Bounds for the bucket id of the entry with the larger UID
    uint32_t min_bucket_id = 0;
    uint32_t max_bucket_id = -1;

    // This is a map of matches between an entry with a larger UID and entries with smaller UIDs
    // Size here will be dynamically adjusted if necessary
    vector<vector<in_bucket_type>> match_staging;
    map<in_bucket_type, uint32_t, FwdYCEntryCompare<conf, table_index - 1>> match_staging_map;

    FxCalculator *fx;

public:

    Penguin<out_bucket_type, conf.interlace_factor> *output_bucket_penguin;
    Penguin<out_gid_type, conf.interlace_factor> *output_gid_penguin;

    MatchQueue(
            Penguin<out_bucket_type, conf.interlace_factor> *output_bucket_penguin_in,
            Penguin<out_gid_type, conf.interlace_factor> *output_gid_penguin_in) :
            output_bucket_penguin(output_bucket_penguin_in),
            output_gid_penguin(output_gid_penguin_in) {
        fx = new FxCalculator(conf.K, table_index + 2);
        match_staging.resize(300);
        for (auto &it : match_staging) {
            it.reserve(conf.max_gid_stub_val);
        }
    }

    ~MatchQueue() {
        delete fx;
    }

    void setBounds(uint32_t min_bucket_id_in, uint32_t max_bucket_id_in) {
        min_bucket_id = min_bucket_id_in;
        max_bucket_id = max_bucket_id_in;
    }

    void addMatch(in_bucket_type a, in_bucket_type b) {
        auto bigger_entry = a;
        auto smaller_entry = b;
        if (bigger_entry.gid < smaller_entry.gid) {
            swap(bigger_entry, smaller_entry);
        }

        if ((bigger_entry.getY() / kBC < min_bucket_id) || (bigger_entry.getY() / kBC > max_bucket_id)) {
            // Don't use matches that are out of bounds.
            return;
        }

        // Find if we have an output queue for this one yet
        int32_t dest_queue = -1;
        if (match_staging_map.contains(bigger_entry)) {
            dest_queue = match_staging_map[bigger_entry];
        } else {
            // Find an open staging vector
            for (uint32_t i = 0; i < match_staging.size(); i++) {
                if (match_staging[i].empty()) {
                    dest_queue = i;
                    break;
                }
            }
            if (dest_queue == -1) {
                // Out of queues! create another
                // Expensive, so hopefully almost never needed
                dest_queue = match_staging.size();
                match_staging.push_back(vector<in_bucket_type>());
            }

            // Cache so we know to put others here as well
            match_staging_map[bigger_entry] = dest_queue;
        }

        if (match_staging[dest_queue].size() >= conf.max_gid_stub_val) {
            throw runtime_error("Too many matches to index with gid stubs!");
        }

        match_staging[dest_queue].push_back(smaller_entry);
    }

    void sendMatchOut(
            in_bucket_type bigger_gid_entry,
            in_bucket_type smaller_gid_entry,
            uint8_t match_id
    ) {
        auto left_entry = bigger_gid_entry;
        auto right_entry = smaller_gid_entry;
        // Flip based on y
        if (left_entry.getY() > right_entry.getY()) {
            swap(left_entry, right_entry);
        }

        // Setup gid entry first
        auto new_gid_entry = FwdGIDEntry<conf, table_index>();
        new_gid_entry.setY(bigger_gid_entry.gid * conf.max_gid_stub_val + match_id);
        new_gid_entry.right_gid = smaller_gid_entry.gid;
        output_gid_penguin->addEntry(new_gid_entry);

        // Now setup bucket entry for next iteration
        const int64_t c_len_in = conf.GetCLen(table_index - 1);
        const int64_t c_len_out = conf.GetCLen(table_index);

        uint8_t input_data[64];
        // Zero out so we can or in easily
        memset(input_data, 0, 64);

        // Copy in R
        bitCopy<conf.K + kExtraBits + c_len_in, 64,
                in_bucket_type::c_len_bytes * 8 - c_len_in, in_bucket_type::c_len_bytes>(input_data, right_entry.c);

        // Copy in L
        bitCopy<conf.K + kExtraBits, 64, in_bucket_type::c_len_bytes * 8 - c_len_in, in_bucket_type::c_len_bytes>(
                input_data, left_entry.c);

        // Copy in F1
        {
            uint64_t f = bswap_64(left_entry.getY());
            bitCopy<0, 64,
                    8 * 8 - (conf.K + kExtraBits), 8>(input_data, (uint8_t *) &f);
        }

        blake3_hasher hasher;
        blake3_hasher_init(&hasher);
        blake3_hasher_update(&hasher,
                             input_data, (conf.K + kExtraBits + c_len_in * 2 + 7) / 8);

        uint8_t hash_bytes[32];
        blake3_hasher_finalize(&hasher, hash_bytes, sizeof(hash_bytes));

        FwdYCEntry<conf, table_index> new_entry;
        new_entry.gid = new_gid_entry.getY();

        if constexpr (table_index < 5) {
            new_entry.setY(bswap_64(*(uint64_t *) hash_bytes) >> (64 - (conf.K + kExtraBits)));
        } else {
            new_entry.setY(bswap_64(*(uint64_t *) hash_bytes) >> (64 - (conf.K)));
        }

        memset(new_entry.c, 0, new_entry.c_len_bytes);
        if constexpr (table_index < 2) {
            // Copy in R
            bitCopy<out_bucket_type::c_len_bytes * 8 - c_len_in, out_bucket_type::c_len_bytes,
                    in_bucket_type::c_len_bytes * 8 - c_len_in, in_bucket_type::c_len_bytes>(new_entry.c,
                                                                                             right_entry.c);

            // Copy in L
            bitCopy<out_bucket_type::c_len_bytes * 8 - c_len_in * 2, out_bucket_type::c_len_bytes,
                    in_bucket_type::c_len_bytes * 8 - c_len_in, in_bucket_type::c_len_bytes>(new_entry.c, left_entry.c);
        } else if constexpr (table_index < 5) {
            // Copy from hash result
            bitCopy<out_bucket_type::c_len_bytes * 8 - c_len_out, out_bucket_type::c_len_bytes,
                    conf.K + kExtraBits, 32>(new_entry.c, hash_bytes);
            assert(*new_entry.c < (1ULL << (c_len_out - (out_bucket_type::c_len_bytes - 1) * 8)));
        }
        uint64_t entry_i = output_bucket_penguin->addEntry(new_entry);

        if (new_entry.gid == 82421925)
        {
            cout << "Found" << endl;
        }

        FwdYCEntry<conf, table_index> test_entry = output_bucket_penguin->readEntry(new_entry.getRow(), entry_i);
        assert(new_entry.getY() == test_entry.getY());
        assert(new_entry.gid == test_entry.gid);
        assert(memcmp(new_entry.c, test_entry.c, test_entry.c_len_bytes) == 0);

        // Compare with stock calculation system
        uint64_t c_len = kVectorLens[table_index + 2] * conf.K;
        uint128_t cl = 0;
        memcpy(&cl, left_entry.c, left_entry.c_len_bytes);
        cl = bswap_128(cl) >> (128 - left_entry.c_len_bytes*8);
        uint128_t cr = 0;
        memcpy(&cr, right_entry.c, left_entry.c_len_bytes);
        cr = bswap_128(cr) >> (128 - left_entry.c_len_bytes*8);
        auto out = fx->CalculateBucket(
                Bits(left_entry.getY(), conf.K+kExtraBits),
                Bits(cl, c_len),
                Bits(cr, c_len));
        uint128_t stock_y = 0;
        uint128_t stock_c = 0;
        if (table_index < 5) {
            stock_y = out.first.Slice(0, conf.K+kExtraBits).GetValue();
            uint8_t buff[32];
            memset(buff, 0, sizeof(buff));
            out.second.ToBytes(buff);
            stock_c = Util::SliceInt128FromBytes(
                    buff, 0, kVectorLens[table_index + 3] * conf.K);
        }
        else
        {
            stock_y = out.first.Slice(0, conf.K).GetValue();
        }
        uint128_t c_other = 0;
        memcpy(&c_other, test_entry.c, test_entry.c_len_bytes);
        c_other = bswap_128(c_other) >> (128 - test_entry.c_len_bytes*8);
        assert(test_entry.getY() == stock_y);
        assert(c_other == stock_c);

    }

    uint64_t flushQueueForEntries(vector<in_bucket_type> entries)
    {
        uint64_t entries_sent = 0;
        // Scan for queues that should get dumped now
        for (auto & bigger_entry : entries) {
            if (!match_staging_map.contains(bigger_entry))
                continue;
            uint32_t queue_id = match_staging_map[bigger_entry];

            if (match_staging[queue_id].size() > 0) {
                // sort by GID
                struct {
                    bool operator()(in_bucket_type &a, in_bucket_type &b) const { return a.gid < b.gid; }
                } customLess;
                sort(match_staging[queue_id].begin(), match_staging[queue_id].end() - 1, customLess);
            }

            for (uint8_t i = 0; i < match_staging[queue_id].size(); i++)
            {
                // Dump this entry!
                sendMatchOut(bigger_entry, match_staging[queue_id][i], i);
                entries_sent++;
            }
            match_staging[queue_id].clear();
            match_staging_map.erase(bigger_entry);
        }
        return entries_sent;
    }
};



void assert_matching(uint64_t lout, uint64_t rout)
{

    int64_t bucket_id_lout = lout/kBC;
    int64_t bucket_id_rout = rout/kBC;
    assert(bucket_id_lout + 1 == bucket_id_rout);
    int64_t b_id_l = (lout%kBC)/kC;
    int64_t c_id_l = (lout%kBC)%kC;
    int64_t b_id_r = (rout%kBC)/kC;
    int64_t c_id_r = (rout%kBC)%kC;
    int64_t m = (kB + b_id_r - b_id_l)%kB;
    assert((uint64_t)m < (1UL<<kExtraBits));
    int64_t a = c_id_r - c_id_l;
    int64_t b = 2*m + (bucket_id_lout%2);
    int64_t c = (b*b)%kC;
    int64_t d = (kC + a)%kC;
    assert(c == d);

}

template <PlotConf conf>
template <uint8_t table_index>
Penguin<FwdYCEntry<conf, table_index>, conf.interlace_factor>*  Plotter<conf>::Context::createTable(
        vector<Penguin<FwdYCEntry<conf, table_index-1>, conf.interlace_factor>*> prev_penguins,
        atomic<uint64_t>& num_matches_out)
{
    auto output_bucket_penguin = new Penguin<FwdYCEntry<conf, table_index>, conf.interlace_factor>();
    auto output_gid_penguin = new Penguin<FwdGIDEntry<conf, table_index>, conf.interlace_factor>();
    this->forward_pass_gid_penguins[table_index] = (void*)output_gid_penguin;

    vector<thread> threads;

    for (auto & cpu_id : cpu_ids)
    {
        threads.push_back(thread( [this, cpu_id, &output_bucket_penguin, &output_gid_penguin, &prev_penguins, &num_matches_out] {

            PinToCpuid(cpu_id);
            uint64_t num_matches_in_thread = 0;
            vector<vector<uint32_t>> right_map(kBC);
            for (auto &i : right_map) {
                i.clear();
            }
            std::vector<uint32_t> right_map_clean;

            YCPenguinUnloader<conf, table_index> penguin_unloader(prev_penguins);
            MatchQueue<conf, table_index> match_queue(
                    output_bucket_penguin,
                    output_gid_penguin
                    );

            constexpr uint32_t batchSize = YCPenguinUnloader<conf, table_index>::temp_buckets_needed * 32;
            uint64_t bucket_id = 0;
            while (bucket_id < YCPenguinUnloader<conf, table_index>::bc_bucket_num)
            {

                uint64_t first_bucket_id_in_batch = this->plotter->coordinator.fetch_add(batchSize);
                uint64_t last_bucket_id_in_batch = bucket_id+batchSize;

                if (last_bucket_id_in_batch > YCPenguinUnloader<conf, table_index>::bc_bucket_num-1)
                {
                    last_bucket_id_in_batch = YCPenguinUnloader<conf, table_index>::bc_bucket_num-1;
                }

                bucket_id = first_bucket_id_in_batch;

                match_queue.setBounds(first_bucket_id_in_batch, last_bucket_id_in_batch);

                if (bucket_id > 0)
                {
                    // Start loading from one bucket back
                    bucket_id--;
                }
                while (bucket_id < last_bucket_id_in_batch)
                {
                    auto left_bucket_entries = penguin_unloader.getBucketEntries(bucket_id);
                    auto right_bucket_entries = penguin_unloader.getBucketEntries(bucket_id + 1);

                    // Process bucket_id and bucket_id + 1
                    uint16_t parity = bucket_id % 2;
                    for (size_t yl : right_map_clean) {
                        right_map[yl].clear();
                    }
                    right_map_clean.clear();
                    for (size_t pos_R = 0; pos_R < right_bucket_entries.size(); pos_R++) {
                        auto right_entry = right_bucket_entries[pos_R];
                        uint64_t r_y = right_entry.getY() % kBC;
                        if (right_map[r_y].empty())
                        {
                            right_map_clean.push_back(r_y);
                        }
                        right_map[r_y].push_back(pos_R);
                    }

                    for (size_t pos_L = 0; pos_L < left_bucket_entries.size(); pos_L++) {
                        auto left_entry = left_bucket_entries[pos_L];
                        for (uint8_t i = 0; i < kExtraBitsPow; i++) {
                            uint16_t r_target = L_targets[parity][left_entry.getY()%kBC][i];
                            for (auto& pos_R : right_map[r_target]) {
                                match_queue.addMatch(left_entry, right_bucket_entries[pos_R]);
                                assert_matching(left_entry.getY(), right_bucket_entries[pos_R].getY());
                            }
                        }
                    }

                    num_matches_in_thread += match_queue.flushQueueForEntries(left_bucket_entries);

                    if (bucket_id == YCPenguinUnloader<conf, table_index>::bc_bucket_num-2)
                    {
                        // We have emptied the last bucket. This means we should also dump the right bucket from the queue
                        num_matches_in_thread += match_queue.flushQueueForEntries(right_bucket_entries);
                    }
                    bucket_id++;
                }
            }

            num_matches_out += num_matches_in_thread;
        }));
    }
    for (auto &it: threads)
    {
        it.join();
    }

    return output_bucket_penguin;
}

template <PlotConf conf, uint8_t table_index>
vector<Penguin<FwdYCEntry<conf, table_index>, conf.interlace_factor>*> DoTable(
        Plotter<conf>* plotter,
        vector<Penguin<FwdYCEntry<conf, table_index-1>, conf.interlace_factor>*> input_yc_penguins)
{
    StatusUpdate::StartSeg("0." + to_string((uint32_t)table_index) + "0M");
    vector<thread> threads;
    plotter->coordinator = 0;
    atomic<uint64_t> num_matches = 0;
    vector<Penguin<FwdYCEntry<conf, table_index>, conf.interlace_factor>*> output_yc_penguin(plotter->contexts.size());
    for (uint32_t i = 0; i < plotter->contexts.size(); i++)
    {
        threads.push_back(thread( [plotter, i, &output_yc_penguin, &input_yc_penguins, &num_matches] {
            output_yc_penguin[i] = plotter->contexts[i]->template createTable<table_index>(input_yc_penguins, num_matches);
        }));
    }
    for (auto & thread : threads)
    {
        thread.join();
    }
    cout << " (" << num_matches << " matches)";
    StatusUpdate::StartSeg("0." + to_string((uint32_t)table_index) + "1S");
    for (auto & penguin : input_yc_penguins)
    {
        delete penguin;
    }

    return output_yc_penguin;
}

template <PlotConf conf>
void Plotter<conf>::create(array<uint8_t, 32> id_in)
{
    id = id_in;

    load_tables();

    vector<thread> threads;
    coordinator = 0;
    vector<Penguin<FwdYCEntry<conf, -1>, conf.interlace_factor>*> first_yc_penguin(contexts.size());
    StatusUpdate::StartSeg("0.-1M");
    for (uint32_t i = 0; i < contexts.size(); i++)
    {
        threads.push_back(thread( [this,
                                   i,
                                   &first_yc_penguin] {
            first_yc_penguin[contexts[i]->numa_node] = contexts[i]->createFirstTable();
        }));
    }
    for (auto & thread : threads)
    {
        thread.join();
    }

    auto i1 = DoTable<conf, 0>(this, first_yc_penguin);
    auto i2 = DoTable<conf, 1>(this, i1);
    auto i3 = DoTable<conf, 2>(this, i2);
    auto i4 = DoTable<conf, 3>(this, i3);
    auto i5 = DoTable<conf, 4>(this, i4);
    auto i6 = DoTable<conf, 5>(this, i5);

    for (uint32_t i = 0; i < contexts.size(); i++)
    {
        contexts[i]->forward_pass_final_yc_penguin = i6[i];
    }

    StatusUpdate::EndSeg();
}

#include "explicit_templates.hpp"