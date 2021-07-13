#include <vector>
#include <atomic>
#include <cstring>
#include <ctime>
#include <condition_variable>

#include "penguin.hpp"
#include "buffer.hpp"
#include "calculate_bucket.hpp"
#include "encoding.hpp"
#include "plotter.hpp"
#include "pos_constants.hpp"
#include "thread_mgr.hpp"
#include "bitcopy.hpp"
#include "status_update.hpp"
#include "penguin_unloaders.hpp"

using namespace std;

template<PlotConf conf>
Penguin<FwdYCEntry<conf, -1>, true>* Plotter<conf>::Context::createFirstTable()
{
    auto output_penguin = new Penguin<FwdYCEntry<conf, -1>, true>();

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
                  auto entry = output_penguin->newEntry(buff[i]);
                  uint64_t c = bswap_64(x+i);
                  memcpy((void*)entry.data->c, ((uint8_t*)&c) + 8 - entry.c_len_bytes, entry.c_len_bytes);
                  entry.data->gid = x+i;
              }
            }}));
    }
    for (auto &it: threads)
    {
        it.join();
    }

    return output_penguin;
}
/*
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
        const int64_t c_len_in = conf.getCLen(table_index - 1);
        const int64_t c_len_out = conf.getCLen(table_index);

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

        new_entry.gid = new_gid_entry.getY();
        output_bucket_penguin->addEntry(new_entry);
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

*/

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

unsigned int countSetBits(int n)
{
    unsigned int count = 0;
    while (n) {
        n &= (n - 1);
        count++;
    }
    return count;
}

template <PlotConf conf, uint8_t table_index>
class EntryMatcher
{
public:

    static constexpr  uint32_t num_entries_per_buffer = 17;
    static constexpr uint32_t num_buckets_per_row = FwdYCEntry<conf, table_index-1>::num_bc_buckets_per_row;

    struct entryWithC
    {
        uint64_t bucket_id : ceillog2(num_buckets_per_row);
        uint128_t gid : FwdYCEntry<conf, table_index-1>::gid_len_bits;
        uint8_t c[FwdYCEntry<conf, table_index-1>::c_len_bytes];
    };

    struct entryWithoutC
    {
        uint64_t bucket_id : ceillog2(num_buckets_per_row);
        uint128_t gid : FwdYCEntry<conf, table_index-1>::gid_len_bits;
    };

    using entry_struct = typename std::conditional<(FwdYCEntry<conf, table_index-1>::c_len_bytes > 0), entryWithC, entryWithoutC>::type;

    array<array<entry_struct, num_entries_per_buffer>, kBC> match_staging;
    array<bitset<num_entries_per_buffer>, kBC> match_staging_used_map;

    void loadRow(vector<Penguin<FwdYCEntry<conf, table_index-1>, true>*> penguins, uint64_t row_id)
    {
        for (auto & it : match_staging_used_map)
        {
            it = 0;
        }
        for (auto & penguin : penguins)
        {
            const uint64_t num_entries = penguin->getCountInRow(row_id);
            for (uint64_t entry_id = 0; entry_id < num_entries; entry_id++)
            {
                auto entry = penguin->readEntry(row_id, entry_id);
                entry_struct tmp_entry;
                tmp_entry.bucket_id = entry.data->packed_y / kBC;
                tmp_entry.gid = entry.data->gid;
                if constexpr(FwdYCEntry<conf, table_index - 1>::c_len_bytes > 0) {
                    memcpy(tmp_entry.c, entry.data->c, FwdYCEntry<conf, table_index - 1>::c_len_bytes);
                }
                uint32_t bc = entry.data->packed_y % kBC;
                // Find destination
                uint32_t pos = (tmp_entry.bucket_id * num_entries_per_buffer) /
                               FwdYCEntry<conf, table_index - 1>::num_bc_buckets_per_row;
                uint8_t attempt_ct = 0;
                while (match_staging_used_map[bc][pos]) {
                    // Store lesser of the two values
                    if (tmp_entry.bucket_id < match_staging[bc][pos].bucket_id) {
                        swap(tmp_entry, match_staging[bc][pos]);
                    }
                    pos = (pos + 1) % num_entries_per_buffer;
                    attempt_ct++;
                    if (attempt_ct == num_entries_per_buffer) {
                        throw runtime_error("Filled up a match staging area!");
                    }
                }
                // Store new value
                match_staging[bc][pos].bucket_id = tmp_entry.bucket_id;
                match_staging[bc][pos].gid = tmp_entry.gid;
                if constexpr(FwdYCEntry<conf, table_index - 1>::c_len_bytes > 0) {
                    memcpy(match_staging[bc][pos].c, tmp_entry.c,
                           FwdYCEntry<conf, table_index - 1>::c_len_bytes);
                }
                // Mark used
                match_staging_used_map[bc][pos] = true;
            }
        }
    }
};

template <PlotConf conf>
template <uint8_t table_index>
Penguin<FwdYCEntry<conf, table_index>, true>*  Plotter<conf>::Context::createTable(
        vector<Penguin<FwdYCEntry<conf, table_index-1>, true>*> prev_penguins,
        atomic<uint64_t>& num_matches_out) {
    PinToCpuid(cpu_ids[0]);
    auto output_bucket_penguin = new Penguin<FwdYCEntry<conf, table_index>, true>();
    auto output_gid_penguin = new Penguin<FwdGIDEntry<conf, table_index>, true>();
    this->forward_pass_gid_penguins[table_index] = (void *) output_gid_penguin;

    vector<EntryMatcher<conf, table_index>> entry_matchers(3);
    atomic<int64_t> latest_active_matcher_buffer = -1;
    atomic<uint64_t> context_coordinator = 0;
    uint64_t last_matcher_buffer = (1ULL << 63);

    vector<thread> threads;

    for (auto &cpu_id : cpu_ids) {
        threads.push_back(thread([this,
                                         cpu_id,
                                         &entry_matchers,
                                         &last_matcher_buffer,
                                         &latest_active_matcher_buffer,
                                         &context_coordinator,
                                         &output_bucket_penguin,
                                         &output_gid_penguin,
                                         &num_matches_out] {

            PinToCpuid(cpu_id);
            uint64_t num_matches_in_thread = 0;

            struct MatchData {
                uint32_t bc;
                typename EntryMatcher<conf, table_index>::entry_struct *entry;
            };
            struct {
                bool operator()(const MatchData lhs, const MatchData &rhs) const {
                    return lhs.entry->gid < rhs.entry->gid;
                }
            } matchCompare;
            vector<MatchData> matches;
            matches.reserve(conf.max_gid_stub_val);

            while (true) {
                uint32_t batch_num = this->plotter->coordinator++;
                uint16_t b_id = batch_num % kB;
                uint32_t matcher_buffer_needed = batch_num / kB;

                // Wait for the necessary buffer to be available
                {
                    uint32_t matcher_buffer_avail = latest_active_matcher_buffer.load();
                    if (matcher_buffer_avail < matcher_buffer_needed) {
                        latest_active_matcher_buffer.wait(matcher_buffer_avail);
                    }
                    if (last_matcher_buffer < matcher_buffer_needed)
                    {
                        break;
                    }
                    // TODO: handle finish
                }
                EntryMatcher<conf, table_index> *matcher = &entry_matchers[matcher_buffer_needed %
                                                                          entry_matchers.size()];

                // Staging area is filled, now search for matches
                for (uint16_t c_id = 0; c_id < kC; c_id++) {
                    uint16_t bc = b_id * kC + c_id;

                    // bc targets for right and left match candidates
                    uint16_t right_targets[2][64];
                    uint16_t left_targets[2][64];
                    for (uint8_t parity = 0; parity < 2; parity++) {
                        for (uint16_t m = 0; m < kExtraBitsPow; m++) {
                            right_targets[parity][m] =
                                    ((b_id + m) % kB) * kC + ((bc + (2 * m + parity) * (2 * m + parity)) % kC);
                            left_targets[parity][m] =
                                    ((b_id - m) % kB) * kC + ((bc - (2 * m + parity) * (2 * m + parity)) % kC);
                        }
                    }

                    const uint32_t num_entries_per_buffer = EntryMatcher<conf, table_index>::num_entries_per_buffer;

                    for (uint8_t cpos = 0; cpos < num_entries_per_buffer; cpos++) {
                        if (!matcher->match_staging_used_map[bc][cpos]) {
                            continue;
                        }
                        // Scan for left and right matches of this value, only using those with lower GID then we have
                        const uint32_t our_bucket = matcher->match_staging[bc][cpos].bucket_id;
                        const bool parity = our_bucket % 2;
                        const uint128_t our_gid = matcher->match_staging[bc][cpos].gid;
                        matches.clear();
                        {
                            // looking for entries right of us
                            uint32_t r_bucket = our_bucket + 1;
                            uint32_t start_pos = (r_bucket * num_entries_per_buffer) /
                                                 FwdYCEntry<conf, table_index - 1>::num_bc_buckets_per_row;
                            for (auto &r_bc : right_targets[parity]) {
                                uint32_t pos = start_pos;
                                uint8_t attempt_ct = 0;
                                while (matcher->match_staging_used_map[r_bc][pos]) {
                                    if (matcher->match_staging[r_bc][pos].bucket_id == r_bucket) {
                                        if (matcher->match_staging[r_bc][pos].gid < our_gid) {
                                            MatchData match;
                                            match.entry = &(matcher->match_staging[r_bc][pos]);
                                            match.bc = r_bc;
                                            matches.push_back(match);
                                        }
                                    }
                                    pos++;
                                    attempt_ct++;
                                    if (attempt_ct == num_entries_per_buffer) {
                                        break;
                                    }
                                }
                            }
                        }
                        {
                            // looking for entries left of us
                            uint32_t l_bucket = our_bucket - 1;
                            uint32_t start_pos = (l_bucket * num_entries_per_buffer) /
                                                 FwdYCEntry<conf, table_index - 1>::num_bc_buckets_per_row;
                            for (auto &l_bc : left_targets[parity]) {
                                uint32_t pos = start_pos;
                                uint8_t attempt_ct = 0;
                                while (matcher->match_staging_used_map[l_bc][pos]) {
                                    if (matcher->match_staging[l_bc][pos].bucket_id == l_bucket) {
                                        if (matcher->match_staging[l_bc][pos].gid < our_gid) {
                                            MatchData match;
                                            match.entry = &(matcher->match_staging[l_bc][pos]);
                                            match.bc = l_bc;
                                            matches.push_back(match);
                                        }
                                    }
                                    pos++;
                                    attempt_ct++;
                                    if (attempt_ct == num_entries_per_buffer) {
                                        break;
                                    }
                                }
                            }
                        }
                        num_matches_in_thread += matches.size();
                        if (matches.size() > conf.max_gid_stub_val) {
                            throw new runtime_error("Too many matches for selected max_gid_stub_val");
                        }
                        // Sort matches
                        sort(matches.begin(), matches.end(), matchCompare);
                        // Process matches
                        uint128_t new_gid = our_gid * conf.max_gid_stub_val;
                        for (auto &match : matches) {
                            bool flip = match.entry->bucket_id > our_bucket;
                            uint32_t left_y = flip ? bc + our_bucket * kBC : match.bc + match.entry->bucket_id * kBC;
                            //uint32_t right_y = flip ? match.bc + match.entry->bucket_id*kBC : bc + our_bucket*kBC;
                            uint8_t *left_c = flip ? matcher->match_staging[bc][cpos].c : match.entry->c;
                            uint8_t *right_c = flip ? match.entry->c : matcher->match_staging[bc][cpos].c;

                            // Setup gid entry first
                            auto new_gid_entry = output_gid_penguin->newEntry(new_gid);
                            new_gid_entry.data->right_gid = match.entry->gid;

                            // Now setup bucket entry for next iteration
                            const int64_t c_len_bits_in = conf.getCLen(table_index - 1);
                            const int64_t c_len_bits_out = conf.getCLen(table_index);
                            const uint64_t c_len_bytes_in = FwdYCEntry<conf, table_index - 1>::c_len_bytes;
                            const uint64_t c_len_bytes_out = FwdYCEntry<conf, table_index>::c_len_bytes;

                            uint8_t input_data[64];
                            // Zero out so we can or in easily
                            memset(input_data, 0, 64);

                            // Copy in R
                            bitCopy<conf.K + kExtraBits + c_len_bits_in, 64,
                                    c_len_bytes_in * 8 - c_len_bits_in, FwdYCEntry<conf, table_index - 1>::c_len_bytes>(
                                    input_data, right_c);

                            // Copy in L
                            bitCopy<conf.K + kExtraBits, 64, c_len_bytes_in * 8 - c_len_bits_in, c_len_bytes_in>(
                                    input_data, left_c);

                            // Copy in F1
                            {
                                uint64_t f = bswap_64(left_y);
                                bitCopy<0, 64,
                                        8 * 8 - (conf.K + kExtraBits), 8>(input_data, (uint8_t *) &f);
                            }

                            blake3_hasher hasher;
                            blake3_hasher_init(&hasher);
                            blake3_hasher_update(&hasher,
                                                 input_data, (conf.K + kExtraBits + c_len_bits_in * 2 + 7) / 8);

                            uint8_t hash_bytes[32];
                            blake3_hasher_finalize(&hasher, hash_bytes, sizeof(hash_bytes));

                            uint64_t new_y = 0;
                            if constexpr (table_index < 5) {
                                new_y = bswap_64(*(uint64_t *) hash_bytes) >> (64 - (conf.K + kExtraBits));
                            } else {
                                new_y = bswap_64(*(uint64_t *) hash_bytes) >> (64 - (conf.K));
                            }

                            auto new_entry = output_bucket_penguin->newEntry(new_y);

                            if constexpr(c_len_bytes_out > 0) {
                                memset(new_entry.data->c, 0, c_len_bytes_out);
                                if constexpr (table_index < 2) {
                                    // Copy in R
                                    bitCopy<c_len_bytes_out * 8 - c_len_bits_in, c_len_bytes_out,
                                            c_len_bytes_in * 8 - c_len_bits_in, c_len_bytes_in>(new_entry.data->c,
                                                                                                right_c);

                                    // Copy in L
                                    bitCopy<c_len_bytes_out * 8 - c_len_bits_in * 2, c_len_bytes_out,
                                            c_len_bytes_in * 8 - c_len_bits_in, c_len_bytes_in>(new_entry.data->c,
                                                                                                left_c);
                                } else if constexpr (table_index < 5) {
                                    // Copy from hash result
                                    bitCopy<c_len_bytes_out * 8 - c_len_bits_out, c_len_bytes_out,
                                            conf.K + kExtraBits, 32>(new_entry.data->c, hash_bytes);
                                    assert(*new_entry.data->c < (1ULL << (c_len_bits_out - (c_len_bytes_out - 1) * 8)));
                                }
                            }
                            new_entry.data->gid = new_gid;

                            new_gid++;
                        }
                    }
                }
            }
            num_matches_out += num_matches_in_thread;
        }));
    }

    // Fill the entry_matchers datastructures ahead of the main matching process
    uint32_t row_id;
    uint32_t matcher_id = 0;
    while ((row_id = plotter->coordinator++) < FwdYCEntry<conf, table_index-1>::num_rows)
    {
        entry_matchers[matcher_id%entry_matchers.size()].loadRow(prev_penguins, row_id);
        latest_active_matcher_buffer = matcher_id;
        latest_active_matcher_buffer.notify_all();

        // Wait for local coordinator to tick over to the latest active buffer
        uint64_t latest_coordinator_val;
        while ((latest_coordinator_val = context_coordinator.load()) < kB*matcher_id)
        {
            context_coordinator.wait(latest_coordinator_val);
        }

        matcher_id++;
    }
    last_matcher_buffer = latest_active_matcher_buffer;
    latest_active_matcher_buffer.notify_all();


    for (auto &it: threads)
    {
        it.join();
    }

    return output_bucket_penguin;
}

template <PlotConf conf, uint8_t table_index>
vector<Penguin<FwdYCEntry<conf, table_index>, true>*> DoTable(
        Plotter<conf>* plotter,
        vector<Penguin<FwdYCEntry<conf, table_index-1>, true>*> input_yc_penguins)
{
    StatusUpdate::StartSeg("0." + to_string((uint32_t)table_index) + "0M");
    vector<thread> threads;
    plotter->coordinator = 0;
    atomic<uint64_t> num_matches = 0;
    vector<Penguin<FwdYCEntry<conf, table_index>, true>*> output_yc_penguin(plotter->contexts.size());
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

    vector<thread> threads;
    coordinator = 0;
    vector<Penguin<FwdYCEntry<conf, -1>, true>*> first_yc_penguin(contexts.size());
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
#include "plotconf.hpp"
#include "penguin_entries.hpp"
#include "penguin_unloaders.hpp"
