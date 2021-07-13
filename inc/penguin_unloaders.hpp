//
// Created by christian on 6/27/21.
//

#ifndef WHISPERPLOT_PENGUIN_UNLOADERS_HPP
#define WHISPERPLOT_PENGUIN_UNLOADERS_HPP

#include "plotconf.hpp"
#include "explicit_templates.hpp"
#include "status_update.hpp"
#include "bitcopy.hpp"
#include "thread_mgr.hpp"
#include "pos_constants.hpp"
#include "plotter.hpp"
#include "encoding.hpp"
#include "calculate_bucket.hpp"
#include "buffer.hpp"
#include "penguin.hpp"
#include <ctime>
#include <cstring>
#include <atomic>
#include <vector>

template<PlotConf conf, class entry_type>
class BucketPenguinUnloader
{
    // These are BC buckets from the chia algo

    static constexpr uint32_t max_entries_per_bc_bucket = 350;

    std::vector<std::vector<entry_type>> temp_buckets;

    int64_t latest_bucket_started_load = -1;
    int64_t latest_bucket_finished_load = -1;
    int64_t latest_row_loaded = -1;

    const std::vector<Penguin<entry_type, conf.interlace_factor>*> penguins;

public:
    static constexpr uint32_t bc_bucket_num = (1ULL << (conf.K + kExtraBits)) / kBC;
    static constexpr uint32_t temp_buckets_needed = (entry_type::getMaxY()/conf.num_rows) / kBC + 12;

    explicit BucketPenguinUnloader(const std::vector<Penguin<entry_type, conf.interlace_factor>*> penguins_in)
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

        if (entry_type::getMaxY() <= last_value_contained)
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

    std::vector<entry_type>& getBucketEntries(uint32_t bucket_id)
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
            uint32_t first_row_needed = entry_type::getRowFromY(first_value_needed);
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

template<PlotConf conf, class entry_type, uint64_t num_entries_per_bucket>
class SortedParkedPenguinUnloader
{

    std::vector<std::vector<entry_type>> temp_buckets;

    std::vector<entry_type> temp_row;

    int64_t latest_bucket_started_load = -1;
    int64_t latest_bucket_finished_load = -1;
    int64_t latest_row_loaded = -1;

    uint64_t latest_row_start_pos_scanned = 0;
    uint64_t latest_row_start_pos = 0;

    const std::vector<Penguin<entry_type, conf.interlace_factor>*> penguins;

public:
    static constexpr uint32_t temp_buckets_needed = (entry_type::max_entries_per_row/num_entries_per_bucket)+1;

    explicit SortedParkedPenguinUnloader(const std::vector<Penguin<entry_type, conf.interlace_factor>*> penguins_in)
            : penguins(penguins_in)
    {
        temp_row.reserve(entry_type::max_entries_per_row);
        temp_buckets.resize(temp_buckets_needed);
        for (auto &i : temp_buckets)
        {
            i.reserve(num_entries_per_bucket);
        }
    }

    uint64_t getRowStartPos(uint32_t row_id)
    {
        assert(row_id >= latest_row_start_pos_scanned);
        while (latest_row_start_pos_scanned < row_id)
        {
            latest_row_start_pos_scanned++;
            for (auto & penguin : penguins)
            {
                latest_row_start_pos += penguin.getCountInRow(latest_row_start_pos_scanned);
            }
        }
        return latest_row_start_pos;
    }

    uint64_t getRowEndPos(uint32_t row_id)
    {
        uint64_t end_pos = getRowStartPos(row_id);
        for (auto & penguin : penguins)
        {
            end_pos += penguin.getCountInRow(row_id);
        }
        return end_pos;
    }

    void loadRow(uint32_t row_id)
    {
        uint64_t first_pos_contained = getRowStartPos(row_id);
        uint64_t first_bucket_needed = first_pos_contained/num_entries_per_bucket;

        uint64_t last_pos_contained = getRowStartPos(row_id)-1;
        for (auto & penguin : penguins)
        {
            last_pos_contained += penguin->getCountInRow(row_id);
        }
        uint64_t last_bucket_needed = last_pos_contained/num_entries_per_bucket;

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

        temp_row.clear();

        // Load in the row
        for (auto& penguin : penguins)
        {
            for (uint32_t entry_id = 0; entry_id < penguin->getCountInRow(row_id); entry_id++)
            {
                temp_row.push_back(penguin->readEntry(row_id, entry_id));
            }
        }

        // Sort the temp row
        struct {
            bool operator()(entry_type a, entry_type b) const { return a.getY() < b.getY();}
        } customLess;
        sort(temp_row.begin(), temp_row.end(), customLess);

        // Dump into temp buckets
        uint32_t current_pos = first_pos_contained;
        for (auto& entry : temp_row)
        {
            current_pos++;
            temp_buckets[current_pos/num_entries_per_bucket].push_back(entry);
        }

        latest_bucket_finished_load = (last_pos_contained+1)/num_entries_per_bucket - 1;
        latest_row_loaded = row_id;

        if (row_id >= conf.num_rows)
        {
            latest_bucket_finished_load = last_bucket_needed;
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

    std::vector<entry_type>& getBucketEntries(uint32_t bucket_id)
    {
        while (latest_bucket_finished_load < (int32_t)(bucket_id)) {
            // Load next row
            uint32_t row_id = latest_row_loaded + 1;

            if (row_id >= conf.num_rows) {
                latest_bucket_finished_load = bucket_id + 1;
                break;
            }

            uint64_t first_pos_needed = bucket_id * num_entries_per_bucket;
            while (getRowEndPos(row_id) < first_pos_needed) {
                row_id++;
            }

            loadRow(row_id);
        }
        return temp_buckets[bucket_id%temp_buckets_needed];
    }
};

#endif //WHISPERPLOT_PENGUIN_UNLOADERS_HPP
