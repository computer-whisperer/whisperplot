//
// Created by christian on 6/5/21.
//

#include <atomic>
#include <thread>
#include <vector>
#include "plotter.hpp"

using namespace std;

template <PlotConf conf>
int32_t Plotter<conf>::find_proofs(uint128_t challenge_in)
{
    F1Calculator f1(conf.K, id.data());
    uint64_t challenge = challenge_in%(1ULL<<conf.K);
    vector<uint64_t> gids;

    uint64_t row_id = FwdYCEntry<conf, 5>::getRowFromY(challenge);
    for (auto & context : contexts)
    {
        for (uint64_t i = 0; i < context->forward_pass_final_yc_penguin->getCountInRow(row_id); i++)
        {
            auto entry = context->forward_pass_final_yc_penguin->readEntry(row_id, i);
            if (entry.getY() == challenge)
            {
                gids.push_back(entry.data->gid);
            }
        }
    }

    uint64_t num_good = 0;
    uint64_t num_fail_matches = 0;
    uint64_t num_fail_value = 0;

    for (auto & gid : gids)
    {
        uint64_t y = 0;
        uint128_t c = 0;
        auto ret = check_match_and_return_values<5>(gid, y, c);

        if (y == challenge)
        {
            num_good++;
        }
        else if (ret)
        {
            num_fail_matches++;
        }
        else
        {
            num_fail_value++;
        }
    }

    return (int32_t)num_good;
}

template <PlotConf conf>
template <uint8_t table_index>
uint8_t Plotter<conf>::check_match_and_return_values(uint128_t gid, uint64_t& y, uint128_t& c)
{

    uint64_t yl = 0;
    uint64_t yr = 0;
    uint128_t cl = 0;
    uint128_t cr = 0;

    vector<FwdGIDEntry<conf, table_index>> candidate_entries;

    uint64_t row_id = FwdGIDEntry<conf, table_index>::getRowFromY(gid);
    for (auto & context : contexts)
    {
        auto penguin = (Penguin<FwdGIDEntry<conf, table_index>, true>*) context->forward_pass_gid_penguins[table_index];
        for (uint64_t i = 0; i < penguin->getCountInRow(row_id); i++)
        {
            auto entry = penguin->readEntry(row_id, i);
            if (entry.getY() == gid)
            {
                candidate_entries.push_back(entry);
            }
        }
    }

    if (candidate_entries.size() > 1)
    {
        throw runtime_error("Found 2 entries with same gid!");
    }

    if (candidate_entries.empty())
    {
        return 1;
    }

    auto entry = candidate_entries[0];

    uint8_t l_ret = 0;
    uint8_t r_ret = 0;
    switch(table_index)
    {
        case 5:
            l_ret = check_match_and_return_values<4>(entry.getY()/conf.max_gid_stub_val, ref(yl), ref(cl));
            r_ret = check_match_and_return_values<4>(entry.data->right_gid, ref(yr), ref(cr));
            break;
        case 4:
            l_ret = check_match_and_return_values<3>(entry.getY()/conf.max_gid_stub_val, ref(yl), ref(cl));
            r_ret = check_match_and_return_values<3>(entry.data->right_gid, ref(yr), ref(cr));
            break;
        case 3:
            l_ret = check_match_and_return_values<2>(entry.getY()/conf.max_gid_stub_val, ref(yl), ref(cl));
            r_ret = check_match_and_return_values<2>(entry.data->right_gid, ref(yr), ref(cr));
            break;
        case 2:
            l_ret = check_match_and_return_values<1>(entry.getY()/conf.max_gid_stub_val, ref(yl), ref(cl));
            r_ret = check_match_and_return_values<1>(entry.data->right_gid, ref(yr), ref(cr));
            break;
        case 1:
            l_ret = check_match_and_return_values<0>(entry.getY()/conf.max_gid_stub_val, ref(yl), ref(cl));
            r_ret = check_match_and_return_values<0>(entry.data->right_gid, ref(yr), ref(cr));
            break;
        case 0:
            F1Calculator f1(conf.K, id.data());
            cl = entry.getY()/conf.max_gid_stub_val;
            cr = entry.data->right_gid;
            yl = f1.CalculateF(Bits(cl, conf.K)).GetValue();
            yr = f1.CalculateF(Bits(cr, conf.K)).GetValue();
            break;
    }
    if (l_ret)
        return l_ret;
    if (r_ret)
        return r_ret;

    FxCalculator fx(conf.K, table_index+2);
    if (yl > yr)
    {
        swap(yl, yr);
        swap(cl, cr);
    }

    int64_t bucket_id_lout = yl/kBC;
    int64_t bucket_id_rout = yr/kBC;
    if (bucket_id_lout + 1 != bucket_id_rout)
    {
        return 1;
    }
    int64_t b_id_l = (yl%kBC)/kC;
    int64_t c_id_l = (yl%kBC)%kC;
    int64_t b_id_r = (yr%kBC)/kC;
    int64_t c_id_r = (yr%kBC)%kC;
    int64_t m = (kB + b_id_r - b_id_l)%kB;
    if ((uint64_t)m >= (1UL<<kExtraBits))
    {
        return 1;
    }
    int64_t a = c_id_r - c_id_l;
    int64_t b = 2*m + (bucket_id_lout%2);
    int64_t c2 = (b*b)%kC;
    int64_t d = (kC + a)%kC;
    if (c2 != d)
    {
        return 1;
    }

    uint64_t c_len = kVectorLens[table_index + 2] * conf.K;
    auto out = fx.CalculateBucket(
            Bits(yl, conf.K+kExtraBits),
            Bits(cl, c_len),
            Bits(cr, c_len));
    if (table_index < 5) {
        y = out.first.Slice(0, conf.K+kExtraBits).GetValue();
        uint8_t buff[32];
        memset(buff, 0, sizeof(buff));
        out.second.ToBytes(buff);
        c = Util::SliceInt128FromBytes(
                buff, 0, kVectorLens[table_index + 3] * conf.K);
    }
    else
    {
        y = out.first.Slice(0, conf.K).GetValue();
    }

    return 0;
}


template <PlotConf conf>
void Plotter<conf>::checkFullPlot()
{
    cout << "Verifying all proofs" << endl;
    uint64_t part_start_seconds = time(nullptr);

    atomic<uint64_t> proofs_found = 0;
    atomic<uint64_t> proofs_verified = 0;
    atomic<uint64_t> proofs_failed_matching = 0;
    atomic<uint64_t> proofs_failed_value = 0;

    vector<thread> threads;
    coordinator = 0;
    for (uint32_t i = 0; i < num_threads; i++)
    {
        threads.push_back(thread( [this,
                                  &proofs_found,
                                  &proofs_verified,
                                  &proofs_failed_matching,
                                  &proofs_failed_value]
        {
            vector<uint128_t> line_points;
            while(true)
            {
                uint64_t row_id = this->coordinator++;

                if (row_id >= conf.num_rows)
                  break;

                for (auto & context : contexts)
                {
                    for (uint64_t i = 0; i < context->forward_pass_final_yc_penguin->getCountInRow(row_id); i++)
                    {
                        auto entry = context->forward_pass_final_yc_penguin->readEntry(row_id, i);
                        proofs_found++;

                        uint64_t gid = entry.data->gid;
                        uint64_t y = entry.getY();

                        uint64_t y_calculated;
                        uint128_t c;
                        uint8_t ret = check_match_and_return_values<5>(gid, ref(y_calculated), ref(c));
                        if (ret)
                        {
                            proofs_failed_matching++;
                        }
                        else if (y != y_calculated)
                        {
                            proofs_failed_value++;
                        }
                        else
                        {
                            proofs_verified++;
                        }
                    }
                }
            }

            }));
    }
    for (auto &it: threads)
    {
        it.join();
    }

    cout << "Found " << proofs_found << " proofs in total" << endl;
    cout << "Found " << proofs_verified << " proofs that are good" << endl;
    cout << "Found " << proofs_failed_matching << " proofs that do not meet matching conditions" << endl;
    cout << "Found " << proofs_failed_value << " proofs that meet all matching conditions, but do not match provided f7" << endl;
    cout << "Verification finished in " << time(nullptr) - part_start_seconds << "s" << endl;
}


/*
template <PlotConf conf>
void Plotter<conf>::check_full_table(uint8_t table_index)
{
    cout << "Verifying all entries in table " << static_cast<uint32_t>(table_index) << " and below" << endl;
    uint64_t part_start_seconds = time(nullptr);


    atomic<uint64_t> entries_verified = 0;
    atomic<uint64_t> entries_failed_matching = 0;

    vector<thread> threads;
    std::atomic<uint64_t> coordinator = 0;
    for (auto & cpu_id : cpu_ids)
    {
        threads.push_back(thread( [this,
                                          &coordinator,
                                          table_index,
                                          &entries_verified,
                                          &entries_failed_matching]
          {
              uint64_t num_entries = (*(phase1_graph_parks[table_index].end()-1))->start_pos + (*(phase1_graph_parks[table_index].end()-1))->size();
              while(true) {
                  uint64_t i = (coordinator)++;
                  if (i >= num_entries)
                  {
                      break;
                  }
                  uint64_t y_calculated;
                  uint128_t c;
                  uint8_t ret = check_match_and_return_values(table_index, i, ref(y_calculated), ref(c));
                  if (ret) {
                      entries_failed_matching++;
                  } else {
                      entries_verified++;
                  }
              }
          }));
    }
    for (auto &it: threads)
    {
        it.join();
    }
    uint64_t num_entries = (*(phase1_graph_parks[table_index].end()-1))->start_pos + (*(phase1_graph_parks[table_index].end()-1))->size();
    cout << "Found " << num_entries << " entries in total" << endl;
    cout << "Found " << entries_verified << " entries that are good" << endl;
    cout << "Found " << entries_failed_matching << " entries that do not meet matching conditions" << endl;
    cout << "Verification finished in " << time(nullptr) - part_start_seconds << "s" << endl;
}
*/

template <PlotConf conf>
void Plotter<conf>::find_many_proofs(uint32_t n)
{
    uint32_t successes = 0;
    uint32_t failures = 0;
    for (uint32_t i = 0; i < n; i++)
    {
        auto r = find_proofs(i);
        if (r < 0)
        {
            failures++;
        }
        else
        {
            successes += r;
        }
    }
    cout << "Found " << successes << "/" << n << " good Proofs" << endl;
    cout << "Found " << failures << "/" << n << " failed Proofs" << endl;
}

/*
template <PlotConf conf>
void Plotter<conf>::check_parks_integrity()
{
    vector<uint128_t> line_points_check;
    uint128_t prev_last_entry = 0;
    for (uint8_t table_index = 0; table_index < 6; table_index++)
    {
        prev_last_entry = 0;
        for (auto& park : phase1_graph_parks[table_index])
        {
            if (park->size() > 0)
            {
                line_points_check.resize(park->size());
                park->readEntries(line_points_check);
                assert(line_points_check[0] > prev_last_entry);
                prev_last_entry = line_points_check[park->size()-1];
            }
        }
    }
    prev_last_entry = 0;
    for (auto& park : phase1_final_parks)
    {
        if (park->size() > 0)
        {
            line_points_check.resize(park->size());
            park->readEntries(line_points_check);
            assert(line_points_check[0] > prev_last_entry);
            prev_last_entry = line_points_check[park->size()-1];
        }
    }
}
*/
#include "explicit_templates.hpp"