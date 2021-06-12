//
// Created by christian on 6/5/21.
//

#include <atomic>
#include <thread>
#include <vector>
#include "plotter.hpp"

using namespace std;

template <uint8_t K, uint32_t num_rows>
void Plotter<K, num_rows>::assert_matching(uint64_t lout, uint64_t rout)
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

template <uint8_t K, uint32_t num_rows>
int32_t Plotter<K, num_rows>::find_proof(uint128_t challenge_in)
{
    F1Calculator f1(K, id);
    uint64_t challenge = challenge_in%(1ULL<<K);
    //cout << "Looking for : " << challenge << endl;
    int32_t proof_num = 1;
    while (true)
    {
        uint64_t pos;
        uint8_t found = 0;
        for (auto& park : phase1_final_parks)
        {
            vector<uint128_t> line_points(park->size());
            park->readEntries(line_points);
            for (auto& line_point : line_points)
            {
                pos = line_point&((1ULL << (K+1))-1);
                uint64_t y = line_point>>(K+1);
                if (y == challenge) {
                    found++;
                    if (found == proof_num)
                    {
                        break;
                    }
                }
            }
            if (found == proof_num) {
                break;
            }
        }
        if (found == proof_num)
        {
            // cout << "Got inputs: ";
            vector<uint64_t> x(64);
            for (uint32_t j = 0; j < 64; j++)
            {
                uint64_t next_pos = pos;
                uint32_t l;
                for (l = 5; l < 6; l--)
                {
                    // Find next_pos
                    uint64_t entries_so_far = 0;
                    for (auto& park : phase1_graph_parks[l])
                    {
                        entries_so_far += park->size();
                        if (next_pos < entries_so_far)
                        {
                            // pos is in this park
                            vector<uint128_t> line_points(park->size());
                            park->readEntries(line_points);
                            uint128_t line_point = line_points[park->size() - (entries_so_far - next_pos)];
                            auto res = Encoding::LinePointToSquare(line_point);
                            uint8_t mask = 1ULL<<(l);
                            if (j & mask)
                            {
                                next_pos = res.first;
                            }
                            else
                            {
                                next_pos = res.second;
                            }
                            break;
                        }
                    }
                }
                x[j] = next_pos;
                /*
                cout << x[j];
                if (j < 63)
                {
                    cout << ", ";
                }*/
            }
            // cout << endl;
            // Check for matching conditions and re-combine
            vector<Bits> input_collations;
            vector<Bits> input_fs;
            for (uint32_t fi = 1; fi <= 7; fi++)
            {
                vector<Bits> output_fs;
                vector<Bits> output_collations;
                if (fi == 1)
                {
                    for (unsigned long i : x)
                    {
                        Bits b = Bits(i, K);
                        output_collations.push_back(b);
                        output_fs.push_back(f1.CalculateF(b));
                    }
                }
                else
                {
                    FxCalculator fx(K, fi);
                    for (uint32_t i = 0; i < input_collations.size(); i += 2)
                    {
                        // Swap inputs if needed
                        auto yl = input_fs[i];
                        auto yr = input_fs[i+1];
                        auto cl = input_collations[i];
                        auto cr = input_collations[i+1];
                        if (yl.GetValue() > yr.GetValue())
                        {
                            swap(yl, yr);
                            swap(cl, cr);
                        }
                        assert_matching(yl.GetValue(), yr.GetValue());
                        auto out = fx.CalculateBucket(yl, cl, cr);
                        output_fs.push_back(out.first);
                        output_collations.push_back(out.second);
                    }
                }
                input_collations = output_collations;
                input_fs = output_fs;
            }
            //cout << "Result of tree: f7(...) = " << input_fs[0].Slice(0, K).GetValue() << endl;
            if (input_fs[0].Slice(0, K).GetValue() == challenge)
            {
                proof_num++;
            }
            else
            {
                return -1;
            }
        }
        else
        {
            return proof_num-1;
        }
    }

}

template <uint8_t K, uint32_t num_rows>
uint8_t Plotter<K, num_rows>::check_match_and_return_values(uint8_t table_index, uint64_t position, uint64_t& y, uint128_t& c)
{

    uint64_t yl = 0;
    uint64_t yr = 0;
    uint128_t cl = 0;
    uint128_t cr = 0;

    bool found_it = false;
    uint128_t line_point = 0;

    // Find entry
    uint64_t entries_so_far = 0;
    for (auto& park : phase1_graph_parks[table_index])
    {
        entries_so_far += park->size();
        if (position < entries_so_far)
        {
            // pos is in this park
            vector<uint128_t> line_points(park->size());
            park->readEntries(line_points);
            line_point = line_points[park->size() - (entries_so_far - position)];
            found_it = true;
            break;
        }
    }

    if (!found_it)
    {
        return 1;
    }
    auto square = Encoding::LinePointToSquare(line_point);

    if(table_index >= 1)
    {
        uint8_t l_ret = check_match_and_return_values(table_index-1, square.first, ref(yl), ref(cl));
        if (l_ret)
        {
            return l_ret;
        }
        uint8_t r_ret = check_match_and_return_values(table_index-1, square.second, ref(yr), ref(cr));
        if (r_ret)
        {
            return r_ret;
        }
    }
    else
    {
        F1Calculator f1(K, id);
        cl = square.first;
        cr = square.second;
        yl = f1.CalculateF(Bits(cl, K)).GetValue();
        yr = f1.CalculateF(Bits(cr, K)).GetValue();
    }

    FxCalculator fx(K, table_index+2);
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

    uint64_t c_len = kVectorLens[table_index + 2] * K;
    auto out = fx.CalculateBucket(
            Bits(yl, K+kExtraBits),
            Bits(cl, c_len),
            Bits(cr, c_len));
    if (table_index < 5) {
        y = out.first.Slice(0, K+kExtraBits).GetValue();
        uint8_t buff[32];
        memset(buff, 0, sizeof(buff));
        out.second.ToBytes(buff);
        c = Util::SliceInt128FromBytes(
                buff, 0, kVectorLens[table_index + 3] * K);
    }
    else
    {
        y = out.first.Slice(0, K).GetValue();
    }

    return 0;
}

template <uint8_t K, uint32_t num_rows>
void Plotter<K, num_rows>::checkFullPlotThread(
        std::atomic<uint64_t>* coordinator,
        std::atomic<uint64_t>* proofs_found_out,
        std::atomic<uint64_t>* proofs_verified_out,
        std::atomic<uint64_t>* proofs_failed_matching_out,
        std::atomic<uint64_t>* proofs_failed_value_out)
{
    uint64_t proofs_found = 0;
    uint64_t proofs_verified = 0;
    uint64_t proofs_failed_matching = 0;
    uint64_t proofs_failed_value = 0;

    vector<uint128_t> line_points;
    while(true)
    {
        uint64_t park_id = (*coordinator)++;

        if (park_id >= phase1_final_parks.size())
            break;

        auto park = phase1_final_parks[park_id];

        if (park->size() > 0)
        {
            line_points.resize(park->size());
            park->readEntries(line_points);

            for (auto & line_point : line_points)
            {
                proofs_found++;

                uint64_t pos = line_point%(1ULL << (K+1));
                uint64_t y = line_point>>(K+1);

                uint64_t y_calculated;
                uint128_t c;
                uint8_t ret = check_match_and_return_values(5, pos, ref(y_calculated), ref(c));
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
    *proofs_verified_out += proofs_verified;
    *proofs_found_out += proofs_found;
    *proofs_failed_matching_out += proofs_failed_matching;
    *proofs_failed_value_out += proofs_failed_value;
}

template <uint8_t K, uint32_t num_rows>
void Plotter<K, num_rows>::check_full_plot()
{
    cout << "Verifying all proofs" << endl;
    uint64_t part_start_seconds = time(nullptr);

    atomic<uint64_t> proofs_found = 0;
    atomic<uint64_t> proofs_verified = 0;
    atomic<uint64_t> proofs_failed_matching = 0;
    atomic<uint64_t> proofs_failed_value = 0;

    vector<thread> threads;
    std::atomic<uint64_t> coordinator = 0;
    for (auto & cpu_id : cpu_ids)
    {
        threads.push_back(thread( [this,
                                   &coordinator,
                                  &proofs_found,
                                  &proofs_verified,
                                  &proofs_failed_matching,
                                  &proofs_failed_value]
                                  {this->checkFullPlotThread(
                                    &coordinator,
                                    &proofs_found,
                                    &proofs_verified,
                                    &proofs_failed_matching,
                                    &proofs_failed_value);}));
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


template <uint8_t K, uint32_t num_rows>
void Plotter<K, num_rows>::check_full_table(uint8_t table_index)
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

template <uint8_t K, uint32_t num_rows>
void Plotter<K, num_rows>::find_many_proofs(uint32_t n)
{
    uint32_t successes = 0;
    uint32_t failures = 0;
    for (uint32_t i = 0; i < n; i++)
    {
        auto r = find_proof(i);
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

template <uint8_t K, uint32_t num_rows>
void Plotter<K, num_rows>::check_parks_integrity()
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

#include "explicit_templates.hpp"