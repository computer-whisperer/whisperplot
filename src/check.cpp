//
// Created by christian on 6/5/21.
//

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
void Plotter<K, num_rows>::check_table1()
{
    // Check that all matching conditions from the first table actually work
    F1Calculator f1(K, id);

    for (auto& park : phase1_graph_parks[0]) {
        vector<uint128_t> line_points(park->size());
        park->readEntries(line_points);
        for (auto &line_point : line_points) {
            auto res = Encoding::LinePointToSquare(line_point);

            Bits x = f1.CalculateF(Bits(res.first, K));
            Bits y = f1.CalculateF(Bits(res.second, K));
            if (x.GetValue() > y.GetValue()) {
                swap(x, y);
            }
            assert_matching(x.GetValue(), y.GetValue());
        }
    }
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
    for (uint8_t table_index = 0; table_index < 7; table_index++)
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