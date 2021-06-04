#include <vector>
#include <atomic>
#include <cstring>
#include <ctime>

#include "penguin.hpp"
#include "buffer.hpp"
#include "calculate_bucket.hpp"
#include "encoding.hpp"
#include "park.hpp"
#include "plotter.hpp"
#include "pos_constants.hpp"
#include "thread_mgr.hpp"

using namespace std;

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

template <uint8_t K>
void Plotter<K>::phase1ThreadA(
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
void Plotter<K>::phase1ThreadB(
        uint32_t cpu_id,
        std::atomic<uint64_t> * coordinator,
        map<uint32_t, vector<uint64_t>> *new_entry_positions,
        map<uint32_t, Penguin<YCPackedEntry<K, table_index - 1>> *> *prev_penguins,
        Penguin<YCPackedEntry<K, table_index>> * new_yc_penguin,
        Penguin<LinePointEntryUIDPackedEntry<K>> * new_line_point_penguin)
{
    PinToCpuid(cpu_id);
  //  F1Calculator f1(K, id);
    vector<vector<uint32_t>> right_map(kBC);
    for (auto& i : right_map) {
        i.clear();
    }
    std::vector<uint32_t> right_map_clean;
    FxCalculator fx(K, table_index+2);

    // These are BC buckets from the chia algo
    constexpr uint32_t bc_buckets_per_sort_row = YCPackedEntry<K, table_index - 1>::num_bc_buckets_per_sort_row;
    constexpr uint32_t bc_bucket_num = (1ULL << (K+kExtraBits))/kBC;
    constexpr uint64_t num_entries_per_bc_bucket = YCPackedEntry<K, table_index - 1>::sizes.num_entries_per_bc_bucket;
    vector<PackedArray<YCPackedEntry<K, table_index - 1>, YCPackedEntry<K, table_index - 1>::sizes.num_entries_per_bc_bucket>> temp_buckets(bc_buckets_per_sort_row*2);
    uint32_t latest_bc_bucket_loaded = 0;
    vector<vector<uint32_t>> entry_positions(bc_buckets_per_sort_row*2);
    for (auto & i : entry_positions)
    {
        i.reserve(num_entries_per_bc_bucket);
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

            if (latest_bc_bucket_loaded <= bucket_id+1) {
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
                        right_entry.sort_row = (bucket_id+1) / bc_buckets_per_sort_row;

                        // Calculate the next F and dump a new entry into the bucket index for the next table
                        uint64_t c_len = kVectorLens[table_index + 2] * K;
                        auto out = fx.CalculateBucket(
                                Bits(left_entry.getY(), K + kExtraBits),
                                Bits(left_entry.c, c_len),
                                Bits(right_entry.c, c_len));

                        YCPackedEntry<K, table_index> new_entry;

                        if (table_index < 5) {
                            new_entry.setY(out.first.GetValue());
                            uint8_t buff[32];
                            memset(buff, 0, sizeof(buff));
                            out.second.ToBytes(buff);
                            new_entry.c = Util::SliceInt128FromBytes(
                                    buff, 0, kVectorLens[table_index + 3] * K);
                        }
                        else
                        {
                            new_entry.setY(out.first.Slice(0, K).GetValue());
                        }
                        uint64_t new_entry_id = new_yc_penguin->addEntry(new_entry);
                        uint64_t sort_row = new_entry.sort_row;

                        // Emit a new line point for sorting into the phase1 output buffer
                        uint64_t x, y;
                        if (table_index == 0) // First table just gets the two x values
                        {
                            x = left_entry.c;
                            y = right_entry.c;
                            assert_matching(left_entry.getY(), right_entry.getY());
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
void Plotter<K>::phase1ThreadC(
        uint32_t cpu_id,
        std::atomic<uint64_t>* coordinator,
        map<uint32_t, vector<uint64_t>> *new_entry_positions,
        map<uint32_t, Penguin<LinePointEntryUIDPackedEntry<K>> *> line_point_bucket_indexes,
        vector<DeltaPark<line_point_delta_len_bits> *> *parks,
        Buffer* buffer)
{
    PinToCpuid(cpu_id);
   // F1Calculator f1(K, id);
    vector<pair<uint32_t, LinePointEntryUIDPackedEntry<K>>> entries;
    entries.reserve(LinePointEntryUIDPackedEntry<K>::max_entries_per_sort_row);
    vector<uint128_t> line_points;
    line_points.reserve(LinePointEntryUIDPackedEntry<K>::max_entries_per_sort_row);
    vector<uint128_t> line_points_check(LinePointEntryUIDPackedEntry<K>::max_entries_per_sort_row);
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
                auto entry = penguin->readEntry(row_id, i);
                entries.push_back(pair<uint32_t, LinePointEntryUIDPackedEntry<K>>(numa_node,entry));
            }
            num_entries += entries_in_numa;
            //penguin->popRow(row_id);
        }
        assert(num_entries < LinePointEntryUIDPackedEntry<K>::max_entries_per_sort_row);

// Sort the bucket
        struct {
            bool operator()(pair<uint32_t, LinePointEntryUIDPackedEntry<K>> a, pair<uint32_t, LinePointEntryUIDPackedEntry<K>> b) const { return a.second.getLinePoint() < b.second.getLinePoint(); }
        } customLess;
        sort(entries.begin(), entries.begin()+num_entries, customLess);
// We have a sorted entries list! Emit to output buffer and record new positions
        auto new_park = new DeltaPark<line_point_delta_len_bits>(num_entries);
        new_park->start_pos = total_entries_so_far;
        uint64_t num_bytes = new_park->GetSpaceNeeded();
        new_park->bind(buffer->data + buffer->GetInsertionOffset(num_bytes));
        line_points.clear();
        for (uint64_t i = 0; i < num_entries; i++)
        {
            line_points.push_back(entries[i].second.getLinePoint());
            assert(entries[i].second.entry_uid < (1ULL<<(K+2)));
            (*new_entry_positions)[entries[i].first][entries[i].second.entry_uid] = i + total_entries_so_far;
        }
        new_park->addEntries(line_points);
        (*parks)[row_id] = new_park;
/*
        // Check that we can recover the data
        new_park->readEntries(line_points_check.data());
        for (uint32_t i = 0; i < num_entries; i++)
        {
            assert(line_points[i] == line_points_check[i]);
        }

        if (table_index == 0)
        {
            for (uint32_t i = 0; i < num_entries; i++)
            {
                auto res = Encoding::LinePointToSquare(line_points[i]);

                Bits x = f1.CalculateF(Bits(res.first, K));
                Bits y = f1.CalculateF(Bits(res.second, K));
                if (x.GetValue() > y.GetValue())
                {
                    swap(x, y);
                }
                assert_matching(x.GetValue(), y.GetValue());
            }
        }
        */
    }
}

template <uint8_t K>
void Plotter<K>::phase1ThreadD(
        uint32_t cpu_id,
        std::atomic<uint64_t>* coordinator,
        map<uint32_t, vector<uint64_t>> *new_entry_positions,
        std::map<uint32_t, Penguin<YCPackedEntry<K, 5>>*> line_point_penguins,
        vector<DeltaPark<finaltable_y_delta_len_bits> *> *parks,
        Buffer* buffer)
{
    PinToCpuid(cpu_id);
    vector<YCPackedEntry<K, 5>> entries(YCPackedEntry<K, 5>::max_entries_per_sort_row);
    vector<uint128_t> line_points(YCPackedEntry<K, 5>::max_entries_per_sort_row);
    uint64_t total_entries_so_far = 0;
    while (true)
    {
        uint64_t row_id = coordinator->fetch_add(1);
        if (row_id >= YCPackedEntry<K, 5>::num_sort_rows)
            break;

        for (auto& [numa_node, penguin] : line_point_penguins) {
            total_entries_so_far += penguin->getCountInRow(row_id);
        }
// Sort all entries in this bucket
        entries.clear();
        uint64_t num_entries = 0;
        for (auto& [numa_node, penguin] : line_point_penguins)
        {
            uint64_t entries_in_numa = penguin->getCountInRow(row_id);
            for (uint64_t i = 0; i < entries_in_numa; i++)
            {
                uint64_t new_pos = (*new_entry_positions)[numa_node][penguin->getUniqueIdentifier(row_id, i)];
                auto e = penguin->readEntry(row_id, i);

                e.setY(e.getY()<<(K+1) | new_pos);
                entries.push_back(e);
            }
            num_entries += entries_in_numa;
            penguin->popRow(row_id);
        }
        //assert(num_entries < YCPackedEntry<K, 5>::max_entries_per_sort_row);

// Sort the bucket
        struct {
            bool operator()(YCPackedEntry<K, 5> a, YCPackedEntry<K, 5> b) const { return a.getY() < b.getY(); }
        } customLess;
        sort(entries.begin(), entries.begin()+num_entries, customLess);

// We have a sorted entries list! Emit to output buffer
        auto new_park = new DeltaPark<finaltable_y_delta_len_bits>(num_entries);
        new_park->start_pos = total_entries_so_far - num_entries;
        uint64_t num_bytes = new_park->GetSpaceNeeded();
        new_park->bind(buffer->data + buffer->GetInsertionOffset(num_bytes));
        for (uint32_t i = 0; i < num_entries; i++)
        {
            line_points[i] = entries[i].getY();
        }
        new_park->addEntries(line_points);
        (*parks)[row_id] = new_park;
    }
}

template <uint8_t K>
template <int8_t table_index>
map<uint32_t, Penguin<YCPackedEntry<K, table_index>>*> Plotter<K>::phase1DoTable(
        map<uint32_t, Penguin<YCPackedEntry<K, table_index - 1>>*> prev_bucket_indexes)
{
    map<uint32_t, Penguin<YCPackedEntry<K, table_index>>*> next_bucket_indexes;
    map<uint32_t, Penguin<LinePointEntryUIDPackedEntry<K>>*> new_line_point_bucket_indexes;
    for (auto& numa_node : GetNUMANodesFromCpuIds(cpu_ids))
    {
        d_new_entry_positions[numa_node].resize((1ULL<<(K+2)));
        next_bucket_indexes[numa_node] = new Penguin<YCPackedEntry<K, table_index>>();
        next_bucket_indexes[numa_node]->pops_required_per_row = 2*YCPackedEntry<K, table_index>::num_bc_buckets_per_sort_row;
        new_line_point_bucket_indexes[numa_node] = new Penguin<LinePointEntryUIDPackedEntry<K>>();
    }

    cout << "Part B"<< (uint32_t)table_index;
    uint64_t start_seconds = time(nullptr);
    vector<thread> threads;
    std::atomic<uint64_t> coordinator = 0;
    for (auto & cpu_id : cpu_ids)
    {
        uint32_t numa_node = numa_node_of_cpu(cpu_id);
        threads.push_back(thread(
                Plotter<K>::phase1ThreadB<table_index>,
                cpu_id,
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
    DeltaPark<line_point_delta_len_bits> test_park(LinePointEntryUIDPackedEntry<K>::max_entries_per_sort_row);
    buffers.push_back(new Buffer(test_park.GetSpaceNeeded() * LinePointEntryUIDPackedEntry<K>::num_sort_rows));
    phase1_graph_parks.push_back(vector<DeltaPark<line_point_delta_len_bits>*>(LinePointEntryUIDPackedEntry<K>::num_sort_rows));

    coordinator = 0;
    cout << "Part C"<< (uint32_t)table_index;
    start_seconds = time(nullptr);
    threads.clear();
    for (auto & cpu_id : cpu_ids)
    {
        threads.push_back(thread(
                Plotter<K>::phase1ThreadC<table_index>,
                cpu_id,
                &coordinator,
                &d_new_entry_positions,
                new_line_point_bucket_indexes,
                &(phase1_graph_parks[table_index]),
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
void Plotter<K>::phase1()
{
    load_tables();

    uint64_t phase_start_seconds = time(nullptr);

    map<uint32_t, Penguin<YCPackedEntry<K, -1>>*> first_bucket_indexes;
    for (auto& numa_node : GetNUMANodesFromCpuIds(cpu_ids))
    {
        first_bucket_indexes[numa_node] = new Penguin<YCPackedEntry<K, -1>>();
        first_bucket_indexes[numa_node]->pops_required_per_row = 2*YCPackedEntry<K, -1>::num_bc_buckets_per_sort_row;
    }

    cout << "Part A";
    uint64_t part_start_seconds = time(nullptr);
    vector<thread> threads;
    std::atomic<uint64_t> coordinator = 0;
    for (auto & cpu_id : cpu_ids)
    {
        uint32_t numa_node = numa_node_of_cpu(cpu_id);
        threads.push_back(thread(
                Plotter<K>::phase1ThreadA,
                cpu_id,
                &coordinator,
                id,
                first_bucket_indexes[numa_node]));
    }

    for (auto &it: threads)
    {
        it.join();
    }
    cout << " (" << time(nullptr) - part_start_seconds << "s)" << endl;

    auto i0 = phase1DoTable<0>(first_bucket_indexes);
    auto i1 = phase1DoTable<1>(i0);
    auto i2 = phase1DoTable<2>(i1);
    auto i3 = phase1DoTable<3>(i2);
    auto i4 = phase1DoTable<4>(i3);
    auto i5 = phase1DoTable<5>(i4);

    DeltaPark<finaltable_y_delta_len_bits> test_park(YCPackedEntry<K, 5>::max_entries_per_sort_row);
    buffers.push_back(new Buffer(test_park.GetSpaceNeeded() * YCPackedEntry<K, 5>::num_sort_rows));
    phase1_final_parks.resize(YCPackedEntry<K, 5>::num_sort_rows);



    cout << "Part D";
    part_start_seconds = time(nullptr);
    threads.clear();
    coordinator = 0;
    for (auto & cpu_id : cpu_ids)
    {
        threads.push_back(thread(
                Plotter<K>::phase1ThreadD,
                cpu_id,
                &coordinator,
                &d_new_entry_positions,
                i5,
                &phase1_final_parks,
                buffers[6]));
    }
    for (auto &it: threads)
    {
        it.join();
    }
    cout << " (" << time(nullptr) - part_start_seconds << "s)" << endl;
    cout << "Phase 1 finished in " << time(nullptr) - phase_start_seconds << "s" << endl;
}


template <uint8_t K>
void Plotter<K>::check()
{
    // Check that all matching conditions from the first table actually work
	F1Calculator f1(K, id);
	/*
    for (auto& park : phase1_graph_parks[0])
    {
        vector<uint128_t> line_points(park->size());
        park->readEntries(line_points);
        for (auto& line_point : line_points)
        {
            auto res = Encoding::LinePointToSquare(line_point);

            Bits x = f1.CalculateF(Bits(res.first, K));
            Bits y = f1.CalculateF(Bits(res.second, K));
            if (x.GetValue() > y.GetValue())
            {
                swap(x, y);
            }
            assert_matching(x.GetValue(), y.GetValue());
        }
    }
    */
	// Test that we can draw good proofs from this
	uint64_t challenge = 0x4819;
	challenge = challenge%(1ULL<<K);
	uint64_t i;
	cout << "Looking for : " << challenge << endl;
    uint64_t pos;
    bool found = false;
    for (auto& park : phase1_final_parks)
    {
        vector<uint128_t> line_points(park->size());
        park->readEntries(line_points);
        for (auto& line_point : line_points)
        {
            pos = line_point&((1ULL << (K+1))-1);
            uint64_t y = line_point>>(K+1);
            if (y == challenge) {
                found = true;
                break;
            }
        }
	}
	if (found)
	{
		cout << "Got inputs: ";
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
			cout << x[j];
			if (j < 63)
			{
				cout << ", ";
			}
		}
		cout << endl;
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
				for (i = 0; i < input_collations.size(); i += 2)
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
		cout << "Result of tree: f7(...) = " << input_fs[0].Slice(0, K).GetValue() << endl;
	}
	else
	{
		cout << "Could not find " << challenge << " in table 7." << endl;
	}
}

#include "explicit_templates.hpp"