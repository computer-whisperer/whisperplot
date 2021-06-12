#ifndef PHASE1_HPP
#define PHASE1_HPP

#include <iostream>
#include <vector>
#include "penguin.hpp"
#include "buffer.hpp"
#include "park.hpp"
#include "plotter.hpp"
#include "packed_array.hpp"

uint8_t getK(std::string plot_fname);

template <uint8_t K, uint32_t num_rows>
class Plotter
{
public:

    static constexpr uint64_t GetMaxY(int8_t table_index)
    {
        if (table_index < 5)
            return 1ULL << (K+kExtraBits);
        else
            return 1ULL << K;
    }

    static constexpr uint128_t GetMaxLinePoint(int8_t table_index)
    {
        uint128_t max_y = (1ULL << K)*1.1 - 1;
        uint128_t max_x = (1ULL << K)*1.1;
        if (table_index == 0)
        {
            max_y = (1ULL << K)-1;
            max_x = 1ULL << K;
        }

        uint64_t a = max_x, b = max_x - 1;
        if (a % 2 == 0)
            a /= 2;
        else
            b /= 2;

        return (uint128_t)a * b + max_y;
    }

    static constexpr uint64_t GetMeanEntryCount(int8_t table_index)
    {
        return 1ULL << K;
    }

    static constexpr uint64_t GetCLen(int8_t table_index)
    {
        switch (table_index)
        {
            case -1:
                return K;
            case 0:
                return K*2;
            case 1:
            case 2:
                return K*4;
            case 3:
                return K*3;
            case 4:
                return K*2;
            default:
                return 0;
        }
    }

    template<int8_t table_index>
    class YCPackedEntry : public PackedEntry<num_rows, GetMaxY(table_index)+1, GetMeanEntryCount(table_index)>
    {
        using parent = PackedEntry<num_rows, GetMaxY(table_index)+1, GetMeanEntryCount(table_index)>;
    public:
        static constexpr uint64_t c_len_bits = GetCLen(table_index);

        struct packed_type {
            uint32_t y : parent::trimmed_y_len_bits;
            uint64_t c : c_len_bits;
        };

        static constexpr uint64_t len_bits = parent::trimmed_y_len_bits + c_len_bits;
        //static constexpr uint64_t len_bits = sizeof(packed_type)*8;
        uint128_t c;
        inline void pack(uint8_t * dest, uint64_t offset) override
        {
            assert(this->y < GetMaxY(table_index)+1);
            /*
            struct packed_type temp;
            temp.y = this->y;
            temp.c = c;
            memcpy(dest + offset/8, &temp, sizeof(temp));
             */

            /*

            auto packed_dat = (struct packed_type*)dest;
            packed_dat->y = this->y;
            packed_dat->c = c;*/
/*
            dest += offset/8;
            uint128_t temp = this->y | (c << parent::trimmed_y_len_bits);
            temp = bswap_128(temp);
            memcpy(dest, &temp, (len_bits+7)/8);
            */
            std::span<uint8_t> s{dest, (offset + len_bits + 7)/8};
            if (parent::trimmed_y_len_bits)
            {
                bitpacker::insert(s, offset, parent::trimmed_y_len_bits,  this->y);
            }
            if (c_len_bits)
            {
                assert(c < (((uint128_t)(1ULL)) << c_len_bits));
                bitpacker::insert(s, offset + parent::trimmed_y_len_bits, c_len_bits,  c);
            }

        }

        inline void unpack(uint8_t * src, uint64_t offset) override
        {
            /*
            struct packed_type temp;
            memcpy(&temp, src + offset/8, sizeof(temp));
            this->y = temp.y;
            c = temp.c;
             */
            /*
            src += offset/8;
            auto packed_dat = (struct packed_type*)src;
            this->y = packed_dat->y;
            c = packed_dat->c;
*/
            /*
            src += offset/8;
            uint128_t temp = 0;
            memcpy(&temp, src, (len_bits+7)/8);
            temp = bswap_128(temp);
            this->y = temp % parent::row_divisor;
            c = temp >> parent::trimmed_y_len_bits;
            // = this->y | (c << parent::trimmed_y_len_bits);*/

            std::span<uint8_t> s{src, (offset + len_bits + 7)/8};
            if (parent::trimmed_y_len_bits)
            {
                this->y = bitpacker::extract<uint128_t>(s, offset, parent::trimmed_y_len_bits);
            }
            if (c_len_bits)
            {
                c = bitpacker::extract<uint128_t>(s, offset + parent::trimmed_y_len_bits, c_len_bits);
            }

        }
    };

    template<int8_t table_index>
    class YCTempPackedEntry : public PackedEntry<GetMaxY(table_index)/kBC, GetMaxY(table_index)+1, GetMeanEntryCount(table_index)>
    {
        using parent = PackedEntry<GetMaxY(table_index)/kBC, GetMaxY(table_index)+1, GetMeanEntryCount(table_index)>;
    public:
        static constexpr uint64_t c_len_bits = GetCLen(table_index);
        static constexpr uint64_t len_bits = parent::trimmed_y_len_bits + c_len_bits;
        uint128_t c;
        inline void pack(uint8_t * dest, uint64_t offset) override
        {
            assert(this->y < GetMaxY(table_index)+1);
            std::span<uint8_t> s{dest, (offset + len_bits + 7)/8};
            if (parent::trimmed_y_len_bits)
            {
                bitpacker::insert(s, offset, parent::trimmed_y_len_bits,  this->y);
            }
            if (c_len_bits)
            {
                assert(c < (((uint128_t)(1ULL)) << c_len_bits));
                bitpacker::insert(s, offset + parent::trimmed_y_len_bits, c_len_bits,  c);
            }
        }

        inline void unpack(uint8_t * src, uint64_t offset) override
        {
            std::span<uint8_t> s{src, (offset + len_bits + 7)/8};
            if (parent::trimmed_y_len_bits)
            {
                this->y = bitpacker::extract<uint128_t>(s, offset, parent::trimmed_y_len_bits);
            }
            if (c_len_bits)
            {
                c = bitpacker::extract<uint128_t>(s, offset + parent::trimmed_y_len_bits, c_len_bits);
            }
        }
    };

    template<int8_t table_index>
    class LinePointUIDPackedEntry : public PackedEntry<num_rows, GetMaxLinePoint(table_index)+1, GetMeanEntryCount(table_index)*2>
    {
        using parent = PackedEntry<num_rows, GetMaxLinePoint(table_index)+1, GetMeanEntryCount(table_index)*2>;
    public:
        static constexpr uint64_t uid_len_bits = K+2;
        static constexpr uint64_t len_bits = parent::trimmed_y_len_bits + uid_len_bits;
        uint64_t uid;
        inline void pack(uint8_t * dest, uint64_t offset) override
        {
            assert(this->y < GetMaxLinePoint(table_index)+1);
            std::span<uint8_t> s{dest, (offset + len_bits + 7)/8};
            if (parent::trimmed_y_len_bits)
            {
                bitpacker::insert(s, offset, parent::trimmed_y_len_bits,  this->y);
            }
            if (uid_len_bits)
            {
                assert(uid < (1ULL << uid_len_bits));
                bitpacker::insert(s, offset + parent::trimmed_y_len_bits, uid_len_bits,  uid);
            }
        }

        inline void unpack(uint8_t * src, uint64_t offset) override
        {
            std::span<uint8_t> s{src, (offset + len_bits + 7)/8};
            if (parent::trimmed_y_len_bits)
            {
                this->y = bitpacker::extract<uint64_t>(s, offset, parent::trimmed_y_len_bits);
            }
            if (uid_len_bits)
            {
                uid = bitpacker::extract<uint64_t>(s, offset + parent::trimmed_y_len_bits, uid_len_bits);
            }
        }
    };

    static constexpr uint8_t line_point_delta_len_bits = K + 10;
    static constexpr uint8_t finaltable_y_delta_len_bits = K + 10;
    static constexpr uint64_t max_entries_per_graph_table = (1ULL << K)*1.1;

    std::vector<Buffer*> buffers;
    std::vector<std::vector<Park*>> phase1_graph_parks;
    std::vector<DeltaPark<finaltable_y_delta_len_bits>*> phase1_final_parks;
    std::string filename;

    inline Plotter(const uint8_t* id_in, uint8_t* memo_in, uint32_t memo_size_in, std::vector<uint32_t> cpu_ids_in, std::string filename_in)
    {
        memcpy(id, id_in, sizeof(id));
        cpu_ids = cpu_ids_in;
        filename = filename_in;
        memo = memo_in;
        memo_size = memo_size_in;
    }

    inline Plotter(std::vector<uint32_t> cpu_ids_in)
    {
        cpu_ids = cpu_ids_in;
    }

    void phase1();
    void phase2();
    void phase3();
    void phase4();
    void read(std::string plot_fname, bool decompress=true);

    int32_t find_proof(uint128_t challenge);
    void find_many_proofs(uint32_t n);
    void check_parks_integrity();
    static void assert_matching(uint64_t lout, uint64_t rout);
    void check_full_plot();

    void check_full_table(uint8_t table_index);

    using p2_final_positions_type = PackedArray<PackedEntry<1, 1ULL << (K+1), 1>, max_entries_per_graph_table>;
    using entries_used_type = BasePackedArray<BooleanPackedEntry, max_entries_per_graph_table, 8>;

    using phase1_new_positions_type = BasePackedArray<PackedEntry<1, (1ULL<<(K+1)), 1>, (1ULL<<(K+1)), 8>;

    template <int8_t table_index>
    using p1_buckets_done_type = BasePackedArray<BooleanPackedEntry, GetMaxY(table_index)/kBC + 1, 8>;


private:
    Buffer* output_buffer;
    uint8_t id[32];
    uint8_t* memo;
    uint32_t memo_size;
    std::vector<std::thread> phase2b_threads;
    std::vector<uint32_t> cpu_ids;
    std::map<uint32_t, phase1_new_positions_type*> phase1_new_entry_positions;
    std::vector<entries_used_type>* entries_used;
    std::vector<p2_final_positions_type>* phase2_final_positions;
    uint64_t final_table_begin_pointers[10];
    uint64_t pointer_table_offset;
    std::vector<std::vector<Park*>> final_parks;

    static void phase1ThreadA(
            uint32_t cpu_id,
            std::atomic<uint64_t>* coordinator,
            const uint8_t* id,
            Penguin<YCPackedEntry<-1>>* new_penguin
    );

    template <int8_t table_index>
    static void phase1ThreadB(
            uint32_t cpu_id,
            const uint8_t* id,
            std::atomic<uint64_t>* coordinator,
            p1_buckets_done_type<table_index-1> * bucket_left_done,
            p1_buckets_done_type<table_index-1> * bucket_right_done,
            std::map<uint32_t, phase1_new_positions_type*>* new_entry_positions,
            std::map<uint32_t, Penguin<YCPackedEntry<table_index - 1>>*>* prev_penguins,
            Penguin<YCPackedEntry<table_index>>* new_yc_penguin,
            Penguin<LinePointUIDPackedEntry<table_index>>* new_line_point_penguin);

    template <int8_t table_index>
    static void phase1ThreadC(
            uint32_t cpu_id,
            const uint8_t* id,
            std::atomic<uint64_t>* coordinator,
            std::map<uint32_t, phase1_new_positions_type*>* new_entry_positions,
            std::map<uint32_t, Penguin<LinePointUIDPackedEntry<table_index>>*> line_point_bucket_index,
            std::vector<Park*> *parks,
            Buffer* buffer
    );

    static void phase1ThreadD(
            uint32_t cpu_id,
            std::atomic<uint64_t>* coordinator,
            std::map<uint32_t, phase1_new_positions_type*>* new_entry_positions,
            std::map<uint32_t, Penguin<YCPackedEntry<5>>*> line_point_penguins,
            std::vector<DeltaPark<finaltable_y_delta_len_bits> *> *parks,
            Buffer* buffer
    );

    template <int8_t table_index>
    std::map<uint32_t, Penguin<YCPackedEntry<table_index>>*> phase1DoTable(
            std::map<uint32_t, Penguin<YCPackedEntry<table_index - 1>>*> prev_bucket_indexes);

    template <int8_t table_index>
    static void phase2ThreadA(
            uint32_t cpu_id,
            std::atomic<uint64_t>* coordinator,
            std::vector<Park*>* parks,
            entries_used_type* current_entries_used,
            entries_used_type* prev_entries_used
    );

    static void phase2ThreadB(
            entries_used_type* entries_used,
            p2_final_positions_type* final_positions);

    template <int8_t table_index>
    void phase2DoTable();

    template <int8_t table_index>
    static void phase3ThreadA(
            uint32_t cpu_id,
            std::atomic<uint64_t>* coordinator,
            std::vector<Park*>* temporary_parks,
            entries_used_type* entries_used,
            std::vector<p2_final_positions_type>* final_positions,
            Buffer* output_buffer,
            std::vector<std::vector<Park*>>* final_parks,
            uint64_t start_offset);

    template <int8_t table_index>
    void phase3DoTable();

    template <int8_t table_index>
    void readGraphTable(std::ifstream& file_stream, bool decompress);

    void readFinalTable(std::ifstream& file_stream);

    void checkFullPlotThread(
            std::atomic<uint64_t>* coordinator,
            std::atomic<uint64_t>* proofs_found_out,
            std::atomic<uint64_t>* proofs_verified_out,
            std::atomic<uint64_t>* proofs_failed_matching_out,
            std::atomic<uint64_t>* proofs_failed_value_out);

    uint8_t check_match_and_return_values(uint8_t table_index, uint64_t position, uint64_t& y, uint128_t& c);
};



#endif
