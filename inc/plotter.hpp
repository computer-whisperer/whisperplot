#ifndef PHASE1_HPP
#define PHASE1_HPP

#include <iostream>
#include <vector>
#include "penguin.hpp"
#include "buffer.hpp"
#include "park.hpp"
#include "plotter.hpp"
#include "packed_array.hpp"

struct PlotConf
{
    uint8_t K;
    uint32_t num_rows;
    uint32_t interlace_factor;
};

uint8_t getK(std::string plot_fname);

template <PlotConf conf>
class Plotter
{
public:

    static constexpr uint64_t max_entries_per_graph_table = (1ULL << conf.K)*1.1;

    static constexpr uint64_t GetMaxY(int8_t table_index)
    {
        if (table_index < 5)
            return 1ULL << (conf.K+kExtraBits);
        else
            return 1ULL << conf.K;
    }

    static constexpr uint128_t GetMaxLinePoint(int8_t table_index)
    {
        uint128_t max_y = (1ULL << conf.K)*1.1 - 1;
        uint128_t max_x = (1ULL << conf.K)*1.1;
        if (table_index == 0)
        {
            max_y = (1ULL << conf.K)-1;
            max_x = 1ULL << conf.K;
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
        return 1ULL << conf.K;
    }

    static constexpr uint64_t GetCLen(int8_t table_index)
    {
        switch (table_index)
        {
            case -1:
                return conf.K;
            case 0:
                return conf.K*2;
            case 1:
            case 2:
                return conf.K*4;
            case 3:
                return conf.K*3;
            case 4:
                return conf.K*2;
            default:
                return 0;
        }
    }

    template<int8_t table_index>
    class YCPackedEntry : public PackedEntry<conf.num_rows, GetMaxY(table_index)+1, GetMeanEntryCount(table_index)>
    {
        using parent = PackedEntry<conf.num_rows, GetMaxY(table_index)+1, GetMeanEntryCount(table_index)>;
    public:
        static constexpr uint64_t c_len_bits = GetCLen(table_index);

        struct packed_type {
            uint32_t y : parent::trimmed_y_len_bits;
            uint64_t c : c_len_bits;
        };

        static constexpr uint64_t packed_y_len_bytes = ((parent::trimmed_y_len_bits+31)/32)*4;
        static constexpr uint64_t packed_c_len_bytes = ((c_len_bits+31)/32)*4;

        static constexpr uint64_t len_bits = (packed_y_len_bytes + packed_c_len_bytes)*8;
        //static constexpr uint64_t len_bits = parent::trimmed_y_len_bits + c_len_bits;
        //static constexpr uint64_t len_bits = sizeof(packed_type)*8;
        uint128_t c;
        inline void pack(uint8_t * dest, uint64_t offset) override
        {
            assert(this->y < GetMaxY(table_index)+1);

            dest += offset/8;

            if (packed_y_len_bytes <= 2)
            {
                *((uint16_t*)dest) = this->y;
            }
            else if (packed_y_len_bytes <= 4)
            {
                *((uint32_t*)dest) = this->y;
            }
            else
            {
                *((uint64_t*)dest) = this->y;
            }



            if (packed_c_len_bytes <= 4)
            {
                *((uint32_t*)(dest + packed_y_len_bytes)) = c;
            }
            else if (packed_c_len_bytes <= 8)
            {
                *((uint64_t*)(dest + packed_y_len_bytes)) = c;
            }
            else
            {
                *((uint128_t*)(dest + packed_y_len_bytes)) = c;
            }

        }

        inline void unpack(uint8_t * src, uint64_t offset) override
        {
            src += offset/8;
            if (packed_y_len_bytes <= 2)
            {
                this->y = *((uint16_t*)src) % parent::row_divisor;
            }
            else if (packed_y_len_bytes <= 4)
            {
                this->y = *((uint32_t*)src) % parent::row_divisor;
            }
            else
            {
                this->y = *((uint64_t*)src) % parent::row_divisor;
            }

            uint128_t c_mask = 0;
            if (c_len_bits == 128)
            {
                c_mask--;
            }
            else
            {
                c_mask = ((uint128_t)1) << c_len_bits;
            }

            if (packed_c_len_bytes <= 4)
            {
                c = *((uint32_t*)(src + packed_y_len_bytes)) % c_mask;
            }
            else if (packed_c_len_bytes <= 8)
            {
                c = *((uint64_t*)(src + packed_y_len_bytes)) % c_mask;
            }
            else
            {
                c = *((uint128_t*)(src + packed_y_len_bytes)) % c_mask;
            }

        }
    };

    template<int8_t table_index>
    class YCTempPackedEntry : public PackedEntry<GetMaxY(table_index)/kBC, GetMaxY(table_index)+1, GetMeanEntryCount(table_index)>
    {
        using parent = PackedEntry<GetMaxY(table_index)/kBC, GetMaxY(table_index)+1, GetMeanEntryCount(table_index)>;
    public:
        static constexpr uint64_t c_len_bits = GetCLen(table_index);
        uint128_t c;

        static constexpr uint64_t packed_y_len_bytes = ((parent::trimmed_y_len_bits+31)/32)*4;
        static constexpr uint64_t packed_c_len_bytes = ((c_len_bits+31)/32)*4;

        static constexpr uint64_t len_bits = (packed_y_len_bytes + packed_c_len_bytes)*8;

        inline void pack(uint8_t * dest, uint64_t offset) override
        {
            assert(this->y < GetMaxY(table_index)+1);

            dest += offset/8;

            if (packed_y_len_bytes <= 2)
            {
                *((uint16_t*)dest) = this->y;
            }
            else if (packed_y_len_bytes <= 4)
            {
                *((uint32_t*)dest) = this->y;
            }
            else
            {
                *((uint64_t*)dest) = this->y;
            }



            if (packed_c_len_bytes <= 4)
            {
                *((uint32_t*)(dest + packed_y_len_bytes)) = c;
            }
            else if (packed_c_len_bytes <= 8)
            {
                *((uint64_t*)(dest + packed_y_len_bytes)) = c;
            }
            else
            {
                *((uint128_t*)(dest + packed_y_len_bytes)) = c;
            }

        }

        inline void unpack(uint8_t * src, uint64_t offset) override
        {
            src += offset/8;
            if (packed_y_len_bytes <= 2)
            {
                this->y = *((uint16_t*)src) % parent::row_divisor;
            }
            else if (packed_y_len_bytes <= 4)
            {
                this->y = *((uint32_t*)src) % parent::row_divisor;
            }
            else
            {
                this->y = *((uint64_t*)src) % parent::row_divisor;
            }

            uint128_t c_mask = 0;
            if (c_len_bits == 128)
            {
                c_mask--;
            }
            else
            {
                c_mask = ((uint128_t)1) << c_len_bits;
            }

            if (packed_c_len_bytes <= 4)
            {
                c = *((uint32_t*)(src + packed_y_len_bytes)) % c_mask;
            }
            else if (packed_c_len_bytes <= 8)
            {
                c = *((uint64_t*)(src + packed_y_len_bytes)) % c_mask;
            }
            else
            {
                c = *((uint128_t*)(src + packed_y_len_bytes)) % c_mask;
            }

        }
    };

    template<int8_t table_index>
    class LinePointUIDPackedEntry : public PackedEntry<conf.num_rows, GetMaxLinePoint(table_index)+1, GetMeanEntryCount(table_index)*2>
    {
        using parent = PackedEntry<conf.num_rows, GetMaxLinePoint(table_index)+1, GetMeanEntryCount(table_index)*2>;
    public:
        static constexpr uint64_t uid_len_bits = conf.K+2;
   //     static constexpr uint64_t len_bits = parent::trimmed_y_len_bits + uid_len_bits;

        static constexpr uint64_t packed_y_len_bytes = ((parent::trimmed_y_len_bits+31)/32)*4;
        static constexpr uint64_t packed_uid_len_bytes = ((uid_len_bits+31)/32)*4;

        static constexpr uint64_t len_bits = (packed_y_len_bytes + packed_uid_len_bytes)*8;

        uint64_t uid;

        inline void pack(uint8_t * dest, uint64_t offset) override
        {
            assert(this->y < GetMaxLinePoint(table_index)+1);

            dest += offset/8;

            if (packed_y_len_bytes <= 2)
            {
                *((uint16_t*)dest) = this->y;
            }
            else if (packed_y_len_bytes <= 4)
            {
                *((uint32_t*)dest) = this->y;
            }
            else if (packed_y_len_bytes <= 8)
            {
                *((uint64_t*)dest) = this->y;
            }
            else
            {
                *((uint128_t*)dest) = this->y;
            }



            if (packed_uid_len_bytes <= 4)
            {
                *((uint32_t*)(dest + packed_y_len_bytes)) = uid;
            }
            else if (packed_uid_len_bytes <= 8)
            {
                *((uint64_t*)(dest + packed_y_len_bytes)) = uid;
            }
            else
            {
                *((uint128_t*)(dest + packed_y_len_bytes)) = uid;
            }

        }

        inline void unpack(uint8_t * src, uint64_t offset) override
        {
            src += offset/8;


            if (packed_y_len_bytes <= 2)
            {
                this->y = *((uint16_t*)src) % parent::row_divisor);
            }
            else if (packed_y_len_bytes <= 4)
            {
                this->y = *((uint32_t*)src) % parent::row_divisor;
            }
            else if (packed_y_len_bytes <= 8)
            {
                this->y = *((uint64_t*)src) % parent::row_divisor;
            }
            else
            {
                this->y = *((uint128_t*)src) % parent::row_divisor;
            }

            uint128_t uid_mask = 0;
            if (uid_len_bits == 128)
            {
                uid_mask--;
            }
            else
            {
                uid_mask = ((uint128_t)1) << uid_len_bits;
            }

            if (packed_uid_len_bytes <= 4)
            {
                uid = *((uint32_t*)(src + packed_y_len_bytes)) % uid_mask;
            }
            else if (packed_uid_len_bytes <= 8)
            {
                uid = *((uint64_t*)(src + packed_y_len_bytes)) % uid_mask;
            }
            else
            {
                uid = *((uint128_t*)(src + packed_y_len_bytes)) % uid_mask;
            }

        }

    };

    template<int8_t table_index>
    class PositionPackedEntry : public PackedEntry<conf.num_rows, max_entries_per_graph_table, GetMeanEntryCount(table_index)>
    {
        using parent = PackedEntry<conf.num_rows, max_entries_per_graph_table, GetMeanEntryCount(table_index)>;

        // Y of this packed entry represents the UID of the other entry that maps to this entry
    public:
        static constexpr uint64_t pos_len_bits = K+1;

        static constexpr uint64_t packed_y_len_bytes = ((parent::trimmed_y_len_bits+31)/32)*4;
        static constexpr uint64_t packed_pos_len_bytes = ((pos_len_bits+31)/32)*4;

        static constexpr uint64_t len_bits = (packed_y_len_bytes + packed_c_len_bytes)*8;
        uint64_t pos;
        inline void pack(uint8_t * dest, uint64_t offset) override
        {
            assert(this->y < GetMaxY(table_index)+1);

            dest += offset/8;

            if (packed_y_len_bytes <= 2)
            {
                *((uint16_t*)dest) = this->y;
            }
            else if (packed_y_len_bytes <= 4)
            {
                *((uint32_t*)dest) = this->y;
            }
            else
            {
                *((uint64_t*)dest) = this->y;
            }


            if (packed_pos_len_bytes <= 2)
            {
                *((uint16_t*)(dest + packed_y_len_bytes)) = pos;
            }
            else if (packed_pos_len_bytes <= 4)
            {
                *((uint32_t*)(dest + packed_y_len_bytes)) = pos;
            }
            else if (packed_pos_len_bytes <= 8)
            {
                *((uint64_t*)(dest + packed_y_len_bytes)) = pos;
            }
            else
            {
                *((uint128_t*)(dest + packed_y_len_bytes)) = pos;
            }

        }

        inline void unpack(uint8_t * src, uint64_t offset) override
        {
            src += offset/8;
            if (packed_y_len_bytes <= 2)
            {
                this->y = *((uint16_t*)src) % parent::row_divisor;
            }
            else if (packed_y_len_bytes <= 4)
            {
                this->y = *((uint32_t*)src) % parent::row_divisor;
            }
            else
            {
                this->y = *((uint64_t*)src) % parent::row_divisor;
            }

            uint64_t pos_mask = 0;
            if (pos_len_bits == 64)
            {
                pos_mask--;
            }
            else
            {
                pos_mask = ((uint64_t)1) << pos_mask;
            }

            if (packed_c_len_bytes <= 2)
            {
                pos = *((uint16_t*)(src + packed_y_len_bytes)) % pos_mask;
            }
            else if (packed_c_len_bytes <= 4)
            {
                pos = *((uint32_t*)(src + packed_y_len_bytes)) % pos_mask;
            }
            else if (packed_c_len_bytes <= 8)
            {
                pos = *((uint64_t*)(src + packed_y_len_bytes)) % pos_mask;
            }
            else
            {
                pos = *((uint128_t*)(src + packed_y_len_bytes)) % pos_mask;
            }

        }
    };

    static constexpr uint8_t line_point_delta_len_bits = conf.K + 10;
    static constexpr uint8_t finaltable_y_delta_len_bits = conf.K + 10;

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

    using p2_final_positions_type = PackedArray<PackedEntry<1, 1ULL << (conf.K+1), 1>, max_entries_per_graph_table>;
    //using entries_used_type = BasePackedArray<BooleanPackedEntry, max_entries_per_graph_table, 8>;

    //using phase1_new_positions_type = BasePackedArray<PackedEntry<1, (1ULL<<(conf.K+1)), 1>, (1ULL<<(conf.K+1)), 32>;

    template <int8_t table_index>
    using p1_buckets_done_type = BasePackedArray<BooleanPackedEntry, GetMaxY(table_index)/kBC + 1, 8>;


private:
    Buffer* output_buffer;
    uint8_t id[32];
    uint8_t* memo;
    uint32_t memo_size;
    std::vector<std::thread> phase2b_threads;
    std::vector<uint32_t> cpu_ids;
    std::map<uint32_t, vector<uint64_t>> phase1_new_entry_positions;
    std::vector<uint8_t*> entries_used;
    std::vector<p2_final_positions_type>* phase2_final_positions;
    uint64_t final_table_begin_pointers[10];
    uint64_t pointer_table_offset;
    std::vector<std::vector<Park*>> final_parks;

    static void phase1ThreadA(
            uint32_t cpu_id,
            std::atomic<uint64_t>* coordinator,
            const uint8_t* id,
            Penguin<YCPackedEntry<-1>, conf.interlace_factor>* new_penguin
    );

    template <int8_t table_index>
    struct phase1TableData {
        std::map<uint32_t, Penguin<YCPackedEntry<table_index>, conf.interlace_factor>*> bucket_penguins;
        std::map<uint32_t, Penguin<PositionPackedEntry<table_index>, conf.interlace_factor>*> position_penguins;
    };

    template <int8_t table_index>
    static void phase1ThreadB(
            uint32_t cpu_id,
            const uint8_t* id,
            std::atomic<uint64_t>* coordinator,
            phase1TableData<table_index-1> prev_table_data,
            Penguin<YCPackedEntry<table_index>, conf.interlace_factor>* new_yc_penguin,
            Penguin<LinePointUIDPackedEntry<table_index>, conf.interlace_factor>* new_line_point_penguin);

    template <int8_t table_index>
    static void phase1ThreadC(
            uint32_t cpu_id,
            const uint8_t* id,
            std::atomic<uint64_t>* coordinator,
            std::map<uint32_t, Penguin<PositionPackedEntry<table_index>, conf.interlace_factor>*> position_penguins,
            std::map<uint32_t, Penguin<LinePointUIDPackedEntry<table_index>, conf.interlace_factor>*> line_point_penguin,
            std::vector<Park*> *parks,
            Buffer* buffer
    );

    static void phase1ThreadD(
            uint32_t cpu_id,
            std::atomic<uint64_t>* coordinator,
            std::map<uint32_t, uint64_t*>& new_entry_positions,
            std::map<uint32_t, Penguin<YCPackedEntry<5>, conf.interlace_factor>*> line_point_penguins,
            std::vector<DeltaPark<finaltable_y_delta_len_bits> *> *parks,
            Buffer* buffer
    );

    template <int8_t table_index>
    phase1TableData<table_index> phase1DoTable(phase1TableData<table_index-1> prev_table_data);

    template <int8_t table_index>
    static void phase2ThreadA(
            uint32_t cpu_id,
            std::atomic<uint64_t>* coordinator,
            std::vector<Park*>* parks,
            uint8_t* current_entries_used,
            uint8_t* prev_entries_used
    );

    static void phase2ThreadB(
            uint8_t* entries_used,
            p2_final_positions_type* final_positions);

    template <int8_t table_index>
    void phase2DoTable();

    template <int8_t table_index>
    static void phase3ThreadA(
            uint32_t cpu_id,
            std::atomic<uint64_t>* coordinator,
            std::vector<Park*>* temporary_parks,
            uint8_t* entries_used,
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
