//
// Created by christian on 6/7/21.
//

#ifndef WHISPERPLOT_PACKED_ARRAY_HPP
#define WHISPERPLOT_PACKED_ARRAY_HPP

#include <atomic>
#include "calculate_bucket.hpp"
#include "pos_constants.hpp"
#include "atomic_bitpacker.hpp"
#include "bitpacker.hpp"

constexpr uint128_t floorlog2(uint128_t x)
{
    return x == 1 ? 0 : 1+floorlog2(x >> 1);
}

constexpr uint128_t ceillog2(uint128_t x)
{
    return x == 1 ? 0 : floorlog2(x - 1) + 1;
}

template <uint64_t num_rows, uint128_t y_num_states, uint64_t mean_total_entries, uint8_t squares>
class PackedEntry
{
    static constexpr uint64_t GetMaxEntries()
    {
        if (num_rows == 1)
        {
            return mean_entries_per_row;
        }
        else if (mean_entries_per_row < 10)
        {
            return mean_entries_per_row*10;
        }
        else if (mean_entries_per_row < 100)
        {
            return mean_entries_per_row*2;
        }
        else
        {
            return mean_entries_per_row*1.4;
        }
    }

    uint64_t row = 0;
protected:
    uint128_t y = 0;
    static constexpr uint64_t getMaxStoredY()
    {
        if (squares > 0)
        {
            return y_num_states;
        }
        else
        {
            return y_num_states / num_rows;
        }
    }
    static constexpr uint64_t getTrimmedYLenBits()
    {
        if (squares > 0)
        {
            return ceillog2(y_num_states);
        }
        else
        {
            return ceillog2(y_num_states / num_rows);
        }
    }
public:
    static const constexpr uint64_t mean_entries_per_row = mean_total_entries/num_rows;

    static const constexpr uint64_t num_rows_v = num_rows;
    static const constexpr uint64_t len_bits = (getTrimmedYLenBits() == 0) ? 1 : getTrimmedYLenBits();
    static const constexpr uint64_t max_entries_per_row = GetMaxEntries();
    PackedEntry() = default;

    inline uint64_t getRow()
    {
        if (squares > 0)
        {
            return getRowFromY(y);
        }
        else
        {
            return row;
        }
    }

    inline void setRow(uint64_t row_in)
    {
        if (squares > 0)
        {
            // Not used
        }
        else
        {
            row = row_in;
        }
    }

    static constexpr uint64_t getRowFromY(uint128_t y_value)
    {
        if (squares > 0)
        {
            uint128_t key_value = y_value;
            uint128_t key_value_max = y_num_states;

            for (uint8_t i = 0; i < squares; i++)
            {
                // clip to at most the top 32 bits so we don't overflow when we square
                uint128_t clip_divisor = key_value_max/(1ULL << 32);
                if (clip_divisor == 0)
                {
                    clip_divisor = 1;
                }
                key_value = key_value/clip_divisor;
                key_value_max = key_value_max/clip_divisor;
                key_value = key_value*(key_value-1);
                key_value_max = key_value_max*(key_value_max-1);
            }
            uint64_t squared_row_divisor = key_value_max/((uint64_t)num_rows);
            return key_value/squared_row_divisor;
        }
        else
        {
            return y_value/(y_num_states / num_rows);
        }
    }

    static constexpr uint128_t getFirstYInRow(uint64_t row_id)
    {
        if (squares > 0)
        {
            throw std::logic_error("Not implemented!");
        }
        else
        {
            return (y_num_states / num_rows)*row_id;
        }
    }

    inline explicit PackedEntry(uint128_t value_in)
    {
        setY(value_in);
    }
    inline uint128_t getY()
    {
        if (squares > 0)
        {
            return y;
        } else
        {
            return y + (y_num_states / num_rows)*row;
        }
    }
    inline void setY(uint128_t value)
    {
        assert(value < y_num_states);
        if (squares > 0)
        {
            y = value;
        } else
        {
            row = value/(y_num_states / num_rows);
            y = value%(y_num_states / num_rows);
        }
    }

    inline virtual void pack(uint8_t * dest, uint64_t offset)
    {
        assert(y < getMaxStoredY());
        std::span<uint8_t> s{dest, (offset + len_bits + 7)/8};
        bitpacker::insert(s, offset, len_bits,  y);
    }

    inline virtual void unpack(uint8_t * src, uint64_t offset)
    {
        std::span<uint8_t> s{src, (offset + len_bits + 7)/8};
        y = bitpacker::extract<uint128_t>(s, offset, len_bits);
    }
};

// Setting the backing store to zero is a sufficient initializer
template <class entry_type, uint64_t reserved_count, uint8_t bit_alignment>
class BasePackedArray
{
public:
    static constexpr uint64_t entry_len_bits = ((entry_type::len_bits + bit_alignment - 1) / bit_alignment) * bit_alignment;
    static constexpr uint64_t reserved_len_bits = entry_len_bits*reserved_count;
    static constexpr uint64_t reserved_len_bytes = (reserved_len_bits + 7) / 8;
    uint8_t data[reserved_len_bytes];

    BasePackedArray()= default;
    BasePackedArray(const BasePackedArray &a)
    {
        for (uint64_t i = 0; i < reserved_len_bytes; i++)
        {
            data[i] = (uint8_t)a.data[i];
        }
    }

    inline void fill(entry_type value)
    {
        for (uint64_t i = 0; i < reserved_count; i++)
        {
            value.pack(data, i*entry_len_bits);
        }
    }

    inline void set(uint64_t index, entry_type value)
    {
        assert(index < reserved_count);
        value.pack(data, index*entry_len_bits);
    }

    inline entry_type get(uint64_t i)
    {
        entry_type value;
        assert(i < reserved_count);
        value.unpack(data, i*entry_len_bits);
        return value;
    }
};

#endif //WHISPERPLOT_PACKED_ARRAY_HPP
