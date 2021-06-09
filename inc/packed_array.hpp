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

template <uint64_t num_rows, uint128_t y_num_states, uint64_t mean_total_entries>
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
            return mean_entries_per_row*1.2;
        }
    }
protected:
    uint128_t y = 0;
public:
    static constexpr uint128_t row_divisor = (y_num_states / num_rows);
    static constexpr uint64_t trimmed_y_len_bits = ceillog2(row_divisor);
    static constexpr uint64_t mean_entries_per_row = mean_total_entries/num_rows;
    uint64_t row = 0;
    static constexpr uint64_t num_rows_v = num_rows;
    static constexpr uint64_t len_bits = (trimmed_y_len_bits == 0) ? 1 : trimmed_y_len_bits;
    static_assert(len_bits > 0);
    static constexpr uint64_t max_entries_per_row = GetMaxEntries();
    PackedEntry() = default;

    inline explicit PackedEntry(uint128_t value_in)
    {
        setY(value_in);
    }
    inline uint128_t getY()
    {
        return y + row_divisor*row;
    }
    inline void setY(uint128_t value)
    {
        row = value/row_divisor;
        y = value%row_divisor;
        assert(value < y_num_states);
    }

    inline void setYtemp(uint128_t value)
    {
        row = value/row_divisor;
        y = value%row_divisor;
    }

    inline virtual void pack(uint8_t * dest, uint64_t offset)
    {
        assert(y < y_num_states);
        std::span<uint8_t> s{dest, (offset + len_bits + 7)/8};
        bitpacker::insert(s, offset, len_bits,  y);
    }

    inline virtual void pack_atomic(std::atomic<uint8_t> * dest, uint64_t offset)
    {
        assert(y < y_num_states);
        std::span<std::atomic<uint8_t>> s{dest, (offset + len_bits + 7)/8};
        bitpacker::insert(s, offset, len_bits,  y);
    }

    inline virtual void unpack(uint8_t * src, uint64_t offset)
    {
        std::span<uint8_t> s{src, (offset + len_bits + 7)/8};
        y = bitpacker::extract<uint128_t>(s, offset, len_bits);
    }

    inline virtual void unpack_atomic(std::atomic<uint8_t> * dest, uint64_t offset)
    {
        std::span<std::atomic<uint8_t>> s{dest, (offset + len_bits + 7)/8};
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

// Setting the backing store to zero is a sufficient initializer
template <class entry_type, uint64_t reserved_count, uint8_t bit_alignment>
class AtomicPackedArray
{
public:
    static constexpr uint64_t entry_len_bits = ((entry_type::len_bits + bit_alignment - 1) / bit_alignment) * bit_alignment;
    static_assert(entry_len_bits > 0);
    static_assert(reserved_count > 0);
    static_assert(bit_alignment > 0);
    static constexpr uint64_t reserved_len_bits = entry_len_bits*reserved_count;
    static constexpr uint64_t reserved_len_bytes = (reserved_len_bits + 7) / 8;
    std::atomic<uint8_t> data[reserved_len_bytes];

    AtomicPackedArray()= default;
    AtomicPackedArray(const AtomicPackedArray &a)
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
            value.pack_atomic(data, i*entry_len_bits);
        }
    }

    inline void set(uint64_t index, entry_type value)
    {
        assert(index < reserved_count);
        value.pack_atomic(data, index*entry_len_bits);
    }

    inline entry_type get(uint64_t i)
    {
        entry_type value;
        assert(i < reserved_count);
        value.unpack_atomic(data, i*entry_len_bits);
        return value;
    }
};


using BooleanPackedEntry = PackedEntry<1, 2, 1>;

template<class T, uint64_t reserved_count>
using PackedArray = BasePackedArray<T, reserved_count, 1>;

#endif //WHISPERPLOT_PACKED_ARRAY_HPP
