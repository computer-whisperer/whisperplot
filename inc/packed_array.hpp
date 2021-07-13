//
// Created by christian on 6/7/21.
//

#ifndef WHISPERPLOT_PACKED_ARRAY_HPP
#define WHISPERPLOT_PACKED_ARRAY_HPP

#include <atomic>
#include "calculate_bucket.hpp"
#include "pos_constants.hpp"


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
