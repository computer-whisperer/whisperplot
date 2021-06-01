#ifndef PARK_HPP
#define PARK_HPP
#include "util.hpp"
#include <cstdint>
#include <vector>
#include <thread>
#include <mutex>
#include <sys/mman.h>

template <uint8_t delta_len>
class TemporaryPark {
    uint128_t start_linepoint;
    uint8_t* data;
    uint64_t num_entries;
    static constexpr uint64_t entry_len_bits = delta_len;

public:
    uint64_t start_pos;
    inline TemporaryPark(uint64_t num_entries_in)
    {
        num_entries = num_entries_in;
    }
    inline void bind(uint8_t* buffer) { data = buffer; }
    inline uint64_t size() { return num_entries; }
    inline uint64_t GetSpaceNeeded()
    {
        uint128_t bits_needed = entry_len_bits * num_entries;
        return (bits_needed + 7) / 8;
    }

    inline uint128_t bswap_128(uint128_t value)
    {
        auto parts = (uint64_t *)&value;
        return (((uint128_t)bswap_64(parts[0]))<<64) | bswap_64(parts[1]);
    }

    inline void addEntries(uint128_t* src)
    {
        start_linepoint = src[0];
        uint128_t prev_linepoint_written = src[0];

        std::span<bitpacker::byte_type> dest{data, GetSpaceNeeded()};

        for (uint32_t i = 0; i < num_entries; i++)
        {
            // Only worry about not overwriting the previous entry
            assert(prev_linepoint_written <= src[i]);
            uint128_t delta = src[i] - prev_linepoint_written;
            prev_linepoint_written = src[i];
            assert(delta < (1ULL << entry_len_bits));
            bitpacker::insert(dest, i*entry_len_bits, entry_len_bits,  delta);
            /*
            uint64_t bits_offset = i*entry_len_bits;
            uint64_t bytes_offset = bits_offset / 8;
            bits_offset = bits_offset % 8;
            uint8_t mask = ~((1<<(8-bits_offset))-1);
            *(uint128_t *)(data + bytes_offset) = bswap_128((delta << ((128 - entry_len_bits) - bits_offset)) | ((uint128_t)(*(data + bytes_offset)&mask))<<120);
             */
        }
    }
    inline void readEntries(uint128_t* dest)
    {
        std::span<bitpacker::byte_type> src{data, GetSpaceNeeded()};
        uint128_t last_linepoint = start_linepoint;
        for (uint32_t i = 0; i < num_entries; i++)
        {
            /*
            uint64_t bits_offset = i*entry_len_bits;
            uint64_t bytes_offset = bits_offset / 8;
            bits_offset = bits_offset % 8;
            uint128_t mask = (1ULL<<entry_len_bits)-1;
            uint128_t delta = bswap_128(*(uint128_t *)(data + bytes_offset))>>((128-bits_offset)-entry_len_bits);
            delta &= mask;
             */
            uint128_t delta = bitpacker::extract<uint128_t>(src, i*entry_len_bits, entry_len_bits);
            dest[i] = last_linepoint + delta;
            last_linepoint += delta;
        }
    }
};

#endif