#ifndef PARK_HPP
#define PARK_HPP
#include "chiapos_util.hpp"
#include "encoding.hpp"
#include <cstdint>
#include <vector>
#include <thread>
#include <mutex>
#include <sys/mman.h>

class Park {
protected:
    uint8_t* data;
    uint64_t num_entries;
    uint128_t start_value;
public:
    uint64_t start_pos;
    inline Park(uint64_t num_entries_in)
    {
        num_entries = num_entries_in;
        start_value = 0;
        start_pos = 0;
        data = nullptr;
    }
    [[nodiscard]] inline uint64_t size() const { return num_entries; }
    inline void bind(uint8_t* buffer) { data = buffer; }
    virtual void addEntries(std::vector<uint128_t>& src){}
    virtual void readEntries(std::vector<uint128_t>& dest){}
};

template <uint8_t delta_len>
class DeltaPark : public Park {

    static constexpr uint64_t entry_len_bits = delta_len;

public:
    explicit DeltaPark(uint64_t num_entries_in) : Park(num_entries_in){}

    inline uint64_t GetSpaceNeeded()
    {
        uint128_t bits_needed = entry_len_bits * num_entries;
        return (bits_needed + 7) / 8;
    }

    inline void addEntries(std::vector<uint128_t>& src) override
    {
        start_value = src[0];
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
        }
    }
    inline void readEntries(std::vector<uint128_t>& dest) override
    {
        std::span<bitpacker::byte_type> src{data, GetSpaceNeeded()};
        uint128_t last_linepoint = start_value;
        for (uint32_t i = 0; i < num_entries; i++)
        {
            uint128_t delta = bitpacker::extract<uint128_t>(src, i*entry_len_bits, entry_len_bits);
            dest[i] = last_linepoint + delta;
            last_linepoint += delta;
        }
    }
};

template <uint32_t value_len, uint32_t delta_len, uint32_t stub_minus_bits, double max_average_delta, double R>
class CompressedPark : public Park {
    uint32_t getStubsSize()
    {
        return ((num_entries - 1) * (delta_len - stub_minus_bits) + 7) / 8;
    }

    uint32_t getMaxDeltasSize()
    {
        return ((num_entries - 1) * max_average_delta+7) / 8;
    }

public:
    explicit CompressedPark(uint64_t num_entries_in) : Park(num_entries_in){}

    inline uint64_t GetSpaceNeeded()
    {
        return value_len + getStubsSize() + getMaxDeltasSize();
    }

    inline void addEntries(std::vector<uint128_t>& src) override
    {
        start_value = src[0];

        // Since we have approx 2^k line_points between 0 and 2^2k, the average
        // space between them when sorted, is k bits. Much more efficient than storing each
        // line point. This is divided into the stub and delta. The stub is the least
        // significant (k-kMinusStubs) bits, and largely random/incompressible. The small
        // delta is the rest, which can be efficiently encoded since it's usually very
        // small.
        std::vector<uint64_t> park_stubs(src.size()-1);
        std::vector<uint8_t> park_deltas(src.size()-1);

        for (uint32_t i = 1; i < src.size(); i++)
        {
            assert(src[i] > src[i-1]);
            uint128_t big_delta = src[i] - src[i-1];
            uint64_t stub = big_delta & ((1ULL << (delta_len - stub_minus_bits)) - 1);
            uint64_t small_delta = big_delta >> (delta_len - stub_minus_bits);
            assert(small_delta < 256);
            park_deltas[i-1] = small_delta;
            park_stubs[i-1] = stub;
        }

        // Parks are fixed size, so we know where to start writing. The deltas will not go over
        // into the next park.
        uint8_t * dest = data;

        // First entry is static
        uint128_t first_line_point = start_value << (128 - value_len);
        Util::IntTo16Bytes(dest, first_line_point);
        dest += (value_len+7)/8;

        // We use ParkBits instead of Bits since it allows storing more data
        ParkBits park_stubs_bits;
        for (uint32_t i = 0; i < src.size()-1; i++)
        {
            park_stubs_bits.AppendValue(park_stubs[i], (delta_len - stub_minus_bits));
        }
        uint32_t stubs_size = getStubsSize();
        uint32_t stubs_valid_size = cdiv(park_stubs_bits.GetSize(), 8);
        park_stubs_bits.ToBytes(dest);
        memset(dest + stubs_valid_size, 0, stubs_size - stubs_valid_size);
        dest += stubs_size;

        // The stubs are random so they don't need encoding. But deltas are more likely to
        // be small, so we can compress them
        uint8_t *deltas_start = dest + 2;
        size_t deltas_size = Encoding::ANSEncodeDeltas(park_deltas, src.size()-1, R, deltas_start);

        if (!deltas_size) {
            // Uncompressed
            deltas_size = src.size()-1;
            Util::IntToTwoBytesLE(dest, deltas_size | 0x8000);
            memcpy(deltas_start, park_deltas.data(), deltas_size);
        } else {
            // Compressed
            Util::IntToTwoBytesLE(dest, deltas_size);
        }

        dest += 2 + deltas_size;

        assert(GetSpaceNeeded() > (uint64_t)(dest - data));

        memset(dest, 0x00, GetSpaceNeeded() - (dest - data));
    }
    inline void readEntries(std::vector<uint128_t>& dest) override
    {
        // This is the checkpoint at the beginning of the park
        uint8_t * src = data;
        dest[0] = Util::SliceInt128FromBytes(data, 0, value_len);
        src += (value_len+7)/8;

        // Reads EPP stubs
        auto* stubs_bin = src;
        src += getStubsSize();

        // Reads EPP deltas
        uint32_t max_deltas_size_bits = getMaxDeltasSize() * 8;
        auto* deltas_bin = new uint8_t[max_deltas_size_bits / 8];

        // Reads the size of the encoded deltas object
        uint16_t encoded_deltas_size = 0;
        memcpy((uint8_t*)&encoded_deltas_size, src, sizeof(uint16_t));
        src += 2;

        std::vector<uint8_t> deltas;

        if (0x8000 & encoded_deltas_size) {
            // Uncompressed
            encoded_deltas_size &= 0x7fff;
            deltas.resize(encoded_deltas_size);
            memcpy(deltas.data(), src, encoded_deltas_size);
        } else {
            // Compressed
            memcpy(deltas_bin, src, encoded_deltas_size);

            // Decodes the deltas
            deltas = Encoding::ANSDecodeDeltas(deltas_bin, encoded_deltas_size, num_entries - 1, R);
        }
        src += encoded_deltas_size;

        uint32_t start_bit = 0;
        uint8_t stub_size = (delta_len - stub_minus_bits);

        for (uint32_t i = 0; i < num_entries-1; i++)
        {
            uint64_t stub = Util::EightBytesToInt(stubs_bin + start_bit / 8);
            stub <<= start_bit % 8;
            stub >>= 64 - stub_size;
            start_bit += stub_size;

            dest[i+1] = dest[i] + (deltas[i]<<stub_size) + stub;
        }
        delete[] deltas_bin;
    }
};

#endif