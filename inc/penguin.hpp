#ifndef BUCKET_PAGE_INDEX_HPP
#define BUCKET_PAGE_INDEX_HPP

#include <atomic>
#include <sys/mman.h>
#include "calculate_bucket.hpp"
#include "pos_constants.hpp"

#include "atomic_bitpacker.hpp"
#include "bitpacker.hpp"

template<uint8_t K, int8_t table_index> struct BucketPageIndexEntry_YCSizes {
    static constexpr uint32_t y_len_bits = K+kExtraBits;
    static constexpr uint32_t c_len_bits = K;
    static constexpr uint32_t num_entries_per_bc_bucket = 320;
};

template<uint8_t K> struct BucketPageIndexEntry_YCSizes<K, 0> {
    static constexpr uint32_t y_len_bits = K+kExtraBits;
    static constexpr uint32_t c_len_bits = K*2;
    static constexpr uint32_t num_entries_per_bc_bucket = 320;
};

template<uint8_t K> struct BucketPageIndexEntry_YCSizes<K, 1> {
    static constexpr uint32_t y_len_bits = K+kExtraBits;
    static constexpr uint32_t c_len_bits = K*4;
    static constexpr uint32_t num_entries_per_bc_bucket = 320;
};

template<uint8_t K> struct BucketPageIndexEntry_YCSizes<K, 2> {
    static constexpr uint32_t y_len_bits = K+kExtraBits;
    static constexpr uint32_t c_len_bits = K*4;
    static constexpr uint32_t num_entries_per_bc_bucket = 320;
};

template<uint8_t K> struct BucketPageIndexEntry_YCSizes<K, 3> {
    static constexpr uint32_t y_len_bits = K+kExtraBits;
    static constexpr uint32_t c_len_bits = K*3;
    static constexpr uint32_t num_entries_per_bc_bucket = 320;
};

template<uint8_t K> struct BucketPageIndexEntry_YCSizes<K, 4> {
    static constexpr uint32_t y_len_bits = K+kExtraBits;
    static constexpr uint32_t c_len_bits = K*2;
    static constexpr uint32_t num_entries_per_bc_bucket = 320;
};

template<uint8_t K> struct BucketPageIndexEntry_YCSizes<K, 5> {
    static constexpr uint32_t y_len_bits = K;
    static constexpr uint32_t c_len_bits = 0;
    static constexpr uint32_t num_entries_per_bc_bucket = 16000;
};

constexpr unsigned floorlog2(unsigned x)
{
    return x == 1 ? 0 : 1+floorlog2(x >> 1);
}

constexpr unsigned ceillog2(unsigned x)
{
    return x == 1 ? 0 : floorlog2(x - 1) + 1;
}

template<class T, uint64_t len_bits_in> struct SimplePackedEntry
{
    T value;

    static constexpr uint64_t len_bits = len_bits_in;

    inline explicit SimplePackedEntry()
    {

    }

    inline explicit SimplePackedEntry(T value_in)
    {
        value = value_in;
    }

    template <class T1>
    inline void pack(T1 * dest, uint64_t offset)
    {
        std::span<T1> s{dest, (offset + len_bits + 7)/8};
        bitpacker::insert(s, offset, len_bits,  value);
    }

    template <class T1>
    inline void unpack(T1 * src, uint64_t offset)
    {
        std::span<T1> s{src, (offset + len_bits + 7)/8};
        value = bitpacker::extract<T>(s, offset, len_bits);
    }
};

template<uint8_t K, int8_t table_index> struct YCPackedEntry
{
    uint64_t y_offset;
    uint64_t sort_row;
    uint128_t c;

    static struct BucketPageIndexEntry_YCSizes<K, table_index> sizes;
    static constexpr uint64_t num_bc_buckets_per_sort_row = 4096;
    static constexpr uint64_t sort_row_divisor = kBC * num_bc_buckets_per_sort_row;
    static constexpr uint64_t num_sort_rows = (1ULL << sizes.y_len_bits) / sort_row_divisor + 1;
    static constexpr uint64_t max_entries_per_sort_row = sizes.num_entries_per_bc_bucket * num_bc_buckets_per_sort_row;
    static constexpr uint64_t max_entries = max_entries_per_sort_row*num_sort_rows;
    static constexpr uint64_t trimmed_y_len_bits = sizes.y_len_bits - floorlog2(num_sort_rows);
    static constexpr uint64_t len_bits = trimmed_y_len_bits + sizes.c_len_bits;

    inline uint64_t getY()
    {
        return y_offset + sort_row * sort_row_divisor;
    }

    inline void setY(uint64_t y)
    {
        sort_row = y / sort_row_divisor;
        y_offset = y % sort_row_divisor;
    }

    template <class T>
    inline void pack(T * dest, uint64_t offset)
    {
        std::span<T> s{dest, (offset + len_bits + 7)/8};
        if (trimmed_y_len_bits)
        {
            bitpacker::insert(s, offset, trimmed_y_len_bits,  y_offset);
        }
        if (sizes.c_len_bits)
        {
            bitpacker::insert(s, offset + trimmed_y_len_bits, sizes.c_len_bits,  c);
        }
    }

    template <class T>
    inline void unpack(T * src, uint64_t offset)
    {
        std::span<T> s{src, (offset + len_bits + 7)/8};
        if (trimmed_y_len_bits)
        {
            y_offset = bitpacker::extract<uint64_t>(s, offset, trimmed_y_len_bits);
        }
        if (sizes.c_len_bits)
        {
            c = bitpacker::extract<uint128_t>(s, offset + trimmed_y_len_bits, sizes.c_len_bits);
        }
    }

};

template<uint8_t K> struct LinePointEntryUIDPackedEntry
{
    uint64_t sort_row;
    uint128_t line_point_offset;
    uint64_t entry_uid;
    static constexpr uint64_t linepoint_len_bits = K*2;
    static constexpr uint64_t entry_uid_len_bits = K+2;

    static constexpr uint64_t sort_row_divisor = (1ULL << (linepoint_len_bits - 13));
    static constexpr uint64_t num_sort_rows = (1ULL << 13) + 1;
    static constexpr uint64_t max_entries_per_sort_row = ((1ULL << (K + 2)) / num_sort_rows) * 1.2;
    static constexpr uint64_t max_entries = max_entries_per_sort_row*num_sort_rows;
    static constexpr uint64_t trimmed_line_point_len_bits = linepoint_len_bits - floorlog2(num_sort_rows);
    static constexpr uint64_t len_bits = trimmed_line_point_len_bits + entry_uid_len_bits;

    inline uint64_t getLinePoint()
    {
        return line_point_offset + sort_row * sort_row_divisor;
    }

    inline void setLinePoint(uint128_t line_point)
    {
        sort_row = line_point / sort_row_divisor;
        line_point_offset = line_point % sort_row_divisor;
    }

    template <class T>
    inline void pack(T * dest, uint64_t offset)
    {
        std::span<T> s{dest, (offset + len_bits + 7)/8};
        bitpacker::insert(s, offset, trimmed_line_point_len_bits,  line_point_offset);
        bitpacker::insert(s, offset+trimmed_line_point_len_bits, entry_uid_len_bits,  entry_uid);
    }

    template <class T>
    inline void unpack(T * src, uint64_t offset)
    {
        std::span<T> s{src, (offset + len_bits + 7)/8};
        line_point_offset = bitpacker::extract<uint128_t>(s, offset, trimmed_line_point_len_bits);
        entry_uid = bitpacker::extract<uint64_t>(s, offset+trimmed_line_point_len_bits, entry_uid_len_bits);
    }
};

// Setting the backing store to zero is a sufficient initializer
template <class E, class E2, class T, uint64_t reserved_count>
struct BasePackedArray
{
    static constexpr uint64_t reserved_len_bits = T::len_bits*reserved_count;
    static constexpr uint64_t reserved_len_bytes = (reserved_len_bits + 7) / 8;
    E2 count;
    E data[reserved_len_bytes];

    BasePackedArray()= default;
    BasePackedArray(const BasePackedArray &a)
    {
        count = (uint64_t)a.count;
        for (uint64_t i = 0; i < (T::len_bits*count+7)/8; i++)
        {
            data[i] = (uint8_t)a.data[i];
        }
    }

    inline uint64_t size()
    {
        return count;
    }

    inline void clear()
    {
        count = 0;
    }

    inline void fill(T value)
    {
        count = reserved_count;
        for (uint64_t i = 0; i < reserved_count; i++)
        {
            value.pack(data, i*T::len_bits);
        }
    }

    inline uint64_t append(T value)
    {
        uint64_t index = count++;
        assert(index < reserved_count);
        assert((index*T::len_bits + T::len_bits + 7)/8 <= reserved_len_bytes);
        value.pack(data, index*T::len_bits);
        return index;
    }

    inline void set(uint64_t index, T value)
    {
        assert(index < reserved_count);
        if (count < index)
        {
            count = index+1;
        }
        value.pack(data, index*T::len_bits);
    }


    inline T read(uint64_t i)
    {
        T value;
        assert(i < count);
        value.unpack(data, i*T::len_bits);
        return value;
    }
};

using BooleanPackedEntry = SimplePackedEntry<uint8_t, 1>;

template<class T, uint64_t reserved_count>
using AtomicPackedArray = BasePackedArray<std::atomic<uint8_t>, std::atomic<uint64_t>, T, reserved_count>;

template<class T, uint64_t reserved_count>
using PackedArray = BasePackedArray<uint8_t, uint64_t, T, reserved_count>;

template <class T> class Penguin
{
    using Row = AtomicPackedArray<T, T::max_entries_per_sort_row>;

    static constexpr uint64_t row_count = T::num_sort_rows;

    std::atomic<uint32_t> pop_counts[row_count];
    Row * rows;

public:
    uint64_t pops_required_per_row = 1;
    inline Penguin()
    {
        rows = (Row*) mmap(nullptr, sizeof(Row)*row_count, PROT_READ|PROT_WRITE, MAP_ANONYMOUS|MAP_SHARED, -1, 0);
        madvise(rows, sizeof(Row)*row_count, MADV_DONTDUMP);

        for (auto & it : pop_counts)
        {
            it = 0;
        }
    }

    inline ~Penguin()
    {
        munmap(rows, sizeof(Row)*row_count);
    }

    inline uint64_t addEntry(T entry)
    {
        assert(entry.sort_row < T::num_sort_rows);
        return rows[entry.sort_row].append(entry);
    }

    inline uint32_t getCountInRow(uint64_t row_id)
    {
        assert(row_id < row_count);
        return rows[row_id].count;
    }

    inline T readEntry(uint64_t row_id, uint32_t entry_id)
    {
        assert(row_id < row_count);
        T entry = rows[row_id].read(entry_id);
        entry.sort_row = row_id;
        return entry;
    }
    inline void popRow(uint64_t row_id)
    {
        assert(row_id < row_count);
        uint32_t pop_num = pop_counts[row_id].fetch_add(1);
        if (pop_num == pops_required_per_row)
        {
            // Every single bucket from this row has been popped twice, free the memory!
            madvise(rows + row_id, sizeof(Row), MADV_FREE);
        }
    }
    inline uint64_t getUniqueIdentifier(uint64_t bucket_id, uint32_t entry_id)
    {
        return bucket_id + entry_id*T::num_sort_rows;
    }

};


#endif