#ifndef BUCKET_PAGE_INDEX_HPP
#define BUCKET_PAGE_INDEX_HPP

#include <atomic>
#include <sys/mman.h>
#include "calculate_bucket.hpp"
#include "pos_constants.hpp"

template<uint8_t K, int8_t table_index> struct BucketPageIndexEntry_YCSizes {
    static constexpr uint32_t y_len_bits = K+kExtraBits;
    static constexpr uint32_t c_len_bits = K;
};

template<uint8_t K> struct BucketPageIndexEntry_YCSizes<K, 0> {
    static constexpr uint32_t y_len_bits = K+kExtraBits;
    static constexpr uint32_t c_len_bits = K*2;
};

template<uint8_t K> struct BucketPageIndexEntry_YCSizes<K, 1> {
    static constexpr uint32_t y_len_bits = K+kExtraBits;
    static constexpr uint32_t c_len_bits = K*4;
};

template<uint8_t K> struct BucketPageIndexEntry_YCSizes<K, 2> {
    static constexpr uint32_t y_len_bits = K+kExtraBits;
    static constexpr uint32_t c_len_bits = K*4;
};

template<uint8_t K> struct BucketPageIndexEntry_YCSizes<K, 3> {
    static constexpr uint32_t y_len_bits = K+kExtraBits;
    static constexpr uint32_t c_len_bits = K*3;
};

template<uint8_t K> struct BucketPageIndexEntry_YCSizes<K, 4> {
    static constexpr uint32_t y_len_bits = K+kExtraBits;
    static constexpr uint32_t c_len_bits = K*2;
};

template<uint8_t K> struct BucketPageIndexEntry_YCSizes<K, 5> {
    static constexpr uint32_t y_len_bits = K;
    static constexpr uint32_t c_len_bits = 0;
};

constexpr unsigned floorlog2(unsigned x)
{
    return x == 1 ? 0 : 1+floorlog2(x >> 1);
}

constexpr unsigned ceillog2(unsigned x)
{
    return x == 1 ? 0 : floorlog2(x - 1) + 1;
}

template<uint8_t K, int8_t table_index> struct YCBucketEntry
{
    uint64_t y;
    uint128_t c;


    static struct BucketPageIndexEntry_YCSizes<K, table_index> sizes;
    static constexpr uint32_t bucket_divisor = kBC;
    static constexpr uint32_t num_buckets = (1ULL << sizes.y_len_bits)/bucket_divisor + 1;
    static constexpr uint32_t max_entries_per_bucket = 350;
    static constexpr uint32_t trimmed_y_len_bits = sizes.y_len_bits - floorlog2(num_buckets);
    static constexpr uint32_t packed_entry_len_bits = trimmed_y_len_bits + sizes.c_len_bits;
    static constexpr uint32_t packed_entry_len_bytes = (packed_entry_len_bits+7)/8;

    inline uint32_t get_bucket_id()
    {
        return y/bucket_divisor;
    }

    inline void pack(uint8_t * dest)
    {
        auto b = Bits();
        if (trimmed_y_len_bits) {
            uint64_t y_trimmed = y % bucket_divisor;
            b += Bits(y_trimmed, trimmed_y_len_bits);
        }
        if (sizes.c_len_bits)
            b += Bits(c, sizes.c_len_bits);
        b.ToBytes(dest);
    }

    inline void unpack(uint8_t * src, uint32_t bucket_id)
    {
        uint32_t bits_offset = 0;
        if (trimmed_y_len_bits)
        {
            uint64_t y_trimmed = Util::SliceInt64FromBytes(src, bits_offset, trimmed_y_len_bits);
            y = y_trimmed + bucket_id*bucket_divisor;
        }
        bits_offset += trimmed_y_len_bits;
        if (sizes.c_len_bits)
            c = Util::SliceInt128FromBytes(src, bits_offset, sizes.c_len_bits);
    }
};

template<uint8_t K> struct LinePointEntryUIDBucketEntry
{
    uint128_t line_point;
    uint64_t entry_uid;
    static constexpr uint64_t linepoint_len_bits = K*2;
    static constexpr uint64_t entry_uid_len_bits = K+2;

    static constexpr uint64_t bucket_divisor = (1ULL << (linepoint_len_bits-13));
    static constexpr uint64_t num_buckets = (1ULL << 13) + 1;
    static constexpr uint64_t max_entries_per_bucket = ((1ULL<<(K+2))/num_buckets)*1.2;
    static constexpr uint64_t trimmed_line_point_len_bits = linepoint_len_bits - floorlog2(num_buckets);
    static constexpr uint64_t packed_entry_len_bits = trimmed_line_point_len_bits + entry_uid_len_bits;
    static constexpr uint64_t packed_entry_len_bytes = (packed_entry_len_bits+7)/8;

    inline uint32_t get_bucket_id()
    {
        return line_point / bucket_divisor;
    }

    inline void pack(uint8_t * dest)
    {
        auto b = Bits();
        uint64_t linepoint_trimmed = line_point % bucket_divisor;
        b += Bits(linepoint_trimmed, trimmed_line_point_len_bits);
        b += Bits(entry_uid, entry_uid_len_bits);
        b.ToBytes(dest);
    }

    inline void unpack(uint8_t * src, uint32_t bucket_id)
    {
        uint64_t linepoint_trimmed = Util::SliceInt64FromBytes(src, 0, trimmed_line_point_len_bits);
        line_point = linepoint_trimmed + bucket_id * bucket_divisor;
        entry_uid = Util::SliceInt64FromBytes(src, trimmed_line_point_len_bits, entry_uid_len_bits);
    }
};

template <class T> class BucketPageIndex
{
    static constexpr uint64_t buckets_per_page = 1000;
    static constexpr uint64_t num_index_cols = T::max_entries_per_bucket;
    static constexpr uint64_t num_index_rows = T::num_buckets/buckets_per_page + 1;

    static constexpr uint64_t entry_len_bytes = T::packed_entry_len_bytes;
    static constexpr uint64_t page_size = entry_len_bytes*buckets_per_page;
    static constexpr uint64_t index_row_size = page_size*num_index_cols;
    static constexpr uint64_t index_size = index_row_size*num_index_rows;

    // Otherwise known as, how many entries are in buckets below me?
    uint64_t bucket_offsets[T::num_buckets];

    std::atomic<uint32_t> bucket_entry_counts[T::num_buckets];
    std::atomic<uint32_t> row_pop_counts[num_index_rows];
    uint8_t * pages;


public:
    uint8_t pops_required_per_bucket = 1;
    BucketPageIndex()
    {
        pages = (uint8_t*) mmap(nullptr, index_size, PROT_READ|PROT_WRITE, MAP_ANONYMOUS|MAP_SHARED, -1, 0);

        for (auto & it : bucket_entry_counts)
        {
            it = 0;
        }
        for (auto & it : row_pop_counts)
        {
            it = 0;
        }
    }

    ~BucketPageIndex()
    {
        munmap(pages, index_size);
    }

    void SumBucketOffsets()
    {
        uint64_t ctr = 0;
        for (uint64_t i = 0; i < T::num_buckets; i++)
        {
            bucket_offsets[i] = ctr;
            ctr += bucket_entry_counts[i];
        }
    }

    uint64_t GetBucketOffset(uint64_t bucket_id)
    {
        return bucket_offsets[bucket_id];
    }

    inline uint64_t InsertEntry(T entry)
    {
        uint64_t bucket_id = entry.get_bucket_id();
        uint32_t entry_id = bucket_entry_counts[bucket_id].fetch_add(1);
        assert(entry_id < T::max_entries_per_bucket);
        assert(bucket_id < T::num_buckets);
        uint8_t* dest = FindEntry(bucket_id, entry_id);
        entry.pack(dest);
        return entry_id;
    }

    inline uint32_t GetCountInBucket(uint64_t bucket_id)
    {
        return bucket_entry_counts[bucket_id];
    }

    inline T ReadEntry(uint64_t bucket_id, uint32_t entry_id)
    {
        assert(entry_id < T::max_entries_per_bucket);
        assert(bucket_id < T::num_buckets);
        uint8_t* src = FindEntry(bucket_id, entry_id);
        T entry;
        entry.unpack(src, bucket_id);
        return entry;
    }
    inline void PopBucket(uint64_t bucket_id)
    {
        uint32_t index_row = bucket_id/buckets_per_page;
        uint32_t pop_num = row_pop_counts[index_row].fetch_add(1);
        if (pop_num == buckets_per_page*pops_required_per_bucket)
        {
            // Every single bucket from this row has been popped twice, free the memory!
            madvise(pages + index_row_size*index_row, index_row_size, MADV_FREE);
        }
    }
    inline uint64_t GetUniqueIdentifier(uint64_t bucket_id, uint32_t entry_id)
    {
        return bucket_id + entry_id*T::num_buckets;
    }

private:
    [[nodiscard]] inline uint8_t * FindEntry(uint64_t bucket_id, uint32_t entry_id) const
    {
        uint32_t index_row = bucket_id/buckets_per_page;
        uint32_t index_col = entry_id;
        uint32_t page_col = bucket_id%buckets_per_page;
        return pages + index_row_size*index_row + page_size*index_col + entry_len_bytes*page_col;
    }

};


#endif