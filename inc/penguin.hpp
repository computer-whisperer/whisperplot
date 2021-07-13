#ifndef BUCKET_PAGE_INDEX_HPP
#define BUCKET_PAGE_INDEX_HPP

#include <atomic>
#include <sys/mman.h>
#include "calculate_bucket.hpp"
#include "pos_constants.hpp"

#include "packed_array.hpp"

template <class entry_type, bool threadsafe>
class Penguin
{
    static constexpr bool use_hugepages = false;

    static constexpr uint64_t hugepage_len = 1ULL << 21;

    static constexpr uint64_t row_len_bytes = (entry_type::max_entries_per_row*entry_type::len_bits+7)/8;

    static constexpr uint64_t row_count = entry_type::num_rows;

    using counter_type = typename std::conditional<threadsafe, std::atomic<uint64_t>, uint64_t>::type;

    counter_type entry_counts[row_count];
    uint8_t* rows[row_count];

    using overlap_entry_id_type = uint64_t;
    static constexpr size_t overlap_row_len = entry_type::max_overlaps_per_row*sizeof(overlap_entry_id_type);
    static constexpr size_t overlap_total_allocation_len =
            use_hugepages ? ((overlap_row_len*row_count*2 + hugepage_len - 1)/hugepage_len)*hugepage_len : overlap_row_len*row_count*2;
    overlap_entry_id_type * overlap_allocation_data;

    counter_type upper_overlap_entry_counts[row_count];
    counter_type lower_overlap_entry_counts[row_count];

    static constexpr uint64_t row_size_bytes = use_hugepages ?  ((row_len_bytes + hugepage_len - 1)/hugepage_len)*hugepage_len : row_len_bytes;

    static uint8_t* do_mmap(uint64_t len_bytes)
    {
        int flags = MAP_ANONYMOUS | MAP_PRIVATE | MAP_NORESERVE;
        if (use_hugepages)
            flags |= MAP_HUGETLB;
        uint8_t* dat = (uint8_t *) mmap(nullptr, len_bytes, PROT_READ | PROT_WRITE, flags, -1, 0);
        if (dat == nullptr)
        {
            throw std::runtime_error("Failed to mmap for penguin");
        }
        if (madvise(dat, len_bytes, MADV_DONTDUMP))
        {
            throw std::runtime_error("Failed to MADV_DONTDUMP for penguin");
        }
        return dat;
    }

public:
    inline Penguin()
    {
        // Allocate memory for overlap rows
        if constexpr(overlap_total_allocation_len > 0)
        {
            overlap_allocation_data = (overlap_entry_id_type*)do_mmap(overlap_total_allocation_len);
        }
        for (uint32_t i = 0; i < row_count; i++) {
            rows[i] = do_mmap(row_size_bytes);
        }
        for (auto& it : entry_counts)
        {
            it = 0;
        }
        for (auto& it : upper_overlap_entry_counts)
        {
            it = 0;
        }
        for (auto& it : lower_overlap_entry_counts)
        {
            it = 0;
        }
    }

    inline ~Penguin()
    {
        if constexpr(overlap_total_allocation_len > 0) {
            munmap(overlap_allocation_data, overlap_total_allocation_len);
        }
        for (uint32_t i = 0; i < row_count; i++)
        {
            if (1 || use_hugepages)
            {
                munmap(rows[i], row_size_bytes);
            }
            else
            {
                madvise(rows[i], row_size_bytes, MADV_FREE);
            }
        }
    }

    inline entry_type newEntry(uint128_t y)
    {
        uint64_t row_id = entry_type::getRowFromY(y);
        assert(row_id < row_count);
        entry_type entry(rows[row_id], entry_counts[row_id]++, row_id);
        entry.setY(y);
        return entry;
    }

    inline uint64_t getCountInRow(uint64_t sort_row_id)
    {
        assert(sort_row_id < row_count);
        return entry_counts[sort_row_id];
    }

    inline entry_type readEntry(uint64_t sort_row_id, uint64_t entry_id)
    {
        assert(sort_row_id < row_count);
        return entry_type(rows[sort_row_id], entry_id, sort_row_id);
    }
    inline void popRow(uint64_t row_id)
    {
        assert(row_id < row_count);
        if (use_hugepages)
        {
            if (munmap(rows[row_id], row_size_bytes))
            {
                std::cerr << "Failed to munmap:" << strerror(errno) << std::endl;
            }
        }
        else
        {
            if (madvise(rows[row_id], row_size_bytes, MADV_FREE))
            {
                std::cerr << "Failed to MADV_FREE:" << strerror(errno) << std::endl;
            }
        }

    }

};


#endif