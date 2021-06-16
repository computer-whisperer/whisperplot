#ifndef BUCKET_PAGE_INDEX_HPP
#define BUCKET_PAGE_INDEX_HPP

#include <atomic>
#include <sys/mman.h>
#include "calculate_bucket.hpp"
#include "pos_constants.hpp"

#include "atomic_bitpacker.hpp"
#include "bitpacker.hpp"
#include "packed_array.hpp"



template <class entry_type, uint32_t interlace_factor>
class Penguin
{
    static constexpr bool use_hugepages = false;

    static constexpr uint64_t hugepage_len = 1ULL << 21;

    using Row = BasePackedArray<entry_type, entry_type::max_entries_per_row*interlace_factor, 8>;

    static constexpr uint64_t sort_row_count = entry_type::num_rows_v;
    static constexpr uint64_t store_row_count = (sort_row_count + interlace_factor - 1)/interlace_factor;

    std::atomic<uint32_t> pop_counts[store_row_count];

    std::atomic<uint64_t> entry_counts[sort_row_count];
    Row* rows[store_row_count];

    static constexpr uint64_t row_size_bytes = use_hugepages ?  ((sizeof(Row) + hugepage_len - 1)/hugepage_len)*hugepage_len : sizeof(Row);

public:
    inline Penguin()
    {
        for (uint32_t i = 0; i < store_row_count; i++) {
            int flags = MAP_ANONYMOUS | MAP_PRIVATE | MAP_NORESERVE;
        if (use_hugepages)
                flags |= MAP_HUGETLB;
            uint8_t *test = (uint8_t *) mmap(nullptr, row_size_bytes, PROT_READ | PROT_WRITE, flags, -1, 0);
            rows[i] = (Row *) test;
            if (test == nullptr)
            {
                throw std::runtime_error("Failed to mmap for penguin");
            }

            if (madvise(test, row_size_bytes, MADV_DONTDUMP))
            {
                throw std::runtime_error("Failed to MADV_DONTDUMP for penguin");
            }
        /*    if (use_hugepages)
            {
                if (madvise(test, row_size_bytes, MADV_HUGEPAGE))
                {
                    throw std::runtime_error("Failed to MADV_HUGEPAGE for penguin");
                }
            }*/
        }
        for (auto& it : entry_counts)
        {
            it = 0;
        }
        if (interlace_factor > 1)
        {
            for (auto& it : pop_counts)
            {
                it = 0;
            }
        }
    }

    inline ~Penguin()
    {
        for (uint32_t i = 0; i < store_row_count; i++)
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

    inline uint64_t addEntry(entry_type entry)
    {
        assert(entry.row < sort_row_count);
        uint64_t idx = entry_counts[entry.row]++;
        uint64_t store_row_id = entry.row / interlace_factor;
        uint32_t interlace = entry.row % interlace_factor;
        assert(idx < entry_type::max_entries_per_row);
        rows[store_row_id]->set(interlace + idx*interlace_factor, entry);
        return idx;
    }

    inline uint64_t getCountInRow(uint64_t sort_row_id)
    {
        assert(sort_row_id < sort_row_count);
        return entry_counts[sort_row_id];
    }

    inline entry_type readEntry(uint64_t sort_row_id, uint64_t entry_id)
    {
        assert(sort_row_id < sort_row_count);
        uint64_t store_row_id = sort_row_id / interlace_factor;
        uint32_t interlace = sort_row_id % interlace_factor;
        entry_type entry = rows[store_row_id]->get(interlace + entry_id*interlace_factor);
        entry.row = sort_row_id;
        return entry;
    }
    inline void popRow(uint64_t row_id)
    {
        assert(row_id < sort_row_count);

        uint64_t store_row_id = row_id/interlace_factor;
        if (interlace_factor > 1)
            pop_counts[store_row_id]++;

        if ((interlace_factor == 1) || (pop_counts[store_row_id] == interlace_factor))
        {
            if (use_hugepages)
            {
                if (munmap(rows[store_row_id], row_size_bytes))
                {
                    std::cerr << "Failed to munmap:" << strerror(errno) << std::endl;
                }
            }
            else
            {
                if (madvise(rows[store_row_id], row_size_bytes, MADV_FREE))
                {
                    std::cerr << "Failed to MADV_FREE:" << strerror(errno) << std::endl;
                }
            }
        }
    }
    inline uint64_t getUniqueIdentifier(uint64_t sort_row_id, uint32_t entry_id)
    {
        return sort_row_id + entry_id * entry_type::num_rows_v;
    }

};


#endif