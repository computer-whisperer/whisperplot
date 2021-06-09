#ifndef BUCKET_PAGE_INDEX_HPP
#define BUCKET_PAGE_INDEX_HPP

#include <atomic>
#include <sys/mman.h>
#include "calculate_bucket.hpp"
#include "pos_constants.hpp"

#include "atomic_bitpacker.hpp"
#include "bitpacker.hpp"
#include "packed_array.hpp"

template <class entry_type>
class Penguin
{
    using Row = BasePackedArray<entry_type, entry_type::max_entries_per_row, 8>;

    static constexpr uint64_t row_count = entry_type::num_rows_v;

    std::atomic<uint64_t> entry_counts[row_count];
    Row * rows;

public:
    inline Penguin()
    {
        rows = (Row*) mmap(nullptr, sizeof(Row)*row_count, PROT_READ|PROT_WRITE, MAP_ANONYMOUS|MAP_SHARED, -1, 0);
        madvise(rows, sizeof(Row)*row_count, MADV_DONTDUMP);
    }

    inline ~Penguin()
    {
        munmap(rows, sizeof(Row)*row_count);
    }

    inline uint64_t addEntry(entry_type entry)
    {
        assert(entry.row < row_count);
        uint64_t idx = entry_counts[entry.row]++;
        assert(idx < entry_type::max_entries_per_row);
        rows[entry.row].set(idx, entry);
        return idx;
    }

    inline uint64_t getCountInRow(uint64_t row_id)
    {
        assert(row_id < row_count);
        return entry_counts[row_id];
    }

    inline entry_type readEntry(uint64_t row_id, uint64_t entry_id)
    {
        assert(row_id < row_count);
        entry_type entry = rows[row_id].get(entry_id);
        entry.row = row_id;
        return entry;
    }
    inline void popRow(uint64_t row_id)
    {
        assert(row_id < row_count);
        madvise(rows + row_id, sizeof(Row), MADV_FREE);
    }
    inline uint64_t getUniqueIdentifier(uint64_t bucket_id, uint32_t entry_id)
    {
        return bucket_id + entry_id*entry_type::num_rows_v;
    }

};


#endif