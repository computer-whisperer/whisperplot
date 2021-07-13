//
// Created by christian on 6/24/21.
//

#ifndef WHISPERPLOT_PLOTCONF_HPP
#define WHISPERPLOT_PLOTCONF_HPP

struct PlotConf {
    const uint8_t K;
    const uint32_t num_rows;
    const uint32_t interlace_factor;

    constexpr PlotConf(uint8_t K_in,
        uint32_t num_rows_in,
        uint32_t interlace_factor_in):
        K(K_in),
        num_rows(num_rows_in),
        interlace_factor(interlace_factor_in)
    {}

    static constexpr uint8_t max_gid_stub_val = 26;

    constexpr uint64_t getMaxEntriesForGraphTable() const
    {
        return (1ULL << K) * 1.1;
    }

    constexpr uint64_t getMeanEntriesForGraphTable() const
    {
        return 1ULL << K;
    }

    constexpr uint64_t getMaxY(int8_t table_index) const
    {
        if (table_index < 5)
            return 1ULL << (K+kExtraBits);
        else
            return 1ULL << K;
    }

    constexpr uint128_t getMaxGID(int8_t table_index) const {
        if (table_index == -1) {
            return (1ULL << K);
        }
        else
        {
            return getMaxGID(table_index - 1) * max_gid_stub_val;
        }
    }

    constexpr uint64_t getCLen(int8_t table_index) const
    {
        switch (table_index)
        {
            case -1:
                return K;
            case 0:
                return K*2;
            case 1:
            case 2:
                return K*4;
            case 3:
                return K*3;
            case 4:
                return K*2;
            default:
                return 0;
        }
    }

};

#include <iostream>
#include <vector>
#include "penguin.hpp"
#include "buffer.hpp"
#include "park.hpp"
#include "packed_array.hpp"
#include "thread_mgr.hpp"

#endif //WHISPERPLOT_PLOTCONF_HPP
