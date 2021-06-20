//
// Created by christian on 6/19/21.
//

#ifndef WHISPERPLOT_BITCOPY_HPP
#define WHISPERPLOT_BITCOPY_HPP

template <uint64_t dest_start_bit, uint64_t dest_num_bytes, uint64_t src_start_bit, uint64_t src_num_bytes>
constexpr void bitCopy(uint8_t* dest, uint8_t* src)
{
    // If any of these occurred, something overflowed
    static_assert(dest_start_bit < (1ULL << 60));
    static_assert(src_start_bit < (1ULL << 60));
    static_assert(dest_num_bytes < (1ULL << 60));
    static_assert(src_num_bytes < (1ULL << 60));
    if constexpr ((dest_num_bytes <= 0) || (src_num_bytes <= 0))
    {
        return;
    }
    else if constexpr (dest_start_bit >= 8)
    {
        if constexpr ((dest_start_bit/8) > dest_num_bytes)
        {
            return;
        }
        else
        {
            return bitCopy<dest_start_bit%8, dest_num_bytes-(dest_start_bit/8),
                    src_start_bit, src_num_bytes>
                    (dest+(dest_start_bit/8), src);
        }
    }
    else if constexpr (src_start_bit >= 8)
    {
        if constexpr ((src_start_bit/8) > src_num_bytes)
        {
            return;
        }
        else
        {
            return bitCopy<dest_start_bit, dest_num_bytes,
                    src_start_bit%8, src_num_bytes - (src_start_bit/8)>
                    (dest, src+(src_start_bit/8));
        }
    } else
    {
        dest[0] |= ((src[0] << src_start_bit)&0xFF) >> dest_start_bit;

        const int8_t delta = static_cast<int8_t>(dest_start_bit) - static_cast<int8_t>(src_start_bit);
        if constexpr (dest_num_bytes > 1)
        {
            if constexpr (delta>0)
            {
                dest[1] |= (src[0] << (8-delta));
            }
        }
        if constexpr (src_num_bytes > 1)
        {
            if constexpr (delta<0)
            {
                dest[0] |= (src[1] >> (8+delta));
            }
        }

        if constexpr ((src_num_bytes > 1)&&(dest_num_bytes > 1)) {
            dest[1] |= (((src[1] >> (8-src_start_bit))&0xFF) << (8-dest_start_bit));
        }

        return bitCopy<dest_start_bit, dest_num_bytes-1, src_start_bit, src_num_bytes-1>(dest+1, src+1);
    }
}


#endif //WHISPERPLOT_BITCOPY_HPP
