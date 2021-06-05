// Copyright 2018 Chia Network Inc

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//    http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
#include "bits.hpp"

#include <stdint.h>

#include <algorithm>
#include <array>
#include <bitset>
#include <iostream>
#include <map>
#include <utility>
#include <vector>


#include "b3/blake3.h"

#include "chacha8.h"
#include "pos_constants.hpp"
#include "chiapos_util.hpp"
#include "calculate_bucket.hpp"


uint16_t L_targets[2][kBC][kExtraBitsPow];

bool table_initialized = false;
void load_tables()
{
    for (uint8_t parity = 0; parity < 2; parity++) {
        for (uint16_t i = 0; i < kBC; i++) {
            uint16_t indJ = i / kC;
            for (uint16_t m = 0; m < kExtraBitsPow; m++) {
                uint16_t yr =
                    ((indJ + m) % kB) * kC + (((2 * m + parity) * (2 * m + parity) + i) % kC);
                L_targets[parity][i][m] = yr;
            }
        }
    }
    table_initialized = true;
}

// Class to evaluate F1
F1Calculator::F1Calculator(uint8_t k, const uint8_t* orig_key)
{
	uint8_t enc_key[32];
	size_t buf_blocks = cdiv(k << kBatchSizes, kF1BlockSizeBits) + 1;
	this->k_ = k;
	this->buf_ = new uint8_t[buf_blocks * kF1BlockSizeBits / 8 + 7];

	// First byte is 1, the index of this table
	enc_key[0] = 1;
	memcpy(enc_key + 1, orig_key, 31);

	// Setup ChaCha8 context with zero-filled IV
	chacha8_keysetup(&this->enc_ctx_, enc_key, 256, NULL);
}

F1Calculator::~F1Calculator()
{
	delete[] buf_;
}


// F1(x) values for x in range [first_x, first_x + n) are placed in res[].
// n must not be more than 1 << kBatchSizes.
void F1Calculator::CalculateBuckets(uint64_t first_x, uint64_t n, uint64_t *res)
{
	uint64_t start = first_x * k_ / kF1BlockSizeBits;
	// 'end' is one past the last keystream block number to be generated
	uint64_t end = cdiv((first_x + n) * k_, kF1BlockSizeBits);
	uint64_t num_blocks = end - start;
	uint32_t start_bit = first_x * k_ % kF1BlockSizeBits;
	uint8_t x_shift = k_ - kExtraBits;

	assert(n <= (1U << kBatchSizes));

	chacha8_get_keystream(&this->enc_ctx_, start, num_blocks, buf_);
	for (uint64_t x = first_x; x < first_x + n; x++) {
		uint64_t y = Util::SliceInt64FromBytes(buf_, start_bit, k_);

		res[x - first_x] = (y << kExtraBits) | (x >> x_shift);

		start_bit += k_;
	}
}


// Class to evaluate F2 .. F7.




