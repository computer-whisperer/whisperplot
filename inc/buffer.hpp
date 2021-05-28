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

#ifndef SRC_CPP_BUFFERS_HPP_
#define SRC_CPP_BUFFERS_HPP_

#include <atomic>
#include <string>
#include <thread>
#include "bits.hpp"

extern std::string buffer_tmpdir;
extern bool buffer_use_tmp;

class Buffer
{
public:
    bool is_file_backed = true;
    bool is_shared;
    uint8_t *data = NULL;
    uint64_t data_len = 0;
    bool remove_on_destroy = false;
    std::atomic<uint64_t>* insert_pos;

    uint64_t entry_len = 0;

    Buffer(const uint64_t size, std::string name="");
    uint64_t GetInsertionOffset(uint64_t len);
    uint64_t PushEntry(Bits bits);
    uint64_t Count();
    ~Buffer();

    void SwapOutAsync();
    void SwapInAsync(bool shared);

    void WaitForSwapIn();
    void WaitForSwapOut();

    uint64_t InsertString(std::string s);
    uint64_t InsertData(void * data, size_t data_len);
    void Truncate(uint64_t new_size);

    private:

	std::thread swapinthread;
	std::thread swapoutthread;
	bool is_swapped = false;
	bool is_swapping = false;

	std::string fname;
	int fd = -1;

    void SwapOut();
    void SwapIn(bool shared);
};


#endif  // SRC_CPP_BUFFERS_HPP_
