//
// Created by christian on 7/17/21.
//
#include <iostream>
#include <string>
#include <filesystem>
#include <vector>
#include <atomic>
#include <cstring>
#include <sched.h>
#include <ctime>
#include <chrono>
#include <condition_variable>
#include <math.h>
#include <sys/mman.h>

using namespace std;


constexpr uint64_t num_entries = 1ULL << 28;

uint64_t* big_buffer;

constexpr uint64_t multiplier = 5;

constexpr bool use_hugepages = false;
constexpr bool use_prefetch = true;

void allocBuffer()
{
    auto last_segment_start = std::chrono::high_resolution_clock::now();
    cout << "Allocating and touching memory...";
    uint64_t len_bytes = sizeof(uint64_t)*multiplier*num_entries;
    int flags = MAP_ANONYMOUS | MAP_PRIVATE | MAP_NORESERVE;
    if (use_hugepages)
        flags |= MAP_HUGETLB;
    big_buffer = (uint64_t *) mmap(nullptr, len_bytes, PROT_READ | PROT_WRITE, flags, -1, 0);
    if (big_buffer == nullptr)
    {
        throw std::runtime_error("Failed to mmap");
    }
    if (madvise(big_buffer, len_bytes, MADV_DONTDUMP))
    {
        throw std::runtime_error("Failed to MADV_DONTDUMP");
    }

    memset(big_buffer, 5, sizeof(uint64_t)*multiplier*num_entries);
    uint64_t ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - last_segment_start).count();
    std::cout << " (" << ms << "ms)" << std::endl;
}

void doTest(uint64_t row_count)
{
    vector<uint64_t> entry_counts(row_count);
    uint64_t row_len = multiplier*num_entries/row_count;
    for (auto & it : entry_counts)
    {
        it = 0;
    }

    srand(332);

    auto last_segment_start = std::chrono::high_resolution_clock::now();

    std::cout << row_count << ", ";

    for (uint64_t i = 0; i < num_entries; i++)
    {
        uint64_t val = rand();
        uint64_t row = val%row_count;
        uint64_t entry = entry_counts[row]++;
        if (entry >= row_len)
        {
            throw std::runtime_error("Too many entries");
        }
        big_buffer[row*row_len + entry] = val;
        if (use_prefetch)
            __builtin_prefetch(&(big_buffer[row*row_len + entry+1]));
    }

    uint64_t ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - last_segment_start).count();
    std::cout << ms << ", " << static_cast<uint64_t>(ms/log2((double)row_count)) << std::endl;

}

int main(int argc, char *argv[]) {
    std::cout.setf( std::ios_base::unitbuf );

    int cpuAffinity = argc > 1 ? atoi(argv[1]) : -1;

    if (cpuAffinity > -1)
    {
        cpu_set_t mask;
        int status;

        CPU_ZERO(&mask);
        CPU_SET(cpuAffinity, &mask);
        status = sched_setaffinity(0, sizeof(mask), &mask);
        if (status != 0)
        {
            perror("sched_setaffinity");
        }
    }

    allocBuffer();
    std::cout << "row_count, time, time/log2(row_count)" << endl;
    for (uint32_t i = 5; i < 18; i++)
    {
        doTest(1ULL << i);
    }
}
