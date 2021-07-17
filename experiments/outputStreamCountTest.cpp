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
#include <condition_variable>
#include <sys/mman.h>

using namespace std;

uint64_t g_seed = 23;

// Compute a pseudorandom integer.
// Output value in range [0, 32767]
int fast_rand(void) {
    g_seed = (214013*g_seed+2531011);
    return (g_seed>>16)&0x7FFF;
}

constexpr uint64_t num_entries = 1ULL << 28;

uint64_t* big_buffer;

constexpr uint64_t multiplier = 4;

constexpr bool use_hugepages = false;

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

    memset(big_buffer, 5, sizeof(uint64_t)*2*num_entries);
    std::cout << " (" << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - last_segment_start).count() << "ms)" << std::endl;
}

void doTest(uint64_t row_count)
{
    vector<uint64_t> entry_counts(row_count);
    uint64_t row_len = multiplier*num_entries/row_count;
    for (auto & it : entry_counts)
    {
        it = 0;
    }

    g_seed = 332;

    auto last_segment_start = std::chrono::high_resolution_clock::now();

    std::cout << row_count << ", ";

    for (uint64_t i = 0; i < num_entries; i++)
    {
        uint64_t val = fast_rand();
        uint64_t row = val%row_count;
        uint64_t entry = entry_counts[row]++;
        if (entry >= row_len)
        {
            throw std::runtime_error("Too many entries");
        }
        big_buffer[row*row_len + entry] = val;
    }

    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - last_segment_start).count() << std::endl;

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
    for (uint32_t i = 0; i < 18; i++)
    {
        doTest(1ULL << i);
    }
}
