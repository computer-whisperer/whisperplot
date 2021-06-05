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

#include <cassert>
#include <chrono>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "chiapos_util.hpp"


#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#include <processthreadsapi.h>
#include "uint128_t.h"
#else
// __uint__128_t is only available in 64 bit architectures and on certain
// compilers.

// Allows printing of uint128_t
std::ostream &operator<<(std::ostream &strm, uint128_t const &v)
{
    strm << "uint128(" << (uint64_t)(v >> 64) << "," << (uint64_t)(v & (((uint128_t)1 << 64) - 1))
         << ")";
    return strm;
}

#endif

Timer::Timer()
{
	wall_clock_time_start_ = std::chrono::steady_clock::now();
#if _WIN32
	::GetProcessTimes(::GetCurrentProcess(), &ft_[3], &ft_[2], &ft_[1], &ft_[0]);
#else
	cpu_time_start_ = clock();
#endif
}

char *TimerGetNow()
{
	auto now = std::chrono::system_clock::now();
	auto tt = std::chrono::system_clock::to_time_t(now);
	return ctime(&tt);  // ctime includes newline
}

void Timer::PrintElapsed(const std::string &name) const
{
	auto end = std::chrono::steady_clock::now();
	auto wall_clock_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
							 end - this->wall_clock_time_start_)
							 .count();

#if _WIN32
	FILETIME nowft_[6];
	nowft_[0] = ft_[0];
	nowft_[1] = ft_[1];

	::GetProcessTimes(::GetCurrentProcess(), &nowft_[5], &nowft_[4], &nowft_[3], &nowft_[2]);
	ULARGE_INTEGER u[4];
	for (size_t i = 0; i < 4; ++i) {
		u[i].LowPart = nowft_[i].dwLowDateTime;
		u[i].HighPart = nowft_[i].dwHighDateTime;
	}
	double user = (u[2].QuadPart - u[0].QuadPart) / 10000.0;
	double kernel = (u[3].QuadPart - u[1].QuadPart) / 10000.0;
	double cpu_time_ms = user + kernel;
#else
	double cpu_time_ms =
		1000.0 * (static_cast<double>(clock()) - this->cpu_time_start_) / CLOCKS_PER_SEC;
#endif

	double cpu_ratio = static_cast<int>(10000 * (cpu_time_ms / wall_clock_ms)) / 100.0;

	std::cout << name << " " << (wall_clock_ms / 1000.0) << " seconds. CPU (" << cpu_ratio
			  << "%) " << TimerGetNow();
}







#if defined(_WIN32) || defined(__x86_64__)
    void Util::CpuID(uint32_t leaf, uint32_t *regs)
    {
#if defined(_WIN32)
        __cpuid((int *)regs, (int)leaf);
#else
        __get_cpuid(leaf, &regs[0], &regs[1], &regs[2], &regs[3]);
#endif /* defined(_WIN32) */
    }

    bool Util::HavePopcnt(void)
    {
        // EAX, EBX, ECX, EDX
        uint32_t regs[4] = {0};

        Util::CpuID(1, regs);
        // Bit 23 of ECX indicates POPCNT instruction support
        return (regs[2] >> 23) & 1;
    }
#endif /* defined(_WIN32) || defined(__x86_64__) */


