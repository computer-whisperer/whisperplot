//
// Created by christian on 5/29/21.
//

#ifndef WHISPERPLOT_THREAD_MGR_HPP
#define WHISPERPLOT_THREAD_MGR_HPP
#include <vector>
#include <numa.h>

inline void PinToCpuid(uint32_t cpuid)
{
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(cpuid, &cpuset);
    pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
}

inline std::vector<uint32_t> GetNUMANodesFromCpuIds(const std::vector<uint32_t>& cpuids)
{
    std::vector<uint32_t> numa_nodes;
    for (auto& cpuid : cpuids)
    {
        uint32_t numa_node = numa_node_of_cpu(cpuid);
        if (std::find(numa_nodes.begin(), numa_nodes.end(), numa_node) == numa_nodes.end())
        {
            numa_nodes.push_back(numa_node);
        }
    }
    return numa_nodes;
}

#endif //WHISPERPLOT_THREAD_MGR_HPP
