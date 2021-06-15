//
// Created by christian on 6/14/21.
//
#include "status_update.hpp"
#include <iostream>

static uint64_t last_segment_start;

void StatusUpdate::Init()
{
    last_segment_start = 0;
}

void StatusUpdate::StartSeg(std::string seg_name)
{
    if (last_segment_start > 0)
    {
        EndSeg();
    }
    std::cout << seg_name;
    last_segment_start = time(nullptr);
}

void StatusUpdate::EndSeg()
{
    std::cout << " (" << time(nullptr) - last_segment_start << "s)" << std::endl;
    last_segment_start = 0;
}