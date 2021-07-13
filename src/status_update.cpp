//
// Created by christian on 6/14/21.
//
#include "status_update.hpp"
#include <iostream>
#include <chrono>

static std::chrono::high_resolution_clock::time_point last_segment_start;
static bool has_started = false;

void StatusUpdate::Init()
{
    has_started = false;
}

void StatusUpdate::StartSeg(std::string seg_name)
{
    if (has_started)
    {
        EndSeg();
    }
    std::cout << seg_name;
    last_segment_start = std::chrono::high_resolution_clock::now();
    has_started = true;
}

void StatusUpdate::EndSeg()
{
    std::cout << " (" << duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - last_segment_start).count() << "ms)" << std::endl;
    has_started = false;
}