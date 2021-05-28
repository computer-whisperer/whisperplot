cmake_minimum_required(VERSION 3.19)
project(whisperplot)

set(CMAKE_CXX_STANDARD 20)

include_directories(inc)

add_executable(whisperplot src/main.cpp src/phase1.cpp src/buffers.cpp)