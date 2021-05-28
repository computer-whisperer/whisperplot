cmake_minimum_required(VERSION 3.19)
project(whisperplot)

set(CMAKE_CXX_STANDARD 20)

set(FSE_LIB chiapos_inc/FiniteStateEntropy/lib)
set(FSE_FILES
        ${FSE_LIB}/fse_compress.c
        ${FSE_LIB}/fse_decompress.c
        ${FSE_LIB}/entropy_common.c
        ${FSE_LIB}/hist.c
        )
add_library(fse ${FSE_FILES})

set(BLAKE3_SRC
        chiapos_inc/b3/blake3.c
        chiapos_inc/b3/blake3_portable.c
        chiapos_inc/b3/blake3_dispatch.c
        chiapos_inc/b3/blake3_avx2_x86-64_unix.S
        chiapos_inc/b3/blake3_avx512_x86-64_unix.S
        chiapos_inc/b3/blake3_sse41_x86-64_unix.S
        )

include_directories(inc chiapos_inc)

add_executable(whisperplot
        src/main.cpp
        src/phase1.cpp
        src/buffer.cpp
        src/calculate_bucket.cpp
        chiapos_inc/chacha8.c
        ${BLAKE3_SRC})

find_package(Threads REQUIRED)
add_library(uint128 STATIC chiapos_inc/uint128_t/uint128_t.cpp)
target_include_directories(uint128 PUBLIC chiapos_inc/uint128_t)


target_link_libraries(whisperplot fse stdc++fs Threads::Threads)