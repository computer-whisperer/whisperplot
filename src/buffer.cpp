#include "../inc/buffer.hpp"
#include <sys/mman.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <cassert>
#include <fcntl.h>
#include <unistd.h>
#include <fstream>
#include <string>
#include <iostream>
#include <cerrno>

using namespace std;

std::string buffer_tmpdir = "/tmp";
bool buffer_use_tmp = false;

Buffer::Buffer(const uint64_t size, string name)
{
    data_len = size;
    is_file_backed = buffer_use_tmp || !name.empty();
    if (name.empty() && is_file_backed)
    {
        name = buffer_tmpdir + "/chiapos_tmp_buff_" + to_string(rand()) + ".buf";
        remove_on_destroy = true;
    }
    else
    {
        remove_on_destroy = false;
    }
    fname = name;
    insert_pos = new std::atomic<uint64_t>(0);
    if (is_file_backed)
    {
        fd = open(fname.c_str(), O_RDWR|O_CREAT);
        assert(fd > 0);
        int res = ftruncate(fd, data_len);
        if (res < 0)
        {
            std::cout << "log(-1) failed: " << std::strerror(errno) << '\n';
        }
    }
    SwapIn(true);
}

uint64_t Buffer::GetInsertionOffset(uint64_t len)
{
  uint64_t offset = insert_pos->fetch_add(len);
  assert((offset+len) <= data_len);
  return offset;
}

uint64_t Buffer::PushEntry(Bits bits)
{
	uint64_t offset = GetInsertionOffset(entry_len);
	bits.ToBytes(data + offset);
	return offset;
}

uint64_t Buffer::InsertString(string s)
{
	uint64_t o = GetInsertionOffset(s.size());
    memcpy(data + o, s.data(), s.size());
	return o;
}

uint64_t Buffer::InsertData(void * data_in, size_t data_len)
{
	uint64_t o = GetInsertionOffset(data_len);
    memcpy(data + o, data_in, data_len);
    return o;
}

uint64_t Buffer::Count()
{
	return *insert_pos/entry_len;
}

Buffer::~Buffer()
{
    if (swapinthread.joinable())
	    swapinthread.join();
    if (swapoutthread.joinable())
	    swapoutthread.join();
	SwapOut();
    if (is_file_backed)
    {
        close(fd);
        if (remove_on_destroy)
            remove(fname.c_str());
    }
}

void Buffer::Truncate(uint64_t new_size)
{
    if (!is_file_backed)
        return;
	SwapOutAsync();
	WaitForSwapOut();
	data_len = new_size;
	int res = ftruncate(fd, new_size);
	if (res < 0)
	{
		std::cout << "log(-1) failed: " << std::strerror(errno) << '\n';
	}
}

void Buffer::SwapOut()
{
	is_swapping = true;
    munmap(data, data_len);
    data = NULL;
    is_swapped = true;
    is_swapping = false;
}

void Buffer::SwapOutAsync()
{
    if (!is_file_backed)
        return;
	if (is_swapping && swapinthread.joinable())
		swapinthread.join();
	if (!is_swapping && !is_swapped)
		swapoutthread = thread(&Buffer::SwapOut, this);
}

void Buffer::SwapIn(bool shared)
{
	is_swapping = true;
    is_shared = shared;
	uint8_t flags = shared ? MAP_SHARED : MAP_PRIVATE;
    if (!is_file_backed)
    {
        flags = MAP_ANONYMOUS|MAP_SHARED;
    }
    data = (uint8_t *) mmap(NULL, data_len, PROT_READ|PROT_WRITE, flags, fd, 0);
    assert(data != MAP_FAILED);
    is_swapped = false;
    is_swapping = false;
}

void Buffer::SwapInAsync(bool shared)
{
    if (!is_file_backed)
        return;
	if (is_swapping && swapoutthread.joinable())
		swapoutthread.join();
	if (!is_swapping && is_swapped)
		swapinthread = thread(&Buffer::SwapIn, this, shared);
}

void Buffer::WaitForSwapIn()
{
    if (!is_file_backed)
        return;
    if (is_swapping && swapoutthread.joinable())
		swapoutthread.join();
	if (is_swapping && is_swapped)
		swapinthread.join();
}

void Buffer::WaitForSwapOut()
{
    if (!is_file_backed)
        return;
    if (is_swapping && swapinthread.joinable())
		swapinthread.join();
	if (is_swapping && !is_swapped)
		swapoutthread.join();
}
