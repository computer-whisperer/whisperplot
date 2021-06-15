#include <string>
#include <iostream>
#include <fstream>
#include "plotter.hpp"
#include "entry_sizes.hpp"
// Read a plot back into memory

using namespace std;

struct plot_header {
    uint8_t magic[19];
    uint8_t id[32];
    uint8_t k;
    uint8_t fmt_desc_len[2];
    uint8_t fmt_desc[50];
};

// Using this method instead of simply seeking will prevent segfaults that would arise when
// continuing the process of looking up qualities.
static void SafeSeek(std::ifstream& disk_file, uint64_t seek_location) {
    disk_file.seekg(seek_location);

    if (disk_file.fail()) {
        std::cout << "goodbit, failbit, badbit, eofbit: "
                  << (disk_file.rdstate() & std::ifstream::goodbit)
                  << (disk_file.rdstate() & std::ifstream::failbit)
                  << (disk_file.rdstate() & std::ifstream::badbit)
                  << (disk_file.rdstate() & std::ifstream::eofbit)
                  << std::endl;
        throw std::runtime_error("badbit or failbit after seeking to " + std::to_string(seek_location));
    }
}

static void SafeRead(std::ifstream& disk_file, uint8_t* target, uint64_t size) {
    int64_t pos = disk_file.tellg();
    disk_file.read(reinterpret_cast<char*>(target), size);

    if (disk_file.fail()) {
        std::cout << "goodbit, failbit, badbit, eofbit: "
                  << (disk_file.rdstate() & std::ifstream::goodbit)
                  << (disk_file.rdstate() & std::ifstream::failbit)
                  << (disk_file.rdstate() & std::ifstream::badbit)
                  << (disk_file.rdstate() & std::ifstream::eofbit)
                  << std::endl;
        throw std::runtime_error("badbit or failbit after reading size " +
                                 std::to_string(size) + " at position " + std::to_string(pos));
    }
}

template <PlotConf conf>
template <int8_t table_index>
void Plotter<conf>::readGraphTable(std::ifstream& file_stream, bool decompress)
{
    uint64_t table_len_bytes = final_table_begin_pointers[table_index+1] - final_table_begin_pointers[table_index];
    Buffer* table_buffer;
    if (decompress)
        table_buffer = new Buffer(table_len_bytes*2);
    else
        table_buffer = new Buffer(table_len_bytes);

    SafeSeek(file_stream, final_table_begin_pointers[table_index]);
    uint8_t* decompress_buff;
    if (decompress)
    {
        decompress_buff = (uint8_t*)malloc(table_len_bytes);
        SafeRead(file_stream, decompress_buff, table_len_bytes);
    } else
    {
        SafeRead(file_stream, table_buffer->data, table_len_bytes);
    }

    buffers.push_back(table_buffer);

    vector<uint128_t> line_points;

    // Generate park list
    uint64_t park_start_pos = 0;
    uint64_t park_size_bytes = EntrySizes::CalculateParkSize(conf.K, table_index+1);
    uint64_t last_park_end = 0;
    while (last_park_end < table_len_bytes)
    {
        Park* p;
        if (table_index == 0) {
            using cp = CompressedPark<conf.K * 2, conf.K, kStubMinusBits, (uint32_t)(kMaxAverageDeltaTable1*100), (uint32_t)(kRValues[table_index]*100)>;
            cp* p2 = new cp(kEntriesPerPark);
            p = p2;
        }
        else {
            using cp = CompressedPark<conf.K * 2, conf.K, kStubMinusBits, (uint32_t)(kMaxAverageDelta*100), (uint32_t)(kRValues[table_index]*100)>;
            cp* p2 = new cp(kEntriesPerPark);
            p = p2;
        }
        if (decompress)
        {
            p->bind(decompress_buff + last_park_end);
            line_points.resize(p->size());
            p->readEntries(line_points);
            delete p;

            using dp = DeltaPark<line_point_delta_len_bits>;
            dp* p2 = new dp(kEntriesPerPark);
            p2->bind(table_buffer->data + table_buffer->GetInsertionOffset(p2->GetSpaceNeeded()));
            p2->addEntries(line_points);
            p = p2;
        } else
        {
            p->bind(table_buffer->data + last_park_end);
        }
        p->start_pos = park_start_pos;
        phase1_graph_parks[table_index].push_back(p);

        last_park_end += park_size_bytes;
        park_start_pos += kEntriesPerPark;
    }
}

template <PlotConf conf>
void Plotter<conf>::readFinalTable(std::ifstream& file_stream)
{
    uint32_t P7_park_size = Util::ByteAlign((conf.K + 1) * kEntriesPerPark) / 8;

    uint64_t p7_len_bytes = final_table_begin_pointers[7] - final_table_begin_pointers[6];
    auto p7_buff = (uint8_t*)malloc(p7_len_bytes);
    SafeSeek(file_stream, final_table_begin_pointers[6]);
    SafeRead(file_stream, p7_buff, p7_len_bytes);

    uint64_t c1_len_bytes = final_table_begin_pointers[8] - final_table_begin_pointers[7];
    auto c1_buff = (uint8_t*)malloc(c1_len_bytes);
    SafeSeek(file_stream, final_table_begin_pointers[7]);
    SafeRead(file_stream, c1_buff, c1_len_bytes);

    file_stream.seekg (0, ios::end);
    uint64_t c3_len_bytes = ((uint64_t)file_stream.tellg()) - final_table_begin_pointers[9];
    auto c3_buff = (uint8_t*)malloc(c3_len_bytes);
    SafeSeek(file_stream, final_table_begin_pointers[9]);
    SafeRead(file_stream, c3_buff, c3_len_bytes-1);

    auto table_buffer = new Buffer(p7_len_bytes*4);
    buffers.push_back(table_buffer);

    uint64_t max_possible_pos = (*(phase1_graph_parks[5].end()-1))->start_pos + (*(phase1_graph_parks[5].end()-1))->size();

    // Iterate over c3, read data from p7, and generate the parks
    uint64_t park_size_bytes = EntrySizes::CalculateC3Size(conf.K);
    uint64_t c3_last_park_end = 0;
    uint128_t last_park_last_val = 0;
    uint64_t entry_i = 0;
    uint64_t c1_index = 0;
    vector<uint128_t> line_points;
    while (c3_last_park_end < c3_len_bytes) {
        uint8_t *park_input_data = c3_buff + c3_last_park_end;
        uint128_t raw_values[kCheckpoint1Interval];

        // Grab value from c1
        uint32_t c1_entry_size_bytes = Util::ByteAlign(conf.K) / 8;
        Bits c1_entry_bits = Bits(c1_buff+(c1_index++)*c1_entry_size_bytes, c1_entry_size_bytes, Util::ByteAlign(conf.K));
        raw_values[0] = c1_entry_bits.Slice(0, conf.K).GetValue();

        // Deconstruct c3 park
        uint64_t encoded_size = Bits(park_input_data, 2, 16).GetValue();
        park_input_data += 2;

        assert(encoded_size < park_size_bytes);

        std::vector<uint8_t> deltas = Encoding::ANSDecodeDeltas(park_input_data, encoded_size, kCheckpoint1Interval-1,
                                                                kC3R);
        park_input_data += encoded_size;
        for (uint32_t i = 0; i < deltas.size(); ++i)
        {
            raw_values[i+1] = raw_values[i] + deltas[i];
        }
        line_points.resize(deltas.size()+1);
        uint32_t i;
        for (i = 0; i < deltas.size()+1; ++i)
        {
            // Get matching pos
            uint64_t park_num = (entry_i)/kEntriesPerPark;
            uint64_t park_offset = (entry_i)%kEntriesPerPark;
            entry_i++;
            std::span<uint8_t> s{p7_buff + park_num*P7_park_size, (park_offset*(conf.K+1) + (conf.K+1) + 7)/8};

            if (park_num*P7_park_size >= p7_len_bytes)
            {
                // Out of P7 entries
                break;
            }

            auto pos = bitpacker::extract<uint64_t>(s, park_offset*(conf.K+1), conf.K+1);
            // sanity check
            assert(pos < max_possible_pos);
            line_points[i] = (raw_values[i] << (conf.K+1)) | pos;

            if ((i > 0) && (line_points[i] < line_points[i-1]))
            {
                // Out of C3 entries
                break;
            }

        }

        assert(line_points[0] > last_park_last_val);
        last_park_last_val = line_points[i-1];

        auto p = new DeltaPark<finaltable_y_delta_len_bits>(i);
        p->bind(table_buffer->data + table_buffer->GetInsertionOffset(p->GetSpaceNeeded()));
        p->start_pos = entry_i - i;
        phase1_final_parks.push_back(p);
        p->addEntries(line_points);
        c3_last_park_end += park_size_bytes;
    }

    free(p7_buff);
    free(c3_buff);
    free(c1_buff);
}

uint8_t getK(string plot_fname)
{
    struct plot_header header{};

    std::ifstream disk_file(plot_fname, std::ios::in | std::ios::binary);

    if (!disk_file.is_open()) {
        throw std::invalid_argument("Invalid file " + plot_fname);
    }
    // 19 bytes  - "Proof of Space Plot" (utf-8)
    // 32 bytes  - unique plot id
    // 1 byte    - k
    // 2 bytes   - format description length
    // x bytes   - format description
    // 2 bytes   - memo length
    // x bytes   - memo

    SafeRead(disk_file, (uint8_t*)&header, sizeof(header));
    if (memcmp(header.magic, "Proof of Space Plot", sizeof(header.magic)) != 0)
        throw std::invalid_argument("Invalid plot header magic");

    uint16_t fmt_desc_len = Util::TwoBytesToInt(header.fmt_desc_len);

    if (fmt_desc_len == kFormatDescription.size() &&
        !memcmp(header.fmt_desc, kFormatDescription.c_str(), fmt_desc_len)) {
        // OK
    } else {
        throw std::invalid_argument("Invalid plot file format");
    }

    disk_file.close();
    return header.k;
}

template <PlotConf conf>
void Plotter<conf>::read(string plot_fname, bool decompress) {
    if (decompress)
    {
        cout << "Importing and decompressing plot ";
    }
    else
    {
        cout << "Importing plot ";
    }
    uint64_t part_start_seconds = time(nullptr);

    struct plot_header header{};
    this->filename = plot_fname;

    std::ifstream disk_file(filename, std::ios::in | std::ios::binary);

    if (!disk_file.is_open()) {
        throw std::invalid_argument("Invalid file " + filename);
    }
    // 19 bytes  - "Proof of Space Plot" (utf-8)
    // 32 bytes  - unique plot id
    // 1 byte    - k
    // 2 bytes   - format description length
    // x bytes   - format description
    // 2 bytes   - memo length
    // x bytes   - memo

    SafeRead(disk_file, (uint8_t*)&header, sizeof(header));
    if (memcmp(header.magic, "Proof of Space Plot", sizeof(header.magic)) != 0)
        throw std::invalid_argument("Invalid plot header magic");

    uint16_t fmt_desc_len = Util::TwoBytesToInt(header.fmt_desc_len);

    if (fmt_desc_len == kFormatDescription.size() &&
        !memcmp(header.fmt_desc, kFormatDescription.c_str(), fmt_desc_len)) {
        // OK
    } else {
        throw std::invalid_argument("Invalid plot file format");
    }

    memcpy(this->id, header.id, sizeof(header.id));
    if (conf.K != header.k)
    {
        cout << "K of file is " << header.k << ", not " << conf.K << endl;
        return;
    }
    assert(conf.K == header.k);
    SafeSeek(disk_file, offsetof(struct plot_header, fmt_desc) + fmt_desc_len);

    uint8_t size_buf[2];
    SafeRead(disk_file, size_buf, 2);
    this->memo_size = Util::TwoBytesToInt(size_buf);
    this->memo = new uint8_t[this->memo_size];
    SafeRead(disk_file, this->memo, this->memo_size);

    uint8_t pointer_buf[8];
    for (uint8_t i =0; i < 10; i++) {
        SafeRead(disk_file, pointer_buf, 8);
        this->final_table_begin_pointers[i] = Util::EightBytesToInt(pointer_buf);
    }

    phase1_graph_parks.resize(6);

    // Read in tables
    readGraphTable<0>(disk_file, decompress);
    readGraphTable<1>(disk_file, decompress);
    readGraphTable<2>(disk_file, decompress);
    readGraphTable<3>(disk_file, decompress);
    readGraphTable<4>(disk_file, decompress);
    readGraphTable<5>(disk_file, decompress);
    readFinalTable(disk_file);

    disk_file.close();

    cout << "(" << time(nullptr) - part_start_seconds << "s)" << endl;
}

#include "explicit_templates.hpp"