#include <iostream>
#include <string>
#include "phase1.hpp"
#include "cxxopts.hpp"
#include "thread_mgr.hpp"

using namespace std;

void HexToBytes(const string &hex, uint8_t *result)
{
    for (uint32_t i = 0; i < hex.length(); i += 2) {
        string byteString = hex.substr(i, 2);
        uint8_t byte = (uint8_t)strtol(byteString.c_str(), NULL, 16);
        result[i / 2] = byte;
    }
}

string Strip0x(const string &hex)
{
    if (hex.size() > 1 && (hex.substr(0, 2) == "0x" || hex.substr(0, 2) == "0X")) {
        return hex.substr(2);
    }
    return hex;
}

using namespace std;

int main(int argc, char *argv[]) {
    std::cout.setf( std::ios_base::unitbuf );
    cxxopts::Options options(
        "whisperplot", "A better plotter!");

    // Default values
    uint8_t k = 22;
    uint8_t num_threads = 1;
    string filename = "plot.dat";
    string tempdir = ".";
    string finaldir = ".";
    string operation = "help";
    string memo = "0102030405";
    string id = "022fb42c08c12de3a6af053880199806532e79515f94e83461612101f9412f9e";

    options.allow_unrecognised_options().add_options()(
            "k, size", "Plot size", cxxopts::value<uint8_t>(k))(
            "r, threads", "Number of threads", cxxopts::value<uint8_t>(num_threads))(
            "t, tempdir", "Temporary directory", cxxopts::value<string>(tempdir))(
            "d, finaldir", "Final directory", cxxopts::value<string>(finaldir))(
            "f, file", "Filename", cxxopts::value<string>(filename))(
            "m, memo", "Memo to insert into the plot", cxxopts::value<string>(memo))(
            "i, id", "Unique 32-byte seed for the plot", cxxopts::value<string>(id))(
            "help", "Print help");
    auto result = options.parse(argc, argv);

    std::vector<uint32_t> cpu_ids;
    for (uint32_t cpuid = 0; cpuid < num_threads; cpuid++)
    {
        cpu_ids.push_back(cpuid);
    }

    cout << "WhisperPlot!" << endl;
    cout << " k = " << static_cast<int>(k) << endl;
    cout << " filename = " << filename << endl;
    cout << " id = " << id << endl;
    cout << " threads = "<< static_cast<int>(num_threads) << endl;
    cout << " numa nodes = "<< static_cast<int>(GetNUMANodesFromCpuIds(cpu_ids).size()) << endl;

    id = Strip0x(id);
    if (id.size() != 64) {
        cout << "Invalid ID, should be 32 bytes (hex)" << endl;
        exit(1);
    }
    memo = Strip0x(memo);
    if (memo.size() % 2 != 0) {
        cout << "Invalid memo, should be only whole bytes (hex)" << endl;
        exit(1);
    }
    std::vector<uint8_t> memo_bytes(memo.size() / 2);
    std::array<uint8_t, 32> id_bytes;

    HexToBytes(memo, memo_bytes.data());
    HexToBytes(id, id_bytes.data());



    if (k == 18) {
        auto res = new Phase1<18>(id_bytes.data(), cpu_ids);
        cout << res->final_parks.size() << endl;
    }
    else if (k == 22) {
        auto res = new Phase1<22>(id_bytes.data(), cpu_ids);
        cout << res->final_parks.size() << endl;
    }
    else if (k == 26) {
        auto res = new Phase1<26>(id_bytes.data(), cpu_ids);
        cout << res->final_parks.size() << endl;
    }
    else if (k == 32) {
        auto res = new Phase1<32>(id_bytes.data(), cpu_ids);
        cout << res->final_parks.size() << endl;
    }
    else if (k == 33){
        auto res = new Phase1<33>(id_bytes.data(), cpu_ids);
        cout << res->final_parks.size() << endl;
    }



    return 0;
}
