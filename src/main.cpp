#include <iostream>
#include <string>
#include <filesystem>
#include "plotter.hpp"
#include "cxxopts.hpp"
#include "thread_mgr.hpp"
#include "bls.hpp"
#include "chiapos_util.hpp"
#include "util.hpp"
using namespace bls;

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

template<uint8_t K>
void doPlot(const uint8_t* id_in, const uint8_t* memo_in, const uint32_t memo_size_in, std::vector<uint32_t> cpu_ids, std::string filename)
{
    auto res = new Plotter<K>(id_in, memo_in, memo_size_in, cpu_ids, filename);
    res->phase1();
    res->check();
    res->phase2();
    res->phase3();
    res->phase4();
    cout << res->phase1_final_parks.size() << endl;
}

uint8_t hex_lookup[] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F'};
void print_data_as_hex(uint8_t * data, uint32_t data_len)
{
    uint32_t i = 0;
    while (i < data_len)
    {
        putchar(hex_lookup[(data[i]&0xF0)>>4]);
        putchar(hex_lookup[(data[i]&0xF)]);
        i++;
        if (i%0x20 == 0)
        {
            putchar('\n');
        }
    }
}


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
    string memo = "0102030405";
    string id = "022fb42c08c12de3a6af053880199806532e79515f94e83461612101f9412f9e";
    string farmer_key;
    string pool_key;
    string test_private_key;

    options.allow_unrecognised_options().add_options()(
            "k, size", "Plot size", cxxopts::value<uint8_t>(k))(
            "r, threads", "Number of threads", cxxopts::value<uint8_t>(num_threads))(
            "t, tempdir", "Temporary directory", cxxopts::value<string>(tempdir))(
            "d, finaldir", "Final directory", cxxopts::value<string>(finaldir))(
            "f, file", "Filename", cxxopts::value<string>(filename))(
            "m, memo", "Memo to insert into the plot", cxxopts::value<string>(memo))(
            "i, id", "Unique 32-byte seed for the plot", cxxopts::value<string>(id))(
            "a, farmerkey", "The farmer key", cxxopts::value<string>(farmer_key))(
            "o, poolkey", "The pool key", cxxopts::value<string>(pool_key))(
           // "tpk, testpk", "Test private key", cxxopts::value<string>(test_private_key))(
            "help", "Print help");
    auto result = options.parse(argc, argv);

    std::vector<uint32_t> cpu_ids;
    for (uint32_t cpuid = 0; cpuid < num_threads; cpuid++)
    {
        cpu_ids.push_back(cpuid);
    }


    std::array<uint8_t, 32> id_bytes;
    std::vector<uint8_t> memo_bytes;
    if (!farmer_key.empty())
    {
        // Generate an ID from keys

        farmer_key = Strip0x(farmer_key);
        std::vector<uint8_t> farmer_key_bytes(farmer_key.size() / 2);
        HexToBytes(farmer_key, farmer_key_bytes.data());

        pool_key = Strip0x(pool_key);
        std::vector<uint8_t> pool_key_bytes(pool_key.size() / 2);
        HexToBytes(pool_key, pool_key_bytes.data());

        // Example seed, used to generate private key. Always use
        // a secure RNG with sufficient entropy to generate a seed (at least 32 bytes).
        vector<uint8_t> seed = {0,  50, 6,  244, 24,  199, 1,  25,  52,  88,  192,
                                19, 18, 12, 89,  6,   220, 18, 102, 58,  209, 82,
                                12, 62, 89, 110, 182, 9,   44, 20,  254, 22};

        PrivateKey sk = AugSchemeMPL().KeyGen(seed);

        G1Element farmer_key_pk = G1Element::FromByteVector(farmer_key_bytes);
        G1Element pool_key_pk = G1Element::FromByteVector(pool_key_bytes);
        PrivateKey local_key_sk = AugSchemeMPL().DeriveChildSk(sk, 12381);
        local_key_sk = AugSchemeMPL().DeriveChildSk(local_key_sk, 8444);
        local_key_sk = AugSchemeMPL().DeriveChildSk(local_key_sk, 4);
        local_key_sk = AugSchemeMPL().DeriveChildSk(local_key_sk, 0);
        G1Element plot_key_pk = local_key_sk.GetG1Element() + farmer_key_pk;
        auto t = pool_key_pk.Serialize();
        auto t2 = plot_key_pk.Serialize();
        t.insert(t.end(), t2.begin(), t2.end());
        std::vector<unsigned char> hash(32);
        bls::Util::Hash256(id_bytes.data(), t.data(), t.size());

        t = pool_key_pk.Serialize();
        t2 = farmer_key_pk.Serialize();
        auto t3 = sk.Serialize();
        t.insert(t.end(), t2.begin(), t2.end());
        t.insert(t.end(), t3.begin(), t3.end());
        memo_bytes = t;
    }
    else
    {
        // use plot ID and memo from args
        id = Strip0x(id);
        if (id.size() != 64) {
            cout << "Invalid ID, should be 32 bytes (hex)" << endl;
            exit(1);
        }
        HexToBytes(id, id_bytes.data());

        memo = Strip0x(memo);
        if (memo.size() % 2 != 0) {
            cout << "Invalid memo, should be only whole bytes (hex)" << endl;
            exit(1);
        }
        memo_bytes.resize(memo.size() / 2);
        HexToBytes(memo, memo_bytes.data());
    }

    cout << "WhisperPlot!" << endl;
    cout << " k = " << static_cast<int>(k) << endl;
    cout << " filename = " << filename << endl;
    cout << " id = 0x";
    print_data_as_hex(id_bytes.data(), id_bytes.size());
    cout << " memo = 0x";
    print_data_as_hex(memo_bytes.data(), memo_bytes.size());
    cout << " threads = "<< static_cast<int>(num_threads) << endl;
    cout << " numa nodes = "<< static_cast<int>(GetNUMANodesFromCpuIds(cpu_ids).size()) << endl;


    string final_filename = filesystem::path(finaldir) / filesystem::path(filename);

    switch (k)
    {
        case 18:
            doPlot<18>(id_bytes.data(), memo_bytes.data(), memo_bytes.size(), cpu_ids, final_filename);
            break;
        case 22:
            doPlot<22>(id_bytes.data(), memo_bytes.data(), memo_bytes.size(), cpu_ids, final_filename);
            break;
        case 26:
            doPlot<26>(id_bytes.data(), memo_bytes.data(), memo_bytes.size(), cpu_ids, final_filename);
            break;
        case 32:
            doPlot<32>(id_bytes.data(), memo_bytes.data(), memo_bytes.size(), cpu_ids, final_filename);
            break;
        default:
            cout << "Unsupported k selected, please choose from 18, 22, 26, or 32." << endl;
    }





    return 0;
}

