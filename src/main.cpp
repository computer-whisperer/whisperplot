#include <iostream>
#include <string>
#include "phase1.hpp"

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

int main() {
    std::cout << "WhisperPlot!" << std::endl;
    string id = "022fb42c08c12de3a6af053880199806532e79515f94e83461612101f9412f9e";
    id = Strip0x(id);
    if (id.size() != 64) {
        cout << "Invalid ID, should be 32 bytes (hex)" << endl;
        exit(1);
    }

    string memo = "MEMO";

    std::vector<uint8_t> memo_bytes(memo.size() / 2);
    std::array<uint8_t, 32> id_bytes;

    HexToBytes(memo, memo_bytes.data());
    HexToBytes(id, id_bytes.data());

    Phase1<32> phase_1(id_bytes.data(), 48);
    return 0;
}
