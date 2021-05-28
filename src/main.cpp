#include <iostream>
#include <string>
#include "phase1.hpp"

string Strip0x(const string &hex)
{
    if (hex.size() > 1 && (hex.substr(0, 2) == "0x" || hex.substr(0, 2) == "0X")) {
        return hex.substr(2);
    }
    return hex;
}

int main() {
    std::cout << "Hello, World!" << std::endl;
    string id = "022fb42c08c12de3a6af053880199806532e79515f94e83461612101f9412f9e";
    id = Strip0x(id);
    if (id.size() != 64) {
        cout << "Invalid ID, should be 32 bytes (hex)" << endl;
        exit(1);
    }

    Phase1<32>::Run();
    return 0;
}
