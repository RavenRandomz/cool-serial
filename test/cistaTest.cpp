#include "cista/cista.h"
#include <gtest/gtest.h>

#include <vector>

TEST(barebones, cista)
{
    struct BasicStruct
    {
        int integer;
        double x;
        double y;
        double z;
    };

    BasicStruct input
    {
        .integer = 123234,
        .x = 0.234234,
        .y = -2352.2342,
        .z = 0.345345
    };



    std::vector<unsigned char> bytes{cista::serialize(input)};

    BasicStruct output{*cista::deserialize<BasicStruct>(bytes)};

    // TODO Refactor into byte print function
    //for (auto& byte : bytes)
    //{
    //    std::cout << ", ";
    //    std::cout << std::hex << int(byte);

    //}

    EXPECT_EQ(output.integer, input.integer);
}

TEST(barebones, cistaCut)
{
    struct BasicStruct
    {
        int integer;
    };

    BasicStruct input
    {
        .integer = 123234,
    };



    std::vector<unsigned char> bytes{cista::serialize(input)};

    BasicStruct output{*cista::deserialize<BasicStruct>(bytes)};

    // TODO: Refactor into byte print function
    //for (auto& byte : bytes)
    //{
    //    std::cout << ", ";
    //    std::cout << std::hex << int(byte);

    //}

    EXPECT_EQ(output.integer, input.integer);
}
