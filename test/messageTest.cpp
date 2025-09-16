#include <gtest/gtest.h>

#include "cool_serial/bytes.hpp"
#include "cool_serial/cool_message.hpp"
#include "cista/cista.h"

#include <iostream>

// Extra namespace prevents test structs and fixtures from causing collisions
// CoolSerial v1.0.0 defines the first 5 bytes then the next bytes are data
namespace coolSerial::coolMessageTest

{
    TEST(CoolMessage, Serialization)
    {
        struct TestStruct
        {
            float x;
            float y;
            double z;
        };

        TestStruct data
        {
            .x = 2.234,
            .y = -345345.234234,
            .z = 34535.345345
        };

        Bytes kDataBytes{cista::serialize(data)};

        // The 3rd header section byte (index 3) must equal 0
        CoolMessage message{std::move(kDataBytes), Byte{0}};

        const Bytes kMessageFrame{message.getFrame()};

        // std::cout << message;

        const Bytes kMessageData{kMessageFrame.begin() + 5, kMessageFrame.end()};
        EXPECT_EQ(kMessageData, kDataBytes);
    }
}



