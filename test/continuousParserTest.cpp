#include <gtest/gtest.h>

#include "cool_serial/byte_queue.hpp"
#include "cool_serial/continuous_parser.hpp"
#include "cool_serial/bytes.hpp"
#include "cool_serial/cool_message.hpp"
#include "cista/cista.h"

#include <iostream>

// Extra namespace prevents test structs and fixtures from causing collisions
// CoolSerial v1.0.0 defines the first 5 bytes then the next bytes are data
namespace coolSerial::coolMessageTest

{
    TEST(ContinuousParser, Deserialization)
    {
        // Construct data
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

        Bytes kMessageFrame{message.getFrame()};

        std::cout << message;

        const Bytes kMessageData{kMessageFrame.begin() + 5, kMessageFrame.end()};
        EXPECT_EQ(kMessageData, kDataBytes);

        ByteQueue testQueue{std::deque<Byte>{kMessageFrame.begin(), kMessageFrame.end()}};


        // Test the actual parser
        ContinuousParser parser{testQueue};
        parser.update();
        CoolMessageData recoveredData{parser.getCurrentMessage()};

        EXPECT_EQ(recoveredData.dataType, Byte{0});
        EXPECT_EQ(recoveredData.data, kMessageData);
    }
}



