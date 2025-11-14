#include <gtest/gtest.h>

#include "cool_serial/byte_queue.hpp"
#include "cool_serial/continuous_parser.hpp"
#include "cool_serial/bytes.hpp"
#include "cool_serial/cool_message.hpp"
#include "cista/cista.h"

#include <iostream>

// Extra namespace prevents test structs and fixtures from causing collisions
// CoolSerial v1.0.0 defines the first 5 bytes then the next bytes are data
namespace coolSerial::coolParserTest
{
    // TODO: Add to util or ByteQueue class
    //
    void addBytesToByteQueue(Bytes& bytes, ByteQueue& byteQueue)
    {
        for(const Byte& kByte : bytes)
        {
            byteQueue.push(kByte);
        }
    }
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

        const Bytes kMessageData{kMessageFrame.begin() + 5, kMessageFrame.end()};
        EXPECT_EQ(kMessageData, kDataBytes);

        ByteQueue testQueue{std::deque<Byte>{kMessageFrame.begin(), kMessageFrame.end()}};


        // Test the actual parser
        ContinuousParser parser{testQueue};
        parser.update();
        CoolMessageData recoveredData{parser.getCurrentMessage()};

        EXPECT_EQ(recoveredData.type, Byte{0});
        EXPECT_EQ(recoveredData.data, kMessageData);
    }

    TEST(ContinuousParser, continousDeserilization)
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

        const Bytes kMessageData{kMessageFrame.begin() + 5, kMessageFrame.end()};
        EXPECT_EQ(kMessageData, kDataBytes);


        // A partial buffer will be fed each time
        // [ 89, 0, 10, 0, 57, db, f9, e, 40, 27, a0, a8, c8, 1b, f5, 10, d, eb, dc, e0, 40, ]
        // Partial header with SOF
        Bytes partialHeader{ 0x84, 0x0, 0x10};
        Bytes restOfHeader{ 0x0, 0x57};
        Bytes halfOfData{ 0xdb, 0xf9, 0xe, 0x40};
        Bytes finalHalfOfData{0x27, 0xa0, 0xa8, 0xc8, 0x1b, 0xf5, 0x10, 0xd, 0xeb, 0xdc, 0xe0, 0x40};

        ByteQueue testQueue{};
        // Test the actual parser
        ContinuousParser parser{testQueue};

        // Test if the message processed is false
        parser.reportMessageProcessed();
        addBytesToByteQueue(partialHeader, testQueue);
        parser.update();
        EXPECT_TRUE(parser.currentMessageProcessed());

        addBytesToByteQueue(restOfHeader, testQueue);
        parser.update();
        EXPECT_TRUE(parser.currentMessageProcessed());

        addBytesToByteQueue(halfOfData, testQueue);
        parser.update();
        EXPECT_TRUE(parser.currentMessageProcessed());

        addBytesToByteQueue(finalHalfOfData, testQueue);
        parser.update();

        EXPECT_FALSE(parser.currentMessageProcessed());
        CoolMessageData recoveredData{parser.getCurrentMessage()};

        EXPECT_EQ(recoveredData.type, Byte{0});
        EXPECT_EQ(recoveredData.data, kMessageData);
    }

    TEST(ContinuousParser, ReallyFragmentedDeserialization)
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

        const Bytes kMessageData{kMessageFrame.begin() + 5, kMessageFrame.end()};
        EXPECT_EQ(kMessageData, kDataBytes);

        ByteQueue testQueue{std::deque<Byte>{kMessageFrame.begin(), kMessageFrame.end()}};


        ByteQueue simulatedQueue{};
        // Test the actual parser
        ContinuousParser parser{simulatedQueue};

        while(!testQueue.empty())
        {
            simulatedQueue.push(testQueue.front());
            testQueue.pop();
            parser.update();
        }

        CoolMessageData recoveredData{parser.getCurrentMessage()};

        EXPECT_EQ(recoveredData.type, Byte{0});
        EXPECT_EQ(recoveredData.data, kMessageData);
    }


}



