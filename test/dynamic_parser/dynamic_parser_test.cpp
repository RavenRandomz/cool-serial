#include <gtest/gtest.h>

#include "cool_serial/byte_queue.hpp"
#include "cool_serial/dynamic_parser/dynamic_parser.hpp"
#include "cool_serial/mock/dynamic_parser/data_found_listener_mock.hpp"
#include "cool_serial/bytes.hpp"
#include "cool_serial/cool_message.hpp"
#include "cista/cista.h"

#include <iostream>

// Extra namespace prevents test structs and fixtures from causing collisions
// CoolSerial v1.0.0 defines the first 5 bytes then the next bytes are data
namespace coolSerial::dynamicParserTest
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

TEST(DynamicParser, Deserialization)
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

    DataFoundListenerMock dataListener{};

    // Test the actual parser
    DynamicParser parser{testQueue,dataListener};

    const CoolMessageData kExpected{0, kMessageData};
    EXPECT_CALL(dataListener, dataFound(kExpected));
    parser.update();
}

TEST(DynamicParser, continousDeserilization)
{
    // A partial buffer will be fed each time
    // [ 89, 0, 10, 0, 57, db, f9, e, 40, 27, a0, a8, c8, 1b, f5, 10, d, eb, dc, e0, 40, ]
    // Partial header with SOF
    Bytes partialHeader{ 0x84, 0x0, 0x10};
    Bytes restOfHeader{ 0x0, 0x57};
    Bytes halfOfData{ 0xdb, 0xf9, 0xe, 0x40};
    Bytes finalHalfOfData{0x27, 0xa0, 0xa8, 0xc8, 0x1b, 0xf5, 0x10, 0xd, 0xeb, 0xdc, 0xe0, 0x40};

    const Bytes kFullData{0xdb, 0xf9, 0xe, 0x40, 0x27, 0xa0, 0xa8, 0xc8, 0x1b, 0xf5, 0x10, 0xd, 0xeb, 0xdc, 0xe0, 0x40};

    EXPECT_EQ(kFullData.size(), halfOfData.size() + finalHalfOfData.size());

    ByteQueue testQueue{};

    DataFoundListenerMock dataListener{};
    const CoolMessageData kExpected{0, kFullData};
    EXPECT_CALL(dataListener, dataFound(kExpected));

    // Test the actual parser
    DynamicParser parser{testQueue, dataListener};

    addBytesToByteQueue(partialHeader, testQueue);
    parser.update();

    addBytesToByteQueue(restOfHeader, testQueue);
    parser.update();

    addBytesToByteQueue(halfOfData, testQueue);
    parser.update();

    addBytesToByteQueue(finalHalfOfData, testQueue);
    parser.update();

    EXPECT_TRUE(testQueue.empty());



}

TEST(DynamicParser, ReallyFragmentedDeserialization)
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
    DataFoundListenerMock dataListener{};
    DynamicParser parser{simulatedQueue, dataListener};

    const CoolMessageData kExpected{0, kMessageData};
    EXPECT_CALL(dataListener, dataFound(kExpected));

    while(!testQueue.empty())
    {
        simulatedQueue.push(testQueue.front());
        testQueue.pop();
        parser.update();
    }
}

/**
 * Check that red herrring data is not extracted
 */
TEST(DynamicParser, extraData)
{
    // A partial buffer will be fed each time
    // [ 89, 0, 10, 0, 57, db, f9, e, 40, 27, a0, a8, c8, 1b, f5, 10, d, eb, dc, e0, 40, ]
    // Partial header with SOF
    Bytes partialHeader{ 0x84, 0x0, 0x10};
    Bytes restOfHeader{ 0x0, 0x57};
    // Everythin after 0x40 is a red herring
    Bytes halfOfData{ 0xdb, 0xf9, 0xe, 0x40};
    Bytes redHerring{ 0x84,0xdb, 0xf9, 0xe, 0x40};
    Bytes finalHalfOfData{0x27, 0xa0, 0xa8, 0xc8, 0x1b, 0xf5, 0x10, 0xd, 0xeb, 0xdc, 0xe0, 0x40};

    const Bytes kFullData{0xdb, 0xf9, 0xe, 0x40, 0x27, 0xa0, 0xa8, 0xc8, 0x1b, 0xf5, 0x10, 0xd, 0xeb, 0xdc, 0xe0, 0x40};

    ByteQueue testQueue{};

    DataFoundListenerMock dataListener{};
    const CoolMessageData kExpected{0, kFullData};
    EXPECT_CALL(dataListener, dataFound(kExpected));

    // Test the actual parser
    DynamicParser parser{testQueue, dataListener};

    addBytesToByteQueue(partialHeader, testQueue);
    parser.update();

    addBytesToByteQueue(restOfHeader, testQueue);
    parser.update();

    addBytesToByteQueue(halfOfData, testQueue);
    parser.update();

    addBytesToByteQueue(finalHalfOfData, testQueue);
    addBytesToByteQueue(redHerring, testQueue);
    parser.update();

    // The queue should have th red herring data
    EXPECT_FALSE(testQueue.empty());
}

/**
 *
 * Check that one and only one message is extracted as well as the extraction of multiple
 * messages
 */
TEST(DynamicParser, MultipleMessage)
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
    DataFoundListenerMock dataListener{};
    DynamicParser parser{simulatedQueue, dataListener};

    const CoolMessageData kExpected{0, kMessageData};
    EXPECT_CALL(dataListener, dataFound(kExpected));

    while(!testQueue.empty())
    {
        simulatedQueue.push(testQueue.front());
        testQueue.pop();
        parser.update();
    }

    TestStruct data2
    {
        .x = 5,
        .y = -345.34,
        .z = 99.39
    };

    Bytes kDataBytes2{cista::serialize(data)};

    // The 3rd header section byte (index 3) must equal 0
    CoolMessage message2{std::move(kDataBytes2), Byte{0}};

    Bytes kMessageFrame2{message.getFrame()};

    const Bytes kMessageData2{kMessageFrame.begin() + 5, kMessageFrame.end()};
    const CoolMessageData kExpected2{0, kMessageData2};

    EXPECT_CALL(dataListener, dataFound(kExpected2));

    ByteQueue testQueue2{std::deque<Byte>{kMessageFrame2.begin(), kMessageFrame2.end()}};

    while(!testQueue2.empty())
    {
        simulatedQueue.push(testQueue2.front());
        testQueue2.pop();
        parser.update();
    }
}
}

