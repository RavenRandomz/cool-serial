#include "cool_serial/byte_queue.hpp"
#include "cool_serial/dynamic_parser/dynamic_data_extract.hpp"
#include "cool_serial/mock/dynamic_parser/data_found_listener_mock.hpp"

#include "gmock/gmock.h"
#include <cool_serial/cool_message.hpp>
#include <gtest/gtest.h>

namespace coolSerial::dynamicDataExtractTest
{
/**
 * Test for no false positives regarding empty queue
 */
TEST(DynamicDataExtract, emptyQueueNoFalsePositive)
{
    ByteQueue testQueue{};

    // Fail if listening function is called
    testing::StrictMock<DataFoundListenerMock> listener{};

    const int kDataSize{5};
    const Byte kDataType{0};
    const DataInfo kDataInfo{kDataSize,kDataType};

    DynamicDataExtract extractor{testQueue, listener};
    extractor.setDataExtractionInfo(kDataInfo);
    extractor.update();
    // The listener should not be called
}

/**
 * Test for successful extraction with a buffer that has an equal amount to the desired
 * size
 */
TEST(DynamicDataExtract, basicExtraction)
{
    const Bytes kTestBytes{1,2,3,4,5};
    const std::uint16_t kDataSize{static_cast<uint16_t>(kTestBytes.size())};
    ByteQueue queue{};
    queue.addBytes(kTestBytes);

    const Byte kDataType{0};
    const DataInfo kDataInfo{kDataSize,kDataType};
    testing::StrictMock<DataFoundListenerMock> listener{};

    DynamicDataExtract extractor{queue, listener};
    extractor.setDataExtractionInfo(kDataInfo);

    const CoolMessageData kExpectedMessageData{kDataType, kTestBytes};
    EXPECT_CALL(listener, dataFound(kExpectedMessageData));
    extractor.update();
}

/**
 * Test for successful extraction with a buffer that has an equal amount to the desired
 * size
 */
TEST(DynamicDataExtract, extraExtraction)
{
    const Bytes kTestBytes{1,2,3,4,5};
    const Bytes kTestQueue{1,2,3,4,5, 2, 4, 5};
    const std::uint16_t kDataSize{static_cast<uint16_t>(kTestBytes.size())};
    ByteQueue queue{};
    queue.addBytes(kTestQueue);

    const Byte kDataType{0};
    const DataInfo kDataInfo{kDataSize,kDataType};
    testing::StrictMock<DataFoundListenerMock> listener{};

    DynamicDataExtract extractor{queue, listener};
    extractor.setDataExtractionInfo(kDataInfo);

    const CoolMessageData kExpectedMessageData{kDataType, kTestBytes};
    EXPECT_CALL(listener, dataFound(kExpectedMessageData));
    extractor.update();
}

/**
 * Test for successful extraction with a buffer that has an equal amount to the desired
 * size
 */
TEST(DynamicDataExtract, partialExtraction)
{
    const Bytes kTestBytes{1,2,3,4,5};
    const Bytes kTestBytes1{1,2};
    const Bytes kTestBytes2{3,4,5};
    const std::uint16_t kDataSize{static_cast<uint16_t>(kTestBytes.size())};
    ByteQueue queue{};
    queue.addBytes(kTestBytes1);

    const Byte kDataType{0};
    const DataInfo kDataInfo{kDataSize,kDataType};
    testing::StrictMock<DataFoundListenerMock> listener{};

    DynamicDataExtract extractor{queue, listener};
    extractor.setDataExtractionInfo(kDataInfo);

    const CoolMessageData kExpectedMessageData{kDataType, kTestBytes};
    extractor.update();

    EXPECT_CALL(listener, dataFound(kExpectedMessageData));
    queue.addBytes(kTestBytes2);
    extractor.update();
}
}
