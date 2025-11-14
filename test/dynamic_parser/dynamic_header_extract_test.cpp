#include <cool_serial/header_section.hpp>
#include <deque>
#include <gtest/gtest.h>

#include "cool_serial/dynamic_parser/dynamic_header_extract.hpp"
#include "cool_serial/mock/dynamic_parser/header_found_listener_mock.hpp"


namespace coolSerial::dynamicHeaderExtractTest
{
/**
 * Test to ensure no false positives
 */
TEST(DynamicHeaderExtract, emptyQueue)
{
    ByteQueue testQueue{};
    // There should be 
    testing::StrictMock<HeaderFoundListenerMock> listener{};
    DynamicHeaderExtract extracter{testQueue, listener};

}

/**
 * Generate a header, attempt to extract proper information
 */
TEST(DynamicHeaderExtract, validHeader)
{
    DataInfo kDataInfo
    {
        .dataLength = 10,
        .dataType = 10
    };
    const HeaderSection kTestHeader{kDataInfo};
    const HeaderData kExpectedHeaderData{kTestHeader.getHeaderData()};

    const auto kHeaderBytes{kTestHeader.getSerialized()};

    ByteQueue testQueue{std::deque<Byte>{kHeaderBytes.begin(), kHeaderBytes.end()}};

    testing::StrictMock<HeaderFoundListenerMock> listener{};
    DynamicHeaderExtract extracter{testQueue, listener};

    EXPECT_CALL(listener, headerFound(kExpectedHeaderData));
    extracter.update();
}

/**
 * Generate a header, attempt to extract proper information in multiple steps
 */
TEST(DynamicHeaderExtract, partialHeader)
{
    DataInfo kDataInfo
    {
        .dataLength = 10,
        .dataType = 0
    };

    const HeaderSection kTestHeader{kDataInfo};
    const HeaderData kExpectedHeaderData{kTestHeader.getHeaderData()};


    const auto kHeaderBytes{kTestHeader.getSerialized()};
    // In this test, the header bytes will be spit into two groups so that
    // the ability to extract a header from two incomplete buffer states
    // without a false positive is tested

    const int kSliceIndex{1};

    // NOTE THAT THE LAST is excluded (kHeaderBytes.begin() + kSliceIndex) is not included)
    ByteQueue testQueue{std::deque<Byte>{kHeaderBytes.begin(), kHeaderBytes.begin() + kSliceIndex}};

    testing::StrictMock<HeaderFoundListenerMock> listener{};
    DynamicHeaderExtract extracter{testQueue, listener};
    extracter.update();

    const Bytes kOtherBytes{kHeaderBytes.begin() + kSliceIndex, kHeaderBytes.end()};
    testQueue.addBytes(kOtherBytes);

    // Only expect call from this point
    EXPECT_CALL(listener, headerFound(kExpectedHeaderData));
    extracter.update();
}
}
