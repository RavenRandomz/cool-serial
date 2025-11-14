#include "cool_serial/byte_queue.hpp"
#include "cool_serial/dynamic_parser/dynamic_sof_search.hpp"
#include "cool_serial/mock/dynamic_parser/start_of_frame_listener_mock.hpp"

#include <gtest/gtest.h>

namespace coolSerial::dynamicSofSearchTest
{
/**
 * The test byte queue will NOT contain the SOF 0x84
 */
TEST(DynamicSofSearch, noSof)
{
    // This is taken from the original continuousParser test since it has a known precalculated value
    Bytes finalHalfOfData{0x27, 0xa0, 0xa8, 0xc8, 0x1b, 0xf5, 0x10, 0xd, 0xeb, 0xdc, 0xe0, 0x40};
    ByteQueue testQueue{std::deque<Byte>{finalHalfOfData.begin(), finalHalfOfData.end()}};

    // Fail if listening function is called
    testing::StrictMock<StartOfFrameFoundListenerMock> listener{};

    DynamicSofSearch sofSearcher{testQueue, listener};
    sofSearcher.update();
}

/**
 * The test byte queue will contain the SOF 0x84
 */
TEST(DynamicSofSearch, sofBegin)
{
    // This is adapt from the original continuousParser test since it has a known precalculated value
    Bytes finalHalfOfData{0x84,0x27, 0xa0, 0xa8, 0xc8, 0x1b, 0xf5, 0x10, 0xd, 0xeb, 0xdc, 0xe0, 0x40};
    ByteQueue testQueue{std::deque<Byte>{finalHalfOfData.begin(), finalHalfOfData.end()}};

    // Fail if listening function is called
    testing::StrictMock<StartOfFrameFoundListenerMock> listener{};

    EXPECT_CALL(listener, startOfFrameFound());
    DynamicSofSearch sofSearcher{testQueue, listener};
    sofSearcher.update();
}

/**
 * The test byte queue will contain the SOF 0x84
 */
TEST(DynamicSofSearch, sofMid)
{
    // This is adapted from the original continuousParser test since it has a known precalculated value
    Bytes finalHalfOfData{0x27, 0xa0, 0xa8, 0x84, 0x1b, 0xf5, 0x10, 0xd, 0xeb, 0xdc, 0xe0, 0x40};
    ByteQueue testQueue{std::deque<Byte>{finalHalfOfData.begin(), finalHalfOfData.end()}};

    // Fail if listening function is called
    testing::StrictMock<StartOfFrameFoundListenerMock> listener{};

    EXPECT_CALL(listener, startOfFrameFound());
    DynamicSofSearch sofSearcher{testQueue, listener};
    sofSearcher.update();
}

/**
 * The test byte queue will contain the SOF 0x84
 */
TEST(DynamicSofSearch, sofEnd)
{
    // This is adapted from the original continuousParser test since it has a known precalculated value
    Bytes finalHalfOfData{0x27, 0xa0, 0xa8, 0x82, 0x1b, 0xf5, 0x10, 0xd, 0xeb, 0xdc, 0xe0, 0x84};
    ByteQueue testQueue{std::deque<Byte>{finalHalfOfData.begin(), finalHalfOfData.end()}};

    // Fail if listening function is called
    testing::StrictMock<StartOfFrameFoundListenerMock> listener{};

    EXPECT_CALL(listener, startOfFrameFound());
    DynamicSofSearch sofSearcher{testQueue, listener};
    sofSearcher.update();
}
}
