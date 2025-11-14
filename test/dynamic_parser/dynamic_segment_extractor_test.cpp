#include "cool_serial/byte_queue.hpp"
#include "cool_serial/dynamic_parser/dynamic_segment_extractor.hpp"
#include "cool_serial/mock/dynamic_parser/segment_found_listener_mock.hpp"

#include "gmock/gmock.h"
#include <gtest/gtest.h>

namespace coolSerial::dynamicSofSearchTest
{
/**
 * Test for no false positives regarding empty queue
 */
TEST(DynamicSegmentExtractor, emptyQueueNoFalsePositive)
{
    ByteQueue testQueue{};

    // Fail if listening function is called
    testing::StrictMock<SegmentFoundListenerMock> listener{};

    const int kSegmentSize{5};

    DynamicSegmentExtractor extractor{testQueue, listener, kSegmentSize};
    extractor.update();
    // The listener should not be called
}

/**
 * Test for successful extraction with a buffer that has an equal amount to the desired
 * size
 */
TEST(DynamicSegmentExtractor, basicExtraction)
{
    const Bytes kTestBytes{1,2,3,4,5};
    const int kSegmentSize{static_cast<int>(kTestBytes.size())};
    ByteQueue queue{};
    queue.addBytes(kTestBytes);

    testing::StrictMock<SegmentFoundListenerMock> listener{};

    DynamicSegmentExtractor extractor{queue, listener, kSegmentSize};

    EXPECT_CALL(listener, segmentFound(kTestBytes));
    extractor.update();
}

/**
 * Test for successful extraction with a buffer that has an equal amount to the desired
 * size but ignoring extra bytes
 */
TEST(DynamicSegmentExtractor, extraByteIgnore)
{
    const Bytes kTestBytes{1,2,3,4,5,6, 7, 8};
    const Bytes kExpectedBytes{1,2,3,4};
    const int kSegmentSize{static_cast<int>(kExpectedBytes.size())};
    ByteQueue queue{};
    queue.addBytes(kTestBytes);

    testing::StrictMock<SegmentFoundListenerMock> listener{};

    DynamicSegmentExtractor extractor{queue, listener, kSegmentSize};

    EXPECT_CALL(listener, segmentFound(kExpectedBytes));
    extractor.update();
}
}
