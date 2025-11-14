#ifndef COOL_SERIAL_DYNAMIC_PARSER_SEGMENT_FOUND_LISTENER_MOCK_HPP
#define COOL_SERIAL_DYNAMIC_PARSER_SEGMENT_FOUND_LISTENER_MOCK_HPP

#include "cool_serial/dynamic_parser/segment_found_listener.hpp"
#include <gmock/gmock.h>

namespace coolSerial
{
class SegmentFoundListenerMock : public SegmentFoundListener
{
public:
    MOCK_METHOD(void, segmentFound, (const Bytes& segment), (override));
};
}

#endif

