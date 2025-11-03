#ifndef COOL_SERIAL_DYNAMIC_PARSER_START_OF_FRAME_FOUND_LISTENER_MOCK_HPP
#define COOL_SERIAL_DYNAMIC_PARSER_START_OF_FRAME_FOUND_LISTENER_MOCK_HPP

#include "cool_serial/dynamic_parser/start_of_frame_found_listener.hpp"
#include <gmock/gmock.h>

namespace coolSerial
{
    class StartOfFrameFoundListenerMock : public StartOfFrameFoundListener
    {
        MOCK_METHOD(void, startofFrameFound, (), (override));
    };
}
#endif
