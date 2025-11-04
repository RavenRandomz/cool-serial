#ifndef COOL_SERIAL_DYNAMIC_PARSER_HEADER_FOUND_LISTENER_MOCK_HPP
#define COOL_SERIAL_DYNAMIC_PARSER_HEADER_FOUND_LISTENER_MOCK_HPP

#include "cool_serial/dynamic_parser/header_found_listener.hpp"
#include <gmock/gmock.h>

namespace coolSerial
{
class HeaderFoundListenerMock : public HeaderFoundListener
{
public:
    
    MOCK_METHOD(void, headerFound, (const HeaderSection& header), (override));
};
}

#endif

