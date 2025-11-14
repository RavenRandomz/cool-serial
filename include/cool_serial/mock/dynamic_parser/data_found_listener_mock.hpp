#ifndef COOL_SERIAL_MOCK_DYNAMIC_PARSER_DATA_FOUND_LISTENER_MOCK_HPP
#define COOL_SERIAL_MOCK_DYNAMIC_PARSER_DATA_FOUND_LISTENER_MOCK_HPP
#include "cool_serial/dynamic_parser/data_found_listener.hpp"

#include <gmock/gmock.h>
namespace coolSerial
{
class DataFoundListenerMock : public DataFoundListener
{
public:
    MOCK_METHOD(void, dataFound, (const CoolMessageData& data), (override));
};
}
#endif
