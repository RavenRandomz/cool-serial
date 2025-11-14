
#ifndef COOL_SERIAL_DATA_HANDLER_I_DATA_HANDLER_MOCK_HPP
#define COOL_SERIAL_DATA_HANDLER_I_DATA_HANDLER_MOCK_HPP
#include "cool_serial/data_handler/idata_handler.hpp"
#include <gmock/gmock.h>

namespace coolSerial
{
    class IDataHandlerMock : public IDataHandler
    {
        MOCK_METHOD(void, handleData, (const Bytes& bytes), (override));
    };
}
#endif

