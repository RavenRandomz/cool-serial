#ifndef COOL_SERIAL_I_DATA_HANDLER_HPP
#define COOL_SERIAL_I_DATA_HANDLER_HPP
#include "bytes.hpp"

namespace coolSerial
{
    class IDataHandler
    {
    public:
        /**
         * Trigger a handler's action
         */
        virtual void handleData(Bytes& bytes) = 0;
        virtual ~IDataHandler() {}
    };
}
#endif
