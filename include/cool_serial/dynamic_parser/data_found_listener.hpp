#ifndef COOL_SERIAL_DYNAMIC_PARSER_DATA_FOUND_LISTENER_HPP
#define COOL_SERIAL_DYNAMIC_PARSER_DATA_FOUND_LISTENER_HPP

#include "cool_serial/cool_message.hpp"

namespace coolSerial
{
    /**
     * Whenever a byte segment is extracted
     * a class that implements this can be passed the data
     */
    class DataFoundListener 
    {
    public:
        virtual void dataFound(const CoolMessageData& data) = 0;
        virtual ~DataFoundListener() {}
    };
}
#endif
