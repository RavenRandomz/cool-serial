#ifndef COOL_SERIAL_DYNAMIC_PARSER_HEADER_FOUND_LISTENER_HPP
#define COOL_SERIAL_DYNAMIC_PARSER_HEADER_FOUND_LISTENER_HPP

#include "cool_serial/header_section.hpp"

namespace coolSerial
{
    class HeaderFoundListener
    {
    public:
        virtual void headerFound(const HeaderSection& header) = 0;
        virtual ~HeaderFoundListener() {}
    };
}
#endif
