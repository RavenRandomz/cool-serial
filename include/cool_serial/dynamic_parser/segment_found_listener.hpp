#ifndef COOL_SERIAL_DYNAMIC_PARSER_SEGMENT_FOUND_LISTENER_HPP
#define COOL_SERIAL_DYNAMIC_PARSER_SEGMENT_FOUND_LISTENER_HPP

#include "cool_serial/bytes.hpp"

namespace coolSerial
{
    /**
     * Whenever a byte segment is extracted
     * a class that implements this can be passed the data
     */
    class SegmentFoundListener 
    {
    public:
        virtual void segmentFound(const Bytes& bytes) = 0;
        virtual ~SegmentFoundListener() {}
    };
}
#endif
