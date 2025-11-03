#ifndef COOL_SERIAL_BYTE_QUEUE
#define COOL_SERIAL_BYTE_QUEUE

#include "byte.hpp"
#include <queue>

namespace coolSerial
{
    class ByteQueue : public std::queue<Byte>
    {
    public:
        bool byteAvailable()
        {
            return !empty();
        }
    };
}
#endif
