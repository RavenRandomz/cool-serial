#ifndef COOL_SERIAL_BYTE_QUEUE
#define COOL_SERIAL_BYTE_QUEUE

#include "byte.hpp"
#include <queue>

namespace coolSerial
{
    using ByteQueue = std::queue<Byte>;
}
#endif
