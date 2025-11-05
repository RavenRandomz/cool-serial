#ifndef COOL_SERIAL_BYTE_QUEUE
#define COOL_SERIAL_BYTE_QUEUE

#include "byte.hpp"
#include "bytes.hpp"

#include <queue>

namespace coolSerial
{
    class ByteQueue : public std::queue<Byte>
    {
    public:
        using std::queue<Byte>::queue;

        bool byteAvailable() const
        {
            return !empty();
        }

        Byte getNextPoppedByte()
        {
            const Byte kNext{front()};
            pop();
            return kNext;
        }

        void addBytes(const Bytes& bytes)
        {
            for(const Byte& kByte : bytes)
            {
                push(kByte);
            }
        }
    };
}
#endif
