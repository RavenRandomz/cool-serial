#include "cool_serial/bytes.hpp"

namespace coolSerial
{
/**
     * This class takes a series of bytes, it is assigned the bytes.
     *
     * The implementer should correspond to a specific message code.
     *
     * CoolSerial only carries the data bytes. In order to decode them
     * another protocol must be established.
     *
     * CoolSerial supports up to 256 different protocols.
     */
class IDataHandler
{
public:
    virtual void handleData(const Bytes& bytes) = 0;
    virtual ~IDataHandler() {}
};
}
