#ifndef COOL_SERIAL_CISTA_PROTOCOL_HPP
#define COOL_SERIAL_CISTA_PROTOCOL_HPP
#include "byte.hpp"
namespace coolSerial
{
    /// Cista appends an extra byte after serializing a struct
    /// This is the standard serialization, not the hashed one
    constexpr Byte kCistaStructEnding {0x00};
}
#endif

