#ifndef COOL_SERIAL_HEADER_SECTION_HPP
#define COOL_SERIAL_HEADER_SECTION_HPP

#include "byte.hpp"
#include <cstdint>

namespace coolSerial
{
    constexpr Byte kStartOfFrame{137}; // Approximate fine structure constant :D
    
    /**
     * Defines the header structure. The SOF byte is not considered part of the
     * header section.
     *
     * @param dataLength - how many bytes make up the data section
     * @param messageType - how the data should be interpreted (user defined)
     * @param headerCrc - CRC8 for the header section (minus the CRC8)
     */
    struct HeaderSection
    {
        std::uint16_t dataLength;
        std::uint8_t messageType;
        Byte headerCrc; 
    };

};
#endif
