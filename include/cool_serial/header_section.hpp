#ifndef COOL_SERIAL_HEADER_SECTION_HPP
#define COOL_SERIAL_HEADER_SECTION_HPP

#include "bytes.hpp"
#include "cista/cista.h"
#include "cppcrc/cppcrc.h"

#include <cstdint>

namespace coolSerial
{
    constexpr Byte kStartOfFrame{137}; // Approximate fine structure constant :D

    /**
     * This is the portion of the header that determines how the data section
     * is interpreted.
     *
     * @param dataLength - how many bytes make up the data section
     * @param dataType - how the data should be interpreted (user defined)
     */
    struct DataInfo
    {
        std::uint16_t dataLength;
        std::uint8_t dataType;
    };

    /**
     * Defines the header structure. The SOF byte is not considered part of the
     * header section.
     *
     * @param dataInfo - information about the data section
     * @param headerCrc - CRC8 for dataInfo (minus the CRC8)
     *
     * This is a struct for easier serialization and deserialization
     */
    struct HeaderData
    {
        DataInfo dataInfo;
        Byte dataInfoCrc; 
    };

    class HeaderSection
    {
    public:
        HeaderSection(const DataInfo& dataInfo)
        {
            Bytes serializedInfo{cista::serialize(dataInfo)};
            const uint8_t dataInfoCrc
            {
                CRC8::CRC8::calc(&serializedInfo[0], serializedInfo.size())
            };

            serializedInfo.push_back(dataInfoCrc);
            serialized_ = std::move(serializedInfo);
        }
    private:
        Bytes serialized_;
    };
};
#endif
