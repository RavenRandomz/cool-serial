#ifndef COOL_SERIAL_HEADER_SECTION_HPP
#define COOL_SERIAL_HEADER_SECTION_HPP

#include "bytes.hpp"
#include "cista_protocol.hpp"

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
        /**
         * Generates the serialized version of the HeaderSection
         * with complete crc.
         *
         * Because the serialized form of the dataInfo needs to be
         * serialized to calculate the crc8, it was decided
         * to generate the crc8
         */
        HeaderSection(const DataInfo& dataInfo)
        {
            Bytes serializedInfo{cista::serialize(dataInfo)};

            // HACK: automatic casting of unsigned char to uint8_t
            // seems to be well-defined, but if there are glitches, look here
            const uint8_t dataInfoCrc
            {
                CRC8::CRC8::calc(&serializedInfo[0], serializedInfo.size())
            };

            // Append Crc
            serializedInfo.push_back(dataInfoCrc);

            // Prevent overflow error by adding required Cista struct ending
            serializedInfo.push_back(kCistaStructEnding);

            serialized_ = std::move(serializedInfo);
        }

        /**
         * Returns serialized info 
         */
        const Bytes getSerialized() const
        {
            return serialized_;
        }
    private:
        Bytes serialized_;
    };
}
#endif
