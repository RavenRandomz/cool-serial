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
    constexpr int kDataInfoBytes{3}; // Refer to Cool Serial Protocol
                                       //
    using HeaderBytes = std::array<Byte, 4>;
    using DataInfoBytes = std::array<Byte, kDataInfoBytes>;

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

        DataInfoBytes serialize() const
        {

            const uint8_t firstLengthByte = (dataLength >> 8) & 0xFF;
            const uint8_t secondLengthByte = (dataLength >> 0) & 0xFF;

            // https://stackoverflow.com/questions/54255885/portable-way-of-splitting-n-byte-integer-into-single-bytes
            // HACK: Copy initialization is required to bypass restriction on narrowing conversion
            // Unsigned ints are well-defined across architectures
            // DO NOT REFACTOR: These constants will be optimized away
            DataInfoBytes bytes{firstLengthByte, secondLengthByte, dataType};
            return bytes; 
        }

        Byte getCrc() const
        {
            const DataInfoBytes kSerialized{serialize()};
            return CRC8::CRC8::calc(kSerialized.data(), kSerialized.size());
        }
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

        /**
         * Check if the CRC8 byte matches the
         * Data Info subsection CRC8
         */
        bool isValid() const
        {
            return dataInfo.getCrc() == dataInfoCrc;
        }
    };

    /**
     * This stores the serialized portion of the header section.
     */
    class HeaderSection
    {
    public:
        /**
         * Generate the serial information
         */
        static HeaderBytes generateSerialized(const DataInfo& dataInfo)
        {
            HeaderBytes bytes{};
            const DataInfoBytes kSerializedDataInfo{dataInfo.serialize()};

            //Access indicies 0 - 2
            for(int i{0} ; i < kDataInfoBytes; ++i )
            {
                bytes[i] = kSerializedDataInfo[i];
            }

            bytes[3] = dataInfo.getCrc();
            return bytes;
        }

        /**
         * Generate the serial information
         */
        static HeaderData deserializeBytes(const HeaderBytes& bytes)
        {

            // https://stackoverflow.com/questions/54255885/portable-way-of-splitting-n-byte-integer-into-single-bytes
            // HACK: Copy initialization is required to bypass restriction on narrowing conversion
            // Unsigned ints are well-defined across architectures
            // DO NOT REFACTOR: These constants will be optimized away
            //

            const uint16_t kDataLength{static_cast<uint16_t>(bytes[0] << 8 | bytes[1] << 0)};
            const uint8_t kDataType{bytes[2]}; 
            const Byte kCrc{bytes[3]};

            return HeaderData
            {
                .dataInfo = 
                {
                    .dataLength = kDataLength,
                    .dataType = kDataType,
                },
                .dataInfoCrc = kCrc
            };
        }

        /**
         * Generates the serialized version of the HeaderSection
         * with complete crc.
         *
         * Because the serialized form of the dataInfo needs to be
         * serialized to calculate the crc8, it was decided
         * to generate the crc8
         */
        HeaderSection(const DataInfo& dataInfo)
            :
            serialized_{generateSerialized(dataInfo)}
        {
        }

        /**
         * Returns serialized info 
         */
        const HeaderBytes getSerialized() const
        {
            return serialized_;
        }

        
    private:
        HeaderBytes serialized_;
    };
}
#endif
