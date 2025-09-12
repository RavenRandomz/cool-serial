#ifndef COOL_SERIAL_COOL_MESSAGE_HPP
#define COOL_SERIAL_COOL_MESSAGE_HPP

#include "bytes.hpp"
#include "header_section.hpp"

#include <cassert>
#include <limits>
#include <iostream>
#include <vector>

namespace coolSerial
{
    /**
     * Contains data to be passed onto a handler 
     */
    struct CoolMessageData
    {
        Byte dataType;
        Bytes data;
    };

    /**
     * Preliminary class 
     *
     * Optimizations would be to create fixed-width templated message class which
     * could be used for each specific message type.
     *
     * For now, this will suffice for testing.
     */
    class CoolMessage
    {
    public:
        /// 137 is a cool approximation of a fundamental universal constant :D
        /// Obscura aside, this is always at the beginning of a message
        static constexpr Byte kStartOfFrame{137};

        /// A CoolSerial packet length is stored by 16 bits
        static constexpr size_t kMaxDataSize{std::numeric_limits<uint16_t>::max()};

        /**
         * This is more performant via move semantics :D
         *
         * Turns your data into a CoolSerial compliant
         * serial message!
         */
        CoolMessage(const Bytes&& data, Byte dataType)
        {
            assert(data.size() <= kMaxDataSize);
            // Add Start of Frame
            messageInfo_.push_back(kStartOfFrame);

            const HeaderSection kHeader{generateHeader(data, dataType)};
            const HeaderBytes kHeaderBytes{kHeader.getSerialized()}; 

            // Append header to beginning of transmission
            messageInfo_.insert(messageInfo_.end(), std::begin(kHeaderBytes), std::end(kHeaderBytes));

            // Append data
            messageInfo_.insert(messageInfo_.end(), std::begin(data), std::end(data));
        }

        /**
         * Returns CoolSerial frame in the form of bytes
         */
        Bytes getFrame() const
        {
            return messageInfo_;
        }

        friend std::ostream& operator<<(std::ostream& os, const CoolMessage& message);

    private:
        static HeaderSection generateHeader(const Bytes& data, Byte dataType)
        {
            const DataInfo kDataInfo
            {
                .dataLength = static_cast<uint16_t>(data.size()),
                .dataType = dataType 
            };

            const HeaderSection kHeader{kDataInfo};
            return kHeader;
        }

        Bytes messageInfo_{};
    };

    inline std::ostream& operator<<(std::ostream& os, const CoolMessage& message)
    {
        os << "[ ";
        for(auto& byte : message.messageInfo_)
        {
            os << std::hex << static_cast<int>(byte);
            os << ", ";
        }
        os << "]";
        return os;
    }
}
#endif
