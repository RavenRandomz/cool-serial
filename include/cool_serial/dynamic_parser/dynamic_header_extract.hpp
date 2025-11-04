#ifndef COOL_SERIAL_DYNAMIC_PARSER_DYNAMIC_HEADER_EXTRACT_HPP
#define COOL_SERIAL_DYNAMIC_PARSER_DYNAMIC_HEADER_EXTRACT_HPP

#include "cool_serial/byte_queue.hpp"
#include "cool_serial/header_section.hpp"
#include "cool_serial/dynamic_parser/header_found_listener.hpp"

namespace coolSerial
{
/**
 * This continuously checks a stream to make sure that the
 * header data can be found. A header for standard CoolSerial is 5 bytes.
 * This class can continuously check the queue until a potential header is
 * extracted.
 *
 * This header that is extracted may or may not be valid.
 *
 * WARNING: This causes the byteQueue to dump each byte. This is a destructive
 * operation
 */
class DynamicHeaderExtract
{
public:
    DynamicHeaderExtract(ByteQueue& queue, HeaderFoundListener& listener)
        :
        queue_{queue},
        listener_{listener}
    {}

    void update()
    {
        while(queue_.byteAvailable())
        {
            // This is to optimize performance
            // So instead of using a queue, it is directly written
            // to the array
            headerBytes_[headerByteIndex_] = queue_.getNextPoppedByte();

            if(headerByteIndex_ == kHeaderByteCountMaxIndex)
            {
                proccessHeaderBytes();
                resetHeaderExtraction();
                break;
            }
            // Prepare next interation for next index
            ++headerByteIndex_;
        }
    }

    void proccessHeaderBytes()
    {
        // Reset for next iteration
        const auto kHeaderData{HeaderSection::deserializeBytes(headerBytes_)};
        listener_.headerFound(kHeaderData);
    }

private:
    static const int kHeaderByteCountMaxIndex{3};

    ByteQueue& queue_; 
    HeaderFoundListener& listener_;

    int headerByteIndex_{0};
    HeaderBytes headerBytes_{};

    void resetHeaderExtraction()
    {
        headerByteIndex_ = 0;
    }
};
}
#endif
