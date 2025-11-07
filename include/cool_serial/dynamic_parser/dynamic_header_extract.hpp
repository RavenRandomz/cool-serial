#ifndef COOL_SERIAL_DYNAMIC_PARSER_DYNAMIC_HEADER_EXTRACT_HPP
#define COOL_SERIAL_DYNAMIC_PARSER_DYNAMIC_HEADER_EXTRACT_HPP

#include "cool_serial/byte_queue.hpp"
#include "cool_serial/header_section.hpp"
#include "cool_serial/dynamic_parser/header_found_listener.hpp"

#include "cool_serial/dynamic_parser/dynamic_segment_extractor.hpp"
#include "cool_serial/dynamic_parser/segment_found_listener.hpp"
#include <cassert>

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
class DynamicHeaderExtract : public SegmentFoundListener
{
public:
    DynamicHeaderExtract(ByteQueue& queue, HeaderFoundListener& listener)
        :
        queue_{queue},
        listener_{listener}
    {}

    /**
     * The segmentExractor will call this function when it has located a segment.
     */
    void segmentFound(const Bytes& bytes) override
    {
        const HeaderBytes kHeaderBytes{
            bytes[0],
            bytes[1],
            bytes[2],
            bytes[3]
        };

        assert(bytes.size() == kHeaderBytes.size() && "Valid buffer required");

        // Reset for next iteration
        const auto kHeaderData{HeaderSection::deserializeBytes(kHeaderBytes)};
        listener_.headerFound(kHeaderData);
        assert(!segmentExtractor_.isMidExtraction() && "Segment extractor should be ready for new cycle");
    }

    void update()
    {
        segmentExtractor_.update();
    }

private:
    static const int kHeaderByteCountMaxIndex{3};

    ByteQueue& queue_; 
    HeaderFoundListener& listener_;

    DynamicSegmentExtractor segmentExtractor_{queue_, *this, kHeaderByteCountMaxIndex + 1};
};
}
#endif
