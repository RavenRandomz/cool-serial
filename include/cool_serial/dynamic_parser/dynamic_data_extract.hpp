#ifndef COOL_SERIAL_DYNAMIC_PARSER_DYNAMIC_DATA_EXTRACT_HPP
#define COOL_SERIAL_DYNAMIC_PARSER_DYNAMIC_DATA_EXTRACT_HPP

#include "cool_serial/byte_queue.hpp"
#include "cool_serial/dynamic_parser/data_found_listener.hpp"
#include "cool_serial/header_section.hpp"

#include "cool_serial/dynamic_parser/dynamic_segment_extractor.hpp"
#include "cool_serial/dynamic_parser/segment_found_listener.hpp"

#include <cinttypes>

namespace coolSerial
{
/**
 * This continuously checks a stream to make sure that the
 * retquired serialized data is found.
 *
 * This informs a listener whenever the data is found.
 *
 * WARNING: This causes the byteQueue to dump each byte. This is a destructive
 * operation
 *
 * WARNING: You will need to synchronize the data extraction by locating
 * a valid header (data section follows the header).
 */
class DynamicDataExtract : public SegmentFoundListener
{
public:
    DynamicDataExtract(ByteQueue& queue, DataFoundListener& listener)
        :
        queue_{queue},
        listener_{listener}
    {}


    void segmentFound(const Bytes& bytes) override
    {
        // Consider keeping object alive to prevent allocation reallocation
        // profile first before optimizing
        const CoolMessageData kExtractedData{extractedDataType_, bytes};
        listener_.dataFound(kExtractedData);
    }

    void update()
    {
        segmentExtractor_.update();
    }

    /**
     * Set what data type apend to extracted data metadata
     * and length the data segment to extract is
     * (what a mouthfull)
     */
    void setDataExtractionInfo(const DataInfo& dataInfo)
    {
        extractedDataType_ = dataInfo.dataType;
        segmentExtractor_.setSegmentSize(dataInfo.dataLength);
    }

private:
    ByteQueue& queue_; 
    DataFoundListener& listener_;

    Byte extractedDataType_{0}; 
    DynamicSegmentExtractor segmentExtractor_{queue_, *this};
};
}
#endif
