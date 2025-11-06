#ifndef COOL_SERIAL_DYNAMIC_PARSER_DYNAMIC_SEGMENT_EXTRACTOR
#define COOL_SERIAL_DYNAMIC_PARSER_DYNAMIC_SEGMENT_EXTRACTOR

#include "cool_serial/byte_queue.hpp"
#include "cool_serial/dynamic_parser/segment_found_listener.hpp"
#include "cassert"

namespace coolSerial
{
/**
 * This continuously checks a stream gathering data until
 * it has obtained a number of bytes equal to the segmentSize
 *
 * update() must be repeatedly called.
 *
 * WARNING: This causes the byteQueue to dump each byte. This is a destructive
 * operation
 */
class DynamicSegmentExtractor
{
public:
    DynamicSegmentExtractor(ByteQueue& queue, SegmentFoundListener& listener, int segmentSize)
        :
        queue_{queue},
        listener_{listener}
    {
        setSegmentSize(segmentSize);
    }

    /**
     * This will take all of the bytes within the queue
     * and attempt to extract a segment of the specified size.
     *
     * If the segment is found, the queue will not be completely emptied.
     */
    void update()
    {
        // The gist is that it sequentially takes
        // a byte into the segment buffer until
        // the segment buffer is full.
        //
        // The listener will be have their segmentFound function
        // called which will pass the information.
        while(queue_.byteAvailable())
        {
            // This is to optimize performance
            // So instead of using a queue, it is directly written
            // to the array
            segmentBytes_[extractedIndex_] = queue_.getNextPoppedByte();

            if(extractedIndex_ == maxSegmentIndex_)
            {
                listener_.segmentFound(segmentBytes_);
                resetExtraction();
                break;
            }
            // Prepare next interation for next index
            ++extractedIndex_;
        }
    }

    /**
     * WARNING: setting segment size when
     * midExtraction is true can lead
     * to unexpected results
     */
    void setSegmentSize(int size)
    {
        assert(size >= 1 && "Invalid segment size");
        maxSegmentIndex_ = size - 1;
    }

    bool isMidExtraction()
    {
        return extractedIndex_ != 0;
    }

    /**
     * The current message buffer will be overwritten
     * with the next extraction.
     */
    void resetExtraction()
    {
        extractedIndex_ = 0;
    }

private:
    ByteQueue& queue_; 
    SegmentFoundListener& listener_;
    int maxSegmentIndex_{0};

    Bytes segmentBytes_{};
    int extractedIndex_{0};
};
}
#endif
