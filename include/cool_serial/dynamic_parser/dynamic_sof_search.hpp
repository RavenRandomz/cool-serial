#ifndef COOL_SERIAL_DYNAMIC_PASER_DYNAMIC_SOF_SEARCH_HPP
#define COOL_SERIAL_DYNAMIC_PASER_DYNAMIC_SOF_SEARCH_HPP

#include "cool_serial/byte_queue.hpp"
#include "cool_serial/cool_message.hpp"
#include "cool_serial/dynamic_parser/start_of_frame_found_listener.hpp"

namespace coolSerial
{
    /**
     * This can continuously search for the coolSerial SOF
     * Upon finding the SOF, it informs a single listener.
     *
     * (SOF: Start of Frame) Look up serial protocols such as RefSerial
     * or CoolSerial, etc. for more info
     *
     * Every time update() is called, it checks all available
     * bytes in the queue until the SOF is located. No more bytes
     * will be popped otherwise.
     *
     * Due to performance reasons, there is only one listener.
     * A listener adapter can be informed which would them pass
     * the call towards a collection of listeners if need be.
     *
     * WARNIG: This will dump the non-SOF bytes in the ByteQueue when searching
     * for the SOF byte.
     */
    class DynamicSofSearch
    {
    public:
        DynamicSofSearch(
            ByteQueue& queue,
            StartOfFrameFoundListener& listener
        ):
            queue_{queue},
            frameFoundListener_{listener}
            
        {
        }

        void update()
        {
            //Consider having an upper limit
            while(queue_.byteAvailable())
            {
                sofListenerInformCheck();
            }
        }
    private:
        ByteQueue& queue_;
        StartOfFrameFoundListener& frameFoundListener_;

        void sofListenerInformCheck()
        {
            const Byte kNext{queue_.getNextPoppedByte()};
            if(kNext == CoolMessage::kStartOfFrame)
            {
                frameFoundListener_.startOfFrameFound();
            }
        }
    };
}
#endif
