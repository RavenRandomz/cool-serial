#ifndef COOL_SERIAL_DYNAMIC_PARSER_DYNAMIC_PARSER_HPP
#define COOL_SERIAL_DYNAMIC_PARSER_DYNAMIC_PARSER_HPP
#include "cool_serial/byte_queue.hpp"

#include "cool_serial/dynamic_parser/dynamic_sof_search.hpp"
#include "cool_serial/dynamic_parser/dynamic_data_extract.hpp"
#include "cool_serial/dynamic_parser/dynamic_header_extract.hpp"
#include "cool_serial/dynamic_parser/data_found_listener.hpp"
#include <cool_serial/dynamic_parser/segment_found_listener.hpp>
#include <cool_serial/dynamic_parser/start_of_frame_found_listener.hpp>

namespace coolSerial
{
/**
 * This is able to parse an incomplete
 * stream of bytes (varying chunks)
 * When a complete Cool Serial message
 * has been recovered.
 *
 * Listeners can be attached.
 */
class DynamicParser : public StartOfFrameFoundListener, public SegmentFoundListener, public DataFoundListener
{
public:
    void update()
    {
    }
private:
     /**
      * Allows other classes to mesh in
      */
    class State
    {
    public:
        virtual void update() = 0;
    };

    class StartOfFrameSearch : public State, public DynamicSofSearch
    {
    public:
        using DynamicSofSearch::DynamicSofSearch;
        void update() override
        {
            DynamicSofSearch::update();
        }
    };

    class HeaderExtract: public State, public DynamicHeaderExtract 
    {
    public:
        using DynamicHeaderExtract::DynamicHeaderExtract;
        void update() override
        {
            DynamicHeaderExtract::update();
        }
    };

    class DataExtract: public State, public DynamicDataExtract
    {
    public:
        using DynamicDataExtract::DynamicDataExtract;
        void update() override
        {
            DynamicDataExtract::update();
        }
    };

    ByteQueue& byteBuffer_;
    DataFoundListener& dataFounListener_;
};
}
#endif
