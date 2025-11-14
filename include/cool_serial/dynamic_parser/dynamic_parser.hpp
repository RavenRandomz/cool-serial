#ifndef COOL_SERIAL_DYNAMIC_PARSER_DYNAMIC_PARSER_HPP
#define COOL_SERIAL_DYNAMIC_PARSER_DYNAMIC_PARSER_HPP
#include "cool_serial/byte_queue.hpp"

#include "cool_serial/dynamic_parser/dynamic_sof_search.hpp"
#include "cool_serial/dynamic_parser/dynamic_data_extract.hpp"
#include "cool_serial/dynamic_parser/dynamic_header_extract.hpp"
#include "cool_serial/dynamic_parser/data_found_listener.hpp"
#include <cool_serial/dynamic_parser/header_found_listener.hpp>
#include <cool_serial/dynamic_parser/segment_found_listener.hpp>
#include <cool_serial/dynamic_parser/start_of_frame_found_listener.hpp>
#include <functional>

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
class DynamicParser : public StartOfFrameFoundListener, public HeaderFoundListener, public DataFoundListener
{
public:
    DynamicParser(ByteQueue& byteBuffer, DataFoundListener& listener):
        dataFoundListener_{listener},
        startOfFrameSearch_{byteBuffer, *this},
        headerExtract_{byteBuffer, *this},
        dataExtract_ {byteBuffer, *this}
    {}
    void update()
    {
        state_.get().update();
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

    class HeaderExtract : public State, public DynamicHeaderExtract 
    {
    public:
        using DynamicHeaderExtract::DynamicHeaderExtract;
        void update() override
        {
            DynamicHeaderExtract::update();
        }
    };

    class DataExtract : public State, public DynamicDataExtract
    {
    public:
        using DynamicDataExtract::DynamicDataExtract;
        void update() override
        {
            DynamicDataExtract::update();
        }
    };

    DataFoundListener& dataFoundListener_;

    StartOfFrameSearch startOfFrameSearch_;
    HeaderExtract headerExtract_;
    DataExtract dataExtract_;

    std::reference_wrapper<State> state_{startOfFrameSearch_};

    // State change logic
    void startOfFrameFound() override
    {
        state_ = headerExtract_;
        // In case the buffer has more of the message
        update();
    }

    void headerFound(const HeaderData& headerData) override
    {
        if (headerData.isValid())
        {
            dataExtract_.setDataExtractionInfo(headerData.dataInfo);
            state_ = dataExtract_;
            // In case the buffer has more of the message
            update();
        }
        else
        {
            // Do not call update, wait until next cycle
            state_ = startOfFrameSearch_;
        }
    }

    void dataFound(const CoolMessageData& data) override
    {
        // Do not call update, wait until next cycle
        // at most one message is extracted every update()
        dataFoundListener_.dataFound(data);
        state_ = startOfFrameSearch_;
    }
};
}
#endif
