#ifndef COOL_SERIAL_CONTINUOUS_PARSE_HPP
#define COOL_SERIAL_CONTINUOUS_PARSE_HPP

#include "byte_queue.hpp"
#include "cool_message.hpp"
#include "header_section.hpp"

namespace coolSerial
{
    /**
     * This class connects to a queue of bytes.
     *
     *
     */
    class ContinuousParser
    {
    public:
        enum class StateType
        {
            startOfFrameSearch,
            headerParse,
            dataParse,
            updateMessage
        };

        class State 
        {
        public:
            virtual void update() = 0;
            virtual ~State() {};
        };

        class UpdateMessage : public State
        {
        public:
            UpdateMessage(ContinuousParser& parser) : parser_{parser}
            {
            }

            void update() override
            {

                parser_.setCurrentMessage
                    (
                        CoolMessageData
                        {
                            .dataType = dataType_,
                            .data = std::move(messageData_)
                        }
                    );
                messageData_ = Bytes{};
            }

            void setData(Bytes bytes)
            {
                messageData_ = std::move(bytes);
            }

            void setDataType(Byte dataType)
            {
                dataType_ = dataType;
            }
        private:
            ContinuousParser& parser_;
            Bytes messageData_;
            Byte dataType_;
        };

        class StartOfFrameSearch : public State
        {
        public:
            StartOfFrameSearch(ContinuousParser& parser) : parser_{parser}
            {}

            void update() override
            {
                //Consider having an upper limit
                while(parser_.byteAvailable())
                {
                    const Byte kNext{parser_.getNextPoppedByte()};
                    if(kNext == CoolMessage::kStartOfFrame)
                    {
                        // We can assume that the next 5 bytes are header bytes.
                        parser_.setState(StateType::headerParse);
                        parser_.update();
                        break;
                    }
                }
            }

        private:
            ContinuousParser& parser_; 
        };

        class DataParse : public State
        {
        public:
            DataParse(ContinuousParser& parser, UpdateMessage& updateMessage)
                :
                parser_{parser},
                updateMessage_{updateMessage}
            {
            }

            void update() override
            {
                while(parser_.byteAvailable())
                {
                    dataBytes_.push_back(parser_.getNextPoppedByte());
                    // Current byte is being proccessed so it doesn't count
                    --bytesRemaining_;

                    if(bytesRemaining_ == 0)
                    {
                        passToMessageUpdater();
                        break;
                    }

                }
            }

            void passToMessageUpdater()
            {
                parser_.setState(StateType::updateMessage);
                updateMessage_.setDataType(dataInfo_.dataType);
                updateMessage_.setData(std::move(dataBytes_));
                // Add new data to prevent crash (potential)
                dataBytes_ = Bytes{};
                parser_.update();
            }

            void setDataInfo(const DataInfo& dataInfo)
            {
                dataInfo_ = dataInfo;
                bytesRemaining_ = dataInfo.dataLength;
            }
        private:
            ContinuousParser& parser_;
            UpdateMessage& updateMessage_;
            DataInfo dataInfo_{};
            int bytesRemaining_{0};
            Bytes dataBytes_{};
        };

        class HeaderParse : public State
        {
        public:
            HeaderParse(ContinuousParser& parser, DataParse& dataParse)
                :
                parser_{parser},
                dataParse_{dataParse}
            {}

            void update() override
            {
                while(parser_.byteAvailable())
                {

                    headerBytes_[headerByteIndex_] = parser_.getNextPoppedByte();
                    std::cout << headerByteIndex_;

                    if(headerByteIndex_ == kHeaderByteCountMaxIndex)
                    {
                        proccessHeaderBytes();
                        break;
                    }
                    // Prepare next interation for next index
                    ++headerByteIndex_;
                }
            }

            void proccessHeaderBytes()
            {
                // Reset for next iteration
                headerByteIndex_ = 0;

                 const auto kHeaderData{HeaderSection::deserializeBytes(headerBytes_)};
                if (kHeaderData.isValid())
                {
                    parser_.setState(StateType::dataParse);
                    dataParse_.setDataInfo(kHeaderData.dataInfo);
                    parser_.update();
                }
                else 
                {
                    // Dump header and associated data until next start of frame
                    parser_.setState(StateType::startOfFrameSearch);
                    // Wait until next update
                }
            }


        private:
            static const int kHeaderByteCountMaxIndex{3};

            int headerByteIndex_{0};
            HeaderBytes headerBytes_{};
            ContinuousParser& parser_; 
            DataParse& dataParse_;
        };

        ContinuousParser(ByteQueue& byteQueue) : byteQueue_{byteQueue}
        {}

        void update()
        {
            state_->update();
        }

        void setState(StateType stateType)
        {
            switch(stateType)
            {
                case StateType::dataParse:
                    state_ = &dataParse_;
                    break;
                case StateType::headerParse:
                    state_ = &headerParse_;
                    break;
                case StateType::updateMessage:
                    state_ = &updateMessage_;
                    break;
                case StateType::startOfFrameSearch:
                    state_ = &startOfFrameSearch_;
                    break;
            }
        }


        void setCurrentMessage(const CoolMessageData& messageData)
        {
            currentMessage_ = messageData;
        }

        CoolMessageData getCurrentMessage() const 
        {
            return currentMessage_;
        }
        bool byteAvailable()
        {
            return !byteQueue_.empty();
        }

        Byte getNextPoppedByte()
        {
            const Byte kNext{byteQueue_.front()};
            byteQueue_.pop();
            return kNext;
        }
    private:

        // State handlers
        StartOfFrameSearch startOfFrameSearch_{*this};
        UpdateMessage updateMessage_{*this};
        DataParse dataParse_{*this, updateMessage_};
        HeaderParse headerParse_{*this, dataParse_};

        ByteQueue& byteQueue_;
        State* state_{&startOfFrameSearch_};
        CoolMessageData currentMessage_{};
    };
}
#endif
