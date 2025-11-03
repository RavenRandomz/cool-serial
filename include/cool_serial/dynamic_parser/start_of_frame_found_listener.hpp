#ifndef COOL_SERIAL_DYNAMIC_PARSER_START_OF_FRAME_FOUND_LISTENER_HPP
#define COOL_SERIAL_DYNAMIC_PARSER_START_OF_FRAME_FOUND_LISTENER_HPP
namespace coolSerial
{
    class StartOfFrameFoundListener
    {
    public:
        virtual void startOfFrameFound() = 0;
        virtual ~StartOfFrameFoundListener() {}
    };
}
#endif
