#include "cool_serial/data_handler/handler_map.hpp"

namespace coolSerial
{
    /**
     * This listens for data which has been found (usually called by a parser).
     *
     * It will refer the data to the proper type depending on the data id.
     *
     * Make sure that all protocols that are being recieved have a handler.
     */
    class DataRouter
    {
    public:
        DataRouter(const HandlerMap& handlerMap):
            handlerMap_{handlerMap}
        {}

        DataRouter() : handlerMap_{}
        {}
    private:
        HandlerMap handlerMap_;
    };
}
