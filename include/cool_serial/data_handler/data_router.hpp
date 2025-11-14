#ifndef COOL_SERIAL_DATA_HANDLER_DATA_ROUTER_HPP
#define COOL_SERIAL_DATA_HANDLER_DATA_ROUTER_HPP

#include "cool_serial/data_handler/handler_map.hpp"
#include "cool_serial/data_handler/idata_handler.hpp"
#include "cool_serial/dynamic_parser/data_found_listener.hpp"

namespace coolSerial
{
    /**
     * This listens for data which has been found (usually called by a parser).
     *
     * It will refer the data to the proper type depending on the data id.
     *
     * Make sure that all protocols that are being recieved have a handler.
     */
    class DataRouter : public DataFoundListener
    {
    public:
        DataRouter(const HandlerMap& handlerMap):
            handlerMap_{handlerMap}
        {}

        DataRouter() : handlerMap_{}
        {}

        /**
         * Change a handler.
         *
         * If the handler has not been set, it will add a new map
         */
        void setMapping(Byte dataType, IDataHandler& handler)
        {
            handlerMap_.at(dataType) = handler;
        }

        void dataFound(const CoolMessageData& data)
        {
            IDataHandler& handler{handlerMap_.at(data.type).get()};
            handler.handleData(data.data);
        }
    private:
        HandlerMap handlerMap_;
    };
}
#endif
