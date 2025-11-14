#ifndef COOL_SERIAL_DATA_HANDLER_HANDLER_MAP_HPP
#define COOL_SERIAL_DATA_HANDLER_HANDLER_MAP_HPP
#include <unordered_map>

#include "cool_serial/data_handler/idata_handler.hpp"
namespace coolSerial
{
    /**
     * This is used to connect IDataHandlers to a byte protocol
     *
     * WARNING: This does NOT check for bounds. Max protocol number is
     * 256. (Ids range from 0 to 255)
     */
    using HandlerMap = std::unordered_map<coolSerial::Byte, std::reference_wrapper<IDataHandler>>;
}
#endif
