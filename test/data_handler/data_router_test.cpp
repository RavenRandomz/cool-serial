#include "cool_serial/data_handler/data_router.hpp"
#include "cool_serial/mock/data_handler/idata_handler_mock.hpp"

#include <cool_serial/cool_message.hpp>
#include <cool_serial/data_handler/handler_map.hpp>
#include <gtest/gtest.h>

namespace coolSerial::dataRouterTest
{
/**
 * This checks for routing of different data wit different
 * types to different handlers.
 *
 * The appropriate handler should be called
 */
TEST(DataRouter, correctRouting)
{
    testing::StrictMock<IDataHandlerMock> handler0{};
    testing::StrictMock<IDataHandlerMock> handler5{};
    HandlerMap testMap
    {
    {0, handler0},
    {5, handler5}
    };
    DataRouter router{testMap};

    Bytes bytes0{1,2,3,4,5};

    CoolMessageData data0{0, bytes0};

    EXPECT_CALL(handler0, handleData(bytes0));
    router.dataFound(data0);
}
}
