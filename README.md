Cool Serial
===========
Error correction is defined by each message type. The message is usually a 
serialized object which is simply deserialized. UART's scope is to merely 
transmit these serialized objects in a minimal yet error-resistant manner.

Cool Serial is a transit protocol. It's main purpose is conveying an arbitrary 
sequence of bytes which are then sent to whatever handles them.

Normally, when using CoolSerial on a practical level, you will need two 
protocls: Cool Serial and another which defines how that data section works.

## The Protocol

### SOF Section
0. Start of Frame: 132 - Z in ascii
### Header Section
1. Length uint\_16 most significant byte
2. Length uint\_16 byte least significant
3. Data type byte : types defined by user
4. CRC8 for bytes 1-3
    * CRC8/SMBUS, specifically:
    * check  - 0xF4
    * poly   - 0x07
    * init   - 0x00
    * refin  - false
    * refout - true
    * xorout - 0x00

### Data Section
5. Data bytes

# cool-serial basic API

This section goes over the important classes and methods used in cool-serial. 

cool-serial provides a C++ implementation of the CoolSerial protocol. It does 
this by providing: the ability to generate a series of bytes which adhere to the 
protocol, a parser which can take an incomplete series of bytes until it is 
used, and a router which will forward intercepted data type.

cool-serial is implemented with unit tests in mind. The mock subdirectory 
provides mock classes for all of the interfaces.

**Note that these examples assume that the code takes place within the 
coolSerial namespace. In your codebase, you must do coolSerial::Bytes or 
coolSerial::CoolMessage**

## Message Classes
These implement data sections

### CoolMessage

```cpp
Bytes bytes {1, 2, 3, 4, 5};
Byte dataType{20};

CoolMessage message{data, dataType};
```
The data which contains the serialized bytes of the protocol which corresponds 
to the dataType must be passed into the constructor

```cpp
byteQueue.addBytes(message.getFrame());
```

In order to get the serialized message bytes, use getFrame(). This is a raw 
coolSerial message. This can be added to a queue.

### CoolMessageData

This simply wraps the byets and a data type together. This is the format that is 
recieved by those who implement IDataFoundListener

```cpp
Bytes dataBytes{1, 2, 3, 4, 5};
Byte dataType{20};

CoolMessageData coolData{.dataType = dataType, .data = dataBytes};
```

### DataInfo
This class stores the data length and data type. It is capable of serializing 
the 1st 3 bytes in a header section after the SOF.

```cpp
DataInfo dataInfo
{
    .dataLength = 10
    .dataType = 1
};
```
There will be rename so that it's dataInfo.length isntead of dataInfo.dataLength 
:)

```cpp
DataInfoBytes bytes {dataInfo.serialize()};
```

Note that DataInfoBytes is an std::array, not a vector. This is for performance 
reasons. Remember that these classes are mainly for implementation of the 
protocol and parsing rather than regular use.

### HeaderSection
The header section class is seldom worked with individually. It mainly involves 
the serialization, deserialization, and verification of headers within the 
implementation.

There are common std::array instances called HeaderBytes and DataInfoBytes. 
These are for performance reasons.

The DataInfoBytes are for specifically the 1st three bytes of the data section.
HeaderBytes contains the entire header (minus the SOF byte)

```cpp
HeaderSection headerSection{dataInfo};
HeaderBytes bytes{headerSection.getSerialized()}; // obtain serialized bytes
HeaderData dataInfo{headerSection.getHeaderData}
```
HeaderSection stores its data as the actual serialized bytes (it's meant to 
contain the header section of a coolSerial message).


### Bytes

This is an alias for std::vector<Byte>;

### Byte

Byte is an alias for unsigned char

## Transmission Classes

These classes are meant to relay data to and from cool-seril

### ByteQueue
Byte Queue adds syntactic sugar to a standard Queue

```cpp
ByteQueue buffer{};

buffer.push(getNextByte());
```
You'll have to interface with a uart library to update it reguarly.


```cpp
const Byte kByte {buffer.getNextPoppedByte()};
```

```
buffer.hasNextByte();
```
Check if it is empty or not.



## Reading Classes

These classes are involved in taking a stream of bytes which may or may not be a 
complete cool serial frame.

### DynamicParser
It uses a state machine design on whether or not to search for SOF, attempt to 
extract a header, or complete data. Update is reguarly called. 

Whenever it extracts a cool serial message, it will forward it to that object to 
handle.

```cpp
ByteQueue buffer{};
DataRouter router{};

// pass by reference
DynamicParser parser{buffer, dataRouter};
```
Dynamic parser requires a ByteQueue in its constructor. In addition, it requires 
a DataFoundListener which it will call void dataFound(CoolMessageData).

```cpp
// Assume the buffer has been filled with data but not enough
parser.update()

// the router will not have its dataFound() called.

// Assume at least the rest of the bytes has been added to the buffer

parser.update()
```

So every "program cyle" the parser will have update() called to process whether 
is in the buffer. Note that this will cause the buffer to empty bytes.
### ContinuousParser

This is a deprecated implementation. It requires the user to reguarly check for 
if the message has been processed or not every update(). This is very wasteful 
as it is unlikely for all of the frame's bytes to be within the ByteBuffer.

It is kept for backwards compatibility reasons. It will be removed in v2.0.0

## Handling Classes

These classes take data which has been found from a parser and send them to 
whatever class is meant to handle them.

### IDataHandler

Any class that implements handleData(const Bytes& data);

It can update a currentClass (currentFloat for example) or pass the deserialized 
data (if it deserialized the data) onto some listeners.

### DataMap

```cpp
testing::StrictMock<IDataHandlerMock> handler0{};
testing::StrictMock<IDataHandlerMock> handler5{};

HandlerMap testMap
{
    {0, handler0},
    {5, handler5}
};
```
Note that in its constructor, it maps a uint\_8 to a 
std::reference\_wrapper<IDataHandler>. So you don't have to pass by address.

### DataRouter

Data Router takes a handler map. It also is a DataFoundListenerWhenever it has 
```
DataRouter router{testMap};
```
DataRouter requires a handler map which tells it which handlers to call 
handleData(const Bytes& data). This is triggered when it's dataFound() function 
has been called. This can be manually called, but 

```
handler.dataFound(data);
```
Whatever handler is associated to the data's type in the test map will be 
called. The data's data (the actual bytes) will be 

# Basic Operations

This provides an overview on how to use the CoolSerial classes


## Writing a message


```cpp
struct TestStruct
{
    float x;
    float y;
    double z;
};

TestStruct data
{
    .x = 2.234,
    .y = -345345.234234,
    .z = 34535.345345
};

Bytes kDataBytes{cista::serialize(data)};

CoolMessage message{std::move(kDataBytes), Byte{0}};

const Bytes kMessageFrame{message.getFrame()};
```
Note that cista is a serialization/deserialization library.

CoolMessage requires a data type (uint8) and the data itself (bytes).
## Sending a message
You will need a uart library to push bytes from a message frame:

for example, assuming the function is uart::tx::sendbyte(Byte); (this will 
depend on which serial communication library you use).
```
for (const kByte& byte : kMessageFrame)
{
    uart::tx::sendByte(byte);
}
```

## Receiving a message
Assuming the function is uart::rx::popByte() and uart::rx::byteAvailable();
```
while (uart::rx::byterAvailable())
{
    // a coolSerial buffer
    buffer.push(uart::rx::popByte());
}
```

From there, if the buffer is passed via reference to a DynamicParser, the 
dynamic parser can send it to a DataFoundListener. For your convenience, 
DataRouter can take maps

## Handling Incoming Messages

As long as a class implements DataFoundListener (in dynamic\_parser 
subdirectory), it can be passed by reference to a DynamicParser's constructor. 
Usually, you would want to use the premade well-tested DataRouter. But if you 
want to get adventurous, you can implement your own classes, too.
