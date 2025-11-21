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

## IDataHandler
## DataRouter

# Basic Operations

This provides an overview on how to use the CoolSerial classes


## Writing a message

## Sending a message
## Receiving a message
## Handling a data

## Handling Incoming Messages



