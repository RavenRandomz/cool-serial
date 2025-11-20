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

### Header Section
0. Start of Frame: 132 - Z in ascii
1. Length uint\_16 most significant byte
2. Length uint\_16 byte least significant
3. Data type byte : types defined by user
4. CRC8 for bytes 1-3
CRC8/SMBUS, specifically:
check  - 0xF4
poly   - 0x07
init   - 0x00
refin  - false
refout - true
xorout - 0x00

Data Section
-----------
5. Data bytes 
