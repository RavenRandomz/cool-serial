Cool Serial
===========
Error correction is defined by each message type. The message is usually a serialized object which is simply deserialized.
UART's scope is to merely transmit these serialized objects in a minimal yet error-resistant manner.

Header Section
---------------
0. Start of Frame: 0xAE
1. Length uint_16 byte no. 1
2. Length uint_16 byte no. 2
3. Message type byte : types defined by user
4. CRC8 for bytes 2-5

Data Section
-----------
6. Data bytes 
