from crc import Calculator, Configuration

config = Configuration(
    width = 8,
    check = 0xf4,
    polynomial = 0x07,
    init_value = 0x00,
    reverse_input = False,
    reverse_output = True,
    xor_output = 0x00
)

calculator = Calculator(config)