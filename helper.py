class BytesIntEncoder:
    @staticmethod
    def encode(s: str) -> int:
        return int.from_bytes(s.encode(), byteorder='big')

    @staticmethod
    def decode(i: int) -> str:
        return i.to_bytes(((i.bit_length() + 7) // 8), byteorder='big').decode()
