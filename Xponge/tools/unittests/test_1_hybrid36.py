"""
    Hybrid-36 reference implementation and tests
"""
digits_upper = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
digits_lower = digits_upper.lower()
digits_upper_values = {digit: value for value, digit in enumerate(digits_upper)}
digits_lower_values = {digit: value for value, digit in enumerate(digits_lower)}

__all__ = ["test_hybrid36_reference_decode",
           "test_hybrid36_reference_cases",
           "test_hybrid36_invalid_literals"]


def encode_pure(digits, value):
    """
        encode value using the given digits
    """
    assert value >= 0
    if value == 0:
        return digits[0]
    base = len(digits)
    result = []
    while value != 0:
        rest = value // base
        result.append(digits[value - rest * base])
        value = rest
    result.reverse()
    return "".join(result)


def decode_pure(digits_values, text):
    """
        decode string using the digit/value mapping
    """
    result = 0
    base = len(digits_values)
    for char in text:
        result = result * base + digits_values[char]
    return result


def hy36encode(width, value):
    """
        encode value as base-10/upper-case base-36/lower-case base-36 hybrid
    """
    i = value
    if i >= 1 - 10 ** (width - 1):
        if i < 10 ** width:
            return ("%%%dd" % width) % i
        i -= 10 ** width
        if i < 26 * 36 ** (width - 1):
            i += 10 * 36 ** (width - 1)
            return encode_pure(digits_upper, i)
        i -= 26 * 36 ** (width - 1)
        if i < 26 * 36 ** (width - 1):
            i += 10 * 36 ** (width - 1)
            return encode_pure(digits_lower, i)
    raise ValueError("value out of range.")


def hy36decode(width, text):
    """
        decode base-10/upper-case base-36/lower-case base-36 hybrid
    """
    if len(text) == width:
        first = text[0]
        if first == "-" or first == " " or first.isdigit():
            try:
                return int(text)
            except ValueError:
                if text == " " * width:
                    return 0
        elif first in digits_upper_values:
            try:
                return decode_pure(digits_upper_values, text) - 10 * 36 ** (width - 1) + 10 ** width
            except KeyError:
                pass
        elif first in digits_lower_values:
            try:
                return decode_pure(digits_lower_values, text) + 16 * 36 ** (width - 1) + 10 ** width
            except KeyError:
                pass
    raise ValueError("invalid number literal.")


def test_hybrid36_reference_decode():
    """
        compare Xponge decode with reference implementation
    """
    import Xponge.load as load
    for width in (4, 5):
        for value in (-999, -78, -6, 0, 12, 345, 6789):
            text = hy36encode(width, value)
            assert load._pdb_hybrid36_decode(width, text) == hy36decode(width, text)


def test_hybrid36_reference_cases():
    """
        check representative boundary cases from the reference exercise
    """
    import Xponge.load as load
    cases = [
        (4, "    ", 0),
        (4, "  -0", 0),
        (4, "-999", -999),
        (4, "  -6", -6),
        (4, "   0", 0),
        (4, "9999", 9999),
        (4, "A000", 10000),
        (4, "A00Z", 10035),
        (4, "A010", 10036),
        (4, "AZZZ", 10000 + 36 ** 3 - 1),
        (4, "ZZZZ", 10000 + 26 * 36 ** 3 - 1),
        (4, "a000", 10000 + 26 * 36 ** 3),
        (4, "azzz", 10000 + 26 * 36 ** 3 + 36 ** 3 - 1),
        (4, "zzzz", 10000 + 2 * 26 * 36 ** 3 - 1),
        (5, "     ", 0),
        (5, "   -0", 0),
        (5, "-9999", -9999),
        (5, "    0", 0),
        (5, "99999", 99999),
        (5, "A0000", 100000),
        (5, "A000Z", 100035),
        (5, "A0010", 100036),
        (5, "AZZZZ", 100000 + 36 ** 4 - 1),
        (5, "ZZZZZ", 100000 + 26 * 36 ** 4 - 1),
        (5, "a0000", 100000 + 26 * 36 ** 4),
        (5, "azzzz", 100000 + 26 * 36 ** 4 + 36 ** 4 - 1),
        (5, "zzzzz", 100000 + 2 * 26 * 36 ** 4 - 1),
    ]
    for width, text, value in cases:
        assert hy36decode(width, text) == value
        assert load._pdb_hybrid36_decode(width, text) == value


def test_hybrid36_invalid_literals():
    """
        invalid literal handling
    """
    import Xponge.load as load
    invalid = [
        (4, ""),
        (4, "    0"),
        (4, " abc"),
        (4, "abc-"),
        (4, "A=BC"),
        (4, "40a0"),
        (4, "40A0"),
        (5, ""),
        (5, "     0"),
        (5, " abcd"),
        (5, "ABCD-"),
        (5, "a=bcd"),
        (5, "410b0"),
        (5, "410B0"),
    ]
    for width, text in invalid:
        try:
            hy36decode(width, text)
        except ValueError as exc:
            assert str(exc) == "invalid number literal."
        else:
            raise RuntimeError("Exception expected.")
        try:
            load._pdb_hybrid36_decode(width, text)
        except ValueError as exc:
            assert str(exc) == "invalid number literal."
        else:
            raise RuntimeError("Exception expected.")
