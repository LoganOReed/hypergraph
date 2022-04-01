def square(num):
    """Squares a given number.

    Takes an integer `num` and returns the square of `num`.

    Parameters
    ----------
    num : int
        The integer that will be squared

    Returns
    -------
    out : int
        The integer that is the square of `num`
    """
    return num * num


def test_square():
    assert square(3) == 9
