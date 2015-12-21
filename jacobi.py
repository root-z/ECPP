def jacobi(a, n):
    """

    Args:
        (a/n)
    Returns:
        Jacobi Symbol
    """
    symbol = 1
    while True:
        if n < 1 or n % 2 == 0:
            raise ValueError("Illegal n")
        if n == 1:
            return symbol
        a = a % n
        if a == 0:
            return 0
        if a == 1:
            return symbol

        if a & 1 == 0:
            if n % 8 in (3, 5):
                symbol = -symbol
            a >>= 1
            continue

        if a % 4 == 3 and n % 4 == 3:
            symbol = -symbol

        a, n = n, a
    return
