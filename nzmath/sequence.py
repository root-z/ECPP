def generator_fibonacci(n=None):
    """
    Generate Fibonacci number up to n-th term if n is assigned
    else infinity.
    """
    a = 0
    b = 1

    if None == n:
        while True:
            yield b
            a += b
            yield a
            b += a
    else:
        count = 0
        while True:
            yield b
            count += 1
            if n <= count:
                break
            a += b

            yield a
            count += 1
            if n <= count:
                break
            b += a


FIBONACCI = {0:0, 1:1}
def fibonacci(n):
    """
    param non-negative integer n
    return the n-th term of the Fibonacci
    effect FIBONACCI[n] = fibonacci(n)
    """
    if n < 0:
        raise ValueError("fibonacci(n)  0 <= n  ?")

    if n in FIBONACCI:
        return FIBONACCI[n]

    m = n >> 1
    if n & 1 == 0:
        f1 = fibonacci(m - 1)
        f2 = fibonacci(m)
        FIBONACCI[n] = (f1 + f1 + f2) * f2
    else:
        f1 = fibonacci(m)
        f2 = fibonacci(m + 1)
        FIBONACCI[n] = f1 ** 2 + f2 ** 2

    return FIBONACCI[n]
