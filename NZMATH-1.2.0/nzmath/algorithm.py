"""
Abstract (higher order) functions
"""

import functools

def digital_method(coefficients, val, add, mul, act, power, zero, one):
    """
    Evaluate a univariate polynomial at val.

    The polynomial is given as an iterator 'coefficients' which yields
    (degree, coefficient) pair in descending order of degree.

    If the polynomial is of R-coefficients, then val should be in an
    R-algebra D.

    All operations 'add', 'mul', 'act', 'power', 'zero', 'one' should
    be explicitly given, where: 
      'add' means addition (D x D -> D),
      'mul' multiplication (D x D -> D),
      'act' action of R (R x D -> D),
      'power' powering (D x Z -> D) and
      'zero', 'one' constants in D.
    """
    result = zero
    val_pow = {0:one, 1:val}
    d = None
    for d_next, c_next in coefficients:
        if d:
            diff = d - d_next
            if diff not in val_pow:
                val_pow[diff] = power(val, diff)
            result = add(mul(result, val_pow[diff]), act(c_next, one))
        else:
            result = act(c_next, one)
        d = d_next
    if d:
        return mul(result, val_pow.get(d, power(val, d)))
    else:
        return result

def digital_method_func(add, mul, act, power, zero, one):
    """
    Return a function which evaluates polynomial (binary function)
    constructed from given functions.
    """
    result = functools.partial(digital_method, add=add, mul=mul,
                               act=act, power=power, zero=zero, one=one)
    result.___doc__ = """
    Evaluate a univariate polynomial at val.

    The polynomial is given as an iterator 'coefficients' which yields
    (degree, coefficient) pair in descending order of degree.

    If the polynomial is of R-coefficients, then val should be in an
    R-algebra D.
    """
    return result

def rl_binary_powering(element, index, mul, square=None, one=None):
    """
     powering by using right-left binary method
     (assume that index >= 0)
    """
    square = _set_square(mul, square)
    sol = one
    mul_part = element
    while True:
        if index & 1:
            try:
                sol = mul(sol, mul_part)
            except: # one is None!
                sol = mul_part
        index >>= 1
        if not index:
            return sol
        else:
            mul_part = square(mul_part)

def lr_binary_powering(element, index, mul, square=None, one=None):
    """
    powering by using left-right binary method.
    (assume that index >= 0)
    """
    import math
    square = _set_square(mul, square)
    if not index:
        return one
    spot = 1 << long(math.log(index, 2))
    sol = one
    while True:
        if spot & index:
            try:
                sol = mul(sol, element)
            except: # one is None!
                sol = element
        spot >>= 1
        if not spot:
            return sol
        sol = square(sol)

def window_powering(element, index, mul, square=None, one=None):
    """
    powering by using small-window method
    window size is selected by average analystic optimization
    (assume that index >= 0)
    """
    import math
    square = _set_square(mul, square)
    if not index:
        return one
    log_n = long(math.log(index, 2))
    # Find the proper window size
    size = 2
    pow_size = 2
    while log_n > (size + 1) * (size + 2) * (pow_size - 1):
        pow_size <<= 1
        size += 1
    # Precomputation
    sqr = square(element)
    pre_table = [element]
    pow_size -= 1
    for i in range(pow_size):
        pre_table.append(mul(pre_table[-1], sqr))

    near_n = 1 << log_n
    spot = near_n
    while spot:
        if not(spot & index):
            sol = square(sol)
            spot >>= 1
        else:
            # Find the window
            f_spot, e_spot = spot, spot >> (size - 1)
            t_size = size
            if not(e_spot):
                e_spot = 1
                t_size = int(math.log(f_spot, 2)) + 1
            while True:
                if e_spot & index:
                    spot = (e_spot >> 1)
                    window = 0
                    sqr = 1
                    while f_spot != e_spot: # Compute value of window
                        if index & e_spot:
                            window += sqr
                        sqr <<= 1
                        e_spot <<= 1
                    window += sqr
                    break
                t_size -= 1
                e_spot <<= 1
            # Compute a part of powering
            if f_spot == near_n:
                sol = pre_table[(window - 1) >> 1]
            else:
                for i in range(t_size):
                    sol = square(sol)
                sol = mul(sol, pre_table[(window - 1) >> 1])
    return sol

def _set_square(mul, square):
    """
    Set square if square is None
    (for the initialization of methods for powering)
    """
    if square == None:
        return lambda x : mul(x, x)
    else:
        return square

def powering_func(mul, square=None, one=None, type=0):
    """
    Return a function which computes powering (binary function)
    constructed from given functions.
    """
    if type == 0:
        result = functools.partial(
            rl_binary_powering, mul=mul, one=one, square=square)
        result.___doc__ = """
        compute powering, that is, element ** index.

        This function uses right-left binary method.
        (assume that index >= 0)
        """
    elif type == 1:
        result = functools.partial(
            lr_binary_powering, mul=mul, one=one, square=square)
        result.___doc__ = """
        compute powering, that is, element ** index.

        This function uses left-right binary method.
        (assume that index >= 0)
        """
    elif type == 2:
        result = functools.partial(
            window_powering, mul=mul, one=one, square=square)
        result.___doc__ = """
        compute powering, that is, element ** index.

        This function uses window method.
        (assume that index >= 0)
        """
    return result
