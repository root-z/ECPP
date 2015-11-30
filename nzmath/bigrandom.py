#bigrandom.py

import sys
import random as _random
import nzmath.arith1 as arith1


def randrange(start, stop=None, step=1):
    """
    Choose a random item from range([start,] stop[, step]).
    (Return long integer.)

    see random.randrange
    """
    if stop is None:
        if abs(start) < sys.maxint:
            return _random.randrange(start)
    elif abs(stop - start) < sys.maxint:
        return _random.randrange(start, stop, step)

    negative_step = False
    if stop is None:
        start, stop = 0, start
    if step < 0:
        start, stop, step = -start, -stop, -step
        negative_step = True

    _validate_for_randrange(start, stop, step)

    if (stop - start) % step != 0:
        v = (stop - start)//step + 1
    else:
        v = (stop - start)//step

    bitlength = arith1.log(v, base=2)
    if v != 2**bitlength:
        bitlength += 1
    randomized = v + 1 # to pass the first test
    while randomized >= v:
        randomized = 0
        for i in range(bitlength):
            if random() < 0.5:
                randomized += 1 << i
    result = randomized * step + start
    if negative_step:
        return -result
    return result

def _validate_for_randrange(start, stop, step):
    """
    Check validity of arguments for randrange.
    """
    if step == 0:
        raise ValueError("zero step for randrange()")
    elif start != long(start):
        raise ValueError("non-integer arg 1 for randrange()")
    elif stop != long(stop):
        raise ValueError("non-integer stop for randrange()")
    elif step != long(step):
        raise ValueError("non-integer step for randrange()")
    if start >= stop:
        raise ValueError("empty range for randrange()")

random = _random.random
seed = _random.seed


def map_choice(mapping, upperbound):
    """
    Return a choice from a set given as the image of the mapping from
    natural numbers (more precisely range(upperbound)).  In other
    words, it is equivalent to
    random.choice([mapping(i) for i in range(upperboud) if mapping(i) != None])
    if upperbound is small enough for the list size limit.

    The mapping can be a partial function, i.e. it may return None for
    some input.  However, if the resulting set is empty, it will end
    up with an infinite loop.
    """
    result = None
    while result is None:
        result = mapping(randrange(upperbound))
    return result


__all__ = ['random', 'randrange', 'seed', 'map_choice']
