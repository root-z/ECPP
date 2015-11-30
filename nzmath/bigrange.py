"""
bigrange

Generators for range like sequences.
"""


def count(n=0):
    """
    Count up infinitely from 'n' (default to 0),

    see itertools.count
    """
    while True:
        yield n
        n += 1


def range(start, stop=None, step=None):
    """
    Generate range like finite integer sequence but can generate more
    than sys.maxint elements.
    """
    if step is None:
        step = 1
    elif not isinstance(step, (int, long)):
        raise ValueError("non-integer step for range()")
    if not isinstance(start, (int, long)):
        raise ValueError("non-integer arg 1 for range()")
    if stop is None:
        start, stop = 0, start
    elif not isinstance(stop, (int, long)):
        raise ValueError("non-integer stop for range()")

    if step > 0:
        n = start
        while n < stop:
            yield n
            n += step
    elif step < 0:
        n = start
        while n > stop:
            yield n
            n += step
    else:
        raise ValueError("zero step for range()")


def arithmetic_progression(init, difference):
    """
    Generate an arithmetic progression start form 'init' and
    'difference' step.
    """
    return _iterate(lambda x: x + difference, init)


def geometric_progression(init, ratio):
    """
    Generate a geometric progression start form 'init' and multiplying
    'ratio'.
    """
    return _iterate(lambda x: x * ratio, init)


def _iterate(func, init):
    """
    Generate (infinitely) a sequence (init, func(init), func(func(init)), ...)
    """
    val = init
    while True:
        yield val
        val = func(val)


def multirange(triples):
    """
    multirange(list of range triples) is an iterator over cartesian
    product of elements of ranges.

    For example, multirange([(1, 10, 3), (1, 10, 4)]) yields
    (1, 1), (1, 5), (1, 9), (4, 1), (4, 5), (4, 9), (7, 1), (7, 5),
    and (7, 9).

    The range triples may be doubles (start, stop) or single (stop,),
    but they have to be always tuples.

    Be cautious that using this module usually means you are trying to
    do brute force looping.
    """
    starts, stops, steps = _split_triples(triples)
    if _empty(starts, stops, steps):
        raise StopIteration("empty range")
    i, values = len(triples) - 1, list(starts)
    yield tuple(values)
    while i >= 0:
        values[i] += steps[i]
        if _range_exhausted(i, values, steps, stops):
            values[i] = starts[i]
            i -= 1
            continue
        yield tuple(values)
        i = len(triples) - 1

def _split_triples(triples):
    """
    Convert a list of triples into three lists each consists of 1st,
    2nd and 3rd element of each triples.

    The range triples may be doubles (start, stop) or single (stop,),
    but they have to be always tuples.
    """
    starts, stops, steps = [], [], []
    for triple in triples:
        if len(triple) == 3:
            start, stop, step = triple
        elif len(triple) == 2:
            start, stop, step = triple[0], triple[1], 1
        elif len(triple) == 1:
            start, stop, step = 0, triple[0], 1
        else:
            raise TypeError("tuple is longer or shorter than expected 1 to 3")
        starts.append(start)
        stops.append(stop)
        steps.append(step)
    return starts, stops, steps

def _empty(starts, stops, steps):
    """
    Report whether there is an empty range in the triples or not.
    """
    for start, stop, step in zip(starts, stops, steps):
        if step * (stop - start) <= 0:
            return True
    return False

def _range_exhausted(i, values, steps, stops):
    """
    Check if a range at i is exhausted.
    """
    return steps[i] * (values[i] - stops[i]) >= 0


def multirange_restrictions(triples, ascending=(), descending=(), strictly_ascending=(), strictly_descending=()):
    """
    multirange_restrictions is an iterator similar to the multirange
    but putting restrictions on each ranges.  A restriction ascending
    is a sequence that specifies the indices where the number emitted
    by the range should be greater than or equal to the number at the
    previous index.  Other restrictions descending, strictly_ascending
    and strictly_descending are similar.

    For example, while multirange([(1, 10, 3), (1, 10, 4)]) yields
    (1, 1), (1, 5), (1, 9), (4, 1), (4, 5), (4, 9), (7, 1), (7, 5),
    and (7, 9),
    multirange_restrictions([(1, 10, 3), (1, 10, 4)], ascending=(1,))
    yields only (1, 1), (1, 5), (1, 9), (4, 5), (4, 9) and (7, 9).

    The range triples may be doubles (start, stop) or single (stop,),
    but they have to be always tuples.
    """
    starts, stops, steps = _split_triples(triples)
    if _empty(starts, stops, steps):
        raise StopIteration("empty range")
    i, values = len(triples) - 1, list(starts)
    while i >= 0:
        if _range_exhausted(i, values, steps, stops):
            values[i] = starts[i]
            i -= 1
        elif (strictly_ascending and _violate_strictly_ascending(values, strictly_ascending) or
              ascending and _violate_ascending(values, ascending) or
              strictly_descending and _violate_strictly_descending(values, strictly_descending) or
              descending and _violate_descending(values, descending)):
            i = min(i + 1, len(triples) - 1)
        else:
            yield tuple(values)
            i = len(triples) - 1
        values[i] += steps[i]

def _violate_ascending(values, ascending):
    """
    Check if values violate ascending restrictions.
    """
    for i in ascending:
        if values[i - 1] > values[i]:
            return True
    return False

def _violate_strictly_ascending(values, ascending):
    """
    Check if values violate strictly-ascending restrictions.
    """
    for i in ascending:
        if values[i - 1] >= values[i]:
            return True
    return False

def _violate_descending(values, descending):
    """
    Check if values violate descending restrictions.
    """
    for i in descending:
        if values[i - 1] < values[i]:
            return True
    return False

def _violate_strictly_descending(values, descending):
    """
    Check if values violate strictly-descending restrictions.
    """
    for i in descending:
        if values[i - 1] <= values[i]:
            return True
    return False
