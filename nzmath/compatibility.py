"""
compatibility between Python version
"""

# __builtins__.set is in >=2.4, sets.Set is in >=2.3.
# Be careful that the compatibility is not perfect.
try:
    set, frozenset
except NameError:
    import sets
    __builtins__["set"] = sets.Set
    __builtins__["frozenset"] = sets.ImmutableSet
    del sets

# __builtins__.cmp is only in <3.0.
# Be careful that the compatibility is not perfect.
try:
    cmp
except NameError:
    cmp = lambda x,y: (x>y) - (x<y)
    __builtins__["cmp"] = cmp
    del cmp

# builtin len() raises OverflowError when the result > sys.maxint.
# I'm not sure this restriction will go away in the future.
# The following function card() ought to be used instead of len()
# for obtaining cardinality of sets or set-like objects.
def card(virtualset):
    """
    Return cardinality of the virtualset.
    """
    if hasattr(virtualset, "card"):
        return virtualset.card()
    return len(virtualset)

__builtins__["card"] = card
