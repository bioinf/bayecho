import random

def bernoulli_test(p):
    """
    Returns True with probability 0 <= p <= 1
    """
    if not 0. <= p <= 1.:
        raise Exception("wrong probability (%s) was passed" % p)
    if random.random() < p:
        return True
    else:
        return False

def loaded_dice(ps):
    cum_ps = [sum(ps[:i+1]) for i in xrange(len(ps))]
    if cum_ps[-1] != 1.0:
        raise Exception("wrong probability sum")
    t = random.random()
    for i in xrange(len(ps)):
        if t < cum_ps[i]:
            return i
