import random

def bernoulli_test(p):
    """
    Returns True with probability 0 <= p <= 1
    """
    if 0. <= p <= 1.:
        if random.random() < p:
            return True
        else:
            return False
    else:
        raise Exception("wrong probability (%s) was passed" % p)
