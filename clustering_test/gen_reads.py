#!/usr/bin/python

import argparse
import random
import utils

DEFAULT_READ_LEN = 10
DEFAULT_COVERAGE_MEAN = 5.
DEFAULT_COVERAGE_STDEV = 1.
MISREAD_IDXS = dict(enumerate("NATGC"))
MISREAD_PROBS = [
                 [1. , 0. , 0. , 0. , 0. ], # from N
                 [.02, .9 , .04, .02, .02], # from A
                 [.02, .04, .9 , .02, .02], # from T
                 [.01, .03, .02, .9 , .04], # from G
                 [.01, .03, .02, .04, .9 ]  # from C
                ]

def check_misread_probs(probs):
    if not all(sum(prob) == 1.0 for prob in probs):
        raise Exception("misread probabilities for each base " \
                        "should add up to 1.0!")

check_misread_probs(MISREAD_PROBS)


descr_str = """Generates sample reads from provided DNA strands (of equal length).
               Warning: misread probabilites are hardcoded, please edit this file
               to change them"""
parser = argparse.ArgumentParser(description=descr_str)

parser.add_argument("-l", "--length", type=int,
                    required=False, default=DEFAULT_READ_LEN,
                    dest="read_len",
                    help="length of read, default: %s" % DEFAULT_READ_LEN)
parser.add_argument("-m", "--coverange_mean", type=float,
                    required=False, default=DEFAULT_COVERAGE_MEAN,
                    dest="cov_mean",
                    help="coverage mean, default: %s" % DEFAULT_COVERAGE_MEAN)
parser.add_argument("-d", "--coverage_stdev", type=float,
                    required=False, default=DEFAULT_COVERAGE_STDEV,
                    dest="cov_stdev",
                    help="coverage stdev, default: %s" % DEFAULT_COVERAGE_STDEV)

args = parser.parse_args()
