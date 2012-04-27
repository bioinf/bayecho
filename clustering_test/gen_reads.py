#!/usr/bin/python

import argparse
import random
import utils
import sys

STRAND_CHARACTERS = "TAGC"
DEFAULT_READ_LEN = 10
DEFAULT_COVERAGE_MEAN = 5.
DEFAULT_COVERAGE_STDEV = 1.
MISREAD_IDX_TO_BASE = dict(enumerate("NATGC"))
MISREAD_BASE_TO_IDX = dict((b, a) for a, b in enumerate("NATGC"))
MISREAD_PROBS = [
                 [1. , 0. , 0. , 0. , 0. ], # from N
                 [.02, .9 , .04, .02, .02], # from A
                 [.02, .04, .9 , .02, .02], # from T
                 [.01, .03, .02, .9 , .04], # from G
                 [.01, .03, .02, .04, .9 ]  # from C
                ]

def check_misread_probs(probs):
    if not all(sum(prob) == 1.0 for prob in probs):
        raise Exception("ERROR: misread probabilities for each base " \
                        "should add up to 1.0")

def gen_read(strand, read_len):
    startpos = random.randint(0, len(strand) - read_len - 1)
    seq = list(strand[startpos:startpos+read_len])
    for i in xrange(len(seq)):
        probs = MISREAD_PROBS[MISREAD_BASE_TO_IDX[seq[i]]]
        seq[i] = MISREAD_IDX_TO_BASE[utils.loaded_dice(probs)]
    return "".join(seq)


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
check_misread_probs(MISREAD_PROBS)

strands = []
s = ""
while True:
    s = sys.stdin.readline().strip()
    if len(s) > 0:
        strands.append(s)
    else:
        break


strands_n = len(strands)
if strands_n == 0:
    exit("ERROR: no strands provided")

strands_len = len(strands[0])
if not all(len(s) == strands_len for s in strands):
    exit("ERROR: strands should have equal lengthes")
if not all(all(c in set(STRAND_CHARACTERS) for c in s) for s in strands):
    exit("ERROR: strands can contain only %s" % STRAND_CHARACTERS)

reads = []
for strand in strands:
    for i in xrange(int(round(random.gauss(args.cov_mean, args.cov_stdev)))):
        reads.append(gen_read(strand, args.read_len))

random.shuffle(reads)

for r in reads:
    print r
