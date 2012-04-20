#!/usr/bin/python

import argparse
import random
import utils

CHARACTERS = "TAGC"
DEFAULT_STRAND_LEN = 50
DEFAULT_STRAND_NUM = 3
DEFAULT_SNP_PROB = 0.03


descr_str = """Generates sample DNA strands with SNPs, ensures that
               strands are different"""
parser = argparse.ArgumentParser(description=descr_str)

parser.add_argument("-n", "--number", type=int,
                    required=False, default=DEFAULT_STRAND_NUM,
                    dest="num",
                    help="number of strands, default: %s" % DEFAULT_STRAND_NUM)
parser.add_argument("-l", "--length", type=int,
                    required=False, default=DEFAULT_STRAND_LEN,
                    dest="len",
                    help="length of strand, default: %s" % DEFAULT_STRAND_LEN)
parser.add_argument("-p", "--prob", type=int,
                    required=False, default=DEFAULT_SNP_PROB,
                    dest="prob",
                    help="probability of SNP, default: %s" % DEFAULT_SNP_PROB)

args = parser.parse_args()


seq = "".join(random.choice(CHARACTERS) for i in xrange(args.len))
all_seqs = [seq]
while len(all_seqs) < args.num:
    new_seq = list(seq[:])
    for i in xrange(args.len):
        if utils.bernoulli_test(args.prob):
            new_seq[i] = random.choice(CHARACTERS)
    new_seq = "".join(new_seq)
    if new_seq != seq:
        all_seqs.append(new_seq)

print "\n".join(all_seqs)
