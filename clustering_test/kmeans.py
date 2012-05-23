#!/usr/bin/python

import random
import utils
from operator import itemgetter

#
# parameters
#
MISREAD_IDX_TO_BASE = dict(enumerate("NATGC"))
MISREAD_BASE_TO_IDX = dict((b, a) for a, b in enumerate("NATGC"))
MISREAD_PROBS = [
                 [.96, .01, .01, .01, .01], # from N
                 [.01, .9 , .04, .03, .02], # from A
                 [.01, .04, .9 , .02, .03], # from T
                 [.01, .03, .02, .9 , .04], # from G
                 [.01, .02, .03, .04, .9 ]  # from C
                 ]

READS = ["___ATGC",
         "GCGCTGC",
         "GCGCTAC",
         "GCGCT__",
         "TATATAT",
         "TCTAT__",
         "TATGTAT"]

K = 2

#
# functions
#
def distance(a, b):
    assert(len(a) == len(b))
    d = 1.0
    for i in xrange(len(a)):
        if a[i] == "_" or b[i] == "_":
            d *= 1.0 # _ can be anything
        else:
            a_idx = MISREAD_BASE_TO_IDX[a[i]]
            b_idx = MISREAD_BASE_TO_IDX[b[i]]
            d *= MISREAD_PROBS[a_idx][b_idx]
    # now d contains probability that a was generated from b (or vice versa)
    # distance is the opposite probability
    return 1-d

def min_distance(centroids, a):
    return min(distance(a, c) for c in centroids)

def min_distance_idx(centroids, a):
    distances = [distance(a, c) for c in centroids]
    return min(enumerate(distances), key=itemgetter(1))[0]

def normalize(xs):
    s = sum(xs)
    return [x / s for x in xs]

## if centroid contains some missed bases, we should complement it
def prepare_centroid(centroid, reads):
    if not any(c == "_" for c in centroid):
        return centroid # centroid doesn't have missed bases, return original one
    centroid_lst = list(centroid)
    probs_nonnorm = [1. - distance(centroid, read) for read in reads]
    for i in xrange(len(centroid)):
        if centroid[i] != "_":
            continue
        # we should drop reads that contain "_" in this position
        probs_nonnorm_copy = probs_nonnorm[:] # copying an array
        for read_idx in xrange(len(reads)):
            if reads[read_idx][i] == "_":
                probs_nonnorm_copy[read_idx] = 0.0
        if sum(probs_nonnorm_copy) == 0.0:
            continue # all reads contain _ in this position, despair
        donor_idx = utils.loaded_dice(normalize(probs_nonnorm_copy))
        centroid_lst[i] = reads[donor_idx][i]
    return "".join(centroid_lst)

def kpp(k, reads):
    centroids = [prepare_centroid(random.choice(reads), reads)]
    while len(centroids) < k:
        sqrd_distances = [min_distance(centroids, read) ** 2 for read in reads]
        sqrd_distances_sum = sum(sqrd_distances)
        probs = [d / sqrd_distances_sum for d in sqrd_distances]
        new_centroid_idx = utils.loaded_dice(probs)
        new_centroid = prepare_centroid(reads[new_centroid_idx], reads)
        if new_centroid in centroids:
            continue # ensure that we don't have similar centroids
        centroids.append(new_centroid)
    return centroids

# Critical part!
def find_centroid(reads):
    assert(len(reads) > 0)
    centroid = []
    for i in xrange(len(reads[0])):
        # AFAIK it's MAP
        probs = [1.0, 1.0, 1.0, 1.0, 1.0] # priors (unnormalized)
        for read in reads:
            base = read[i]
            if base == "_": continue
            base_idx = MISREAD_BASE_TO_IDX[base]
            for centroid_base_idx in xrange(len(probs)):
                observed = MISREAD_PROBS[base_idx][centroid_base_idx]
                probs[centroid_base_idx] *= observed
        most_likely_idx = max(enumerate(probs), key=itemgetter(1))[0]
        centroid.append(MISREAD_IDX_TO_BASE[most_likely_idx])
    return "".join(centroid)

def main():
    #
    # init
    #
    centroids = kpp(K, READS)

    #
    # main loop
    #
    prev_centroids = centroids
    first_pass = True
    iteration = 1
    while (first_pass or prev_centroids != centroids) and iteration < 1000:
        first_pass = False
        prev_centroids = centroids

        assignments = [min_distance_idx(centroids, read) for read in READS]
        groups = [[] for _ in centroids]
        for idx, assignment in enumerate(assignments):
            groups[assignment].append(READS[idx])
        centroids = [find_centroid(reads) for reads in groups]
        centroids = sorted(centroids) # handy for debugging

        print "iteration %s completed" % iteration
        iteration += 1

    print "centroids:"
    print centroids
    print "groups:"
    for g in groups: print g

if __name__ == "__main__":
    main()
