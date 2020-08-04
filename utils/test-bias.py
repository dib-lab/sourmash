#!/usr/bin/env python
# code by David Koslicki
# see https://github.com/dib-lab/sourmash/issues/187

import sourmash
import numpy as np
import matplotlib.pyplot as plt
import random
import screed

random.seed(1)

def kmers(seq, kk):
    for i in range(len(seq) - kk + 1):
        kmer = seq[i:i+kk]
        rc = screed.rc(kmer)
        if kmer < rc:
            yield kmer
        else:
            yield rc

n = 10000  # sequence length
ksize = 10  # k-mer length
h = 5000  # number of hashes in sketch
i_range = range(1, 50000, 100)  # range of intersection sizes
true_jaccards = np.zeros(len(i_range))
estimate_jaccards = np.zeros(len(i_range))
it = 0
for i_size in i_range:
    # Append a common string to two different random strings (true jaccard will be ~ i_size/n)
    common_string = ''.join(np.random.choice(['A', 'C', 'T', 'G'], i_size))
    seq1 = ''.join(np.random.choice(['A', 'C', 'T', 'G'], n)) + common_string
    seq2 = ''.join(np.random.choice(['A', 'C', 'T', 'G'], n)) + common_string

#    with open('xyz/{}-1.fa'.format(i_size), 'wt') as fp:
#        fp.write('>a\n{}\n'.format(seq1))
#    with open('xyz/{}-2.fa'.format(i_size), 'wt') as fp:
#        fp.write('>a\n{}\n'.format(seq2))

    # Calculate exact Jaccard index
    kmers1 = set(kmers(seq1, ksize))
    kmers2 = set(kmers(seq2, ksize))

    true_jaccard = len(kmers1.intersection(kmers2)) / float(len(kmers1.union(kmers2)))
    true_jaccards[it] = true_jaccard

    # Calculate sourmash estimate of Jaccard index
    E1 = sourmash.MinHash(n=h, ksize=ksize)
    E2 = sourmash.MinHash(n=h, ksize=ksize)
    E1.add_sequence(seq1)
    E2.add_sequence(seq2)
    estimate_jaccard = E1.jaccard(E2)
    estimate_jaccards[it] = estimate_jaccard
    it += 1

differences = true_jaccards - estimate_jaccards
sorted_true = sorted(true_jaccards)
sorted_estimates = np.array([x for (y, x) in sorted(zip(true_jaccards, estimate_jaccards), key=lambda pair: pair[0])])
sorted_differences = sorted_true - sorted_estimates
plt.figure()
plt.plot(sorted_true, sorted_differences)
axes = plt.gca()
axes.set_ylim([np.min(plt.yticks()[0])*1.5, np.max(plt.yticks()[0])*1.5])
plt.title('True - estimate Jaccard index')
plt.text(0, 0, 'Underestimate', rotation=90, horizontalalignment='center', verticalalignment='bottom', multialignment='center', color='b', fontsize=14)
plt.text(0, 0, 'Overestimate', rotation=90, horizontalalignment='center', verticalalignment='top', multialignment='center', color='r', fontsize=14)
plt.axhline(0, color='black', linestyle='dashed', linewidth=2)
plt.ylabel('Difference')
plt.xlabel('True Jaccard index')
plt.savefig('Differences.png')
print('created Differences.png')

plt.figure()
n, bins, patches = plt.hist(differences, 50, normed=1, facecolor='green', alpha=0.75)
plt.axvline(0, color='b', linestyle='dashed', linewidth=2)
plt.title('Histogram of (true - estimate) Jaccard index\n Mean: %f' % np.mean(differences))
plt.text(0, max(plt.yticks()[0])-1, 'Underestimate', rotation=0, horizontalalignment='left', verticalalignment='top', multialignment='left', color='b', fontsize=14)
plt.text(plt.xticks()[0][1], max(plt.yticks()[0])-1, 'Overestimate', rotation=0, horizontalalignment='left', verticalalignment='top', multialignment='left', color='r', fontsize=14)
plt.xlabel('Difference')
plt.savefig('Histogram.png')

print('created Histogram.png')
