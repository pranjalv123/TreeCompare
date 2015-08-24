import dendropy
import sys
import pylab as pl
import argparse

def compare(reference, test, metric):
    rb = set(reference.encode_bipartitions())
    tb = set(test.encode_bipartitions())

    l = []
    s = []
    
    for i in tb:
        if i not in rb:
            l.append(i)
            s.append( metric(reference, test, i))

    return l, s
    

def small_clade_size_metric(r, t, b):
    n = bin(b.leafset_bitmask).count('1')
    return min(n, len(t.taxon_namespace) - n)

def rf_matric(r, t, b):
    return 1


if __name__ == "__main__":
    t = dendropy.Tree.get_from_path(sys.argv[1], 'newick')
    r = dendropy.Tree.get_from_path(sys.argv[2], 'newick')
    l, s = compare(r, t, small_clade_size_metric)
    pl.hist(s, bins = 1000)
    count = 0
    for i in zip(l, s):
        if i[1] > 100:
            print i
            count += 1
    print count
    pl.show()
