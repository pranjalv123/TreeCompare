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

def rf_metric(r, t, b):
    return 1

def no_normalization(r, t):
    return 1


def ref_normalization(r, t):
    return len(r.edges())

def test_normalization(r, t):
    return len(t.edges())

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare trees")

    parser.add_argument('-m', '--metric', dest='metric', default='rf')
    parser.add_argument('-g', '--show-hist', dest='hist', action='store_true')
    parser.add_argument('-n', '--normalize', dest='normalize', default='ref')
    parser.add_argument('reference')
    parser.add_argument('test')

    args = vars(parser.parse_args())

    tn = dendropy.TaxonNamespace()
    
    t = dendropy.Tree.get_from_path(args['test'], 'newick', taxon_namespace = tn)
    r = dendropy.Tree.get_from_path(args['reference'], 'newick', taxon_namespace = tn)

    l, s = compare(r, t, globals()[args['metric'] + '_metric'])
    print float(sum(s))/float(globals()[(args['normalize'] + '_normalization')](r, t))
    if args['hist']:
        pl.hist(s, bins = 1000)
        count = 0
        for i in zip(l, s):
            if i[1] > 100:
                print i
                count += 1
        print count
        pl.show()
