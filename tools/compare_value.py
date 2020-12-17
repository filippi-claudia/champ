#!/usr/bin/env python
import itertools
import argparse


def read_file(filename):
    f = open(filename)
    # data = itertools.ifilter(predicate, f)
    data = f.readlines()
    f.close()
    return data


def compare_values(data, keyword, ref_values, kw_col=0, val_col=None):
    keyword = keyword.split()
    if val_col is None:
        val_col = kw_col+len(keyword)+1

    for l in data[::-1]:
        l = l.split()
        if l[kw_col:kw_col+len(keyword)] == keyword:
            print(float(l[val_col]), float(ref_values))
            return float(l[val_col]) == float(ref_values)
    print('Warning : keyword not found')
    return False


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="compare values in a file")
    parser.add_argument("filename", help="name of the output file")
    parser.add_argument("keyword", help="name of the variable")
    parser.add_argument("values",  help="reference value")
    args = parser.parse_args()

    data = read_file(args.filename)
    assert(compare_values(data, args.keyword, args.values))
