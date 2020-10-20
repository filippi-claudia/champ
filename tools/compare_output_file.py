#!/usr/bin/env python
import itertools
import argparse


def predicate(line):
    if 'END' in line:
        return False
    return True


def read_file(filename):
    f = open(filename)
    # data = itertools.ifilter(predicate, f)
    data = f.readlines()
    f.close()
    return data


def compare_file(data1, data2):
    for l1, l2 in zip(data1, data2):
        if 'END' in l1:
            continue
        if l1 != l2:
            print('+', l1)
            print('-', l2)
            print('')
            # return False
    return True


def compare_values(data, keyword, ref_values, kw_col=0, val_col=1):
    for l in data:
        if data[kw_col] == keyword:
            assert(float(data[val_col]) == ref_values)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="compare two files")
    parser.add_argument("filename1", help="name of the first file")
    parser.add_argument("filename2", help="name of the second file")
    parser.add_argument("-kw", "--keyword", default=None,
                        help="name of the variable to test")
    parser.add_argument(
        "-v", "--values", default=None, help="reference value")
    args = parser.parse_args()

    data1 = read_file(args.filename1)
    data2 = read_file(args.filename2)

    assert(compare_file(data1, data2))
