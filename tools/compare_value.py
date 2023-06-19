#!/usr/bin/env python
import itertools
import argparse
import math

def read_file(filename):
    f = open(filename)
    # data = itertools.ifilter(predicate, f)
    data = f.readlines()
    f.close()
    return data


def compare_values(data, keyword, ref_values, tolerance=0.0, kw_col=0, val_col=None):
    keyword = keyword.split()
    if val_col is None:
        val_col = kw_col+len(keyword)+1

    for l in data[::-1]:
        l = l.split()
        if l[kw_col:kw_col+len(keyword)] == keyword:
            standard_deviation = float(l[val_col+2])
            print(' Obtained value   :: ', float(l[val_col]), '+-' , float(l[val_col+2]) )
            print(' Reference value  :: ', float(ref_values), '+-' , float(tolerance))
            print('The results should be within +- 2*sqrt(sd**2 + sd**2_ref)')
            return ( abs(float(l[val_col]) - float(ref_values)) <= (2.0*math.sqrt( float(tolerance)*float(tolerance) + standard_deviation*standard_deviation ) ))
    print('Warning : keyword not found')
    return False


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="compare values in a file")
    parser.add_argument("filename", help="name of the output file")
    parser.add_argument("keyword", help="name of the variable")
    parser.add_argument("values",  help="reference value")
    parser.add_argument("tolerance", nargs="?", default=0.0, help="allowed tolerance in the value")
    parser.add_argument(
        "--no_assert", action='store_true', help="skip the assertio and only prints warning")
    args = parser.parse_args()

    data = read_file(args.filename)
    if args.no_assert:
        test = compare_values(data, args.keyword, args.values, args.tolerance)
        if not test:
            print('Value Error : Reference and test values are different')
    else:
        # Strict comparison :: use assert
        assert(compare_values(data, args.keyword, args.values, args.tolerance))
        print ("-------------------------------------------------------------")
        #compare_values(data, args.keyword, args.values, args.tolerance) # delete after testing. Ravindra
