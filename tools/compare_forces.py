#!/usr/bin/env python3
import argparse
import math

def read_file(filename):
    data = [[float(i) for i in line.split()][1:] for line in open(filename)]
    return data

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="compare values in a file")
    parser.add_argument("filename1", help="name of the output file")
    parser.add_argument("filename2", help="name of the variable")
    args = parser.parse_args()

    data1 = read_file(args.filename1)
    data2 = read_file(args.filename2)
    N_atoms = len(data1)
    correct = True
    for ic in range(N_atoms):
        for k in range(3):
            Ediff = data1[ic][k] - data2[ic][k]
            error = 2*math.sqrt(data1[ic][k+3]**2 + data2[ic][k+3]**2)
            if not abs(Ediff) < abs(error):
                correct = False

    assert correct, "Analytic forces not within errorbar"
