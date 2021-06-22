import numpy as np
from pytest import approx


linecount = 0
first_match = []

with open('vmc_only.out', 'r') as f:
    lines = (line.rstrip() for line in f)
    lines = list(line for line in lines if line)

    for i in lines:
        linecount =  linecount + 1
        if (i.split()[0]) == 'enow':
            first_match.append(linecount)

    first_energy = (lines[first_match[0]].split()[0])
    print ("First Energy: " , first_energy)

    # Check the first energy is same or not
    assert (float(first_energy) == approx(-26.06833, rel=1e-6, abs=1e-6))
