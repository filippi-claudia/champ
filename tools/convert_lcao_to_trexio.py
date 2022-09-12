'''
This script can be used to convert the lcao
(all Cartesian but in the XXX,YYY,ZZZ etc order) and bfinfo files into the
newer lcao and basis pointers files which have AO ordering consisten with trexio.

Three files are required to run this script:
1) Geometry file in .xyz format or .geom (older format)
2) bfinfo file (Containing nsx,npx,npy,npz,ndxx... and index of columns of basis grid file)
3) lcao file ( in the XXX,YYY,ZZZ etc ordering)

Author: Ravindra Shinde
Email : r.l.shinde@utwente.nl
Date : September 06, 2022
'''




from os.path import splitext
import argparse
import numpy as np
np.set_printoptions(threshold=np.inf)
from collections import Counter
from collections.abc import Iterable
from itertools import count as icount

# Instantiate the parser
parser = argparse.ArgumentParser(description='Python Converter for conversion of all-Cartesian CHAMP-ordered MOs to trexio ordered MOs for the CHAMP code')

# Required positional argument
parser.add_argument("--lcao", "-s", "--orb", dest='filename_lcao', type=str, required = True,
                    help='Required: Filename (including extension) of .lcao or .orb file in the CHAMP AO-ordered format')

# Required positional argument
parser.add_argument("--bfinfo", "-r", "--bf", dest='filename_bfinfo', type=str, required = True,
                    help='Required: Filename (including extension) of .bfinfo file in the CHAMP format')

# Required positional argument
parser.add_argument("--geom", "-g", "--xyz", dest='filename_geom', type=str, required = True,
                    help='Required: Filename of .xyz or .geom file in the old or new format')

args = parser.parse_args()

print ("Filenames parsed are :")
print (' lcao file old   ::         \t {}'.format(args.filename_lcao))
print (' bfinfo file old ::         \t {}'.format(args.filename_bfinfo))
print (' geom file old   ::         \t {}'.format(args.filename_geom))


def flatten(x):
    if isinstance(x, Iterable):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]

### Read the lcao file first

with open(args.filename_lcao) as f1:
    lines = f1.readlines()

    iorb = 0; mocoeffs = []
    for line in lines:
        # Skip the comments
        if line.startswith('#'):
            continue
        if line.startswith('end'):
            continue
        # Read the number of basis and number of coeffs
        if line.startswith('lcao'):
            ncoeff = int(line.split()[1])
            nbasis = int(line.split()[2])
            continue

        temp = [float(i) for i in list(line.split()) ]
        mocoeffs.append(temp)
        iorb += 1




### Read the geometry file
with open(args.filename_geom) as f2:
    lines = f2.readlines()

    if any("&atoms" in s for s in lines):
        print ("")
        print ("Older format of Geometry file detected")
        dict_atom_type = {}
        atom_type = []
        for line_num, line in enumerate(lines):
            # print ("debug print  all lines ", line_num, line)
            # Skip the comments
            if line.startswith('#'):
                continue

            if line.startswith('&atoms'):
                nctype = int(line.split()[2])
                natoms = int(line.split()[4])
                continue
            if line.startswith('&atom_types'):
                ntokens = len(line.split()) - 1
                for i in range(0,ntokens,2):
                    dict_atom_type[int(line.split()[i+1])] = line.split()[i+2]
                continue

            if line.startswith('geometry'):
                coord_block_start = line_num + 1

        print ('Number of atoms             \t {}'.format(natoms))
        print ('Number of types of atoms    \t {}'.format(nctype))
        print ('Dictionary of atom types    \t {}'.format(dict_atom_type))

        nucleus_coord = np.zeros((natoms,3),dtype=float)
        for i in range(coord_block_start, coord_block_start+natoms):
            nucleus_coord[i-coord_block_start][0] = lines[i].split()[0]
            nucleus_coord[i-coord_block_start][1] = lines[i].split()[1]
            nucleus_coord[i-coord_block_start][2] = lines[i].split()[2]
            atom_type.append(lines[i].split()[3])

        atom_type_symbol = []
        for i in atom_type:
            atom_type_symbol.append(dict_atom_type[int(i)])

        print ('Atom Symbols                \t {}'.format(atom_type_symbol))
        print ('Atom types                  \t {}'.format(atom_type))

    else:
        print ("Newer format of geometry file detected")
        dict_atom_type = {}
        atom_type_symbol = []
        count = 0
        for line_num, line in enumerate(lines):
            if line.startswith('#'):
                continue
            # Read the first line
            ntokens = len(line.split())
            if ntokens == 1:
                natoms = int(line.split()[0])
                nucleus_coord = np.zeros((natoms,3),dtype=float)
                print ('Number of atoms             \t {}'.format(natoms))

            if ntokens == 4:
                nucleus_coord[count][0] = line.split()[1]
                nucleus_coord[count][1] = line.split()[2]
                nucleus_coord[count][2] = line.split()[3]
                atom_type_symbol.append(line.split()[0])
                count += 1


        unique_elements, indices = np.unique(atom_type_symbol, return_index=True)
        unique_atom_type = np.array(atom_type_symbol)[indices]
        nctype = len(unique_elements)

        # assignment of atom types
        for i in range(len(unique_atom_type)):
            dict_atom_type[i+1] = unique_atom_type[i]

        atom_type = []
        for i in atom_type_symbol:
            atom_type.append( list(unique_elements).index(i) + 1 )

        print ('Number of atoms             \t {}'.format(natoms))
        print ('Number of types of atoms    \t {}'.format(nctype))
        # print ('Dictionary of atom types    \t {}'.format(dict_atom_type))

        print ('Atom Symbols                \t {}'.format(atom_type_symbol))
        # print ('Atom types                  \t {}'.format(atom_type))


## utility functions
def get_key(dictionary, val):
    for key, value in dictionary.items():
        if val == value:
            return key


### Read the bfinfo file
with open(args.filename_bfinfo) as f3:
    lines = f3.readlines()

    for line_num, line in enumerate(lines):
        # Skip the comments
        if line.startswith('#'):
            continue
        if line.startswith('qmc_bf_info'):
            coord_block_start = line_num + 1
        if line.startswith('end'):
            continue

    atom_type_symbol = np.array(atom_type_symbol)
    Atoms, unique_atom_index = np.unique(atom_type_symbol, return_index=True)
    Atoms, count_each_type_atoms = np.unique(atom_type_symbol, return_counts=True)
    same_ordered_unique_atoms =  Atoms[np.argsort(unique_atom_index)]

    # print("debug ", Atoms)
    # print("same order ", same_ordered_unique_atoms)
    print ('Unique atoms             \t {}'.format(Atoms))
    print ('Indices of unique atoms  \t {}'.format(unique_atom_index))
    print ('count each atom type     \t {}'.format(count_each_type_atoms))


    num_unique_atoms = len(unique_atom_index)

    dict_num_per_shell = {} #only the odd numbered rows of data
    for i in range(coord_block_start, coord_block_start+2*num_unique_atoms,2):
        dict_num_per_shell[i] = lines[i].split()

    dict_radial_pointers = {} #only the even numbered rows of data
    for i in range(coord_block_start+1, coord_block_start+2*num_unique_atoms,2):
        dict_radial_pointers[i] = lines[i].split()


print (" ")
print (" ")
### Writing part starts here!

## Writing the new geometry file begins here ----------------
new_filename_geom = "champ_v3_trexio_order_" + splitext(args.filename_geom)[0] + ".xyz"
if new_filename_geom is not None:
    if isinstance(new_filename_geom, str):
        ## Write down a symmetry file in the new champ v2.0 format
        with open(new_filename_geom, 'w') as file:

            file.write("{} \n".format(natoms))
            # header line printed below
            file.write("# Converted from the old .geometry file to CHAMP v2 .xyz file \n")

            for element in range(natoms):
                file.write("{:5s} {} {} {} \n".format(atom_type_symbol[element], nucleus_coord[element][0], nucleus_coord[element][1], nucleus_coord[element][2]))

            file.write("\n")
        file.close()
    else:
        raise ValueError
# all the geometry file information written to the file
print ("New geometry file created   \t {} ".format(new_filename_geom))
## Writing the new geometry file ends here ----------------


## Writing the new bfinfo file begins here ----------------
## Get the bfinfo file with cartesian ordering of shell pointers
# Cartesian ordering
order = [[0],
         [0, 1, 2],
         [0, 1, 2, 3, 4, 5],
         [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
         [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]]

## v2.0 champ shell ordering (to be replaced by trexio order)
shells = {}
shells[0] = ['S']
shells[1] = ['X','Y','Z']
shells[2] = ['XX','XY','XZ','YY','YZ','ZZ']
shells[3] = ['XXX','XXY','XXZ','XYY','XYZ','XZZ','YYY','YYZ','YZZ','ZZZ']
shells[4] = ['XXXX','XXXY','XXXZ','XXYY','XXYZ','XXZZ','XYYY','XYYZ','XYZZ','XZZZ','YYYY','YYYZ','YYZZ','YZZZ','ZZZZ']

dict_l_of_shells = {1:0, 3:1, 6:2, 10:3}

## Write the new bfinfo file begins here ----------------
new_filename_bfinfo = "champ_v3_" + splitext(args.filename_bfinfo)[0] + "_basis_pointers.bfinfo"
if new_filename_bfinfo is not None:
    if isinstance(new_filename_bfinfo, str):
        ## Write down a bfinfo file in the new champ v3.0 format
        with open(new_filename_bfinfo, 'w') as file:

            # qmc bfinfo line printed below
            file.write("# Format of the new basis information file champ_v3 \n")
            file.write("# num_ao_per_center, n(s), n(p), n(d), n(f), n(g) \n")
            file.write("# Index of Slm (Range 1 to 35) \n")
            file.write("# Index of column from numerical basis file  \n")
            file.write("qmc_bf_info 1 \n")

            ## Calculate the needed quantities

            for vals in dict_radial_pointers.values():
                pointers = [int(i) for i in vals]
                basis_pointers = sorted(pointers)
                unique_shells, shell_count = np.unique(basis_pointers, return_counts=True)
                # sorted_pointers can be stored as it is now.

                num_shells_per_l = np.zeros(5, dtype=int)

                basis_per_atom = np.sum(shell_count)

                _, temp_counter = np.unique(shell_count, return_counts=True)

                list_slm_index = []
                for l in range(len(temp_counter)):
                    num_shells_per_l[l] = temp_counter[l]
                    for nshell in range(num_shells_per_l[l]):
                        if l == 0:
                            list_slm_index.append(1)
                        elif l == 1:
                            list_slm_index.append([2,3,4])
                        elif l == 2:
                            list_slm_index.append([5,6,7,8,9,10])
                        elif l == 3:
                            list_slm_index.append([11,12,13,14,15,16,17,18,19,20])
                        elif l == 4:
                            list_slm_index.append([21,22,23,24,25,26,27,28,29,30,31,32,33,34,35])

                slm_index = flatten(list_slm_index)

                # Things to write to the file
                # line 1) basis_per_atom, num_shells_per_l
                # line 2) slm_index
                # line 3) basis_pointers

                # pointers to the basis functions
                file.write(f"{basis_per_atom} ")

                for num in num_shells_per_l:
                    file.write(f"{num} ")
                file.write(f"\n")

                # Write the number of types of shells for each unique atom
                for slm in slm_index:
                    file.write(f"{slm} ")
                file.write(f"\n")

                # Write the pointers to the basis functions
                for pointer in basis_pointers:
                    file.write(f"{pointer} ")
                file.write(f"\n")
            file.write("end\n")
        file.close()

    else:
        raise ValueError
print ("New bfinfo file created   \t {} ".format(new_filename_bfinfo))
# all the bfinfo file information written to the file










### Operations on MO coefficients done here
## Calculate the needed quantities

dict_unique_atom_l = {}
dict_shellcount_per_l = {}

count_over_unique = 0
for vals in dict_radial_pointers.values():
    pointers = [int(i) for i in vals]
    basis_pointers = sorted(pointers)
    unique_shells, shell_count = np.unique(basis_pointers, return_counts=True)

    # sorted_pointers can be stored as it is now.

    num_shells_per_l = np.zeros(5, dtype=int)

    basis_per_atom = np.sum(shell_count)

    _, temp_counter = np.unique(shell_count, return_counts=True)

    list_slm_index = []
    list_shell_ang_mom = []
    for l in range(len(temp_counter)):
        num_shells_per_l[l] = temp_counter[l]
        for nshell in range(num_shells_per_l[l]):
            list_shell_ang_mom.append(l)


    dict_unique_atom_l[same_ordered_unique_atoms[count_over_unique]] = list_shell_ang_mom
    # print ("dict ang mom", same_ordered_unique_atoms[count_over_unique], dict_unique_atom_l )
    dict_shellcount_per_l[same_ordered_unique_atoms[count_over_unique]] = list(num_shells_per_l)
    count_over_unique += 1


list_ang_mom_l_all = [[] for _ in range(natoms)]
list_num_shell_per_l_all = [[] for _ in range(natoms)]

for i in range(natoms):
    list_ang_mom_l_all[i] = dict_unique_atom_l[atom_type_symbol[i]]
    list_num_shell_per_l_all[i] = dict_shellcount_per_l[atom_type_symbol[i]]

# print ("shell ang mom all atoms  ", list_ang_mom_l_all )
# print ("dict_shellcount_per_l", list_num_shell_per_l_all)

###


# This part is for reshuffling to make the AO basis in the CHAMP's own ordering
index_dict = {}; shell_representation = {}
new_shell_representation = []
counter = 0; basis_per_atom = []
ind = 0; champ_ao_ordering = []
for atom_index in range(natoms):
    bfcounter = 1; basis_per_atom_counter = 0
    pindex = 0; dindex = 0; findex = 0; gindex = 0
    for l in list_ang_mom_l_all[atom_index]:
        # run a small loop to reshuffle the shell ordering
        if l == 0:
            new_shell_representation.append(shells[l][0])
            champ_ao_ordering.append(ind)
            ind += 1

        local_p = np.zeros((3,list_num_shell_per_l_all[atom_index][1]),dtype='U1')
        local_ind_p = np.zeros((3,list_num_shell_per_l_all[atom_index][1]),dtype=int)
        if l == 1:
            pindex += 1; ind = champ_ao_ordering[-1] + 1
            for j in range(list_num_shell_per_l_all[atom_index][1]):
                #loop over all 3 p orbitals
                for k in order[l]:
                    local_p[k,j] = shells[l][k]
                    local_ind_p[k,j] = ind
                    ind += 1

            if pindex == list_num_shell_per_l_all[atom_index][1]:
                new_shell_representation.extend(list(local_p.flatten()))
                champ_ao_ordering.extend(list(local_ind_p.flatten()))



        local_d = np.zeros((6,list_num_shell_per_l_all[atom_index][2]),dtype='U2')
        local_ind_d = np.zeros((6,list_num_shell_per_l_all[atom_index][2]),dtype=int)
        if l == 2:
            dindex += 1; ind = champ_ao_ordering[-1] + 1
            for j in range(list_num_shell_per_l_all[atom_index][2]):
                #loop over all 6 d orbitals
                for k in order[l]:
                    local_d[k,j] = shells[l][k]
                    local_ind_d[k,j] = ind
                    ind += 1

            if dindex == list_num_shell_per_l_all[atom_index][2]:
                new_shell_representation.extend(list(local_d.flatten()))
                champ_ao_ordering.extend(list(local_ind_d.flatten()))


        local_f = np.zeros((10,list_num_shell_per_l_all[atom_index][3]),dtype='U3')
        local_ind_f = np.zeros((10,list_num_shell_per_l_all[atom_index][3]),dtype=int)
        if l == 3:
            findex += 1; ind = champ_ao_ordering[-1] + 1
            for j in range(list_num_shell_per_l_all[atom_index][3]):
                #loop over all 10 f orbitals
                for k in order[l]:
                    local_f[k,j] = shells[l][k]
                    local_ind_f[k,j] = ind
                    ind += 1

            if findex == list_num_shell_per_l_all[atom_index][3]:
                new_shell_representation.extend(list(local_f.flatten()))
                champ_ao_ordering.extend(list(local_ind_f.flatten()))

        local_g = np.zeros((15,list_num_shell_per_l_all[atom_index][4]),dtype='U4')
        local_ind_g = np.zeros((15,list_num_shell_per_l_all[atom_index][4]),dtype=int)
        if l == 4:
            gindex += 1; ind = champ_ao_ordering[-1] + 1
            for j in range(list_num_shell_per_l_all[atom_index][4]):
                #loop over all 15 g orbitals
                for k in order[l]:
                    local_g[k,j] = shells[l][k]
                    local_ind_g[k,j] = ind
                    ind += 1

            if gindex == list_num_shell_per_l_all[atom_index][4]:
                new_shell_representation.extend(list(local_g.flatten()))
                champ_ao_ordering.extend(list(local_ind_g.flatten()))


        # Get number of AO basis per atom
        for k in order[l]:
            shell_representation[counter] = shells[l][k]
            index_dict[counter] =  counter
            counter += 1
            basis_per_atom_counter += 1
        bfcounter += 1
    basis_per_atom.append(basis_per_atom_counter)

# print ("champ ao ordering: ", champ_ao_ordering)
# print ("new_shell_representation: ", new_shell_representation)
# print ("old_shell_representation: ", shell_representation.values())
# print ("basis per atom: ", basis_per_atom)



mocoeffs = np.array(mocoeffs)
## Reorder orbitals according to the ordering of the CHAMP ordering
reordered_mo_array =  mocoeffs[:,np.argsort(champ_ao_ordering)]


# write the transformed molecular coefficients to the new .lcao file
new_filename_orbitals = "champ_v3_trexio_order_" + splitext(args.filename_lcao)[0] + '.lcao'
if new_filename_orbitals is not None:
    if isinstance(new_filename_orbitals, str):
        ## Write down an orbitals file in the new champ v2.0 format
        with open(new_filename_orbitals, 'w') as file:

            # header line printed below
            file.write("# File created using the champ v3 converter in the trexio order \n")
            file.write("lcao " + str(ncoeff) + " " + str(nbasis) + " 1 " + "\n" )
            # np.savetxt(file, reordered_mo_array, fmt='%.8f')
            np.savetxt(file, reordered_mo_array)
            file.write("end\n")
        file.close()
    else:
        raise ValueError
# all the lcao file information written to the file
print ("New orbital file created   \t {} ".format(new_filename_orbitals))

print ('{}'.format("All files have been converted successfully"))
# The end