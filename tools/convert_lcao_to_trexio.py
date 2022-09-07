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
        print ('Dictionary of atom types    \t {}'.format(dict_atom_type))

        print ('Atom Symbols                \t {}'.format(atom_type_symbol))
        print ('Atom types                  \t {}'.format(atom_type))


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

    # unique_atoms, indices = np.unique(atom_type_symbol, return_index=True)
    _, indices = np.unique(atom_type_symbol, return_index=True)
    unique_atoms = np.array(atom_type_symbol)[np.sort(indices)]
    _, count_each_type_atoms = np.unique(atom_type_symbol, return_counts=True)
    # unique_atoms = unique_atoms[np.sort(indices)]
    print ('unique elements             \t {}'.format(unique_atoms))
    print ('indices                     \t {}'.format(np.sort(indices)))

    print ('count each atom type        \t {}'.format(count_each_type_atoms))
    num_unique_atoms = len(unique_atoms)

    dict_num_per_shell = {} #only the odd numbered rows of data
    for i in range(coord_block_start, coord_block_start+2*num_unique_atoms,2):
        dict_num_per_shell[i] = lines[i].split()

    print ("dict num per shell", dict_num_per_shell)


    dict_radial_pointers = {} #only the even numbered rows of data
    for i in range(coord_block_start+1, coord_block_start+2*num_unique_atoms,2):
        dict_radial_pointers[i] = lines[i].split()

    print ("radial pointers", dict_radial_pointers)


print (" ")
print (" ")
### Writing part starts here!

## Get the bfinfo file with cartesian ordering of shell pointers
# Cartesian ordering
order = [[0],
         [0, 1, 2],
         [0, 1, 2, 3, 4, 5],
         [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]]

## NEW champ shell ordering
shells = {}
shells[0] = ['S']
shells[1] = ['X','Y','Z']
shells[2] = ['XX','XY','XZ','YY','YZ','ZZ']
shells[3] = ['XXX','XXY','XXZ','XYY','XYZ','XZZ','YYY','YYZ','YZZ','ZZZ']

## OLD champ shell ordering
old_shells = {}
old_shells[0] = ['S']
old_shells[1] = ['X','Y','Z']
old_shells[2] = ['XX','YY','ZZ','XY','XZ','YZ']
old_shells[3] = ['XXX','YYY','ZZZ','XXY','XXZ','YYX','YYZ','ZZX','ZZY','XYZ']


dict_l_of_shells = {1:0, 3:1, 6:2, 10:3}
index_unique_atom = 0
basis_shell_ang_mom_unique = {}
dict_shell_counter = {}

## Write the new bfinfo file begins here ----------------
new_filename_bfinfo = "champ_v3_trexio_order_" + splitext(args.filename_bfinfo)[0] + ".bfinfo"
if new_filename_bfinfo is not None:
    if isinstance(new_filename_bfinfo, str):
        ## Write down a symmetry file in the new champ v2.0 format
        with open(new_filename_bfinfo, 'w') as file:

            # qmc bfinfo line printed below
            file.write("# Format of the new basis information file champ_v3 \n")
            file.write("# num_ao_per_center, n(s), n(p), n(d), n(f), n(g) \n")
            file.write("# Index of Slm (Range 1 to 35) \n")
            file.write("# Index of column from numerical basis file  \n")
            file.write("qmc_bf_info 1 \n")

            ## Calculate the needed quantities

            index_unique_atom = 0
            basis_shell_ang_mom_unique = {}
            dict_shell_counter = {}

            for vals in dict_radial_pointers.values():
                pointers = [int(i) for i in vals]
                basis_pointers = sorted(pointers)
                unique_shells, shell_count = np.unique(basis_pointers, return_counts=True)
                # sorted_pointers can be stored as it is now.

                num_shells_per_l = np.zeros(5, dtype=int)

                basis_per_atom = np.sum(shell_count)

                _, temp_counter = np.unique(shell_count, return_counts=True)

                for l in range(len(temp_counter)):
                    num_shells_per_l[l] = temp_counter[l]

                slm_index = []
                for ind in range(len(shell_count)):
                    val = shell_count[ind]
                    count = 1
                    for j in range(val):
                        if val == 1:
                            slm_index.append(count)
                            count += 1
                        else:
                            slm_index.append(count+val-2)
                            count += 1

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

## Write the new geometry file begins here ----------------
new_filename_geom = "champ_v3_trexio_order_" + splitext(args.filename_geom)[0] + ".xyz"
if new_filename_geom is not None:
    if isinstance(new_filename_geom, str):
        ## Write down a symmetry file in the new champ v2.0 format
        with open(new_filename_geom, 'w') as file:

            file.write("{} \n".format(natoms))
            # header line printed below
            file.write("# Converted from the old .geometry file to CHAMP v2 .xyz file \n")

            for element in range(natoms):
                file.write("{:5s} {: 0.6f} {: 0.6f} {: 0.6f} \n".format(atom_type_symbol[element], nucleus_coord[element][0], nucleus_coord[element][1], nucleus_coord[element][2]))

            file.write("\n")
        file.close()
    else:
        raise ValueError
# all the geometry file information written to the file
print ("New geometry file created   \t {} ".format(new_filename_geom))




### Operations on MO coefficients done here

## CHAMP old AO ordering (where d is moved to s orbital)
# ao_order_old= ['S','XX','X','Y','Z','YY','ZZ','XY','XZ','YZ','XXX','YYY','ZZZ','XXY','XXZ','YYX','YYZ','ZZX','ZZY','XYZ']
## Cartesian Ordering CHAMP
# ao_order_new = ['S',[all 'X'],[all 'Y'],[all 'Z'],'XX','XY','XZ','YY','YZ','ZZ','XXX','XXY','XXZ','XYY','XYZ','XZZ','YYY','YYZ','YZZ','ZZZ']


# ## Reorder orbitals according to the ordering of the CHAMP ordering
# reordered_mo_array = dict_mo["coefficient"][:,champ_ao_ordering]


# print ("reconstruting the old shell representation first atom")
# old_shell_reprensentation_per_atom = []
# for element in atom_type_symbol:
#     temp_old_shell = []
#     for ind, num in enumerate(dict_shell_counter[element]):
#         # S with XX of D
#         if ind == 0:
#             for j in range(num):
#                 temp_old_shell.append(old_shells[ind][0])

#             for jj in range(dict_shell_counter[element][2]):
#                 temp_old_shell.append(old_shells[2][0]+str(jj))

#         # P
#         if ind == 1:
#             temp = np.zeros((3,num),dtype='U2')
#             for j in range(num):
#                 for k in order[1]:
#                     temp[k,j] = old_shells[1][k] + str(j)

#             temp_old_shell.extend(list(temp.flatten()))

#         # D without XX
#         if ind == 2:
#             temp = np.zeros((6,num),dtype='U3')
#             for j in range(num):
#                 for k in order[2][1:]:  #exclude XX
#                     temp[k,j] = old_shells[2][k] + str(j)

#             temp_old_shell.extend(list(temp.flatten()))

#         # F
#         if ind == 3:
#             temp = np.zeros((10,num),dtype='U4')
#             for j in range(num):
#                 for k in order[3]:
#                     temp[k,j] = old_shells[3][k] + str(j)

#             temp_old_shell.extend(list(temp.flatten()))

#     temp_old_shell = list(filter(None, temp_old_shell))
#     old_shell_reprensentation_per_atom.append(temp_old_shell)



# # Create a copy of mocoeffs; only the d coeffs will change; rest will be shuffled
# transformed_mocoeffs = np.asarray(mocoeffs)

# # Get the list of indices of starting point for each atom
# basis_start_index = np.cumsum(basis_per_atom)
# basis_start_index = np.insert(basis_start_index, 0, 0, axis=0)
# # Ignore the last number from the above list


# summ = 0; final_list_indices = []
# for index, element in enumerate(atom_type_symbol):
#     temporary = []
#     B = np.array(distinct_shell_reprensentation_per_atom[index])
#     A = np.array(old_shell_reprensentation_per_atom[index])

#     #take care of repeated number corresponding to S shells
#     c = Counter(B)
#     iters = {k: icount(1) for k, v in c.items() if v > 1}
#     indexed_B = [x+str(next(iters[x])) if x in iters else x for x in B]
#     indexed_B =np.array(indexed_B)


#     c = Counter(A)
#     iters = {k: icount(1) for k, v in c.items() if v > 1}
#     indexed_A = [x+str(next(iters[x])) if x in iters else x for x in A]
#     indexed_A =np.array(indexed_A)

#     # print ("indexed A ", indexed_A)
#     # print ("indexed B ", indexed_B)

#     xsorted = np.argsort(indexed_B)
#     res = xsorted[np.searchsorted(indexed_B[xsorted], indexed_A)]

#     temporary = [str(a + summ) for a in res]
#     final_list_indices.extend(temporary)
#     summ = summ + basis_per_atom[index]


#     ## Do the molecular coefficient transformation here.
#     # Select the XX, YY, ZZ coefficients from the old file for each atom
#     # and then apply the following transformation::
#     CS = np.sqrt(5.0)
#     CD = 1.0/np.sqrt(3.0)
#     # a=oldcoeff('XX')
#     # b=oldcoeff('YY')
#     # c=oldcoeff('ZZ')
#     # newcoeff('XX') = a/CS - b/2.0 + c / (2.0*CD);
#     # newcoeff('YY') = a/CS - b/2.0 - c / (2.0*CD);
#     # newcoeff('ZZ') = a/CS + b;

#     # make sure that you loop over all the d coeffs
#     # print ("to search in list A", A)
#     # print ("basis per atom ", basis_per_atom)
#     for num_d in range(dict_shell_counter[element][2]):
#         index_xx = np.where( A == 'XX' + str(num_d))[0][0] + basis_start_index[index]
#         index_yy = np.where( A == 'YY' + str(num_d))[0][0] + basis_start_index[index]
#         index_zz = np.where( A == 'ZZ' + str(num_d))[0][0] + basis_start_index[index]
#         # print ("index what i want xx,yy,zz", index_xx,index_yy,index_zz)
#         for iorb in range(ncoeff):
#             a = mocoeffs[iorb][index_xx]
#             b = mocoeffs[iorb][index_yy]
#             c = mocoeffs[iorb][index_zz]
#             transformed_mocoeffs[iorb][index_xx] = a/CS - b/2.0 + c/(2.0*CD)
#             transformed_mocoeffs[iorb][index_yy] = a/CS - b/2.0 - c/(2.0*CD)
#             transformed_mocoeffs[iorb][index_zz] = a/CS + b




# flat_list = [item for sublist in old_shell_reprensentation_per_atom for item in sublist]
# #Convert to numpy array
# final_list_indices = np.asarray(final_list_indices)
# #Convert to numpy int array
# final_list_indices = final_list_indices.astype(int)
# # print ("list of indices final ", final_list_indices)
# ## Rearrange the transofrmed mocoeffs array with index array from final_list_indices
# argindex = final_list_indices.argsort()
# for iorb in range(ncoeff):
#     transformed_mocoeffs[iorb] = transformed_mocoeffs[iorb][argindex]

# # write the transformed molecular coefficients to the new .lcao file
# new_filename_orbitals = "champ_v3_trexio_order_" + splitext(args.filename_lcao)[0] + '_orbitals.lcao'
# if new_filename_orbitals is not None:
#     if isinstance(new_filename_orbitals, str):
#         ## Write down an orbitals file in the new champ v2.0 format
#         with open(new_filename_orbitals, 'w') as file:

#             # header line printed below
#             file.write("# File created using the old champ spherical to cartesian converter \n")
#             file.write("lcao " + str(ncoeff) + " " + str(nbasis) + " 1 " + "\n" )
#             np.savetxt(file, transformed_mocoeffs, fmt='%.8f')
#             file.write("end\n")
#         file.close()
#     else:
#         raise ValueError
# # all the lcao file information written to the file
# print ("New orbital file created   \t {} ".format(new_filename_orbitals))

print ('{}'.format("All files have been converted successfully"))
# The end