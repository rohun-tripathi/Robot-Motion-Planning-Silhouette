import sys

import numpy as np

import shared as SH


def processInput(inputfile):
    origin_list = []
    ellipse_list = []
    inp = open(inputfile, "r").readlines()  # Write the code to validate the input files
    for index, l in enumerate(inp):
        parts = l.strip().split()
        if index == 0:
            SH.num = len(inp) - 1  # Number of ellipses = Number of lines - line for dim
            SH.dim = int(parts[0])  # dim is for the dimensions being used
            continue
        else:
            if SH.dim == 0:
                print "error in input as dim = 0"
                sys.exit()
            origin = []  # first detail in line is origin
            for i in range(SH.dim):
                origin.append(float(parts[i]))
            origin_list.append(origin)  # store the center values for the nth ellipse
            carry = SH.dim;  # carry tracks how many have been read and where to read from next

            Umatrix = []  # holds the U matrix in 2d list form
            for i in range(SH.dim):
                Ulist = []  # List for input
                for j in range(SH.dim):
                    Ulist.append(float(parts[j * SH.dim + i + carry]))
                Umatrix.append(Ulist)

            Sigmalist = []  # holds the sigma matrix
            carry = SH.dim + SH.dim * SH.dim  # Carry gets updated after reads
            for i in range(SH.dim):  # Calculations for A of matrix
                term = float(parts[i + carry])
                inver = 1 / (term * term)
                Sigmalist.append(inver)
                if index == 1:
                    SH.primA.append(term)

            eye = np.identity(SH.dim, dtype=float)
            for i in range(SH.dim):
                eye[i][i] = Sigmalist[i]  # now eye stores the SIGMA matrix
            U = np.array(Umatrix, float)
            trans = U.T
            # print "U and trans = " , U, trans		#multiply the three
            mult1 = np.dot(U, eye)
            mult2 = np.dot(mult1, trans)  # A = U * eye * U.T
            # print "A for this iteration = ", mult2
            ellipse_list.append(mult2)
    return origin_list, ellipse_list
