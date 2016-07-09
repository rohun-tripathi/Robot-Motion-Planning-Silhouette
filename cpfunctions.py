import pprint as pprint
import sys
from copy import deepcopy

import numpy as np
import numpy.linalg as linalg

import shared as SH


def ReduceSingleEllipsoid(arrayYZ, CPslicevector, travaxis, ellipse, origin, debug=False):
    topcorner = ellipse[0][0]
    planeState = CPslicevector[travaxis]
    highercenter = origin[0]
    # denom = 1 - aii (a - ci ) ^ 2
    denominator = 1 - topcorner * (planeState - highercenter) * (planeState - highercenter)

    if debug: print "topcorner, planeState, highercenter == ", topcorner, planeState, highercenter
    if debug: print "denominator == ", denominator, "\n"

    firstcenter = np.array(origin[1:])
    trans = firstcenter.T
    oldmatrix = np.array([x[1:] for x in ellipse[1:]])
    if debug: print "Firstcenter, it's trans, oldmatrix == ", firstcenter, trans, oldmatrix

    value = np.dot(np.dot(trans, oldmatrix), firstcenter)
    denominator = denominator - value

    if debug: print value
    if debug: print "denominator == ", denominator, "\n"

    firstrowA = np.array(ellipse[0][1:])
    trans = firstrowA.T
    denominator = denominator + 2 * (planeState - highercenter) * np.dot(trans, firstcenter)

    if debug: print "denominator == ", denominator, "\n"

    center2nd = np.array(arrayYZ)
    trans = center2nd.T
    # numerator = 1 - np.dot(trans, center2nd) #This is the mistake

    arrayofD = []

    reduceEllipseDim = len(ellipse[0]) - 1  # n-1

    if len(center2nd) != reduceEllipseDim:
        print "Error! len(center2nd) != reduceEllipseDim\nExiting"
        sys.exit()

    for irow in range(reduceEllipseDim):
        for icol in range(reduceEllipseDim):
            arrayofD.append(center2nd[irow] * center2nd[icol])

    if debug: print "arrayofD == ", arrayofD

    # begins the preparation for equation. AX = Y
    # X in the equations is the terms of the "B" metrix (for the ellipse in the lower dim)  len -> (n-1) * (n-1)
    # A is the coefficients of the terms in the "B" matrix
    # Y is the constant term of each equation

    A = [];
    Y = []

    for irow in range(reduceEllipseDim):
        for icol in range(reduceEllipseDim):
            apq = ellipse[1 + irow][1 + icol]
            Y.append(apq)  # Accounting for the lowering in dimensions
            mult = np.multiply(apq, arrayofD)
            mult[(irow) * reduceEllipseDim + icol] += denominator
            A.append(mult.tolist())

    if debug:
        print "A == "
        pprint.pprint(A)

    if debug:
        print "here"

    X = np.linalg.solve(A, Y)

    if debug:
        print "X == "
        pprint.pprint(X)

    # produce the newmatrix

    newmatrix = []
    for irow in range(reduceEllipseDim):
        temp = []
        for icol in range(reduceEllipseDim):
            temp.append(X[irow * reduceEllipseDim + icol])
        newmatrix.append(temp)

    # tempellise = deepcopy(ellipse)

    # newmatrix = [x[1:] for x in tempellise[1:]  ]
    # if debug: print "Original	 newmatrix is == ",newmatrix
    # for col in range( len(newmatrix)):
    # 	for row in range( len(newmatrix[0]) ):
    # 		newmatrix[row][col] = newmatrix[row][col] * numerator/denominator
    # center2nd = arrayYZ

    # if debug: print "Obtained newmatrix is == ",newmatrix
    # if debug: print center2nd

    return newmatrix, center2nd


def ReduceEllipsoids(considerlist, considerYZ, CPslicevector, travaxis, ellipselist, originlist, debug=False):
    if debug: print "considerlist == ", considerlist
    if debug: print "considerYZ == ", considerYZ
    RecursionEllipses = []
    RecursionOrigins = []

    # The first (outer) Ellipse
    arrayYZ = CPslicevector[(travaxis + 1):]
    ellipse = deepcopy(ellipselist[0])
    origin = originlist[0]

    if debug: print "Original ellipse and origin and CPslicevector == ", ellipse, origin, CPslicevector

    newmatrix, center2nd = ReduceSingleEllipsoid(arrayYZ, CPslicevector, travaxis, ellipse, origin, True)
    RecursionEllipses.append(newmatrix[:])
    RecursionOrigins.append(center2nd[:])

    print considerlist

    for index, term in enumerate(considerlist):
        if term == 1:
            arrayYZ = considerYZ[index]

            ellipse = ellipselist[index + 1]  # Accounting for the outer ellipse
            origin = originlist[index + 1]

            newmatrix, center2nd = ReduceSingleEllipsoid(arrayYZ, CPslicevector, travaxis, ellipse, origin, True)

            RecursionEllipses.append(newmatrix[:])
            RecursionOrigins.append(center2nd[:])
    print "RecursionEllipses, RecursionOrigins == ", RecursionEllipses, RecursionOrigins
    # raw_input("Press Enter to continue.... ")
    return RecursionEllipses, RecursionOrigins


def retrieveCriticalPointsFromHigherLevel(critical, presentslice, nextslice, debug=False):
    CritAtThisSlice = []
    if debug: print "In RecurCheck and critical, presentslice, nextslice == ", critical, presentslice, nextslice
    for item in critical:
        # The < and <= essentially means that the obstacle cannot be sticking to the first point in the arena
        if presentslice < item[0][0] <= nextslice:
            CritAtThisSlice.append(item)
    return CritAtThisSlice


def retrieveActiveCPsForSlice(CriticalX, CriticalYZ, presentslice, nextslice, debug=False):
    activeCPsForSlice = []
    for tupl, rest in zip(CriticalX, CriticalYZ):
        if presentslice < tupl[0] < nextslice:
            criticalPoint = [tupl[0]] + rest[0]
            activeCPsForSlice.append([criticalPoint, "start"])
        elif presentslice < tupl[1] < nextslice:
            criticalPoint = [tupl[1]] + rest[1]
            activeCPsForSlice.append([criticalPoint, "end"])
    return activeCPsForSlice


def retrieveCPValuesAlongLowerDim(activeCPsForSlice):
    CPValuesForRemaningDim = []
    for cp in activeCPsForSlice:
        CPValuesForRemaningDim.append([cp[0][1:], cp[1]])
    return CPValuesForRemaningDim


def ellipseUnderConsider(considerlist, CriticalYZ, CriticalX, nextslice, travaxis, debug=False):
    considerYZ = []
    for index, listitem in enumerate(considerlist):
        if listitem == 0:
            considerYZ.append([])
            continue
        else:
            templist = [0 for x in range(travaxis + 1, SH.dim)]
            # Xvalues = (x - x0) / (x1 - x0)
            Xvalues = (nextslice - CriticalX[index][0]) / (CriticalX[index][1] - CriticalX[index][0])

            for dimen in range(SH.dim - (travaxis + 1)):
                templist[dimen] = CriticalYZ[index][1][dimen] - CriticalYZ[index][0][dimen]  # Val = (y1 - y0)
                templist[dimen] *= Xvalues  # Val = (x - x0) * (y1 - y0) / (x1 - x0)
                templist[dimen] += CriticalYZ[index][0][dimen]  # Val = (x - x0) * (y1 - y0) / (x1 - x0) + y0

            considerYZ.append(templist)
    return considerYZ


############################
# CPcalculate - Calculates the Critical Points, for all but the primA ellipse
# The returned list - Criticalpts, is a list of length n-1, where n is the ellipses being considered for this slice
# Create separate lists for the travaxis and other axes, both of length n-1
############################
def cpCalculate(originlist, ellipselist, debug=False):

    CriticalPointPairs = []
    for ellipse, origin in zip(ellipselist[1:], originlist[1:]):

        ellipseInverse = linalg.pinv(np.matrix(ellipse))

        # Calculate root of ellinver[0][0]. denominator -> a - c1 = sqrt(A^{-1}[0][0])
        mult = ellipseInverse.item((0, 0))
        denominator = np.sqrt(mult)

        # Calculate ellipseInverse * e1. e1 = [1, 0, 0, ..... 0]
        ellipseInvRow = np.array(ellipseInverse[0][:])

        # The two critical values come using the +/- for the sqrt value
        cpPairForEllipse = []
        for PosnegValue in [-1, 1]:
            cpVector = np.array(origin) + PosnegValue * (1 / denominator) * ellipseInvRow
            cpPairForEllipse.append(np.array(cpVector)[0].tolist())
        CriticalPointPairs.append(cpPairForEllipse)

    return CriticalPointPairs


def axisRange(Criticalpt, debug=False):
    # This function returns the range of the values of the travaxis and the other axis separately for easier calculation
    traversalAxisRange = []
    remainingAxisRange = []
    for Cpoint in Criticalpt:
        if len(Cpoint) != 2:
            raise ValueError("No Two Critical Points Found.", Cpoint, Criticalpt)

        # The smaller value, in respect to the X axis is first in the list
        traversalAxisRange.append([Cpoint[0][0], Cpoint[1][0]])
        remainingAxisRange.append([Cpoint[0][(1):], Cpoint[1][(1):]])

    if debug: print "Travaxis and otheraxis = ", traversalAxisRange, remainingAxisRange
    return traversalAxisRange, remainingAxisRange


####################
# Considerlist - holds the list of the ellipses that need to be considered for the intersect function
# The returned list is of dimension n. If the vaule is 1, the ellipse is to be considered, else not
# I think, the < and < are fine. For equality, there will be a critical slice at that slice
####################
def ellipseSliceintersect(CriticalX, value):
    considerlist = []
    for index, xcrit in enumerate(CriticalX):
        if min(xcrit) < value < max(xcrit):  # Or should it be <= and <= in both the inequalities??
            considerlist.append(1)
        else:
            considerlist.append(0)
    return considerlist


    # def ReduceSingleEllipsoid(arrayYZ, CPslicevector, travaxis, ellipse, origin, debug = False):
    # 	center2nd = np.array( arrayYZ )
    # 	trans = center2nd.T
    # 	numerator = 1 - np.dot(trans, center2nd)

    # 	if debug: print "d and d.T == ", center2nd, trans
    # 	if debug: print "numerator == ", numerator, "\n"

    # 	topcorner = ellipse[0][0]
    # 	planeState = CPslicevector[travaxis]
    # 	highercenter = origin[0]
    # 	denominator = 1 - topcorner * (planeState - highercenter)  * (planeState - highercenter)	# denom = 1 - aii (a - ci ) ^ 2

    # 	if debug: print "topcorner, planeState, highercenter == ", topcorner, planeState, highercenter
    # 	if debug: print "denominator == ", denominator, "\n"

    # 	firstcenter = np.array( origin[1:] )
    # 	trans = firstcenter.T
    # 	oldmatrix = np.array( [x[1:] for x in ellipse[1:]  ])
    # 	if debug: print "Firstcenter, it's trans, oldmatrix == ", firstcenter, trans, oldmatrix

    # 	value = np.dot( np.dot( trans, oldmatrix ), firstcenter )
    # 	denominator = denominator - value

    # 	if debug: print value
    # 	if debug: print "denominator == ", denominator, "\n"

    # 	firstrowA = np.array( ellipse[0][1:] )
    # 	trans = firstrowA.T
    # 	denominator = denominator + 2 * (planeState - highercenter) * np.dot( trans, firstcenter )

    # 	if debug: print "denominator == ", denominator, "\n"

    # 	tempellise = deepcopy(ellipse)

    # 	newmatrix = [x[1:] for x in tempellise[1:]  ]
    # 	if debug: print "Original	 newmatrix is == ",newmatrix
    # 	for col in range( len(newmatrix)):
    # 		for row in range( len(newmatrix[0]) ):
    # 			newmatrix[row][col] = newmatrix[row][col] * numerator/denominator
    # 	center2nd = arrayYZ

    # 	if debug: print "Obtained newmatrix is == ",newmatrix
    # 	if debug: print center2nd

    # 	return newmatrix, center2nd
