import pprint as pprint
import sys
from copy import deepcopy

import numpy as np
import numpy.linalg as linalg

import shared as SH

import rdfunctions as rd

import CriticalPointsDetail

LARGE_VALUE = sys.maxsize


def extract_matrix(X, reduced_dim):
    result_matrix = []
    for row_iterator in range(reduced_dim):
        temp = []
        for col_iterator in range(reduced_dim):
            temp.append(X[row_iterator * reduced_dim + col_iterator])
        result_matrix.append(temp)

    return result_matrix


def reduce_ellipsoid(array_along_y_z, critical_slice_vector, traversal_axis, ellipsoid, origin):

    reduced_dim = len(ellipsoid[0]) - 1
    if len(array_along_y_z) != reduced_dim:
        raise Exception("Error! len(reduced_center) != reduced_dim\nExiting")

    top_corner = ellipsoid[0][0]
    critical_point = critical_slice_vector[traversal_axis]
    top_most_center = origin[0]

    initial_matrix = np.array([x[1:] for x in ellipsoid[1:]])

    # denominator = 1 - aii (a - ci ) ^ 2
    denominator = 1 - top_corner * (critical_point - top_most_center) * (critical_point - top_most_center)

    initial_center = np.array(origin[1:])
    value = np.dot(np.dot(initial_center.T, initial_matrix), initial_center)
    denominator -= value

    first_row = np.array(ellipsoid[0][1:])
    denominator += 2 * (critical_point - top_most_center) * np.dot(first_row.T, initial_center)

    reduced_center = np.array(array_along_y_z)
    # transpose_center = reduced_center.T
    # This was the mistake
    # numerator = 1 - np.dot(transpose_center, reduced_center)

    stretched_center = []

    for row_iterator in range(reduced_dim):
        for col_iterator in range(reduced_dim):
            stretched_center.append(reduced_center[row_iterator] * reduced_center[col_iterator])

    # begins the preparation for equation. AX = Y
    # X in the equations is the terms of the "B" matrix (for the ellipse in the lower dim)  len -> (n-1) * (n-1)
    # A is the coefficients of the terms in the "B" matrix
    # Y is the constant term of each equation

    A = []
    Y = []

    for row_iterator in range(reduced_dim):
        for col_iterator in range(reduced_dim):
            initial_matrix_term = ellipsoid[1 + row_iterator][1 + col_iterator]
            Y.append(initial_matrix_term)  # Accounting for the lowering in dimensions
            result = np.multiply(initial_matrix_term, stretched_center)
            result[row_iterator * reduced_dim + col_iterator] += denominator
            A.append(result.tolist())

    X = np.linalg.solve(A, Y)
    return extract_matrix(X, reduced_dim), reduced_center


def ReduceEllipsoids(considerlist, considerYZ, CPslicevector, travaxis, ellipselist, originlist, debug=False):
    recursion_ellipsoids = []
    recursion_origins = []

    # The first (outer) Ellipse
    arrayYZ = CPslicevector[(travaxis + 1):]
    ellipse = deepcopy(ellipselist[0])
    origin = originlist[0]

    reduced_matrix, center2nd = reduce_ellipsoid(arrayYZ, CPslicevector, travaxis, ellipse, origin)
    recursion_ellipsoids.append(reduced_matrix[:])
    recursion_origins.append(center2nd[:])

    for index, term in enumerate(considerlist):
        if term == 1:
            arrayYZ = considerYZ[index]

            ellipse = ellipselist[index + 1]  # Accounting for the outer ellipse
            origin = originlist[index + 1]

            reduced_matrix, center2nd = reduce_ellipsoid(arrayYZ, CPslicevector, travaxis, ellipse, origin)

            recursion_ellipsoids.append(reduced_matrix[:])
            recursion_origins.append(center2nd[:])
    return recursion_ellipsoids, recursion_origins


def retrieveCriticalPointsFromHigherLevel(critical, presentslice, nextslice, debug=False):
    CritAtThisSlice = []
    if debug: print "In RecurCheck and critical, presentslice, nextslice == ", critical, presentslice, nextslice
    for item in critical:
        # The < and <= essentially means that the obstacle cannot be sticking to the first point in the arena
        if presentslice < item[0][0] <= nextslice:
            CritAtThisSlice.append(item)
    return CritAtThisSlice


def retrieveActiveCPsForSlice(criticalPointPairList, presentSlice, nextSlice):
    activeCPsForSlice = []
    for cpPair in criticalPointPairList:
        firstCP = cpPair.getFirstCP()
        secondCP = cpPair.getSecondCP()

        if cpPair.getValidityOfFirst() and presentSlice < firstCP[0] < nextSlice:
            activeCPsForSlice.append([firstCP[:], "start"])
        elif cpPair.getValidityOfSecond() and presentSlice < secondCP[0] < nextSlice:
            activeCPsForSlice.append([secondCP[:], "end"])
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
def cpCalculate(originList, ellipseList, travAxis):

    CriticalPointPairs = []
    for ellipse, origin in zip(ellipseList[1:], originList[1:]):

        ellipseInverse = linalg.pinv(np.matrix(ellipse))

        # Calculate root of ellinver[0][0]. denominator -> a - c1 = sqrt(A^{-1}[0][0])
        multiply = ellipseInverse.item((0, 0))
        denominator = np.sqrt(multiply)

        # Calculate ellipseInverse * e1. e1 = [1, 0, 0, ..... 0]
        ellipseInvRow = np.array(ellipseInverse[0][:])

        # The two critical values come using the +/- for the sqrt value
        cpPair = populateCriticalPointValueAndValidity(denominator, ellipseInvRow, ellipseList, origin, originList,
                                                       travAxis)
        CriticalPointPairs.append(cpPair)

    return CriticalPointPairs


def populateCriticalPointValueAndValidity(denominator, ellipseInvRow, ellipseList, origin, originList, travAxis):
    cpPair = CriticalPointsDetail.CriticalPointPair()

    cpVector = np.array(origin) + -1 * (1 / denominator) * ellipseInvRow
    if rd.checkValidityOfPoint(cpVector[0], originList, ellipseList, travAxis):
        cpPair.setFirstCPReturnSelf(np.array(cpVector)[0].tolist()).setValidityOfFirst(True)
    else:
        cpPair.setFirstCPReturnSelf(np.array(cpVector)[0].tolist()).setValidityOfFirst(False)

    cpVector = np.array(origin) + 1 * (1 / denominator) * ellipseInvRow
    if rd.checkValidityOfPoint(cpVector[0], originList, ellipseList, travAxis):
        cpPair.setSecondCPReturnSelf(np.array(cpVector)[0].tolist()).setValidityOfSecond(True)
    else:
        cpPair.setSecondCPReturnSelf(np.array(cpVector)[0].tolist()).setValidityOfSecond(False)
    return cpPair


def axisRange(criticalPointPairList):
    # This function returns the range of the values of the travaxis and the other axis separately for easier calculation
    traversalAxisRange = []
    remainingAxisRange = []
    for cpPair in criticalPointPairList:

        firstCP = cpPair.getFirstCP()
        secondCP = cpPair.getSecondCP()
        traversalAxisRange.append([firstCP[0], secondCP[0]])
        remainingAxisRange.append([firstCP[1:], secondCP[1:]])

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
