import sys
from copy import deepcopy

import auxilary
import cpfunctions as CriticalPointFunctions
import rdfunctions as rd
import shared as SHARED

CONSTANT_TO_ACCOUNT_FOR_PRIMARY_ELLIPSE_SKIPPED_FROM_LIST = 1


def createRoad(traversalAxis, criticalPoints, parentVector, originList, ellipseMatrixList, debug=False):
    stateList, firstLink, sliceVector, returnvec, pastVector, presentVector = initialize \
        (traversalAxis, parentVector, originList, ellipseMatrixList, debug)
    errorOut = open("debuginfo.txt", "w")

    CriticalPointPairs = CriticalPointFunctions.cpCalculate(traversalAxis, originList, ellipseMatrixList, debug)
    CriticalAlongTraversal, CriticalAlongOthers = CriticalPointFunctions.axisRange(CriticalPointPairs, debug)

    # Begin the iteration
    startPoint = sliceVector[traversalAxis]  # major todo
    for iteration in range(SHARED.iterate):
        presentVector = []
        nextSlice = getNextSlice(startPoint, iteration)  # nextSlice is a value

        if iteration == 0 or iteration == SHARED.iterate - 1:
            pastVector, returnvec = processForFirstOrLastSlice(iteration, nextSlice, pastVector, presentVector,
                                                               returnvec, sliceVector, traversalAxis, debug)
            continue

        booleanValuesForEllipsesToConsiderExceptPrimary = CriticalPointFunctions.ellipseSliceintersect(CriticalAlongTraversal, nextSlice)

        # The returned list is of dimension n-1.
        otherAxesValueForEllsInConsider = CriticalPointFunctions.ellipseUnderConsider\
            (booleanValuesForEllipsesToConsiderExceptPrimary, CriticalAlongOthers,
             CriticalAlongTraversal, nextSlice, traversalAxis, debug)

        sliceVector[traversalAxis] = nextSlice  # No need for the old slice information?
        # Technically if there were any CPs they should have been dealt with by this stage - todo - to be implemented
        #############################################
        # The Intersect and creation of presentVector
        #############################################

        # first for the first ellipse
        for axis2 in range(1, 2):  # (travaxis+1, SH.dim): todo - scale to n
            findIntersectionPointsBetweenEllipseAndSlice(axis2, ellipseMatrixList[0], originList[0], presentVector,
                                                         traversalAxis, sliceVector[:])

        for index, ellipseFlag in enumerate(booleanValuesForEllipsesToConsiderExceptPrimary):
            if ellipseFlag == 0: continue

            # We have to find the correct ellipse and remember to account for the first ellipse being primA
            if debug: print "Working for obstacle, index no. ", index
            ellipse = ellipseMatrixList[index + CONSTANT_TO_ACCOUNT_FOR_PRIMARY_ELLIPSE_SKIPPED_FROM_LIST]
            origin = originList[index + CONSTANT_TO_ACCOUNT_FOR_PRIMARY_ELLIPSE_SKIPPED_FROM_LIST]

            if not len(otherAxesValueForEllsInConsider[index]) == SHARED.dim - (traversalAxis + 1):
                raise IndexError("Error! The \"rest of slice vector\" and vector from otherAxesValueForEllsInConsider should have same len.",
                                 otherAxesValueForEllsInConsider[index], sliceVector, SHARED.dim - (traversalAxis + 1))

            vector = extractCorrectVectorFromOtherAxesValue(otherAxesValueForEllsInConsider, index, sliceVector, traversalAxis, debug)

            for axis2 in range(1, 2):  # (travaxis+1, SH.dim): todo
                findIntersectionPointsBetweenEllipseAndSlice(axis2, ellipse, origin, presentVector, traversalAxis,
                                                             vector)

        pastVector = auxilary.linkPresentAndPastVector(pastVector, presentVector, [], debug)

        # Now Begins Work of Critical points
        nextSlice = getNextSlice(startPoint, iteration + 1)  # nextSlice is a value of next iteration
        # ActiveCP and restCP have the value for any critcal points found at this stage. The active is along the travaxis and restCP is for recursion
        activeCP, restCP = CriticalPointFunctions.RecursionPoints(CriticalAlongTraversal, CriticalAlongOthers,
                                                                  sliceVector[traversalAxis],
                                                                  nextSlice, False)
        # To find which critical points from higher level have to be attached at this stage
        CritAtThisSlice = CriticalPointFunctions.RecurCheck(criticalPoints, sliceVector[traversalAxis], nextSlice,
                                                            False)

        if len(activeCP) + len(CritAtThisSlice) > 1:
            print "System Limitation. We are just working with one critical point between two slices presently. A solution could be to increase the iterate variable, so slices come closer"
            print "It could be be the critical points from higher dimension and this one are crowding together."
            sys.exit(0)
        if len(activeCP) + len(CritAtThisSlice) == 0: continue

        # We work on the CritAtThisSlice one first
        # if len(CritAtThisSlice) > 0:						#If count is greater than 0, first recursion worked
        # 	print activeCP, restCP, CritAtThisSlice, sliceVector
        # 	raw_input("reached here")

        if (SHARED.dim - traversalAxis) == 2:  # Reached the 2D case - Base Case
            CPvector = sliceVector[:]  # CPvector will be manipulated to represent the critical val

            # print "CritAtThisSlice at sliceVector == ", CritAtThisSlice, sliceVector, activeCP
            try:
                point = CritAtThisSlice[0][0]  # System Limitation, only account for one critical point
            except IndexError:
                point = [activeCP[0][0]] + restCP[0][0]
            for index, term in enumerate(point):
                CPvector[traversalAxis + index] = term  # CPvector will be added to the Graph after this
            VectorNum = rd.addToVertices(CPvector)
            try:
                startendflag = CritAtThisSlice[0][1]
            except IndexError:
                startendflag = restCP[0][1]
            if startendflag == "start":
                returnvec[1].append(VectorNum)
            elif startendflag == "end":
                returnvec[2].append(VectorNum)
            else:
                print "There is an error in the travaxis == ", traversalAxis, "The Criticalpt here does not have appended text (start/end)"
                sys.exit(0)

            auxilary.linkPresentAndPastVector(pastVector, [VectorNum], [], False)  # Add this VectorNum to the network
            pastVector.append(VectorNum)  # Added to pastvectors as will be needed in next slice
            continue  # Done for this iteration


        else:  # Not 2D. We have to prepare to call the lower dimension
            nextSlice = activeCP[0][0]  # The travaxis value for the recursion slice
            booleanValuesForEllipsesToConsiderExceptPrimary = CriticalPointFunctions.ellipseSliceintersect(CriticalAlongTraversal, nextSlice)
            otherAxesValueForEllsInConsider = CriticalPointFunctions.ellipseUnderConsider(booleanValuesForEllipsesToConsiderExceptPrimary, CriticalAlongOthers,
                                                                     CriticalAlongTraversal, nextSlice,
                                                                     traversalAxis, False)
            sliceVector[traversalAxis] = nextSlice

            # Get the Ellipselist and the originlist for the lower dimensions
            RecursionEll, RecursionOri = CriticalPointFunctions.ReduceEllipsoids(booleanValuesForEllipsesToConsiderExceptPrimary, otherAxesValueForEllsInConsider,
                                                                                 sliceVector,
                                                                                 traversalAxis,
                                                                                 deepcopy(ellipseMatrixList),
                                                                                 deepcopy(originList), debug)

            # obtainedvec is the otherside of returnvec
            obtainedvec = createRoad(traversalAxis + 1, restCP, sliceVector, RecursionOri, RecursionEll, debug=False)

            presentVector = obtainedvec[0][:]
            presentVector.extend(obtainedvec[2][:])  # For Doubts regarding the obtainedvec/returnvec refer to README.md
            auxilary.endRecur_link(pastVector, presentVector, False)

            pastVector = obtainedvec[0][:]
            pastVector.extend(obtainedvec[1][:])  # Preparing the real pastVector for the next slice.

    for edgeindex, edge in enumerate(SHARED.adjmatrix):
        print >> errorOut, edgeindex, "edge == ", edge
    for coorindex, coor in enumerate(SHARED.adjcoordinates):
        print >> errorOut, coorindex, "coor == ", coor
    return returnvec


# ------------------------------------ Private Methods-----------------------------------#

# The statelist gives the status of whether a particular slice is inside or outside of the ellipse in question
# State dictates where a slice
# interacts with the nth ellipse
# Called inout in earlier versions

def processForFirstOrLastSlice(iteration, nextSlice, pastvector, presentVector, returnvec, sliceVector,
                               traversalAxis, debug):
    if debug: print "in Slice " + str(iteration) + " sliceVector == ", sliceVector
    sliceVector[traversalAxis] = nextSlice  # No need for the old information
    if debug: print "in Slice 0, sliceVector == ", sliceVector
    returnvec, presentVector = LastOrFirstSlice(returnvec, presentVector, sliceVector)
    if iteration == SHARED.iterate - 1:
        pastvector = auxilary.complete_link(pastvector, presentVector, False)
    else:
        pastvector = presentVector[:]
    return pastvector, returnvec


# we have to extract the correct vector using the slice vector and the list in consider YZ
def extractCorrectVectorFromOtherAxesValue(considerYZ, index, sliceVector, traversalAxis, debug):
    vector = sliceVector[:]
    for i in range(SHARED.dim - (traversalAxis + 1)):
        vector[(traversalAxis + 1) + i] = considerYZ[index][i]
    if debug: print "For the sliceVector ", sliceVector, "gotten the following vector to run intersect along, ", vector
    return vector


def findIntersectionPointsBetweenEllipseAndSlice(axis2, ellipse, origin, presentVector, traversalAxis, vector):
    solution, valid = auxilary.getIntersectionPointAlongtravAxis(ellipse, origin, vector[(traversalAxis):],
                                                                 axis2, False)
    if valid == 1:
        if solution[0] == solution[1]:
            print "Check2, Error! This is not supposed to happen, both the solutions are equal."
            sys.exit(0)
        for sol in solution:
            vector[traversalAxis + 1] = sol
            VectorNum = rd.addToVertices(vector)
            presentVector.append(VectorNum)
            # returnvec[0].append( VectorNum )


def initialize(traversalAxis, parentVector, originList, ellipseMatrixList, debug=False):
    sliceEllipseStateList = [0 for x in range(0, len(ellipseMatrixList))]
    firstLinkInSlice = 0

    # Calculate startPoint
    sol, valid = auxilary.getIntersectionPointAlongtravAxis \
        (ellipseMatrixList[0], originList[0], parentVector[traversalAxis:], traversalAxis, False)

    if valid != 1:
        raise ValueError("Error in initialize function in calculation of start point of slicing  "
                         "This is a serious error and never should have come up. Exiting")
    elif debug:
        print "In initialize function in calculation of start point of slicing, Solutions == ", sol

    sliceVector = parentVector[:]  # In the first call, parent vector is [0 for x in range(SH.dim)]
    sliceVector[traversalAxis] = min(sol[0], sol[1]);  # Start point of recursion

    if debug:
        print "In intialize, the debug is on so the info is as follows"
        print "firstlink == ", firstLinkInSlice
        print "slicevector == ", sliceVector
        print "startpoint == ", sliceVector[traversalAxis]

    returnvec = [[], [], []]
    pastvector = []
    presentvector = []
    return sliceEllipseStateList, firstLinkInSlice, sliceVector, returnvec, pastvector, presentvector


def getNextSlice(startpt, iteration):
    # Distance moved on along the traversal axis in this iteration
    travaxisDist = startpt + ((startpt) * (-2) * (iteration)) / (SHARED.iterate - 1)
    return travaxisDist


def LastOrFirstSlice(returnvec, presentvector, slicevector):
    VectorNum = rd.addToVertices(slicevector)
    returnvec[0].append(VectorNum)
    presentvector = [VectorNum]
    return returnvec, presentvector
