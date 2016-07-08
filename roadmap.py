import sys
from copy import deepcopy

import auxilary
import cpfunctions as CriticalPointFunctions
import rdfunctions as rd
import shared as SHARED

from CreateRoadContext import RoadContext

CONSTANT_TO_ACCOUNT_FOR_PRIMARY_ELLIPSE_SKIPPED_FROM_LIST = 1


def createRoad(context, debug=False):

    context = initialize(context)
    #stateList

    originList = RoadContext.getOriginList(context)
    ellipseMatrixList = RoadContext.getEllipseList(context)
    CriticalPointPairs = CriticalPointFunctions.cpCalculate(originList, ellipseMatrixList, debug)
    CriticalAlongTraversal, CriticalAlongOthers = CriticalPointFunctions.axisRange(CriticalPointPairs, debug)

    traversalAxis = RoadContext.getTraversalAxis(context)
    sliceVector = RoadContext.getSliceVector(context)
    sliceVectorList = RoadContext.getSliceVectorList(context)
    pastVector = [] if len(sliceVectorList) < 2 else sliceVectorList[len(sliceVectorList) - 2]
    returnvec = RoadContext.getReturnvec(context)

    startPoint = sliceVector[traversalAxis]
    for iteration in range(SHARED.iterate):
        presentVector = []
        nextSliceValue = getNextSlice(startPoint, iteration)

        if iteration == 0 or iteration == SHARED.iterate - 1:
            pastVector, returnvec = processForFirstOrLastSlice(iteration, nextSliceValue, pastVector, presentVector,
                                                               returnvec, sliceVector, traversalAxis, debug)
            continue

        booleanValuesForEllipsesToConsiderExceptPrimary = CriticalPointFunctions.ellipseSliceintersect(CriticalAlongTraversal, nextSliceValue)

        otherAxesValueForEllsToConsiderExceptPrimary = CriticalPointFunctions.ellipseUnderConsider\
            (booleanValuesForEllipsesToConsiderExceptPrimary, CriticalAlongOthers,
             CriticalAlongTraversal, nextSliceValue, traversalAxis, debug)

        sliceVector[traversalAxis] = nextSliceValue  # No need for the old slice information?
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

            if not len(otherAxesValueForEllsToConsiderExceptPrimary[index]) == SHARED.dim - (traversalAxis + 1):
                raise IndexError("Error! The \"rest of slice vector\" and vector from otherAxesValueForEllsToConsiderExceptPrimary should have same len.",
                                 otherAxesValueForEllsToConsiderExceptPrimary[index], sliceVector, SHARED.dim - (traversalAxis + 1))

            vector = extractCorrectVectorFromOtherAxesValue(otherAxesValueForEllsToConsiderExceptPrimary, index, sliceVector, traversalAxis, debug)

            for axis2 in range(1, 2):  # (travaxis+1, SH.dim): todo
                findIntersectionPointsBetweenEllipseAndSlice(axis2, ellipse, origin, presentVector, traversalAxis,
                                                             vector)

        pastVector = auxilary.linkPresentAndPastVector(pastVector, presentVector, [], debug)

        # Now Begins Work of Critical points
        nextSliceValue = getNextSlice(startPoint, iteration + 1)  # nextSliceValue is a value of next iteration
        activeCPsForSlice = CriticalPointFunctions.retrieveActiveCPsForSlice(CriticalAlongTraversal, CriticalAlongOthers,
                                                                  sliceVector[traversalAxis],
                                                                  nextSliceValue)
        # To find which critical points from higher level have to be attached at this stage
        criticalPoints = RoadContext.getCriticalPoints(context)
        CritAtThisSlice = CriticalPointFunctions.RecurCheck(criticalPoints, sliceVector[traversalAxis], nextSliceValue,
                                                            False)

        if len(activeCPsForSlice) + len(CritAtThisSlice) > 1:
            print "System Limitation. We are just working with one critical point between two slices presently. A solution could be to increase the iterate variable, so slices come closer"
            print "It could be be the critical points from higher dimension and this one are crowding together."
            sys.exit(0)
        if len(activeCPsForSlice) + len(CritAtThisSlice) == 0: continue

        if (SHARED.dim - traversalAxis) == 2:  # Reached the 2D case - Base Case
            CPvector = sliceVector[:]  # CPvector will be manipulated to represent the critical val

            # print "CritAtThisSlice at sliceVector == ", CritAtThisSlice, sliceVector, activeCP
            try:
                point = CritAtThisSlice[0][0]  # System Limitation, only account for one critical point
            except IndexError:
                point = activeCPsForSlice[0][0]
            for index, term in enumerate(point):
                CPvector[traversalAxis + index] = term  # CPvector will be added to the Graph after this
            VectorNum = rd.addToVertices(CPvector)
            try:
                startendflag = CritAtThisSlice[0][1]
            except IndexError:
                startendflag = activeCPsForSlice[0][1]
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
            nextSliceValue = activeCPsForSlice[0][0][0]  # The travaxis value for the recursion slice
            booleanValuesForEllipsesToConsiderExceptPrimary = CriticalPointFunctions.ellipseSliceintersect(CriticalAlongTraversal, nextSliceValue)
            otherAxesValueForEllsToConsiderExceptPrimary = CriticalPointFunctions.ellipseUnderConsider(booleanValuesForEllipsesToConsiderExceptPrimary, CriticalAlongOthers,
                                                                     CriticalAlongTraversal, nextSliceValue,
                                                                     traversalAxis, False)
            sliceVector[traversalAxis] = nextSliceValue

            # Get the Ellipselist and the originlist for the lower dimensions
            RecursionEll, RecursionOri = CriticalPointFunctions.ReduceEllipsoids(booleanValuesForEllipsesToConsiderExceptPrimary, otherAxesValueForEllsToConsiderExceptPrimary,
                                                                                 sliceVector,
                                                                                 traversalAxis,
                                                                                 deepcopy(ellipseMatrixList),
                                                                                 deepcopy(originList), debug)

            # obtainedvec is the otherside of returnvec
            cpValuesAlongRemainingDim = CriticalPointFunctions.retrieveCPValuesAlongLowerDim(activeCPsForSlice)
            recursionContext = RoadContext()
            recursionContext.setCriticalPointsReturnSelf(cpValuesAlongRemainingDim).setTraversalAxisReturnSelf(traversalAxis + 1).\
                setEllipseListReturnSelf(RecursionEll).setOriginListReturnSelf(RecursionOri).setParentVectorReturnSelf(sliceVector)
            obtainedvec = createRoad(recursionContext, debug)

            presentVector = obtainedvec[0][:]
            presentVector.extend(obtainedvec[2][:])  # For Doubts regarding the obtainedvec/returnvec refer to README.md
            auxilary.endRecur_link(pastVector, presentVector, False)

            pastVector = obtainedvec[0][:]
            pastVector.extend(obtainedvec[1][:])  # Preparing the real pastVector for the next slice.

    for edgeindex, edge in enumerate(SHARED.adjmatrix):
        print >> SHARED.errorOut, edgeindex, "edge == ", edge
    for coorindex, coor in enumerate(SHARED.adjcoordinates):
        print >> SHARED.errorOut, coorindex, "coor == ", coor
    return returnvec


# ------------------------------------ Private Methods-----------------------------------#

# The statelist gives the status of whether a particular slice is inside or outside of the ellipse in question
# State dictates where a slice
# interacts with the nth ellipse
# Called inout in earlier versions

def processForFirstOrLastSlice(iteration, nextSlice, pastvector, presentVector, returnvec, sliceVector,
                               traversalAxis, debug):
    sliceVector[traversalAxis] = nextSlice  # No need for the old information
    returnvec, presentVector = LastOrFirstSlice(returnvec, sliceVector)
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


def initialize(context):
    traversalAxis = RoadContext.getTraversalAxis(context)
    parentVector = RoadContext.getParentVector(context)
    originList = RoadContext.getOriginList(context)
    ellipseMatrixList = RoadContext.getEllipseList(context)

    # Solve for startPoint, if solvable
    sol, valid = auxilary.getIntersectionPointAlongtravAxis \
        (ellipseMatrixList[0], originList[0], parentVector[traversalAxis:], traversalAxis, False)

    if valid != 1: raise ValueError("Error in initialize function in calculation of start point of slicing  "
                         "This is a serious error and never should have come up. Exiting")

    sliceVector = parentVector[:]  # In the first call, parent vector is [0 for x in range(SH.dim)]
    sliceVector[traversalAxis] = min(sol[0], sol[1]);  # Start point of recursion

    RoadContext.setSliceEllipseStateListReturnSelf(context, [0 for x in range(0, len(ellipseMatrixList))]) \
        .setSliceVectorReturnSelf(sliceVector).setSliceVectorListReturnSelf([]).\
        setReturnvecListReturnSelf([[], [], []])
    return context


def getNextSlice(startpt, iteration):
    # Distance moved on along the traversal axis in this iteration
    travaxisDist = startpt + ((startpt) * (-2) * (iteration)) / (SHARED.iterate - 1)
    return travaxisDist


def LastOrFirstSlice(returnvec, slicevector):
    VectorNum = rd.addToVertices(slicevector)
    returnvec[0].append(VectorNum)
    presentvector = [VectorNum]
    return returnvec, presentvector
