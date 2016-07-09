import sys
from copy import deepcopy

import auxilary
import cpfunctions as CriticalPointFunctions
import rdfunctions as rd
import shared as SHARED

from RecursionResult import VectorsToLink
from CreateRoadContext import RoadContext

CONSTANT_TO_ACCOUNT_FOR_PRIMARY_ELLIPSE_SKIPPED_FROM_LIST = 1


#Consider the solution for having a two lists passed down in the context. One has the start set of the CP being passed
#Other has the set relevant to called recursion.
#reason : will remove this travaxis and travaxis + 1 bullshit
#the whole point is needed when adding it to the graph

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
    returnvec = VectorsToLink()

    startPoint = sliceVector[traversalAxis]
    for iteration in range(SHARED.iterate):
        presentVector = resetVector()
        nextSliceValue = getNextSlice(startPoint, iteration)

        if iteration == 0 or iteration == SHARED.iterate - 1:
            pastVector, returnvec = processForFirstOrLastSlice(iteration, nextSliceValue, pastVector,
                                                               returnvec, sliceVector, traversalAxis)
            continue

        booleanValuesForEllipsesToConsiderExceptPrimary = CriticalPointFunctions.ellipseSliceintersect(CriticalAlongTraversal, nextSliceValue)

        otherAxesValueForEllsToConsiderExceptPrimary = CriticalPointFunctions.ellipseUnderConsider\
            (booleanValuesForEllipsesToConsiderExceptPrimary, CriticalAlongOthers,
             CriticalAlongTraversal, nextSliceValue, traversalAxis, debug)

        sliceVector[traversalAxis] = nextSliceValue  # No need for the old slice information?

        #todo - to be implemented, or does it?
        # Technically if there were any CPs they should have been dealt with by this stage -
        #############################################

        checkIntersectionAlongOtherAxisForEllipsoidsForSlice(booleanValuesForEllipsesToConsiderExceptPrimary, debug,
                                                             ellipseMatrixList, originList,
                                                             otherAxesValueForEllsToConsiderExceptPrimary,
                                                             presentVector, sliceVector, traversalAxis)

        pastVector = auxilary.linkPresentAndPastVector(pastVector, presentVector, [], debug)

        # Now Begins Work of Critical points
        nextSliceValue = getNextSlice(startPoint, iteration + 1)  # nextSliceValue is a value of next iteration
        activeCPsForSlice = CriticalPointFunctions.retrieveActiveCPsForSlice(CriticalAlongTraversal, CriticalAlongOthers,
                                                                  sliceVector[traversalAxis],
                                                                  nextSliceValue)
        # To find which critical points from higher level have to be attached at this stage
        criticalPoints = RoadContext.getCriticalPoints(context)
        criticalPointsFromHigherLevel = CriticalPointFunctions.retrieveCriticalPointsFromHigherLevel(criticalPoints, sliceVector[traversalAxis], nextSliceValue,
                                                                                                     False)
        criticalPointsForSlice = []
        criticalPointsForSlice.extend(criticalPointsFromHigherLevel)
        criticalPointsForSlice.extend(activeCPsForSlice)

        if len(criticalPointsForSlice) == 0: continue
        if len(criticalPointsForSlice) > 1:
            print "System Limitation. We are just working with one critical point between two slices presently. A solution could be to increase the iterate variable, so slices come closer"
            print "It could be be the critical points from higher dimension and this one are crowding together."
            sys.exit(0)

        if (SHARED.dim - traversalAxis) == 2:  # 2D case - Base Case
            pastVector = connectCPInBase2DCase(criticalPointsForSlice, pastVector, returnvec, sliceVector, traversalAxis)
            continue


        else:  # Not 2D. We have to prepare to call the lower dimension
            nextSliceValue = criticalPointsForSlice[0][0][0]  # The travaxis value for the recursion slice
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

            cpValuesAlongRemainingDim = CriticalPointFunctions.retrieveCPValuesAlongLowerDim(criticalPointsForSlice)
            recursionContext = RoadContext()
            recursionContext.setCriticalPointsReturnSelf(cpValuesAlongRemainingDim).setTraversalAxisReturnSelf(traversalAxis + 1).\
                setEllipseListReturnSelf(RecursionEll).setOriginListReturnSelf(RecursionOri).setParentVectorReturnSelf(sliceVector)
            vectorFromLowerRecursion = createRoad(recursionContext, debug)

            presentVector = vectorFromLowerRecursion.getIntersectionPointsAlongOtherAxis()[:]
            presentVector.extend(vectorFromLowerRecursion.getCriticalPointsToLinkToPast()[:])
            auxilary.endRecur_link(pastVector, presentVector, False)

            pastVector = vectorFromLowerRecursion.getIntersectionPointsAlongOtherAxis()[:]
            pastVector.extend(vectorFromLowerRecursion.getCriticalPointsToLinkToFuture()[:])  # Preparing the real pastVector for the next slice.

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

def processForFirstOrLastSlice(iteration, nextSlice, pastvector, returnvec, sliceVector,
                               traversalAxis):
    sliceVector[traversalAxis] = nextSlice  # No need for the old information

    VectorNum = rd.addToVertices(sliceVector[:])
    returnvec.getIntersectionPointsAlongOtherAxis().append(VectorNum)  # done
    presentVector = [VectorNum]

    if iteration == SHARED.iterate - 1:
        pastvector = auxilary.complete_link(pastvector, presentVector, False)
    else:
        pastvector = presentVector[:]
    return pastvector, returnvec


def checkIntersectionAlongOtherAxisForEllipsoidsForSlice(booleanValuesForEllipsesToConsiderExceptPrimary, debug,
                                                         ellipseMatrixList, originList,
                                                         otherAxesValueForEllsToConsiderExceptPrimary, presentVector,
                                                         sliceVector, traversalAxis):
    # first for the first ellipse
    axisToCheckIntersectionAlong = 1
    findIntersectionPointsBetweenEllipseAndSlice(axisToCheckIntersectionAlong, ellipseMatrixList[0], originList[0],
                                                 presentVector, traversalAxis, sliceVector[:])

    for index, ellipseFlag in enumerate(booleanValuesForEllipsesToConsiderExceptPrimary):
        if ellipseFlag == 0: continue

        # We have to find the correct ellipse and remember to account for the first ellipse being primA
        ellipse = ellipseMatrixList[index + CONSTANT_TO_ACCOUNT_FOR_PRIMARY_ELLIPSE_SKIPPED_FROM_LIST]
        origin = originList[index + CONSTANT_TO_ACCOUNT_FOR_PRIMARY_ELLIPSE_SKIPPED_FROM_LIST]

        if not len(otherAxesValueForEllsToConsiderExceptPrimary[index]) == SHARED.dim - (traversalAxis + 1):
            raise IndexError(
                "Error! The \"rest of slice vector\" and vector from otherAxesValueForEllsToConsiderExceptPrimary should have same len.",
                otherAxesValueForEllsToConsiderExceptPrimary[index], sliceVector, SHARED.dim - (traversalAxis + 1))

        vector = extractCorrectVectorFromOtherAxesValue(otherAxesValueForEllsToConsiderExceptPrimary, index,
                                                        sliceVector, traversalAxis, debug)

        findIntersectionPointsBetweenEllipseAndSlice(axisToCheckIntersectionAlong, ellipse, origin, presentVector,
                                                     traversalAxis, vector)


def connectCPInBase2DCase(criticalPointsForSlice, pastVector, returnvec, sliceVector, traversalAxis):
    cpVector = sliceVector[:]
    point = criticalPointsForSlice[0][0]  # System Limitation, only account for one critical point

    for index, term in enumerate(point):
        cpVector[traversalAxis + index] = term
    VectorNum = rd.addToVertices(cpVector)

    startendflag = criticalPointsForSlice[0][1]
    if startendflag == "start":
        returnvec.getCriticalPointsToLinkToFuture().append(VectorNum)
    elif startendflag == "end":
        returnvec.getCriticalPointsToLinkToPast().append(VectorNum)

    presentVector = [VectorNum]
    auxilary.linkPresentAndPastVector(pastVector, presentVector, [])  # Add this VectorNum to the network
    pastVector.append(VectorNum)  # Added to pastvectors as will be needed in next slice
    return pastVector



# we have to extract the correct vector using the slice vector and the list in consider YZ
def extractCorrectVectorFromOtherAxesValue(considerYZ, index, sliceVector, traversalAxis, debug):
    vector = sliceVector[:]
    for i in range(SHARED.dim - (traversalAxis + 1)):
        vector[(traversalAxis + 1) + i] = considerYZ[index][i]
    return vector


def findIntersectionPointsBetweenEllipseAndSlice(axisToCheckIntersectionAlong, ellipse, origin, presentVector, traversalAxis, vector):
    solution, valid = auxilary.getIntersectionPointAlongtravAxis(ellipse, origin, vector[(traversalAxis):],
                                                                 axisToCheckIntersectionAlong, False)
    if valid == 1:
        if solution[0] == solution[1]:
            print "Check2, Error! This is not supposed to happen, both the solutions are equal."
            sys.exit(0)
        for sol in solution:
            vector[traversalAxis + 1] = sol
            VectorNum = rd.addToVertices(vector)
            presentVector.append(VectorNum)


def initialize(context):
    traversalAxis = RoadContext.getTraversalAxis(context)
    parentVector = RoadContext.getParentVector(context)
    originList = RoadContext.getOriginList(context)
    ellipseMatrixList = RoadContext.getEllipseList(context)

    # Solve for startPoint, if solvable
    sol, valid = auxilary.getIntersectionPointAlongtravAxis \
        (ellipseMatrixList[0], originList[0], parentVector[traversalAxis:], 0, False)

    if valid != 1: raise ValueError("Error in initialize function in calculation of start point of slicing  "
                         "This is a serious error and never should have come up. Exiting")

    sliceVector = parentVector[:]  # In the first call, parent vector is [0 for x in range(SH.dim)]
    sliceVector[traversalAxis] = min(sol[0], sol[1]);  # Start point of recursion

    RoadContext.setSliceEllipseStateListReturnSelf(context, [0 for x in range(0, len(ellipseMatrixList))]) \
        .setSliceVectorReturnSelf(sliceVector).setSliceVectorListReturnSelf([])
    return context


def resetVector():
    return []


def getNextSlice(startpt, iteration):
    # Distance moved on along the traversal axis in this iteration
    travaxisDist = startpt + ((startpt) * (-2) * (iteration)) / (SHARED.iterate - 1)
    return travaxisDist

