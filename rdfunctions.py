import shared as SH
import sys
import numpy as np

from CreateRoadContext import RoadContext


def addToRoadmap(vector1, vector2, Distance,
                 String=''):  # String might have details like the function that called the addition to the adjacency matrix
    SH.adjmatrix.append([vector1, vector2, Distance, String])

# The value in adjcoordinates is added only here.
# This should be done ONLY after checking the the value is a valid one.


def addToVerticesUsingContext(vector, context):
    originList = RoadContext.getOriginList(context)
    ellipseMatrixList = RoadContext.getEllipseList(context)
    travAxis = RoadContext.getTraversalAxis(context)
    return addToVertices(vector, originList, ellipseMatrixList, travAxis)


def addToVertices(vector, originList, ellipseMatrixList, travAxis):
    if not checkValidityOfPoint(vector[travAxis:], originList, ellipseMatrixList, travAxis):
        return -1
    SH.valcount += 1
    SH.adjcoordinates.append(vector[:])
    return SH.valcount


def checkValidityOfPoint(pointVector, originList, ellMatrixList, travAxis):
    numEllipsoids = len(originList)
    dimension = SH.dim - travAxis

    for i in range(numEllipsoids):
        shiftedPoint = []
        temp = []
        for j in range(0, dimension):
            shiftedPoint.append(pointVector[j] - originList[i][j])
        temp.append(shiftedPoint)

        pointTranspose = np.array(temp)
        point = pointTranspose.T

        multiply = ellMatrixList[i]
        multiply = np.dot(pointTranspose, multiply)
        multiply = np.dot(multiply, point)

        if pointLiesOutsidePrimary(i, multiply[0][0]):
            return False
        if pointLiesInsideObstacle(i, multiply[0][0]):
            return False

    return True

# ------------------------------------ Private Methods-----------------------------------#


def pointLiesInsideObstacle(index, value):
    return index != 0 and abs(value - 1.0) > 1e-10 and value < 1.0

def pointLiesOutsidePrimary(index, value):
    return index == 0 and abs(value - 1.0) > 1e-10 and value > 1.0
