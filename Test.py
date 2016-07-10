# !usr/bin/env python

import processAndValidateInput  # Actually use this
import plot as myplt
import roadmap as roadmap
import shared as SH
import CreateRoadContext
import rdfunctions


def testPlottingFor4D():
    # initialize the global variables
    SH.init()
    debug = False

    inputFile = "input4D.txt"

    originList, ellipseList = processAndValidateInput.processInput(inputFile)
    context = CreateRoadContext.RoadContext()
    context.setEllipseListReturnSelf(ellipseList).setParentVectorReturnSelf([0 for x in range(SH.dim)]). \
        setOriginListReturnSelf(originList)

    roadmap.createRoad(context, debug)

    myplt.project4DTo3DAndDisplay(originList, ellipseList)


def testValidityFunction():
    # initialize the global variables
    SH.init()
    debug = False

    inputFile = "input4D.txt"

    originList, ellipseList = processAndValidateInput.processInput(inputFile)
    basicVector = [1,2,3,4]
    for i in range (50):
        basicVector[1] += 1
        value = rdfunctions.checkValidityOfPoint(basicVector, originList, ellipseList, 0)
        print basicVector, value
    return

try:
    testValidityFunction()
except Exception:
    raise