# !usr/bin/env python

import processAndValidateInput  # Actually use this
import plot as myplt
import roadmap as roadmap
import shared as SH
import CreateRoadContext

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

try:
    testPlottingFor4D()
except Exception:
    raise