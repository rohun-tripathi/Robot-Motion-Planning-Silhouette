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

    inputFile = "input/input4D.txt"

    originList, ellipseList = processAndValidateInput.process_file_input(inputFile)
    context = CreateRoadContext.RoadContext()
    context.setEllipseListReturnSelf(ellipseList).setParentVectorReturnSelf([0 for x in range(SH.dim)]). \
        setOriginListReturnSelf(originList)

    roadmap.createRoad(context, debug)

    myplt.projectHigherToLowerAndDisplay(originList, ellipseList)

def test3DOneAngularObstacle():
    # !usr/bin/env python

    import processAndValidateInput  # Actually use this
    import plot as myplt
    import roadmap as roadmap
    import shared as SH
    import CreateRoadContext

    # initialize the global variables
    SH.init()
    debug = False

    inputFile = "input/input3DOneAngular.txt"
    # inputFile = "input3DThreeEllipsoidsOneAtAnAngle.txt"
    # inputFile = "input3D.txt"


    originList, ellipseList = processAndValidateInput.process_file_input(inputFile)

    context = CreateRoadContext.RoadContext()
    context.setEllipseListReturnSelf(ellipseList).setParentVectorReturnSelf([0 for x in range(SH.dim)]). \
        setOriginListReturnSelf(originList)

    roadmap.createRoad(context, debug)

    # tree = main(valcount, adj)
    # fig = plt.figure(i)


    myplt.plot3DEllipseAndGraph(SH.num, originList, ellipseList, SH.adjmatrix, SH.adjcoordinates)


def testValidityFunction():
    # initialize the global variables
    SH.init()
    debug = False

    inputFile = "input/input4D.txt"

    originList, ellipseList = processAndValidateInput.process_file_input(inputFile)
    basicVector = [1,2,3,4]
    for i in range (50):
        basicVector[1] += 1
        value = rdfunctions.checkValidityOfPoint(basicVector, originList, ellipseList, 0)
        print basicVector, value
    return

try:
    test3DOneAngularObstacle()
except Exception:
    raise