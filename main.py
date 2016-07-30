# !usr/bin/env python

import processAndValidateInput
import plot as myplt
import roadmap as roadmap
import shared as SH
import CreateRoadContext
from djikstra import Digraph

SH.init()
diGraph = Digraph()
debug = False
outputFile = open("outputFile.txt", "w")

inputFile = "input3DOneAngular.txt"
# inputFile = "input3DThreeEllipsoidsOneAtAnAngle.txt"
# inputFile = "input3D.txt"


originList, ellipseList = processAndValidateInput.processInput(inputFile)

context = CreateRoadContext.RoadContext()
context.setEllipseListReturnSelf(ellipseList).setParentVectorReturnSelf([0 for x in range(SH.dim)]).\
    setOriginListReturnSelf(originList).setCriticalPointsReturnSelf([[[5,0,0], "start"]])

criticalPoints = roadmap.createRoad(context, debug)

for point in SH.adjmatrix:
    Digraph.addEdge(diGraph, point[0], point[1], point[2]);

dist, path = Digraph.min_path(diGraph, 0, len(SH.adjcoordinates) - 1);
print >> outputFile, dist, "\n\n\n\n", path

if len(path) < 2 :
    print "Length of the path is less than 2 nodes. This is an error situation. Path == ", path
    sys.exit(1)
tree = []
for index in range(1, len(path)):
    tree.append([path[index-1], path[index]])

myplt.plot3DEllipseAndGraph(SH.num, originList, ellipseList, tree, SH.adjcoordinates)
# myplt.plot3DEllipseAndGraph(SH.num, originList, ellipseList, SH.adjmatrix, SH.adjcoordinates)
# myplt.projectHigherToLowerAndDisplay(originList, ellipseList)
