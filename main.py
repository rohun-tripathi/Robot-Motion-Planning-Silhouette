# !usr/bin/env python
import processAndValidateInput
import plot as myplt
import roadmap as roadmap
import shared as SH
import CreateRoadContext
from djikstra import Digraph


def plot_point_as_well():
    START = "start"
    END = "end"

    def findEnd():
        pass

    SH.init()
    diGraph = Digraph()
    debug = False
    outputFile = open("file_output/outputFile.txt", "w")

    inputFile = "input/input3D.txt"

    originList, ellipseList = processAndValidateInput.process_file_input(inputFile)

    endCriticalPoint = [0, 5.3, 0]
    context = CreateRoadContext.RoadContext()
    context.setEllipseListReturnSelf(ellipseList).setParentVectorReturnSelf([0 for x in range(SH.dim)]). \
        setOriginListReturnSelf(originList)
    #    setOriginListReturnSelf(originList).setCriticalPointsReturnSelf([[endCriticalPoint[:], END]])

    criticalPoints = roadmap.createRoad(context, debug)

    for point in SH.adjmatrix:
        Digraph.addEdge(diGraph, point[0], point[1], point[2])

    start = 0
    end = 491
    dist, path = Digraph.min_path(diGraph, start, end)
    # print >> outputFile, dist, "\n\n\n\n", path

    for index, coor in enumerate(SH.adjcoordinates):
        print >> outputFile, index, coor

    if len(path) < 2:
        print "Length of the path is less than 2 nodes. This is an error situation. Path == ", path
        sys.exit(1)
    tree = []
    for index in range(1, len(path)):
        tree.append([path[index - 1], path[index]])

    myplt.plot3DEllipseAndGraph(SH.num, originList, ellipseList, tree, SH.adjcoordinates)


def test_point():
    SH.init()
    inputFile = "input/input3D.txt"
    originList, ellipseList = processAndValidateInput.process_file_input(inputFile)
    create_road_map_and_plot(ellipseList, originList)


def run_silhouette_method(input_lines, file_name):
    SH.init()
    originList, ellipseList = processAndValidateInput.process_input(input_lines)
    create_road_map_and_plot(ellipseList, originList, file_name)


def create_road_map_and_plot(ellipseList, originList, file_name="foo.png"):
    context = CreateRoadContext.RoadContext()
    context.setEllipseListReturnSelf(ellipseList).setParentVectorReturnSelf([0 for x in range(SH.dim)]). \
        setOriginListReturnSelf(originList)
    criticalPoints = roadmap.createRoad(context)
    myplt.plot3DEllipseAndGraph(SH.num, originList, ellipseList, SH.adjmatrix, SH.adjcoordinates, file_name)


if __name__ == "__main__":
    test_point()