# !usr/bin/env python

import processAndValidateInput  # Actually use this
import plot as myplt
import roadmap as roadmap
import shared as SH


def printInitialValues():
    if debug:
        print "originlist = ", originList
        print "ellipselist = ", ellipseList
        print "num = ", SH.num
        print "dim = ", SH.dim
        print "adjmatrix = ", SH.adjmatrix
        print "primA = ", SH.primA
        print "adjcoordinates = ", SH.adjcoordinates


# This code assumes that the input format is correct
inputFile = "input3D.txt"

# initialize the global variables
SH.init()
debug = False

# processinput(originlist, ellipselist)
originList, ellipseList = processAndValidateInput.processInput([], [], inputFile)
printInitialValues()

# Start the process CreateRoad(traversalAxis, CriticalPoints, parentvector ,originlist, ellipselist )

roadmap.createRoad(0, [], [0 for x in range(SH.dim)], originList, ellipseList, False)

# tree = main(valcount, adj)
# fig = plt.figure(i)

# createbasis(num)

# plt.show()
# plt.close(fig)
# del fig
# plotthis3()


# START WITH THIS COMMENT
myplt.plotthis(SH.num, originList, ellipseList, SH.adjmatrix, SH.adjcoordinates)