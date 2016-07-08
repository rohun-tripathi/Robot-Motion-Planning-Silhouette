# !usr/bin/env python

import processAndValidateInput  # Actually use this
import plot as myplt
import roadmap as roadmap
import shared as SH
import CreateRoadContext

# initialize the global variables
SH.init()
debug = False

inputFile = "input3D_2.txt"

originList, ellipseList = processAndValidateInput.processInput([], [], inputFile)

context = CreateRoadContext.RoadContext()
context.setEllipseListReturnSelf(ellipseList).setParentVectorReturnSelf([0 for x in range(SH.dim)]).\
    setOriginListReturnSelf(originList)

roadmap.createRoad(context, debug)

# tree = main(valcount, adj)
# fig = plt.figure(i)

# createbasis(num)

# plt.show()
# plt.close(fig)
# del fig
# plotthis3()


# START WITH THIS COMMENT
myplt.plotthis(SH.num, originList, ellipseList, SH.adjmatrix, SH.adjcoordinates)