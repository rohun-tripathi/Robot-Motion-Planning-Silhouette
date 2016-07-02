
#!usr/bin/env python
import numpy as np
import numpy.linalg as linalg
import sys, time, random
from optparse import OptionParser

import shared as SH
import auxilary as aux
import graph as graph
import plot as myplt
import complement #Actually use this
import roadmap as roadmap

###This code assumes that the input format is correct
debug = False
SH.init()				#initialize the global variables

(options, args) = OptionParser().parse_args()

originlist, ellipselist = complement.processinput([], [], args [0])			#processinput(originlist, ellipselist):

if debug == True:
	print "originlist = ",  originlist
	print "ellipselist = ",  ellipselist
	print "num = ",  SH.num
	print "dim = ",  SH.dim
	print "adjmatrix = ",  SH.adjmatrix
	print "primA = ", SH.primA
roadmap.CreateRoad(0, [], [0 for x in range(SH.dim)] ,originlist, ellipselist, False)		#CreateRoad(traversalAxis, CriticalPoints, parentvector ,originlist, ellipselist )	



# tree = main(valcount, adj)
# fig = plt.figure(i)

# createbasis(num)

# plt.show()
# plt.close(fig)
# del fig
#plotthis3()


#START WITH THIS COMMENT
myplt.plotthis(SH.num,originlist,ellipselist,SH.adjmatrix,SH.adjcoordinates)