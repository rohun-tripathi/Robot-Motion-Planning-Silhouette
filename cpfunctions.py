import numpy as np
import numpy.linalg as linalg
import sys, time, random

import shared as SH
import auxilary as aux
import graph as graph
import plot as myplt
import complement as cmpl #Actually use this


def EllUnderConsider(considerlist , CriticalYZ, CriticalX, nextslice, travaxis, debug = False):
	if debug == True: print considerlist , CriticalYZ, CriticalX, nextslice, travaxis
	considerYZ = []
	for index, listitem in enumerate(considerlist):
		if listitem == 0:
			considerYZ.append([])
			continue
		else:
			templist = [0 for x in range(travaxis + 1, SH.dim) ]
			Xvalues = (nextslice - CriticalX[index][0]) / (CriticalX[index][1] - CriticalX[index][0]) 	#Xvalues = (x - x0) / (x1 - x0)

			for dimen in range(SH.dim - (travaxis + 1) ):
				if debug == True: print CriticalYZ[index], CriticalYZ[index][1]
				templist[dimen] = CriticalYZ[index][1][dimen]  - CriticalYZ[index][0][dimen]	#Val = (y1 - y0)
				templist[dimen] *= Xvalues										#Val = (x - x0) * (y1 - y0) / (x1 - x0)	
				templist[dimen] += CriticalYZ[index][0][dimen]					#Val = (x - x0) * (y1 - y0) / (x1 - x0) + y0
			considerYZ.append(templist)
	return considerYZ

def CPcalculate(travaxis, originlist, ellipselist, debug = False):
	#I need to this for each ellipse
	Criticalpt = []
	for ellipse, origin in zip( ellipselist[1:], originlist[1:]):
		ellinver = linalg.pinv( np.matrix(ellipse))
		if debug == True: print "ellinver == " ,ellinver

		#Equation we are using is point = center + 1/(sqrt(ellinver[0][0])) * ellinver * E1
		#Here E1 is a one colum matrix with only first value = 1 and others 0

		#Calculate ellinver[][]
		mult = ellinver.item((0, 0)) # works only for [0][0]
		if debug == True: print "mult == " ,mult	
		criticalval = np.sqrt(mult)

		if debug == True: print "criticalval == " ,criticalval

		#Calculate ellinver * E1
		#The length depends n the traversal axis presently. if at the zeroth level then the complete thing. Else [travaxis:]
		ellinvrow = np.array(ellinver[0][:])
		if debug == True: print "Ellinovrow == " ,ellinvrow
		#Calculate p
		posneg = [-1,1]
		#The two cricat values come using the +/- for the sqrt value
		CPellise = []
		for PosnegValue in posneg:
			cpvector = PosnegValue * 1/criticalval * ellinvrow + np.array(origin)
			#The main equation ^ of this part of the code. Equation of "Criticality" :P
			if debug == True: print "cpvector == " ,cpvector
			CPellise.append(np.array(cpvector)[0].tolist())
		if debug == True: print "CPellise == " ,CPellise
		Criticalpt.append(CPellise)
		if debug == True: print "Criticalpt == ", Criticalpt
	return Criticalpt

def axisrange(Criticalpt, debug = False):
	#This fucntion returns the range of the values of the travaxis and the other axis separately for easier calculation
	travaxisrange = []
	otheraxisrange = []
	for Cpoint in Criticalpt:
		travaxisrange.append( [Cpoint[0][0], Cpoint[1][0]] )					#The smaller value, in respect to the X axis is first in the list
		otheraxisrange.append ( [Cpoint[0][ (1) :], Cpoint[1][ (1) :]] )
	if debug == True: print "Travaxis and otheraxis = ", travaxisrange, otheraxisrange
	return travaxisrange, otheraxisrange

#Considerlist - holds the list of the ellipses that need to be considered for the intersect function
#The returned list is of dimension n. If the vaule is 1, the ellipse is to be considered, else not
def EllSliceintersect(CriticalX, value):
	considerlist = []
	for index, xcrit in enumerate(CriticalX):
		if min(xcrit) < value < max(xcrit):   #Or should it be <= and <= in both the inequalities??
			considerlist.append(1)
		else:
			considerlist.append(0)
	return considerlist
