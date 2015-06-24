import numpy as np
import numpy.linalg as linalg
import sys, time, random

import shared as SH
import auxilary as aux
import graph as graph
import plot as myplt
import complement as cmpl #Actually use this
import cpfunctions

def initRecursion():
	returnvec = [[],[],[]]
	pastvector = []
	presentvector = []
	return returnvec, pastvector, presentvector

def initialize (travaxis, parentvector, originlist, ellipselist, debug = False):
	statelist = [0 for x in range(0,len(ellipselist))]	#The statelist gives the status of whether a particular slice is inside or outside of the ellipse in question
											# State dictates where a slice
											# interacts with the nth ellipse
											# Called inout in earlier versions
	firstlink = 0							#I think this the value that is to be used to distinguish a change in intersection points due to critical points and due to first point in the slicing
	
	#Code for startpoint
	sol, valid = aux.intersect(ellipselist[0], originlist[0], parentvector[(travaxis):], travaxis, False)
	if valid != 1:
		print "Error in rdfunctions in initialize function in calculation of start point of slicing"
		print "This is a serious error and never should have come up. Exiting"
		sys.exit(0)
	if valid == 1 and debug == True:
		print "In rdfunctions in initialize function in calculation of start point of slicing"
		print sol
		time.sleep(2)

	startpoint = min(sol[0], sol[1]);		#Start point of recursion

	slicevector = parentvector[:]			#In the first call, parent vector is [0 for x in range(SH.dim)]
	slicevector[travaxis] = startpoint;		#Start point of recursion

	if debug == True:
		print "In rdfunctions, in intialize, the debug is on so the info is as follows"
		print "statelist == ", statelist
		print "firstlink == ", firstlink
		print "slicevector == ", slicevector
		print "startpoint == ", startpoint

	return statelist, firstlink, slicevector, startpoint


def LastorFirstSlice(returnvec, presentvector, slicevector):
	VectorNum = AddToVertices(slicevector)
	returnvec[0].append(VectorNum)
	presentvector = [VectorNum]
	return returnvec, presentvector

def AddToRoadmap(vector1, vector2, Distance, String = ''):	#String might have details like the function that called the addition to the adjacency matrix
	SH.adjmatrix.append( [vector1, vector2, Distance, String] )

def AddToVertices(vector1):
	SH.valcount += 1
	SH.adjcoordinates.append(vector1[:])
	return SH.valcount




def getNextSlice(startpt, iteration):
	travaxisDist = startpt + ((startpt) * (-2) * (iteration) )/(SH.iterate-1) 	#Distance moved on along the traversal axis in this iteration
	return travaxisDist