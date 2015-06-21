import numpy as np
import numpy.linalg as linalg
import sys, time, random

import shared as SH
import auxilary as aux
import graph as graph
import plot as myplt
import complement as cmpl #Actually use this
import cpfunctions as cp
import rdfunctions as rd


#The main method
def CreateRoad(travaxis, critical, parentvector ,originlist, ellipselist, debug = False):
	#Start of initialize	
	
	errorout = open("debuginfo.txt", "w")

	statelist, firstlink, slicevector, startpt = rd.initialize(travaxis, parentvector, originlist, ellipselist, False)
	returnvec, pastvector, presentvector = rd.initRecursion()

	if debug == True:
		print "Originlist = ",  originlist
		print "ellipselist = ",  ellipselist
		print "num = ",  SH.num
		print "dim = ",  SH.dim
		print "slicevector = ",  slicevector
		print "adjcoordinates = ",  SH.adjcoordinates

	#CPcalculate - Calculates the Critical Points, for all but the primA ellipse
	#The returned list - Criticalpts, is a list of length n-1, where n is the ellipses being considered for this slice
	#Create separate lists for the travaxis and other axes, both of length n-1
	Criticalpts = cp.CPcalculate(travaxis, originlist, ellipselist, False)
	CriticalX, 	CriticalYZ = cp.axisrange(Criticalpts, True	)

	#Begin the iteration
	for iteration in range(SH.iterate):
		presentvector = []
		nextslice = rd.getNextSlice(startpt, iteration)		#nextslice is a value
		
		if iteration == 0 or iteration == SH.iterate -1 :
			## ## All this goes to a different function

			if debug == True: print "in Slice "+ str(iteration) + " slicevector == ", slicevector
			slicevector[travaxis] = nextslice 			#No need for the old information
			if debug == True: print "in Slice 0, slicevector == ", slicevector
			returnvec, presentvector = rd.LastorFirstSlice(returnvec, presentvector, slicevector)
			if iteration == SH.iterate -1 :
				pastvector = aux.complete_link(pastvector, presentvector, False)
			else:
				pastvector = presentvector[:]

			continue

		# cplist = cp.recursionpoints(CriticalX, slicevector[travaxis], nextslice)# Holds the list of the Criticalpts that need to taken care of by recursion
		# I will code the cplists code later, for christsake
		# When I do, code for  tha 2D case first

		#Considerlist - holds the list of the ellipses that need to be considered for the intersect function
		#The returned list is of dimension n-1. If the vaule is 1, the ellipse is to be considered, else not
		considerlist = cp.EllSliceintersect(CriticalX, nextslice)

		slicevector[travaxis] = nextslice 			#No need for the old slice information?
				#Technically if there were any CPs they should have been dealt with by this stage
		
		#The returned list is of dimension n-1.
		considerYZ = cp.EllUnderConsider(considerlist , CriticalYZ, CriticalX,  nextslice, travaxis, False)

#############################################
		#The Intersect and creation of presentvector
#############################################

		#first for the first ellipse
		for axis2 in range(travaxis+1, travaxis + 2):# (travaxis+1, SH.dim):
			solution, valid = aux.intersect(ellipselist[0], originlist[0], slicevector, axis2)
			vector = slicevector[:]

			if valid ==1 :		#Later check also that it is not inside any of the other ellipsoids
				if solution[0] == solution[1] and firstlink!= 1:
					print "Error, this is not supposed to happen, both the solutions are equal. But letting it be for Now. iteration == ", iteration
					print solution, slicevector
					#sys.exit(0)
				for sol in solution:
					vector[axis2] = sol
					VectorNum = rd.AddToVertices(vector[:])
					presentvector.append( VectorNum )
					returnvec[0].append( VectorNum )
				#Some processing needs to be done
			

		for index, listitem in enumerate(considerlist):
			if listitem == 0 :
				continue
			#else: Now the rest of the code for this loop

			# We have to find the correct ellipse and remember to account for the first ellipse being primA 			
			print "Working for obstacle, index no. ", index
			ellipse = ellipselist[index + 1]
			origin = originlist[index + 1]

			# also we have to extract the correct vector using the slice vector and the list in consider YZ
			vector = slicevector[:]

			#The rest of slice vector and the vector from considerYZ shud have same len, check this
			if not len(considerYZ[index]) == SH.dim - (travaxis + 1):
				print "Check2, Error! This is not supposed to happen, both the matrix should have equal dimensions."
				sys.exit(0)
			
			for i in range(SH.dim - (travaxis + 1) ):
				vector[	(travaxis + 1) + i] = considerYZ[index][i]
			print "For the slicevector ", slicevector, "gotten the following vector to run intersect along, ", vector
			
			#Run Intersect
			for axis2 in range(travaxis+1, travaxis + 2):# (travaxis+1, SH.dim):
				solution, valid = aux.intersect(ellipse, origin, vector, axis2, False)
				if valid ==1 :
					if solution[0] == solution[1]:
						print "Check2, Error! This is not supposed to happen, both the solutions are equal."
						sys.exit(0)
					for sol in solution:
						vector[axis2] = sol
						VectorNum = rd.AddToVertices(vector)
						presentvector.append( VectorNum )
						returnvec[0].append( VectorNum )
					print "got solution == ", solution
###########################################

		pastvector = aux.link(pastvector, presentvector, [], False)



	for edgeindex, edge in enumerate( SH.adjmatrix) :
		print >> errorout, edgeindex, "edge == ", edge
	for coorindex, coor in enumerate( SH.adjcoordinates) :
		print >> errorout, coorindex, "coor == ", coor	
	return returnvec