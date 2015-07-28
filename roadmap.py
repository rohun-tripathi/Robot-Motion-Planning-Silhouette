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

from copy import copy, deepcopy


#The main method
def CreateRoad(travaxis, critical, parentvector ,originlist, ellipselist, debug = False):
	#Start of initialize	
	errorout = open("debuginfo.txt", "w")
	print originlist, ellipselist
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
	CriticalX, 	CriticalYZ = cp.axisrange(Criticalpts, False	)

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

		
		#Considerlist - holds the list of the ellipses that need to be considered for the intersect function
		#The returned list is of dimension n-1. If the vaule is 1, the ellipse is to be considered, else not
		considerlist = cp.EllSliceintersect(CriticalX, nextslice)

		#The returned list is of dimension n-1.
		considerYZ = cp.EllUnderConsider(considerlist , CriticalYZ, CriticalX,  nextslice, travaxis, False)

		slicevector[travaxis] = nextslice 			#No need for the old slice information?
				#Technically if there were any CPs they should have been dealt with by this stage
#############################################
		#The Intersect and creation of presentvector
#############################################

		#first for the first ellipse
		for axis2 in range(1,2):# (travaxis+1, SH.dim):
			solution, valid = aux.intersect(ellipselist[0], originlist[0], slicevector[(travaxis):], axis2)
			#print >> errorout, "Sent == ", ellipselist[0], originlist[0], slicevector[(travaxis):], axis2
			vector = slicevector[:]

			if valid ==1 :		#Later check also that it is not inside any of the other ellipsoids
				if solution[0] == solution[1] and firstlink!= 1:
					print "Error, this is not supposed to happen, both the solutions are equal. But letting it be for Now. iteration == ", iteration
					print solution, slicevector
					#sys.exit(0)
				for sol in solution:
					print >> errorout, "For axis == ", axis2, " for vector == ",vector, " found solution == ", sol
			
					vector[travaxis + 1] = sol
					
					VectorNum = rd.AddToVertices(vector[:])
					presentvector.append( VectorNum )
					#returnvec[0].append( VectorNum )
				#Some processing needs to be done
			

		for index, listitem in enumerate(considerlist):
			if listitem == 0 :
				continue
			#else: Now the rest of the code for this loop

			# We have to find the correct ellipse and remember to account for the first ellipse being primA 			
			if debug == True :print "Working for obstacle, index no. ", index
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
			if debug == True :print "For the slicevector ", slicevector, "gotten the following vector to run intersect along, ", vector
			
			#Run Intersect
			for axis2 in range(1,2):# (travaxis+1, SH.dim):
				solution, valid = aux.intersect(ellipse, origin, vector[(travaxis):], axis2, False)
				if valid ==1 :
					if solution[0] == solution[1]:
						print "Check2, Error! This is not supposed to happen, both the solutions are equal."
						sys.exit(0)
					for sol in solution:
						vector[travaxis + 1] = sol
						VectorNum = rd.AddToVertices(vector)
						presentvector.append( VectorNum )
						#returnvec[0].append( VectorNum )
					if debug == True :print "got solution == ", solution
###########################################
#Intersect Done

		#Last Step of this iteration
		pastvector = aux.link(pastvector, presentvector, [], False)

		#Now Begins Work of Critical points
		nextslice = rd.getNextSlice(startpt, iteration + 1)		#nextslice is a value of next iteration
		#ActiveCP and restCP have the value for any critcal points found at this stage. The active is along the travaxis and restCP is for recursion
		activeCP, restCP = cp.RecursionPoints(CriticalX, CriticalYZ, slicevector[travaxis], nextslice, False)
		#To find which critical points from higher level have to be attached at this stage
		CritAtThisSlice = cp.RecurCheck (critical, slicevector[travaxis], nextslice, False)
		
		if len(activeCP) + len(CritAtThisSlice)> 1:
			print "System Limitation. We are just working with one critical point between two slices presently. A solution could be to increase the iterate variable, so slices come closer"
			print "It could be be the critical points from higher dimension and this one are crowding together."
			sys.exit(0)
		if len(activeCP) + len(CritAtThisSlice)== 0: continue
		
		#We work on the CritAtThisSlice one first
		# if len(CritAtThisSlice) > 0:						#If count is greater than 0, first recursion worked
		# 	print activeCP, restCP, CritAtThisSlice, slicevector
		# 	raw_input("reached here")

		if (SH.dim - travaxis) == 2: 						#Reached the 2D case - Base Case
			CPvector = slicevector[:]						#CPvector will be manipulated to represent the critical val
			
			# print "CritAtThisSlice at slicevector == ", CritAtThisSlice, slicevector, activeCP
			try:
				point = CritAtThisSlice[0][0]					#System Limitation, only account for one critical point
			except IndexError:
				point = [activeCP[0][0]] + restCP[0][0]
			for index, term in enumerate(point):
				CPvector[travaxis + index] = term			#CPvector will be added to the Graph after this
			VectorNum = rd.AddToVertices(CPvector)
			try:
				startendflag = CritAtThisSlice[0][1]
			except IndexError:
				startendflag = restCP[0][1]
			if startendflag == "start": returnvec[1].append(VectorNum)
			elif startendflag == "end": returnvec[2].append(VectorNum)
			else:
				print "There is an error in the travaxis == ", travaxis, "The Criticalpt here does not have appended text (start/end)" 
				sys.exit(0)

			aux.link(pastvector, [VectorNum], [], False)	#Add this VectorNum to the network
			pastvector.append(VectorNum)					#Added to pastvectors as will be needed in next slice	
			continue										#Done for this iteration


		else:												#Not 2D. We have to prepare to call the lower dimension
			nextslice = activeCP[0][0]						#The travaxis value for the recursion slice
			considerlist = cp.EllSliceintersect(CriticalX, nextslice)
			considerYZ = cp.EllUnderConsider(considerlist , CriticalYZ, CriticalX,  nextslice, travaxis, False)
			slicevector[travaxis] = nextslice 			

			#Get the Ellipselist and the originlist for the lower dimensions
			RecursionEll, RecursionOri = cp.ReduceEllipsoids(considerlist, considerYZ, slicevector, travaxis, deepcopy(ellipselist), deepcopy(originlist), False)

			#obtainedvec is the otherside of returnvec
			obtainedvec = CreateRoad(travaxis+1, restCP, slicevector ,RecursionOri, RecursionEll, debug = False)
			
			presentvector = obtainedvec[0][:]
			presentvector.extend(obtainedvec[2][:])		#For Doubts regarding the obtainedvec/returnvec refer to README.md
			aux.endRecur_link(pastvector, presentvector, False)

			pastvector = obtainedvec[0][:]
			pastvector.extend(obtainedvec[1][:])		#Preparing the real pastvector for the next slice. 


	for edgeindex, edge in enumerate( SH.adjmatrix) :
		print >> errorout, edgeindex, "edge == ", edge
	for coorindex, coor in enumerate( SH.adjcoordinates) :
		print >> errorout, coorindex, "coor == ", coor	
	return returnvec