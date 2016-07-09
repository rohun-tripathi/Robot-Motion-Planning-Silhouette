import math
import sys

import numpy as np
import numpy.linalg as linalg

import rdfunctions
import shared as SH


#####################
# links the the present vector points of the last slice of an obstacle to the previous slice
# Mimics behaviour of complete_link, which is not other wise captured in last point of obstacle
# diffList is the list of all points that have not been connected, and need to be connected 
def endRecur_link(P, V, debug=False):  # P - past and V - Presentvector
    linked = []
    for i in range(0, len(V)):
        inf = 0;
        val = 23456789  # high value
        for j in range(0, len(P)):
            one = np.array(SH.adjcoordinates[V[i]])
            two = np.array(SH.adjcoordinates[P[j]])
            three = linalg.norm(one - two)

            if val >= three:
                val = three
                inf = P[j]
        if debug == True: print "For the ", i, "th presentvector, the inf== ", inf
        rdfunctions.addToRoadmap(V[i], inf, val, "link")  # the weight is val
        linked.append(inf)
    diffList = [x for x in P if x not in linked]
    if debug == True: print "diffList, linked, pastvector == ", diffList, linked, P
    linkPresentAndPastVector(V, diffList, [], debug)


#####################
# function Complete_link - 
# links the the present vector points in V and past vector points in P
# creates adjacency list adj. adjval stores the point in N dimension for each point index
# returns the index of vectors in V corresponding to numbers in adjval 
# The functionality pf link has been understood and the others are just variations as per needs
def complete_link(P, V, debug=False):  # past and now vector
    if len(V) > 1:
        print "In Complete Link, Pastvector and Present-vector== ", P, V
        print "The length of Pastvector is greator than 1. \n Should this have happened? Check it.\nExiting"
        sys.exit(0)

    for i in range(0, len(V)):
        inf = 0;
        val = 23456789  # high value
        for j in range(0, len(P)):
            one = np.array(SH.adjcoordinates[V[i]])
            two = np.array(SH.adjcoordinates[P[j]])
            three = linalg.norm(one - two)

            if debug == True: print "In complete link, For the ", i, "the presentvector, the inf== ", inf
            rdfunctions.addToRoadmap(V[i], P[j], val, "complete_link")  # the weight is val
    return V


#####################
# function link - 
# links the the present vector points in V and past vector points in P
# creates adjacency list adj. adjval stores the point in N dimension for each point index
# returns the index of vectors in V corresponding to numbers in adjval 
# The functionality pf link has been understood and the others are just variations as per needs
def linkPresentAndPastVector(P, V, CV, debug=False):  # past and now vector
    if debug == True: print "Pastvector and Present-vector and CV== ", P, V, CV
    # move the two below to the same function
    for i in range(0, len(V)):
        inf = 0
        val = 23456789  # high value
        for j in range(0, len(P)):
            one = np.array(SH.adjcoordinates[V[i]])
            two = np.array(SH.adjcoordinates[P[j]])
            three = linalg.norm(one - two)

            if val >= three:
                val = three
                inf = P[j]
        if debug == True: print "For the ", i, "th presentvector, the inf== ", inf
        rdfunctions.addToRoadmap(V[i], inf, val, "link")  # the weight is val

    for i in range(0, len(CV)):
        inf = 0;
        val = 23456789  # high value
        for j in range(0, len(P)):
            one = np.array(SH.adjcoordinates[CV[i]])
            two = np.array(SH.adjcoordinates[P[j]])
            three = linalg.norm(one - two)

            if val >= three:
                val = three
                inf = P[j]
        if debug == True: print "For the ", i, "th CV, the inf== ", inf
        rdfunctions.addToRoadmap(CV[i], inf, val, "link")  # the weight is val
    V.extend(CV)
    return V


#############################
# Find the intersect point along the nth axis
# A is the ellipse
# X is the array with the points X = [x1, x2 ....]
# O is the array with the origin
# n is one of the axes after the travaxis

def getIntersectionPointAlongtravAxis(A, O, X, n, debug=False):
    # This is important
    # n = 0
    # This marks the change in strategy from last model to this one
    # In the last model, the traversal went downward, in this the ellipses go downward, they are calculated for the lower dimentsion
    # This means that the axis of traversal stays "0" for the present ellipse
    rank = len(X)

    if debug == True: print "intersect called, axis = ", n
    if debug == True: print "A and X and O == ", A, X, O
    if debug == True: print "Rank == ", rank

    templist = [[0, 0] for i in range(0, rank)]  # will be used for the first multiplication

    for i in range(0, rank):  # This will be the column in the A matrix
        prod = 0;
        for j in range(0, rank):  # Go thru the calculations for jth matrix
            if j == n:
                templist[i][1] = A[j][i]  # The templist[x][1] is the coeff of the variable being solved
            # prod -= A[j][i] * O[j]
            else:
                prod += A[j][i] * (X[j] - O[j])  # Product is the constant from this multiplication
            if debug == True: print "Prod == ", prod
        templist[i][0] = prod
        if debug == True: print "Templist == ", templist
    # y^2
    component = [0, 0,
                 0]  # represents a, b and c component of the equation. It is independent of the dimensions of the ellipse
    for i in range(0, rank):
        if i == n:
            component[0] = templist[i][1]  # adds to the coeff of X^2, the variable being solved for
            component[1] += templist[i][0]  # adds to the coeff of X, the variable being solved for
        else:
            component[2] += templist[i][0] * (X[i] - O[i])  # preparing the c
            component[1] += templist[i][1] * (X[i] - O[i])  # adds to the coeff of X, the variable being solved for

    if debug == True: print "Component == ", component

    component[2] -= 1  # account for the ellipse equation

    solution, valid = solver(component)  # Solves ax^2 + bX + c = 0
    if debug == True: print "Solution and valid after function solver", solution, valid
    if debug == True: print "Substracting with the center, value of O[n-1], which is == ", O[n - 1]
    if debug == True: print "Where n is == ", n
    for x in range(0, len(solution)):
        solution[x] += O[n]  # This is important. this is because we solved for (Xn - On), not just for Xn.

    if debug == True: print "The returned values, solution and valid == ", solution, valid
    if debug == True: gap = raw_input("Done with Intersect, Press Enter...")
    return solution, valid


#######################
# This function solves the algebraic equation ax^2 + bx + c = 0
# Returns the solution and a valid flag to indicate if the solution is valid or not
def solver(component):
    solution = []

    a = component[0];
    b = component[1];
    c = component[2]
    d = (b ** 2) - (4 * a * c)  # discriminant
    if d < 0:
        return solution, 0  # No Solution
    # find two solutions
    sol1 = (-b - math.sqrt(d)) / (2 * a)
    sol2 = (-b + math.sqrt(d)) / (2 * a)

    solution.append(sol1)
    solution.append(sol2)
    return solution, 1



    # def link_special(P, V, adjval, adj):			#past and now vector

    # 	#print "p and V== ", P, V
    # 	for i in range(0,len(V)):
    # 		inf = 0; val = 23456789	#high value
    # 		for j in range(0,len(P)):
    # 			one = np.array( adjval[V[i]])
    # 			two = np.array( adjval[P[j]] )
    # 			three = linalg.norm( one - two )
    # 		 	#print "1, P[j] , 2 ,3, val == ",one ,P[j], two,  three, val, "\n"
    # 		 	#print "adjval == ", adjval
    # 		 	if val >= three :
    # 		 		val = three
    # 		 		inf = P[j]
    # 		#print "inf== ", inf
    # 		adj.append([inf,V[i][:], val, "link_special"])	#ULTA append!

    # 	return adj

    # def complete_link(P, V, CV, valcount, adjval, adj):			#past and now vector
    # 	point = []

    # 	#print "p and V and CV== ", P, V, CV
    # 	for i in range(0,len(V)):
    # 		inf = 0; val = 23456789	#high value

    # 		tempi = V[i][:]
    # 		adjval.append(tempi)

    # 		for j in range(0,len(P)):
    # 			one = np.array(V[i])
    # 			two = np.array( adjval[P[j]] )
    # 			three = linalg.norm( one - two )
    # 		 	#print "1, P[j] , 2 ,3, val == ",one ,P[j], two,  three, val, "\n"
    # 		 	#print "adjval == ", adjval
    # 			adj.append([valcount,inf, val, "complete_link"])	#the weight is val
    # 		point.append(valcount)
    # 		valcount += 1
    # 	for i in range(0,len(CV)):
    # 		inf = 0; val = 23456789	#high value

    # 		tempi = CV[i][:]
    # 		adjval.append(tempi)

    # 		for j in range(0,len(P)):
    # 			one = np.array(CV[i])
    # 			two = np.array( adjval[P[j]] )
    # 			three = linalg.norm( one - two )
    # 		 	#print "1, P[j] , 2 ,3, val == ",one ,P[j], two,  three, val, "\n"
    # 		 	#print "adjval == ", adjval
    # 		 	if val >= three :
    # 		 		val = three
    # 		 		inf = P[j]
    # 			adj.append([valcount,inf, val, "complete_link"])	#the weight is val
    # 		point.append(valcount)
    # 		valcount += 1
    # 	return point, valcount, adjval, adj
