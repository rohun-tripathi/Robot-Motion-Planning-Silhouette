import numpy as np
import sys
import numpy.linalg as linalg
import math

# Find the intersect point along the nth axis
# A is the ellipse and X is the array with the points X = [x1, x2 ....]
def intersect(A,O,X,n):
	
	#print "A== " ,A, X
	rank=len(X)
	templist = [[0,0] for i in range(0,rank)] #will be used for the first multiplication

	for i in range(0, rank):
		prod = 0;
		for j in range(0,rank):
			if j == n:
				templist[i][1]=A[j][i]
				prod -= A[j][i] * O[j]
			else:
				prod += A[j][i] * (X[j] - O[j] )
		templist[i][0]=prod
	#y^2
	component = [0,0,0]			#represents a, b and c component
	for i in range(0, rank):
		if i == n:
			component[0] = templist[i][1]
			component[1] += templist[i][0]
		else:
			component[2] += templist[i][0] * X[i]
			component[1] += templist[i][1] * X[i]
	component[2] -= 1 			#account for the ellipse equation	
	#print "component == " ,component
	solution,valid = solver(component)
	for x in range(0,len(solution)):
		solution[x] += O[n-1]
	return solution,valid

# Find the critical point along the cth axis for the nth axis (c is x and n is y)
# A is the ellipse and X is the array with the points X = [x1, x2 ....]
def cpfind(A,O,X,x,y):
	rank=len(X)
	templist = [[0,0,0] for i in range(0,rank)] #order: constant, x, y

	for i in range(0, rank):
		prod = 0;
		for j in range(0,rank):
			if j == x:
				templist[i][1] += A[j][i]
				templist[i][0] -= A[j][i] * (O[x])
			elif j == y:
				templist[i][2] += A[j][i]
				templist[i][0] -= A[j][i] * (O[y])
			else:
				templist[i][0] += A[j][i] * (X[j]-O[j])
		
	#print "templist== ", templist

	component = [0,0,0,0,0,0]			#represents a, b and c component
		
	for i in range(0,rank):
		if i == x:
			component[3] += templist[i][0]
		elif i == y:
			component[1] += templist[i][0]
		else:
			component[5] += templist[i][0] * (X[i]-O[i])
	#doing for n
	for i in range(0,rank):
		if i == x:
			component[4] += templist[i][2]
		elif i == y:
			component[0] += templist[i][2]
		else:
			component[1] += templist[i][2] * (X[i]-O[i])
	#doing for c
	for i in range(0,rank):
		if i == x:
			component[2] += templist[i][1]
		elif i == y:
			component[4] += templist[i][1]
		else:
			component[3] += templist[i][1] * (X[i]-O[i])

	p = component[0]
	r = component[1]
	delt = component[2]
	q = component[3]
	s = component[4]
	t = component[5] - 1

	for i in range(0,len(templist)):
		t -= templist[i][0] * O[i]
	for i in range(0,len(templist)):
		r -= templist[i][1] *  O[i]
	for i in range(0,len(templist)):
		q -= templist[i][2] *  O[i]

	comp = [0,0,0]
	comp[0] = (s * s) - (4 * delt * p)
	comp[2] = q *q - 4* delt * t
	comp[1] = 2*q* s - 4* delt * r
	#print "componenet, comp== ",component, comp
	solution,valid = solver(comp)
	# for x in range(0,len(solution)):
	# 	solution[x] += O[n-1]
	
	return solution,valid


# function link - 
# links the the present vector points in V and past vector points in P
# creates adjacency list adj. adjval stores the point in N dimension for each point index
# returns the index of vectors in V corresponding to numbers in adjval 

def link_special(P, V, adjval, adj):			#past and now vector
	
	#print "p and V== ", P, V
	for i in range(0,len(V)):		
		inf = 0; val = 23456789	#high value
		for j in range(0,len(P)):
			one = np.array( adjval[V[i]])
			two = np.array( adjval[P[j]] )
			three = linalg.norm( one - two )
		 	#print "1, P[j] , 2 ,3, val == ",one ,P[j], two,  three, val, "\n"
		 	#print "adjval == ", adjval
		 	if val >= three :
		 		val = three
		 		inf = P[j]
		#print "inf== ", inf
		adj.append([inf,V[i], val])	#ULTA append!
		
	return adj

def link(P, V, CV, valcount, adjval, adj):			#past and now vector
	point = []

	#print "p and V and CV== ", P, V, CV
	for i in range(0,len(V)):		
		inf = 0; val = 23456789	#high value
		for j in range(0,len(P)):
			one = np.array(V[i])
			two = np.array( adjval[P[j]] )
			three = linalg.norm( one - two )
		 	#print "1, P[j] , 2 ,3, val == ",one ,P[j], two,  three, val, "\n"
		 	#print "adjval == ", adjval
		 	if val >= three :
		 		val = three
		 		inf = P[j]
		tempi = V[i][:]
		#print "inf== ", inf
		adjval.append(tempi)
		adj.append([valcount,inf, val])	#the weight is val
		point.append(valcount)
		valcount += 1
	for i in range(0,len(CV)):		
		inf = 0; val = 23456789	#high value
		for j in range(0,len(P)):
			one = np.array(CV[i])
			two = np.array( adjval[P[j]] )
			three = linalg.norm( one - two )
		 	#print "1, P[j] , 2 ,3, val == ",one ,P[j], two,  three, val, "\n"
		 	#print "adjval == ", adjval
		 	if val >= three :
		 		val = three
		 		inf = P[j]
		tempi = CV[i][:]
		#print "inf== ", inf
		adjval.append(tempi)
		adj.append([valcount,inf, val])	#the weight is val
		point.append(valcount)
		valcount += 1
	return point, valcount, adjval, adj

def complete_link(P, V, CV, valcount, adjval, adj):			#past and now vector
	point = []

	#print "p and V and CV== ", P, V, CV
	for i in range(0,len(V)):		
		inf = 0; val = 23456789	#high value
		
		tempi = V[i][:]
		adjval.append(tempi)
		
		for j in range(0,len(P)):
			one = np.array(V[i])
			two = np.array( adjval[P[j]] )
			three = linalg.norm( one - two )
		 	#print "1, P[j] , 2 ,3, val == ",one ,P[j], two,  three, val, "\n"
		 	#print "adjval == ", adjval
			adj.append([valcount,inf, val])	#the weight is val
		point.append(valcount)
		valcount += 1
	for i in range(0,len(CV)):		
		inf = 0; val = 23456789	#high value
		
		tempi = CV[i][:]
		adjval.append(tempi)
		
		for j in range(0,len(P)):
			one = np.array(CV[i])
			two = np.array( adjval[P[j]] )
			three = linalg.norm( one - two )
		 	#print "1, P[j] , 2 ,3, val == ",one ,P[j], two,  three, val, "\n"
		 	#print "adjval == ", adjval
		 	if val >= three :
		 		val = three
		 		inf = P[j]
			adj.append([valcount,inf, val])	#the weight is val
		point.append(valcount)
		valcount += 1
	return point, valcount, adjval, adj

def solver(component):
	solution = []

	a=  component[0]; b = component[1]; c = component[2]
	d = (b**2) - (4*a*c) 	#discriminant
	if d<0:
		return solution,0
	# find two solutions
	sol1 = (-b-math.sqrt(d))/(2*a)
	sol2 = (-b+math.sqrt(d))/(2*a)

	solution.append(sol1)
	solution.append(sol2)
	return solution,1