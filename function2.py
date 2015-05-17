#!usr/bin/env python
import numpy as np
import numpy.linalg as linalg
import sys
import fileinput
from auxilary import *
from graph import *
from plot import *

ellips = []		#stores the origin (center) values for diff ellipses
ellarr = []		#stores the A matrix for diff ellipses
adj = []; adjval = []; valcount = 0; num = 0;dim = 0;primA = []; inout = [];

def method(travaxis,vector, maxminus, crit, presslice):	#traversal axis
	global valcount,num,dim,adj,adjval;
	basecase = 0
	if (len(vector) - travaxis == 2):
		basecase = 1 		#indeed the base case

	metreturn = []
	pastvector = []
	
	startpt = maxminus			#-primA[travaxis]; 	
	iterate = 20
	vector[travaxis] = startpt # don't have to do  this as vector has come as already set
	tempv = vector[:]
	adjval.append(tempv)
	metreturn.append(valcount); #VERY IMP
	pastvector.append(valcount); valcount+=1;		#DO NOT APPEND BEFORE SENDING!!!!!
	cpt = 1 #to start things of for now.
	
	count = 0
	while 1: 	
		count += 1;
		if (count==iterate): 
			break;
		#put checks for base case later
		vector[travaxis] = startpt + ((startpt) * (-2) * count )/(iterate-1) 
		print crit, vector[travaxis]

		if crit == None:
			pass
		elif (vector[travaxis]<crit[travaxis] and crit[travaxis] < startpt + ((startpt)*(-2)*(count+1))/(iterate-1)):
			vector[travaxis] = crit[travaxis]
			count -= 1
		V = []
		itercpt = 0	
		atleastone = 0		#tracks if the plane intersects aatleast one ellipse
		#print "ellarr==", ellarr
		for index,ellipse in enumerate(ellarr):
			origin = ellips[index]
			solution, valid = intersect(ellipse, origin,vector, (travaxis+1) +1)		#always check maxima along the next axis
			#print "vect, solution == ", vector, solution
			if valid == 1:
				#put checks to check that each point is actully valid and not entering another ellipse!
				inout[index] = 1;
				atleastone =1;

				if solution[1] != solution[0]:
					itercpt += 2;					
					vector2 = vector[:]; 			vector2[travaxis+1] = solution[0]
					vector3 = vector[:]; 			vector3[travaxis+1] = solution[1]
					V.append(vector2); 	V.append(vector3);
				else:
					itercpt += 1;					
					vector2 = vector[:]; 			vector2[travaxis+1] = solution[0]
					V.append(vector2);
			else:
				pass
				#print "Invalid for count== ", count, " vector== ", vector, "ellipse== ", index
		#print "count == ", count
		if count == iterate-1:
			#print "here, V = ",V,pastvector
			pastvector, valcount , adjval, adj = link(pastvector, V, valcount, adjval, adj)
			metreturn.append(valcount-1);
			return metreturn
		if atleastone == 0:
			continue;			#it means it does not need to connect this point
		
		critical = None
		temvector = pastvector
		if (itercpt > cpt or itercpt < cpt):	
			#print "\nChange in cpt, travaxis== ",travaxis
			cpt = itercpt
			if presslice == 1:
				presslice = 0

			else:

				if dim - travaxis == 2: #check for 2d
					print "yahan h hum"
					pastvector, valcount, adjval, adj = link(pastvector, V, valcount, adjval, adj)
					metreturn.append(valcount-1)
					print "Done"
					crit = None
					continue
				solution, valid = cpfind(ellipse,vector,travaxis,travaxis+1)
				#print "solution, vector and valid == ", solution,vector, valid
				
				tempv1 = vector[:]
				v2 = vector[:]
				#print tempv1, solution, valid
				if valid ==1:
					tempv1[travaxis] = solution[0]		#for the first
					one = np.array(tempv1)
					two = np.array( v2 )
					three = linalg.norm( one - two )
					#print "one and two", one, two
					tempv1[travaxis] = solution[1]
					one  = np.array(tempv1)
				 	four = linalg.norm( one - two )
				 	if three<four:
				 		tempv1[travaxis] = solution[0]
				 	#else leave it as it is
					max2minus = 0;			#need to calculate maxminus
					solution, valid = intersect(ellarr[0], tempv1,ellips[0], (travaxis+1) +1)		#always check maxima along the next axis
					if(solution[0]  < 0):
						max2minus = solution[0]
					else:
						max2minus = solution[1]
					critical = tempv1[:]
					
					temvector = method(travaxis+1 , tempv1, max2minus, critical, 1)
					past2vector, valcount, adjval, adj = link(temvector, V, valcount, adjval, adj)
					adj = link_special(temvector,pastvector,adjval,adj)
					pastvector = past2vector[:]
					continue
					
		#print "Before link, P, V", pastvector,V
		pastvector, valcount, adjval, adj = link(temvector, V, valcount, adjval, adj)
		
		
	

dim = 0
num = 0
inp1 = open("take_into_account.txt","r")#inp = open("input.txt","r")
inp = inp1.readlines()
for index, l in enumerate(inp):
	line = l.strip()
	parts = line.split()
	if len(parts) == 2:
		#dim is for the dimensions being used
		#num is the number of ellipses being checked
		num = int(parts[0])
		dim = int(parts[1])
	else:
		if dim == 0:
			print "error in input as dim =0"
			sys.exit()		
		origin = []		#first detail in line is origin
		for i in range(0,dim):
			origin.append(float(parts[i]))
		ellips.append(origin) #store the center values for the nth ellipse
		
		carry = dim;	#carry tracks how many have been read and where to read from next
		
		Siglist = [] #holds the U matrix in 2d list form
		for i in range(0,dim):
			origin = []	#reusing this variable/list
			for j in range(0,dim):
				origin.append(float(parts[j*dim + i + carry]))
			Siglist.append(origin)

		Sigmalist = [] 	#holds the sigma matrix
		carry = dim + dim*dim
		for i in range(0,dim):
			term = float(parts[i + carry])
			inver = 1/(term*term)
			Sigmalist.append(inver)
			if index == 1:
				primA.append(term)
		
		eye = np.identity(dim,dtype = float)
		for i in range(0,dim):
			eye[i][i] = Sigmalist[i] 		#now eye stores the SIGMA matrix
		U = np.array(Siglist,float)
		trans = U.T
		#print "U and trans = " , U, trans		#multiply the three
		mult1 = np.dot(U,eye)
		mult2 = np.dot(mult1,trans)		#A = U * eye * U.T
		#print "A for this iteration = ", mult2
		ellarr.append(mult2)

travaxis = 0
vector = []
for x in range(0,dim):
	vector.append(0)
for x in range(0,len(ellarr)):
	inout.append(0)
maxminus = -primA[travaxis];
critical = None;
presslice = 1
method(travaxis,vector , maxminus, critical, presslice)

#print "\n\nadj== ", adj
tree = main(valcount, adj)
print tree
#for x in range(0,adj):

plotthis(num,ellips,ellarr,adj,adjval)