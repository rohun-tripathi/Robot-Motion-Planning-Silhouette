## This the introduction, a readme of sorts
#ellips is a list of lists. List of the origins of the nth ellipse
#ellarr stores the A matrix for diff ellipses

#primA stores the terms for the primary axis. The one to traverse along on the first iteration

#Presslice is to avoid any kind of recursion in first slice due to change in Critical points, the starting point

#imports
#!usr/bin/env python
import numpy as np
import numpy.linalg as linalg
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import sys
import time
import random
import fileinput

from auxilary import *
from graph import *
from plot import *
from complement import *




ellips = []		#stores the origin (center) values for diff ellipses
ellarr = []		#stores the A matrix for diff ellipses
adj = []; adjval = []; valcount = 0; num = 0;

#	valcount keeps the count of the value of the poitn being added to the adjacency tree
#	the adj is the list of edges that being created
#	the adjval gives the dimension coordinates (for plotting for the nth point)

dim = 0;primA = []; inout = [];
#primA stores the terms for the primary axis. The one to traverse along on the first iteration
#The Inout gives the status of whether a particular slice is inside or outside of the ellipse in question

iterate = 100 	#Something like space parts to complete traversal

#The main method

def method(travaxis,vector, startpt, crit, presslice):	#traversal axis
	global valcount,num,dim,adj,adjval;
	

	
	vector[travaxis] = startpt	
	inplace = vector[:]		#inplace used only in next line
	adjval.append(inplace)	#first point appended
	
	metreturn = []					#Will be returned as "pastvector" for next slice
	metreturn.append(valcount); 	#Appended one end. Append the other end in the last iteration
	
	pastvector = []			#From the previous iteration, slice
	pastvector.append(valcount); valcount+=1;		#DO NOT APPEND BEFORE SENDING!!!!!
	cpt = 1 				#Critical point		
	
	count = 0				#loop over variable "iterate"
	while 1: 	
		count += 1;
		if (count==iterate): 
			break;

		#The different slices at the different positions along the travaxis
		vector[travaxis] = startpt + ((startpt) * (-2) * count )/(iterate-1) 	
		V = []				#Holds vector for this slice. will be connected to pastvector
		CV = []
		itercpt = 0			#Critical Points in this iteration
		atleastone = 0		#tracks if the plane intersects atleast one ellipse and if linking has to be done

		if crit != None:	#If the critical points from higher recursion is not ziltch
			if (vector[travaxis] >= crit[travaxis] and crit[travaxis] > vector[travaxis] - ((startpt)*(-2))/(iterate-1)):
				#print "Found ", crit[travaxis], "between ", vector[travaxis], vector[travaxis] - ((startpt)*(-2))/(iterate-1)
				vector[travaxis] = crit[travaxis]
				crit_temp = crit[:];
				CV.append(crit_temp);
				itercpt += 1;
				
				count -= 1 					#If the crit is not ziltch, then goto the slice in question and also reduce this count coz it has to be done later
		
		ellinques = -1						#assuming at one instant, only one new ellipse is added to the network to be encompassed
		axis2 = travaxis;
		for fg in range(0,2):				#This 0,2 was because of limit to 3 dimensions in this case
			axis2 = axis2 + 1
			if not axis2 < dim:				#Axis2 should be less tham dim
				break;


			for index,ellipse in enumerate(ellarr):		#looping over each ellipse for the given slice
				origin = ellips[index]					#print "\norigin and vector = ", origin, vector
				solution, valid = intersect(ellipse, origin, vector, axis2)		#maxima along next axis
				if valid == 1:
					atleastone = 1;			#used later, it means it does or does not need to link this point

					#check that each point is actully valid and not entering another ellipse!
					#STill have not written that code
					
					if inout[index] == 2:		#Comments on inout in complement	
						inout[index] = 1
					if inout[index] == 0:
						inout[index] = 2
						ellinques = index
					
					if solution[1] != solution[0]:	#two point as compared to one in the caseof a critcial slice!
						itercpt += 2;				#iteration cpts

						vector2 = vector[:]; 			vector2[axis2] = solution[0]
						vector3 = vector[:]; 			vector3[axis2] = solution[1]
						V.append(vector2); 	
						V.append(vector3); 	
						# if vector2 != crit:
						# 	fvalid = f(dim,num,vector2,ellips,ellarr)	#f(..) always returns 1
						# 	if fvalid == 1:
						# 		V.append(vector2); 	
						# if vector3 != crit:
						# 	fvalid = f(dim,num,vector3,ellips,ellarr)	#f(..) always returns 1
						# 	if fvalid == 1:
						# 		V.append(vector3);

					else:							#Solutions equal in casen of a critcial slice!
						itercpt += 1;				#one in the caseof a critcial slice!
						vector2 = vector[:]; 			vector2[axis2] = solution[0]
						if vector2 != crit:
							CV.append(vector2); 			#NO Doubt about validity	
															#CV used here to disambiguate from the other points (V) added to roadmap this iteration

				else:								#If NOT valid
					if inout[index] == 1 or inout[index] == 2:	#Comments on inout in complement, State 1 or 2 means inside.
						inout[index] = 3
						ellinques = index
					else:
						inout[index] = 0
					#print "Invalid for count== ", count, " vector== ", vector, "ellipse== ", index
		
		axis2 =travaxis + 1;			#Back in main while
		
		if itercpt != cpt: print "\nCHANGE : solution, vector, V, axis2, itercpt, presslice", solution, vector, V, axis2, itercpt, presslice
		#This ^ only happens if some intersections of the particular slice with atleast one ellipse
			
		
		if count == iterate-1:						#Last iteration. Return to higher resursion
			pastvector, valcount , adjval, adj = complete_link(pastvector, V, CV, valcount, adjval, adj)
			metreturn.append(valcount-1);			#Appended to metreturn the last point of previous iteration
			return metreturn

		if atleastone == 0:			#it means it does not need to connect in this iteration, no linking
			continue;				#This is the main while(1) loop
		
		##########################
		#The next part is of linking in different manners. Simple linking (link)
		#reverse linking (link_special) for Critical slice cases
		# In the case of the critical slice the approppriate recursion is called unless the system is in the base case(2D)
		# In base case connected by point

		critical = None
		temvector = pastvector 		#used to connect points from lower recursion 
		if (itercpt != cpt):	
			#print "\nChange in cpt, travaxis== ",travaxis
			cpt = itercpt
			if presslice == 1:		#first slice in recursion, the starting point
				presslice = 0
			
			else:
				presslice = 1
				if dim - travaxis == 2: #check for 2d
					print "\nIN 2d, itercpt, vector, pastvector == ",itercpt, vector, pastvector
					pastvector, valcount, adjval, adj = link(pastvector, V,CV ,valcount, adjval, adj)
					metreturn.append(valcount-1)		#last me cp hi append hoga
					#print "\n\nappended to metreturn , presslice == ", valcount, presslice
					#print "2D linking done\n"
					crit = None
					continue
				else:
					if crit == None:
						solution, valid = cpfind(ellarr[ellinques],ellips[ellinques],vector,travaxis,axis2)	#crit plane creator
						
						tempv1 = vector[:]
						v2 = vector[:]
						
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
						 		tempv1[travaxis] = solution[0]			#Find the closer criticl point for respective ellipse to the vector
						 												#else leave it as it is
							
							max2minus = 0;								#need to calculate maxminus
							solution, valid = intersect(ellarr[0],ellips[0], tempv1, (axis2))			#always check maxima along the next axis
							max2minus = solution[1]						#find the point on the primry to start recursion
							if(solution[0]  < 0):					
								max2minus = solution[0]
							
							critical = tempv1[:]
							critical[axis2] = ellips[ellinques][axis2]			#crit creation
							
							temvector = method(axis2 , tempv1, max2minus, critical, 1)		#RECURSION
							past2vector, valcount, adjval, adj = link(temvector, V, CV, valcount, adjval, adj)
							adj = link_special(temvector,pastvector,adjval,adj)
							pastvector = past2vector[:]
					else:
						tempv1 = crit[:]
						
						max2minus = 0;			#need to calculate maxminus
						solution, valid = intersect(ellarr[0],ellips[0], tempv1, (axis2))		#always check maxima along the next axis
						if(solution[0]  < 0):
							max2minus = solution[0]
						else:
							max2minus = solution[1]
						
						temvector = method(axis2 , tempv1, max2minus, crit, 1)
						past2vector, valcount, adjval, adj = link(temvector, V, CV, valcount, adjval, adj)
						adj = link_special(temvector,pastvector,adjval,adj)
						pastvector = past2vector[:]
				crit = None 
				continue
						
		#print "Before link, P, V", pastvector,V
		pastvector, valcount, adjval, adj = link(temvector, V,CV, valcount, adjval, adj)
		#time.sleep(.2)

#---------------------------------------------------------
#Code region to calculate the Gram Smidt Orthogonalization
#----------------------------------------------------------

def modulus(v1):
	mod = float(np.sqrt(float(np.dot(v1,v1))))
	return [x/mod for x in v1]

def gs_cofficient(v1, v2):			#Gram Smidt Coefficient
	return np.dot(v2, v1) / np.dot(v1, v1)
 
def multiply(cofficient, v):
	return [x * cofficient for x in v]
 
def proj(v1, v2):
	return multiply(gs_cofficient(v1, v2) , v1)
 
def gs(X):
	mod_vec = 0.0
	Y = []
	for i in range(len(X)):
		temp_vec = X[i]
		for inY in Y :
			proj_vec = proj(inY, X[i])
			temp_vec = map(lambda x, y : x - y, temp_vec, proj_vec)
		mod_vec = modulus(temp_vec)
		Y.append(mod_vec)
	#print "The final matrix returned from gs Y = ", Y
	return Y
 
#--------------------------------
#this function plots in 3D
#--------------------------------

def createbasis(num):
	# plot
	for turn in range(0,4):
		#print "Enter four points which describes the vector perpendicular to the 3D hyperplane"
		send = []
		filename = ""
		for n1 in range(0,4):
			if n1 == turn:
				filename += str(1) + " "
			else:
				filename += str(0) + " "
				send.append(n1)
		print filename
		time.sleep(.2)
		vector = []
		subspace = []
		for temp in filename.split():
			term = temp.strip()
			vector.append(float(term))
		
		subspace.append(vector)
		#print "Points you entered : " , vector

		#we have put the first one
		#for random
		for i in range(0,3):
			vector = []
			for j in range(0,4):
				term = random.randrange(0,5+1)	#random number between 0 and 6 including them	
				vector.append(float(term))
			subspace.append(vector)
		basis = []
		basis = gs(subspace)
		temp_turn = 200+20+ turn
		plotthis2(num,basis,temp_turn,send)

def plotthis2(num,basis,turn,send):
	#---------------------------------
	#P matrix as described by sir in notes is called Proj
	#---------------------------------
	Projtemp = []  							#contains vectors in plane first and then orthogonal
	for i in range(0,len(basis)):
		Projtemp.append(basis[(i+1) % dim])
	#print "Projtemp = ", Projtemp
	
	Proj = np.array(Projtemp)
	Projtrans = Proj.T
	#print "Proj = ", Proj
	
	ax = fig.add_subplot(turn, projection='3d')
	for k in range(0,num):
	
		center = ellips[k]
		#print "center = ", ellips[k]
		
		A = ellarr[k]
		#print "A from ellipse list = ", A

		#------------------------
		#the projected matrix is Proj*A*Projtrans
		#------------------------

		projected = np.dot(np.dot(Proj,A),Projtrans)
		projcenter = np.dot(Projtrans,center)
		print "Projected and projcenter = " , projected, projcenter
		
		#print "center and projcenter = " , center,  projcenter
		#now we have to select the corner nine points from this projected matrix

		center = []
		A = []
		for i in range(0,3):
			row = []
			for j in range(0,3):
				row.append(projected[i][j])
			A.append(row)
		for i in range(0,3):
			center.append(projcenter[i])

		#This is the part same as before

		U, s, rotation = linalg.svd(A)
		radii = 1.0/np.sqrt(s)

		u = np.linspace(0.0, 2.0 * np.pi, 100)
		v = np.linspace(0.0, np.pi, 100)
		x = radii[0] * np.outer(np.cos(u), np.sin(v))
		y = radii[1] * np.outer(np.sin(u), np.sin(v))
		z = radii[2] * np.outer(np.ones_like(u), np.cos(v))

		for i in range(len(x)):
		    for j in range(len(x)):
		        [x[i,j],y[i,j],z[i,j]] = np.dot([x[i,j],y[i,j],z[i,j]], rotation) + center

		if k==0:
			ax.plot_wireframe(x, y, z,  rstride=4, cstride=4, color='b', alpha=0.2)
		else:
			ax.plot_wireframe(x, y, z,  rstride=4, cstride=4, color='r', alpha=0.2)
	
	adjval_temp = []
	for n1 in adjval:
		n2  = np.array(n1)
		adjval_temp.append(np.dot(Proj,n2))

	for lin in adj:
		print lin
		n = []; m= []; p = [];
		for i in range(0,2):
			n.append(adjval_temp[lin[i]][send[0]])
		for i in range(0,2):
			m.append(adjval_temp[lin[i]][send[1]])
		for i in range(0,2):
			p.append(adjval_temp[lin[i]][send[2]])

		ax.plot_wireframe(n,m,p, color="green", linewidth=2.0, linestyle="-")	
	
def plotthis3():
	# plot
	print "Enter three points which describes the vector perpendicular to the 3D hyperplane"
	filename = str("1")+" "+str("1")+" "+str("1")+" "
	vector = []
	subspace = []
	for temp in filename.split():
		term = temp.strip()
		vector.append(float(term))
	
	subspace.append(vector)
	#print "Points you entered : " , vector

	#we have put the first one
	#for random
	for i in range(0,dim-1):
		vector = []
		for j in range(0,dim):
			term = random.randrange(0,5+1)	#random number between 0 and 6 including them	
			vector.append(float(term))
		subspace.append(vector)
	basis = []
	basis = gs(subspace)
	print basis
	basis = [[1 ,0 ,0],[0, 1, 0], [0,0,1]]
	
	
	#---------------------------------
	#P matrix is called Proj
	#---------------------------------
	
	Projtemp = []  							#contains vectors in plane first and then orthogonal
	for i in range(0,len(basis)):
		Projtemp.append(basis[(i+1) % dim])
	#print "Projtemp = ", Projtemp
	
	Proj = np.array(Projtemp)
	Projtrans = Proj.T
	#print "Proj = ", Proj
	
	fig = plt.figure()
	ax = fig.add_subplot(111)
	for k in range(0,num):
	
		center = ellips[k]
		#print "center = ", ellips[k]
		
		A = ellarr[k]
		#print "A from ellipse list = ", A

		#------------------------
		#the projected matrix is Proj*A*Projtrans
		#------------------------

		projected = np.dot(np.dot(Proj,A),Projtrans)
		projcenter = np.dot(Projtrans,center)
		print "Projected and projcenter = " , projected, projcenter
		
		#print "center and projcenter = " , center,  projcenter
		#now we have to select the corner nine points from this projected matrix

		center = []
		A = []
		for i in range(0,dim-1):
			row = []
			for j in range(0,dim-1):
				row.append(projected[i][j])
			A.append(row)
		for i in range(0,dim-1):
			center.append(projcenter[i])

		#This is the part same as before

		U, s, rotation = linalg.svd(A)
		radii = 1.0/np.sqrt(s)

		u = np.linspace(0.0, 2.0 * np.pi, 100)
		v = np.linspace(0.0, np.pi, 100)
		x = radii[0] * np.outer(np.cos(u), np.sin(v))
		y = radii[1] * np.outer(np.sin(u), np.sin(v))
		#z = radii[2] * np.outer(np.ones_like(u), np.cos(v))

		for i in range(len(x)):
		    for j in range(len(x)):
		        [x[i,j],y[i,j]] = np.dot([x[i,j],y[i,j]] , rotation) + center

		if k==0:
			ax.plot(x, y, c = 'b' , label='2D projection',linewidth=0.3, linestyle="-")
		else:
			ax.plot(x, y, c = 'r')
	for lin in adj:
		print lin
		n = []; m= [];
		for i in range(0,2):
			n.append(adjval[lin[i]][1])
		for i in range(0,2):

			m.append(adjval[lin[i]][0])
		ax.plot(n,m ,color="green", linewidth=2.0, linestyle="-")	
	

	plt.show()
	plt.close(fig)
	del fig




###This code assumes that the input format is correct and as in the adjoint README file

dim = 0
num = 0
inpfile = open("input3D.txt","r")#inp = open("input.txt","r")
inp = inpfile.readlines()
for index, l in enumerate(inp):
	parts = l.strip().split()
	if index == 0:
		num = int(parts[0])				#num is the number of ellipses being checked
		dim = int(parts[1]) 			#dim is for the dimensions being used
	else:
		if dim == 0:
			print "error in input as dim = 0"
			sys.exit()		
		origin = []						#first detail in line is origin
		for i in range(0,dim):
			origin.append(float(parts[i]))
		ellips.append(origin) 			#store the center values for the nth ellipse
		
		carry = dim;					#carry tracks how many have been read and where to read from next
		
		Umatrix = [] 					#holds the U matrix in 2d list form
		for i in range(0,dim):
			Ulist = []					#List for input
			for j in range(0,dim):
				Ulist.append(float(parts[j*dim + i + carry]))
			Umatrix.append(Ulist)

		Sigmalist = [] 					#holds the sigma matrix
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
		U = np.array(Umatrix,float)
		trans = U.T
		#print "U and trans = " , U, trans		#multiply the three
		mult1 = np.dot(U,eye)
		mult2 = np.dot(mult1,trans)		#A = U * eye * U.T
		#print "A for this iteration = ", mult2
		ellarr.append(mult2)

for x in range(0,len(ellarr)):
	inout.append(0)

travaxis = 0
vector = [0 for x in range(dim)]
maxminus = -primA[travaxis];
critical = None
presslice = 1

method(travaxis, vector, maxminus, critical, presslice)

#print "\n\nadj== ", adj
#tree = main(valcount, adj)
#print tree
#for x in range(0,adj):
print "adj == ",adj;

adjout = open("adjout.txt","w")
for term in adj:
	print >> adjout, term;

# fig = plt.figure(i)

# createbasis(num)

# plt.show()
# plt.close(fig)
# del fig
#plotthis3()

plotthis(num,ellips,ellarr,adj,adjval)