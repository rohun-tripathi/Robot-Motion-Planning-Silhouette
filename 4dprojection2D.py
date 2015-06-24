#!usr/bin/env python

#modifying this one to read for and show projection along on
import numpy as np
import sys
import numpy.linalg as linalg
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random

ellips = []		#stores the origin values for diff ellipses
ellarr = []		#stores the A matrix for diff ellipses

dim_reduce = 1; #this will store the value of the number of dimensions being reduced

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
def plotthis(num):
	
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
		for i in range(0,dim-dim_reduce):
			row = []
			for j in range(0,dim-dim_reduce):
				row.append(projected[i][j])
			A.append(row)
		for i in range(0,dim-dim_reduce):
			center.append(projcenter[i])

		#This is the part same as before

		U, s, rotation = linalg.svd(A)
		radii = 1.0/np.sqrt(s)

		u = np.linspace(0.0, 2.0 * np.pi, 50)
		v = np.linspace(0.0, np.pi, 50)
		x = radii[0] * np.outer(np.cos(u), np.sin(v))
		y = radii[1] * np.outer(np.sin(u), np.sin(v))
		#z = radii[2] * np.outer(np.ones_like(u), np.cos(v))

		for i in range(len(x)):
		    for j in range(len(x)):
		        [x[i,j],y[i,j]] = np.dot([x[i,j],y[i,j]] , rotation) + center

		if k==0:
			ax.plot(x, y, c = 'b' , label='2D projection')
		else:
			ax.plot(x, y, c = 'r')

	

	plt.show()
	plt.close(fig)
	del fig


#-------------------------------------
#Checks if the point provided by the user is in the valid region
#-------------------------------------

def f(dim,num,X):
	for i in range(0,num):
		#print "i = " , i
		
		arrlist =[]
		temp = []
		for j in range(0,dim):
			arrlist.append(X[j] - ellips[i][j])	#X is point and ellips stores origins
		temp.append(arrlist)
		#print "temp = " , temp
		arrtrans = np.array(temp)
		arr = arrtrans.T
		#print "arr and arrtrans = " , arr, arrtrans

		#final mult
		mult2 = ellarr[i]

		mult3 = np.dot(arrtrans, mult2)
		mult4 = np.dot(mult3,arr)

		#print "mult3 and 4 = ", mult3, mult4

		if i==0 and mult4[0][0] > 1:
			print "\n\n Final output = 0\n"
			sys.exit(0)
		elif i!=0 and mult4[0][0] < 1:
			print "\n\nFinal output = 0\n"
			sys.exit(0)
	print "\n\nFinal output = 1\n"


#to create a 2D matrix with just one row like below to use for column
#a = np.array([5,4])[np.newaxis]
#or use:
#a = np.array([[5, 4]])
dim = 0
num = 0

#inp = open("input.txt","r")
inp = open("input4D.txt","r")
	
for l in inp.readlines():
	line = l.strip()
	parts = line.split()
	if len(parts) == 2:
		#dim is for the dimensions being used
		#num is the number of ellipses being checked
		num = int(parts[0])
		dim = int(parts[1])
		
		#print "dim = " , dim
	else:
		if dim == 0:
			print "error in input as dim =0"
		
		#first detail in line is origin
		origin = []
		for i in range(0,dim):
			origin.append(float(parts[i]))
		ellips.append(origin) 

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
		eye = np.identity(dim,dtype = float)
		for i in range(0,dim):
			eye[i][i] = Sigmalist[i] 
		
		#print "Sig = " , eye

		#getting the U and the U(Trans)

		U = np.array(Siglist,float)
		trans = U.T

		#print "U and trans = " , U, trans

		#multiply the three
		#----------------------------
		#creation of A happens here
		#----------------------------

		mult1 = np.dot(U,eye) 			#here eye is NOT an identity matrix but infact the sigma matrix
		mult2 = np.dot(mult1,trans)		#A

		#print "A for this iteration = ", mult2
		ellarr.append(mult2)


plotthis(num)

#commenting this for now
#X = [2.8,7,8]	#this is the list of containing point info
#f(dim,num,X)


