import numpy as np
import numpy.linalg as linalg
import sys, time

import shared as SH
import auxilary as aux

def processinput(originlist, ellipselist, debug=False):
	if debug == True: print "In ProcessInput, Check"
	inp = open("input3D.txt","r").readlines()
	for index, l in enumerate(inp):
		parts = l.strip().split()
		if index == 0:
			SH.num = len(inp) - 1				#Number of ellipses = Number of lines - line for dim
			SH.dim = int(parts[0]) 			#dim is for the dimensions being used
			continue
		else:
			if SH.dim == 0:
				print "error in input as dim = 0"
				sys.exit()		
			origin = []						#first detail in line is origin
			for i in range(SH.dim):
				origin.append(float(parts[i]))
			originlist.append(origin) 			#store the center values for the nth ellipse
			carry = SH.dim;					#carry tracks how many have been read and where to read from next
			
			Umatrix = [] 					#holds the U matrix in 2d list form
			for i in range(SH.dim):
				Ulist = []					#List for input
				for j in range(SH.dim):
					Ulist.append(float(parts[j*SH.dim + i + carry]))
				Umatrix.append(Ulist)

			Sigmalist = [] 					#holds the sigma matrix
			carry = SH.dim + SH.dim*SH.dim 			#Carry gets updated after reads
			for i in range(SH.dim):			#Calculations for A of matrix
				term = float(parts[i + carry])
				inver = 1/(term*term)
				Sigmalist.append(inver)
				if index == 1:
					SH.primA.append(term)
			
			eye = np.identity(SH.dim,dtype = float)
			for i in range(SH.dim):
				eye[i][i] = Sigmalist[i] 		#now eye stores the SIGMA matrix
			U = np.array(Umatrix,float)
			trans = U.T
			#print "U and trans = " , U, trans		#multiply the three
			mult1 = np.dot(U,eye)
			mult2 = np.dot(mult1,trans)		#A = U * eye * U.T
			#print "A for this iteration = ", mult2
			ellipselist.append(mult2)
	return originlist, ellipselist




def f(dim,num,X,ellips,ellarr):
	return 1
	for i in range(0,num):
		arrlist =[]
		temp = []
		for j in range(0,dim):
			arrlist.append(X[j] - ellips[i][j])	#X is point and ellips stores origins
		temp.append(arrlist)
		#print "temp = " , temp
		arrtrans = np.array(temp)
		arr = arrtrans.T
		#print "arr and arrtrans = " , arr, arrtrans
		mult2 = ellarr[i]
		mult3 = np.dot(arrtrans, mult2)
		mult4 = np.dot(mult3,arr)
		#print "mult3 and 4 = ", mult3, mult4
		print "index== " , i, "output== ", mult4[0][0]
		if i==0 and mult4[0][0] > 1:
			print "index== " , i, "Final output = 0"
			return 0
			#sys.exit(0)
		elif i!=0 and mult4[0][0] < 1:
			print "index== " , i, "Final output = 0"
			return 0
			#sys.exit(0)
		
	return 1
