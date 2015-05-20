#!usr/bin/env python
import numpy as np
import numpy.linalg as linalg
import sys

from auxilary import *
from graph import *
from plot import *
import time



#ELLINQUES:
# This variable is used temporarily and under the assumption that only one new ellipse will be encompassed in a slice
# Can use inout array instead if there multiple such cases

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
