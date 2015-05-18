import numpy as np
import sys
import numpy.linalg as linalg
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import math

def plotthis(num,ellips,ellarr, tree, adjval):		#module from the net
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	for k in range(0,num):#fornow
		center = ellips[k]			#print "center = ", center
		A = ellarr[k] 				#print "a = ", A
	
		U, s, rotation = linalg.svd(A)
		radii = 1.0/np.sqrt(s)

		u = np.linspace(0.0, 2* np.pi, 100)
		v = np.linspace(0.0, 2*np.pi, 100)
		x = radii[0] * np.outer(np.cos(u), np.sin(v))
		y = radii[1] * np.outer(np.sin(u), np.sin(v))
		z = radii[2] * np.outer(np.ones_like(u), np.cos(v))

		for i in range(len(x)):
		    for j in range(len(x)):
		        [x[i,j],y[i,j],z[i,j]] = np.dot([x[i,j],y[i,j],z[i,j]], rotation) + center

		if k==0:		#Different colors
			ax.plot_wireframe(x, y, z,  rstride=4, cstride=4, color='b', alpha=0.2)
		else:
			ax.plot_wireframe(x, y, z,  rstride=4, cstride=4, color='r', alpha=0.2)
		
		# print "x = ", x, len(x[0]) 		# print "y = ", y, len(y) 		# print "Z = ", z, len(z)
	for lin in tree:
		n = []; m= []; p = [];
		for i in range(0,2):
			n.append(adjval[lin[i]][0])
		for i in range(0,2):
			m.append(adjval[lin[i]][1])
		for i in range(0,2):
			p.append(adjval[lin[i]][2])

		ax.plot_wireframe(n,m,p, color="green", linewidth=2.0, linestyle="-")


	plt.show()
	plt.close(fig) 	#print "finished plotthis"
	del fig
