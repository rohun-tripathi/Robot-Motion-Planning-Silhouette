def init():
	global originlist, ellipselist, adjmatrix, adjcoordinates, valcount,num,dim,iterate, primA

	adjmatrix = []		# the adjmatrix is the list of edges that being created
	adjcoordinates = [];		# the adjval gives the dimension coordinates (for plotting for the nth point)
	valcount = -1; 		# valcount keeps the number of the value of the point being added to the adjacency tree
	primA = []			#primA stores the terms for the primary axis. The one to traverse along on the first iteration

	iterate = 100 	#Number of slices to cover the ellipsoid in nD space

	num = 0;
	dim = 0;

	