#STill don't know presslice

def method(travaxis,vector, startpt, crit, presslice):	#traversal axis
	global valcount,num,dim,adj,adjval;
	
	metreturn = []			#Will be returned as "pastvector" for next slice
	pastvector = []
	
	vector[travaxis] = startpt
	inplace = vector[:]		#inplace used only in next line
	adjval.append(inplace)
	metreturn.append(valcount); #VERY IMP
	pastvector.append(valcount); valcount+=1;		#DO NOT APPEND BEFORE SENDING!!!!!
	cpt = 1 				#Critical point		
	
	count = 0				#loop over variable "iterate"
	while 1: 	
		count += 1;
		if (count==iterate): 
			break;

		vector[travaxis] = startpt + ((startpt) * (-2) * count )/(iterate-1) 
		V = []				#Holds vector for this slice. will be connected to pastvector
		CV = []
		itercpt = 0			#Critical Points in this iteration
		atleastone = 0		#tracks if the plane intersects atleast one ellipse

		if crit == None:
			pass
		elif (vector[travaxis] >= crit[travaxis] and crit[travaxis] > vector[travaxis] - ((startpt)*(-2))/(iterate-1)):
			#print "found ", crit[travaxis], "between ", vector[travaxis], vector[travaxis] - ((startpt)*(-2))/(iterate-1)
			vector[travaxis] = crit[travaxis]
			crit_temp = crit[:];
			CV.append(crit_temp); itercpt += 1;
			#time.sleep(2)
			count -= 1
		
		ellinques = -1	#assuming at one instant, only one new ellipse is added to the network to be encompassed
		axis2 = travaxis;
		for fg in range(0,2):
			axis2 = axis2 + 1
			if axis2 < dim:
				pass
			else:
				break;
			for index,ellipse in enumerate(ellarr):
				origin = ellips[index]
				#print "\norigin and vector = ", origin, vector
				solution, valid = intersect(ellipse, origin, vector, axis2)		#maxima along next axis
				if valid == 1:
					#check that each point is actully valid and not entering another ellipse!
					if inout[index] == 2:		#Comments on inout in complement	
						inout[index] = 1
					if inout[index] == 0:
						inout[index] = 2
						ellinques = index
					atleastone = 1;

					if solution[1] != solution[0]:
						itercpt += 2;					
						vector2 = vector[:]; 			vector2[axis2] = solution[0]
						vector3 = vector[:]; 			vector3[axis2] = solution[1]
						if vector2 != crit:
							fvalid = f(dim,num,vector2,ellips,ellarr)
							if fvalid == 1:
								V.append(vector2); 	
						if vector3 != crit:
							fvalid = f(dim,num,vector3,ellips,ellarr)
							if fvalid == 1:
								V.append(vector3);
					else:
						itercpt += 1;	
						#print "Solution when equal == ", solution, vector
						#time.sleep(2)				
						vector2 = vector[:]; 			vector2[axis2] = solution[0]
						if vector2 != crit:
							CV.append(vector2); 			#NO Doubt about validity		
				else:
					if inout[index] == 1 or inout[index] == 2:	#Comments on inout in complement
						inout[index] = 3
						ellinques = index
					else:
						inout[index] = 0
					pass
					#print "Invalid for count== ", count, " vector== ", vector, "ellipse== ", index
		#print "\nsolution, vector, V, axis2, itercpt, presslice", solution, vector, V, axis2, itercpt, presslice
		#time.sleep(.5)
		
		axis2 =travaxis + 1;
		if itercpt != cpt:
			print "\nCHANGE : solution, vector, V, axis2, itercpt, presslice", solution, vector, V, axis2, itercpt, presslice
			#time.sleep(2)
		
		if count == iterate-1:
			pastvector, valcount , adjval, adj = complete_link(pastvector, V, CV, valcount, adjval, adj)
			metreturn.append(valcount-1);
			#print "\n\nappended to metreturn == ", valcount
	
			#print "\n\nSlicing done as travaxis = ",travaxis
			return metreturn

		if atleastone == 0:
			continue;			#it means it does not need to connect this point
		


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
						#print "solution from cpfind == ", solution, vector, valid
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
						 		tempv1[travaxis] = solution[0]			#Find the one closer to the vector
						 	#else leave it as it is
							max2minus = 0;			#need to calculate maxminus
							solution, valid = intersect(ellarr[0],ellips[0], tempv1, (axis2))		#always check maxima along the next axis
							if(solution[0]  < 0):
								max2minus = solution[0]
							else:
								max2minus = solution[1]		#find the point on the primry to start recursion
							critical = tempv1[:]
							critical[axis2] = ellips[ellinques][axis2]			#crit creation
							
							#print "solution for calling == ", tempv1, max2minus, critical
							#time.sleep(2)

							temvector = method(axis2 , tempv1, max2minus, critical, 1)
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
