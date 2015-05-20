#!usr/bin/env python
import sys

rank = []
parent = []
tree = []

def union(edge,weight):
	x = edge[0]; 	y = edge[1];
	a = find(x);	b = find(y);
	if a == b:	#check for disjoint
		return weight
	else:
		if rank[a] < rank[b]:
			parent[a] = b
		if rank[a] > rank[b]:
			parent[b] = a
		if rank[a] == rank[b]:
			parent[b] = a
			rank[a] += 1
		tree.append(edge)
		return weight + float(edge[2])

def find(x):
	if parent[x] == x:
		return x;
	parent[x] = find(parent[x]);	#path compression
	return int(parent[x]);

def main(nodenum, edges):

	weight = 0.0
	#print "edges = " , edges
	for counter in range(0,nodenum):#for each node
		parent.append(int(counter))
		rank.append(int(0))

	tempedges = sorted(edges, key = lambda x: int(x[2])) #sorting by weight
	for counter in range(0,len(tempedges)):	#call for each edge
		weight = union(tempedges[counter],weight)	#sending one edge
	print "\n\n\nTotal weight = " , weight
	print "\n\nEdges in MST / Planned path = \n"
	
	for lin in tree:
		print lin[0],lin[1],lin[2]
	return tree