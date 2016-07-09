Read this README to gain an idea about the functiosn of different files : 

4
0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 20.0 20.0 20.0 20.0
0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 10.0 10.0 10.0 10.0

The Recursion Points function can be improved. It does not account for the obstacle that starts and ends between two slices.

Lists that are screwing with me:
StateList - [0 for x in range( len(ellipselist) )]
This makes it for ALL the ellipsoids that collide with slice in question

The retrunvec list in the createroad function. This function has three lists of its own.
The first is the normal intersection points of the recursion that have to be conencted to the pastvector when the recursion returns.
The second and third are for the critical points. The first of those are only to be added to the future slices and not this one.
The third is the critical points that are the end of an ellipse and are not to added to future vectors, just pastvectors.

1. auxilary.py : this contains auxilary functios used by the other python files
2. auxilary_backup is just a copy of the above file. Last checked it was the same file on 27th Feb 4:32

3. trial.py 


##############################
INPUT FILE FORMAT
The input files are named as such. The general formulation is:
First line : the first integer is the number of obstacles and the second integer is the number of dimensions
#################

##Inout Explained
# This is a global variable that dictates tha state of the slice wrt to different ellipses

# State 1 : means the slice is intersecting and had intersected in the last slice.
# State 2 : means the slice has just started intersecting
# State 3 : means the slice has just stopped intersecting
# State 0 : means the slice is NOT intersecting and had NOT interseted in the last slice 

The complement.py has this function of name "f" but all it does is return 1 and the system is working at this stage. So don't worry about it.


## This the introduction, a readme of sorts
#ellips is a list of lists. List of the origins of the nth ellipse
#ellarr stores the A matrix for diff ellipses

#primA stores the terms for the primary axis. The one to traverse along on the first iteration

#Presslice is to avoid any kind of recursion in first slice due to change in Critical points, the starting point

#ELLINQUES:
# This variable is used temporarily and under the assumption that only one new ellipse will be encompassed in a slice
# Can use inout array instead if there multiple such cases
