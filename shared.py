def init():
    global errorOut, originlist, ellipselist, adjmatrix, adjcoordinates, valcount, num, dim, iterate, primA

    # the adjmatrix is the list of edges that being created
    adjmatrix = []

    # the adjcoordinates gives the dimension coordinates (for plotting for the nth point)
    adjcoordinates = []

    # valcount keeps the number of the value of the point being added to the adjacency tree
    valcount = -1

    # primA stores the terms for the primary axis. The one to traverse along on the first iteration
    primA = []

    iterate = 100  # Number of slices to cover the ellipsoid in nD space

    num = 0
    dim = 0

    errorOut = open("file_output/debuginfo.txt", "w")
