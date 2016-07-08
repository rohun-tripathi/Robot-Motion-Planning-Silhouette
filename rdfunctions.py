import shared as SH


def addToRoadmap(vector1, vector2, Distance,
                 String=''):  # String might have details like the function that called the addition to the adjacency matrix
    SH.adjmatrix.append([vector1, vector2, Distance, String])


def addToVertices(vector):
    SH.valcount += 1
    SH.adjcoordinates.append(vector[:])
    return SH.valcount