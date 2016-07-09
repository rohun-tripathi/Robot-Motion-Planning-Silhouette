# noinspection PyUnresolvedReferences
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import numpy.linalg as linalg
import matplotlib.pyplot as plt
import random
import shared as SH


def plotthis(num, ellips, ellarr, tree, adjval):  # module from the net
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # This part plots the ellipses
    for k in range(0, num):
        center = ellips[k]  # print "center = ", center
        A = ellarr[k]  # print "a = ", A

        U, s, rotation = linalg.svd(A)
        radii = 1.0 / np.sqrt(s)

        u = np.linspace(0.0, 2 * np.pi, 100)
        v = np.linspace(0.0, 2 * np.pi, 100)
        x = radii[0] * np.outer(np.cos(u), np.sin(v))
        y = radii[1] * np.outer(np.sin(u), np.sin(v))
        z = radii[2] * np.outer(np.ones_like(u), np.cos(v))

        for i in range(len(x)):
            for j in range(len(x)):
                [x[i, j], y[i, j], z[i, j]] = np.dot([x[i, j], y[i, j], z[i, j]], rotation) + center

        if k == 0:  # Different colors
            ax.plot_wireframe(x, y, z, rstride=4, cstride=4, color='b', alpha=0.2)
        else:
            ax.plot_wireframe(x, y, z, rstride=4, cstride=4, color='r', alpha=0.2)

            # print "x = ", x, len(x[0]) 		# print "y = ", y, len(y) 		# print "Z = ", z, len(z)

    # This part plots the tree in the ellipses
    for index, lin in enumerate(tree):
        n = [];
        m = [];
        p = [];
        for i in range(0, 2):
            n.append(adjval[lin[i]][0])
        for i in range(0, 2):
            m.append(adjval[lin[i]][1])
        for i in range(0, 2):
            p.append(adjval[lin[i]][2])
        # a = [n,m,p]
        # pprint.pprint(a)
        ax.plot(n, m, p, color="green", linewidth=2.0, linestyle="-")

    plt.show()
    plt.close(fig)  # print "finished plotthis"
    del fig


# ---------------------------------------------------------
# Code region to calculate the Gram Smidt Orthogonalization
# ----------------------------------------------------------

def modulus(v1):
    mod = float(np.sqrt(float(np.dot(v1, v1))))
    return [x / mod for x in v1]


def gs_cofficient(v1, v2):  # Gram Smidt Coefficient
    return np.dot(v2, v1) / np.dot(v1, v1)


def multiply(cofficient, v):
    return [x * cofficient for x in v]


def proj(v1, v2):
    return multiply(gs_cofficient(v1, v2), v1)


def gs(X):
    mod_vec = 0.0
    Y = []
    for i in range(len(X)):
        temp_vec = X[i]
        for inY in Y:
            proj_vec = proj(inY, X[i])
            temp_vec = map(lambda x, y: x - y, temp_vec, proj_vec)
        mod_vec = modulus(temp_vec)
        Y.append(mod_vec)
        # print "The final matrix returned from gs Y = ", Y
    return Y


# --------------------------------
# this function plots in 3D
# --------------------------------

def getBasisAndSendVector(positionOnPlot):
    send = []
    filename = ""
    for n1 in range(0, 4):
        if n1 == positionOnPlot:
            filename += str(1) + " "
        else:
            filename += str(0) + " "
            send.append(n1)
    vector = []
    subspace = []
    for temp in filename.split():
        term = temp.strip()
        vector.append(float(term))
    subspace.append(vector)
    # we have put the first one
    # for random
    for i in range(0, 3):
        vector = []
        for j in range(0, 4):
            term = random.randrange(0, 5 + 1)  # random number between 0 and 6 including them
            vector.append(float(term))
        subspace.append(vector)
    basis = gs(subspace)
    return send, basis


def getBasisFromPosition(positionOnPlot):
    basisForIteration = []
    baseBasisVector = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]
    for index, term in enumerate(baseBasisVector):
        if index == positionOnPlot: continue
        basisForIteration.append(term[:])
    return basisForIteration


def project4DTo3DAndDisplay(originList, ellMatrixList):
    fig = plt.figure()
    for positionOnPlot in range(0, 4):
        basis = getBasisFromPosition(positionOnPlot)
        fig = plot4DUsingBasis(basis, positionOnPlot, originList, ellMatrixList, fig)
    plt.show()
    plt.close(fig)
    del fig


def plot4DUsingBasis(basis, positionOnPlot, originList, ellMatrixList, fig):
    plotConfigPlusPosition = 200 + 20 + positionOnPlot
    ax = fig.add_subplot(plotConfigPlusPosition, projection='3d')

    # ---------------------------------
    # P matrix as described by sir in notes is called projVector
    # ---------------------------------
    projVector = np.array(basis)
    projTranspose = projVector.T

    renderEllipsoids(ax, ellMatrixList, originList, projTranspose, projVector)
    renderGraph(ax, projVector)

    return fig


def renderGraph(ax, projVector):
    adjvalProjected = []
    for n1 in SH.adjcoordinates:
        n2 = np.array(n1)
        adjvalProjected.append(np.dot(projVector, n2))  # Todo check that the dimension reduces.
    for lin in SH.adjmatrix:
        print lin
        n = [];
        m = [];
        p = [];
        for i in range(0, 2):
            n.append(adjvalProjected[lin[i]][0])
        for i in range(0, 2):
            m.append(adjvalProjected[lin[i]][1])
        for i in range(0, 2):
            p.append(adjvalProjected[lin[i]][2])

        ax.plot_wireframe(n, m, p, color="green", linewidth=1.5, linestyle="-")


def renderEllipsoids(ax, ellMatrixList, originList, projTranspose, projVector):
    for k in range(0, len(originList)):

        center = originList[k]
        ellipse = ellMatrixList[k]

        # ------------------------
        # the projectedMatrix matrix is projVector*A*projTranspose
        # ------------------------

        A = np.dot(np.dot(projVector, ellipse), projTranspose)
        center = np.dot(projVector, center)

        # This is the part same as before
        U, s, rotation = linalg.svd(A)
        radii = 1.0 / np.sqrt(s)

        u = np.linspace(0.0, 2 * np.pi, 100)
        v = np.linspace(0.0, 2 * np.pi, 100)
        x = radii[0] * np.outer(np.cos(u), np.sin(v))
        y = radii[1] * np.outer(np.sin(u), np.sin(v))
        z = radii[2] * np.outer(np.ones_like(u), np.cos(v))

        for i in range(len(x)):
            for j in range(len(x)):
                [x[i, j], y[i, j], z[i, j]] = np.dot([x[i, j], y[i, j], z[i, j]], rotation) + center

        if k == 0:
            ax.plot_wireframe(x, y, z, rstride=4, cstride=4, color='b', alpha=0.2)
        else:
            ax.plot_wireframe(x, y, z, rstride=4, cstride=4, color='r', alpha=0.2)


def plot2DProjection():
    # plot
    print "Enter three points which describes the vector perpendicular to the 3D hyperplane"
    filename = str("1") + " " + str("1") + " " + str("1") + " "
    vector = []
    subspace = []
    for temp in filename.split():
        term = temp.strip()
        vector.append(float(term))

    subspace.append(vector)
    # print "Points you entered : " , vector

    # we have put the first one
    # for random
    for i in range(0, dim - 1):
        vector = []
        for j in range(0, dim):
            term = random.randrange(0, 5 + 1)  # random number between 0 and 6 including them
            vector.append(float(term))
        subspace.append(vector)
    basis = []
    basis = gs(subspace)
    print basis
    basis = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

    # ---------------------------------
    # P matrix is called Proj
    # ---------------------------------

    Projtemp = []  # contains vectors in plane first and then orthogonal
    for i in range(0, len(basis)):
        Projtemp.append(basis[(i + 1) % dim])
        # print "Projtemp = ", Projtemp

    Proj = np.array(Projtemp)
    Projtrans = Proj.T
    # print "Proj = ", Proj

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for k in range(0, num):

        center = ellips[k]
        # print "center = ", ellips[k]

        A = ellarr[k]
        # print "A from ellipse list = ", A

        # ------------------------
        # the projected matrix is Proj*A*Projtrans
        # ------------------------

        projected = np.dot(np.dot(Proj, A), Projtrans)
        projcenter = np.dot(Projtrans, center)
        print "Projected and projcenter = ", projected, projcenter

        # print "center and projcenter = " , center,  projcenter
        # now we have to select the corner nine points from this projected matrix

        center = []
        A = []
        for i in range(0, dim - 1):
            row = []
            for j in range(0, dim - 1):
                row.append(projected[i][j])
            A.append(row)
        for i in range(0, dim - 1):
            center.append(projcenter[i])

            # This is the part same as before

        U, s, rotation = linalg.svd(A)
        radii = 1.0 / np.sqrt(s)

        u = np.linspace(0.0, 2.0 * np.pi, 100)
        v = np.linspace(0.0, np.pi, 100)
        x = radii[0] * np.outer(np.cos(u), np.sin(v))
        y = radii[1] * np.outer(np.sin(u), np.sin(v))
        # z = radii[2] * np.outer(np.ones_like(u), np.cos(v))

        for i in range(len(x)):
            for j in range(len(x)):
                [x[i, j], y[i, j]] = np.dot([x[i, j], y[i, j]], rotation) + center

        if k == 0:
            ax.plot(x, y, c='b', label='2D projection', linewidth=0.3, linestyle="-")
        else:
            ax.plot(x, y, c='r')
    for lin in adj:
        print lin
        n = [];
        m = [];
        for i in range(0, 2):
            n.append(adjval[lin[i]][1])
        for i in range(0, 2):
            m.append(adjval[lin[i]][0])
        ax.plot(n, m, color="green", linewidth=2.0, linestyle="-")

    plt.show()
    plt.close(fig)
    del fig
