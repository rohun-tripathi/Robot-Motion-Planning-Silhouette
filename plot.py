# noinspection PyUnresolvedReferences
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import numpy.linalg as linalg
import matplotlib.pyplot as plt


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

def createbasis(num):
    # plot
    for turn in range(0, 4):
        # print "Enter four points which describes the vector perpendicular to the 3D hyperplane"
        send = []
        filename = ""
        for n1 in range(0, 4):
            if n1 == turn:
                filename += str(1) + " "
            else:
                filename += str(0) + " "
                send.append(n1)
        print filename
        vector = []
        subspace = []
        for temp in filename.split():
            term = temp.strip()
            vector.append(float(term))

        subspace.append(vector)
        # print "Points you entered : " , vector

        # we have put the first one
        # for random
        for i in range(0, 3):
            vector = []
            for j in range(0, 4):
                term = random.randrange(0, 5 + 1)  # random number between 0 and 6 including them
                vector.append(float(term))
            subspace.append(vector)
        basis = []
        basis = gs(subspace)
        temp_turn = 200 + 20 + turn
        plotthis2(num, basis, temp_turn, send)


def plotthis2(num, basis, turn, send):
    # ---------------------------------
    # P matrix as described by sir in notes is called Proj
    # ---------------------------------
    Projtemp = []  # contains vectors in plane first and then orthogonal
    for i in range(0, len(basis)):
        Projtemp.append(basis[(i + 1) % dim])
        # print "Projtemp = ", Projtemp

    Proj = np.array(Projtemp)
    Projtrans = Proj.T
    # print "Proj = ", Proj

    ax = fig.add_subplot(turn, projection='3d')
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
        for i in range(0, 3):
            row = []
            for j in range(0, 3):
                row.append(projected[i][j])
            A.append(row)
        for i in range(0, 3):
            center.append(projcenter[i])

            # This is the part same as before

        U, s, rotation = linalg.svd(A)
        radii = 1.0 / np.sqrt(s)

        u = np.linspace(0.0, 2.0 * np.pi, 100)
        v = np.linspace(0.0, np.pi, 100)
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

    adjval_temp = []
    for n1 in adjval:
        n2 = np.array(n1)
        adjval_temp.append(np.dot(Proj, n2))

    for lin in adj:
        print lin
        n = [];
        m = [];
        p = [];
        for i in range(0, 2):
            n.append(adjval_temp[lin[i]][send[0]])
        for i in range(0, 2):
            m.append(adjval_temp[lin[i]][send[1]])
        for i in range(0, 2):
            p.append(adjval_temp[lin[i]][send[2]])

        ax.plot_wireframe(n, m, p, color="green", linewidth=2.0, linestyle="-")


def plotthis3():
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
