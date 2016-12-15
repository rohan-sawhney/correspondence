import numpy as np
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse.linalg import spsolve

# PyAMG algebraic multigrid
from pyamg import *


def meanCurvatureFlow(mesh, stepSize):
    """
    Convenience function to perform one step of curvature flow. Modifies
    positions directly on the input mesh
    """
    curvatureFlower = cy.MeanCurvatureFlower(mesh, hStep = stepSize)
    curvatureFlower.flowStep()


class MeanCurvatureFlower(object):
    """
    Pronounced 'MeanCurvatureFlow-er', unrelated to plants.

    An object which performs mean curvature flow on a mesh.
    First the object is constructed, then flowStep() is called
    repeatedly to perform 1 or more steps of the flow.

    The parameter hStep gives the size of the step to take.
    """


    def __init__(self, mesh, hStep, updateEachIteration=False):

        self.mesh = mesh
        self.meshInd = self.mesh.enumerateVertices()

        self.hStep = hStep
        self.updateEachIteration = updateEachIteration

        # Store the solver object (= stiffness matrix). Lets us be much faster
        # when updateEachIteration=False
        self.currSolver = None

        # If we're going to be using the same matrix the whole time, lets go ahead
        # and calculate it now so the first call to flowStep() isn't dramatically
        # slower than the rest.
        self.updateStiffnessMatrix()

    def flowStep(self):
        """
        Perform a step of mean curature flow, updating the positions of the mesh vertices
        """

        print("Taking mean curvature flow step of size h="+str(self.hStep))

        if (self.currSolver is None) or self.updateEachIteration:
            self.updateStiffnessMatrix()

        # Solve 3 problems, one each in x,y,z
        nVerts = len(self.mesh.verts)
        for i in range(3):

            #print("  Solving system in " + (['x','y','z'][i]))

            b = np.zeros((nVerts,1))
            for vert in self.mesh.verts:
                b[self.meshInd[vert]] = vert.pos[i]

            sol = self.currSolver.solve(b, tol=1e-8)

            for vert in self.mesh.verts:
                vert.pos[i] = sol[self.meshInd[vert]]

    def updateStiffnessMatrix(self):
        """
        Recomputes the stiffness matrix used in the Poisson problem
        """
        print("Updating mean curvature flow stiffness matrix")

        nVerts = len(self.mesh.verts)

        ## Build the poisson matrix
        A = lil_matrix((nVerts, nVerts))

        # Iterate over the interior vertices, adding terms for the neighbors of each
        for vert in self.mesh.verts:
            i = self.meshInd[vert]
            A[i,i] = 1.0
            for (neighEdge, neighVert) in vert.adjacentEdgeVertexPairs():
                w = neighEdge.cotanLaplace
                j = self.meshInd[neighVert]
                A[i,j] -= self.hStep * w
                A[i,i] += self.hStep * w


        # Solve the system of equations
        A = A.tocsr()
        self.currSolver = ruge_stuben_solver(A)
