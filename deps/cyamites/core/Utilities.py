import numpy as np
import scipy.sparse as sparse
from math import pi

# Normalizes a numpy vector
# This methods modifies its argument in place, but also returns a reference to that
# array for chaining.
# Works on both single vectors and nx3 arrays of vectors (perfomed in-place).
# If zeroError=False, then this function while silently return a same-sized 0
# for low-norm vectors. If zeroError=True it will throw an exception
#
# TODO I feel like this method just shouldn't exist. It creates many bugs and little
# value.... Maybe make an inPlace=False default arg?
def normalize(vec, zeroError=False):


    # Used for testing zeroError
    eps = 0.00000000001

    # Use separate tests for 1D vs 2D arrays (TODO is there a nicer way to do this?)
    if(len(vec.shape) == 1):

        norm = np.linalg.norm(vec)
        if(norm < 0.0000001):
            if(zeroError):
                raise ArithmeticError("Cannot normalize function with norm near 0")
            else:
                vec[0] = 0
                vec[1] = 0
                vec[2] = 0
                return vec
        vec[0] /= norm
        vec[1] /= norm
        vec[2] /= norm
        return vec

    elif(len(vec.shape) == 2):

        # Compute norms for each vector
        norms = np.sqrt( vec[:,0]**2 + vec[:,1]**2 + vec[:,2]**2 )

        # Check for norm zero, if checking is enabled
        if(zeroError and np.any(norms < 0.00000000001)):
            raise ArithmeticError("Cannot normalize function with norm near 0")

        # Normalize in place
        # oldSettings = np.seterr(invalid='ignore')    # Silence warnings since we check above if the user cares
        vec[:,0] /= norms
        vec[:,1] /= norms
        vec[:,2] /= norms
        # np.seterr(**oldSettings)

    else:
        raise ValueError("I don't know how to normalize a vector array with > 2 dimensions")

    return vec


# Normalizes a numpy vector.
# This method returns a new (normalized) vector
# Works on both single vectors and nx3 arrays of vectors (perfomed in-place).
# If zeroError=False, then this function while silently return a same-sized 0
# for low-norm vectors. If zeroError=True it will throw an exception
def normalized(vec, zeroClip=False):

    # Used for testing zeroError
    eps = 0.00000000001

    # Use separate tests for 1D vs 2D arrays (TODO is there a nicer way to do this?)
    if(len(vec.shape) == 1):

        norm = np.linalg.norm(vec)
        if(zeroClip and norm < eps):
            return np.zeros_like(vec)
        return vec / norm

    elif(len(vec.shape) == 2):

        # Compute norms for each vector
        norms = np.sqrt( vec[:,0]**2 + vec[:,1]**2 + vec[:,2]**2 )

        # Check for norm zero, if checking is enabled
        if(zeroClip and  np.any(norms < 0.00000000001)):
            raise ArithmeticError("Cannot normalize function with norm near 0")

        # Normalize in place
        # oldSettings = np.seterr(invalid='ignore')    # Silence warnings since we check above if the user cares
        vec = vec.copy()
        vec[:,0] /= norms
        vec[:,1] /= norms
        vec[:,2] /= norms
        # np.seterr(**oldSettings)

    else:
        raise ValueError("I don't know how to normalize a vector array with > 2 dimensions")

    return vec


# An alias for np.linal.norm, because typing that is ugly
def norm(vec, *args, **kwargs):
    return np.linalg.norm(vec, *args, **kwargs)


# A quicker cross method when calling on a single vector
def cross(u, v):
    return np.array((
        u[1]*v[2] - u[2]*v[1],
        u[2]*v[0] - u[0]*v[2],
        u[0]*v[1] - u[1]*v[0]
        ))

def dot(u,v):
    return np.dot(u,v)

def clamp(val, lower = -float('inf'), upper = float('inf')):
    if val > upper:
        val = upper
    if val < lower:
        val = lower
    return val

def regAngle(theta):
    """
    Returns the argument mapped in to (-pi,pi]
    """
    while theta > pi: theta = theta - 2*pi
    while theta <= -pi: theta = theta + 2*pi
    return theta

def circlePairs(lst):
    """
    Iterate through a list returning [i],[(i+1)%N] circular pairs, including the
    (last,first) pair
    """
    i = iter(lst)
    first = prev = item = i.next()
    for item in i:
        yield prev, item
        prev = item
    yield item, first

def Vector3D(x,y,z):
    return np.array([float(x),float(y),float(z)])

def printVec3(v,iDig=None):
    if iDig is not None:
        return ("({:."+str(iDig)+"f}, {:."+str(iDig)+"f}, {:."+str(iDig)+"f})").format(v[0], v[1], v[2])
    else:
        return "({:.5f}, {:.5f}, {:.5f})".format(v[0], v[1], v[2])

def checkSymmetric(mat, throwError=True):
    '''
    Verifies that a matrix is symmetric and errors out if not.
    '''

    # Doing this check in a scale-independent way is not trivial. The solution
    # here is to first normalize the matrix such that the mean element magnitude is
    # 1.0, then perform a test against an aboslute epsilon. This is not a perfect solution
    eps = 0.00000001
    errorFound = False

    # Handle a sparse matrix
    if sparse.issparse(mat):
        # Handle a dense matrix

        dokMat = sparse.dok_matrix(mat)  # Convert to a better sparse format for this

        sumVal = 0.0
        nNonzero = 0
        for (i,j) in dokMat.iterkeys():
            sumVal += abs(dokMat[i,j])
            nNonzero += 1

        # Normalize
        meanVal = sumVal / nNonzero
        dokMatN = dokMat / meanVal

        # Check everything
        for (i,j) in dokMat.iterkeys():

                if abs(dokMatN[i,j] - dokMatN[j,i]) > eps:
                    errorFound = True
                    print("\n Matrix is non-symmetric with unequal values on the off-diagonal (" + str(i) + "," + str(j) + ")")
                    print("     values: " + str(dokMat[i,j]) + "  ,  " + str(dokMat[j,i]))


    # Handle a dense matrix
    else:

        # Normalize
        sumVal = np.sum(np.abs(mat))
        nNonzero = np.count_nonzero(mat)
        meanVal = sumVal / nNonzero
        matN = mat / meanVal

        # Check everything
        for i in range(mat.shape[0]):
            for j in range(i,math.shape[1]):

                if abs(matN[i,j] - matN[j,i]) > eps:
                    errorFound = True
                    print("\n Matrix is non-symmetric with unequal values on the off-diagonal (" + str(i) + "," + str(j) + ")")
                    print("     values: " + str(mat[i,j]) + "  ,  " + str(mat[j,i]))


    if errorFound and throwError:
        raise Error("ERROR: checkSymmetric test failed")

def checkHermitian(mat, throwError=True):
    '''
    Verifies that a matrix is Hermitian (conjugate symmetric) and errors out if not.
    '''

    # Doing this check in a scale-independent way is not trivial. The solution
    # here is to first normalize the matrix such that the mean element magnitude is
    # 1.0, then perform a test against an aboslute epsilon. This is not a perfect solution
    eps = 0.00000001
    errorFound = False

    # Handle a sparse matrix
    if sparse.issparse(mat):
        # Handle a dense matrix

        dokMat = sparse.dok_matrix(mat)  # Convert to a better sparse format for this

        sumVal = 0.0
        nNonzero = 0
        for (i,j) in dokMat.iterkeys():
            sumVal += abs(dokMat[i,j])
            nNonzero += 1

        # Normalize
        meanVal = sumVal / nNonzero
        dokMatN = dokMat / meanVal

        # Check everything
        for (i,j) in dokMat.iterkeys():

                if abs(dokMatN[i,j] - dokMatN[j,i].conjugate()) > eps:
                    errorFound = True
                    if i == j:
                        print("\n Matrix is non-hermitian with complex values on the diagonal (" + str(i) + "," + str(j) + ")")
                        print("     values: " + str(dokMat[i,j]) + "  ,  " + str(dokMat[j,i]))
                    else:
                        print("\n Matrix is non-hermitian with non-conjugate values on the off-diagonal (" + str(i) + "," + str(j) + ")")
                        print("     values: " + str(dokMat[i,j]) + "  ,  " + str(dokMat[j,i]))


    # Handle a dense matrix
    else:

        # Normalize
        sumVal = np.sum(np.abs(mat))
        nNonzero = np.count_nonzero(mat)
        meanVal = sumVal / nNonzero
        matN = mat / meanVal

        # Check everything
        for i in range(mat.shape[0]):
            for j in range(i,math.shape[1]):

                if abs(matN[i,j] - matN[j,i].conjugate()) > eps:
                    errorFound = True
                    if i == j:
                        print("\n Matrix is non-hermitian with complex values on the diagonal (" + str(i) + "," + str(j) + ")")
                        print("     values: " + str(mat[i,j]) + "  ,  " + str(mat[j,i]))
                    else:
                        print("\n Matrix is non-hermitian with non-conjugate values on the off-diagonal (" + str(i) + "," + str(j) + ")")
                        print("     values: " + str(mat[i,j]) + "  ,  " + str(mat[j,i]))


    if errorFound and throwError:
        raise Error("ERROR: checkHermitian test failed")
