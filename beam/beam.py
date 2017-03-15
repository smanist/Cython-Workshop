# Standard libraries
import time
import numpy as np
import matplotlib.pyplot as plt

# User libraries
import npUtils
import pyUtils
# import pyUtils_alter as pyUtils
import cyUtils
import pcUtils

class Assembler(object):
    """
    This class provides the functions for
    (1) Assembling the element stiffness matrices and load vectors into their global counterparts
    (2) Solving the resulting linear system using given matrix library
    """
    def __init__(self, side, band, block, utils):
        self.side  = side   # Equals to the total number of DOFs of the matrix
        self.band  = band   # Equals to the number of DOFs per element
        self.block = block  # Equals to the number of DOFs per node
        self.utils = utils  # The library for matrix operations
        self.reset()

    def reset(self):
        # Set the matrix and vectors to zero
        self.K = np.zeros((self.side, self.side))  # The stiffness matrix
        self.F = np.zeros((self.side,))            # The loading vector
        self.u = np.zeros((self.side,))            # The solution vector

    def assembleMat(self, coef, k0, k1):
        # Assemble the element stiffness matrices into the global stiffness matrix
        self.utils.AddMat(self.K, coef, k0, k1, self.block, self.band)

    def assembleVec(self, coef, f0, f1):
        # Assemble the element load vector into the global load vector
        self.utils.AddVec(self.F, coef, f0, f1, self.block, self.band)

    def applyBC(self, bc, penalty):        
        # Apply boundary conditions to the global stiffness matrix
        ref = np.amax(np.diag(self.K))
        for i in bc:
            self.K[i, i] += ref*penalty

    def solve(self):
        # Solve for the linear system
        self.utils.Decomp(self.K, self.band)
        self.utils.Solve(self.K, self.F, self.u, self.band)
        return self.u

class ElemEB(object):
    """
    A two-node element for Euler-Bernoulli beam.
    Each node has following degrees-of-freedoms (DOF's): Deflection (0) and rotation (1).
    The hermite functions are used for the shape functions.
    The bending stiffness and load distribution are assumed linear over the element.
    See Formulation.ipynb for detailed derivation.
    """
    def __init__(self, L=1.0):
        self.eNode = 2  # Number of node of this element
        self.eDof  = 2  # Number of DOF's of this element

        # The reference element
        l1 = 1.0/L
        l2 = l1**2
        l3 = l1**3
        # Stiffness matrix and load vector contributed by LEFT node
        self._K0 = np.array([[  6.0*l3,  4.0*l2,  -6.0*l3,  2.0*l2],
                             [  4.0*l2,  3.0*l1,  -4.0*l2,  1.0*l1],
                             [ -6.0*l3, -4.0*l2,   6.0*l3, -2.0*l2],
                             [  2.0*l2,  1.0*l1,  -2.0*l2,  1.0*l1]])
        self._F0 = np.array([7.0/20.0, L/20.0, 3.0/20.0, -L/30.0])*L
        # Stiffness matrix and load vector contributed by RIGHT node
        self._K1 = np.array([[  6.0*l3,  2.0*l2,  -6.0*l3,  4.0*l2],
                             [  2.0*l2,  1.0*l1,  -2.0*l2,  1.0*l1],
                             [ -6.0*l3, -2.0*l2,   6.0*l3, -4.0*l2],
                             [  4.0*l2,  1.0*l1,  -4.0*l2,  3.0*l1]])
        self._F1 = np.array([3.0/20.0, L/30.0, 7.0/20.0, -L/20.0])*L

        # Boundary types
        self.bc_type = {"cc" : [0, 1],  # Clamped, deflection=rotation=0
                        "ss" : [0],     # Simply-supported, deflection=0
                        "rr" : [1],     # "Roller"-supported, rotation=0
                        "ff" : []}      # Free

class FESolver(object):
    """
    Solves the bending of an Euler-Bernoulli beam due to transverse loading.
    Geometry: Length L with uniformly-distributed N elements, variable thickness
    Formulation: See ElemEB
    """
    def __init__(self, **kwargs):
        self.def_opt = {
            "boundary"  : [],    # Boundary conditions at two ends
            "penalty"   : 1e6,   # Penalty factor for applying boundary conditions
            "length"    : 0.0,   # Length of the beam
            "load"      : None,  # Load distribution function
            "nElem"     : 0,     # Number of elements
            "modulus"   : 0.0,   # Young's modulus of the beam
            "thickness" : [],    # Thickness distribution, either a constant,
                                 # or a list [[position], [value]] (to be interpolated)
            "utils"     : None,  # The library for matrix operations
        }
        self._initParams(**kwargs)
        
    def __call__(self, **kwargs):
        t1 = time.time()
        self._updateParams(**kwargs)
        t2 = time.time()
        self._setupSystem()
        t3 = time.time()
        self._solveSystem()
        t4 = time.time()

        print("Update parameters: {0:5.4f} ms".format((t2-t1)*1e3))
        print("Setup System:      {0:5.4f} ms".format((t3-t2)*1e3))
        print("Solve System:      {0:5.4f} ms".format((t4-t3)*1e3))

    def plotSolution(self, fig=None, sty='b-'):
        # The deflection
        disp = self.u[::2]

        # Assign figure object
        if fig is None:
            fig = plt.figure()
        else:
            plt.figure(fig.number)
        
        # plot solution
        plt.plot(self.coor, disp, sty)
        
        return fig

    def _initParams(self, **kwargs):
        # Initialize the list of parameters
        for key in self.def_opt:
            setattr(self, key, self.def_opt[key])
        for key in kwargs:
            if key in self.def_opt.keys():
                setattr(self, key, kwargs[key])
        
        # Generate mesh
        self.coor = np.linspace(0, self.length, self.nElem+1)

        # Generate element
        self.elem = ElemEB(L = self.length/self.nElem)
        
        # Allocate data array
        self.nNode = self.nElem+1
        self.nDof  = self.elem.eNode*self.nNode
        self._band = self.elem.eNode*self.elem.eDof
        self.u     = np.zeros((self.nDof,))
        
        # Generate the assembler
        self.KF = Assembler(self.nDof, self._band, self.elem.eNode, self.utils)

        # Process boundary conditions
        bc_node = [0, self.nNode-1]
        self._bc = []
        for i, bc in enumerate(bc_node):
            for j in self.elem.bc_type[self.boundary[i]]:
                self._bc.append(self.elem.eDof*bc+j)

    def _updateParams(self, **kwargs):
        # Update parameters if defined
        for key in kwargs:
            if key in self.def_opt:
                setattr(self, key, kwargs[key])

        # Check validity of thickness and loading data
        if len(self.thickness) == 1:
            self.thickness = [[0.0, self.length],
                              [self.thickness[0], self.thickness[0]]]

        # The elastic constants
        self._h  = np.interp(self.coor, self.thickness[0], self.thickness[1])
        self._EI = self.modulus * self._h**3 / 12.0

        # The distributed loading
        self._p  = self.load(self.coor)

    def _setupSystem(self):
        # Assemble the stiffness matrix and load vector
        self.KF.reset()
        self.KF.assembleMat(self._EI, self.elem._K0, self.elem._K1)
        self.KF.assembleVec(self._p,  self.elem._F0, self.elem._F1)

        # Apply boundary conditions
        self.KF.applyBC(self._bc, self.penalty)
    
    def _solveSystem(self):
        self.u = self.KF.solve()

if __name__ == "__main__":
    # Parameters
    E  = 1e8
    L  = 1.0
    h0 = 2e-2
    h1 = 1e-2
    P0 = 1.0
    P1 = 2.0
    I0 = E*h0**3/12

    # Creating the solver object
    opts = {
        "boundary"  : ["cc", "ff"],
        "length"    : L,
        "nElem"     : 1000,
        "modulus"   : E,
        # "utils"     : npUtils,
        # "utils"     : pyUtils,
        # "utils"     : cyUtils,
        "utils"     : pcUtils,
        }
    sol = FESolver(**opts)

    # Constant loading, constant thickness
    f1 = lambda x: P0/I0/24.0 * x**2 * (x**2 - 4.0*L*x + 6.0*L*L)
    sol(load = lambda x: P0*np.ones_like(x),
        thickness = [h0])
    g = sol.plotSolution(sty='b-')
    plt.plot(sol.coor, f1(sol.coor), 'k:', label='Case 1')

    # Linear loading, constant thickness
    p2 =  (P0+2*P1)/12/I0 * L*L
    p3 = -(P0+P1)/12/I0 * L
    p4 =   P0/24/I0
    p5 = -(P0-P1)/120/I0/L
    f2 = lambda x: x*x*(p2 + x*(p3 + x*(p4 + p5*x)))
    sol(load = lambda x: P0+(P1-P0)/L*x,
        thickness = [h0])
    sol.plotSolution(fig=g, sty='g-')
    plt.plot(sol.coor, f2(sol.coor), 'k:', label='Case 2')

    # Linear loading, linear thickness
    f3 = lambda x: ( -365.982 + 546.804*x - 182.907*x**2 + 4.*x**3 + 528*np.log(2-x) - 408*x*np.log(2-x) + 72*x*x*np.log(2-x) ) * 250000/E/(2-x)
    sol(load = lambda x: P0+(P1-P0)/L*x,
        thickness = [[0, L], [h0, h1]])
    sol.plotSolution(fig=g, sty='r-')
    plt.plot(sol.coor, f3(sol.coor), 'k:', label='Case 3')

    plt.grid()
    plt.legend(loc=0)
    plt.xlabel('x')
    plt.ylabel('Deflection')
    
    plt.show()
