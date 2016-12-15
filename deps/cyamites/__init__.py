# Imports from core
from core.InputOutput import *
from core.TriSoupMesh import TriSoupMesh
from core.HalfEdgeMesh import *
from core.Solvers import *

# Imports from viewer
#from viewer.MeshDisplayQT import MeshDisplayQT
from viewer.MeshDisplay import MeshDisplay


# Import from Geometry Processing
from geometry_processing.MeanCurvatureFlow import *


#. Register a signal handler for crtl-c
import signal, sys
def signal_handler(signal, frame):
        print('\nExiting due to SIGINT. Have a nice day.')
        sys.exit(0)
signal.signal(signal.SIGINT, signal_handler)

