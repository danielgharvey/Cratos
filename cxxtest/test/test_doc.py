#
# Import and execute the Python test driver for the user guide examples
#

# Imports
try:
    import pyutilib.th as unittest
    pyutilib_available=True
except:
    pyutilib_available=False
import os
from os.path import dirname, abspath, abspath, basename
import sys

if pyutilib_available:
    currdir = dirname(abspath(__file__))+os.sep
    datadir = os.sep.join([dirname(dirname(abspath(__file__))),'doc','examples'])+os.sep

    os.chdir(datadir)
    sys.path.insert(0, datadir)

    from test_examples import *

# Execute the tests
if __name__ == '__main__':
    unittest.main()
