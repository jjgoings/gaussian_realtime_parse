from spectra import Spectra 
from realtime import RealTime

# Sample script to print out debug information for l512
# As an example, we use 'test_rabi.log' found in the 'test' subdirectory.

# Read in test_rabi (it knows its the log file) and parse it
rabi = RealTime('./test/test_rabi')

# Print out sanity checks, energy conservation, etc to terminal.
rabi.test()

# For STO-3G H2+ ONLY, we can test Rabi oscillations. It will dump out some 
# relevant information to 'rabi-analysis.csv', in case you want to plot, etc.
rabi.test_rabi()

