from rpy2.robjects.packages import importr
# import R's "base" package
base = importr('base')

# import R's utility package
utils = importr('utils')

# select a mirror for R packages
utils.chooseCRANmirror(ind=1) # select the first mirror in the list

# R vector of ints, R vector of floats, R vector of strings, embedded R, global environment
from rpy2.robjects.vectors import IntVector,FloatVector,StrVector
from rpy2.robjects import r,globalenv
