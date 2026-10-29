using CairoMakie
using ColorSchemes

#------------------------------------------------------------------------------
# First we define the plots as insets
# Then we will use the multiple dispatch to create a function with default axis
# The other options should be automatically followed through;
#------------------------------------------------------------------------------

include("finite_size_eigstat.loaddata.jl")
include("finite_size_eigstat.dependencies.jl")
include("finite_size_eigstat.base.jl")
include("finite_size_eigstat.inset.jl")

#include("finite_size_eigstat_inset.jl")

