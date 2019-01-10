module ModularForms

using Nemo

export eisenstein_series_poly
export big_oh
export delta_poly

include("delta.jl")
include("eis_series.jl")
include("eisenstein.jl")
include("generators.jl")
include("vm_basis.jl")
include("hecke.jl")

end # module
