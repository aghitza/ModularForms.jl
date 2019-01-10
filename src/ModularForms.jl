module ModularForms

using Nemo

export eisenstein_series_poly
export eisenstein_series_qexp
export big_oh
export delta_poly
export prime_range
export victor_miller_basis

include("delta.jl")
include("eis_series.jl")
include("generators.jl")
include("vm_basis.jl")
include("hecke.jl")
include("prime_range.jl")

end # module
