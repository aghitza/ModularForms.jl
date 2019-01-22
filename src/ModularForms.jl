module ModularForms

using Nemo

export eisenstein_series_poly
export eisenstein_series_qexp
export poly_to_power_series
export delta_poly
export prime_range
export victor_miller_basis
export Delta16, Delta18, Delta20, Delta22, Delta26
export eta_quotient

include("delta.jl")
include("poly_to_power_series.jl")
include("eis_series.jl")
include("generators.jl")
include("vm_basis.jl")
include("hecke.jl")
include("prime_range.jl")
include("eta_quotient.jl")

end # module
