include("delta-test.jl")
include("eisenstein-test.jl")
include("generators-test.jl")
include("victormiller-test.jl")
include("hecke-test.jl")

function test_all()
	test_delta()
	test_eis_series()
	test_generators()
	test_vm_basis()
	test_hecke()
end
