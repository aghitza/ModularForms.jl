include("../src/delta.jl")

function test_delta_length()
	print("delta.length...")
	
	@test length(delta_qexp(5)) == 5
	@test length(delta_qexp(500)) == 500

	println("PASS")
end

function test_delta()
	test_delta_length()

	println("")
end
