include("../src/delta.jl")

function test_big_oh_correctness()
	print("big_oh.correctness...")

	R, q = PolynomialRing(QQ, "q")
	f1 = R(0)
	f2 = q + 3q^2 + 5q^3 - 4q^4
	f3 = 5 - q^2 + 4q^5
	
	g1 = poly_to_power_series(f1, 3) 
	S = parent(g1)
	q = gen(S)
  h1 = 0 + O(q^3)
	@test g1 - h1 == 0

	g2 = poly_to_power_series(f2, 3) 
	S = parent(g2)
	q = gen(S)
  h2 = q + 3q^2 + O(q^3)
	@test g2 - h2 == 0
	
	g3 = poly_to_power_series(f2, 7) 
	S = parent(g3)
	q = gen(S)
  h3 = q + 3q^2 + 5q^3 - 4q^4 + O(q^7)
	@test g3 - h3 == 0
	
	g4 = poly_to_power_series(f3, 5) 
	S = parent(g4)
	q = gen(S)
  h4 = 5 - q^2 + O(q^5)
	@test g4 - h4 == 0
	
	g5 = poly_to_power_series(f3, 6) 
	S = parent(g5)
	q = gen(S)
  h5 = 5 - q^2 + 4q^5 + O(q^6)
	@test g5 - h5 == 0

	println("PASS")
end

function test_big_oh()
	test_big_oh_correctness()

	println("")
end
