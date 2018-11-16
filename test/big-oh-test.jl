include("../src/delta.jl")

function test_big_oh_correctness()
	print("big_oh.correctness...")

	R, q = PolynomialRing(QQ, "q")
	f1 = R(0)
	f2 = q + 3q^2 + 5q^3 - 4q^4
	f3 = 5 - q^2 + 4q^5
	
	g1 = big_oh(f1, 3) 
	S = parent(g1)
	q = gen(S)
	@test g1 == 0 + O(q^3)

	g2 = big_oh(f2, 3) 
	S = parent(g2)
	q = gen(S)
	@test g2 == q + 3q^2 + O(q^3)
	
	g3 = big_oh(f2, 7) 
	S = parent(g3)
	q = gen(S)
	@test g3 == q + 3q^2 + 5q^3 - 4q^4 + O(q^7)
	
	g4 = big_oh(f3, 5) 
	S = parent(g4)
	q = gen(S)
	@test g4 == 5 - q^2 + O(q^5)
	
	g5 = big_oh(f3, 6) 
	S = parent(g5)
	q = gen(S)
	@test g5 == 5 - q^2 + 4q^5 + O(q^6)

	println("PASS")
end

function test_big_oh()
	test_big_oh_correctness()

	println("")
end
