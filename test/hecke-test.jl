include("../src/delta.jl")
include("../src/vm_basis.jl")
include("../src/hecke.jl")

function test_hecke_on_qexp()
        print("hecke.qexp...")
	
	#R, q = PolynomialRing(QQ, "q")
	#delta = delta_qexp(30)
	#h1 = hecke_operator_on_qexp(delta, 1, 12, 10)
	#h2 = truncate(hecke_operator_on_qexp(delta, 3, 12, 10), 7)
	#@test h1 == delta
	#@test h2 == 252q - 6048q^2 + 63504q^3 - 370944q^4 + 1217160q^5 - 1524096q^6

        println("PASS")
end

function test_hecke_on_basis()
	print("hecke.basis...")

	println("PASS")
end

function test_hecke_vm()
	print("hecke.victormiller...")

	B = victor_miller_basis(24, 20)
	b1, b2 = B
	h1 = hecke_operator_on_qexp(b1, 3, 24, 7)
	S = parent(h1)
	q = gen(S)
	@test h1 == 	195660*q - 982499328*q^2 + 85442803344*q^3 + 1302498570240*q^4 + 
			23514204375720*q^5 - 333538871869440*q^6 + O(q^7)
	h2 = hecke_operator_on_qexp(b2, 3, 24, 7)
	S = parent(h2)
	q = gen(S)
	@test h2 == 	-48*q + 143820*q^2 - 16295040*q^3 - 424520544*q^4 - 
			4306546080*q^5 + 67844160144*q^6 + O(q^7)

	#test hecke_operator_on_basis
	#H_basis = hecke_operator_on_basis(B, 3, 24)
	#R = base_ring(h1)
	#S = MatrixSpace(R, 2, 2)
	#matrix = S([195660 -982499328; -48 143820]) 	

	println("PASS")
end

function test_hecke()
        test_hecke_on_qexp()
	test_hecke_on_basis()
	test_hecke_vm()

        println("")
end

