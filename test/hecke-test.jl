include("../src/delta.jl")
include("../src/hecke.jl")

function test_hecke_on_qexp()
        print("hecke.qexp...")
	
	R, q = PolynomialRing(QQ, "q")
	delta = delta_qexp(30)
	h1 = hecke_operator_on_qexp(delta, 1, 12, 10)
	h2 = truncate(hecke_operator_on_qexp(delta, 3, 12, 10), 7)
	@test h1 == delta
	#@test h2 == 252q - 6048q^2 + 63504q^3 - 370944q^4 + 1217160q^5 - 1524096q^6

        println("PASS")
end

function test_hecke_on_basis()
	print("hecke.basis...")

	println("PASS")
end

function test_hecke()
        test_hecke_on_qexp()
	test_hecke_on_basis()

        println("")
end

