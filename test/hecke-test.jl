include("../src/hecke.jl")

function test_hecke_on_qexp()
        print("hecke.qexp...")

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

