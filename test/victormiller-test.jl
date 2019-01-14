include("../src/vm_basis.jl")

function test_vm_basis_length()
   print("vm_basis.length...")

   @test length(victor_miller_basis(12)) == dim_Sk(12)
   @test length(victor_miller_basis(32)) == dim_Sk(32)

   println("PASS")
end

function test_vm_basis_type()
   print("vm_basis.type...")

   @test typeof(victor_miller_basis(12)[1]) == fmpz_rel_series
   @test typeof(victor_miller_basis(32)[2]) == fmpz_rel_series 

   println("PASS")
end

function test_vm_basis_correctness()
   print("vm_basis.correctness...")

   vm1 = victor_miller_basis(12, 6) 
   S = parent(vm1[1])
   q = gen(S)
   @test vm1[1] - (q-24q^2+252q^3-1472q^4+4830q^5+O(q^6)) == 0

   vm2 = victor_miller_basis(24, 6)
   S = parent(vm2[1])
   q = gen(S)
   @test vm2[1] - (q+195660q^3+12080128q^4+44656110q^5+O(q^6)) == 0 
   @test vm2[2] - (q^2-48q^3+1080q^4-15040q^5+O(q^6)) == 0

   vm3 = victor_miller_basis(32, 6)
   S = parent(vm3[1])
   q = gen(S) 
   @test vm3[1] - (q+50220q^3+87866368q^4+18647219790q^5+O(q^6)) == 0
   @test vm3[2] - (q^2+432q^3+39960q^4-1418560q^5+O(q^6)) == 0

   println("PASS")
end

function test_vm_basis()
   test_vm_basis_length()
   test_vm_basis_type()
   test_vm_basis_correctness()

   println("")
end

