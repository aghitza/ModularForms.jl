include("../src/eta_product.jl")

function test_eta_product_length()
   print("eta_product.length...")

   @test eta_product([[1,24]]).length == 10
   @test eta_product([[1,24]], 100).length == 100

   println("PASS")
end

function test_eta_product_correctness()
   print("eta_product.correctness...")

   eta1 = eta_product([[1,24]],100)
   @test eta1 - delta_qexp(100) == 0
   eta2 = eta_product([[1,6],[3,6]],20)
   S = parent(eta2)
   q = gen(S)
   @test eta2 - (q-6q^2+9q^3+4q^4+6q^5-54q^6-40q^7+168q^8+81q^9-36q^10-564q^11+36q^12+638q^13+240q^14+54q^15-1136q^16+882q^17-486q^18-556q^19+O(q^20)) == 0
   eta3 = eta_product([[1,2],[2,2],[3,2],[6,2]],20)
   S = parent(eta3)
   q = gen(S)
   @test eta3 - (q-2q^2-3q^3+4q^4+6q^5+6q^6-16q^7-8q^8+9q^9-12q^10+12q^11-12q^12+38q^13+32q^14-18q^15+16q^16-126q^17-18q^18+20q^19+O(q^20)) == 0

   println("PASS")
end

function test_eta_product()
   test_eta_product_length()
   test_eta_product_correctness()

   println("")
end
