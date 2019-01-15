include("../src/eta_product.jl")
include("../src/eta_product_eigenforms.jl")

function test_eta_product_length()
   print("eta_product.length...")

   @test eta_product([[1,24]]).length == 10
   @test eta_product([[1,24]], 100).length == 100

   println("PASS")
end

function test_eta_product_correctness()
   print("eta_product.correctness...")

   for g in keys(ETA_QUOTIENT_EIGENFORM)
     eta = eta_product(g, 20)
     @test eta - ETA_QUOTIENT_EIGENFORM[g] == 0
   end

   eta1 = eta_product([[1,24]],100)
   @test eta1 - delta_qexp(100) == 0

   println("PASS")
end

function test_eta_product()
   test_eta_product_length()
   test_eta_product_correctness()

   println("")
end
