export eta_product

#Return the eta-product relative to the given collection of integers g as a 
#power series up to precision prec. The output is a cusp form of integral weight.
#The input g is an array of pairs [[t_1,r_1], [t_2,r_2], ..., [t_s,r_s]] where 
#each t_j is a positive integer and each r_j a nonnegative integer. 
#An error is thrown if the sum of t_j*r_j for j in (1,...,s) does not equal 24, 
#since the function only applies to eta-products that are cusp forms.  
function eta_product(g, prec=10)

   R, q = PolynomialRing(ZZ, "q")
   poly = R(1)		#product
   sum = 0 		#error handling
   s = length(g)
   for i in 1:s
      t = g[i][1]
      r = g[i][2]
      poly = poly * eta_qexp(1, prec, q^t)^r	#check prec
      sum += t*r 
   end

   if sum != 24
      error("sum of t_j*r_j for all j must equal 24 if it is a cusp form")
   end

   poly = q * poly
   power_series = poly_to_power_series(poly, prec)

   return power_series

end
