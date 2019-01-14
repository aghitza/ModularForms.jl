include("poly_to_power_series.jl") 
export delta_poly, delta_qexp, delta_k_qexp

#This file contains the functions delta_poly and delta_qexp which return the 
#q-expansion of delta (the weight 12 level 1 cusp form) as a polynomial or 
#powerseries (respectively). Moreover, this file contains the function 
#delta_k_qexp, which returns the q-expansion of the normalised generator for 
#(one of the six) one-dimensional spaces of cusp forms of weight k and level 1. 



#Return the q-expansion of the (normalised) cusp form of weight 12 to precision prec 
#(default: 10) as a polynomial over K (default: ZZ) in the variable var (default: "q"). 
function delta_poly(prec::Int=10, var::String="q", K=ZZ)

   R, q = PolynomialRing(ZZ, var)
   delta = q*eta_qexp(24, prec-1, q) 
   RR, q = PolynomialRing(K, var)

   return RR(delta)
end



#Return the q-expansion of the (normalised) cusp form of weight 12 to precision prec 
#(default: 10) as a power series over K (default: ZZ) in the variable var (default: "q").
function delta_qexp(prec::Int=10, var::String="q", K=ZZ) 
	
   delta = delta_poly(prec, var, K)
   power_series = poly_to_power_series(delta, prec)

   return power_series 
end



#Return the q-expansion of the unique normalised eigenform of weight k and level 1
#for k in [12, 16, 18, 20, 22, 26] as a power series to precision prec (default: 10) 
#in ZZ[[var]] (default: "q"). These eigenforms are the normalised generators for the
#six one-dimensional spaces of cusp forms of level 1. 
function delta_k_qexp(k::Int, prec::Int=10, var::String="q")

   if k == 12
      return delta_qexp(prec, var, ZZ)
   elseif k in [16, 18, 20, 22, 26]
      return delta_qexp(prec, var) * eisenstein_series_qexp(k-12, prec, ZZ, var, "integral") 
   else
      error("k must be one of 12, 16, 18, 20, 22, 26")
   end

end
