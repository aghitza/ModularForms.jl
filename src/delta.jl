include("big_oh.jl") 
export delta_poly, delta_qexp

#This file contains the functions delta_poly and delta_qexp which return the q expansion 
#of delta (the weight 12 level 1 cusp form) as a polynomial or powerseries (respectively).
#The algorithm uses the function eta_qexp from Nemo. 



#Return the q-expansion of the (normalised) cusp form of weight 12 to precision prec 
#(default: 10) as a polynomial over K (default: ZZ) in the variable var (default: "q"). 
function delta_poly(prec=10, var="q", K=ZZ)

	R, q = PolynomialRing(K, var)
	delta = q*eta_qexp(24, prec-1, q) 

	return delta
end



#Return the q-expansion of the (normalised) cusp form of weight 12 to precision prec 
#(default: 10) as a power series over K (default: ZZ) in the variable var (default: "q").
function delta_qexp(prec=10, var="q", K=ZZ) 
	
	delta = delta_poly(prec, var, K)
	power_series = big_oh(delta, prec, var)

	return power_series 
end
