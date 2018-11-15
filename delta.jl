#This file contains the delta function by using eta_qexp from Nemo
#Need "using Nemo" to make it work


#Return the q-expansion of the (normalised) cusp form of weight 12 (level 1) to 
#precision prec (default=10) with coefficients in the ring ZZ. 
function delta_qexp(prec=10, var="q")
	R, q = PolynomialRing(ZZ, var)
	delta = q*eta_qexp(24, prec-1, q) 
	return delta
end

