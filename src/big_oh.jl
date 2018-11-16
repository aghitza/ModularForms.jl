#Construct a relative power series to precision prec (default: 10) from a polynomial f 
#over the same ring as f in the variable var (default: "q")
function big_oh(f, prec=10, var="q")
	K = base_ring(f)
	R, q = PowerSeriesRing(K, prec, var)

	power_series = R(0*q)
	max_power = min(length(f), prec)
	max_power = max_power - 1
	for i in 0:max_power
		coef = coeff(f,i)
		power_series += R(coef*q^i)
	end
	
	return power_series
end
