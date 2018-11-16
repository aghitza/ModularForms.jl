#Construct a relative power series to precision prec (default: 10) from a polynomial f 
#over the same ring as f in the variable var (default: "q").
#We assume that the polynomial f is correct up to precision prec. 
function big_oh(f, prec=10, var="q")

	#initialise 
	K = base_ring(f)
	R, q = PowerSeriesRing(K, prec, var)
	power_series = R(0)

	if f == 0
		return power_series
	end

	max_power = min(length(f), prec)
	max_power = max_power - 1
	for i in 0:max_power
		coef = coeff(f,i)
		power_series += R(coef*q^i)
	end
	
	return power_series
end
