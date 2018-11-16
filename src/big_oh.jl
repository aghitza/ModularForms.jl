#Construct a PowerSeries element to precision prec from a Polynomial element f, over the 
#given ring (which must be the same as the ring corresponding to f)
function big_oh(f, prec, K)				#add var
	R, q = PowerSeriesRing(K, prec, "q")

	power_series = R(0*q)
	max_power = length(f) - 1 	#highest appearing power (OR must be prec?)
	for i in 1:max_power
		coef = coeff(f,i)
		power_series += R(coef*q^i)
	end
	
	println(typeof(power_series))
	return power_series
end
