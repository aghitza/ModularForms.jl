#Return the q-expansion of the normalised weight k Eisenstein series on the modular group 
#to precision prec (default: 10) over QQ in the variable var (default: "q"). 
#The normalisation is such that the linear coefficient is 1. 
function eisenstein_series_qexp(k, prec=10, var="q") 

	#initialise
	R, q = PowerSeriesRing(QQ, prec, var)   	#need to change prec?
	qexp = R(0)

	#error handling
	if k%2 != 0 || k<2
		error("k must be an even positive integer")
	end
	if prec < 0
		error("prec must be an even nonnegative integer")
	elseif prec == 0
		return qexp
	end

	#initialise a0
	a0 = - bernoulli(k) // 2k	
	
	qexp += a0
	for n in 1:prec-1
		qexp += sigma(fmpz(n),k-1)*q^n
	end

	return qexp
end
