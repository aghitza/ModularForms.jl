#Return the q-expansion up to precision prec (default: 10) of the weight k Eisenstein 
#series as a polynomial over K (default: QQ) in the variable var (default: "q"). 
#Three normalizations are available: "linear" (default), "constant", and 
#"integral". 
#If normalization is "linear" then the linear coefficient is 1. If it is "constant"
#then the series will be normalized to have constant term 1. If the normalization is
#"integral" then the series will be normalized to have integer coefficients and no 
#common factors. 
#NOTE: The output polynomial will be over QQ if the normalization is "linear" or 
#"constant" to prevent errors.  
function eisenstein_series_poly(k, prec=10, K=QQ, var="q", normalization="linear")

	#output will be in QQ if normalization is "linear" or "constant"
	if normalization != "integral"
		K = QQ
	end

	R, q = PolynomialRing(K, var)
	qexp = R(0)

	#error handling 
        if k%2 != 0 || k < 2
                error("k must be an even positive integer")
        end
        if prec < 0
                error("prec must be an even nonnegative integer")
        elseif prec == 0
                return qexp
        end
	
	#initialise a0
	a0 = - bernoulli(k) // 2k 
	
	qexp += a0.num
	for n in 1:prec-1
		qexp += a0.den*sigma(fmpz(n),k-1)*q^n
	end

	if normalization == "linear" 
		return qexp*1//a0.den
	elseif normalization == "constant"
		return qexp*1//a0.num
	elseif normalization == "integral"
		return qexp
	else
		error("normalization must be one of 'linear', 'constant', 'integral'")
	end

end
