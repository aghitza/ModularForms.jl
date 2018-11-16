#Return the q-expansion up to precision prec (default: 10) of the weight k Eisenstein 
#series as a polynomial over QQ in the variable var (default: "q"), normalised such that 
#the linear coefficient is 1 
function eisenstein_series_poly(k, prec=10, var="q")

	R, q = PolynomialRing(QQ, var) 	#QQ as ring to prevent bug at the end
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
	a0 = - bernoulli(k) // 2k 	#bernoulli is fmpq so need //
	
	qexp += a0
	for n in 1:prec-1
		qexp += sigma(fmpz(n),k-1)*q^n
	end

	return qexp
end
