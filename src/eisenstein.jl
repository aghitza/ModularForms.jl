#Return the q-expansion up to precision prec of the weight k Eisenstein series as a FLINT 
#Fmpz_poly object, normalised such that coefficients are integers with no common factor 
function eisenstein_series_poly(k, prec)

	R, q = PolynomialRing(QQ, "q") 	#QQ as ring to prevent bug at the end
	qexp = R(0)

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

	#sum of sigma function times q^n
	sigmasum = 0
	for n in 2:prec-1
		sigmasum += sigma(fmpz(n),k-1)*q^n
	end

	qexp += a0 + q + sigmasum 		
	
	return qexp
	
end
