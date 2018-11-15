#Note: in order to use the function one needs to run "using Nemo"


#Return the q-expansion up to precision prec of the weight k Eisenstein series as a FLINT 
#Fmpz_poly object, normalised such that coefficients are integers with no common factor 
function eisenstein_series_poly(k, prec)

        R, q = PolynomialRing(ZZ, "q")
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
	for n in 2:prec
		sigmasum += sigma(n,k-1)*q^n
	end

	#qexp += -a0 + q + sigmasum 		#there is a bug here
	
	return qexp
	
end




#Compute the sum of t-th powers of positive divisors of n (t>=0, n>=1)
function sigma(n, t)
	
	#error 
	if t<0 
		error("t must be a nonnegative integer")
	elseif n<1
		error("n must be a positive integer")
	end

	sum = 0
	for d in 1:n
		#println(sum)
		if n%d == 0
			sum += d^t
		end
	end

	return sum

end
