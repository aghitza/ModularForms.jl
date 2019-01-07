#This file contains the function hecke_operator_on_qexp. 
#NOTE: the function now only works for cusp forms.


#Return an array consisting of all divisors of n if n is nonzero, 
#raise an error otherwise
function divisors(n)

	#error handling
	if n == 0
		error("n must be nonzero")
	end

	n = abs(n)
	array = [div for div in 1:n if n%div == 0]
	return array

end


#Compute the image of the q-expansion f of a modular form under the Hecke operator Tn 
#of weight k. Return a power series to precision prec. 
function hecke_operator_on_qexp(f, n, k, prec=nothing)

	max_prec = Int(ceil(f.prec / Int(n)))	#check if works also for n nonprime 
	if prec == nothing
		prec = max_prec
	elseif prec > max_prec
		error("desired precision is too high given precision of f")
	end

	R, q = PowerSeriesRing(base_ring(f), prec, "q")
	T_p = R(0)	

	l = k-1
	for m in 1:prec-1
		sum = 0
		array = divisors(gcd(n,m)) 	#check if there is a quicker method
		for i in 1:length(array)
			d = array[i]
			if (m*n)%(d*d) == 0
				index = (m*n)//(d*d)
				sum += d^l * coeff(f, Int(index))
			end
		end
		T_p += sum*q^m
	end
	
	return T_p

end
