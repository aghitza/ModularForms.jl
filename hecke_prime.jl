#This file contains the function hecke_operator_on_qexp but then only working for primes

 
#Compute the image of the q expansion f of a modular form under the Hecke operator Tp of
#weight k, for p prime
function hecke_operator_prime(f, p, k, prec)

        T_p = 0
        for m in 1:prec							#check 
		a_mp = coeff(f, m*p)
                coef = a_mp
		if m%p == 0
			coef += p^(k-1) * coeff(f, Int(m/p))
		end
	T_p += coef*q^m
	end
        
	return T_p

end

