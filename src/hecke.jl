#This file contains the functions hecke_operator_on_basis and hecke_operator_on_qexp.
#NOTE: the functions now only work for cusp forms.  
include("vm_basis.jl") 


#Compute the matrix of the Hecke operator Tn of weight k relative to the given basis B of 
#q expansions for a space of cusp forms. 
#NOTE: input B needs to be an array, output is a matrix. 
#NOTE: now only works for space of cusp forms, not for space of all modular forms. 
function hecke_operator_on_basis(B, n, k)
	
	#error handling
	if n < 1
		error("n (= $n) must be a positive integer")
	end 
	if isa(B, Array) == false
		error("B (= $B) must be an array")
	end

	#construct dxd matrix
	ring = base_ring(B[1]) 
	d = dim_Sk(k)
	S = MatrixSpace(ring, d, d) 
	matrix = S()


	#NOTE: Sage first normalises B, and works with list first and only transforms it 
	#into a matrix in the end, quicker? 


	#compute Tn(f) to precision d+1 for each element f of B
	for j in 1:d
		f = B[j]
		T_f = hecke_operator_on_qexp(f, n, k, d+1)
		for i in 1:d
			#Tn for the jth element of B corresponds to the jth row 
			matrix[j,i] = coeff(T_f,i)
		end
	end

	return matrix
end



#Compute the image of the q expansion f of a modular form under the Hecke operator Tp of
#weight k, for p prime.
#Return a power series to precision prec (default: 10).
#NOTE: now only works for cusp forms. 
function hecke_operator_prime(f, p, k, prec=10)

	R, q = PowerSeriesRing(base_ring(f), prec, "q")
        T_p = R(0)

	#NOTE: we start with m = 1, since a0=0 for cusp forms, and
	#we end with m = prec-1 since we want ... + O(^prec)
        for m in 1:prec-1
                a_mp = coeff(f, m*p)
                coef = a_mp
                if m%p == 0
                        coef += p^(k-1) * coeff(f, Int(m/p))
                end
        T_p += coef*q^m
        end

        return T_p
end



#Compute the image of the q expansion f of a modular form under the nth Hecke 
#operator Tn of weight k, where n is a prime p to some power. 
#Return a power series to precision prec (default: 10). 
#NOTE: now only works for cusp forms.
function hecke_operator_prime_power(f, p, power, k, prec=10)
	
	T_p = hecke_operator_prime(f, p, k, prec)

	if power == 1
		return T_p
	elseif power == 2
		T_1 = f
		T_p_power = hecke_operator_prime(T_p, p, k, prec) - p^(k-1)*T_1
		return truncate(T_p_power, prec)
	end

	T_p_power_minus_1_over_T_p = hecke_operator_prime_power(T_p, p, power-1, k, prec)
	T_p_power_minus_2 = hecke_operator_prime_power(f, p, power-2, k, prec)

	T_p_power = T_p_power_minus_1_over_T_p - p^(k-1)*T_p_power_minus_2
	return truncate(T_p_power, prec)
end



#Compute the image of the q expansion f of a modular form under the Hecke operator Tn of 
#weight k.
#Return a power series to precision prec (default: 10). 
function hecke_operator_on_qexp(f, n, k, prec=10)	
	
	if isprime(fmpz(n))
		return hecke_operator_prime(f, n, k, prec)
	end

	#construct array with pairs corresponding to the factorisation of n
	n = fmpz(n)
	factorisation = factor(n)
	array = [fact for fact in factorisation]
	T_n = 1
	for i in 1:length(array)
		p = array[i].first
		p = Int(p)
		power = array[i].second 
		T_p_power = hecke_operator_prime_power(f, p, power, k, prec)
		T_n = T_n*T_p_power 
	end

	return T_n
end
