#This file contains the functions hecke_operator_on_basis and hecke_operator_on_qexp 
include("hecke_prime.jl")
include("vm_basis.jl") 


#Compute the matrix of the Hecke operator Tn of weight k relative to the given basis B of 
#q expansions for a space of cusp forms 
#Note: input B needs to be an array, output is a matrix 
#Note: now only works for space of cusp forms, not for space of all modular forms 
function hecke_operator_on_basis(B, n, k)
	
	#error handling
	if n < 1
		error("n (= $n) must be a positive integer")
	end 
	if isa(B, Array) == false
		error("B (= $B) must be an array")
	end


	#construct dxd matrix
	ring = parent(B[1]) 		#retrieve the relevant ring (or field)
	d = dim_Sk(k)
	S = MatrixSpace(ring, d, d) 
	matrix = S()


	#NOTE: Sage first normalises B, and works with list first and only transforms it 
	#into a matrix in the end, quicker? 


	#compute Tn(f) to precision d for each element f of B
	for j in 1:d
		f = B[j]
		T_f = hecke_operator_on_qexp(f, n, k, d)
		for i in 1:d
			#Tn for the jth element of B corresponds to the jth row 
			matrix[j,i] = coeff(T_f,i)
		end
	end

	return matrix
end



#Compute the image of the q expansion f of a modular form under the Hecke operator Tn of 
#weight k
function hecke_operator_on_qexp(f, n, k, prec)	
	
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
		power = array[i].second 
		T_n = T_n*(hecke_operator_prime(f, p, k, prec)^power)
	end

	return T_n
end





#Example: 
#include("vm_basis.jl") 
#B = victor_miller_basis(24,20)
#b1, b2 = B
#hecke_operator_on_qexp(b1, 3, 24, 20)
#hecke_operator_on_qexp(b2, 3, 24, 20)
#hecke_operator_on_basis(B, 3, 24)
