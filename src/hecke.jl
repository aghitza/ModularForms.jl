#This file contains the functions hecke_operator_on_basis and hecke_operator_on_qexp. 
include("vm_basis.jl")


#Compute the matrix of the Hecke operator Tn of weight k relative to the given basis B 
#of q-expansions for a space of cusp forms. 
#NOTE: input B needs to be an array, output is a matrix. 
function hecke_operator_on_basis(B, n, k)

	#error handling
	if n < 1
		error("n (=$n) must be a positive integer")
	end
	if isa(B, Array) == false
		error("B (=$B) must be an array")
	end

	#construct dxd matrix
	ring = base_ring(B[1])
	d = length(B)			#check if this is always correct
	S = MatrixSpace(ring, d, d)
	matrix = S()

	#compute Tn(f) to precision d+1 for each element f of B
	for j in 1:d
		f = B[j]
		T_f = hecke_operator_on_qexp(f, n, k, d+1)
		#check if cusp form
		if d == dim_Sk(k)
			for i in 1:d
				#Tn for the jth element of B corresponds to the jth row
				matrix[j,i] = coeff(T_f,i)
			end
		else
			for i in 0:d-1
				#Tn for the jth element of B corresponds to the jth row
				matrix[j,i+1] = coeff(T_f,i)
			end		
		end
	end

	return matrix

end


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
	for m in 0:prec-1		#start with 0 to deal with all modular forms
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
