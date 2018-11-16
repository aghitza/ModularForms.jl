#This file contains the function sigma2. However, Nemo has a built-in version of the sigma 
#function, so we will not use sigma2. 

#Compute the sum of t-th powers of positive divisors of n (t>=0, n>=1)
function sigma2(n, t)

	#error
	if t<0
		error("t must be a nonnegative integer")
	elseif n<1
		error("n must be a positive integer")
	end

	sum = 0
	for d in 1:n
		if n%d == 0
			sum += d^t
		end
	end

	return sum
end
