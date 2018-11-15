include("eisenstein.jl") #need sigma function

#Return the q-expansion of the normalised weight k Eisenstein series on the modular group 
#to precision prec in the ring K (default is QQ). 
#The normalisation is such that the linear coefficient is 1. 
function eisenstein_series_qexp(k, prec, K=QQ, var="q") 

	#error handling
	if k%2 != 0 || k<2
		error("k must be an even positive integer")
	end

	a0 = - bernoulli(k) // 2k
	
	#the Sage code does this
	#try 
	#	a0fac = K(1//a0.den)	
	#catch				#check if this is correct
	#	if a0.den == 0 		
	#		error("The numerator of -B_k/2k must be invertible in the ring K")
	
	
	#now only works for default values QQ and "q"
	R, q = PolynomialRing(QQ, "q")
	#R, q = PowerSeriesRing(QQ, prec, "q") 	#works as well, prec is random
	qexp = fmpz(0)
	for n in 1:prec-1
		qexp += sigma(k-1,n)*q^n
	end
	
	qexp += a0

	return qexp

end
