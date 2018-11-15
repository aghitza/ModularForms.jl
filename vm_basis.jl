include("eis_series.jl")

#This file contains a function to compute the Victor Miller basis for a given weight k to 
#any desired precision
#Using proof in book of Stein


#Return the Victor Miller basis of the space of cusp forms of weight k and level 1 
#(i.e. S_k) to precision prec as an array whose entries are power series in ZZ[[var]]
function victor_miller_basis(k, prec, var = "q")

	#error handling
	if prec < 0 
		error("prec must be an even nonnegative integer")
	elseif prec == 0
		return []
	end


	#simple case
	if k == 0
                return [1]
        elseif k%2 != 0 || k<2
                return [] 
        end

	#requirements for integers a,b>=0
	e = mod(k,12) 		#e = 4a + 6b
	if e == 0
		a,b = 0
	#note that a<=3, b<=2 
	#cases (a=3 b=0) and (a=0 b=2) are excluded
	elseif e == 2	#case a=2, b=1
		a = 2
		b = 1
	elseif e == 4 	#case a=1, b=0
		a = 1
		b = 0
	elseif e == 6	#case a=0, b=1
		a = 0
		b = 1
	elseif e == 8 	#case a=2, b=0
		a = 2
		b = 0
	elseif e == 10	#case a=1, b=1 
		a = 1
		b = 1
	end


	E4 = eisenstein_series_qexp(k, prec)
        E6 = eisenstein_series_qexp(k, prec)
        F4 = (-8//bernoulli(4))*E4
        F6 = (-12//bernoulli(6))*E6

	
	#want a (d x prec) matrix where d = dim(Sk)
	#missing part

end
