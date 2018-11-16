include("delta.jl")
include("eis_series.jl") 
include("eisenstein.jl") #temporary to use eisenstein_series_poly

#This file contains functions to compute the dimension of spaces of modular forms or cusp 
#forms, and a function to compute the Victor Miller basis for a given weight k to any desired 
#precision
#Using "Modular Forms: A Computational Approach" by William A. Stein


#Algorithm uses Corollary 2.15 and 2.16 from William A. Stein
#Return the dimension of the space of cusp forms of given weight k
function dim_Sk(k)
	
	#case 1: k is odd or smaller than 14
	if k%2 != 0 || k<14
		dim = 0	
	#case 2: k congruent to 2 (mod 12)
	elseif mod(k,12) == mod(2,12)
		dim = floor((k-12)/12)
	#case 3: k not congruent to 2 (mod 12)
	else 
		dim = floor((k-12)/12) + 1
	end

	return Int64(dim) 	#otherwise returns a float
end


#The following function can be used to test dim_Sk
#Algorithm uses Corollary 2.16 from William A. Stein
#Return the dimension of the space of modular forms of given weight k
function dim_Mk(k)

        #case 1: k is odd or negative 
        if k%2 != 0 || k<0
                dim = 0 
        #case 2: k congruent to 2 (mod 12)
        elseif mod(k,12) == mod(2,12)
                dim = floor(k/12)
        #case 3: k not congruent to 2 (mod 12)
        else 
                dim = floor(k/12) + 1
	end

        return Int64(dim) 	#otherwise returns a float
end



#Function to get the right f_j from g_i in order to construct the Victor Miller basis
function get_f(g, dim)
        for i in 2:dim                          #or start with 1?
                #println("i = $i")
                for j in 1:i-1
                        #println("j = $j")
                        #println("a$i(g$j) = $coeff(g[j],i)")
                        g[j] = g[j] - coeff(g[j],i)*g[i]
                        #println("g$j = $g[j]")
                end
        end
        return g
end



#Algorithm uses the proof of Lemma 2.20 from William A. Stein
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
		a = 0
		b = 0
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


	#E4 = eisenstein_series_qexp(k, prec)
        #E6 = eisenstein_series_qexp(k, prec)
        E4 = eisenstein_series_poly(4, prec)  	#temp
	E6 = eisenstein_series_poly(6, prec) 	#temp
	F4 = (-8//bernoulli(4))*E4
        F6 = (-12//bernoulli(6))*E6
	
	#construct a dx1 matrix g where d = dim(Sk)
	d = dim_Sk(k)
	#R, q = PowerSeriesRing(ZZ, prec, "q") 	#not sure about ZZ, prec
	if d == 0
		return []
	end
	g = [truncate((delta_qexp(prec)^j)*(F6^(2(d-j)+b))*(F4^a),prec) for j in 1:d]
	
	#for j in 1:d
		#problem with PowerSeriesRing
		#g[j,1] = (delta_qexp(prec)^j)*(F6^(2(d-j)+b))*(F4^a)
	#end

	#construct f_i from the g_j
	basis_f = get_f(g,d) 

	return basis_f

end
