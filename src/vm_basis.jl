include("delta.jl")
include("eis_series.jl")
include("big_oh.jl")

#This file contains functions to compute the dimension of spaces of modular forms or cusp 
#forms, and a function to compute the Victor Miller basis for a given weight k to any desired 
#precision (both in terms of polynomials and in terms of power series).
#For the algorithms we make use of "Modular Forms: A Computational Approach" by William A. 
#Stein. 



#Return the dimension of the space of cusp forms of given weight k.
#Algorithm uses Corollary 2.15 and 2.16 from William A. Stein. 
function dim_Sk(k)
	
   d = dim_Mk(k)

   if d == 0
      return 0
   end

   return d - 1 
end



#Return the dimension of the space of modular forms of given weight k.
#Algorithm uses Corollary 2.16 from William A. Stein. 
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
                for j in 1:i-1
                        g[j] = g[j] - coeff(g[j],i)*g[i]
                end
        end

        return g
end



#Return the Victor Miller basis of the space of cusp forms of weight k and level 1 
#(i.e. S_k) to precision prec (default: 10) as an array whose entries are polynomials 
#in ZZ[[var]].
#Algorithm uses the proof of Lemma 2.20 from William A. Stein. 
#NOTE: now returns polynomials over QQ instead of ZZ. 
function victor_miller_basis_poly(k, prec=10, var="q")

	#error handling
	if prec <= 0 
		error("prec must be a positive integer")
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
	#cases (a=3, b=0) and (a=0, b=2) are excluded
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

	F4 = eisenstein_series_poly(4, prec)
	F6 = eisenstein_series_poly(6, prec)
	
	#construct a dx1 matrix g where d = dim(Sk)
	d = dim_Sk(k)
	if d == 0
		return []
	end

	g = [truncate((delta_poly(prec)^j)*(F6^(2(d-j)+b))*(F4^a),prec) for j in 1:d]

	#construct f_i from the g_j
	basis_f = get_f(g,d) 

	return basis_f
end



#Return the Victor Miller basis of the space of cusp forms of weight k and level 1 
#(i.e. S_k) to precision prec (default: 10) as an array whose entries are power series 
#in ZZ[[var]].
#Algorithm uses the proof of Lemma 2.20 from William A. Stein. 
#NOTE: now returns a power series over QQ instead of ZZ.
function victor_miller_basis(k, prec=10, var="q")
	
	vm_basis = victor_miller_basis_poly(k, prec, var)
	power_series = [poly_to_power_series(poly, prec) for poly in vm_basis]

	return power_series
end
