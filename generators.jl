include("delta.jl")
include("eisenstein.jl")

#This file contains functions Delta16, Delta18, Delta20, Delta22, and Delta26 of the 
#corresponding 1D spaces of cusp forms (S_k)
#Delta_k is the unique normalised eigenform of weight k and level 1 for k in 
#{12,16,18,20,22,26}


#Return the q expansion for Delta16 to precision prec
function Delta16(prec=10, var="q")

	#cusp form of weight 16 to precision prec
	delta_E4 = delta_qexp(prec, var)*eisenstein_series_poly(4, prec)

	return (-8//bernoulli(4))*truncate(delta_E4, prec)
end


#Return the q expansion for Delta18 to precision prec
function Delta18(prec=10, var="q")

        #cusp form of weight 18 to precision prec
        delta_E6 = delta_qexp(prec, var)*eisenstein_series_poly(6, prec)

        return (-12//bernoulli(6))*truncate(delta_E6, prec)
end


#Return the q expansion for Delta20 to precision prec
function Delta20(prec=10, var="q")

        #cusp form of weight 20 to precision prec
	E4_E4 = (eisenstein_series_poly(4, prec))^2
        delta_E4_E4 = delta_qexp(prec, var)*truncate(E4_E4, prec)

        return ((-8//bernoulli(4))^2)*truncate(delta_E4_E4, prec)
end


#Return the q expansion for Delta22 to precision prec
function Delta22(prec=10, var="q")

        #cusp form of weight 22 to precision prec
	E4_E6 = eisenstein_series_poly(4, prec)*eisenstein_series_poly(6, prec)
        delta_E4_E6 = delta_qexp(prec, var)*truncate(E4_E6, prec)

        return (-8//bernoulli(4))*(-12//bernoulli(6))*truncate(delta_E4_E6, prec)
end


#Return the q expansion for Delta26 to precision prec
function Delta26(prec=10, var="q")
	
	#cusp form of weight 26 to precision prec
	E4_E4 = (eisenstein_series_poly(4, prec))^2
	E4_E4_E6 = truncate(E4_E4, prec)*eisenstein_series_poly(6, prec)
	delta_E4_E4_E6 = delta_qexp(prec, var)*truncate(E4_E4_E6, prec)

	return ((-8//bernoulli(4))^2)*(-12//bernoulli(6))*truncate(delta_E4_E4_E6, prec)
end


