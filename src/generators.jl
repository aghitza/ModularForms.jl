include("delta.jl")
include("eis_series.jl")
include("big_oh.jl")

#This file contains functions Delta16, Delta18, Delta20, Delta22, and Delta26 of the 
#corresponding 1D spaces of cusp forms (S_k). Each function returns a power series.
#Delta_k is the unique normalised eigenform of weight k and level 1 for k in 
#{12,16,18,20,22,26}.


#Return the q expansion for Delta16 to precision prec (default: 10) in var (default: "q").
function Delta16(prec=10, var="q")

	#cusp form of weight 16 to precision prec
	delta_E4 = delta_poly(prec, var)*eisenstein_series_poly(4, prec)

	delta16 = (-8//bernoulli(4))*truncate(delta_E4, prec)
	delta16_series = big_oh(delta16, prec, var)

	return delta16_series
end


#Return the q expansion for Delta18 to precision prec (default: 10) in var (default: "q").
function Delta18(prec=10, var="q")

        #cusp form of weight 18 to precision prec
        delta_E6 = delta_poly(prec, var)*eisenstein_series_poly(6, prec)

        delta18 = (-12//bernoulli(6))*truncate(delta_E6, prec)
	delta18_series = big_oh(delta18, prec, var)

	return delta18_series
end


#Return the q expansion for Delta20 to precision prec (default: 10) in var (default: "q").
function Delta20(prec=10, var="q")

        #cusp form of weight 20 to precision prec
	E4_E4 = (eisenstein_series_poly(4, prec))^2
        delta_E4_E4 = delta_poly(prec, var)*truncate(E4_E4, prec)

	delta20 = ((-8//bernoulli(4))^2)*truncate(delta_E4_E4, prec)
	delta20_series = big_oh(delta20, prec, var)

	return delta20_series
end


#Return the q expansion for Delta22 to precision prec (default: 10) in var (default: "q"). 
function Delta22(prec=10, var="q")

        #cusp form of weight 22 to precision prec
	E4_E6 = eisenstein_series_poly(4, prec)*eisenstein_series_poly(6, prec)
        delta_E4_E6 = delta_poly(prec, var)*truncate(E4_E6, prec)

	delta22 = (-8//bernoulli(4))*(-12//bernoulli(6))*truncate(delta_E4_E6, prec)
	delta22_series = big_oh(delta22, prec, var)

	return delta22_series
end


#Return the q expansion for Delta26 to precision prec (default: 10) in var (default: "q"). 
function Delta26(prec=10, var="q")
	
	#cusp form of weight 26 to precision prec
	E4_E4 = (eisenstein_series_poly(4, prec))^2
	E4_E4_E6 = truncate(E4_E4, prec)*eisenstein_series_poly(6, prec)
	delta_E4_E4_E6 = delta_poly(prec, var)*truncate(E4_E4_E6, prec)

	delta26 = ((-8//bernoulli(4))^2)*(-12//bernoulli(6))*truncate(delta_E4_E4_E6, prec)
	delta26_series = big_oh(delta26, prec, var)

	return delta26_series
end
