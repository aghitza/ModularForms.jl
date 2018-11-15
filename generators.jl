include("delta.jl")
include("eis_series.jl")

#This file contains functions Delta16, Delta18, Delta20, Delta22, and Delta26 of the 
#corresponding 1D spaces of cusp forms (S_k)
#Delta_k is the unique normalised eigenform of weight k and level 1 for k in 
#{12,16,18,20,22,26}


#Return the q expansion for Delta16 to precision prec
function Delta16(prec=10, var="q")

	#cusp form of weight 16 to precision prec
	delta_times_E4 = delta_qexp(prec, var)*eisenstein_series_qexp(4, prec)

	return (-8/bernoulli(4))*delta_times_E4


#Return the q expansion for Delta18 to precision prec
function Delta18(prec=10, var="q")

        #cusp form of weight 18 to precision prec
        delta_times_E6 = delta_qexp(prec, var)*eisenstein_series_qexp(6, prec)

        return (-12/bernoulli(6))*delta_times_E6


#Return the q expansion for Delta20 to precision prec
function Delta20(prec=10, var="q")

        #cusp form of weight 20 to precision prec
        delta_times_E4_twice = delta_qexp(prec, var)*(eisenstein_series_qexp(4, prec)^2)

        return ((-8/bernoulli(4))^2)*delta_times_E4_twice		#must do //2?
