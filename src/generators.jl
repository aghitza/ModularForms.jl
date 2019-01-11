include("delta.jl")
include("eis_series.jl")
include("big_oh.jl")

#This file contains functions Delta16, Delta18, Delta20, Delta22, and Delta26 of the 
#corresponding 1D spaces of cusp forms (S_k). Each function returns a power series.
#Delta_k is the unique normalised eigenform of weight k and level 1 for k in 
#{12,16,18,20,22,26}.


#Return the q expansion for Delta16 to precision prec (default: 10) in var (default: "q").
function Delta16(prec=10, var="q")
  return delta_qexp(prec, var) * eisenstein_series_qexp(4, prec, ZZ, var, "integral")
end


#Return the q expansion for Delta18 to precision prec (default: 10) in var (default: "q").
function Delta18(prec=10, var="q")
  return delta_qexp(prec, var) * eisenstein_series_qexp(6, prec, ZZ, var, "integral")
end


#Return the q expansion for Delta20 to precision prec (default: 10) in var (default: "q").
function Delta20(prec=10, var="q")
  return delta_qexp(prec, var) * eisenstein_series_qexp(8, prec, ZZ, var, "integral")
end


#Return the q expansion for Delta22 to precision prec (default: 10) in var (default: "q"). 
function Delta22(prec=10, var="q")
  return delta_qexp(prec, var) * eisenstein_series_qexp(10, prec, ZZ, var, "integral")
end


#Return the q expansion for Delta26 to precision prec (default: 10) in var (default: "q"). 
function Delta26(prec=10, var="q")
  return delta_qexp(prec, var) * eisenstein_series_qexp(14, prec, ZZ, var, "integral")
end
