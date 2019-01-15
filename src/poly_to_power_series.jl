export poly_to_power_series

#Convert the polynomial f to a relative power series with 
precision
#prec (default: 10).
#The power series has the same coefficient ring and variable name as f.
#We assume that the polynomial f is correct up to precision prec.
function poly_to_power_series(f::PolyElem, K, prec::Int=10)

   #initialise
   varname = string(gen(parent(f)))
   R, q = PowerSeriesRing(K, prec, varname)

   if f == 0
      return R(0)
   end

   d = min(degree(f), prec-1)
   c = [coeff(f, i) for i in 0:d]
   power_series = R(c, d+1, prec, 0)

   return power_series
end


function poly_to_power_series(f::PolyElem, prec::Int=10)
   return poly_to_power_series(f, base_ring(f), prec)
end
