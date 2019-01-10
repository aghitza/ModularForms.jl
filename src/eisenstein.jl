#Return the q-expansion up to precision prec (default: 10) of the weight
#k Eisenstein series as a polynomial over ZZ of degree prec-1 in the
#variable var (default: "q"). 
#The polynomial is normalized to have integer coefficients and no 
#common factors. 
#The algorithm is taken from the implementation of eisenstein_series_poly
#in SageMath
function eisenstein_series_poly(k, prec=10, var="q")
  a0 = -bernoulli(k) // 2k
  val = fill(ZZ(a0.den), prec)
  val[1] = ZZ(a0.num)
  expt = k - 1
  for p in prime_range(prec-1)
    ppow = p
    mult = ZZ(p)^expt
    term = mult * mult
    last = mult
    while ppow < prec
      ind = ppow
      term_m1 = term - 1
      last_m1 = last - 1
      while ind < prec
        val[ind+1] *= term_m1
        val[ind+1], r = fdivrem(val[ind+1], last_m1)
        ind += ppow
      end
      ppow *= p
      last = term
      term *= mult
    end
  end

  R, q = PolynomialRing(ZZ, var)
  return R(val)
end
