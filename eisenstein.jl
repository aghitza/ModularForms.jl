#Return the q-expansion up to precision prec of the weight k Eisenstein series as a FLINT 
#Fmpz_poly object, normalised such that coefficients are integers with no common factor 

function eisenstein_series_poly(k, prec)

        qexp = 0

        if k%2 != 0 || k < 2
                error("k must be an even positive integer")
        end

        if prec < 0
                error("prec must be an even nonnegative integer")
        elseif prec ==0
                return qexp
        end

end


