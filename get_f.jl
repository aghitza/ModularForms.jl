#Temporary function to get the right f_j from g_i 

#R, q = PowerSeriesRing(ZZ,prec,var)

function get_f(g, dim)
        for i in 2:dim				#or start with 1?
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


