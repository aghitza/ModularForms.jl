var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "ModularForms.jl Documentation",
    "title": "ModularForms.jl Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "#ModularForms.jl-Documentation-1",
    "page": "ModularForms.jl Documentation",
    "title": "ModularForms.jl Documentation",
    "category": "section",
    "text": ""
},

{
    "location": "#ModularForms.delta_qexp",
    "page": "ModularForms.jl Documentation",
    "title": "ModularForms.delta_qexp",
    "category": "function",
    "text": "delta_qexp(prec=10, var=\"q\", K=ZZ)\n\nReturn the q-expansion of the normalized cusp form of weight 12 to  precision prec as a power series over K in the variable var. \n\nArguments\n\nprec::Integer=10: precision of the output \nvar::String=\"q\": variable name\nK=ZZ: base ring of the output\n\nExamples still missing\n\n\n\n\n\n"
},

{
    "location": "#ModularForms.delta_k_qexp",
    "page": "ModularForms.jl Documentation",
    "title": "ModularForms.delta_k_qexp",
    "category": "function",
    "text": "delta_k_qexp(k, prec=10, var=\"q\")\n\nReturn the q-expansion of the unique normalized eigenform of weight  k and level 1 for k in {12, 16, 18, 20, 22, 26} as a power series  to precision prec in var over ZZ.\n\nThese eigenforms are the normalized generators for the six  one-dimensional spaces of cusp forms of level 1.  \n\nArguments\n\nk::Integer: weight of the unique normalized eigenform\nprec::Integer=10: precision of the output \nvar::String=\"q\": variable name\n\nExamples still missing\n\n\n\n\n\n"
},

{
    "location": "#ModularForms.eisenstein_series_qexp",
    "page": "ModularForms.jl Documentation",
    "title": "ModularForms.eisenstein_series_qexp",
    "category": "function",
    "text": "eisenstein_series_qexp(k, prec=10, K=QQ, var=\"q\", normalization=\"linear\")\n\nReturn the q-expansion of the normalized weight k Eisenstein series  on the modular group to precision prec as a power series in the ring  K in variable var, using the given normalization. \n\nThree normalizations are available: \"linear\" (default), \"constant\",  and \"integral\". If the normalization is \"linear\" then the linear  coefficient is 1. If it is \"constant\" then the series will be  normalized to have constant term 1. If the normalization is \"integral\" then the series will be normalized to have integer coefficients and no  common factors.  Note: To prevent errors, the output will be in the ring QQ if the  normalization is \"linear\" or \"constant\".\n\nArguments\n\nk::Integer: even positive integer, weight of the Eisenstein series\nprec::Integer=10: precision of the output \nK=ZZ: base ring of the output \nvar::String=\"q\": variable name\nnormalization::String=\"linear\": normalization to use\n\nExamples still missing\n\n\n\n\n\n"
},

{
    "location": "#ModularForms.eta_quotient",
    "page": "ModularForms.jl Documentation",
    "title": "ModularForms.eta_quotient",
    "category": "function",
    "text": "eta_quotient(g, prec=10)\n\nReturn the eta-quotient relative to the given collection of integers g as a power series up to precision prec. \n\nThe input g is an array of pairs [[t1,r1], [t2,r2], ..., [ts,rs]]  where each tj is a positive integer and each rj a nonnegative integer.  An error is thrown if the sum of tj*rj for j in {1, ..., s} does not equal 24, as the function only applies to eta-quotients that are cusp forms. The output is a cusp form of integral weight.  \n\nArguments\n\ng::Array{Array}: an array of pairs, collection of integers\nprec::Integer=10: precision of the output\n\nExamples still missing\n\n\n\n\n\n"
},

{
    "location": "#ModularForms.victor_miller_basis",
    "page": "ModularForms.jl Documentation",
    "title": "ModularForms.victor_miller_basis",
    "category": "function",
    "text": "victor_miller_basis_poly(k, prec=10, cusp_only=false, var=\"q\")\n\nReturn the Victor Miller basis for modular forms of weight k and level 1 to precision prec as an array whose entries are power series in ZZ[[var]]. If cusp_only is true, then return only a basis for the cuspidal subspace.\n\nThe algorithm uses the proof of Lemma 2.20 from William A. Stein.\n\nArguments\n\nk::Integer: weight\nprec::Integer=10: precision of the output\ncusp_only::Boolean=false\nvar::String=\"q\": variable name\n\nExamples still missing\n\n\n\n\n\n"
},

{
    "location": "#ModularForms.hecke_operator_on_qexp",
    "page": "ModularForms.jl Documentation",
    "title": "ModularForms.hecke_operator_on_qexp",
    "category": "function",
    "text": "hecke_operator_on_qexp(f, n, k, prec=nothing)\n\nCompute the image of the q-expansion f of a modular form under the Hecke operator T_n of weight k. Return a power series to precision prec. \n\nArguments\n\nf::RelSeriesElem: q-expansion\nn::Integer: integer >=1 \nk::Integer: weight \nprec::Integer=nothing: precision of the output\n\nExamples still missing\n\n\n\n\n\n"
},

{
    "location": "#ModularForms.hecke_operator_on_basis",
    "page": "ModularForms.jl Documentation",
    "title": "ModularForms.hecke_operator_on_basis",
    "category": "function",
    "text": "hecke_operator_on_basis(B, n, k)\n\nCompute the matrix of the Hecke operator T_n of weight k relative to the given basis B of q-expansions for a space of  modular forms. \n\nArguments\n\nB::Array: array of q-expansions\nn::Integer: integer >=1 \nk::Integer: weight \n\nExamples still missing\n\n\n\n\n\n"
},

{
    "location": "#ModularForms.prime_range",
    "page": "ModularForms.jl Documentation",
    "title": "ModularForms.prime_range",
    "category": "function",
    "text": "prime_range(n)\n\nReturn an array consisting of all primes up to and including n.\n\n\n\n\n\n"
},

{
    "location": "#ModularForms.poly_to_power_series",
    "page": "ModularForms.jl Documentation",
    "title": "ModularForms.poly_to_power_series",
    "category": "function",
    "text": "poly_to_power_series(f, K, prec=10)\n\nConvert the polynomial f over K to a relative power series  to precision prec as a polynomial.\n\nThe power series has the same coefficient ring and variable  name as f. The input polynomial f must be correct up to  precision prec. \n\nArguments\n\n\' f::PolyElem: polynomial over K, correct up to prec\n\nK: base ring of the output\nprec::Integer=10: precision of the output\n\nExamples still missing\n\n\n\n\n\n"
},

{
    "location": "#Functions-1",
    "page": "ModularForms.jl Documentation",
    "title": "Functions",
    "category": "section",
    "text": "delta_qexp\ndelta_k_qexp\neisenstein_series_qexp\neta_quotient\nvictor_miller_basis\nhecke_operator_on_qexp\nhecke_operator_on_basis\nprime_range\npoly_to_power_series"
},

{
    "location": "#Index-1",
    "page": "ModularForms.jl Documentation",
    "title": "Index",
    "category": "section",
    "text": ""
},

]}
