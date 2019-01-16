# List of normalised generators for the six one-dimensional
# spaces of cusp forms of level 1.

using Nemo

R, q = PowerSeriesRing(ZZ, 20, "q")

GENERATOR = Dict(
  12 => q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 - 6048*q^6 - 16744*q^7 + 84480*q^8 - 113643*q^9 - 115920*q^10 + 534612*q^11 - 370944*q^12 - 577738*q^13 + 401856*q^14 + 1217160*q^15 + 987136*q^16 - 6905934*q^17 + 2727432*q^18 + 10661420*q^19 + O(q^20),
  16 => q + 216*q^2 - 3348*q^3 + 13888*q^4 + 52110*q^5 - 723168*q^6 + 2822456*q^7 - 4078080*q^8 - 3139803*q^9 + 11255760*q^10 + 20586852*q^11 - 46497024*q^12 - 190073338*q^13 + 609650496*q^14 - 174464280*q^15 - 1335947264*q^16 + 1646527986*q^17 - 678197448*q^18 + 1563257180*q^19 + O(q^20),
  18 => q - 528*q^2 - 4284*q^3 + 147712*q^4 - 1025850*q^5 + 2261952*q^6 + 3225992*q^7 - 8785920*q^8 - 110787507*q^9 + 541648800*q^10 - 753618228*q^11 - 632798208*q^12 + 2541064526*q^13 - 1703323776*q^14 + 4394741400*q^15 - 14721941504*q^16 - 5429742318*q^17 + 58495803696*q^18 + 1487499860*q^19 + O(q^20),
  20 => q + 456*q^2 + 50652*q^3 - 316352*q^4 - 2377410*q^5 + 23097312*q^6 - 16917544*q^7 - 383331840*q^8 + 1403363637*q^9 - 1084098960*q^10 - 16212108*q^11 - 16023861504*q^12 + 50421615062*q^13 - 7714400064*q^14 - 120420571320*q^15 - 8939761664*q^16 + 225070099506*q^17 + 639933818472*q^18 - 1710278572660*q^19 + O(q^20),
  22 => q - 288*q^2 - 128844*q^3 - 2014208*q^4 + 21640950*q^5 + 37107072*q^6 - 768078808*q^7 + 1184071680*q^8 + 6140423133*q^9 - 6232593600*q^10 - 94724929188*q^11 + 259518615552*q^12 - 80621789794*q^13 + 221206696704*q^14 - 2788306561800*q^15 + 3883087691776*q^16 + 3052282930002*q^17 - 1768441862304*q^18 - 7920788351740*q^19 + O(q^20),
  26 => q - 48*q^2 - 195804*q^3 - 33552128*q^4 - 741989850*q^5 + 9398592*q^6 + 39080597192*q^7 + 3221114880*q^8 - 808949403027*q^9 + 35615512800*q^10 + 8419515299052*q^11 + 6569640870912*q^12 - 81651045335314*q^13 - 1875868665216*q^14 + 145284580589400*q^15 + 1125667983917056*q^16 - 2519900028948078*q^17 + 38829571345296*q^18 - 6082056370308940*q^19 + O(q^20))
