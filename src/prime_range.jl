function prime_range(n::Int)
  primes = fill(true, n)
  if n == 0
    return primes
  end
  primes[1] = false
  for p in 2:n
    primes[p] || continue
    for j in 2:div(n, p)
      primes[p*j] = false
    end
  end
  return findall(primes)
end
