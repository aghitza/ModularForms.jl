function prime_range(n)
  primes = fill(true, n)
  primes[1] = false
  for p in 2:n
    primes[p] || continue
    for j in 2:div(n, p)
      primes[p*j] = false
    end
  end
  findall(primes)
end
