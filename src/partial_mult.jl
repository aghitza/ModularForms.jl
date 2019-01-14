function partial_mult(f, g, n, i=0, s=0)
  x = gen(parent(f))
  quo = 0
  rem = s
  for j in 0:i
    fs = slice(f, j*n, (j+1)*n-1)
    gs = slice(g, (i-j)*n, (i-j+1)*n-1)
    p = fs*gs
    q, r = divrem(p, x^n)
    quo += q
    rem += r
  end
  return (rem, quo)
end

function slice(f, a, b)
  x = gen(parent(f))
  q, f = divrem(f, x^(b+1))
  f, r = divrem(f, x^a)
  return f
end
