# Regla de los Trapecios: utilizada para aproximar la integral I = Intregal(a..b){f(x)dx}.
# Entrada: los extremos a y b de la integral, un entero positivo n y f(x).
# Salida: una aproximacion AI de I.
def trapecios(a, b, n, &f)
  h = (b - a).fdiv(n)
  ai0 = f.call(a) + f.call(b)
  ai = 0

  (1..n - 1).each do |k|
    x = a + k * h
    ai += f.call(x)
  end

  ai = (ai0 + 2 * ai) * h.fdiv(2)

  puts "Un valor aproximado de la integral para #{n} subintervalos es #{ai}."
  return ai
end

# Regla de Simpson (1 / 3): utilizada para aproximar la integral I = Intregal(a..b){f(x)dx}.
# Entrada: los extremos a y b de la integral, un entero positivo par n = 2m y f(x).
# Salida: una aproximacion AI de I.
def simpson(a, b, n, &f)
  h = (b - a).fdiv(n)
  ai0 = f.call(a) + f.call(b)
  ai1 = 0
  ai2 = 0

  (1..n - 1).each do |i|
    x = a + i * h
    if i.odd?
      ai1 += f.call(x)
    else
      ai2 += f.call(x)
    end
  end

  ai = (ai0 + 4 * ai1 + 2 * ai2) * h.fdiv(3)

  puts "Un valor aproximado de la integral para #{n} subintervalos es #{ai}."
  return ai
end

# Metodo de Romberg: utilizada para aproximar la integral I = Intregal(a..b){f(x)dx}.
# Entrada: los extremos a y b de la integral, un entero positivo n y f(x).
# Salida: una aproximacion de I.
def romberg(a, b, n, &f)
  n = n - 1
  r = Array.new(n + 1) { Array.new(n + 1) }

  r[0][0] = 0.5 * (b - a) * (f.call(a) + f.call(b))

  (1..n).each do |i|
    h = (b - a).fdiv(2 ** i)
    sumatoria = (1..2 ** (i - 1)).inject(0) { |sum, k| sum + f.call(a + (2 * k - 1) * h) }
    r[i][0] = 0.5 * r[i - 1][0] + h * sumatoria
  end

  (1..n).each do |i|
    (1..i).each do |j|
      r[i][j] = r[i][j - 1] + (1.fdiv(4 ** j - 1)) * (r[i][j - 1] - r[i - 1][j - 1])
    end
  end

  puts "Un valor aproximado de la integral es #{r[n][n]}."
  return r[n][n]
end
