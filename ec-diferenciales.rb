# Metodo de Euler: utilizado para encontrar una aproximacion de la solucion del PVI, con solucion unica
#                    y' = f(t, y)
#                    y(t0) = y0
# evaluado en t1.
# Entrada: los valores iniciales t0 y y0, el valor t1, un entero m y la funcion f(t, y).
# Salida: una aproximacion de la solucion y(t) evaluado en t1.
def euler(t0, y0, m, t1, &f)
  h = (t1 - t0).fdiv(m)

  t = Array.new(m + 1)
  y = Array.new(m + 1)

  t[0] = t0
  y[0] = y0

  (1..m).each do |i|
    t[i] = t[i - 1] + h
    y[i] = y[i - 1] + h * f.call(t[i - 1], y[i - 1])
  end

  puts "Una aproximacion de la solucion del PVI dado es y(#{t[m]}) = #{y[m]}."
  return y[m]
end

# Metodo del punto medio: utilizado para encontrar una aproximacion de la solucion del PVI, con solucion unica
#                    y' = f(t, y)
#                    y(t0) = y0
# evaluado en t1.
# Entrada: los valores iniciales t0 y y0, el valor t1, un entero m y la funcion f(t, y).
# Salida: una aproximacion de la solucion y(t) evaluado en t1.
def punto_medio(t0, y0, m, t1, &f)
  h = (t1 - t0).fdiv(m)

  t = Array.new(m + 1)
  y = Array.new(m + 1)

  t[0] = t0
  y[0] = y0

  (1..m).each do |i|
    t[i] = t[i - 1] + h
    y[i] = y[i - 1] + h * f.call(t[i - 1] + 0.5 * h, y[i - 1] + 0.5 * h * f.call(t[i - 1], y[i - 1]))
  end

  puts "Una aproximacion de la solucion del PVI dado es y(#{t[m]}) = #{y[m]}."
  return y[m]
end

# Metodo de Heun: utilizado para encontrar una aproximacion de la solucion del PVI, con solucion unica
#                    y' = f(t, y)
#                    y(t0) = y0
# evaluado en t1.
# Entrada: los valores iniciales t0 y y0, el valor t1, un entero m y la funcion f(t, y).
# Salida: una aproximacion de la solucion y(t) evaluado en t1.
def heun(t0, y0, m, t1, &f)
  h = (t1 - t0).fdiv(m)

  t = Array.new(m + 1)
  y = Array.new(m + 1)

  t[0] = t0
  y[0] = y0

  (1..m).each do |i|
    t[i] = t[i - 1] + h

    k1 = h * f.call(t[i - 1], y[i - 1])
    k2 = h * f.call(t[i], y[i - 1] + k1)
    y[i] = y[i - 1] + 0.5 * (k1 + k2)
  end

  puts "Una aproximacion de la solucion del PVI dado es y(#{t[m]}) = #{y[m]}."
  return y[m]
end

# Metodo de Ralston: utilizado para encontrar una aproximacion de la solucion del PVI, con solucion unica
#                    y' = f(t, y)
#                    y(t0) = y0
# evaluado en t1.
# Entrada: los valores iniciales t0 y y0, el valor t1, un entero m y la funcion f(t, y).
# Salida: una aproximacion de la solucion y(t) evaluado en t1.
def ralston(t0, y0, m, t1, &f)
  h = (t1 - t0).fdiv(m)

  t = Array.new(m + 1)
  y = Array.new(m + 1)

  t[0] = t0
  y[0] = y0

  (1..m).each do |i|
    t[i] = t[i - 1] + h

    k1 = h * f.call(t[i - 1], y[i - 1])
    k2 = h * f.call(t[i - 1] + h * 2.fdiv(3), y[i - 1] + k1 * 2.fdiv(3))
    y[i] = y[i - 1] + 0.25 * (k1 + 3 * k2)
  end

  puts "Una aproximacion de la solucion del PVI dado es y(#{t[m]}) = #{y[m]}."
  return y[m]
end

# Metodo de Runge-Kutta: utilizado para encontrar una aproximacion de la solucion del PVI, con solucion unica
#                    y' = f(t, y)
#                    y(t0) = y0
# evaluado en t1.
# Entrada: los valores iniciales t0 y y0, el valor t1, un entero m y la funcion f(t, y).
# Salida: una aproximacion de la solucion y(t) evaluado en t1.
def runge_kutta(t0, y0, m, t1, &f)
  h = (t1 - t0).fdiv(m)

  t = Array.new(m + 1)
  y = Array.new(m + 1)

  t[0] = t0
  y[0] = y0

  (1..m).each do |i|
    t[i] = t[i - 1] + h

    k1 = f.call(t[i - 1], y[i - 1])
    k2 = f.call(t[i - 1] + h.fdiv(2), y[i - 1] + h.fdiv(2) * k1)
    k3 = f.call(t[i - 1] + h.fdiv(2), y[i - 1] + h.fdiv(2) * k2)
    k4 = f.call(t[i - 1] + h, y[i - 1] + h * k3)
    y[i] = y[i - 1] + h.fdiv(6) * (k1 + 2 * k2 + 2 * k3 + k4)
  end

  puts "Una aproximacion de la solucion del PVI dado es y(#{t[m]}) = #{y[m]}."
  return y[m]
end
