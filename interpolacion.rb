# Diferencias divididas progresivas: utilizado para obtener los coeficientes b0,b1,...,bn de la forma de
# Newton del polinomio interpolante usando diferencias divididas progresivas, conocidos n + 1 puntos
# (x0,y0),(x1,y1),...,(xn,yn) donde x0,x1,...,xn son numeros distintos.
# Entrada: los valores x0,x1,...,xn y y0,y1,...,yn.
# Salida: los coeficientes b0,b1,...,bn de la forma progresiva de Newton del polinomio interpolante
# pn(x) = b0 + b1 * (x - x0) + b2 * (x - x0) * (x - x1) + ... + bn * (x - x0) * (x - x1) * ... * (x - xn-1)
def diferencias_divididas(x, y)
  n = x.size - 1
  fx = Array.new(n + 1) { Array.new(n + 1) }

  (0..n).each do |i|
    fx[i][i] = y[i]
  end

  (1..n).each do |i|
    (0..n - i).each do |j|
      fx[j][j + i] = (fx[j + 1][j + i] - fx[j][j + i - 1]).fdiv(x[j + i] - x[j])
    end
  end

  # Los coefiecientes b[i] son los valores fx[0][i].
  puts "Los coeficientes del polinomio interpolante progresivo de Newton son #{fx[0]}."
  return fx[0]
end

# Método de Lagrange: utilizado para obtener el valor del polinomio interpolador de Lagrange evaluado en un
# punto px, conocidos n + 1 puntos (x0,y0),(x1,y1),...,(xn,yn) donde x0,x1,...,xn son numeros distintos.
# Entrada: el punto px, los valores x0,x1,...,xn y y0,y1,...,yn.
# Salida: el valor py del polinomio interpolador de Lagrange evaluado en px.
def lagrange(px, x, y)
  n = x.size - 1
  
  py = (0..n).inject(0) do |sum, i|
    range = (0..n).to_a
    range.delete(i)
    sum + y[i] * range.inject(1) { |prd, j| prd * (px - x[j]).fdiv(x[i] - x[j]) }
  end

  puts "El valor py del polinomio interpolador de Lagrange evaluado en #{px} es #{py}."
  return py
end
