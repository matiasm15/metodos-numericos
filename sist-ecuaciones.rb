# Sustitucion regresiva: utilizado para encontrar una solucion aproximada X de un sistema triangular superior
# AX = b con A = (aij) de nxn donde A es invertible.
# Entrada: el orden n del sistema, los coeficientes aij con i = 1,2,...,n y j = 1,2,...,n y los terminos 
# independientes bi con i = 1,2,...,n.
# Salida: una solucion aproximada X = (x1,x2,...,xn).
def sustitucion_regresiva(n, a, b)
  n = n - 1

  x = Array.new(n + 1)
  x[n] = b[n].fdiv(a[n][n])

  (0..n - 1).reverse_each do |i|
    sumatoria = (i + 1..n).inject(0) { |sum, k| sum + a[i][k] * x[k] }
    x[i] = (b[i] - sumatoria).fdiv(a[i][i])
  end

  puts "Una solucion aproximada del sistema es X = #{x}"
  return x
end

# Eliminacion gaussiana con sustitucion regresiva: utiizado para obtener una solucion aproximada X de un sistema 
# de la forma
#                 E1: a11 * x1 + ... + a1n * xn = b1
#                 E2: a21 * x1 + ... + a2n * xn = b2
#                 En: an1 * x1 + ... + ann * xn = bn
# Entrada: el orden n del sistema, los coeficientes aij con i = 1,2,...,n y j = 1,2,...,n+1 de la matriz aumentada
# (A:b) con a(i,n+1) = bi con i = 1,2,...,n.
# Salida: una solucion aproximada X = (x1,x2,...,xn) del sistema dado o un mensaje.
def eliminacion_gaussianda(n, a)
  n = n - 1

  for j in (0..n - 1)
    k = (j..n).select { |i| !a[i][j].zero? }.first

    if k.nil?
      puts "El sistema no tiene solucion unica."
      return
    end

    if !k.eql?(j)
      a[k], a[j] = a[j], a[k]
    end

    for i in (j + 1..n)
      m = a[i][j].fdiv(a[j][j])
      resta = a[i].map.with_index { |x, index| x - m * a[j][index] }
      a[i].replace(resta)
    end
  end

  if a[n][n].zero?
    puts "El sistema no tiene solucion unica."
    return
  end

  x = Array.new(n + 1)
  x[n] = a[n][n + 1].fdiv(a[n][n])

  (0..n - 1).reverse_each do |i|
    sumatoria = (i + 1..n).inject(0) { |sum, k| sum + a[i][k] * x[k] }
    x[i] = (a[i][n + 1] - sumatoria).fdiv(a[i][i])
  end

  puts "Una solucion aproximada del sistema es X = #{x}"
  return x
end

# Factorizacion directa de Choleski: utilizada para factorizar una matriz A = (aij) de nxn donde A es real,
# simetrica y definida positiva en la forma A = L * L^T, donde L es triangular inferior.
# Entrada: la dimension n y los elementos aij de la matriz A.
# Salida: los elementos lij de la matriz L.
def choleski(n, a)
  n = n - 1

  l = Array.new(n + 1) { Array.new(n + 1) { 0 } }
  l[0][0] = Math::sqrt(a[0][0])

  (1..n).each do |i|
    l[i][0] = a[i][0].fdiv(l[0][0])
  end

  l[1][1] = Math::sqrt(a[1][1] - l[1][0] ** 2)

  (2..n).each do |i|
    (1..i - 1).each do |j|
      sumatoria = (0..j - 1).inject(0) { |sum, k| sum + l[i][k] * l[j][k] }
      l[i][j] = (a[i][j] - sumatoria).fdiv(l[j][j])
    end

    sumatoria = (0..i - 1).inject(0) { |sum, k| sum + l[i][k] ** 2 }
    l[i][i] = Math::sqrt(a[i][i] - sumatoria)
  end

  puts "Las componentes lij de la matriz L son #{l}."
  return l
end

# Metodo de Jacobi: utilizado para encontrar una solucion aproximada X de un sistema AX = b, con A = (aij) 
# de nxn donde A es real e invertible, b != 0, y aii != 0 para todo i = 1,2,...,n.
# Entrada: el orden n del sistema, las componentes (no nulas) aij con i,j = 1,2,...,n de la matriz A, las 
# componentes bi con i = 1,2,...,n del vector de terminos independientes, las componentes x0i con i = 1,2,...,n 
# de una aproximacion inicial XO = X^(O), una tolerancia error y un numero maximo de iteraciones n_max.
# Salida: una solucion aproximada X = (x1,x2,...,xn) o un mensaje.
def jacobi(n, a, b, x0, error, n_max)
  n = n - 1

  x = Array.new(n + 1)
  for k in (0..n_max)
    (0..n).each do |i|
      sumatoria = (0..n).inject(0) do |sum, j|
        if j.eql?(i)
          sum
        else
          sum + a[i][j] * x0[j]
        end
      end

      x[i] = (b[i] - sumatoria).fdiv(a[i][i])
    end

    resta = x.map.with_index { |xi, i| xi - x0[i] }
    modulo = Math::sqrt(resta.inject(0) { |sum, i| sum + i ** 2 })
    if modulo < error
      puts "Una solucion aproximada es X = #{x}."
      return x
    end

    x0.replace(x)
  end

  puts "Se alcanzo el numero maximo de iteraciones n_max pero no la tolerancia."
end

# Metodo de Gauss-Seidel: utilizado para encontrar una solucion aproximada X de un sistema AX = b, con A = (aij) 
# de nxn donde A es real e invertible, b != 0, y aii != 0 para todo i = 1,2,...,n.
# Entrada: el orden n del sistema, las componentes (no nulas) aij con i,j = 1,2,...,n de la matriz A, las 
# componentes bi con i = 1,2,...,n del vector de terminos independientes, las componentes x0i con i = 1,2,...,n 
# de una aproximacion inicial XO = X^(O), una tolerancia error y un numero maximo de iteraciones n_max.
# Salida: una solucion aproximada X = (x1,x2,...,xn) o un mensaje.
def gauss_seidel(n, a, b, x0, error, n_max)
  n = n - 1

  x = Array.new(n + 1)
  for k in (0..n_max)
    sumatoria = (1..n).inject(0) { |sum, j| sum + a[0][j] * x0[j] }
    x[0] = (b[0] - sumatoria).fdiv(a[0][0])

    (1..n - 1).each do |i|
      sumatoria_1 = (0..i - 1).inject(0) { |sum, j| sum + a[i][j] * x[j] }
      sumatoria_2 = (i + 1..n).inject(0) { |sum, j| sum + a[i][j] * x0[j] }
      x[i] = (b[i] - sumatoria_1 - sumatoria_2).fdiv(a[i][i])
    end

    sumatoria = (0..n - 1).inject(0) { |sum, j| sum + a[n][j] * x[j] }
    x[n] = (b[n] - sumatoria).fdiv(a[n][n])

    resta = x.map.with_index { |xi, i| xi - x0[i] }
    modulo = Math::sqrt(resta.inject(0) { |sum, i| sum + i ** 2 })
    if modulo < error
      puts "Una solucion aproximada es X = #{x}."
      return x
    end

    x0.replace(x)
  end

  puts "Se alcanzo el numero maximo de iteraciones n_max pero no la tolerancia."
end

# Metodo SOR: utilizado para encontrar una solucion aproximada X de un sistema AX = b, con A = (aij) 
# de nxn donde A es real e invertible, b != 0, y aii != 0 para todo i = 1,2,...,n dado un valor del
# parametro w con 0 < w < 2.
# Entrada: el orden n del sistema, las componentes (no nulas) aij con i,j = 1,2,...,n de la matriz A, las 
# componentes bi con i = 1,2,...,n del vector de terminos independientes, las componentes x0i con i = 1,2,...,n 
# de una aproximacion inicial XO = X^(O), un valor del parametro w, una tolerancia error y un numero maximo de
# iteraciones n_max.
# Salida: una solucion aproximada X = (x1,x2,...,xn) o un mensaje.
def sor(n, a, b, x0, w, error, n_max)
  n = n - 1

  x = Array.new(n + 1)
  for k in (0..n_max)
    sumatoria = (1..n).inject(0) { |sum, j| sum + a[0][j] * x0[j] }
    x[0] = (1 - w) * x0[0] + w * (b[0] - sumatoria).fdiv(a[0][0])

    (1..n - 1).each do |i|
      sumatoria_1 = (0..i - 1).inject(0) { |sum, j| sum + a[i][j] * x[j] }
      sumatoria_2 = (i + 1..n).inject(0) { |sum, j| sum + a[i][j] * x0[j] }
      x[i] = (1 - w) * x0[i] + w * (b[i] - sumatoria_1 - sumatoria_2).fdiv(a[i][i])
    end

    sumatoria = (0..n - 1).inject(0) { |sum, j| sum + a[n][j] * x[j] }
    x[n] = (1 - w) * x0[n] + w * (b[n] - sumatoria).fdiv(a[n][n])

    resta = x.map.with_index { |xi, i| xi - x0[i] }
    modulo = Math::sqrt(resta.inject(0) { |sum, i| sum + i ** 2 })
    if modulo < error
      puts "Una solucion aproximada es X = #{x}."
      return x
    end

    x0.replace(x)
  end

  puts "Se alcanzo el numero maximo de iteraciones n_max pero no la tolerancia."
end

# Metodo de Bairstow: utilizado para encontrar un factor cuadratico x^2 - u*x - v de un polinomio con coeficientes
# reales p(x) = a0 + a1*x + a2*x^2 + ... + an*x^n con an != 0 y n >= 2.
# Entrada: el grado n del polinomio p(x), los coeficientes a0,a1,...,an del polinomio p(x), unas aproximaciones iniciales
# u0 y v0 de u y v respectivamente, una tolerancia error y un numero maximo de iteraciones n_max.
# Salida: un factor cuadratico x^2 - u*x - v del polinomio p(x) o un mensaje.
def bairstow(n, a, u0, v0, error, n_max)
  b = Array.new(n + 1)
  c = Array.new(n + 1)

  b[n] = a[n]
  c[n] = 0
  c[n - 1] = a[n]
  for i in (1..n_max)
    b[n - 1] = a[n - 1] + u0 * b[n]
    (0..n - 2).reverse_each do |k|
      b[k] = a[k] + u0 * b[k + 1] + v0 * b[k + 2]
      c[k] = b[k + 1] + u0 * c[k + 1] + v0 * c[k + 2]
    end

    j = c[0] * c[2] - c[1] ** 2
    u1 = u0 + (c[1] * b[1] - c[2] * b[0]).fdiv(j)
    v1 = v0 + (c[1] * b[0] - c[0] * b[1]).fdiv(j)

    if b[1].abs < error and b[0].abs < error
      puts "Un factor cuadratico aproximado del polinomio dado, y los correspondientes valores r y s son x^2 - #{u1} * x - #{v1}, r = #{b[1]} y s = #{b[0]}."
      return u1, v1
    end

    u0 = u1
    v0 = v1
  end

  puts "Se alcanzo el numero maximo de iteraciones n_max pero no la tolerancia."
end
