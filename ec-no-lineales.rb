# Biseccion: utilizado para encontrar una aproximacion de una raiz perteneciente a (a,b) de una ecuacion f(x) = 0,
# donde f es una funcion continua en [a,b] y f(a)f(b) < 0.
# Entrada: los extremos a, b del intervalo, una tolerancia error, un numero maximo de iteraciones n_max y f(x).
# Salida: una raiz aproximada o un mensaje.
def biseccion(a, b, error, n_max, &f)
  for i in (1..n_max)
    c = (a + b).fdiv(2)

    if f.call(c).zero? or ((b - a).fdiv(2) < error)
      puts "Una raiz aproximada de la funcion es #{c}"
      puts "f(#{c}) = #{f.call(c)}"
      return c
    end

    if (f.call(a) * f.call(c)) < 0
      b = c
    else
      a = c
    end
  end

  puts "Se alcanzo el numero maximo de iteraciones n_max pero no la tolerancia."
end

# Punto fijo: utilizado para encontrar una aproximacion de un punto fijo de una funcion f, dada una aproximacion
# inicial x0.
# Entrada: una aproximacion inicial x0, una tolerancia error, un numero maximo de iteraciones n_max y f(x).
# Salida: un punto fijo aproximado o un mensaje.
def punto_fijo(x0, error, n_max, &f)
  for i in (1..n_max)
    c = f.call(x0)

    if (c - x0).abs < error
      puts "Un punto fijo aproximado de la funcion dada es #{c}"
      puts "f(#{c}) = #{f.call(c)}"
      return c
    end

    x0 = c
  end

  puts "Se alcanzo el numero maximo de iteraciones n_max pero no la tolerancia."
end

# Newton-Raphson: utilizado para encontrar una aproximacion de una raiz de una ecuacion f(x) = 0, conocida una
# aproximacion inicial x0.
# Entrada: una aproximacion inicial x0, una tolerancia error, un numero maximo de iteraciones n_max, f(x) y f'(x).
# Salida: una raiz aproximada o un mensaje.
def newton_raphson(x0, error, n_max, f, f_derivada)
  for i in (1..n_max)
    e = f.call(x0)
    d = f_derivada.call(x0)

    if d.zero?
      puts "No se puede continuar con el metodo."
      return
    end

    c = x0 - e.fdiv(d)

    if (f.call(c).abs < error) or (e.fdiv(d).abs < error)
      puts "Una raiz aproximada de la funcion es #{c}"
      puts "f(#{c}) = #{f.call(c)}"
      return c
    end

    x0 = c
  end

  puts "Se alcanzo el numero maximo de iteraciones n_max pero no la tolerancia."
end

# Horner: utilizado para evaluar un polinomio con coeficientes reales p(x) = a0 + a1*x + a2*x^2 + ... + an*x^n y su
# derivada en un numero real z.
# Entrada: el grado n y los coeficientes a0,a1,...,an del polinomio p(x), y el numero real z.
# Salida: b0 = p(z) y c = p'(z).
def horner(n, a, z)
  b = Array.new(n + 1)
  b[n] = a[n]
  c = a[n]

  (1..n - 1).reverse_each do |j|
    b[j] = a[j] + z * b[j + 1]
    c = b[j] + z * c
  end

  b[0] = a[0] + z * b[1]

  puts "p(#{z}) = #{b[0]}"
  puts "p'(#{z}) = #{c}"
  return b[0], c
end

# Newton-Raphson combinado con Horner: utilizado para encontrar un cero aproximado del polinomio
# p(x) = a0 + a1*x + a2*x^2 + ... + an*x^n.
# Entrada: el grado n y los coeficientes a0,a1,...,an del polinomio p(x), una aproximación inicial x0, una
# toleracion error y un numero maximo de iteraciones n_max.
# Salida: un cero aproximado del polinomio p(x) o un mensaje.
def newton_raphson_horner(n, a, x0, error, n_max)
  b = Array.new(n + 1)

  for i in (1..n_max)
    b[n] = a[n]
    c = a[n]

    (1..n - 1).reverse_each do |j|
      b[j] = a[j] + x0 * b[j + 1]
      c = b[j] + x0 * c
    end

    b[0] = a[0] + x0 * b[1]

    if c.zero?
      puts "No se puede continuar con el metodo porque se anulo p'(#{x0})."
      return
    end

    x1 = x0 - b[0].fdiv(c)

    if (b[0].abs < error) or (b[0].fdiv(c).abs < error)
      puts "Una raiz aproximada de p(x) = 0 es #{x1}"
      return x1
    end

    x0 = x1
  end

  puts "Se alcanzo el numero maximo de iteraciones n_max pero no la tolerancia."
end

# Secante: utilizado para encontrar una aproximacion de una raiz de una ecuacion f(x) = 0 conocidas
# dos aproximaciones iniciales x0 y x1.
# Entrada: dos aproximaciones iniciales x0 y xl, una tolerancia error, un numero maximo de iteraciones n_max
# y f(x).
# Salida: una raiz aproximada o un mensaje.
def secante(x0, x1, error, n_max, &f)
  y0 = f.call(x0)
  y1 = f.call(x1)
  for i in (2..n_max)
    if (y1 - y0).zero?
      puts "No se puede aplicar el metodo, porque el denominador en la formula de la secante se anulo."
      return
    end

    x2 = x1 - y1 * (x1 - x0).fdiv(y1 - y0)

    if (x2 - x1).abs < error
      puts "Una aproximacion de una raiz de la ecuacion dada es #{x2}."
      puts "f(#{x2}) = #{f.call(x2)}"
      return x2
    end

    x0 = x1
    y0 = y1
    x1 = x2
    y1 = f.call(x1)
  end

  puts "Se alcanzo el numero maximo de iteraciones n_max pero no la tolerancia."
end
