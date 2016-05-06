# Ecuacion: dx/dt + k*x = sen(t) con x(0) = 1.
# Utiliza esquema adelantado.
def solucion(k, delta_t, n)
  x_delta_t_anterior = 1 - delta_t * k

  resultado = [
    [delta_t, x_delta_t_anterior]
  ]

  for i in (2..n)
    t = i * delta_t
    x_delta_t_anterior = delta_t * Math::sin(t) - x_delta_t_anterior * (delta_t * k - 1)
    resultado << [t, x_delta_t_anterior]
  end

  return resultado
end
