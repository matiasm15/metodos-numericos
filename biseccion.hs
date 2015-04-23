
biseccion f a b error limite | f a == 0 = a
                        	   | f b == 0 = b
                         	   | otherwise = biseccionRecursivo 1 f a b error limite

biseccionRecursivo n f a b error limite | (mult == 0) || (error >= abs (b - a)) || (n >= limite) = c
                                        | mult > 0 = biseccionRecursivo (n + 1) f c b error limite
                                        | mult < 0 = biseccionRecursivo (n + 1) f a c error limite
  where mult = f a * f c
        c = (a + b) / 2
