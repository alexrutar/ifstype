from .numeric import AlgebraicNumberFactory, Rational
from sympy import sqrt, Rational as sy_R, pi

import logging
import timeit
import cProfile

from .ifs_examples import *
from .algorithms import *

from .interval import *

def main():
    logging.basicConfig(filename='example.log',level=logging.DEBUG)
    #  ifs5()


    draw_infinite_type(ifs2m())


    #  expr = ((sy_R(9)+sqrt(sy_R(69)))/sy_R(18))**(sy_R(1,3))+((sy_R(9)-sqrt(sy_R(69)))/sy_R(18))**(sy_R(1,3))
    #  anf = AlgebraicNumberFactory(expr)
    #  x = anf.alpha()
    #  print(x==x)
    


if __name__ == "__main__":
    #  cProfile.run('main()')
    main()
