import logging
from .algorithms import *

def main():
    logging.basicConfig(filename='example.log',level=logging.DEBUG)
    #  transition_types(alphas=[Rational(1,3),Rational(1,9),Rational(1,27)])
    #  transition_types()
    #  test_uq_subdiv()
    example_draw()
    #  normalized_deltas()



if __name__ == "__main__":
    main()
