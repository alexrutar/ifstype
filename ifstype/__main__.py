import logging
from .algorithms import *
from .interval import *

def main():
    logging.basicConfig(filename='example.log',level=logging.DEBUG)

    example_draw()
    #  test_gens()
    #  normalized_deltas()



if __name__ == "__main__":
    main()
