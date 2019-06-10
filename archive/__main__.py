from sympy import *
from generation_rules import *
import itertools

from ifs import IFS
from interval import open_overlap, Interval, NetIntervals


def main():
    ifs1 = IFS(
            [(Rational(1,2),0,Rational(1,3)),
                (Rational(1,4),Rational(1,4),Rational(1,3)),
                (Rational(1,2),Rational(1,2),Rational(1,3))],
            GenRule1)
    ifs2 = IFS(
            [(Rational(1,3),0,Rational(1,4)),
                (Rational(1,5),Rational(4,15),Rational(1,4)),
                (Rational(1,3),Rational(7,15),Rational(1,4)),
                (Rational(1,5),Rational(4,5),Rational(1,4))],
            GenRule2)

    #  print(list(ifs2.new_gen(Rational(1,25))))
    all_nbhd_types = set()
    for alpha in [Rational(1,3),Rational(1,5),Rational(1,9),Rational(1,15),Rational(1,25),Rational(37,51)]:
        print(alpha)
        net = ifs2.net_interval(alpha)
        for nb in net.nbhd_types:
            all_nbhd_types.add(nb)
            print(nb)
        print()
    dct = {nb:idx for idx,nb in enumerate(all_nbhd_types)}
    print(dct)



if __name__ == "__main__":
    main()
