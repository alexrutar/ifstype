from ifstype.examples import *
from ifstype import run_ifs, verify_wft, Generations
from ifstype.exact import Fraction
import cProfile


if __name__ == "__main__":
    #  ifs = wft_3()
    #  gens = Generations(ifs)
    #  tr_g = run_ifs(ifs, "output",with_gens=True,scale='wide',edge_labels='index')
    cProfile.run('Generations(wft_3()).compute_graph()',filename='profile/output.pstats')

    #  from ifstype.exact.symbolic import *



    #  syr = SymbolicRing(("p1","p2","p3"))
    #  p1,p2,p3 = syr.get_symbols()
    #  n = p1+p2+p3**2
    #  syr.set_eval({"p1":Rational(1,2),"p2":Rational(1,10),"p3":Rational(3,4)})

    #  print(0 == p1**2+2*p1*p2+p2**2- (p1+p2)**2)


