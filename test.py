from ifstype.examples import *
from ifstype import run_ifs, verify_fnc, Generations, run_ifs_gens
from ifstype.exact import Fraction
from ifstype.graph import SAdjacencyMatrix
import numpy as np
import cProfile


def test_adj_mat():
    a = Fraction(1,4)
    b = Fraction(1,3)
    # TODO: weird issue printing sums with 1... in SymbolicRing
    tup_arr = [[(a,a),(a,)],[(1,a),(b,)]]

    adj = SAdjacencyMatrix(tup_arr)
    print(adj)
    print(adj.compute_s_val())
    #  print(compute_s_val(arr))
    #  print(adj)

if __name__ == "__main__":
    #  ifs = wft_5()
    #  tr_g = run_ifs(ifs, f"fnc_basic",with_gens=True,scale="relative",depth=2000)

    #  ifs = wft_3()
    #  tr_g = run_ifs(ifs, f"fnc_ex",with_gens=True,scale="relative",depth=2000)

    #  ifs = eft_5()
    #  tr_g = run_ifs(ifs, f"equi_ex",with_gens=False,scale="relative",depth=2000)
    ifs = wft_5()
    #  tr_g = run
    tr_g = run_ifs(ifs,"out",scale='wide',with_gens=True,verbose=True,depth=300)

