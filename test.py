from ifstype.examples import *
from ifstype import run_ifs, verify_fnc, Generations, run_ifs_gens
from ifstype.exact import Fraction
from ifstype.graph import SAdjacencyMatrix
from ifstype.ifs import TransitionMatrix
from ifstype.info import GraphArtist, GraphWriter
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

def test_fnc():
    ifs = fnc_8()
    return run_ifs(ifs,"out",scale='wide',with_gens=True,verbose=True,depth=5000)

def test_reduce():
    tr_g = test_fnc()
    ga = GraphArtist(tr_g)
    gw = GraphWriter(tr_g)

    ga.graph("graph_0.pdf")
    gw.info_to_file("info_0.txt")
    tr_g.reduce_transition_graph()
    ga.graph("graph_1.pdf")
    gw.info_to_file("info_1.txt")



if __name__ == "__main__":
    #  test_reduce()
    #  test_reduce()
    test_fnc()
    #  ifs = wft_5()
    #  tr_g = run_ifs(ifs, f"fnc_basic",with_gens=True,scale="relative",depth=2000)

    #  ifs = wft_3()
    #  tr_g = run_ifs(ifs, f"fnc_ex",with_gens=True,scale="relative",depth=2000)

    #  ifs = eft_5()
    #  tr_g = run_ifs(ifs, f"equi_ex",with_gens=False,scale="relative",depth=2000)

