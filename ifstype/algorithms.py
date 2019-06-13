from sympy import Rational, log, Max, sqrt
import itertools
from collections import defaultdict

from .interval import NetInterval,Interval
from .ifs import IFS,CtrFunc,C
from .gens import Generations
from .draw import Visual
from .select import *

def ifs0():
    # doesn't work yet
    r = (sqrt(5)-1)/2
    return Generations(IFS(
        (CtrFunc(r,Rational(0)),Rational(1,2)),
        (CtrFunc(r,Rational(1)-r),Rational(1,2))),
        finite_type=True)
def ifs1():
    # finite type III
    return Generations(IFS(
        (CtrFunc(Rational(1,2),0), Rational(1,3)),
        (CtrFunc(Rational(1,4),Rational(1,4)), Rational(1,3)),
        (CtrFunc(Rational(1,2),Rational(1,2)), Rational(1,3))),
        finite_type=True)

def ifs2():
    # finite type IV
    return Generations(IFS(
        (CtrFunc(Rational(1,3),0), Rational(1,4)),
        (CtrFunc(Rational(1,5),Rational(4,15)), Rational(1,4)),
        (CtrFunc(Rational(1,3),Rational(7,15)), Rational(1,4)),
        (CtrFunc(Rational(1,5),Rational(4,5)), Rational(1,4))),
        finite_type=True)

def ifs3():
    # not finite type
    return Generations(IFS(
        (CtrFunc(Rational(1,3),0), Rational(1,4)),
        (CtrFunc(Rational(1,4),Rational(4,15)), Rational(1,4)),
        (CtrFunc(Rational(1,3),Rational(7,15)), Rational(1,4)),
        (CtrFunc(Rational(1,5),Rational(4,5)), Rational(1,4))),
        finite_type=False)

def ifs4():
    # finite type IV
    return Generations(IFS(
        (CtrFunc(Rational(1,3),0),Rational(1,4)),
        (CtrFunc(Rational(1,4),Rational(1,4)),Rational(1,4)),
        (CtrFunc(Rational(1,3),Rational(1,2)),Rational(1,4)),
        (CtrFunc(Rational(1,4),Rational(3,4)),Rational(1,4))),
        finite_type=True)

def test_uq_subdiv():
    gn = ifs2()
    diagram = Visual(gn,"example.pdf",1,scale=2)
    for alpha in [Rational(1,2),Rational(1,4),Rational(1,8),Rational(1,16)]:
        iv_net = gn.gen(alpha)
        diagram.interval(iv_net)
        diagram.net(iv_net)
    diagram.nb_set()
    diagram.show()

def example_draw():
    "An example drawing illustrating the interval and net methods"
    gn = ifs2()
    diagram = Visual(gn,"example.pdf",1,scale=3)
    for alpha in [C.n_base,Rational(1),Rational(1,3),Rational(1,5),Rational(1,9),Rational(1,15)]:
        iv_net = gn.gen(alpha)
        diagram.interval(iv_net)
        diagram.net(iv_net)

    diagram.nb_set()
    diagram.show()
    #  iv = NetInterval(Rational(1,3),Rational(1,3),Rational(7,15))
    #  print(f"Transition Type of interval {iv} : {gn.ttype(iv)}")


def global_min_delta(alphas=map(lambda n: Rational(1,5)**n,range(4))):
    "Compute the smallest value of delta among the choices of alpha: computes the entire generation, so it's pretty slow"
    gn = ifs2()
    diagram = Visual(gn,"example.pdf",1,scale=10)
    for alpha in alphas:
        iv_net = gn.gen(alpha)
        iv = iv_net.argmin_delta()
        print("{} : {}".format(iv,float(iv.delta/alpha)))
        diagram.interval(iv_net)
        diagram.net(iv_net,highlight=[iv])
    diagram.nb_set()
    diagram.show()


def local_min_delta(alphas=map(lambda n: Rational(1,5)**n,range(50))):
    "Compute the smallest value of delta; only looks in the current minimum for the next smallest"
    gn = ifs2()
    for (mid,iv_net) in gn.gen_from_select(alphas,select.min_delta):
        print("{}".format(gn.nb_set_type(mid)))
    gn.print_nb_set_types()


def select_point(point=Rational(3,7),alphas=map(lambda n: Rational(1,5)**n,range(50))):
    "test the select algorithm checking a point. To make a minimum delta choice, use select.min_delta"
    gn = ifs2()
    for (mid,iv_net) in gn.gen_from_select(alphas,xval(point)):
        print("{}".format(gn.nb_set_type(mid)))
    gn.print_nb_set_types() 

def test_neighbour_size(alphas=map(lambda n: Rational(1,3)**n,range(30))):
    gn = ifs2()
    cur_max = 0
    for (mid,iv_net) in gn.gen_from_select(alphas,random):
        left,right = iv_net.adjacent(mid)
        if left.is_empty():
            ld = 0
        else:
            ld = abs(log(left.delta/mid.delta))
        if right.is_empty():
            rd = 0
        else:
            rd = abs(log(right.delta/mid.delta))
        print("{} : {} {}".format(mid.alpha, float(ld), float(rd)))
        cur_max = Max(cur_max,ld,rd)
    print(cur_max)

def transition_types(start=Rational(1,5),alphas=[Rational(1,9),Rational(1,15),Rational(1,25)]):
    gn = ifs2()
    out = defaultdict(set)
    info = defaultdict(set)
    diagram = Visual(gn,"example.pdf",1,scale=3)
    iv_net = gn.gen(start)
    diagram.interval(iv_net)
    diagram.net(iv_net)
    for idx,alpha in map(lambda x:(x[0]+1,x[1]),enumerate(alphas)):
        diagram.interval(gn.gen(alpha))
        diagram.net(gn.gen(alpha))
        for net_iv in iv_net:
            new = tuple(gn.nb_set_type(chld) for chld in gn.children(alpha,net_iv)) # tuple of children
            chr_vec = (gn.nb_set_type(net_iv), net_iv.delta, idx)
            out[chr_vec].add(new)
            info[chr_vec].add((new,net_iv,alpha))

    for k,v in out.items():
        print(f"{k} {v}")
        #  if len(v) > 1:
            #  print(f"\n{k}:\n -{info[k]}")


    diagram.nb_set()
    diagram.show()

def normalized_deltas(alphas=[Rational(1,9),Rational(1,15),Rational(1,3),Rational(1,5),Rational(1,27),Rational(1,25),Rational(1,625)]):
    gn = ifs2()
    vals = set()
    for iv_net in map(gn.gen, alphas):
        for net_iv in iv_net:
            vals.add(net_iv.delta/iv_net.alpha)
    print(vals)
    print(len(vals))
