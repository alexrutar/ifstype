from .ifs_examples import *
from .generations import FiniteType, InfiniteType
from .draw import Visual

def draw_finite_type(ifs,filename="example.pdf"):
    "An example drawing illustrating the interval and net methods"
    gn = FiniteType(ifs)
    diagram = Visual(gn,filename,1,scale=3)
    for alpha in gn.ifs.transition_gens(stop=gn.new_transition_stop):
        iv_net = gn.gen(alpha)
        diagram.interval(iv_net)
        diagram.net(iv_net)

    diagram.nb_set()
    diagram.show()

def draw_infinite_type(ifs, stop=Rational(1,25),filename="example.pdf"):
    "An example drawing illustrating the interval and net methods"
    gn = InfiniteType(ifs)
    diagram = Visual(gn,filename,1,scale=3)
    for alpha in gn.ifs.transition_gens(stop=stop):
        iv_net = gn.gen(alpha)
        diagram.interval(iv_net)
        diagram.net(iv_net)

    diagram.nb_set()
    diagram.show()

