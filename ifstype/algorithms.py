from .generations import FiniteType, InfiniteType
from .draw import Visual
from .rational import Rational
from .interval import Interval

from bisect import bisect, bisect_right


def run_finite_type(ifs,file_tag):
    "An example drawing illustrating the interval and net methods"
    gn = FiniteType(ifs)
    diagram = Visual(gn,file_tag + "_generations.pdf",1,scale=3)
    with open(file_tag +"_info.txt",'w+') as f:
        f.write("IFS info:\n")
        f.write(str(ifs))
        f.write("\n\nNeighbour sets:\n")
        f.write(str(gn.nb_mgr))
    for alpha in gn.ifs.transition_gens(stop=gn.new_transition_stop):
        iv_net = gn.gen(alpha)
        diagram.interval(iv_net)
        diagram.net(iv_net)

    diagram.nb_set()
    diagram.show()
    gn.transition_graph.draw(file_tag + "_graph.pdf")

def run_infinite_type(ifs, stop=Rational(1,25),filename="example.pdf"):
    "An example drawing illustrating the interval and net methods"
    gn = InfiniteType(ifs)
    diagram = Visual(gn,filename,2,scale=3)
    for alpha in gn.ifs.transition_gens(stop=stop):
        iv_net = gn.gen(alpha)
        diagram.interval(iv_net)
        #  diagram.net(iv_net)

    diagram.nb_set()
    diagram.show()

def get_alpha_density(ifs, count=10):
    gn = InfiniteType(ifs)
    for alpha in gn.ifs.transition_gens(count=count):
        fun_net = gn.gen(alpha)
        mx = 0
        for idx,iv in enumerate(fun_net.net):
            rep = iv.a + alpha
            ct = bisect_right(fun_net.net,Interval.closed_infty(rep))+1
            mx = max(mx, ct-idx)
        print(mx)

