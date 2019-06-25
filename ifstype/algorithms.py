from .generations import FiniteType, InfiniteType
from .draw import Visual
from .rational import Rational

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
    diagram = Visual(gn,filename,1,scale=3)
    for alpha in gn.ifs.transition_gens(stop=stop):
        iv_net = gn.gen(alpha)
        diagram.interval(iv_net)
        diagram.net(iv_net)

    diagram.nb_set()
    diagram.show()

