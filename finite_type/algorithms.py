from graph_tool.all import *
import pathlib

from .ifs_tree import Generations, NetInterval
from .ifs import IFS, CtrFunc
from .numerics.rational import Rational

def run_ifs(ifs,foldername,with_gens=False):
    """Run the ifs and place output files in foldername"""
    p = pathlib.Path(f"{foldername}/")
    p.mkdir(parents=True, exist_ok=True) # create directory

    ifs_gen = Generations(ifs)
    if with_gens:
        net_ivs = ifs_gen.compute_graph(compute_net_iv=True)
        ifs_gen.draw_gens(net_ivs,f"{foldername}/gens.pdf",scale='wide')
    else:
        ifs_gen.compute_graph()
    ifs_gen.draw_graph(f"{foldername}/graph.pdf",edge_labels='index')
    ifs_gen.write_info(f"{foldername}/info.txt")
    return ifs_gen.tr_graph
