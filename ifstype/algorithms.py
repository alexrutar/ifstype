from graph_tool.all import *
import pathlib

from .generations import Generations
from .info import GraphArtist, GraphWriter

class _Printer:
    def __init__(self, verbose):
        self.verbose = verbose

    def write(self, string):
        if self.verbose:
            print(string)

def run_ifs(ifs,foldername,with_gens=False,scale='wide',edge_labels='index',depth=2000,verbose=True):
    """Run the ifs and place output files in foldername"""
    p = pathlib.Path(f"{foldername}/")
    p.mkdir(parents=True, exist_ok=True) # create directory
    Pr = _Printer(verbose)

    Pr.write("\nComputing neighbour sets...")
    tr_g = Generations(ifs).compute_graph(depth=depth)
    if tr_g.is_wft:
        Pr.write(f"  IFS is weak finite type with {tr_g.g.num_vertices()} neighbour sets.\n")
    else:
        Pr.write(f"  Computed {tr_g.g.num_vertices()} neighbour sets without terminating; cannot guarantee weak finite type. Re-run with greater depth if desired.\n")
    draw = GraphArtist(tr_g)
    write = GraphWriter(tr_g)

    # possibly draw generations
    if with_gens:
        draw.gens(f"{foldername}/gens.pdf",scale=scale)
        Pr.write(f"- Net intervals saved to '{foldername}/gens.pdf'.")

    draw.graph(f"{foldername}/graph.pdf",edge_labels=edge_labels)
    Pr.write(f"- Graph saved to '{foldername}/graph.pdf'.")
    write.info_to_file(f"{foldername}/info.txt")
    Pr.write(f"- Info saved to '{foldername}/info.txt'.")

    return tr_g

def verify_wft(ifs,depth=1000):
    """Run the ifs, and determine if it is finite type or not."""
    ifs_gen = Generations(ifs)
    print(depth)
    return ifs_gen.verify_wft(depth=depth)[0]
