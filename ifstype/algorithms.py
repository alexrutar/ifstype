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

def run_ifs_gens(ifs,filename,explore,scale='wide',verbose=True,depth=2000):
    """Compute variations on the "generations" diagram for the given IFS.
    Run the ifs and place output files in foldername"""
    Pr = _Printer(verbose)

    Pr.write("\nComputing neighbour sets...")
    tr_g = Generations(ifs).compute_graph(depth=depth)
    if tr_g.is_wft:
        Pr.write(f"  IFS is weak finite type with {tr_g.g.num_vertices()} neighbour sets.\n")
    else:
        Pr.write(f"  Computed {tr_g.g.num_vertices()} neighbour sets without"
                " terminating; cannot guarantee weak finite type. Re-run with"
                " greater depth if desired.\n")
    draw = GraphArtist(tr_g)
    fl = f"{filename}.pdf"
    draw.gens(fl,scale=scale,depth=explore)
    Pr.write(f"- Net intervals saved to '{fl}'.")

    return tr_g

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
        Pr.write(f"  Computed {tr_g.g.num_vertices()} neighbour sets without"
                " terminating; cannot guarantee weak finite type. Re-run with"
                " greater depth if desired.\n")
    draw = GraphArtist(tr_g)
    write = GraphWriter(tr_g)

    # possibly draw generations
    if with_gens:
        gens_fname = f"{foldername}/gens.pdf"
        draw.gens(gens_fname,scale=scale)
        Pr.write(f"- Net intervals saved to '{gens_fname}'.")

    graph_fname = f"{foldername}/graph.pdf"
    info_fname = f"{foldername}/info.txt"

    draw.graph(graph_fname,edge_labels=edge_labels)
    write.info_to_file(info_fname)

    Pr.write(f"- Graph saved to '{graph_fname}'.")
    Pr.write(f"- Info saved to '{info_fname}'.")

    return tr_g

def verify_fnc(ifs,depth=1000):
    """Run the ifs, and determine if it is finite type or not."""
    ifs_gen = Generations(ifs)
    return ifs_gen.verify_wft(depth=depth)[0]
