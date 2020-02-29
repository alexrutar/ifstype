""":mod:`ifstype.algorithms`
=============================

This module implements various convenience functions which generate technical
information about IFSs as well as visualizations, and collects that information
in organized output files.

Public module attributes:

* :function:`run_ifs`
* :function:`run_ifs_gens`
* :function:`verify_fnc`

"""
from graph_tool.all import *
import pathlib

from .generations import Generations
from .info import GraphArtist, GraphWriter

class _Printer:
    """Generic printer class to display information when the printer is
    intialized in the verbose setting."""
    def __init__(self, verbose):
        self.verbose = verbose

    def write(self, string):
        if self.verbose:
            print(string)

    def in_prog(self):
        self.write("\nComputing neighbour sets...")

    def is_fnc(self,n):
        self.write("  IFS satisfies the finite neighbour condition with"
                   f" {n} neighbour sets.\n")
    def not_fnc(self,n):
        self.write(f"  Computed {n} neighbour sets without"
                   " terminating; cannot guarantee finite neighbour condition."
                   " Re-run with greater depth if desired.\n")
    def saved_gens(self,fname):
        self.write(f"- Net intervals saved to '{fname}'.")

    def saved_graph(self,fname):
        self.write(f"- Graph saved to '{fname}'.")

    def saved_info(self,fname):
        self.write(f"- Info saved to '{fname}'.")

def run_ifs_gens(ifs,filename,explore,scale='wide',verbose=True,depth=2000):
    """Compute variations on the "generations" diagram for the given IFS.
    Run the ifs and place output files in foldername"""
    Pr = _Printer(verbose)

    Pr.in_prog()
    tr_g = Generations(ifs).compute_graph(depth=depth)
    num_nbsets_computed = tr_g.g.num_vertices()
    if tr_g.is_fnc:
        Pr.is_fnc(num_nbsets_computed)
    else:
        Pr.not_fnc(num_nbsets_computed)

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

    Pr.in_prog()
    tr_g = Generations(ifs).compute_graph(depth=depth)
    num_nbsets_computed = tr_g.g.num_vertices()
    if tr_g.is_fnc:
        Pr.is_fnc(num_nbsets_computed)
    else:
        Pr.not_fnc(num_nbsets_computed)

    draw = GraphArtist(tr_g)
    write = GraphWriter(tr_g)

    # possibly draw generations
    if with_gens:
        gens_fname = f"{foldername}/gens.pdf"
        draw.gens(gens_fname,scale=scale)
        Pr.saved_gens(gens_fname)

    graph_fname = f"{foldername}/graph.pdf"
    info_fname = f"{foldername}/info.txt"

    draw.graph(graph_fname,edge_labels=edge_labels)
    write.info_to_file(info_fname)

    Pr.saved_graph(graph_fname)
    Pr.saved_info(info_fname)

    return tr_g

def verify_fnc(ifs,depth=2000):
    """Run the ifs, and determine if it is finite type or not."""
    ifs_gen = Generations(ifs)
    return ifs_gen.verify_fnc(depth=depth)[0]
