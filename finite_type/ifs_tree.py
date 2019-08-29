import itertools
from sortedcontainers import SortedSet
from collections import defaultdict
import numpy as np

from .numerics.rational import Constants as C
from .numerics.interval import Interval

from .transition import Neighbour, NeighbourSet, TransitionGraph, TransitionMatrix
from .ifs import CtrFunc
from .draw import IFSDrawDispatch


class NetInterval(Interval):
    """A special interval type representing a net interval of generation alpha.
    This contains the interval information, as well as the neighbour set."""

    def __new__(cls,a,b,alpha,nb_set):
        """NetInterval class, based on interval."""
        self = super().__new__(cls,a,b)
        self.nb_set = nb_set
        self.alpha = alpha
        return self

    def __hash__(self):
        # net intervals should also distinguish on alpha
        return hash((self.a,self.b,self.alpha))

    def transition_gen(self):
        return self.nb_set.lmax*self.delta

    def normalization_func(self):
        return CtrFunc(self.delta,self.a)

    def containing_funcs(self):
        return (self.normalization_func().compose(f) for f in self.nb_set)

    @classmethod
    def base(cls):
        return cls(C.n_0,C.n_1,C.n_base,NeighbourSet.base())

    @classmethod
    def from_funcs(cls,a,b,alpha,funcs):
        self = super().__new__(cls,a,b)
        self.nb_set = NeighbourSet(Neighbour.from_f(f,self) for f in funcs)
        self.alpha = alpha
        return self

    # representation
    def __str__(self):
        return f"NetIv({self.alpha})[{self.a},{self.b}]"
    def __repr__(self):
        return f"NetInterval(left={self.left},right={self.right},alpha={self.alpha},nb_set={self.nb_set})"

class Generations:
    """Creates a graph representation in which the nodes are neighbour sets, with possibly corresponding net intervals.
    Continually explores along edges that are not yet in the graph, up to specified depth."""
    def __init__(self, ifs, depth=None,with_transition=True):
        self.ifs = ifs
        self.tr_graph = TransitionGraph()
        self.with_transition = with_transition
        self.net_ivs = set()
        self.draw_dispatch = IFSDrawDispatch(self)

    # main transition graph construction method
    def compute_graph(self, depth=2000, start=None, compute_net_iv=False, transition_graph=None):
        """If depth is none, explores until fully explored.
        Arguments:
        - depth (default 2000): the minimum number of net intervals for which children have been computed before stopping
        - start (default [0,1]): the net interval in which to start computing the graph
        - compute_net_iv (default False): return the net intervals computed when determining the graph
        - destination (default self.tr_graph): a TransitionGraph instance to use when computing the graph

        Returns:
        - net_ivs: a (possibly empty) iterable of net intervals
        """
        print("\nComputing neighbour sets...")
        if transition_graph is None:
            transition_graph = self.tr_graph

        if start is None:
            start = NetInterval.base()

        transition_graph.add_nb_set(start.nb_set) # add root vertex
        computed_children = set() # set to check for nb_sets with children which have already been computed

        net_iv_save = defaultdict(set)
        to_update = [start]
        while(len(to_update)>0 and len(computed_children)<depth):
            print(f"  finished: {len(computed_children)}")
            new = []
            for net_iv in to_update:
                if compute_net_iv:
                    net_iv_save[net_iv.nb_set].add(net_iv)
                if net_iv.nb_set not in computed_children:
                    children, t_indices, lengths,  transitions = self.children(net_iv,with_transition=self.with_transition,normalized=not compute_net_iv)
                    new.extend(children)
                    for ch, t_index, length, transition in zip(children, t_indices, lengths, transitions):
                        # transitions has wrong length
                        transition_graph.add_nb_set(ch.nb_set,prec=net_iv.nb_set,transition=transition,length=length, t_index=t_index)
                    computed_children.add(net_iv.nb_set)


            to_update = new

        if len(to_update) == 0:
            transition_graph.is_wft = True
            zero_degrees = transition_graph.remove_terminal_vertices()

            print("\nIFS is weak finite type.")

            # and remove them from the saved net intervals
            if compute_net_iv:
                for zd_neighbour in zero_degrees:
                    del net_iv_save[zd_neighbour]
                

        else:
            print(f"\nIFS hit depth {depth} with uncomputed neighbour sets.\n  Re-run `compute_graph` with option depth=N for greater depth.")
        return set(itertools.chain.from_iterable(net_iv_save.values()))


    # -----------------------------------------------------------------
    # child computation methods
    # -----------------------------------------------------------------
    def child_intervals(self, new_nbs):
        new_eps = set(itertools.chain.from_iterable(f.endpoints() for f in new_nbs))
        child_endpoints = [C.n_0] + sorted(ep for ep in new_eps if 0 < ep and ep < 1) + [C.n_1]
        return (Interval.closed(a,b) for a,b in zip(child_endpoints,child_endpoints[1:]))

    def nb_children(self, nb_set,with_transition=True):
        """Construct the un-normalized children of the nb_set."""
        # maximal neighbours, i.e. with |L|=lmax
        max_nbs = nb_set.maximal_nbs()
        if with_transition:
            # set elements are triples (p, orig_f, ext_f) where p is the probability associated with the function, orig_f is the initial, and ext_f is the extension
            non_max_nbs = set((C.n_1,nb,nb) for nb in nb_set.nonmaximal_nbs()) # non-maximal neighbours, which extend with probability 1
            new_nbs = set(self.ifs.extend_with_prb(max_nbs)) # extensions of the maximal neighbours, which also tracks the letter by which it was extended

            # children intervals, with endpoints that come from some new neighbour
            ch_ivls = self.child_intervals(nb[2] for nb in new_nbs)

            # transition_triples consist of lookup dictionaries for each child, where the key is the (orig_f, target_f) and the value is the probability of the extension
            # the target_neighbour is computed by normalizing the extension function against the new child interval
            transition_triples = [(ch_ivl,defaultdict(int,{(orig_f,Neighbour.from_f(ext_f,ch_ivl)):p for p,orig_f,ext_f in non_max_nbs.union(set(f for f in new_nbs if f[2].interval().supset(ch_ivl)))})) for ch_ivl in ch_ivls]
            transition_triples = [(ch,triple) for ch,triple in transition_triples if len(triple) > 0]

            # get the new neighbour set by reading the target_neighbour in the transition triples
            net_intervals = [NetInterval(ch_ivl.a,ch_ivl.b,nb_set.lmax,NeighbourSet(nbt[1] for nbt in triple.keys())) for ch_ivl,triple in transition_triples]
            # get the transitions by reading the probabilities from the transition_triples dictionaries
            transitions = [TransitionMatrix([[triple[(nb_par,nb_ch)] for nb_ch in net_iv.nb_set] for nb_par in nb_set]) for net_iv, (_,triple) in zip(net_intervals, transition_triples)]
            return net_intervals, transitions

        else:
            non_max_nbs = set(nb_set.nonmaximal_nbs())
            new_nbs = set(self.ifs.extend(max_nbs)) # extensions of the maximal neighbours
            return [NetInterval.from_funcs(ch_ivl.a,ch_ivl.b,nb_set.lmax,itertools.chain(non_max_nbs,(f for f in new_nbs if f.interval().supset(ch_ivl)))) for ch_ivl in self.child_intervals(new_nbs)], None


    def children(self, net_iv, with_transition=False, normalized = False):
        """Construct the children of a net interval"""
        normal_net_ivs, transitions = self.nb_children(net_iv.nb_set,with_transition=with_transition)
        lengths = [nt.delta for nt in normal_net_ivs]

        if normalized:
            net_intervals = [NetInterval(C.n_0,C.n_1, C.n_base, nt.nb_set) for nt in normal_net_ivs]
        else:
            net_intervals = [NetInterval(net_iv.a+nt.a*net_iv.delta,net_iv.a+nt.b*net_iv.delta,nt.alpha*net_iv.delta, nt.nb_set) for nt in normal_net_ivs]
        t_indices = [nt.a for nt in normal_net_ivs]
        if with_transition:
            return net_intervals, t_indices, lengths, transitions
        else:
            return net_intervals, t_indices, lengths, None


    # -----------------------------------------------------------------
    # visualization / information methods
    # -----------------------------------------------------------------
    def draw_gens(self,net_iv_save,filename,**kwargs):
        print("\nDrawing interval levels...")
        self.draw_dispatch.draw_gens(net_iv_save,filename,**kwargs)
        print(f"  Done! Levels saved to '{filename}'.")

    def draw_graph(self,filename,**kwargs):
        print("\nDrawing transition graph...")
        self.draw_dispatch.draw_graph(filename,**kwargs)
        print(f"  Done! Graph saved to '{filename}'.")

    def write_info(self,filename):
        # write the neighbour sets
        print("\nWriting graph info to file...")
        with open(filename,'w') as outfile:
            outfile.write(str(self.ifs))
            outfile.write("\n\nNeighbour Sets:\n")
            outfile.write("\n".join(f"{self.tr_graph.get_identifier(nb_set)} : {nb_set}" for nb_set in self.tr_graph.all_nb_sets()))
            outfile.write("\n\nEdge Information:\n")
            for e in self.tr_graph.g.edges():
                outfile.write("\n")
                outfile.write(self.tr_graph.edge_info_string(e))
        print(f"  Done! Info saved to '{filename}'.")



