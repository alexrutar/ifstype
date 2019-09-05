import itertools
from collections import defaultdict
import numpy as np

from .exact import Constants as C, Interval
from .ifs import AffineFunc, Neighbour, NeighbourSet, NetInterval, TransitionMatrix
from .graph import TransitionGraph


class Generations:
    """
    Constructor for iterated function systems, based on a fixed IFS.
    Has methods to compute children and neighbour sets based on parents.
    """
    def __init__(self, ifs):
        self.ifs = ifs

    # main transition graph construction method
    def compute_graph(self, depth=2000, root=None):
        """Arguments:
        - depth (default 2000): the minimum number of net intervals for which children have been computed before stopping.
        - root: the root net interval from which the transition graph should be computed.

        Returns:
        - tr_graph: a TransitionGraph object
        """
        if root is None:
            root = NetInterval.base()
        transition_graph = TransitionGraph(root, self.ifs)

        transition_graph.add_nb_set(root.nb_set) # add root vertex
        computed_children = set() # set to check for nb_sets with children which have already been computed
        to_update = {root.nb_set}

        while(len(to_update)>0 and len(computed_children)<depth):
            all_children_params = list(itertools.chain.from_iterable([ch+(nb_set,) for ch in self.children(nb_set,with_transition=True)] for nb_set in to_update if nb_set not in computed_children))
            for ch_nb_set, t_index, length, transition, parent_nb_set in all_children_params:
                transition_graph.add_nb_set(ch_nb_set,parent=parent_nb_set,transition=transition,length=length, t_index=t_index)
            computed_children.update(cp[4] for cp in all_children_params)
            # TODO: maybe faster with list here
            to_update = set(cp[0] for cp in all_children_params)

        # if IFS is weak finite type:
        if len(to_update) == 0:
            transition_graph.is_wft = True
            transition_graph.remove_terminal_vertices()

        return transition_graph

    def verify_wft(self, depth=2000):
        """A pared down version of compute_graph, where we only verify if the IFS is weak finite type (without creating transition graph).
        - returns (bool, depth):
        -   bool is if WFT guarantee
        -   depth is number of neighbour set children computed (effective depth): this is lower bound if infinite, and upper bound if finite (not exact in either case).
        Does not save any information within the ifs_gens object
        """
        computed_children = set() # set to check for nb_sets with children which have already been computed
        to_update = {NeighbourSet.base()}

        while(len(to_update)>0 and len(computed_children)<depth):
            # compute all children where the nb_set is not yet computed, and add them to the update list
            # TODO: maybe faster to do this with a list?
            new = set(itertools.chain.from_iterable(self.children(nb_set,with_transition=False) for nb_set in to_update if nb_set not in computed_children))
            # add the update list to the list of computed children
            computed_children.update(to_update)
            to_update = new

        return len(to_update) == 0, len(computed_children)

    # -----------------------------------------------------------------
    # child computation methods
    # -----------------------------------------------------------------
    def child_intervals(self, new_nbs):
        new_eps = set(itertools.chain.from_iterable(f.endpoints() for f in new_nbs))
        child_endpoints = [C.n_0] + sorted(ep for ep in new_eps if 0 < ep and ep < 1) + [C.n_1]
        return (Interval.closed(a,b) for a,b in zip(child_endpoints,child_endpoints[1:]))

    def nb_children(self, nb_set,with_transition=True):
        """Construct the un-normalized children of the nb_set."""
        max_nbs = nb_set.maximal_nbs() # maximal neighbours, i.e. with |L|=lmax
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
            holding = [(ch_ivl,tuple(itertools.chain(non_max_nbs,(f for f in new_nbs if f.interval().supset(ch_ivl))))) for ch_ivl in self.child_intervals(new_nbs)]
            return [NetInterval.from_funcs(ch_ivl.a,ch_ivl.b,nb_set.lmax,containing_funcs).nb_set for ch_ivl, containing_funcs in holding if len(containing_funcs)>0]


    def children(self, nb_set, with_transition=True):
        """Construct the children of a net interval"""
        if with_transition:
            normal_net_ivs, transitions = self.nb_children(nb_set,with_transition=True)
            # (nb_set, t_index, length, transition) in order for each child
            return [(nt.nb_set, nt.a, nt.delta, tr) for nt, tr in zip(normal_net_ivs, transitions)]
        else:
            return self.nb_children(nb_set,with_transition=False)
