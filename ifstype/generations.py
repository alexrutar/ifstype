""":mod:`ifstype.generations`
=============================

This module implements the :class:`ifstype.generations.Generations` which governs the computation of the transition graph from an IFS and other child dependency relationships.

Functions and classes defined in this module:

* :class:`Generations`

"""
import itertools
from collections import defaultdict
import numpy as np
import typing
import numbers

from .exact import Constants as C, Interval
from .ifs import AffineFunc, Neighbour, NeighbourSet, NetInterval, TransitionMatrix, IFS
from .graph import TransitionGraph, EdgeInfo


class Generations:
    """Compute children relationships and the transition graph based on a fixed iterated function system.

    Initialization:

    * :meth:`__init__`

    Methods to compute the transition graph:

    * :meth:`compute_graph`

    Methods to determine IFS properties:

    * :meth:`verify_wft`

    Methods to compute children:

    * :meth:`child_intervals`
    * :meth:`children`
    * :meth:`children_with_transition`

    """
    def __init__(self, ifs:IFS) -> None:
        """Initialize the Generations object.

        :param ifs: the iterated function system

        """
        self.ifs = ifs

    # main transition graph construction method
    def compute_graph(self, depth:int=2000, root:typing.Optional[NetInterval]=None) -> TransitionGraph:
        """Compute the transition graph based on the stored IFS.

        This method will compute the children for at minimum `depth` net intervals before stopping; if this method does not
        terminate earlier, the transition graph may be incomplete.

        To verify completeness of thetransition graph, see :attr:`ifstype.Graph.is_wft`

        See also :meth:`verify_wft`.

        :param depth: the minimum number of net intervals for which children have been computed before stopping.
        :param root: the root net interval from which the transition graph should be computed.
        :return: the transition graph

        """
        if root is None:
            root = NetInterval()
        transition_graph = TransitionGraph(root, self.ifs)

        transition_graph.add_nb_set(root.nb_set) # add root vertex
        computed_children = set() # set to check for nb_sets with children which have already been computed
        to_update = {root.nb_set}

        while(len(to_update)>0 and len(computed_children)<depth):
            all_children_params = list(itertools.chain.from_iterable([ch + (nb_set,) for ch in self.children_with_transition(nb_set)] for nb_set in to_update if nb_set not in computed_children))
            for ch_nb_set, e_info, parent_nb_set in all_children_params:
                transition_graph.add_nb_set(ch_nb_set,parent=parent_nb_set,edge_info=e_info)
            computed_children.update(cp[2] for cp in all_children_params)
            to_update = set(cp[0] for cp in all_children_params)

        # if IFS is weak finite type:
        if len(to_update) == 0:
            transition_graph.is_wft = True
            transition_graph.remove_terminal_vertices()

        return transition_graph

    def verify_wft(self, depth:int=2000) -> typing.Tuple[bool,int]:
        """Attempt to verify if the stored IFS is weak finite type.

        This method is a pared down version of :meth:`compute_graph`, where only the graph dependency structure is computed (without creating the transition graph).
        The method will compute the children for at minimum `depth` distinct neighbour sets.
        It returns a verification pair ``(bool, computed)``.

        See also :meth:`compute_graph`.

        * If ``bool`` is true, then the IFS is weak finite type with at most ``computed`` distinct neighbour sets (there may be fewer).
        * If ``bool`` is false, then the IFS may or may not be weak finite type, but there are at least ``computed`` distinct neighbour sets.

        :param depth: the minimum number of net intervals for which children have been computed before stopping.
        :return: the verification pair

        """
        computed_children = set() # nb_sets with children which have already been computed
        to_update = {NeighbourSet()} # nb_sets to compute children next iteration

        while(len(to_update)>0 and len(computed_children)<depth):
            # compute all children where the nb_set is not yet computed, and add them to the update list
            # TODO: maybe faster to do this with a list?
            new = set(itertools.chain.from_iterable(self.children(nb_set) for nb_set in to_update if nb_set not in computed_children))
            # add the update list to the list of computed children
            computed_children.update(to_update)
            to_update = new

        return len(to_update) == 0, len(computed_children)

    # -----------------------------------------------------------------
    # child computation methods
    # -----------------------------------------------------------------
    def child_intervals(self, new_nbs:typing.Iterable[Neighbour]) -> typing.Sequence[Interval]:
        """Compute all possible intervals formed by endpoints below a net interval [0,1] with from an iterable of neighbours `new_nbs`.

        :param new_nbs: an iterable of neighbours
        :return: a list of intervals in ascending order

        """
        new_eps = set(itertools.chain.from_iterable((f(C.n_0),f(C.n_1)) for f in new_nbs))
        child_endpoints = [C.n_0] + sorted(ep for ep in new_eps if 0 < ep and ep < 1) + [C.n_1]
        return [Interval(a,b) for a,b in zip(child_endpoints,child_endpoints[1:])]

    def children_with_transition(self, nb_set:NeighbourSet) -> typing.Sequence[typing.Tuple[NeighbourSet, EdgeInfo]]:
        """Compute the transition tuple for every child of the given neighbour set.

        This creates a sequence of :class:`ifstype.graph.EdgeInfo`, one for each child of the neighbour se.
        These quantities are invariant of the specific net interval which has given neighbour set and are used as edge attributes of the transition graph.
        See :meth:`compute_graph` for a usage case, and :meth:`children` for an analgous method without computing the edge info

        :param nb_set: the neighbour set for which to compute transition tuples
        :return: ordered sequence of pairs of the child neighbour set and the edge info, one for each child

        """
        max_nbs = nb_set.maximal_nbs() # maximal neighbours, i.e. with |L|=lmax
        # set elements are triples (p, orig_f, ext_f) where p is the probability associated with the function, orig_f is the initial, and ext_f is the extension
        non_max_nbs = set((C.n_1,nb,nb) for nb in nb_set.nonmaximal_nbs()) # non-maximal neighbours, which extend with probability 1
        new_nbs = set(self.ifs.extend(max_nbs,with_prob=True)) # extensions of the maximal neighbours, which also tracks the letter by which it was extended

        # children intervals, with endpoints that come from some new neighbour
        ch_ivls = self.child_intervals(nb[2] for nb in new_nbs)

        # transition_triples consist of lookup dictionaries for each child, where the key is the (orig_f, target_f) and the value is the probability of the extension
        # the target_neighbour is computed by normalizing the extension function against the new child interval
        transition_triples = [(ch_ivl,defaultdict(int,{(orig_f,Neighbour.from_aff(ext_f,ch_ivl)):p \
                for p,orig_f,ext_f in non_max_nbs.union(set(f for f in new_nbs if f[2].interval().supset(ch_ivl)))})) \
                for ch_ivl in ch_ivls]
        transition_triples = [(ch,triple) for ch,triple in transition_triples if len(triple) > 0]

        # get the new neighbour set by reading the target_neighbour in the transition triples
        normal_net_ivs = [NetInterval(ch_ivl.a,ch_ivl.b,nb_set.lmax,NeighbourSet(nbt[1] for nbt in triple.keys())) for ch_ivl,triple in transition_triples]
        # get the transitions by reading the probabilities from the transition_triples dictionaries
        transitions = [TransitionMatrix([[triple[(nb_par,nb_ch)] for nb_ch in net_iv.nb_set] for nb_par in nb_set]) for net_iv, (_,triple) in zip(normal_net_ivs, transition_triples)]

        return [(nt.nb_set, EdgeInfo(nt.a, nt.delta, tr)) for nt, tr in zip(normal_net_ivs, transitions)]


    def children(self, nb_set:NeighbourSet) -> typing.Sequence[NeighbourSet]:
        """Construct the neighbour sets which are the children of a net interval.

        See also :meth:`children_with_transition` if more information is desired

        :param nb_set: the neighbour set to compute the children of
        :return: sequence of neighbour sets

        """
        max_nbs = nb_set.maximal_nbs() # maximal neighbours, i.e. with |L|=lmax
        non_max_nbs = set(nb_set.nonmaximal_nbs())
        new_nbs = set(self.ifs.extend(max_nbs)) # extensions of the maximal neighbours
        holding = [(ch_ivl,tuple(itertools.chain(non_max_nbs,(f for f in new_nbs if f.interval().supset(ch_ivl))))) for ch_ivl in self.child_intervals(new_nbs)]
        return [NetInterval.from_funcs(ch_ivl.a,ch_ivl.b,nb_set.lmax,containing_funcs).nb_set for ch_ivl, containing_funcs in holding if len(containing_funcs)>0]
