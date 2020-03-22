""":mod:`ifstype.generations`
=============================

This module implements the :class:`ifstype.generations.Generations` class which
governs the computation of the transition graph from an IFS and other child
dependency relationships.

Public module attributes:

* :class:`Generations`
* :class:`LengthFunction`

"""
import itertools
from collections import defaultdict
import numpy as np
from typing import Sequence, Iterable, Optional, Tuple
from numbers import Real

from .exact import Constants as C, Interval

from .ifs import (
    AffineFunc, Neighbour, NeighbourSet, NetInterval,
    TransitionMatrix, IFS
)

from .graph import TransitionGraph, EdgeInfo


class LengthFunction:
    """Class with existing length functions for general usage.

    Implemented length functions:

    * :meth:`measure`
    * :meth:`transition`

    Note that each method must take as arguments two parameters, where the
    first is the net interval forresponding to the child and the second is the
    neighbour set of the parent.
    """
    @staticmethod
    def measure(nt:NetInterval, nb_set:NeighbourSet) -> Real:
        """Measure length function.

        :param nt: normalized child
        :param nb_set: parent neighbour set

        """
        return nt.delta
    @staticmethod
    def transition(nt:NetInterval, nb_set:NeighbourSet) -> Real:
        """Transition length function.

        :param nt: normalized child
        :param nb_set: parent neighbour set

        """
        return nt.delta*nt.nb_set.lmax/nb_set.lmax

class Generations:
    """Compute children relationships and the transition graph based on a fixed
    iterated function system.

    Initialization:

    * :meth:`__init__`

    Methods to compute the transition graph:

    * :meth:`compute_graph`

    Methods to determine IFS properties:

    * :meth:`verify_fnc`

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
    def compute_graph(
            self,
            depth: int=2000,
            root: Optional[NetInterval]=None,
            non_red_nbs: Optional[set] = None
        ) -> TransitionGraph:
        """Compute the transition graph based on the stored IFS.

        This method will compute the children for at minimum `depth` net
        intervals before stopping; if this method does not terminate earlier,
        the transition graph may be incomplete.

        To verify completeness of thetransition graph, see
        :attr:`ifstype.Graph.is_fnc`

        See also :meth:`verify_fnc`.

        :param depth: the minimum number of net intervals for which children
                      have been computed before stopping.
        :param root: the root net interval from which the transition graph
                     should be computed.
        :param non_red_nbs: a (previously computed) set of reduced
                                   neighbours from which the neighbours
                                   most be drawn
        :return: the transition graph

        """
        if root is None:
            root = NetInterval()
        transition_graph = TransitionGraph(root, self.ifs)

        transition_graph.add_nb_set(root.nb_set)
        computed_children = set()
        ch_to_compute = {root.nb_set}

        while(len(ch_to_compute)>0 and len(computed_children)<depth):
            all_children_params = list(itertools.chain.from_iterable(
                [
                    ch + (nb_set,) for ch
                    in self.children_with_transition(
                        nb_set,
                        non_red_nbs=non_red_nbs)]
                    for nb_set
                in ch_to_compute
                if nb_set not in computed_children))

            for ch_nb_set, e_info, parent_nb_set in all_children_params:
                transition_graph.add_nb_set(
                    ch_nb_set,
                    parent=parent_nb_set,
                    edge_info=e_info)

            computed_children.update(cp[2] for cp in all_children_params)
            ch_to_compute = set(cp[0] for cp in all_children_params)

        # if IFS is weak finite type:
        if len(ch_to_compute) == 0:
            transition_graph.is_fnc = True
            if non_red_nbs is not None:
                transition_graph.collapse()
            else:
                transition_graph.remove_terminal_vertices()

        return transition_graph

    def verify_fnc(self, depth:int=2000) -> Tuple[bool,int]:
        """Attempt to verify if the stored IFS satisfies the finite neighbour
        condition.

        This method is a pared down version of :meth:`compute_graph`, where
        only the graph dependency structure is computed (without creating the
        transition graph). The method will compute the children for at minimum
        `depth` distinct neighbour sets. It returns a verification pair
        ``(bool, computed)``.

        See also :meth:`compute_graph`.

        * If ``bool`` is true, then the IFS satisfies the finite neighbour
          condition with at most ``computed`` distinct neighbour sets
          (there may be fewer).
        * If ``bool`` is false, then the IFS may or may not be weak finite
          type, but there are at least ``computed`` distinct neighbour sets.

        :param depth: the minimum number of net intervals for which children
                      have been computed before stopping.
        :return: the verification pair

        """
        computed_children = set()
        ch_to_compute = {NeighbourSet()}

        while(len(ch_to_compute)>0 and len(computed_children)<depth):
            # TODO: maybe faster to do this with a list?
            new = set(
                itertools.chain.from_iterable(self.children(nb_set) for nb_set
                in ch_to_compute if nb_set not in computed_children))
            computed_children.update(ch_to_compute)
            ch_to_compute = new

        return len(to_update) == 0, len(computed_children)

    # -----------------------------------------------------------------
    # child computation methods
    # -----------------------------------------------------------------
    def child_intervals(
            self,
            new_nbs: Iterable[Neighbour]
        ) -> Sequence[Interval]:
        """Compute all possible normalized intervals formed by endpoints below
        a net interval [0,1] with from an iterable of neighbours `new_nbs`.

        :param new_nbs: an iterable of neighbours
        :return: a list of intervals in ascending order

        """
        new_eps = set(
            itertools.chain.from_iterable(
                (f(C.n_0),f(C.n_1)) for f in new_nbs))
        child_endpoints = [C.n_0] \
                + sorted(ep for ep in new_eps if 0 < ep and ep < 1) \
                + [C.n_1]
        return [
            Interval(a,b) for a,b
            in zip(child_endpoints,child_endpoints[1:])]

    def children_with_transition(
            self,
            nb_set: NeighbourSet,
            non_red_nbs=None,
            length_func=LengthFunction.transition
            # TODO: actually implement "length functions", think about proper
            # arguments
        ) -> Sequence[Tuple[NeighbourSet, EdgeInfo]]:
        """Compute the transition tuple for every child of the given neighbour
        set.

        This creates a sequence of :class:`ifstype.graph.EdgeInfo`, one for
        each child of the neighbour set. These quantities are invariant of the
        specific net interval which has given neighbour set and are used as
        edge attributes of the transition graph.
        See :meth:`compute_graph` for a usage case, and :meth:`children` for an
        analgous method without computing the edge info

        :param nb_set: the neighbour set for which to compute transition tuples
        :param length_func: the length function used to assign a length to each
                            edge
        :return: ordered sequence of pairs of the child neighbour set and the
                 edge info, one for each child

        """
        max_nbs = nb_set.maximal_nbs()
        # since we want to compute the transition matrix, we have triples
        # (p, orig_f, ext_f) where
        #   - p is the probability associated with the extension function
        #   - orig_f is the initial function
        #   - ext_f is the extension function
        non_max_nbs = set((C.n_1,nb,nb) for nb in nb_set.nonmaximal_nbs())
        new_nbs = self.ifs.extend(max_nbs,with_prob=True)

        ch_ivls = self.child_intervals(nb[2] for nb in new_nbs)

        # To compute the transition matrix, transition_pairs consists of the
        # child interval and a corresponding defaultdicts (with default 0) for
        # each child, where the key is the (orig_f, nb)f) and the value is the
        # probability of the extension
        transition_pairs = [
            (ch_ivl, defaultdict(
                int,
                {
                    (orig_f,Neighbour.from_aff(ext_f,ch_ivl)):p
                    for p,orig_f,ext_f
                    in non_max_nbs.union(set(
                        f for f in new_nbs
                        if f[2].interval().supset(ch_ivl)))}))
            for ch_ivl in ch_ivls]

        # remove neighbours which are not reduced
        if non_red_nbs is not None:
            transition_pairs = [
                (ch, defaultdict(
                    int,
                    {k:p for k,p in lookup.items()
                        if k[1] not in non_red_nbs}))
                for ch,lookup in transition_pairs]

        # remove net intervals with no neighbours
        transition_pairs = [
            (ch,lookup) for ch,lookup
            in transition_pairs
            if len(lookup) > 0]

        ch_net_ivs = [
            NetInterval(
                ch_ivl.a,
                ch_ivl.b,
                nb_set.lmax,
                NeighbourSet(nbt[1] for nbt in lookup.keys()))
            for ch_ivl,lookup in transition_pairs]

        transitions = [
                TransitionMatrix(
                    [
                        [lookup[(nb_par,nb_ch)] for nb_ch
                            in net_iv.nb_set.sorted_iter()]
                        for nb_par in nb_set.sorted_iter()])
                for net_iv, (_,lookup) in zip(ch_net_ivs, transition_pairs)]

        return [(
                nt.nb_set,
                EdgeInfo(nt.a, nt.delta, length_func(nt,nb_set), tr))
            for nt, tr in zip(ch_net_ivs, transitions)]


    def children(self, nb_set:NeighbourSet) -> Sequence[NeighbourSet]:
        """Construct the neighbour sets which are the children of a net
        interval.

        See also :meth:`children_with_transition`.

        :param nb_set: the neighbour set to compute the children of
        :return: sequence of neighbour sets

        """
        max_nbs = nb_set.maximal_nbs()
        non_max_nbs = set(nb_set.nonmaximal_nbs())
        new_nbs = self.ifs.extend(max_nbs)
        holding = [
            (ch_ivl, tuple(
                itertools.chain(
                    non_max_nbs,
                    (f for f in new_nbs if f.interval().supset(ch_ivl)))))
            for ch_ivl in self.child_intervals(new_nbs)]
        return [
            NetInterval.from_funcs(
                ch_ivl.a,
                ch_ivl.b,
                nb_set.lmax,
                containing_funcs).nb_set for ch_ivl, containing_funcs
            in holding if len(containing_funcs)>0]




