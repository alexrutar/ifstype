from .gens_base import BaseGenerations, Gen, GenKey
from .neighbour import FiniteNbMgr, InfiniteNbMgr
from .rational import Constants as C
from .interval import Interval,NetInterval, View
from .graph import TransitionGraph


class FiniteType(BaseGenerations):
    def __init__(self, ifs, existing_nb_sets=None):
        nb_mgr = FiniteNbMgr(existing_nb_sets)
        super().__init__(ifs, nb_mgr)

        # initialize the neighbour set
        to_update = [NetInterval(C.n_0,C.n_1,C.n_base)]
        gft = C.n_base
        # transition tree is a list of intervals along with the views required to observe all types
        # if you want to make a graph from a transition tree, just iterate through the transition tree and draw the types

        # initial transition generation
        self.transition_tree = [GenKey(C.n_base,View(Interval(C.n_0,C.n_1)))]

        while(len(to_update)>0):
            new = []
            for net_iv in to_update:
                new_nb = self.nb_set(net_iv)

                if new_nb not in self.nb_mgr:
                    transition = super().ttype(net_iv)
                    self.nb_mgr.add(new_nb, transition)
                    ch = self.im_children(net_iv) # recomputing is free because of hashing

                    # track this as a relevant area to draw in the tree
                    self.transition_tree.append(GenKey(ch.alpha, View(net_iv)))


                    gft = min(gft, ch.alpha)
                    new.extend(ch)
            to_update = new

        self.transition_graph = TransitionGraph(self.nb_mgr)
        self.new_transition_stop = gft

    def ttype(self, net_iv):
        return self.nb_mgr.transitions[self.nb_set(net_iv)]

class InfiniteType(BaseGenerations):
    def __init__(self, ifs, existing_nb_sets=None):
        super().__init__(ifs,InfiniteNbMgr(existing_nb_sets))
        self.new_transition_stop = 0

    # wrapper for the base method, since we might need to add new neighbour types to the manager
    def gen(self, *args, **kwargs):
        new_gen = super().gen(*args, **kwargs)
        self.nb_mgr.update(new_gen.nb_set(net_iv) for net_iv in new_gen)
        return new_gen
