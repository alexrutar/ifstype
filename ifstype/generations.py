from .gens_base import BaseGenerations, Gen
from .neighbour import FiniteNbMgr, InfiniteNbMgr
from .rational import Constants as C
from .interval import NetInterval
from .graph import TransitionGraph


class FiniteType(BaseGenerations):
    def __init__(self, ifs, existing_nb_sets=None):
        nb_mgr = FiniteNbMgr(existing_nb_sets)
        super().__init__(ifs, nb_mgr)
        self.transition_graph = TransitionGraph(nb_mgr)

        # initialize the neighbour set
        to_update = [NetInterval(C.n_0,C.n_1,C.n_base)]
        gft = C.n_base

        while(len(to_update)>0):
            new = []
            for net_iv in to_update:
                new_nb = self.nb_set(net_iv)
                if new_nb not in self.nb_mgr:
                    self.nb_mgr.add(new_nb)
                    ch = self.im_children(net_iv)
                    self.transition_graph.add(new_nb, (self.nb_set(nv) for nv in ch))
                    gft = min(gft, ch.alpha)
                    new.extend(ch)
            to_update = new

        self.transition_graph.set_labels(self.nb_mgr.nb_set_type)
        self.new_transition_stop = gft


class InfiniteType(BaseGenerations):
    def __init__(self, ifs, existing_nb_sets=None):
        super().__init__(ifs,InfiniteNbMgr(existing_nb_sets))
        self.new_transition_stop = 0

    # wrapper for the base method, since we might need to add new neighbour types to the manager
    def gen(self, *args, **kwargs):
        new_gen = super().gen(*args, **kwargs)
        self.nb_mgr.update(new_gen.nb_set(net_iv) for net_iv in new_gen)
        return new_gen
