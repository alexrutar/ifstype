from .gens_base import BaseGenerations, Gen
from .neighbour import FiniteNbMgr, InfiniteNbMgr
from .numeric import Constants as C

class FiniteType(BaseGenerations):
    def __init__(self, ifs, existing_nb_sets=None):
        super().__init__(ifs)

        self.nb_mgr = FiniteNbMgr(existing_nb_sets)
        # initialize the neighbour set
        to_update = [C.net_iv_base]
        gft = C.n_base
        while(len(to_update)>0):
            new = []
            for net_iv in to_update:
                new_nb = self.nb_set(net_iv)
                if new_nb not in self.nb_mgr:
                    self.nb_mgr.add(new_nb)
                    ch = self.im_children(net_iv)
                    gft = min(gft, ch.alpha)
                    new.extend(ch)
            to_update = new

        self.new_transition_stop = gft


class InfiniteType(BaseGenerations):
    def __init__(self, ifs, existing_nb_sets=None):
        super().__init__(ifs)

        self.nb_mgr = InfiniteNbMgr(existing_nb_sets)
        self.new_transition_stop = 0
