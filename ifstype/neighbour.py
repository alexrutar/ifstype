import typing
from sortedcontainers import SortedSet
from sympy import Rational

from .interval import NetInterval
from .ifs import C

class Neighbour(typing.NamedTuple):
    "A Neighbour object, with operations defined with respect to a net_iv"
    a: Rational
    L: Rational

    @classmethod
    def from_f(cls,f,net_iv):
        delta = net_iv.delta()
        return cls(
                (net_iv.a-f(C.n_0))/delta,
                f.r/delta)

    def to_f(self,net_iv):
        delta = net_iv.delta()
        return CtrFunc(self.L*delta,delta*self.a+net_iv.a)

class NeighbourSet(tuple):
    def __new__(self,neighbours):
        return super().__new__(self,sorted(neighbours))

    def __str__(self):
        return str(tuple((nb.a,nb.L) for nb in self))


# a NeighbourManager is any class which supports the functions
# - 
class InfiniteNbMgr:
    def __init__(self, existing_nb_sets=None):
        if existing_nb_sets is None:
            self.nb_set_types = {NeighbourSet((Neighbour(C.n_0,C.n_1),)):0}
            self.num_nb_set = 1
        else:
            self.nb_set_types = {}
            self.update(existing_nb_sets)
        self.is_complete = False

    def __str__(self):
        return "\n".join(f"{v} : {k}" for k,v in self)

    def __iter__(self):
        return iter(self.nb_set_types.keys())

    def nb_set_type(self, nb):
        "Return a unique identification string for the corresponding neighbour type"
        return self.nb_set_types[nb]

    def update(self, nb_iter):
        for nb in set(nb_iter):
            self.add(nb)

    def add(self, nb):
        if nb not in self.nb_set_types.keys():
            self.nb_set_types[nb] = self.num_nb_set
            self.num_nb_set += 1

class FiniteNbMgr:
    def __init__(self,existing_nb_sets=None):
        """Takes a generation object as an argument in order to compute all the neighbour types.
        """
        self.is_complete = True
        if existing_nb_sets is None:
            self.nb_set_types = SortedSet([])
        else:
            self.nb_set_types = SortedSet(existing_nb_sets)

    def create_nb_set(self,gn):
        to_update = [NetInterval(C.n_1,C.n_0,C.n_1)]
        while(len(to_update)>0):
            new = []
            for net_iv in to_update:
                new_nb = gn.nb_set(net_iv)
                if new_nb not in self.nb_set_types:
                    self.nb_set_types.add(new_nb)
                    new.extend(gn.im_children(net_iv))
            to_update = new

    def __iter__(self):
        return iter(self.nb_set_types)

    def __str__(self):
        return "\n".join(f"{v} : {k}" for k,v in self)

    def nb_set_type(self, nb):
        return self.nb_set_types.index(nb)

    # should never need to add or update new elements!
    def update(self, nb_iter):
        pass

    def add(self,nb):
        self.nb_set_types.add(nb)
