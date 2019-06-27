import typing
from sortedcontainers import SortedSet
from .rational import Rational, Constants as C

class Neighbour(typing.NamedTuple):
    "A Neighbour object, with operations defined with respect to a net_iv"
    a: Rational
    L: Rational

    @classmethod
    def from_f(cls,f,net_iv):
        delta = net_iv.delta
        return cls(
                (net_iv.a-f(C.n_0))/delta,
                f.r/delta)

    def to_f(self,net_iv):
        delta = net_iv.delta
        return CtrFunc(self.L*delta,delta*self.a+net_iv.a)

class NeighbourSet(tuple):
    def __new__(self,neighbours):
        return super().__new__(self,sorted(neighbours))

    def __str__(self):
        return ", ".join(f"({nb.a},{nb.L})" for nb in self)


# a NeighbourManager is any class which supports the functions
# - 
class InfiniteNbMgr:
    def __init__(self):
        self.nb_set_types = {NeighbourSet((Neighbour(C.n_0,C.n_1),)):0}
        self.num_nb_set = 1
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
        self.nb_set_types = SortedSet([])

        # transition / matrices are dictionaries
        self.transitions={}
        self.transition_matrix={}

    def __iter__(self):
        return iter(self.nb_set_types)

    def __str__(self):
        return "\n".join(f"{self.nb_set_type(nb)} : {nb}" for nb in self)

    def num_nb_sets(self):
        return len(self.nb_set_types)

    def nb_set_type(self, nb):
        return self.nb_set_types.index(nb)

    def update(self, nb_iter):
        raise NotImplementedError

    def add(self,nb,ttype):
        self.nb_set_types.add(nb)
        self.transitions[nb] = ttype

class TransitionType:
    # computes and stores the transition type and transition matrix to each subtype
    def __init__(self, transition_gen):
        self.ttype = tuple(transition_gen)

    def shifts(self):
        return (tty[0] for tty in self.ttype)

    def deltas(self):
        return (tty[1] for tty in self.ttype)

    def nb_sets(self):
        return (tty[2] for tty in self.ttype)
