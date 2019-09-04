"""
Core functions and classes concerning the technical operation of Iterated Function Systems.
"""
import itertools
import numpy as np
import typing
import inspect

from .exact import Rational, Constants as C, Interval
from .exact.symbolic import SymbolicRing, SymbolicMatrix

class CtrFunc(typing.NamedTuple):
    """A CtrFunc class on parameters (r,d) represents a contraction function f(x)=r*x+d.
    """
    r: Rational
    d: Rational

    @classmethod
    def id(cls):
        """The identity function."""
        return cls(C.n_1,C.n_0)

    def sign(self):
        return self.r>0

    def inverse(self):
        """Compute the inverse function of the linear map."""
        return CtrFunc(1/self.r,-self.d/self.r)

    def compose(self, ct_f):
        """Compute the right composition self(ct_f(x))."""
        return CtrFunc(self.r*ct_f.r, self.d+self.r*ct_f.d)

    def fixed_point(self):
        """Compute the point which the contraction function fixes."""
        return self.d/(C.n_1-self.r)

    def endpoints(self):
        return (self.d,self.d+self.r)

    def interval(self,iv=Interval.closed(C.n_0,C.n_1)):
        """Compute the interval self(iv) corresponding to the function."""
        left = self(iv.a)
        right = self(iv.b)
        if left <= right:
            return Interval(a=left,b=right,has_left=iv.has_left,has_right=iv.has_right)
        else:
            return Interval(a=right,b=left,has_left=iv.has_right,has_right=iv.has_left)

    def normalize(self, interval):
        "Normalize with respect to the given interval."
        cur_iv = self.interval(iv=interval)
        new_iv = (cur_iv - interval.a)/interval.delta
        if self.r > C.n_0:
            return CtrFunc(new_iv.delta, new_iv.a)
        else:
            return CtrFunc(-new_iv.delta, new_iv.b)

    def __call__(self, x):
        return self.r*x+self.d
    def __repr__(self):
        return f"CtrFunc(r={self.r},d={self.d})"
    def __str__(self):
        return f"f:x*({self.r}) + ({self.d})"

def ifs_family(fn):
    """Decorator to construct families of iterated function systems parametrized by some set of functions.
    When decorated with @ifs_family, call functions by first specifying probs if necessary, followed by keyword arguments.
    """
    def wrapper(**kwargs):
        fn_kwargs = {k:v.default for k,v in inspect.signature(fn).parameters.items()}
        func_params = {**fn_kwargs,**kwargs}
        try:
            probs = func_params['probs']
        except KeyError:
            raise KeyError("Function decorated by 'ifs_family' has no default keyword 'probs'.")
        if probs is None:
            return IFS.uniform_p(*fn(**func_params))
        else:
            return IFS(fn(**func_params),probs)
    return wrapper


class IFS:
    # an IFS is essentially a factory for words and ctr funcs
    def __init__(self, funcs, probabilities=None):
        """funcs is a list of functions"""
        # check params
        assert all(C.n_0<abs(f.r)<C.n_1 for f in funcs), "IFS contraction factors must have 0<|r|<1"
        assert all(C.n_0<=p for p in probabilities) and sum(probabilities) == C.n_1, "IFS probabilities must be non-negative and sum to 1"

        self.syr = SymbolicRing((f"p{i+1}" for i in range(len(funcs))))
        self.p = self.syr.get_symbols() # use symbolic probabilities

        if probabilities is not None:
            sorted_f_pairs = sorted(zip(funcs,probabilities),key=lambda x:(x[0].d,x[0].r))
            self.f = tuple(f[0] for f in sorted_f_pairs)
            probs = tuple(f[1] for f in sorted_f_pairs)
            self.set_probs(probs)
        else:
            self.f = tuple(sorted(funcs,key=lambda x:(x[0].d,x[0].r)))

        self.r = tuple(f.r for f in self.f)
        self.d = tuple(f.d for f in self.f)
        self.normalize()

    def set_probs(self,probs):
        assert len(self.f) == len(probs), "funcs and probabilities must be the same size"
        self.syr.set_eval({f"p{i+1}":p for i,p in enumerate(probs)})

    def __str__(self):
        return f"IFS(funcs={self.f})"

    @classmethod
    def uniform_p(cls, *funcs):
        return cls(funcs,[Rational(1,len(funcs)) for _ in funcs])

    def normalize(self):
        cvx_hull = self.convex_hull()
        self.f = [f.normalize(cvx_hull) for f in self.f]
        self.d = [f.d for f in self.f]

    def convex_hull(self):
        """Convex Hull computation, as adapted from József Vass' paper, which can be found at https://arxiv.org/abs/1502.03788 Section 3.2."""
        def value(tup):
            return sum(1 for a in tup if self.r[a] < 0) % 2

        def ct_from_tuple(tup):
            "Create the contraction function associated to the tuple argument"
            f = CtrFunc.id()
            for letter in tup:
                f = f.compose(self.f[letter])
            return f

        n=2 if all(r>0 for r in self.r) else 4

        all_x = ((t for t in itertools.product(range(len(self.f)),repeat=n) if value(t) == 0) for n in range(1,n+1))
        all_b = (itertools.product(range(len(self.f)),repeat=n) for n in range(n-1,-1,-1))
        all_addresses = itertools.chain.from_iterable(itertools.product(b_vals,x_vals) for b_vals, x_vals in zip(all_b,all_x))
        ext = [ct_from_tuple(b)(ct_from_tuple(x).fixed_point()) for b,x in all_addresses]

        return Interval.closed(min(ext),max(ext))


    def extend(self, ctr_f_itbl):
        """Return a generator of extensions of ctr_f by all possible words."""
        return itertools.chain.from_iterable((ctr_f.compose(f) for f in self.f) for ctr_f in ctr_f_itbl)

    def extend_with_prb(self, ctr_f_itbl):
        """Return a generator of extensions of ctr_f by all possible words along with the probabilities associated to those words."""
        return itertools.chain.from_iterable(((p,ctr_f,ctr_f.compose(f)) for p,f in zip(self.p,self.f)) for ctr_f in ctr_f_itbl)

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

class Neighbour(CtrFunc):
    @classmethod
    def from_f(cls,f,iv):
        func = CtrFunc(iv.delta,iv.a).inverse().compose(f)
        return cls(func.r,func.d)

    def to_f(cls, net_iv):
        return net_iv.normalization_func().compose(self)

    @property
    def L(self):
        return self.r

    @property
    def a(self):
        return -self.d

class NeighbourSet(tuple):
    """
    A neigbour set is just a sorted tuple of neighbours (contraction functions).
    """
    def __new__(cls,nb_itbl):
        self = super().__new__(cls,sorted(set(nb_itbl)))
        self.lmax = max(abs(nb.L) for nb in self)
        return self

    def __str__(self):
        return ", ".join(f"({nb.d},{nb.L})" for nb in self)

    @classmethod
    def base(cls):
        return cls((Neighbour(d=0,r=1),))

    def maximal_nbs(self):
        """Compute the neigbours in self.nb_set of maximal size"""
        return (nb for nb in self if abs(nb.L) == self.lmax)

    def nonmaximal_nbs(self):
        """Compute the neigbours in self.nb_set not of maximal size"""
        return (nb for nb in self if abs(nb.L) != self.lmax)

class TransitionMatrix(SymbolicMatrix):
    """TransitionMatrix class to represent the transition matrix associated to an edge in the transition graph.
    """
    def pos_row(self):
        return all(any(x!=0 for x in row) for row in self.matrix)

    def spectral_radius(self):
        return np.abs(np.linalg.eigvals(np.array(self))).max()
