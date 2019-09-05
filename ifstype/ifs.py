"""
.. module:: ifs
   :synopsis: Core functionality for IFS simulation.

.. moduleauthor:: Alex Rutar <github.com/alexrutar>
"""

import itertools
import numpy as np
import typing
import inspect
import functools
from quicktions import Fraction

from .exact import Constants as C, Interval
from .exact.symbolic import SymbolicRing, SymbolicMatrix

class AffineFunc(typing.NamedTuple):
    """An AffineFunc is a real-valued affine function ``f(x)=r*x+d``.

    An AffineFunc is immutable.

    Core attributes:

    * :attr:`AffineFunc.r`
    * :attr:`AffineFunc.d`

    Methods for creating the function:

    * :meth:`AffineFunc.id`

    Methods to query attributes:

    * :meth:`AffineFunc.endpoints`
    * :meth:`AffineFunc.fixed_point`
    * :meth:`AffineFunc.interval`
    * :meth:`AffineFunc.sign`

    Methods for behaviour as a function:

    * :meth:`AffineFunc.__call__`
    * :meth:`AffineFunc.compose`
    * :meth:`AffineFunc.inverse`

    Methods for miscellany:

    * :meth:`AffineFunc.__repr__`
    * :meth:`AffineFunc.__str__`

    """
    r: Fraction
    d: Fraction

    @classmethod
    def id(cls):
        """Create the affine function ``f`` such that ``f(x)=x``, the identity affine function.

        :return: identity affine function
        """
        return cls(C.n_1,C.n_0)

    def sign(self):
        """Determine the sign of the affine function.

        :return: True if :attr:`r` is strictly positive, else False
        """
        return self.r>0

    def fixed_point(self):
        """Compute the point which the affine function fixes.

        >>> ctf = AffineFunc(Fraction(2),Fraction(3))
        >>> ctf.fixed_point()
        Fraction(-3, 1)
        >>> ctf(ctf.fixed_point())
        Fraction(-3, 1)

        The affine factor must have absolute value not equal to 1.

        :return: ``p`` such that ``f(p) == p``.
        """
        return self.d/(C.n_1-self.r)

    def endpoints(self):
        """Compute the endpoints of the interval ``f([0,1])`` corresponding to the affine function.

        :return: ``(f(0),f(1))``
        """
        return (self.d,self.d+self.r)

    def interval(self,initial_iv=Interval.closed(C.n_0,C.n_1)):
        """Compute the set image ``f(initial_iv)={f(x):x in initial_iv}`` corresponding to the function.

        :param initial_iv: interval to compute image of
        :return: image of the interval
        """
        left = self(initial_iv.a)
        right = self(initial_iv.b)
        if left <= right:
            return Interval(a=left,b=right,has_left=initial_iv.has_left,has_right=initial_iv.has_right)
        else:
            return Interval(a=right,b=left,has_left=initial_iv.has_right,has_right=initial_iv.has_left)


    def __call__(self, x):
        """Use the class as a python function.

        :param x: call value
        :return: ``f(x)``
        """
        return self.r*x+self.d

    def inverse(self):
        """Compute the inverse function of the linear map.

        Requires that :attr:`r` is non-zero.

        :return: a new affine function instance ``g(x)`` such that ``g(f(x)=f(g(x))=x``
        """
        return AffineFunc(1/self.r,-self.d/self.r)

    def compose(self, g):
        """Compute the affine function resulting from right composition with the affine function g.

        :param g: any affine function
        :return: a new affine function instance ``f(g(x))``
        """
        return AffineFunc(self.r*g.r, self.d+self.r*g.d)


    def __repr__(self):
        """Return string representation of affine function.


        :return: string representation
        """
        return f"AffineFunc(r={self.r},d={self.d})"

    def __str__(self):
        """Return user-readable string representation of affine function.

        :return: string representation.
        """
        return f"f:x*({self.r}) + ({self.d})"

# add proper docstrings for sphinx-autodoc
AffineFunc.r.__doc__ = """\
The affine factor.
"""
AffineFunc.d.__doc__ = """\
The translation factor.
"""
AffineFunc.__new__.__doc__ = """\
Creates a new affine function instance representing the function f(x) = :attr:`r` x + :attr:`d`.

Note that the parameters :attr:`r` and :attr:`d` must be hashable numeric types.

:param r: the linear term
:param d: the constant term
:return: affine function instance
"""


def ifs_family(ifs_func):
    """Convenience decorator to construct families of iterated function systems parametrized by some set of values.

    Used to decorate functions of the form

    .. code-block::

       def ifs(probs=def_p, a_1=def_1, ..., a_n=def_n):
           return [AffineFunc(...), ..., AffineFunc(...)

    where ``a_1,...,a_n`` are arbitrary parameter names, arbitrary default values ``def_1,...,def_n`` for the parameters, and ``def_p`` default argument for probabilities (see :meth:`IFS.set_probs`).
    The arguments must all be specified as keyword arguments.
    The decorated function can be called with the same keyword arguments and returns an :class:`IFS` instance from the corresponding probabilities and contraction functions.

    :param ifs_func: the function being decorated
    :return: decorated function
    """
    @functools.wraps(ifs_func)
    def wrapper(**kwargs):
        fn_kwargs = {k:v.default for k,v in inspect.signature(ifs_func).parameters.items()}
        func_params = {**fn_kwargs,**kwargs}
        try:
            probs = func_params['probs']
        except KeyError:
            raise KeyError("Function decorated by 'ifs_family' has no default keyword 'probs'.")
        if probs is None:
            return IFS.uniform_p(*ifs_func(**func_params))
        else:
            return IFS(ifs_func(**func_params),probs)
    return wrapper


class IFS:
    """A class representing an iterated function system such that the convex hull of the invariant compact set is [0,1].

    After initialization, the instance is guaranteed to have instance variables

    * :var f: normalized contraction functions sorted first by shift, then by contraction factor

    """
    def __init__(self, funcs, probabilities=None):
        """Initialize iterated function system instance.

        The `funcs` argument is a finite iterable of AffineFunc instances with linear factor having absolute value strictly between 0 and 1.

        If `probabilities` is specified, it must be a list of real values strictly between 0 and 1 with sum 1, with the same length as `funcs`.

        :param funcs: an iterable of :class:`AffineFunc` instances
        :param probabilities: an optional list of probabilities
        """
        # check params
        funcs = list(funcs)
        assert all(C.n_0<abs(f.r)<C.n_1 for f in funcs), "IFS affine factors must have 0<|r|<1"

        self.syr = SymbolicRing((f"p{i+1}" for i in range(len(funcs))))
        self.p = self.syr.get_symbols() # use symbolic probabilities

        if probabilities is not None:
            probabilities = list(probabilities)
            assert all(C.n_0<=p for p in probabilities) and sum(probabilities) == C.n_1, "IFS probabilities must be non-negative and sum to 1"
            sorted_f_pairs = sorted(zip(funcs,probabilities),key=lambda x:(x[0].d,x[0].r))
            self.f = tuple(f[0] for f in sorted_f_pairs)
            probs = tuple(f[1] for f in sorted_f_pairs)
            self.set_probs(probs)

        else:
            self.f = tuple(sorted(funcs,key=lambda x:(x[0].d,x[0].r)))

        self.normalize()

    def set_probs(self,probs):
        assert len(self.f) == len(probs), "funcs and probabilities must be the same size"
        self.syr.set_eval({f"p{i+1}":p for i,p in enumerate(probs)})

    def __str__(self):
        return f"IFS(funcs={self.f})"

    @classmethod
    def uniform_p(cls, *funcs):
        return cls(funcs,[Fraction(1,len(funcs)) for _ in funcs])

    def normalize(self):
        def conv_normalize(aff, interval):
            cur_iv = aff.interval(initial_iv=interval)
            new_iv = (cur_iv - interval.a)/interval.delta
            if aff.r > C.n_0:
                return AffineFunc(new_iv.delta, new_iv.a)
            else:
                return AffineFunc(-new_iv.delta, new_iv.b)

        cvx_hull = self.invariant_convex_hull()
        self.f = [conv_normalize(f,cvx_hull) for f in self.f]
        self.d = [f.d for f in self.f]

    def invariant_convex_hull(self):
        """Convex Hull computation, as adapted from József Vass' paper, which can be found at https://arxiv.org/abs/1502.03788 Section 3.2."""
        def value(tup):
            return sum(1 for a in tup if self.f[a].r < 0) % 2

        def ct_from_tuple(tup):
            "Create the affine function associated to the tuple argument"
            f = AffineFunc.id()
            for letter in tup:
                f = f.compose(self.f[letter])
            return f

        n=2 if all(f.r>0 for f in self.f) else 4

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
        return AffineFunc(self.delta,self.a)

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

class Neighbour(AffineFunc):
    @classmethod
    def from_f(cls,f,iv):
        func = AffineFunc(iv.delta,iv.a).inverse().compose(f)
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
    A neigbour set is just a sorted tuple of neighbours (affine functions).
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
