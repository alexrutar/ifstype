"""
.. module:: ifs
   :synopsis: Core functionality for IFS simulation.

.. moduleauthor:: Alex Rutar <github.com/alexrutar>
"""
import numbers
import itertools
import numpy as np
import inspect
import functools
from quicktions import Fraction
import typing
import attr

from .exact import Constants as C, Interval
from .exact.symbolic import SymbolicRing, SymbolicMatrix, SymbolicElement

@attr.s(frozen=True, slots=True)
class AffineFunc:
    """An AffineFunc is a real-valued affine function ``f(x)=r*x+d``.

    An AffineFunc is an immutable storage class.
    When called without arguments, it defaults to the identity function.

    Core attributes:

    * :attr:`AffineFunc.r`
    * :attr:`AffineFunc.d`

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
    r: numbers.Real = attr.ib(default=C.n_1) #: The affine factor
    d: numbers.Real = attr.ib(default=C.n_0) #: The translation factor

    def sign(self) -> bool:
        """Determine the sign of the affine function.

        :return: True if :attr:`r` is strictly positive, else False
        """
        return self.r>0

    def fixed_point(self) -> numbers.Real:
        """Compute the point which the affine function fixes.

        >>> aff = AffineFunc(Fraction(2),Fraction(3))
        >>> aff.fixed_point()
        Fraction(-3, 1)
        >>> aff(aff.fixed_point())
        Fraction(-3, 1)

        The affine factor must have absolute value not equal to 1.

        :return: ``p`` such that ``f(p) == p``.
        """
        return self.d/(C.n_1-self.r)

    def endpoints(self) -> typing.Tuple[numbers.Real,numbers.Real]:
        """Compute the endpoints of the interval ``f([0,1])`` corresponding to the affine function.

        :return: ``(f(0),f(1))``
        """
        return (self.d,self.d+self.r)

    def interval(self,initial_iv=Interval()) -> Interval:
        """Compute the set image ``f(initial_iv)={f(x):x in initial_iv}`` corresponding to the function.

        :param initial_iv: interval to compute image of
        :return: image of the interval
        """
        left = self(initial_iv.a)
        right = self(initial_iv.b)
        if left <= right:
            return Interval(left,right)
        else:
            return Interval(right,left)


    def __call__(self, x:numbers.Real) -> numbers.Real:
        """Use the class as a python function.

        :param x: call value
        :return: result of applying function to x
        """
        return self.r*x+self.d

    def inverse(self) -> 'AffineFunc':
        """Compute the inverse function of the linear map.

        Requires that :attr:`r` is non-zero.

        :return: a new affine function instance ``g(x)`` such that ``g(f(x)=f(g(x))=x``
        """
        return self.__class__(1/self.r,-self.d/self.r)

    def compose(self, g:'AffineFunc') -> 'AffineFunc':
        """Compute the affine function resulting from right composition with the affine function g.

        :param g: any affine function
        :return: a new affine function instance ``f(g(x))``
        """
        return self.__class__(self.r*g.r, self.d+self.r*g.d)

    def __str__(self) -> str:
        """Return user-readable string representation of affine function.

        :return: string representation.
        """
        return f"f:x*({self.r}) + ({self.d})"


class IFS:
    """A class representing an iterated function system such that the convex hull of the invariant compact set is [0,1].

    After initialization, the instance has the following attributes:

    * :attr:`funcs`
    * :attr:`probs`

    The following methods are also defined:

    * :meth:`__init__`
    * :meth:`__str__`
    * :meth:`extend`
    * :meth:`invariant_convex_hull`

    """
    def __init__(self, funcs:typing.Sequence[AffineFunc], probabilities:typing.Optional[typing.Sequence[numbers.Real]]=None) -> None:
        """Initialize iterated function system instance.

        The `funcs` argument is a list of AffineFunc instances with linear factor having absolute value strictly between 0 and 1.
        If `probabilities` is specified, it must be a list of real values strictly between 0 and 1 with sum 1, with the same length as `funcs`.

        .. warning:: The `funcs` argument may not be equal to the :attr:`funcs` attribute, since during initialization `funcs` is sorted and normalized to have invariant convex hull [0,1].

        :param funcs: a sequence of :class:`AffineFunc` instances
        :param probabilities: an optional list of probabilities
        """
        # check params
        assert all(C.n_0<abs(f.r)<C.n_1 for f in funcs), "IFS affine factors must have 0<|r|<1"

        # use symbolic probabilities
        self.syr = SymbolicRing((f"p{i+1}" for i in range(len(funcs))))
        self._probs = self.syr.get_symbols()

        if probabilities is not None:
            # add probabilities, resorting if necessary
            assert all(C.n_0<=p for p in probabilities) and sum(probabilities) == C.n_1, "IFS probabilities must be non-negative and sum to 1"
            sorted_f_pairs = sorted(zip(funcs,probabilities),key=lambda x:(x[0].d,x[0].r))
            self._funcs = tuple(f[0] for f in sorted_f_pairs)
            self.probs = tuple(f[1] for f in sorted_f_pairs)

        else:
            self._funcs = tuple(sorted(funcs,key=lambda x:(x[0].d,x[0].r)))

        self._normalize()

    @classmethod
    def uniform_p(cls, funcs:typing.Sequence[AffineFunc]) -> 'IFS':
        return cls(funcs,[Fraction(1,len(funcs)) for _ in funcs])

    @property
    def funcs(self) -> typing.Sequence[AffineFunc]:
        """Tuple of contraction functions associated to the IFS.
        """
        return self._funcs

    @property
    def probs(self) -> typing.Sequence[SymbolicElement]:
        """Tuple of symbolic probabilities `ifstype.exact.SymbolicElement`, one for each associated to each contraction function.

        :setter: Associate numeric values to the probabilities which sum to 1.
        """
        return self._probs

    @probs.setter
    def probs(self, probs: typing.Sequence[numbers.Real]) -> None:
        assert len(self.funcs) == len(probs), "funcs and probabilities must be the same size"
        self.syr.set_eval({f"p{i+1}":p for i,p in enumerate(probs)})

    def __str__(self) -> str:
        """Return user-readable string representation of the iterated function system.

        :return: string representation
        """
        return f"IFS(funcs={self.funcs})"

    def _normalize(self) -> None:
        """Normalize the contraction functions of the instance such that the invariant convex hull is [0,1]

        :return: None
        """
        def conv_normalize(aff, interval):
            cur_iv = aff.interval(initial_iv=interval)
            new_iv = Interval((cur_iv.a- interval.a)/interval.delta,(cur_iv.b- interval.a)/interval.delta)
            if aff.r > C.n_0:
                return AffineFunc(new_iv.delta, new_iv.a)
            else:
                return AffineFunc(-new_iv.delta, new_iv.b)

        cvx_hull = self.invariant_convex_hull()
        self._funcs = [conv_normalize(f,cvx_hull) for f in self.funcs]

    def invariant_convex_hull(self) -> Interval:
        """Compute the convex hull of the invariant set K.

        The algorithm used is adapted from from József Vass' paper, Section 3.2, in https://arxiv.org/abs/1502.03788.

        :return: closed interval convex hull
        """
        def value(tup):
            return sum(1 for a in tup if self.funcs[a].r < 0) % 2

        def ct_from_tuple(tup):
            "Create the affine function associated to the tuple argument"
            f = AffineFunc()
            for letter in tup:
                f = f.compose(self.funcs[letter])
            return f

        n=2 if all(f.r>0 for f in self.funcs) else 4

        all_x = ((t for t in itertools.product(range(len(self.funcs)),repeat=n) if value(t) == 0) for n in range(1,n+1))
        all_b = (itertools.product(range(len(self.funcs)),repeat=n) for n in range(n-1,-1,-1))
        all_addresses = itertools.chain.from_iterable(itertools.product(b_vals,x_vals) for b_vals, x_vals in zip(all_b,all_x))
        ext = [ct_from_tuple(b)(ct_from_tuple(x).fixed_point()) for b,x in all_addresses]

        return Interval(min(ext),max(ext))

    def extend(self, aff_iterable: typing.Iterable[AffineFunc], with_prob:bool=False) -> typing.Union[typing.Iterable[typing.Tuple[SymbolicElement,AffineFunc,AffineFunc]],typing.Iterable[AffineFunc]]:
        """Return a generator of extensions of elements of aff_iterable by all possible functions in the IFS.

        Chains the generator where for each function g of aff_iterable and f of the instance variable f, we compute g.compose(f).

        If with_prob is True, the elements of the iterables are tuples ``(p,aff_init,aff_res)`` where ``p`` is the symbolic probability associated
        with the function, ``aff_init`` is the initial affine function, and ``aff_res`` is the resulting affine function after composing with some element of :attr:`funcs`

        :param aff_iterable: an iterable of AffineFunc instances
        :param with_prob: boolean to return additional information.
        :return: an iterable of all extensions
        """
        if with_prob:
            return itertools.chain.from_iterable(((p,ctr_f,ctr_f.compose(f)) for p,f in zip(self.probs,self.funcs)) for ctr_f in aff_iterable)
        else:
            return itertools.chain.from_iterable((ctr_f.compose(f) for f in self.funcs) for ctr_f in aff_iterable)

def ifs_family(ifs_func:typing.Callable[..., typing.Sequence[AffineFunc]]) -> typing.Callable[..., IFS]:
    """Convenience decorator to construct families of iterated function systems parametrized by some set of values.

    Used to decorate functions of the form

    .. code-block::

       def ifs(probs=def_p, a_1=def_1, ..., a_n=def_n):
           return [AffineFunc(...), ..., AffineFunc(...)]

    where ``a_1,...,a_n`` are arbitrary parameter names, arbitrary default values ``def_1,...,def_n`` for the parameters, and ``def_p`` default argument for probabilities (see :attr:`IFS.probs`).
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
            return IFS.uniform_p(ifs_func(**func_params))
        else:
            return IFS(ifs_func(**func_params),probs)
    return wrapper



@attr.s(frozen=True,slots=True)
class Neighbour(AffineFunc):
    @classmethod
    def from_aff(cls,f:AffineFunc,iv:Interval) -> 'Neighbour':
        func = AffineFunc(iv.delta,iv.a).inverse().compose(f)
        return cls(func.r,func.d)

    @property
    def L(self) -> numbers.Real:
        return self.r

    @property
    def a(self) -> numbers.Real:
        return -self.d

class NeighbourSet(tuple):
    """
    A neigbour set is just a sorted tuple of neighbours (affine functions).
    """
    def __new__(cls,nb_itbl:typing.Iterable[Neighbour]) -> 'NeighbourSet':
        self = super().__new__(cls,sorted(set(nb_itbl)))
        self.lmax = max(abs(nb.L) for nb in self)
        return self

    def __str__(self) -> str:
        return ", ".join(f"({nb.d},{nb.L})" for nb in self)

    @classmethod
    def base(cls) -> 'NeighbourSet':
        return cls((Neighbour(d=0,r=1),))

    def maximal_nbs(self) -> typing.Iterable[Neighbour]:
        """Compute the neigbours in self.nb_set of maximal size"""
        return (nb for nb in self if abs(nb.L) == self.lmax)

    def nonmaximal_nbs(self) -> typing.Iterable[Neighbour]:
        """Compute the neigbours in self.nb_set not of maximal size"""
        return (nb for nb in self if abs(nb.L) != self.lmax)

@attr.s(frozen=True,slots=True)
class NetInterval(Interval):
    """A special interval type representing a net interval of generation alpha.
    This contains the interval information, as well as the neighbour set.

    In addition to the attributes and methods inherited from the :class:`ifstype.exact.Interval` base class, a net interval also has the following attributes and methods:

    Attributes:

    * :attr:`alpha`
    * :attr:`nb_set`

    Methods for creation:

    * :meth:`from_funcs`

    Methods for computation with net intervals: 

    * :meth:`transition_gen`
    * :meth:`normalization_func`
    * :meth:`containing_funcs`

    Methods for distinctive string representation:

    * :meth:`__repr__`
    * :meth:`__str__`

    """
    alpha: numbers.Real = attr.ib(default=C.n_base) #: the generation
    nb_set: NeighbourSet= attr.ib(default=NeighbourSet.base()) #: the neighbour set

    @classmethod
    def from_funcs(cls,a:numbers.Real,b:numbers.Real,alpha:numbers.Real,funcs:typing.Iterable[AffineFunc]) -> 'NetInterval':
        """Compute the neighbour set of the net interval given an iterable of (not necessarily distinct) affine functions S such that S([0,1]) contains [`a`,`b`].

        :return: the net interval
        """
        iv = Interval(a,b)
        nb_set = NeighbourSet(Neighbour.from_aff(f,iv) for f in funcs)
        return cls(a,b,alpha,nb_set)

    def transition_gen(self) -> numbers.Real:
        """Compute the largest generation less that :attr:`alpha` for which the net interval is no longer a valid.

        :return: the transition generation
        """
        return self.nb_set.lmax*self.delta

    def normalization_func(self) -> AffineFunc:
        """Compute the unique affine function S such that S([0,1]) is equal to the net interval.
        
        :return: the affine function
        """
        return AffineFunc(self.delta,self.a)

    def containing_funcs(self) -> typing.Iterable[AffineFunc]:
        """Compute an iterable of functions which correspond to the neighbours in the neighbour set of the net interval.
        Each affine function generates a distinct neighbour of the net interval.

        :return: an iterable of affine functions
        """
        return (self.normalization_func().compose(f) for f in self.nb_set)

    # representation
    def __str__(self) -> str:
        return f"NetIv({self.alpha})[{self.a},{self.b}]"
    def __repr__(self) -> str:
        return f"NetInterval(left={self.left},right={self.right},alpha={self.alpha},nb_set={self.nb_set})"

class TransitionMatrix(SymbolicMatrix):
    """TransitionMatrix class to represent the transition matrix associated to an edge in the transition graph.
    """
    def pos_row(self) -> bool:
        return all(any(x!=0 for x in row) for row in self.matrix)

    def spectral_radius(self) -> numbers.Real:
        return np.abs(np.linalg.eigvals(np.array(self))).max()

