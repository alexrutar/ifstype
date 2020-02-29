""":mod:`ifstype.ifs`
=====================

This module contains core classes and methods for use in the simulation
of any iterated function system.

Public module attributes:

* :class:`AffineFunc`
* :class:`IFS`
* :func:`ifs_family`
* :class:`Neighbour`
* :class:`NeighbourSet`
* :class:`NetInterval`
* :class:`TransitionMatrix`

"""
from numbers import Real
import itertools
import numpy as np
import inspect
import functools
from typing import Union, Tuple, Iterable, Sequence, AbstractSet, Callable
import attr

from .exact import (
    Constants as C, Interval, Fraction,
    SymbolicRing, SymbolicMatrix, SymbolicElement
)

@attr.s(frozen=True, slots=True)
class AffineFunc:
    """An AffineFunc is a real-valued affine function ``f(x)=r*x+d``.

    This class is an immutable storage class.
    When called without arguments, it defaults to the identity function.

    Core attributes:

    * :attr:`AffineFunc.r`
    * :attr:`AffineFunc.d`

    Methods for function properties:

    * :meth:`AffineFunc.fixed_point`
    * :meth:`AffineFunc.interval`

    Methods for behaviour as a function:

    * :meth:`AffineFunc.__call__`
    * :meth:`AffineFunc.compose`
    * :meth:`AffineFunc.inverse`

    """
    r: Real = attr.ib(default=C.n_1) #: The affine factor
    d: Real = attr.ib(default=C.n_0) #: The translation factor

    def fixed_point(self) -> Real:
        """Compute the point which the affine function fixes.

        >>> aff = AffineFunc(Fraction(2),Fraction(3))
        >>> aff.fixed_point()
        Fraction(-3, 1)
        >>> aff(aff.fixed_point())
        Fraction(-3, 1)

        The affine factor must have absolute value not equal to 1.

        :raises ValueError: if linear coefficient is equal to 1
        :return: ``p`` such that ``f(p) == p``.
        """
        try:
            return self.d/(C.n_1-self.r)
        except ZeroDivisionError:
            raise ValueError("Linear coefficient must not be equal to 1")

    def interval(self,initial_iv: Interval=Interval()) -> Interval:
        """Compute the set image ``f(initial_iv)={f(x):x in initial_iv}``
        corresponding to the function.

        >>> aff = AffineFunc(-Fraction(1,2),Fraction(3,4))
        >>> aff.interval()
        Interval(a=Fraction(1, 4), b=Fraction(3, 4))

        :param initial_iv: interval to compute image of
        :return: image of the interval
        """
        left = self(initial_iv.a)
        right = self(initial_iv.b)
        if left <= right:
            return Interval(left,right)
        else:
            return Interval(right,left)


    def __call__(self, x: Real) -> Real:
        """Use the class as a python function.

        >>> aff = AffineFunc(Fraction(2),Fraction(3))
        >>> aff(Fraction(1,2))
        Fraction(4, 1)

        :param x: call value
        :return: result of applying function to x
        """
        return self.r*x+self.d

    def compose(self, g: 'AffineFunc') -> 'AffineFunc':
        """Compute the affine function resulting from right composition with
        the affine function g.

        >>> aff = AffineFunc(Fraction(2),Fraction(1,3))
        >>> afg = AffineFunc(Fraction(2),Fraction(3))
        >>> aff.compose(afg)
        AffineFunc(r=Fraction(4, 1), d=Fraction(19, 3))

        :param g: any affine function
        :return: a new affine function instance ``f(g(x))``
        """
        return self.__class__(self.r*g.r, self.d+self.r*g.d)

    def inverse(self) -> 'AffineFunc':
        """Compute the inverse function of the linear map.

        Requires that :attr:`r` is non-zero.

        >>> aff = AffineFunc(Fraction(2),Fraction(1,3))
        AffineFunc(r=Fraction(1, 2), d=Fraction(-1, 6))

        :return: a new affine function instance ``g(x)`` such that
                 ``g(f(x)=f(g(x))=x``

        :raises ZeroDivisionError: if :attr:`r` is zero
        """
        try:
            return self.__class__(1/self.r,-self.d/self.r)
        except ZeroDivisionError:
            raise ZeroDivisionError("Cannot compute inverse of affine "
                                    "function with linear coefficient 0")


class IFS:
    """A class representing an iterated function system such that the convex
    hull of the invariant compact set is [0,1].

    After initialization, the instance has the following attributes:

    * :attr:`funcs`
    * :attr:`probs`

    The following methods are also defined:

    * :meth:`__init__`
    * :meth:`__str__`
    * :meth:`extend`
    * :meth:`invariant_convex_hull`

    """
    def __init__(self, funcs: Sequence[AffineFunc]) -> None:
        """Initialize iterated function system instance.

        The `funcs` argument is a list of AffineFunc instances with linear
        factor having absolute value strictly between 0 and 1.

        .. warning:: The ``funcs`` to the :attr:`IFS.funcs` attribute, since
                     `funcs` is normalized to have invariant convex hull [0,1].

        :raises ValueError: if the linear coefficients r do not have absolute
                            value strictly between 0 and 1

        :param funcs: a sequence of :class:`AffineFunc` instances
        """
        # check params
        if not all(C.n_0<abs(f.r)<C.n_1 for f in funcs):
            raise ValueError("funcs linear coefficients r must have 0<|r|<1")

        # use symbolic probabilities
        self._syr = SymbolicRing((f"p{i+1}" for i in range(len(funcs))))
        self._funcs = tuple(funcs)
        # return list of probabilities in order
        self._probs = self._syr.get_symbols()

        self._normalize()

    @property
    def funcs(self) -> Sequence[AffineFunc]:
        """Tuple of contraction functions associated to the IFS."""
        return self._funcs

    @property
    def probs(self) -> Sequence[SymbolicElement]:
        """Tuple of symbolic probabilities `ifstype.exact.SymbolicElement`, one
        for each associated to each contraction function.

        If probabilities is set with a sequence of real numbers, they must be
        strictly greater than 0 and sum to 1

        :setter: Associate numeric values to the probabilities which sum to 1.
        """
        return self._probs

    @probs.setter
    def probs(self, probs: Sequence[Real]) -> None:
        if all(C.n_0<=p for p in probs) and sum(probs) != C.n_1:
            raise ValueError(
                "IFS probabilities must be non-negative and sum to 1")

        if len(self.funcs) != len(probs):
            raise ValueError("funcs and probabilities must be the same size")

        self._syr.set_eval({f"p{i+1}":p for i,p in enumerate(probs)})

    def __str__(self) -> str:
        """Return user-readable string representation of the iterated function
        system.

        :return: string representation
        """
        return f"IFS(funcs={self.funcs})"

    def _normalize(self) -> None:
        """Normalize the contraction functions of the instance such that the
        invariant convex hull is [0,1]

        :return: None
        """
        def conv_normalize(aff,interval):
            iv_func = AffineFunc(interval.b-interval.a,interval.a)
            return iv_func.inverse().compose(aff).compose(iv_func)

        cvx_hull = self.invariant_convex_hull()
        self._funcs = tuple(conv_normalize(f,cvx_hull) for f in self.funcs)

    def invariant_convex_hull(self) -> Interval:
        """Compute the convex hull of the invariant set K.

        The algorithm used is adapted from from József Vass' paper, Section
        3.2, in https://arxiv.org/abs/1502.03788.

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

        all_x = [
            (t for t
             in itertools.product(range(len(self.funcs)),repeat=n)
             if value(t) == 0)
                 for n in range(1,n+1)]
        all_b = [
            itertools.product(range(len(self.funcs)),repeat=n) for n
            in range(n-1,-1,-1)]

        all_addresses = itertools.chain.from_iterable(
            itertools.product(b_vals,x_vals) for b_vals, x_vals
            in zip(all_b,all_x))

        ext = [
            ct_from_tuple(b)(ct_from_tuple(x).fixed_point()) for b,x
            in all_addresses]

        return Interval(min(ext),max(ext))

    def extend(
            self,
            aff_iterable: Iterable[AffineFunc],
            with_prob: bool=False
            ) -> Union[
                AbstractSet[Tuple[SymbolicElement,AffineFunc,AffineFunc]],
                AbstractSet[AffineFunc]]:
        """Return a generator of extensions of elements of aff_iterable by all
        possible functions in the IFS.

        Chains the generator where for each function g of aff_iterable and f of
        the instance variable f, we compute g.compose(f).

        If with_prob is True, the elements of the iterables are tuples
        ``(p,aff_init,aff_res)`` where ``p`` is the symbolic probability
        associated with the function, ``aff_init`` is the initial affine
        function, and ``aff_res`` is the resulting affine function after
        composing with some element of :attr:`funcs`.

        :param aff_iterable: an iterable of AffineFunc instances
        :param with_prob: boolean to return additional information.
        :return: an iterable of all extensions
        """
        if with_prob:
            return set(itertools.chain.from_iterable(
                (
                    (p,ctr_f,ctr_f.compose(f)) for p,f
                    in zip(self.probs,self.funcs))
                for ctr_f in aff_iterable))
        else:
            return set(itertools.chain.from_iterable(
                (ctr_f.compose(f) for f in self.funcs) for ctr_f
                in aff_iterable))


def ifs_family(
        ifs_func: Callable[..., Sequence[AffineFunc]]
    ) -> Callable[..., IFS]:
    """Convenience decorator to construct families of iterated function systems
    parametrized by some set of values.

    Used to decorate functions of the form::

        def ifs(probs=def_p, a_1=def_1, ..., a_n=def_n):
            return [AffineFunc(...), ..., AffineFunc(...)]

    where ``a_1,...,a_n`` are arbitrary parameter names, arbitrary default
    values ``def_1,...,def_n`` for the parameters, and ``def_p`` default
    argument for probabilities (see :attr:`.IFS.probs`).
    The arguments must all be specified as keyword arguments.
    The decorated function has the same keyword arguments and returns an
    :class:`IFS` instance from the corresponding probabilities and contraction
    functions.

    .. warning:: The ``funcs`` return value of the undecorated function may not
                 be equal to the :attr:`IFS.funcs` attribute, since `funcs` is
                 sorted and normalized to have invariant convex hull [0,1].

    :param ifs_func: the function being decorated
    :return: decorated function

    """
    @functools.wraps(ifs_func)
    def wrapper(**kwargs):
        fn_kwargs = {
            k:v.default for k,v
            in inspect.signature(ifs_func).parameters.items()}
        func_params = {**fn_kwargs,**kwargs}
        try:
            probs = func_params['probs']
        except KeyError:
            raise ValueError(
                f"Function '{ifs_func.__name__}' decorated by 'ifs_family' "
                f"has no default keyword 'probs'.")

        funcs = ifs_func(**func_params)


        if probs is not None:
            if len(funcs) != len(probs):
                raise ValueError(
                    f"List of funcs returned by '{ifs_func.__name__}' "
                    f"decorated by 'ifs_family' has different length than "
                    f"keyword 'probs'.")

            # add probabilities, resorting if necessary
            sorted_f_pairs = sorted(
                zip(funcs,probs),key=lambda x:(x[0].d,x[0].r))
            funcs = [f[0] for f in sorted_f_pairs]
            probs = [f[1] for f in sorted_f_pairs]
            ifs = IFS(funcs)
            ifs.probs = probs
            return ifs

        else:
            funcs = sorted(funcs,key=lambda x:(x.d,x.r))
            return IFS(funcs)

    return wrapper


@attr.s(frozen=True,slots=True)
class Neighbour(AffineFunc):
    """A neighbour is a special type of normalized affine function, used to
    represent a neighbour of a net interval.

    This class is an immutable storage class.

    Methods for creation:

    * :meth:`from_aff`

    Convenience attributes:

    * :attr:`a`
    * :attr:`L`

    """

    @property
    def a(self) -> Real:
        """The a descriptor of the neighbour."""
        return self.d

    @property
    def L(self) -> Real:
        """The L descriptor of the neighbour."""
        return self.r

    @classmethod
    def from_aff(cls,aff: AffineFunc, interval: Interval) -> 'Neighbour':
        """Create the neighbour corresponding to `aff` by normalizing against
        `interval`.

        :param aff: the affine function
        :param interval: the interval
        :return: the corresponding neighbour
        """
        func = AffineFunc(interval.delta,interval.a).inverse().compose(aff)
        return cls(func.r,func.d)

@attr.s(frozen=True, slots=True)
class NeighbourSet:
    """A neigbour set represents a set of unique neighbours.

    This class is an immutable storage class.

    Attributes:

    * :attr:`neighbours`
    * :attr:`lmax`

    Iteration and containment:

    * :meth:`__iter__`
    * :meth:`__contains__`

    Canonical string representation:

    * :meth:`__str__`

    Convenince methods:

    * :meth:`maximal_nbs`
    * :meth:`nonmaximal_nbs`

    """
    neighbours: AbstractSet[Neighbour] = attr.ib(
        converter=frozenset,
        default=(Neighbour(),)) #: the set of neighbours

    @property
    def lmax(self):
        """The maximum of the absolute value of the affine factor among
        neighbours.
        """
        return max(abs(nb.L) for nb in self)

    def __iter__(self) -> Iterable[Neighbour]:
        """An iterable returning the neighbours.
        The neighbours are not iterated in any particular order.

        :return: the neighbour iterable
        """
        return iter(self.neighbours)

    def __contains__(self,nb: Neighbour) -> bool:
        """Check if a neighbour is contained in the neighbour set

        :param nb: the neighbour
        :return: True if and only if the neighbour is in the neighbour set.
        """
        return nb in self.neighbours

    def __str__(self) -> str:
        """Return a canonical string representation of the neighbour set.
        Assuming that the string representations for the coefficients of the
        neighbour functions are unique, the strings for neighbours are the same
        if and only if the neighbour sets are equal.

        :return: the string representation
        """
        return ", ".join(f"({nb.d},{nb.L})" for nb in sorted(self))

    def maximal_nbs(self) -> Iterable[Neighbour]:
        """Compute the neigbours in the neighbour set of maximal size.

        See also :meth:`nonmaximal_nbs`.

        :return: an iterable of maximal neighbours
        """
        return (nb for nb in self if abs(nb.L) == self.lmax)

    def nonmaximal_nbs(self) -> Iterable[Neighbour]:
        """Compute the neigbours in the neighbour set not of maximal size.

        See also :meth:`maximal_nbs`.

        :return: an iterable of neighbours which are not maximal.
        """
        return (nb for nb in self if abs(nb.L) != self.lmax)


@attr.s(frozen=True, slots=True)
class NetInterval(Interval):
    """A special interval type representing a net interval of generation alpha.
    This contains the interval information, as well as the neighbour set.

    In addition to the attributes and methods inherited from the
    :class:`ifstype.exact.Interval` base class, a net interval also has the
    following attributes and methods:

    Attributes:

    * :attr:`alpha`
    * :attr:`nb_set`

    Methods for creation:

    * :meth:`from_funcs`

    Methods for computation with net intervals:

    * :meth:`transition_gen`
    * :meth:`normalization_func`
    * :meth:`containing_funcs`

    """
    alpha: Real = attr.ib(default=C.n_base) #: the generation
    nb_set: NeighbourSet= attr.ib(default=NeighbourSet()) #: the neighbour set

    @classmethod
    def from_funcs(
            cls,
            a: Real,
            b: Real,
            alpha: Real,
            funcs: Iterable[AffineFunc]
        ) -> 'NetInterval':
        """Compute the neighbour set of the net interval given an iterable of
        (not necessarily distinct) affine functions S such that S([0,1])
        contains [`a`,`b`].

        :return: the net interval
        """
        iv = Interval(a,b)
        nb_set = NeighbourSet(Neighbour.from_aff(f,iv) for f in funcs)
        return cls(a,b,alpha,nb_set)

    def transition_gen(self) -> Real:
        """Compute the largest generation less that :attr:`alpha` for which the
        net interval is no longer a valid.

        :return: the transition generation
        """
        return self.nb_set.lmax*self.delta

    def normalization_func(self) -> AffineFunc:
        """Compute the unique affine function S such that S([0,1]) is equal to
        the net interval.

        :return: the affine function
        """
        return AffineFunc(self.delta,self.a)

    def containing_funcs(self) -> Iterable[AffineFunc]:
        """Compute an iterable of functions which correspond to the neighbours
        in the neighbour set of the net interval. Each affine function
        generates a distinct neighbour of the net interval.

        :return: an iterable of affine functions
        """
        return (self.normalization_func().compose(f) for f in self.nb_set)


class TransitionMatrix(SymbolicMatrix):
    """A transition matrix is a special class used to represent the transition
    matrix associated to an edge in the transition graph.

    Methods not inherited from :class:`ifstype.exact.SymbolicMatric`:

    * :meth:`pow_row`
    * :meth:`spectral_radius`
    """
    def pos_row(self) -> bool:
        """Check if the matrix has the positive row property, in other words
        that each row contains a non-zero entry.

        :return: True if and only if it has the positive row property.
        """
        return all(any(x!=0 for x in row) for row in self.matrix)

    def spectral_radius(self) -> Real:
        """Compute the approximate spectral radius of the transition matrix,
        assuming that probabilities have been set in the corresponding IFS.

        See :attr:`.IFS.probs` for details on how to set probability values.

        :return: real-valued spectral radius
        """
        return np.abs(np.linalg.eigvals(np.array(self))).max()

