import itertools
import typing

from .numerics.rational import Rational, Constants as C
from .numerics.interval import Interval


class CtrFunc(typing.NamedTuple):
    r: Rational
    d: Rational

    @classmethod
    def id(cls):
        return cls(C.n_1,C.n_0)

    def inverse(self):
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


class IFS:
    # an IFS is essentially a factory for words and ctr funcs
    def __init__(self, funcs, probabilities):
        """funcs is a list of functions"""
        # check params
        assert len(funcs) == len(probabilities), "funcs and probabilities must be the same size"
        assert all(C.n_0<abs(f.r)<C.n_1 for f in funcs), "IFS contraction factors must have 0<|r|<1"
        assert all(C.n_0<=p for p in probabilities) and sum(probabilities) == C.n_1, "IFS probabilities must be non-negative and sum to 1"

        sorted_f_pairs = sorted(zip(funcs,probabilities),key=lambda x:(x[0].d,x[0].r))

        self.idx = tuple(range(len(funcs)))

        self.p = [f[1] for f in sorted_f_pairs]

        self.f = [f[0] for f in sorted_f_pairs]
        self.r = [f.r for f in self.f]
        self.d = [f.d for f in self.f]
        self.normalize()

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

    def __str__(self):
        return "IFS with normalized contraction functions\n- " + "\n- ".join(f"{f}" for f in self.f) + "\nand probabilities\n- " + "\n- ".join(f"{p}" for p in self.p)

    def extend(self, ctr_f_itbl):
        """Return a generator of extensions of ctr_f by all possible words."""
        return itertools.chain.from_iterable((ctr_f.compose(f) for f in self.f) for ctr_f in ctr_f_itbl)

    def extend_with_prb(self, ctr_f_itbl):
        """Return a generator of extensions of ctr_f by all possible words along with the probabilities associated to those words."""
        return itertools.chain.from_iterable(((p,ctr_f,ctr_f.compose(f)) for p,f in zip(self.p,self.f)) for ctr_f in ctr_f_itbl)
