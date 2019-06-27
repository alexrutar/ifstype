import operator
from functools import reduce
import itertools
import typing
from sortedcontainers import SortedSet,SortedList

from .rational import Rational, Constants as C
from .interval import Interval

def app(word,target):
    return reduce(operator.mul, [target[i] for i in word],1)

class Word:
    def __init__(self, ab, rm, p, f):
        self.rm = rm
        self.ab = tuple(ab)
        self.p = p
        self.len = len(ab)
        self.f = f
        self.r = f.r

    def interval(self):
        "Create the interval associated with word."
        return self.f.interval()

    def copy(self):
        return Word(self.ab, self.rm, self.p, self.f)

    def __hash__(self):
        tup = (self.ab, self.f, self.rm, self.p)
        return hash(tup)

    def __str__(self):
        return str(self.ab)

    def __repr__(self):
        tup = (self.ab, self.f, self.rm, self.p)
        return str(tup)

    @classmethod
    def empty(cls):
        return cls("",None,C.n_1,CtrFunc.id())

    def extend(self, a, ct_f, pnew):
        new_f = self.f.compose(ct_f)
        new_p = self.p * pnew
        new_ab = self.ab + (a,)
        if self.len == 0:
            new_rm = 1
        else:
            new_rm = self.r
        return Word(new_ab, new_rm, new_p, new_f)

    # compute / evaluate generations
    def gen(self):
        return Interval.open_closed(self.r, self.rm)

    def is_gen(self,alpha):
        "True if word is of generation alpha"
        return alpha in self.gen()
    def pre_gen(self,alpha):
        "True if word is not yet of generation alpha"
        return self.r >= alpha
    def post_gen(self,alpha):
        "True if word has proper prefix of generation alpha"
        return self.rm < alpha

    def diff(self, extension):
        "extension has self as a prefix. Get the part of extension not in self"
        if extension.len == self.len:
            return Word.empty()
        else:
            new_rm = extension.rm/self.r
            new_p = extension.p/self.p
            new_ab = extension.ab[len(self.ab):]
            new_f = CtrFunc(extension.r/self.r, (extension.f.a-self.f.a)/self.r)
            return Word(new_ab, new_rm, new_p, new_f)


class CtrFunc(typing.NamedTuple):
    r: Rational
    a: Rational

    @classmethod
    def id(cls):
        return cls(C.n_1,C.n_0)

    def fixed_point(self):
        return self.a/(C.n_1-self.r)

    def compose(self, ct_f):
        return CtrFunc(self.r*ct_f.r, self.a+self.r*ct_f.a)

    def interval(self,iv=Interval.closed(C.n_0,C.n_1)):
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
        return self.r*x+self.a
    def __repr__(self):
        return f"f:{self.r}*x + {self.a}"
    def __str__(self):
        return f"f:{self.r}*x + {self.a}"


class IFS:
    # an IFS is essentially a factory for words and ctr funcs
    def __init__(self, funcs, probabilities):
        """funcs is a (CtrFunc, p) pair where 0<p<1"""
        # check params
        assert all(C.n_0<abs(f.r) and abs(f.r)<C.n_1 for f in funcs), "IFS contraction factors must have 0<|r|<1"
        assert all(C.n_0<=p and p<C.n_1 for p in probabilities) and sum(probabilities) == C.n_1, "IFS probabilities must be non-negative and sum to 1"

        sorted_f_pairs = sorted(zip(funcs,probabilities),key=lambda x:x[0].a)

        if any(f.r<C.n_0 for f in funcs):
            print("Warning: convex hull normalization doesn't work with negative factors! Assuming convex hull is [0,1]")
            self.f = [f[0] for f in sorted_f_pairs]
        else:
            cvx_hull = self.convex_hull(funcs)
            self.f = [f[0].normalize(cvx_hull) for f in sorted_f_pairs]

        self.p = [f[1] for f in sorted_f_pairs]

        self.r = [f.r for f in self.f]
        self.a = [f.a for f in self.f]

        self.idx = tuple(range(len(funcs)))

        ab = [abs(r) for r in self.r]
        self.rmin = min(ab)
        self.rmax = max(ab)

    def __str__(self):
        return "IFS with normalized contraction functions\n- " + "\n- ".join(f"{f}" for f in self.f) + "\nand probabilities\n- " + "\n- ".join(f"{p}" for p in self.p)

    @staticmethod
    def convex_hull(funcs):
        fixed_points = {f.fixed_point() for f in funcs}
        iterates = set(itertools.chain.from_iterable({f(p) for p in fixed_points} for f in funcs))
        return Interval.closed(min(iterates), max(iterates))

    @classmethod
    def uniform_p(cls, *funcs):
        return cls(funcs,[Rational(1,len(funcs)) for _ in funcs])

    def transition_gens(self,stop=0,count=None):
        """
        Compute the numbers which look like a product of |r_i| that are >= stop in increasing order.
        Defaults to an infinite generator.
        TODO: implement start point? (products that are >= start)
        """
        abs_r = set(abs(r) for r in self.r)
        sorted_r=SortedSet([C.n_base])
        itbl = itertools.count(0)
        ct = 0
        for n,j in ((i,self.rmax**i) for i in itbl):
            sorted_r.update(x for x in (reduce(operator.mul,tup,1) for tup in itertools.product(abs_r,repeat=n)) if x>=stop)
            while True:
                try:
                    cur = sorted_r.pop()
                except IndexError:
                    return
                yield cur
                ct += 1
                if count is not None and ct >= count:
                    return
                if cur == j:
                    break

    def ct_from_word(self,*ab):
        "Create the contraction function associated to the tuple argument"
        f = CtrFunc.id()
        for letter in ab:
            f = f.compose(self.f[letter])
        return f

    def new(self,*ab):
        "Create a new word from a tuple"
        if len(ab) == 0:
            return Word.empty()
        else:
            rm = app(ab[:-1], self.r)
            p = app(ab, self.p)
            return Word(ab,rm,p,self.ct_from_word(*ab))

    def extend(self, word, letter):
        "Extend word by letter"
        return word.extend(letter,self.f[letter],self.p[letter])
