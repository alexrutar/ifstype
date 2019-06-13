from sympy import Rational
import operator
import functools
import itertools
import typing
from sortedcontainers import SortedList
from functools import reduce

from .interval import Interval

def app(word,target):
    return functools.reduce(operator.mul, [target[i] for i in word],1)

class Word:
    def __init__(self, ab, rm, p, f):
        self.rm = rm
        self.ab = tuple(ab)
        self.p = p
        self.len = len(ab)
        self.f = f
        self.r = f.r

    def interval(self,*args,**kwargs):
        "Create the interval associated with word"
        return self.f.interval(*args,**kwargs)

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
        return cls("",None,Rational(1),CtrFunc.id())

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
        return Interval(self.r,self.rm)

    def is_gen(self,alpha):
        "True if word is of generation alpha"
        return self.r < alpha and self.rm >= alpha
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
        return cls(Rational(1),Rational(0))

    def compose(self, ct_f):
        return CtrFunc(self.r*ct_f.r, self.a+self.r*ct_f.a)

    def interval(self,iv=Interval(0,1)):
        left = self(iv.a)
        right = self(iv.b)
        if left <= right:
            return Interval(left,right)
        else:
            return Interval(right,left)

    def __call__(self, x):
        return self.r*x+self.a
    def __repr__(self):
        return "f:x*{}+{}".format(self.r,self.a)
    def __str__(self):
        return "f:x*{}+{}".format(self.r,self.a)

# define a class containing constants
def constant(f):
    def fset(self, value):
        raise AttributeError("Cannot change constant values")
    def fget(self):
        return f()
    return property(fget, fset)

class _Const:
    @constant
    def n_base():
        return Rational(2)
    @constant
    def n_0():
        return Rational(0)
    @constant
    def n_1():
        return Rational(1)
    @constant
    def iv_0_1():
        return Interval(Rational(0),Rational(1))

C = _Const()

class IFS:
    # an IFS is essentially a factory for words and ctr funcs
    def __init__(self, *funcs):
        "funcs is a pair (ct_f,pr)"
        funcs = sorted(funcs,key=lambda x:x[0].a)

        self.f = [f[0] for f in funcs]
        self.p = [f[1] for f in funcs]

        self.r = [f.r for f in self.f]
        self.a = [f.a for f in self.f]

        self.idx = tuple(range(len(funcs)))
        self.rmin = min(abs(r) for r in self.r)
        self.rmax = max(abs(r) for r in self.r)

    def transition_gens(self,start=1):
        abs_r = sorted(set(abs(r)) for r in self.r)
        jumps = map(lambda n:self.rmax**n,itertools.count(0))
        sorted_r=SortedSet([])
        for n,j in enumerate(jumps):
            sorted_r.update(reduce(operator.mul,itertools.product(abs_r,repeat=n),1))
            cur = sorted_r.pop()


    def __str__(self):
        return str((self.f,self.p))

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
