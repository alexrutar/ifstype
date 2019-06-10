from sympy import *
from interval import Interval, NetIntervals, Word
import operator
import functools
import itertools


def app(word,target):
    return functools.reduce(operator.mul, [target[i] for i in word],1)

class IFS:
    """Past IFS: makes generational assumptions, normalizes wrt generation, etc"""
    # funcs is a [(contraction,shift,probability), ...]
    # assumes support is [0,1]
    # a word is a (tup,r) pair
    def __init__(self, funcs, gen_rule):
        self.gen_rule = gen_rule(self)
        self.gens = {}
        funcs.sort(key=lambda x:x[1])

        self.r = [a[0] for a in funcs]
        self.a = [a[1] for a in funcs]
        self.p = [a[2] for a in funcs]

        self.idx = tuple(range(len(funcs)))
        self.ct = len(funcs)
        self.rmin = min(self.r)

    def new_gen(self,n):
        if n not in self.gens.keys():
            self.gens[n] = list(self.gen_rule.new(n))
        return self.gens[n]

    def get_gen(self,w):
        return self.gen_rule.gen(word)

    # methods on words
    def word_new(self,*ab):
        if len(ab) == 0:
            return Word.empty()
        else:
            rm = app(ab[:-1], self.r)
            return Word(ab,rm*self.r[ab[-1]],rm)
    def word_extend(self, word, letter):
        return word.extend(letter,self.r[letter])
    def word_compose(self,word):
        "Create the contraction function associated to word"
        ctr = 1
        shift = 0
        for c in reversed(word.ab):
            ctr *= self.r[c]
            shift = shift*self.r[c] + self.a[c]
        return lambda x: ctr*x+shift
    def word_interval(self,word,interval=Interval(0,1)):
        "Compute the interval associated to a word"
        f = self.word_compose(word)
        left = f(interval.a)
        right = f(interval.b)
        if interval.a < interval.b:
            return Interval(left,right)
        else:
            return Interval(right,left)
    def net_interval(self,n):
        self.net = NetIntervals(self.new_gen(n),self)
        return self.net
