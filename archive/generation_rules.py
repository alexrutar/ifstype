# These are generation rules, which are all functions which take as arguments
# guaranteed methods:
#   new: returns an iterable of words of generation n
#   gen: computes the generation of a given word
import itertools
from interval import Word
from sympy import *
class GenRule0:
    def __init__(self, ifs):
        self.ifs = ifs

    def new(self,n,start=None):
        return map(lambda w: self.ifs.word_new(*w),itertools.product(self.ifs.idx,repeat=n))

    def gen(self,word):
        return word.len


class GenRule1:
    def __init__(self, ifs):
        self.ifs = ifs

    def new(self,n,start=None):
        """Compute all words of generation n
        start is an iterable yielding a covering set of words
        This works since, for any infinite word w_inf, there is some index i such that word_inf[:i] is of generation n"""
        if start is None:
            start = GenRule0(self.ifs).new(n-1)
        genl = self.ifs.rmin**n
        start = list(start)
        while(len(start) > 0):
            rep = []
            for i,word in itertools.product(self.ifs.idx,start):
                new = self.ifs.word_extend(word,i)
                if new.rm > genl and new.r <= genl:
                    yield new
                else:
                    rep.append(new)
            start = rep

    def gen(self,word):
        "Compute the generation of a word (possibly None)"
        gen = floor(log(word.r)/log(self.ifs.rmin))
        genp = floor(log(word.rm)/log(self.ifs.rmin))
        if genp < gen:
            return gen
        else:
            return None

class GenRule2:
    def __init__(self, ifs):
        self.ifs = ifs

    def new(self,alpha,start=[Word.empty()]):
        """Compute all words of generation n
        start is an iterable yielding a covering set of words
        This works since, for any infinite word w_inf, there is some index i such that word_inf[:i] is of generation n"""
        start = list(start)
        if 0<alpha and alpha<1:
            while(len(start) > 0):
                rep = []
                for i,word in itertools.product(self.ifs.idx,start):
                    new = self.ifs.word_extend(word,i)
                    if new.rm > alpha and new.r <= alpha:
                        yield new
                    else:
                        rep.append(new)
                start = rep
        else:
            yield Word.empty()

    def gen(self,word):
        "Returns an interval (a,b) such that for any alpha in (a,b), word is of generation alpha"
        return Interval(word.r,word.rm)
