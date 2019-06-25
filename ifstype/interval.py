import typing
import operator
import itertools
from sympy import oo # TODO: write custom infinity class
from bisect import bisect

from .rational import Rational, Constants as C

# methods which mutate a list in place
# union
def _iv_union(lst,iv):
    """Compute the union of lst with iv.
    Mutates lst in place."""
    if iv.is_empty:
        return
    idx_l = bisect(lst,iv.cmp_left)
    left_iv = lst[idx_l-1] if idx_l >= 1 else None
    idx_r = bisect(lst,iv.cmp_right)
    outside_right_iv = lst[idx_r] if idx_r < len(lst) else None
    inside_right_iv = lst[idx_r-1] if 1 <= idx_r else None

    del_start = idx_l
    del_num = idx_r-idx_l
    left = iv.left
    right = iv.right
    
    # check if need to include left IV
    if left_iv is not None and left_iv.cmp_right >= iv.cmp_left:
        del_start -= 1
        left = left_iv.left
        del_num += 1

    # check if need to include right IV
    if outside_right_iv is not None and outside_right_iv.cmp_left <= iv.cmp_right:
        right = outside_right_iv.right
        del_num += 1
    elif inside_right_iv is not None and inside_right_iv.cmp_right >= iv.cmp_right:
        right = inside_right_iv.right


    for _ in range(del_num):
        del lst[del_start]

    lst.insert(del_start,Interval.from_lr(left,right))

# intersection
def _iv_intersect(tup,iv):
    """Compute the intersection of tup with iv.
    Outputs a list! (note-different than _iv_union for implementation reasons - we usually don't repeatedly call intersect, unlike union, so it is better to only make one copy.)"""
    if iv.is_empty:
        return []
    idx_l = bisect(tup,iv.cmp_left)
    idx_r = bisect(tup,iv.cmp_right)
    # intersection has size 0
    if idx_r <= 0:
        return []
    # intersection has size one
    elif idx_l == idx_r:
        return [tup[idx_r-1] & iv]
    else:
        # may need to look at one extra interval to the left, if it exists
        ivl = Interval.empty()
        if idx_l >= 1:
            ivl = tup[idx_l-1] & iv

        new = list(tup[idx_l:idx_r-1])
        new.append(tup[idx_r-1] & iv)

        if not ivl.is_empty:
            new.insert(0,ivl)
        return new

def _iv_subset(iv,tup):
    idx = bisect(tup,iv)
    return (idx >= 1 and iv.subset(tup[idx-1])) or (idx < len(tup) and iv.subset(tup[idx]))


class IntervalSet(tuple):
    """An interval class when you also need support for arbitrary unions and complements.
    iv_gen is an interval generator to intialize, where IntervalSet is computed as the union of the intervals in the generator.
    """
    def __new__(cls, iv_gen=None):
        base = []
        if iv_gen is not None:
            for iv in iv_gen:
                _iv_union(base,iv)
        return super().__new__(cls,base)

    # -----------------------------------------------
    # standard interval properties
    # -----------------------------------------------
    @property
    def is_empty(self):
        return len(self) == 0

    @property
    def is_point(self):
        return len(self) == 1 and self[0].is_point

    @property
    def delta(self):
        return sum(iv.delta for iv in self)

    # -----------------------------------------------
    # inclusion methods
    # -----------------------------------------------
    def __contains__(self, item):
        idx = bisect(self, (item,False))
        return (idx >= 1 and item in self[idx-1]) or (idx < len(self) and item in self[idx])
    def subset(self, other):
        return all(_iv_subset(iv, other) for iv in self)
    def contains_interval(self, iv):
        return _iv_subset(iv, self)
    def supset(self, other):
        return other.subset(self)
    def proper_subset(self,other):
        return self.subset(other) and self != other
    def proper_supset(self,other):
        return other.subset(self) and self != other

    # -----------------------------------------------
    # relative construction methods
    # -----------------------------------------------
    def __and__(self,other):
        out = []
        for iv in other:
            out.extend(_iv_intersect(self,iv))
        return super().__new__(self.__class__,out)

    def __or__(self,other):
        lst = list(self)
        for iv in other:
            _iv_union(lst, iv)
        return super().__new__(self.__class__,lst)
    def interior(self):
        return super().__new__(self.__class,(iv.interior() for iv in self))
    def closure(self):
        return super().__new__(self.__class,(iv.closure() for iv in self))


    # -----------------------------------------------
    # misc methods
    # -----------------------------------------------
    def __str__(self):
        return r"{"+" ∪ ".join(str(iv) for iv in self)+r"}"

    def __repr__(self):
        return f"IntervalSet" + super().__repr__()



def check_empty(f):
    def wrapper(self,cnst):
        if self.is_empty:
            return Interval()
        else:
            return f(self,cnst)
    return wrapper

class Interval(tuple):
    "A fast interval class"
    a = property(operator.itemgetter(0))
    b = property(operator.itemgetter(2))
    cmp_left = property(operator.itemgetter(0,1))
    cmp_right = property(operator.itemgetter(2,3))
    def __new__(cls,a=None,b=None,has_left=True,has_right=True,**kwargs):
        # as arguments, has_left=False means open, has_right=True=False means open
        # internally, to maintain ordering, False on the left represents closed, and False on the right represents open

        if a is None or b is None or a>b or (a==b and (has_left,has_right) != (True,True)):
            return super().__new__(cls,(None,True,None,False,*kwargs.values()))
        else:
            return super().__new__(cls,(a,not has_left, b, has_right,*kwargs.values()))
    # -----------------------------------------------
    # construction methods
    # -----------------------------------------------
    @classmethod
    def from_lr(cls,left=(None,False),right=(None,False)):
        return cls(a=left[0],b=right[0],has_left=left[1],has_right=right[1])

    @classmethod
    def open_infty(cls,a=None):
        return cls(a=a,has_left=False,b=oo,has_right=False)

    @classmethod
    def closed_infty(cls,a=None):
        return cls(a=a,has_left=True,b=oo,has_right=False)

    @classmethod
    def infty_open(cls,b=None):
        return cls(a=-oo,has_left=False,b=b,has_right=False)

    @classmethod
    def Linfty(cls,b=None,has_right=True):
        return cls(a=-oo,has_left=False,b=b,has_right=True)

    @classmethod
    def empty(cls):
        "Can also call with Interval()"
        return cls(None,None)

    @classmethod
    def open_closed(cls,a=None,b=None):
        return cls(a,b,has_left=False,has_right=True)

    @classmethod
    def closed_open(cls,a=None,b=None):
        return cls(a,b,has_left=True,has_right=False)

    @classmethod
    def open(cls,a=None,b=None):
        return cls(a,b,has_left=False,has_right=False)

    @classmethod
    def closed(cls,a=None,b=None):
        return cls(a,b,has_left=True,has_right=True)

    @classmethod
    def point(cls,a=None):
        return cls(a,a,has_left=True,has_right=True)

    def _same(self,a,b):
        """Create a new instance of the same interval type but with new endpoints"""
        _,has_left = self.left
        _,has_right = self.right
        return Interval(a=a,b=b,has_left=has_left,has_right=has_right)
    
    # -----------------------------------------------
    # access methods
    # -----------------------------------------------
    @property
    def has_left(self):
        _,has_left = self.cmp_left
        return not has_left

    @property
    def has_right(self):
        _,has_right = self.cmp_right
        return has_right

    @property
    def left(self):
        a,b = self.cmp_left
        return (a,not b)

    @property
    def right(self):
        return self.cmp_right

    # -----------------------------------------------
    # standard interval properties
    # -----------------------------------------------
    @property
    def is_empty(self):
        return self.a is None and self.b is None

    @property
    def is_point(self):
        return not self.is_empty and self.a == self.b

    @property
    def delta(self):
        if self.is_empty:
            return 0
        else:
            return self.b-self.a

    # -----------------------------------------------
    # relative construction methods
    # -----------------------------------------------
    def interior(self):
        return Interval.open(self.a,self.b)
    def closure(self):
        return Interval.closed(self.a,self.b)

    def __and__(self,other: 'Interval'):
        "intersection"
        if self.is_empty or other.is_empty:
            return Interval()
        else:
            a,l = max(self.cmp_left,other.cmp_left)
            b,r = min(self.cmp_right,other.cmp_right)
            return Interval(a=a,b=b,has_left=not l,has_right = r)

    # -----------------------------------------------
    # inclusion methods
    # -----------------------------------------------
    def __contains__(self, item: Rational):
        return (not self.is_empty) and self.cmp_left <= (item,False) and self.cmp_right >= (item,True)
    def subset(self, other):
        return self.is_empty or ((not other.is_empty) and self.cmp_left >= other.cmp_left and other.cmp_right >= self.cmp_right)
    def supset(self, other):
        return other.subset(self)
    def proper_subset(self,other):
        return self.subset(other) and self != other
    def proper_supset(self,other):
        return other.subset(self) and self != other


    # -----------------------------------------------
    # misc methods
    # -----------------------------------------------
    def __str__(self):
        if self.is_empty:
            return "IvEmpty"
        elif self.is_point:
            return f"IvPoint{{{self.a}}}"

        a,has_left = self.left
        b,has_right = self.right
        lchr = '[' if has_left else '('
        rchr = ']' if has_right else ')'
        return f"Iv{lchr}{a},{b}{rchr}"

    def __repr__(self):
        return f"Interval(left={self.left},right={self.right})"

    def __iter__(self):
        return iter((self.cmp_left,self.cmp_right))

    # -----------------------------------------------
    # math operations
    # -----------------------------------------------
    @check_empty
    def __mul__(self, cnst: Rational):
        if cnst < 0:
            return self._same(self.b * cnst, self.a * cnst)
        else:
            return self._same(self.a * cnst, self.b * cnst)

    @check_empty
    def __truediv__(self,cnst: Rational):
        if cnst < C.n_0:
            return self._same(self.b / cnst, self.a / cnst)
        else:
            return self._same(self.a / cnst, self.b / cnst)

    @check_empty
    def __add__(self,cnst: Rational):
        return self._same(self.a + cnst, self.b + cnst)

    @check_empty
    def __sub__(self,cnst: Rational):
        return self._same(self.a - cnst, self.b - cnst)



class NetInterval(Interval):
    """A special interval type representing a net interval of generation alpha.
    This interval is always closed, and if you perform any operations on it, you just get a plain interval back."""
    alpha = property(operator.itemgetter(4))
    def __new__(cls,a,b,alpha):
        self = super().__new__(cls,a=a,b=b,alpha=alpha)
        return self

    # representation
    def __str__(self):
        return f"NetIv({self.alpha})[{self.a},{self.b}]"
    def __repr__(self):
        return f"NetInterval(left={self.left},right={self.right},alpha={self.alpha})"

class View(IntervalSet):
    def __new__(cls, *ivs):
        return super().__new__(cls,iv_gen=ivs)

