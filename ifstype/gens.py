import itertools
import typing
import pdb
import ctypes
from sympy import Matrix, Rational
from collections import defaultdict
from sortedcontainers import SortedList,SortedDict,SortedSet
import copy

from .interval import Interval, NetInterval
from .ifs import Word, C
from .neighbour import NeighbourSet, Neighbour, InfiniteNbMgr, FiniteNbMgr

#  GLOBAL = 0

#  def GLOBAL_val():
    #  if GLOBAL != 0:
        #  return ctypes.cast(GLOBAL, ctypes.py_object).value

class Gen:
    """A Gen is essentially a forgetful generation with masking.
    One can determine the complete behaviour of the children when you know what the functions associated to the words are.
    In this sense, we only save the words as reference for the functions.
    If forgetful is set to true, only one word per function is saved.
    """
    def __init__(self,genkey,words,forgetful=False):
        self.alpha = genkey.alpha
        self.interval = genkey.interval

        # warning: open intersection is not sufficient, since contributions to endpoints is also important!
        self._words = set(w for w in words if not (w.interval() & self.interval).is_empty()) # only keep words which intersect with interval

        # creates a list of functions and a dictionary showing which words generate a given function
        if not forgetful:
            self.f_dict = defaultdict(set)
            for w in self._words:
                self.f_dict[w.f].add(w)
        else:
            self.f_dict = {w.f:(w,) for w in self._words}
            self._words = tuple(w[0] for w in w.f.values())

        endpoints = SortedSet(filter(
            lambda ep: ep in self.interval,
            itertools.chain.from_iterable(map(lambda f:f.interval(), self.fs()))))

        self.net = tuple(NetInterval(self.alpha,a,b) for a,b in zip(endpoints,endpoints[1:]))

    def words(self):
        return iter(self._words)
    def __iter__(self):
        return iter(self.net)

    def nb_set(self,net_iv):
        "Return the local neighbour set of a net_iv"
        if net_iv not in self.net:
            raise ValueError(f"{net_iv} is not a net interval of generation {self.alpha} with view {self.interval}")
        return NeighbourSet(Neighbour.from_f(f,net_iv) for f in self.fs() if f.interval().open_intersect(net_iv))

    def intervals(self):
        return SortedSet(f.interval() for f in self.fs())

    def fs(self):
        return self.f_dict.keys()

    def adjacent(self, net_iv):
        "Compute the adjacent intervals (left,right) if they exist (returns empty if they don't)"
        idx = self.net.index(net_iv)
        try:
            left = self.net[idx-1]
        except IndexError:
            left = Interval.empty()
        try:
            right = self.net[idx+1]
        except IndexError:
            right = Interval.empty()
        return (left,right)

    def random_net_iv(self,n=1):
        "Get a random net interval"
        return random.sample(self.net, n)

    def f_words(self,f):
        "Get the the words which generate a given function"
        try:
            return self.f_dict[f]
        except KeyError:
            raise ValueError(f"{f} is not a function of generation {self.alpha} with view {self.interval}")


    def argmin_delta(self):
        "Returns the net interval with minimum size"
        return min(self.net, key=lambda x:x.delta())

    def word_contributions(self,net_f):
        "Returns an iterable of (function,word) pairs which contribute to net_iv (for recursive calculation purposes)"
        return ((nb.to_f(net_iv),word) for nb in self.nb_set(net_iv) for word in self.f_words(nb.to_f(net_iv)))


class GenKey(typing.NamedTuple):
    """A key used in Generations._gen to compare previous computations, equipped with a comparison relation (not a total ordering!).
    """
    alpha: Rational
    interval: Interval

    def __str__(self):
        return str((self.alpha,self.interval))

    def __hash__(self):
        return hash((self.alpha,self.interval))

    def __ge__(self,other):
        return self.alpha >= other.alpha and self.interval.supset(other.interval)
    def __le__(self,other):
        return self.alpha <= other.alpha and self.interval.subset(other.interval)

    def __gt__(self,other):
        return (self.alpha >= other.alpha and self.interval.proper_supset(other.interval)) or \
                (self.alpha > other.alpha and self.interval.supset(other.interval))
    def __lt__(self,other):
        return (self.alpha <= other.alpha and self.interval.proper_subset(other.interval)) or \
                (self.alpha < other.alpha and self.interval.subset(other.interval))
    def __eq__(self,other):
        return self.alpha == other.alpha and self.interval == other.interval
    def __neq__(self,other):
        return self.alpha != other.alpha or self.interval != other.interval


class Generations:
    def __init__(self, ifs, finite_type=False, existing_nb_sets=None):
        self.ifs = ifs
        # set default generation
        gk = GenKey(C.n_base,Interval(C.n_0,C.n_1))
        self._gens = {gk : Gen(gk, [Word.empty()])}

        # either finite type or infinite type
        if finite_type:
            self.nb_mgr = FiniteNbMgr(existing_nb_sets)
            # initialize the neighbour set
            to_update = [NetInterval(C.n_base,C.n_0,C.n_1)]
            while(len(to_update)>0):
                new = []
                for net_iv in to_update:
                    new_nb = self.nb_set(net_iv)
                    if new_nb not in self.nb_mgr:
                        self.nb_mgr.add(new_nb)
                        ch = self.im_children(net_iv)
                        new.extend(ch)
                to_update = new

        else:
            self.nb_mgr = InfiniteNbMgr(existing_nb_sets)

    def im_children(self, net_iv):
        "Compute the immediate children of net_iv"
        gamma = max(abs(nb.L) for nb in self.nb_set(net_iv))*net_iv.delta()
        return self.children(gamma, net_iv)

    def ttype(self, net_iv):
        ch_gen = self.im_children(net_iv)
        return tuple(((ch.a - net_iv.a)/net_iv.delta(), ch.delta()/net_iv.delta(), self.nb_set(ch)) for ch in ch_gen)

    def children(self, alpha, net_iv):
        "Compute the children (in generation alpha) of an interval as a restricted interval net"
        out = self.gen(alpha, interval=Interval(net_iv.a,net_iv.b))
        return out

    def nb_set(self, net_iv):
        "Compute the neighbour set of a net_iv, wrapper for Gen.nb_set"
        return self.gen(net_iv.alpha, interval=Interval(net_iv.a,net_iv.b)).nb_set(net_iv) # only need the restricted view at net_iv

    # neighbour set functions
    def all_nb_sets(self):
        return str(self.nb_mgr)

    def nb_set_type(self, net_iv):
        "Look up the associated string with the nb_set_type"
        return self.nb_mgr.nb_set_type(self.nb_set(net_iv))

    # main construction methods
    def gen(self, alpha,interval=None):
        """Compute the NetInterval object of generation alpha corresponding to the specified interval, and caches the result in self._gens
        We use as a starting point the set of all words of generation beta where beta >= alpha is minimal, over an interval containing interval
        """
        # compute the starting point, which is the smallest interval of generation beta > alpha
        if interval is None:
            interval = Interval(0,1)

        genkey = GenKey(alpha,interval)
        prev = min((k for k in self._gens.keys() if k >= genkey),
                key=lambda k:(k.alpha,k.interval.delta()))

        if prev == genkey:
            # key already exists
            return self._gens[prev]
        elif prev.alpha == genkey.alpha:
            # alpha already exists: just need to restrict existing gen
            # TODO: implement this faster as an internal restriction method, i.e. create as function from previous gen
            # use self.__class__() to initialize new instance
            new_gen = Gen(genkey, self._gens[prev].words())
        else:
            # alpha doesn't exist, extend from words of next best (first by closest alpha, then by smallest size interval containing it)
            # TODO: after faster restriction, restrict containing interval first, then extend
            new_gen = Gen(genkey, itertools.chain.from_iterable(self.extend_to_gen(genkey.alpha,w) for w in self._gens[prev].words()))

        # add new neighbour types to the governing list (in finite case, perhaps update does nothing)
        self.nb_mgr.update(new_gen.nb_set(net_iv) for net_iv in new_gen)

        self._gens[genkey] = new_gen
        return new_gen

    def gen_from_select(self, alphas, select):
        """
        Compute an interval selection for each value of alpha.
        Given a list of alphas and starting with Interval(0,1), choose a child of the current net interval in the next generation
        alpha, and repeat.

        arguments:
        - alphas: an iterable of Sympy rational values
        - select: any function f:Gen -> NetInterval where NetInterval is in Gen

        returns:
        - a generator which yields pairs (NetInterval, IntervalNet) from generations in alphas
        """
        left = Interval.empty()
        middle = Interval(0,1)
        right = Interval.empty()
        gen = self.gen(1)
        for alpha in alphas:
            # this choice guarantees that we have the left and right intervals in the children
            tmp_gen = self.children(alpha, Interval(left.a,right.b))
            middle = select(self.children(alpha, middle))
            left,right = tmp_gen.adjacent(middle)
            gen = self.gen(alpha, interval=left|middle|right)
            yield (middle,gen)


    # TODO: implement version with pruning for restricted views
    def extend_to_gen(self,alpha,word):
        """Returns an iterable of all words of generation alpha having word as a prefix.

        arguments:
        - alpha: a Sympy rational value to extend word to
        - word: a word of generation before alpha

        returns:
        - a generator which yields words of generation alpha
        """
        if word.is_gen(alpha):
            yield word
        else:
            to_extend = [word]
            while(len(to_extend) > 0):
                itbl = itertools.product(self.ifs.idx,to_extend)
                to_extend = []
                for i,word in itbl:
                    new = self.ifs.extend(word,i)
                    if new.is_gen(alpha):
                        yield new
                    else:
                        to_extend.append(new)


    # relationships and transitions
    def transition(self,beta,net_iv):
        "Compute the transition $T_{beta to net_iv.alpha}(net_iv)$ matrix for net interval net_iv in generation alpha with respect to parent in generation beta > alpha"
        child = net_iv
        child_nb_set = self.nb_set(child)
        child_net = self.gen(child.alpha)

        parent = self.parent(beta,net_iv)
        parent_nb_set = self.nb_set(beta,parent)
        parent_net = self.gen(beta)

        transition = [[C.n_0]*len(parent_nb_set) for _ in child_nb_set]

        for f, word in parent_net.word_contributions(parent):
            col = parent_nb_set.index(iv)
            for new in self.extend_to_gen(alpha,word): # extend those words to be of generation alpha
                nb = Neighbour.from_f(new.f,child)
                if nb in child_nb_set: # check if the (normalized) extension is an element of the neighbour set of the child
                    row = child_nb_set.index(nm)
                    newp = word.diff(new).p
                    transition[row][col] += newp

        return Matrix(transition)


    def parent(self, alpha, net_iv):
        "Compute the parent (in generation alpha) of an interval"
        for par in self.gen(alpha).net:
            if net_iv.subset(par):
                return par
        return None

