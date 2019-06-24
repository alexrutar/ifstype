import itertools
import typing
import operator
from sympy import Matrix
from collections import defaultdict
from sortedcontainers import SortedList,SortedDict,SortedSet

from .interval import Interval, NetInterval, View, IntervalSet
from .ifs import Word
from .neighbour import NeighbourSet, Neighbour, InfiniteNbMgr, FiniteNbMgr
from .numeric import Constants as C, Rational


# TODO: Write a class GeometricGen, which is essentially the same as Gen, but it doesn't save the words at all
# subsequent generations can be computed directly based on the IFS dynamics (if it is finite type)
# use this for FiniteType computations, where you can generate the transition maps using the standard object then go to the faster version after computing everything in advance.
class Gen:
    """A Gen is essentially a forgetful generation with masking.
    One can determine the complete behaviour of the children when you know what the functions associated to the words are.
    In this sense, we only save the words as reference for the functions.
    If forgetful is set to true, only one word per function is saved.
    TODO: allow passing forgetful at generations level.
    """
    def __init__(self,genkey,words,forgetful=False,full_K=False):
        self.alpha = genkey.alpha
        # computing the iteration is pretty heavy, there's probably a faster way to do it
        # maybe maintain a representative set? especially in the finite type case, if loading from saved transition maps

        # warning: open intersection is not sufficient, since contributions to endpoints is also important! TODO why?
        self._words = set(w for w in words if not (View(w.interval()) & genkey.view).is_empty) # only keep words which intersect with the view

        if full_K:
            self.view = genkey.view
            self.iteration=Interval.closed(0,1)
        else:
            self.iteration = IntervalSet(w.interval() for w in self._words) # the generation n net
            self.view = genkey.view & self.iteration

        # creates a list of functions and a dictionary showing which words generate a given function
        if not forgetful:
            self.f_dict = defaultdict(set)
            for w in self._words:
                self.f_dict[w.f].add(w)
        else:
            self.f_dict = {w.f:(w,) for w in self._words}
            self._words = tuple(w[0] for w in w.f.values())

        endpoints = SortedSet(ep for ep in itertools.chain.from_iterable((f(C.n_0),f(C.n_1)) for f in self.fs()) if ep in self.view)

        self.net = tuple(NetInterval(a,b,self.alpha) for a,b in zip(endpoints,endpoints[1:]) if self.iteration.contains_interval(Interval.open(a,b)))

    def __str__(self):
        return ", ".join(str(n) for n in self.net)
    def words(self):
        return iter(self._words)

    def __iter__(self):
        return iter(self.net)

    def nb_set(self,net_iv):
        "Return the local neighbour set of a net_iv"
        if net_iv not in self.net:
            raise ValueError(f"{net_iv} is not a net interval of generation {self.alpha} with view {self.view}")
        return NeighbourSet(Neighbour.from_f(f,net_iv) for f in self.fs() if not (f.interval().interior() & net_iv).is_empty)

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
            raise ValueError(f"{f} is not a function of generation {self.alpha} with view {self.view}")


    def argmin_delta(self):
        "Returns the net interval with minimum size"
        return min(self.net, key=operator.attrgetter('delta'))

    def word_contributions(self,net_f):
        "Returns an iterable of (function,word) pairs which contribute to net_iv (for recursive calculation purposes)"
        return ((nb.to_f(net_iv),word) for nb in self.nb_set(net_iv) for word in self.f_words(nb.to_f(net_iv)))


class GenKey(typing.NamedTuple):
    """A key used in Generations._gen to compare previous computations, equipped with a comparison relation (not a total ordering!).
    """
    alpha: Rational
    view: View

    def __str__(self):
        return str((self.alpha,self.view))

    def __hash__(self):
        return hash((self.alpha,self.view))

    def __ge__(self,other):
        return self.alpha >= other.alpha and self.view.supset(other.view)
    def __le__(self,other):
        return self.alpha <= other.alpha and self.view.subset(other.view)

    def __gt__(self,other):
        return (self.alpha >= other.alpha and self.view.proper_supset(other.view)) or \
                (self.alpha > other.alpha and self.view.supset(other.view))
    def __lt__(self,other):
        return (self.alpha <= other.alpha and self.view.proper_subset(other.view)) or \
                (self.alpha < other.alpha and self.view.subset(other.view))
    def __eq__(self,other):
        return self.alpha == other.alpha and self.view == other.view
    def __neq__(self,other):
        return self.alpha != other.alpha or self.view != other.view


class BaseGenerations:
    """
    The keyword arguments are used to provide known facts about the IFS, which can simplify or speed up the computations.
    - finite_type : whether or not the IFS is `finite type`
    - existing_nb_sets : a pre-computed list of neighbour sets, to be used with finite_type=True, to avoid recomputing the neighbour sets
    - full_K : the invariant compact set is the total interval [0,1]
    """
    def __init__(self, ifs, nb_mgr=None):
        self.ifs = ifs
        gk = GenKey(C.n_base,View(Interval(C.n_0,C.n_1)))
        self._gens = {gk : Gen(gk, [Word.empty()])}
        if nb_mgr is None:
            raise ValueError("Must provide a nb_mgr!")
        else:
            self.nb_mgr = nb_mgr

    def im_alpha(self, net_iv):
        """Compute a generation alpha containing the immediate children of `net_iv`.

        :param net_iv: A valid NetInterval of any generation.
        :return: A Gen object.

        """
        return max(abs(nb.L) for nb in self.nb_set(net_iv))*net_iv.delta

    def im_children(self, net_iv):
        """Compute the immediate children of `net_iv`.

        :param net_iv: A valid NetInterval of any generation.
        :return: A Gen object.

        """
        return self.gen(self.im_alpha(net_iv), view=View(net_iv))

    def parent(self, alpha, net_iv):
        "Compute the parent (in generation alpha) of an interval"
        for par in self.gen(alpha).net:
            if net_iv.subset(par):
                return par
        return None

    def ttype(self, net_iv):
        """Compute the transition type of `net_iv`.

        :param net_iv: A valid NetInterval of any generation.
        :return: A transition type tuple.

        """
        ch_gen = self.im_children(net_iv)
        return tuple(((ch.a - net_iv.a)/net_iv.delta, ch.delta/net_iv.delta, self.nb_set(ch)) for ch in ch_gen)

    # neighbour set functions
    def nb_set(self, net_iv):
        """Compute the neighbour set of `net_iv`.
        
        :param net_iv: A valid NetInterval of any generation.
        :return: A NeighbourSet object.

        """
        return self.gen(net_iv.alpha, view=View(net_iv)).nb_set(net_iv) # only need the restricted view at net_iv

    def nb_set_type(self, net_iv):
        "Return a unique integer (indexed from 0) corresponding to the neighbour set (for human readibility)."
        return self.nb_mgr.nb_set_type(self.nb_set(net_iv))

    # main construction methods
    def gen(self, alpha,view=None):
        """Compute the Gen object of generation alpha corresponding to the specified interval, and caches the result in self._gens
        We use as a starting point the set of all words of generation beta where beta >= alpha is minimal, over an interval containing interval
        """
        if view is None:
            view = View(Interval(0,1))

        genkey = GenKey(alpha,view)
        # compute the starting point, which is the smallest interval of generation beta > alpha
        prev = min((k for k in self._gens.keys() if k >= genkey),
                key=lambda k:(k.alpha,k.view.delta))

        if prev == genkey:
            # key already exists
            return self._gens[prev]
        elif prev.alpha == genkey.alpha:
            # alpha already exists: just need to restrict existing gen
            # TODO: implement this faster as an internal restriction method, i.e. create as function from previous gen
            new_gen = Gen(genkey, self._gens[prev].words())
        else:
            # alpha doesn't exist, extend from words of next best (first by closest alpha, then by smallest size interval containing it)
            # TODO: after faster restriction, restrict containing interval first, then extend
            new_gen = Gen(genkey, itertools.chain.from_iterable(self.extend_to_gen(genkey.alpha,w) for w in self._gens[prev].words()))

        self._gens[genkey] = new_gen
        return new_gen

    def gen_from_select(self, select, stop=None):
        """
        Yield a sequence of pairs (NetIv,Gen) with NetIv in Gen.
        If possible, Gen also contains the neighbouring net intervals
        Compute an interval selection for each value of alpha.
        Given a list of alphas and starting with Interval(0,1), choose a child of the current net interval in the next generation
        alpha, and repeat.

        :param select: any function f:Gen -> NetInterval where NetInterval is in Gen

        returns:
        - a generator which yields pairs (NetInterval, Gen) from generations in alphas
        """
        left = Interval.empty()
        middle = NetInterval(C.n_0,C.n_1,C.n_base)
        right = Interval.empty()
        gen = self.gen(1)
        if stop is None:
            itbl = iter(int, 1) # infinite iterable
        else:
            itbl = range(stop)
        for _ in itbl:
            alpha = self.im_alpha(middle)
            # this choice guarantees that we have the left and right intervals in the children
            tmp_gen = self.gen(alpha,view=View(left,middle,right))

            middle = select(self.gen(alpha,view=View(middle)))
            left,right = tmp_gen.adjacent(middle)
            gen = self.gen(alpha, view=View(left,middle,right))
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

