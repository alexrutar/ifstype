import itertools
from sympy import Matrix
import logging

from .interval import IntervalNet, Interval
from .ifs import Word


class Generations:
    def __init__(self, ifs, nb_set_types={(Interval(0,1),):'0'}):
        self.ifs = ifs
        self._gens = {(1,Interval(0,1)) : IntervalNet(1,[Word.empty()],interval=Interval(0,1))} # dict {alpha: IntervalNet}
        self.nb_set_types=nb_set_types
        self.num_nb_set = 1

    def print_nb_set_types(self):
        for k,v in self.nb_set_types.items():
            print("{} : {}".format(v,k))

    # main construction methods
    def gen(self,alpha, interval=Interval(0,1)):
        """Compute the NetInterval object of generation alpha corresponding to the specified interval, and caches the result in self._gens
        We use as a starting point the set of all words of generation beta where beta >= alpha is minimal, over an interval containing interval
        """
        # compute the starting point, which is the smallest interval of generation beta >
        key = (alpha,interval)
        prev = min((k for k in self._gens.keys() if k[0] >= alpha and interval.subset(k[1])), key=lambda k:(k[0],k[1].delta))
        if prev == key:
            return self._gens[prev]
        elif prev[0] == key[0]:
            # perhaps implement fast restricted view: enough to just drop the elements which are no longer relevant
            self._gens[key] = IntervalNet(alpha, self._gens[prev].words, interval=interval)
            return self._gens[key]
        else:
            start = self._gens[prev].words

            # compute the interval net corresponding to all words extend from "start", with respect to the IFS
            self._gens[key] = IntervalNet(alpha, itertools.chain.from_iterable(self.extend_to_gen(alpha,w) for w in start), interval=interval)

            # add new neighbour types to the governing list
            for nb in self._gens[key].nb_set_types:
                if nb not in self.nb_set_types.keys():
                    self.nb_set_types[nb] = str(self.num_nb_set)
                    self.num_nb_set += 1

            return self._gens[key]

    def gen_from_select(self, alphas, select):
        """
        Compute an interval selection for each value of alpha.
        Given a list of alphas and starting with Interval(0,1), choose a child of the current net interval in the next generation
        alpha, and repeat.

        arguments:
        - alphas: an iterable of Sympy rational values
        - select: any function f:IntervalNet -> NetInterval where NetInterval is in IntervalNet

        returns:
        - a generator which yields pairs (NetInterval, IntervalNet) from generations in alphas
        """
        left = Interval.empty()
        middle = Interval(0,1)
        right = Interval.empty()
        gen = self.gen(1)
        for alpha in alphas:
            # this choice guarantees that we have the left and right intervals in the children
            tmp_gen = self.children(alpha, left|middle|right)
            middle = select(self.children(alpha, middle))
            left,right = tmp_gen.adjacent(middle)
            gen = self.gen(alpha, interval=left|middle|right)
            yield (middle,gen)


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

    def nb_set(self, net_iv):
        "Compute the neighbour set of a net_iv"
        net = self.gen(net_iv.alpha, interval=net_iv) # you only need the restricted view at net_iv
        if net_iv not in net.net:
            raise ValueError("Interval is not a net interval of generation {}".format(net_iv.alpha))
        else:
            return net.nb_set(net_iv)

    def nb_set_type(self, net_iv):
        return self.nb_set_types[self.nb_set(net_iv)]

    # relationships and transitions
    def transition(self,beta,net_iv):
        "Compute the transition $T_{beta to net_iv.alpha}(net_iv)$ matrix for net interval net_iv in generation alpha with respect to parent in generation beta > alpha"
        child = net_iv
        child_nb_set = self.nb_set(child)
        child_net = self.gen(child.alpha)

        parent = self.parent(beta,net_iv)
        parent_nb_set = self.nb_set(beta,parent)
        parent_net = self.gen(beta)

        transition = [[0]*len(parent_nb_set) for _ in child_nb_set]

        track_writes = [[set() for _ in parent_nb_set] for _ in child_nb_set]

        logging.info("-"*40 + "\nTransition: {} -> {} on {}".format(beta,alpha,net_iv))

        for iv,word in parent_net.word_contributions(parent):
            col = parent_nb_set.index(iv)
            logging.info("Start word: {}".format(word))
            for new in self.extend_to_gen(alpha,word): # extend those words to be of generation alpha
                nm = child_net.norm(new.interval(),child)
                if nm in child_nb_set: # check if the (normalized) extension is an element of the neighbour set of the child
                    row = child_nb_set.index(nm) # the index is the index of the child
                    newp = word.diff(new).p
                    transition[row][col] += newp # get the p value of the extension
                    track_writes[row][col].add(newp)
                    logging.info("Extension: {}, {}, {}".format(new, new.r,new.rm))
                    logging.info("Write: ({},{}) = {}".format(row,col,newp))
        return Matrix(transition)

    def children(self, alpha, net_iv):
        "Compute the children (in generation alpha) of an interval as a restricted interval net"
        return self.gen(alpha, interval=net_iv)

    def parent(self, alpha, net_iv):
        "Compute the parent (in generation alpha) of an interval"
        for par in self.gen(alpha).net:
            if net_iv.subset(par):
                return par
        return None

