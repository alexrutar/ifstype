""":mod:`ifstype.graph`
=======================

This module implements the :class:`ifstype.graph.TransitionGraph` class which
implements the transition graph of an IFS and related computations which depend
only on the transition graph.

This module also implements some related transition graph classes.

Public module attributes:

* :class:`TransitionGraph`
* :class:`LocalDim`
* :class:`EdgeInfo`
* :class:`SAdjacencyMatrix`

"""

import graph_tool as gt
import numpy as np
from collections import defaultdict
from typing import Tuple, List, NamedTuple
from functools import reduce
import itertools
from operator import mul
from numbers import Real, Complex

from .exact import (
    Interval, Constants as C, Fraction, Exact,
    SymbolicRing, SymbolicMatrix, SymbolicElement
)

from .ifs import AffineFunc, NetInterval, TransitionMatrix


class LocalDim:
    """Represent local dimensions for easy printing and visualization of what
    the values are.

    Special methods:

    * :meth:`__init__`
    * :meth:`__float__`
    * :meth:`__str__`
    """
    def __init__(self, spr:Real, ln:Real) -> None:
        """Initialize the local dimension with two parameters.

        :param spr: the spectral radius of the transition matrix of the path
        :param ln: the length of the path
        """
        self.spr = spr
        self.length = ln

    def __float__(self) -> float:
        """Return the actual float approximation.

        :return: float approximation
        """
        return float(np.log(float(self.spr))/np.log(float(self.length)))

    def __str__(self) -> str:
        """Return a human-readable string representation

        :return: string representation
        """


class EdgeInfo(NamedTuple):
    """Represent the information intrinsic to a fixed edge (other than the
    source and the target.
    """
    t_index: Exact
    measure: Exact
    length: Exact
    transition: TransitionMatrix

    def new_matrix(self,new_tr_matrix):
        return self.__class__(self.t_index, self.measure, self.length, new_tr_matrix)

    def __mul__(self,other):
        return self.__class__(
                self.t_index + self.measure*other.t_index,
                self.measure*other.measure,
                self.length*other.length,
                self.transition*other.transition)


class SAdjacencyMatrix(SymbolicMatrix):
    """Represent s-adjacency matrices, which are symbolic adjacency matrices
    where the entries are sums of elements of the form r^s, where r is fixed
    and s can vary.

    The matrix is stored symbolically, but particular numpy.ndarray instances
    can be generated automatically.

    >>> s_adj = SAdjacencyMatrix([[(1,),(2,3)],[(),(2,)]])
    >>> print(s_adj)
    [[1, (2)^s + (3)^s],
     [0, (2)^s        ]]
    >>> s_adj.set_s_val(0.5)
    >>> s_adj.spectral_radius()
    1.4142135623730951
    >>> s_adj.compute_s_val()
    (0.99993896484375, 1)


    Initialization:

    * :meth:`__init__`

    Methods:

    * :meth:`set_s_val`
    * :meth:`spectral_radius`
    * :meth:`compute_s_val`

    """
    def __init__(self, mat_values:List[List[Tuple[Real,...]]]) -> None:
        """Initialize the s-adjacency matrix. The only parameter
        ``mat_values`` is a double-nested list of tuples. A tuple
        (r1,...,rn) in position (i,j) represents the entry r1^s + ... + rn^s in
        position (i,j) of the s-adjacency matrix. Empty tuples are treated as
        0.

        :param mat_values: the matrix values

        :raises ValueError: if ``mat_values`` is not a square matrix

        """
        n = len(mat_values)
        if any(len(sub_lst) != n for sub_lst in mat_values):
            raise ValueError("s-adjacency matrix must be a square matrix")

        # extract all possible values in all tuples that are not 0,1 and create
        # a symbolic ring on those elements
        self._vals = {e for lst in mat_values for tup in lst for e in tup
                if e != 0 and e != 1}
        self._syr = SymbolicRing(self._symb_lookup(e) for e in self._vals)

        # evaluate the tuples into corresponding SymbolicRing elements
        new_mat_values_gen = (
            tuple(sum(self._term_lookup(e) for e in tup) for tup in lst)
            for lst in mat_values)

        super().__init__(new_mat_values_gen)

    def _symb_lookup(self,e:Real) -> str:
        """Get the string representation corresponding to a given element.

        :param e: an element which is a term in ``mat_values``
        :return: the string representation

        """
        # get the "SymbolicElement" string corresponding to an eval element
        return f"({e})^s"

    def _term_lookup(self,e:Real) -> SymbolicElement:
        """Get the SymbolicElement corresponding to a given element.

        :param e: an element which is a term in ``mat_values``
        :return: the symbolic element

        """
        if e == 0:
            return 0
        elif e == 1:
            return 1
        else:
            return self._syr.term(self._symb_lookup(e))

    def set_s_val(self,s:Real) -> None:
        """Set the current evaluation value of s to be a fixed real number
        strictly greater than 0

        :param s: evaluation strictly greater than 0

        :raises ValueError: if s <= 0

        """
        if not 0 < s:
            raise ValueError("Invalid s-value not in range 0<s")
        dct={self._symb_lookup(e):float(e)**s for e in self._vals}
        self._syr.set_eval(dct)


    def spectral_radius(self) -> Real:
        """Compute the spectral radius at the current fixed s value.

        .. warning:: The s value must be set using :meth:`set_s_val` before
                     calling this function.
        """
        return np.abs(np.linalg.eigvals(np.array(self))).max()


    def compute_s_val(self,tol:Real=10**(-4)) -> Real:
        """Compute a value s (within :param tol:) such that the spectral
        radius of the s-adjacency matrix is 1. If all entries have value
        between 0 and 1, this value is unique. Returns a lower bound and upper
        bound on s.

        :param tol: error tolerance
        :return: pair (lower s, upper s)
        """
        s_below = 0
        s_above = 1

        # binary search
        while s_above - s_below > tol:
            s_mp = (s_above + s_below)/2
            self.set_s_val(s_mp)
            if self.spectral_radius() < 1:
                s_above = s_mp
            else:
                s_below = s_mp

        return (s_below, s_above)


class TransitionGraph:
    """
    Guarantees: vertex 0 is the root.
    """
    def __init__(self,root,ifs):
        # public attributes
        self.is_fnc = False
        self.g = gt.Graph()
        self.root = root # root net interval
        self.ifs = ifs # associated IFS

        # graph properties
        self._edge_info = self.g.new_edge_property("object")
        self._nb_set_from_vtx = self.g.new_vertex_property("object") # correspondence of vertex objects to neighbour sets

        # additional associations
        self._vtx_lookup = {} # correspondence of neighbour set objects to vertices
        self._edge_lookup = {} # correspondence of edge indices to edge descriptors

    def all_nb_sets(self):
        "Iterable for all neighbour sets, which also allows inclusion checking."
        return self._vtx_lookup.keys()

    # -------------------------------------------------
    # graph properties
    # -------------------------------------------------
    def get_identifier(self, nb_set):
        "Get a unique identifier corresponding to the neighbour set."
        return str(self._vtx_lookup[nb_set])

    def get_nb_set(self, vtx):
        "Get the neighbour set corresponding to a given vertex."
        return self._nb_set_from_vtx[vtx]

    def has_pos_row(self):
        "Check if the transition graph has the positive row property; in other words, that ever transition matrix has a positive entry in every row."
        # TODO: this is wrong / not correct: need to check positive rows along all paths
        return all(e_inf.transition.pos_row() for e_inf in self._edge_info)

    def edge_info(self,e,by_label=False):
        """Return the EdgeInfo object associated with the edge."""
        if by_label:
            return self._edge_info[self._edge_lookup[e]]
        else:
            return self._edge_info[e]


    # -------------------------------------------------
    # compute local dimensions
    # -------------------------------------------------
    def _is_valid_path(self, path, by_label=False):
        """Check if a sequences if edges is in fact a valid path."""
        return all(e1.target() == e2.source() for e1,e2 in zip(path,path[1:]))

    def local_dim(self, loop, by_label=True):
        """Compute the local dimension associated with a loop."""
        if by_label:
            loop = [self._edge_lookup[label] for label in loop]
            assert self._is_valid_path(loop), f"{loop} is not a valid path"
            assert loop[-1].target() == loop[0].source(), f"{loop} start vertex and end vertex are distinct"

        mat = reduce(mul,(self._edge_info[e].transition for e in loop))
        L = reduce(mul, (self._edge_info[e].measure for e in loop))

        return LocalDim(mat.spectral_radius(),L)

    def vertex_local_dims(self,vtx,search_depth=1):
        """Returns a sorted list of all possible local dimensions for a cycle (no repeated vertices) containing vtx."""
        cycles = gt.topology.all_paths(self.g,vtx,vtx,edges=True)
        cycle_prods = itertools.product(cycles, repeat=search_depth)
        return sorted(float(self.local_dim(list(itertools.chain.from_iterable(loop_tups)),by_label=False)) for loop_tups in cycle_prods)
    
    def essential_local_dims(self,search_depth=1):
        ldims = [self.vertex_local_dims(vtx,search_depth=search_depth) for vtx in self.essential_class().get_vertices()]
        return Interval.closed(min(l[0] for l in ldims), max(l[-1] for l in ldims))
    

    # -------------------------------------------------
    # compute adjacency matrix and Hausdorff dimensions
    # -------------------------------------------------

    def adjacency_matrix(self):
        """Compute the weighted adjacency matrix with respect to the edge
        length function"""
        # generate the adjacency matrix of the essential class
        ess = self.essential_class()
        vtx_arr = ess.get_vertices()
        idx_lookup = {vtx:i for i,vtx in enumerate(vtx_arr)}
        evals = [[() for _ in vtx_arr] for _ in vtx_arr]
        for e in ess.edges():
            sr = idx_lookup[int(e.source())]
            tg = idx_lookup[int(e.target())]
            evals[sr][tg] = evals[sr][tg]+ (self.edge_info(e).length,)

        return SAdjacencyMatrix(evals)

    def hausdorff_dim(self,tol=10**(-4)):
        return self.adjacency_matrix().compute_s_val()



    # -------------------------------------------------
    # functions to compute fixed net intervals
    #  from symbolic representation (a sequence of edges)
    # -------------------------------------------------
    def net_ivs_below(self,start=None):
        if start is None:
            start = self.root

        explored_vtxs = set()
        out = set()
        to_explore = [start]
        while(len(to_explore) > 0):
            out.update(to_explore)
            all_children = [self.children(net_iv) for net_iv in to_explore if net_iv.nb_set not in explored_vtxs]
            explored_vtxs.update(net_iv.nb_set for net_iv in to_explore)
            to_explore = list(itertools.chain.from_iterable(all_children))
        return out

    def net_ivs_below_depth(self,depth,start=None):
        if start is None:
            start = self.root

        out = set()
        to_explore = [start]
        while(depth>=0):
            depth -= 1
            out.update(to_explore)
            all_children = [self.children(net_iv) for net_iv in to_explore]
            to_explore = list(itertools.chain.from_iterable(all_children))
        return out

    def _extend_iv_by_edge(self,iv,edge):
        """Given an interval iv, compute the resulting interval after stepping along the edge"""
        return Interval(
                iv.a+self._edge_info[edge].t_index*iv.delta,
                iv.a+(self._edge_info[edge].t_index+self._edge_info[edge].measure)*iv.delta)


    def children(self, net_iv):
        """Compute the children of net_iv using the internal graph properties."""
        vtx = self._vtx_lookup[net_iv.nb_set]
        ivs_and_nbs = ((self._extend_iv_by_edge(net_iv,e),self._nb_set_from_vtx[e.target()]) for e in vtx.out_edges())
        alpha = net_iv.nb_set.lmax*net_iv.delta
        return [NetInterval(iv.a,iv.b,alpha,nb) for iv,nb in ivs_and_nbs]
        

    def net_iv_from_edges(self, edge_seq, by_label=True):
        # first compute where the interval ends up
        if len(edge_seq) == 0:
            return self.root

        if by_label:
            edge_seq = [self._edge_lookup[label] for label in edge_seq]
        iv = reduce(self._extend_iv_by_edge, edge_seq, self.root)

        nb_set = self._nb_set_from_vtx[edge_seq[-1].target()]
        alpha = self._nb_set_from_vtx[edge_seq[-1].source()].lmax*reduce(mul,(self._edge_info[e].measure for e in edge_seq[:-1]), self.root.delta)
        return NetInterval(iv.a,iv.b,alpha,nb_set)


    # -------------------------------------------------
    # compute topological aspects of the graph
    # -------------------------------------------------
    def essential_class(self):
        """Returns a graph_view consisting of the essential class of the graph."""
        components, _, att = gt.topology.label_components(self.g,attractors=True) # get the components and the attractor
        attractor = np.array([(att == True)[comp] == 1 for comp in components.a])
        return gt.GraphView(self.g,vfilt=attractor)

    def components(self):
        """Create a vertex property map `vprop` such that for any v,
        - vprop[v] == -1 if v is an isolated component with no self loop
        - vprop[v] == 0 if v is in the attractor (essential class if satisfies finite neighbour condition)
        - vprop[b] == n for 1 <= n <= N where N is the number of non-trivial components of the graph.
        """
        components, hist, att = gt.topology.label_components(self.g,attractors=True) # get the components and the attractor
        
        # boolean array for attractor of the graph
        attractor = np.array([(att == True)[comp] == 1 for comp in components.a])

        # boolean array for isolated vertices of the graph
        isolated= np.array([(hist == 1)[components[vtx]] and not vtx in vtx.out_neighbours() for vtx in self.g.vertices()],dtype=bool)

        # set isolated components to be type -1
        np.place(components.a, isolated, -1)
        
        # set attractor to be type 0
        np.place(components.a, attractor, 0)

        # set all other loop classes to be type 1 through n
        uniques = np.unique(components.a)
        lookup_map = {elem:idx for idx,elem in enumerate(uniques[uniques > 0],1)}
        components.a = np.array([lookup_map.get(x,x) for x in components.a])

        return components


    # -------------------------------------------------
    # construction methods
    # -------------------------------------------------
    def add_nb_set(self, nb_set, edge_info=None, parent=None):
        if nb_set not in self._vtx_lookup.keys():
            v = self.g.add_vertex()
            self._nb_set_from_vtx[v] = nb_set
            self._vtx_lookup[nb_set] = v
        if parent in self.all_nb_sets():
            e = self.g.add_edge(self._vtx_lookup[parent],self._vtx_lookup[nb_set])
            self._edge_lookup[self.g.edge_index[e]] = e # register edge by index
            self._edge_info[e] = edge_info

    def _filter_pos_row(self):
        """Construct a list of (vtx,i,nb) triples where vtx is a vertex and nb
        is a neighbour such that row i of every outgoing transition matrix is a
        row of zeros.
        """
        non_reduced_list = []
        #  dct = defaultdict(list)
        for vtx in self.g.vertices():
            # list of outgoing transition matrices
            tr_mats = [self.edge_info(e).transition.matrix for e in vtx.out_edges()]
            for i,nb in enumerate(self.get_nb_set(vtx).sorted_iter()):
                # check if there are no offspring in any child
                if all(x == 0 for mat in tr_mats for x in mat[i]):
                    non_reduced_list.append((vtx,i,nb))
        
        return non_reduced_list

    def _prune_vtx(self,vtx_trip):
        """Prune the out edges corresponding to vtx_trip, where
        vtx_trip = (vtx,i,nb) is a triple corresponding to the vertex, the
        row in the outgoing transition matrices that must be removed, and the
        neighbour set corresponding to the row.

        This process removes row i from the outgoing transition matrix, removes
        the neighbour from the neighbour set of vtx, and then removes the
        column i from all the incoming edges to the vertex set.
        If that column makes the transition matrix have size 0, the
        corresponding edge is removed

        """
        vtx,i,nb = vtx_trip

        # update the neighbour set
        new_nb_set = self.get_nb_set(vtx).remove_nb(nb)
        self._nb_set_from_vtx[vtx] = new_nb_set

        # remove the row from the outgoing transition matrices
        for e in vtx.out_edges():
            new_mat = self.edge_info(e).transition.remove_row(i)
            self._edge_info[e] = self.edge_info(e).new_matrix(new_mat)

        # remove the column from the incoming transition matrices
        # if the matrix becomes empty, delete the edge
        for e in vtx.in_edges():
            new_mat = self.edge_info(e).transition.remove_column(i)
            #  if new_mat.is_empty():
                #  self.g.remove_edge(e)
            #  else:
            self._edge_info[e] = self.edge_info(e).new_matrix(new_mat)


    def non_red_nbs(self):
        #  to_reduce = self._filter_pos_row()
        #  non_reduced_neighbours = set()
        #  while(len(to_reduce) > 0):
            #  for red in to_reduce:
                #  # add the neighbour to the list of reduced neighbours
                #  non_reduced_neighbours.add(red[2])
                #  self._prune_vtx(red)
            #  for e in self.g.edges():
                #  if self._edge_info[e].transition.is_empty():
                    #  self.g.remove_edge(e)

            #  to_reduce = self._filter_pos_row()

        #  return non_reduced_neighbours.union(self._terminal_nb_sets)
        return set()

    def _rebuild_lookups(self):
        #  reindex the edges in a sane order and rebuild the lookups
        edges = [(e.source(),e.target(),self._edge_info[e]) for e in self.g.edges()]
        self.g.clear_edges()
        self.g.reindex_edges()

        for src,trg,info in edges:
            new_e = self.g.add_edge(src,trg)
            self._edge_info[new_e] = info

        self._edge_lookup = {self.g.edge_index[e]:e for e in self.g.edges()}
        self._vtx_lookup = {self._nb_set_from_vtx[v]:v for v in self.g.vertices()}

    def _collapse_one(self):
        for vtx in self.g.vertices():
            out = list(vtx.out_edges())
            if len(out) == 1:
                out_edge = out[0]
                for in_edge in vtx.in_edges():
                    new_e = self.g.add_edge(in_edge.source(),out_edge.target())
                    self._edge_info[new_e] = self._edge_info[in_edge]*self._edge_info[out_edge]

                self.g.remove_vertex(vtx)
                break

    def collapse(self):
        n_prev = self.g.num_vertices()
        n = 0
        while True:
            self._collapse_one()
            n_new = self.g.num_vertices()
            if n_new == n_prev:
                break
            else:
                n_prev = n_new
        self._rebuild_lookups()

    def remove_terminal_vertices(self):
        """Remove vertices with out degree 0, repair the _vtx_lookup, and return the corresponding list of neighbour sets removed.
        Warning: this changes edge and vertex numbering."""
        # TODO: use graph views / filters instead of mutating the graph
        # use gt.GraphView(self.g,vfilt=zero_degs), should be a get_vertex option for out_degree 0
        # then update zero_degs each time? does this change the GraphView?
        # then can prune directly from the GraphView on the last run?
        self._terminal_nb_sets = set()
        zero_degs = gt.util.find_vertex(self.g,"out",0)
        while(len(zero_degs) > 0):
            for vtx in zero_degs:
                self._terminal_nb_sets.update(self._nb_set_from_vtx[vtx])
            self.g.remove_vertex(zero_degs)
            zero_degs = gt.util.find_vertex(self.g,"out",0)

        self._rebuild_lookups()


    # -------------------------------------------------
    # saving and loading transition graphs
    # -------------------------------------------------
    # TODO: create internal property maps
    # TODO: write methods to re-load the properties that are not internal property maps
    # TODO: write save methods to filename, and load methods from filename, creating the additional properties upon reloading
    # learn how to implement a python pickle protocol? should be able to include the graph in this too ...
    def save(self, filename):
        # internalize graph properties
        self.g.graph_properties["edge_info"] = self._edge_info
        self.g.graph_properties["nb_sets"] = self._nb_set_from_vtx

