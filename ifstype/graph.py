import graph_tool as gt
import numpy as np
import typing
from functools import reduce
import itertools
from operator import mul

from .exact import Interval, Constants as C, Fraction
from .ifs import AffineFunc, NetInterval, TransitionMatrix


class LocalDim:
    """Class to represent local dimensions, for easy printing / visualization of what the values are."""
    def __init__(self, spr, ln):
        self.spr = spr
        self.length = ln

    def __float__(self):
        return float(np.log(float(self.spr))/np.log(float(self.length)))

    def __str__(self):
        return f"log({self.spr})/log({self.length})"

class EdgeInfo(typing.NamedTuple):
    """Named tuple to represent the edge information intrinsic to an edge (other than the source and target)"""
    t_index: Fraction
    length: Fraction
    transition: TransitionMatrix

class TransitionGraph:
    """
    Guarantees: vertex 0 is the root.
    """
    def __init__(self,root,ifs):
        # public attributes
        self.is_wft = False
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
        L = reduce(mul, (self._edge_info[e].length for e in loop))

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
    # functions to compute fixed net intervals, from symbolic representation (a sequence of edges)
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

    def _extend_iv_by_edge(self,iv,edge):
        """Given an interval iv, compute the resulting interval after stepping along the edge"""
        return Interval(
                iv.a+self._edge_info[edge].t_index*iv.delta,
                iv.a+(self._edge_info[edge].t_index+self._edge_info[edge].length)*iv.delta)


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
        alpha = self._nb_set_from_vtx[edge_seq[-1].source()].lmax*reduce(mul,(self._edge_info[e].length for e in edge_seq[:-1]), self.root.delta)
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
        - vprop[v] == 0 if v is in the attractor (essential class if weak finite type)
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

    def _rebuild_lookups(self):
        self._edge_lookup = {self.g.edge_index[e]:e for e in self.g.edges()}
        self._vtx_lookup = {self._nb_set_from_vtx[v]:v for v in self.g.vertices()}
                
    def remove_terminal_vertices(self):
        """Remove vertices with out degree 0, repair the _vtx_lookup, and return the corresponding list of neighbour sets removed.
        Warning: this changes edge and vertex numbering."""
        # TODO: use graph views / filters instead of mutating the graph
        # use gt.GraphView(self.g,vfilt=zero_degs), should be a get_vertex option for out_degree 0
        # then update zero_degs each time? does this change the GraphView?
        # then can prune directly from the GraphView on the last run?
        zero_degs = gt.util.find_vertex(self.g,"out",0)
        while(len(zero_degs) > 0):
            self.g.remove_vertex(zero_degs)
            zero_degs = gt.util.find_vertex(self.g,"out",0)

        # rebuild vertex and edge lookup
        self._rebuild_lookups()

        # new code:

        #  temp_g = self.g
        #  while True:
            #  zero_deg_filter = (temp_g.get_out_degrees(temp_g.get_vertices()) == 0)
            #  temp_g = GraphView(temp_g,vfilt=zero_deg_filter)

        # just need to prune graph into GraphView

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

