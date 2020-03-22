
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
        print(f"Vertex set {self._nb_set_from_vtx[vtx]} -> {new_nb_set}")
        self._nb_set_from_vtx[vtx] = new_nb_set

        # remove the row from the outgoing transition matrices
        for e in vtx.out_edges():
            new_mat = self.edge_info(e).transition.remove_row(i)
            self._edge_info[e] = self.edge_info(e).new_matrix(new_mat)

        # remove the column from the incoming transition matrices
        # if the matrix becomes empty, delete the edge
        for e in vtx.in_edges():
            new_mat = self.edge_info(e).transition.remove_column(i)
            if new_mat.is_empty():
                self.g.remove_edge(e)
            else:
                self._edge_info[e] = self.edge_info(e).new_matrix(new_mat)


    def _merge_vertices(self,v1,v2):
        """Merge vertices in vtx_lst by combining all in-edges and all out-edges,
        and then removing all trivial loops and duplicated edges (edges with same edge_info)"""
        for e in v1.in_edges():
            # create new edge with same source and new traget
            new_e = self.g.add_edge(e.source(),v2)
            # copy the edge property as well
            self._edge_info[new_e] = self._edge_info[e]

        self.g.remove_vertex(v1)

    def _merge_vertices(self, vtx_list, nb_set):
        # pretty complicated algorithm since we need to remove degree 1 vts
        # first (recursively), then merge the independent vertices
        new_v = self.g.new_vertex()
        self._nb_set_from_vtx[new_v] = nb_set

        # create new vertex with all edges
        for vtx in vtx_list:
            for e in vtx.in_edges():
                new_e = self.g.add_edge(e.source(), new_v)
                self._edge_info[new_e] = self._edge_info[e]
            for e in vtx.out_edges():
                new_e = self.g.add_edge(new_v,e.target())
                self._edge_info[new_e] = self._edge_info[e]

        # remove existing vertices
        self.g.remove_vertex(vtx_list)
        unique_edges = {EdgeInfo}

        # remove doubled edges / identity self-loops
            

    def _rebuild_lookups(self):
        """After all the edges have been pruned, identify vertices which may now have the same (reduced) vertex set.
        Then re-index the edges so that they are labelled contiguously, fix the vertex and edge lookups so that they are in order
        """

        # merge vertices which have the same neighbour set
        reverse = defaultdict(list)
        for v in self.g.vertices():
            reverse[self._nb_set_from_vtx[v]].append(v)
        for k,v in reverse.items():
            print(f"{k} : {[int(l) for l in v]}")

        for nb_set, vtx_lst in reverse.items():
            if len(vtx_lst) >= 2:
                self._merge_vertices(vtx_lst, nb_set)

        # fix reverse lookup dictionaries to reflect changes
        self._edge_lookup = {self.g.edge_index[e]:e for e in self.g.edges()}
        self._vtx_lookup = {self._nb_set_from_vtx[v]:v for v in self.g.vertices()}
        #  print(self._vtx_lookup)

    def reduce_transition_graph(self):
        # reduce the graph
        to_reduce = self._filter_pos_row()
        while(len(to_reduce) > 0):
            for red in to_reduce:
                self._prune_vtx(red)
            to_reduce = self._filter_pos_row()

        zero_degs = gt.util.find_vertex(self.g,"out",0)
        self.g.remove_vertex(zero_degs)


        # identify vertices which may now have the same neighbour set
        # rebuild vertex and edge lookup
        self._rebuild_lookups()

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
        self._edge_lookup = {self.g.edge_index[e]:e for e in self.g.edges()}
        self._vtx_lookup = {self._nb_set_from_vtx[v]:v for v in self.g.vertices()}
        # TODO: write new code
        # new code:

        #  temp_g = self.g
        #  while True:
            #  zero_deg_filter = (temp_g.get_out_degrees(temp_g.get_vertices()) == 0)
            #  temp_g = GraphView(temp_g,vfilt=zero_deg_filter)

        # just need to prune graph into GraphView
