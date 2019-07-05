# -*- coding: utf-8 -*-

# ****************************************************************************
#       Copyright (C) 2019 Wenjie Fang <fwjmath@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.graphs.graph import Graph

class BijectiveMap:
    r"""
    Combinatorial maps for bijective use
    
    This module implements combinatorial maps via canonical labeling with
    rotation systems, and also the basic operations needed to implement
    algorithmically map manipulation for bijections. Plotting is totally
    delegated to other part of sagemath, which may not be adequate to represent
    all the complexity of maps (for instance double edges and loops).
    """

    def __init__(self):
        r"""
        Constructor, returns a map with only the root and two points
        """
        
        # Available corner index
        self._accucnt = 2
        # Number of edges
        self._size = 1
        # Involution indicating edges
        # If the edge is horizontal, the corners linked to an edge are those
        # below on left and above on right.
        self._edges = {0: 1, 1: 0}
        # Permutation indicating corners around each vertex in counter-clockwise
        self._vertices = {0: 0, 1: 1}
        # Inverse permutation for vertices
        self._vertices_inv = {0: 0, 1: 1}
        # Permutation for corners around each face in counter-clockwise
        # !!! Lazy update !!!
        self._faces = {0: 1, 1: 0}
        # Flag indicating if the map is canonical
        # Definition of canonical: 
        # -- all corner labels are distinct and between 0 and size*2 - 1
        # -- self._accucnt == self.size * 2
        # -- faces data is updated
        self._canonical = True
        # Flag indicating if the map is connected
        self._connected = True
        return

    def __check_corner(self, c):
        r"""
        Internal function. Checks if the given corner is valid
        """
        if c not in self._edges:
            raise ValueError("Given corner not in the map")
        return
        
    def __update_vertices_inv(self):
        r"""
        Internal function, updates the inverse permutation for vertices
        """
        self._vertices_inv = {}
        for x in self._vertices:
            self._vertices_inv[self._vertices[x]] = x
        return

    def debug_print(self):
        r"""
        Debug function, used to get all the information
        """
        print(self._edges)
        print(self._vertices)
        return
    
    @staticmethod
    def __from_raw_direct(accucnt, size, edges, vertices):
        r"""
        Internal function, takes on data without checking
        """
        m = BijectiveMap()
        m._accucnt = accucnt
        m._size = size
        m._edges = edges
        m._vertices = vertices
        m.__update_vertices_inv()
        m._faces = {}
        m._canonical = False
        m._connected = False
        return m
    
    @staticmethod
    def from_data(edges, vertices, root):
        r"""
        Returns a combinatorial map with the given data. Checks for sanity.
        
        INPUT:
    
        -- ``edges``: describes how corners are connected by edges
        
        -- ``vertices``: describes how corners surround vertices
        
        -- ``root``: indicates the root corner
        """
        corners = {}
        counter = 0
        # Check that edges is an involution
        for key in edges:
            if key != edges[edges[key]]:
                raise ValueError("The data of edges is not an involution.")
            corners[key] = counter
            counter += 1
        
        # Check that vertices has keys and images all in corners
        for key in vertices:
            if (key not in corners) or (vertices[key] not in corners):
                raise ValueError("Invalid values in the data of vertices")
        
        # Give the root corner a label 0
        if root not in corners:
            raise ValueError("Root not presented in the data")
        olabel = corners[root]
        victim = None
        for key in corners:
            if corners[key] == olabel:
                victim = key
                break
        corners[victim] = olabel
        corners[root] = 0
        
        # Relabel the data and canonicalize
        newedges = {}
        newvertices = {}
        for key in edges:
            newedges[corners[key]] = corners[edges[key]]
        for key in vertices:
            newvertices[corners[key]] = corners[vertices[key]]
        m = BijectiveMap.__from_raw_direct(counter, counter // 2, 
                                                newedges, newvertices)
        m.__canonicalize()
        return m
    
    def __perform_operation(self):
        r"""
        Internal function. Operations on maps must call it to indicate pending
        changes.
        """
        self._canonical = False
    
    def is_canonical(self):
        r"""
        Returns whether the map is with canonicalized labels
        
        OUTPUT:
        
        Whether the map has been canonicalized
        """
        return self._canonical

    def get_size(self):
        r"""
        Returns the size of the map (the number of edges)
        
        OUTPUT:
        
        Size of the map
        """
        return self._size

    def opposite_corner(self, c):
        r"""
        Returns the corner sharing the same edge.

        INPUT:

        -- ``c``: a corner in the map

        OUTPUT:

        Returns the other corner corresponding to the same edge
        """
        self.__check_corner(c)
        return self._edges[c]

    def next_corner_vertex(self, c):
        r"""
        Returns the next corner on the same vertex in counter-clockwise order.

        INPUT:

        -- ``c``: a corner in the map

        OUTPUT:

        Returns the next corner on the same vertex in counter-clockwise order.
        """
        self.__check_corner(c)
        return self._vertices[c]

    def next_corner_face(self, c):
        r"""
        Returns the next corner on the same face in counter-clockwise order.

        INPUT:

        -- ``c``: a corner in the map

        OUTPUT:

        Returns the next corner on the same face in counter-clockwise order.
        """
        self.__check_corner(c)
        return self._vertices[self._edges[c]]
    
    def corners_around_vertex(self, c):
        r"""
        Returns the corners sharing the same vertex as the given corner.
        
        INPUT:
        
        - ``c``: a given corner
        
        OUTPUT:
        
        A list of corners in counter-clockwise order, starting with ``c``
        """
        self.__check_corner(c)
        l = [c]
        current = self._vertices[c]
        while current != c:
            l.append(current)
            current = self._vertices[current]
        return l
    
    def traversal(self, func):
        r"""
        Traverse the map. The traversal is done in vertex-BFS, starting from the
        root corner.
        
        INPUT:
        
        -- ``func``: a function taking two variables, one is the corner, the
        other is the component number.
        """
        visited = {}
        queue = []
        components = 0
        
        def visit(c, comp):
            func(c, comp)
            visited[c] = True
        
        keys = self._edges.keys()
        if 0 not in keys:
            raise ValueError("Internal error: root corner 0 is absent")
        keys.insert(0, 0)
            
        for key in keys:
            if key not in visited:
                components += 1
                queue = self.corners_around_vertex(key)
                for c in queue:
                    visit(c, components)
            else:
                queue = []
            while queue:
                target = self._edges[queue.pop()]
                if target in visited:
                    continue
                if target in queue:
                    continue
                tovisit = self.corners_around_vertex(target)
                for c in tovisit:
                    visit(c, components)
                queue = tovisit + queue
                
        return components
    
    def __canonicalize(self):
        r"""
        Internal function, to keep corner labels in a canonical form
        """
        if self._canonical:
            return
        relabel = {}
        counter = [0]
        queue = []
        
        def getnewlabel(c, comp):
            relabel[c] = counter[0]
            counter[0] += 1
        
        comps = self.traversal(getnewlabel)

        self._accucnt = counter[0]
        self._size = counter[0] // 2
        newedges = {}
        newvertices = {}
        for key in self._edges:
            newedges[relabel[key]] = relabel[self._edges[key]]
        for key in self._vertices:
            newvertices[relabel[key]] = relabel[self._vertices[key]]
        self._edges = newedges
        self._vertices = newvertices
        self.__update_vertices_inv()
        self._faces = {}
        for c in self._edges:
            self._faces[self._vertices[self._edges[c]]] = c
        self._canonical = True
        self._connected = (1 == comps)
        return
        
    def is_connected(self):
        r"""
        Returns whether the map is connected
        """
        if not self._canonical:
            self.__canonicalize()
        return self._connected
    
    def __rename_label(self, c, cp):
        r"""
        Internal function. Renames a label
        """
        # rename c to cp in value
        self._edges[self._edges[c]] = cp
        cprev = self._vertices_inv[c]
        csucc = self._vertices[c]
        self._vertices[cprev] = cp
        self._vertices_inv[csucc] = cp
        cur = c
        if self.is_canonical():
            while self._faces[cur] != c:
                cur = self._faces[cur]
            self._faces[cur] = cp
        # rename c to cp in domain
        self._edges[cp] = self._edges[c]
        self._vertices[cp] = self._vertices[c]
        self._vertices_inv[cp] = self._vertices_inv[c]
        del self._edges[c]
        del self._vertices[c]
        del self._vertices_inv[c]
        if self.is_canonical():
            self._faces[cp] = self._faces[c]
            del self._faces[c]
        
    def reroot(self, c):
        r"""
        Reroot the map at the given corner.
        
        INPUT:
        
        -- ``c``: The corner that we want to relabel as root
        """
        #if not self.is_canonical():
        #    raise RuntimeError("Cannot reroot non-canonical map")
        self.__check_corner(c)
            
        if 0 == c:
            return
        
        self.__rename_label(c, -1)
        self.__rename_label(0, c)
        self.__rename_label(-1, 0)
        return
        
    def __get_components(self):
        r"""
        Internal function. Return a list of lists of corners in the same
        component
        """
        curcomp = [0]
        comps = [[]]
        
        def adder(c, comp):
            if curcomp[0] < comp:
                curcomp[0] = comp
                comps.append([])
            comps[comp].append(c)
        
        self.traversal(adder)
        return comps[1:]
        
    def get_connected_components(self, rootlist=None):
        r"""
        Returns a list of combinatorial maps that are connected components of
        the current map, with roots given by a list of suggested roots. 0 is
        always considered as a root. The list is ordered by increasing label of
        the roots in the original map. All maps returned are canonical.
    
        INPUT:
    
        -- ``rootlist``: A suggested list of roots. If one components contains
        more than two root candidates, than ties are broken by the label, the
        smaller one will be the root.
        
        OUTPUT:
        
        A list of maps as connected components in the original map, ordered by
        increasing label, with roots suggested by the given list
        """
        if rootlist == None:
            rootlist = [0]
        compsedges = self.__get_components()
        comps = []
        for eset in compsedges:
            edges = {}
            vertices = {}
            for x in eset:
                edges[x] = self._edges[x]
                vertices[x] = self._vertices[x]
                curroots = [x for x in rootlist if x in edges]
                if curroots:
                    root = min(curroots)
                else:
                    root = min(edges)
            comps.append(BijectiveMap.from_data(edges, vertices, root))
        return comps

    def add_vertex(self, c):
        r"""
        Operation to add a new edge linked with a new vertex at the current
        corner. The new corner is inserted before the given one in
        counter-clockwise order. Returns the corner on the new vertex.

        INPUT:
        
        -- ``c``: the corner where to add this new vertex
        
        OUTPUT:
        
        Returns the corner label on the new vertex (of degree 1)
        """
        self.__check_corner(c)
        self.__perform_operation()
        c1 = self._accucnt
        c2 = c1 + 1
        self._accucnt += 2
        self._size += 1
        self._edges[c1] = c2
        self._edges[c2] = c1
        cprec = self._vertices_inv[c]
        self._vertices[cprec] = c1
        self._vertices[c1] = c
        self._vertices_inv[c] = c1
        self._vertices_inv[c1] = cprec
        self._vertices[c2] = c2
        self._vertices_inv[c2] = c2
        return c2

    def smooth_out(self, c):
        r"""
        Operation to remove a vertex of degree 2. Return if the removal is done.
        The reason it is not done is when the given corner is not for a vertex
        of degree 2.
        
        INPUT:
        
        -- ``c``: a corner supposedly of a vertex of degree 2
        
        OUTPUT:
        
        Returns if the removal is performed
        """
        self.__check_corner(c)
        c1 = self._vertices[c]
        if self._vertices[c1] != c:
            return False
        if c == 0:
            cp = self._edges[c]
            self.reroot(cp)
            return self.smooth_out(cp)
        if c1 == 0:
            cp = self._edges[c1]
            self.reroot(cp)
            return self.smooth_out(cp)
        self.__perform_operation()
        self._size -= 1
        cc1 = self._edges[c1]
        cc2 = self._edges[c2]
        self._edges[cc1] = cc2
        self._edges[cc2] = cc1
        del self._edges[c1]
        del self._edges[c2]
        del self._vertices[c1]
        del self._vertices[c2]
        return True

    def join_corners(self, c1, c2):
        r"""
        Operation to join two given corners on different vertices. Returns if
        the operation is done (only obstacle is that the corners being on the
        same vertex.
        
        INPUT:
        -- ``c1``, ``c2``: two corners to join
        
        OUTPUT:
        Returns whether the joining operation is performed.
        """
        self.__check_corner(c1)
        self.__check_corner(c2)
        if c2 in self.corners_around_vertex(c1):
            return False
        self.__perform_operation()
        cc1 = self._vertices_inv[c1]
        cc2 = self._vertices_inv[c2]
        self._vertices[cc1] = c2
        self._vertices_inv[c2] = cc1
        self._vertices[cc2] = c1
        self._vertices_inv[c1] = cc2
        return True
    
    def add_edge(self, c1, c2):
        r"""
        Operation to add a new edge between two given corners. If the corners
        are the same, then a loop is added at this corner. The new corners are
        inserted before the given ones in counter-clockwise order. Returns the
        corners of the new edges.
        
        INPUT:
        
        -- ``c1``, ``c2``: two corners in the map, can be identical
        
        OUTPUT:
        
        Returns the new corners of the new edge, in the order inserted for
        ``c1`` and ``c2``. In the case of identical corners, returns in the
        counter-clockwise order
        """
        self.__check_corner(c1)
        self.__check_corner(c2)
        self.__perform_operation()
        cnew = self.add_vertex(c1)
        self.join_corners(c2, cnew)
        return (self._edges[cnew], cnew)

    def delete_edge(self, c):
        r"""
        Operation to delete a given edge. As the empty map is not allowed in
        this package, error will be thrown if there is only one edge. Returns a
        pair of corners indicating the original positions of the deleted edge.
        One of them may be ``None`` if the corresponding vertex is of degree 1.
        

        INPUT:
        
        -- ``c``: a corner indicating the edge to be deleted

        OUTPUT:

        Returns the current corners of the two ends of the deleted edge
        """
        self.__check_corner(c)
        cp = self._edges[c]
        if 1 == self._size:
            raise ValueError("Cannot delete edge when it is the only one")
        if 0 == c:
            cn = self._vertices[c]
            self.reroot(cn)
            return self.delete_edge(cn)
        if 0 == cp:
            cpn = self._vertices[cp]
            self.reroot(cpn)
            return self.delete_edge(cpn)
        self.__perform_operation()
        self._size -= 1
        del self._edges[c]
        del self._edges[cp]
        retc1, retc2 = None, None
        if self._vertices[c] != c:
            cpred = self._vertices_inv[c]
            csucc = self._vertices[c]
            self._vertices[cpred] = csucc
            self._vertices_inv[csucc] = cpred
            retc1 = csucc
        del self._vertices[c]
        del self._vertices_inv[c]
        if self._vertices[cp] != cp:
            cppred = self._vertices_inv[cp]
            cpsucc = self._vertices[cp]
            self._vertices[cppred] = cpsucc
            self._vertices_inv[cpsucc] = cppred
            retc2 = cpsucc
        del self._vertices[cp]
        del self._vertices_inv[cp]
        # Special case of a loop with neighboring corners
        if retc1 == cp:
            retc1 = retc2
        return (retc1, retc2)

    def contract(self, c):
        r"""
        Operation of contracting the edge given by ``c``. Fails when there is
        only one edge, or when the edge to be contracted is a loop
        
        INPUT:
        
        - ``c``: the corner indicating the edge to contract
        """
        self.__check_corner(c)
        if 1 == self._size:
            raise ValueError("Cannot contract edge when it is the only one.")
        if self._edges[c] in self.corners_around_vertex(c):
            raise ValueError("Cannot contract a loop. Use deletion.")
        self.__perform_operation()
        c1, c2 = self.delete_edge(c)
        if (None == c1) or (None == c2):
            return
        self.join_corners(c1, c2)

    def to_graph(self, embedded=True):
        r"""
        Convert the current map to a graph (optionally with embedding
        information) in sage. Canonicalization will be applied first. Due to
        limitation in Sagemath, multiple edges are not supported when
        embedding is needed, and an error will be raised.

        INPUT:

        -- ``embedded``: indicating whether embedding information is added

        OUTPUT:

        A graph object in Sagemath, if succeeded

        NOTE:
        
        Currently the embedding information in sage does not allow us to
        distinguish different edges in multiedge case. Nothing is guaranteed
        here. Judging from the implementation in Sagemath, an improvement may
        not be possible in the near future.
        """
        self.__canonicalize()
        vlist = []
        cornerv = {}
        vcnt = 0
        for c in self._edges:
            if c not in cornerv:
                curvlist = self.corners_around_vertex(c)
                vlist.append(curvlist)
                for x in curvlist:
                    cornerv[x] = vcnt
                vcnt += 1
        gdict = {}
        donec = {}
        for vn in range(len(vlist)):
            edict = {}
            for c in vlist[vn]:
                if c in donec:
                    continue
                cp = self._edges[c]
                cpv = cornerv[cp]
                if cpv not in edict:
                    edict[cpv] = [(c, cp)]
                else:
                    edict[cpv].append((c, cp))
                donec[c] = True
                donec[cp] = True
            gdict[vn] = edict
        g = Graph(gdict, format='dict_of_dicts')
        if embedded:
            # check multiedges
            for v1 in gdict:
                for v2 in gdict[v1]:
                    if len(gdict[v1][v2]) != 1:
                        raise ValueError("Multiedge exists. Not for Sagemath.")
            embdict = {}
            for vn in range(len(vlist)):
                embdict[vn] = map(lambda c: cornerv[self._edges[c]], vlist[vn])
            g.set_embedding(embdict)
        return g
    
    def plot(self):
        r"""
        Returns a plot of the current map. Plotting ability is limited by that
        of Sagemath.
        """
        return self.to_graph().plot(layout="planar", edge_labels=True)
    
