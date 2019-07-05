# BijectiveMap
A Sagemath class for bijective operations of combinatorial maps

It is based on rotation systems. Every operation is done by corners. When an edge is drawed horizontaly, the related corners are those on the lower-left and upper-right. The root corner is always 0.

The following functions are provided:

- `from_data`: construct a map with dictionaries of edge and vertex permutations

- `opposite_corner`: returns the other corner associated to the given one

- `next_corner_vertex`: returns the next corner of the same vertex in counter-clockwise order

- `next_corner_face`: returns the next corner of the same face in counter-clockwise order

- `corners_around_vertex`: returns a list of corners of the same vertex, in counter-clockwise order

- `traversal`: performs traversal with a given function, which can read the current corner and the component number of it

- `is_connected`: returns whether the map is connected

- `reroot`: reroots to the given corner

- `get_connected_components`: returns a list of BijectiveMap corresponding to components

- `add_vertex`: adds a new vertex with a new edges linking to a given corner. The corner on the new vertex is returned, and the other new corner is inserted before the given one in counter-clockwise order.

- `smooth_out`: removes a vertex of degree 2

- `join_corners`: joins two given corners (not on the same vertex)

- `add_edges`: adds an edge between two corners, returns the new corners

- `delete_edges`: deletes an edge, returns the current corners of the two ends of the deleted edge, can be None in the case of degree 1

- `contract`: contracts a non-loop edge

- `to_graph`: convert to a graph object with embeddings in Sagemath

- `plot`: plots the map with Sagemath's functionality. Due to the limitation of Sagemath, it does not work on non-planar maps or disconnected maps.