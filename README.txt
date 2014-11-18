Team:
Tian Lan
Xide Xia

Description:
We construct Mesh basically using the design we started for Assignment 4A. We template Mesh on N, E, T so that the Mesh supports user-specified values: node_value_type, edge_value_type, and triangle_value_type. These user-specified values will be stored inside the internal_node_value, internal_edge_value, and internal_triangle_value respectively. We use a graph and a vector of internal triangles to represent a mesh. The graph stores all internal nodes and internal edges values. We construct a new Node class and a new Edge class in Mesh. In this way, users will only access to the user-specified values by calling Mesh::Node.value()/ Mesh::Edge.value(). The triangle vector stores information about the unique ids of nodes and edges that form a triangle, and the value associated with the triangle. We also make some slight modifications such as using triangle.node(size_type i) to access one of the three nodes of a triangle instead of node1(), node2(), and node3().

User can run the mesh with different data directly without choosing the initial conditions. The system will recognize the input pattern and choose the corresponding initial conditions.


