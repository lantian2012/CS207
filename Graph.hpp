#ifndef CS207_GRAPH_HPP
#define CS207_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

#include "CS207/Util.hpp"
#include "Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  
  struct internal_node;  //for storing nodes in graph
  
  
  

 public:

  /////////////////////////////
  // PUBLIC TYPE DEFINITIONS //
  /////////////////////////////

  /** Type of this graph. */
  typedef Graph graph_type;
  
  /** Synonym for the type of value stored in Node */
  typedef V node_value_type;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  typedef Node node_type;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  typedef Edge edge_type;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  typedef unsigned size_type;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  typedef NodeIterator node_iterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  typedef EdgeIterator edge_iterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  typedef IncidentIterator incident_iterator;

  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    size_ = 0;
    edgesize_ = 0;
    
  }
  /** Default destructor */
  ~Graph() = default;

  /////////////
  // General //
  /////////////

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return size_;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes.clear();
    size_ = 0;
    edgesize_ = 0;
  }

  /////////////////
  // GRAPH NODES //
  /////////////////

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node:private totally_ordered<Node> {
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    Node() {
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return fetch().point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return idx_;
    }

    /** Test whether this node and @a x are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& x) const {
      // HW0: YOUR CODE HERE
      if (x.graph_ == graph_ && x.idx_ == idx_)
        return true;
      //(void) x;          // Quiet compiler warning
      return false;
    }

    /** Test whether this node is less than @a x in the global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& x) const {
      // HW0: YOUR CODE HERE
      return (idx_ < x.idx_);
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Get the value associated with this Node
     */
    node_value_type& value(){
      return fetch().value;
    }

    /** Get the value associated with this Node
     */
    const node_value_type& value() const {
      return fetch().value;
    }

    // Get the size 
    size_type degree() const{
      return fetch().neighbors.size();
    }
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_, idx_, 0);
    }
    incident_iterator edge_end() const{
      return IncidentIterator(graph_, idx_, degree());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    graph_type* graph_;  //pointer to the associated graph
    size_type idx_;      //index of the node
    // a private constructor for graph to construct a Node instance
    Node(const graph_type* graph, size_type idx)
      : graph_(const_cast<graph_type*>(graph)), idx_(idx){
      }
    /**
     Fetch internal node stored in graph
     */
    internal_node& fetch() const {
      if (idx_ >= graph_->size())
        assert(false);
      return graph_->nodes[idx_];
    }
  };

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new size() == old size() + 1
   * @post result_node.index() == old size()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value_in = node_value_type()) {
    // HW0: YOUR CODE HERE
    internal_node new_node;
    new_node.idx = size_;
    new_node.point = position;
    new_node.value = value_in;
    nodes.push_back(new_node);
    size_ ++;
    return Node(this, size_-1); 
  }

  /** Determine if this Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW1: YOUR CODE HERE
    
    return (n < size());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return Node(this, i);
  }

  /////////////////
  // GRAPH EDGES //
  /////////////////

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge:private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      // fetch the node index stored in the graph first
      // then fetch the node from graph
      return graph_->node(node1_);      
       
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      // fetch the node index stored in the graph first
      // then fetch the node from graph
      return graph_->node(node2_);
    }

    /** Test whether this edge and @a x are equal.
     *
     * Equal edges are from the same graph and have the same nodes.
     */
    bool operator==(const Edge& x) const {
      // HW0: YOUR CODE HERE
      //(void) x;          // Quiet compiler warning
      if (graph_ == x.graph_ && node1_ == x.node1_ && node2_ == x.node2_)
        return true;
      return false;
    }

    /** Test whether this edge is less than @a x in the global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The edge ordering relation must obey trichotomy: For any two edges x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Edge& x) const {
      // HW0: YOUR CODE HERE
      return (node1_ < x.node1_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    graph_type* graph_;  //pointer to the associated graph
    size_type node1_;      //the uid of node 1  
    size_type node2_;      //the uid of node 2  
    
    //private constructor for graph to construct edge instance
    Edge(const graph_type* graph, size_type node1, size_type node2)
      : graph_(const_cast<graph_type*>(graph)), node1_(node1), node2_(node2){
      }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edgesize_;
  }
  Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    size_type node1_uid = a.index();
    size_type node2_uid = b.index();
    if (has_edge(a, b)){
      if (node1_uid < node2_uid)
        return Edge(this, node1_uid, node2_uid);
      else
        return Edge(this, node2_uid, node1_uid);
    }
    nodes[node1_uid].neighbors.push_back(node2_uid);
    nodes[node2_uid].neighbors.push_back(node1_uid);
    edgesize_++;
    if (node1_uid < node2_uid)
      return Edge(this, node1_uid, node2_uid);
    else
      return Edge(this, node2_uid, node1_uid);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return true if, for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW1: YOUR CODE HERE
    size_type node1_uid = a.index();
    size_type node2_uid = b.index();
    auto it = std::find(nodes[node1_uid].neighbors.begin(), nodes[node1_uid].neighbors.end(), node2_uid);
    if (it == nodes[node1_uid].neighbors.end())
      return false;
    return true;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    unsigned int count = 0;
    i = i+1;
    unsigned int j = 0;
    unsigned int k = 0;
    while(count < i){
      if (k < nodes[j].neighbors.size() && nodes[j].neighbors[k] > j){
        count++;
      }
      k++;
      if (k >= nodes[j].neighbors.size()){
        k = 0;
        j++;
      }
    }
    if (k == 0){
      j--;
      return Edge(this, j, nodes[j].neighbors.back());
    }
    else{
      k--;
      return Edge(this, j, nodes[j].neighbors[k]);
    }
  }

  ///////////////
  // Iterators //
  ///////////////

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator:private totally_ordered<NodeIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Node value_type;
    /** Type of pointers to elements. */
    typedef Node* pointer;
    /** Type of references to elements. */
    typedef Node& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    //return the Node the interator is pointing to
    Node operator*() const{
      return graph_->node(uid_);
    }

    //return the NodeIterator that points to the
    //next position
    //@pre self.uid_ < graph_->size()
    NodeIterator& operator++(){
      uid_++;
      return *this;
    }
    //True if pointing to the same node in the same graph
    bool operator==(const NodeIterator& x) const{
      if (graph_ == x.graph_ && uid_ == x.uid_)
        return true;
      return false;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    graph_type* graph_;
    size_type uid_;
    
 
    //private constructor for graph to construct NodeIterator instance
    NodeIterator(const graph_type* graph, size_type uid)
      : graph_(const_cast<graph_type*>(graph)), uid_(uid){
      }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  //return the node_iterator that points to the first node
  //if no nodes exist, still points to the initial position
  node_iterator node_begin() const{
    return NodeIterator(this, 0);
  }
  //return the node_iterator that points to the position
  //after the last node
  //if no nodes exist, points to the initial position
  node_iterator node_end() const{
    return NodeIterator(this, nodes.size());
  }

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator:private totally_ordered<EdgeIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:


    Edge operator*() const{
      return Edge(graph_, node_, graph_->nodes[node_].neighbors[nbidx_]);
    }


    EdgeIterator& operator++(){
      //cycles until next edge is found or reaches the end
      do{
        if (nbidx_ < (graph_->nodes[node_].neighbors.size()-1)){
          nbidx_++;
        }
        else{
          nbidx_ = 0;
          node_++;
          while(node_ != graph_->nodes.size() && graph_->nodes[node_].neighbors.size() == 0)
           node_++;
        }
      } while (node_ != graph_->nodes.size() && graph_->nodes[node_].neighbors[nbidx_] < node_); 
      return *this;
    }
    
    bool operator==(const EdgeIterator& x) const{
      return (graph_ == x.graph_ && node_ == x.node_ && nbidx_ == x.nbidx_);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    graph_type* graph_;
    size_type node_;
    size_type nbidx_;
    //private constructor for graph to construct NodeIterator instance
    EdgeIterator(const graph_type* graph, size_type node, size_type nbidx)
      : graph_(const_cast<graph_type*>(graph)), node_(node), nbidx_(nbidx){
      }
    
  };

  // HW1 #3: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const{
    unsigned int i = 0;
    while (i < nodes.size()){
      if (nodes[i].neighbors.size() != 0)
        return EdgeIterator(this, i, 0);
    }
    return EdgeIterator(this, 0, 0);
  }
  edge_iterator edge_end() const{
    unsigned int j = 0;
    for (unsigned int i=0; i<nodes.size(); ++i){
      if (nodes[i].neighbors.size() != 0)
        j = i+1;
    }
    return EdgeIterator(this, j, 0);
  }


  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator:private totally_ordered<IncidentIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    //return the Edge pointing to
    //@post result.node1().index() = uid_;  result.node2().index()=idx_
    Edge operator*() const{
      return Edge(graph_, uid_, graph_->nodes[uid_].neighbors[idx_]);
    }
    IncidentIterator& operator++(){
      if (idx_ == graph_->nodes[uid_].neighbors.size())
        return *this;
      else{
        idx_++;
        return *this;
      }
    }
    bool operator==(const IncidentIterator& x) const{
      return (graph_ == x.graph_ && uid_ == x.uid_ && idx_ == x.idx_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_;
    size_type uid_;    //the uid of the Node that spawns the incident_iterator
    size_type idx_;    //the index of the edge pointing to
    //private constructor for graph to construct NodeIterator instance
    IncidentIterator(const graph_type* graph, size_type uid, size_type idx)
      : graph_(const_cast<graph_type*>(graph)), uid_(uid), idx_(idx){
      }
  };

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  
  struct internal_node {
    Point point;
    size_type idx;
    std::vector<size_type> neighbors;   //the uid of adjacent nodes
    node_value_type value;        //the additional value stored in nodes
  };

  
  std::vector<internal_node> nodes;  //store node information in graph
  size_type size_;
  size_type edgesize_;

};

#endif
