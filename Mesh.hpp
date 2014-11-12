#pragma once
#include "Graph.hpp"
#include "Point.hpp"
/** @file Mesh.hpp
 * @brief A Mesh is composed of nodes, edges, and triangles such that:
 *  -- All triangles have three nodes and three edges.
 *  -- All edges belong to at least one triangle and at most two triangles.
 */

 /** Water column characteristics */
struct QVar {
  double h;   // Height of column
  double hx;  // Height times average x velocity of column
  double hy;  // Height times average y velocity of column

  /** Default constructor.
   *
   * A default water column is 1 unit high with no velocity. */
  QVar()
      : h(1), hx(0), hy(0) {
  }
  /** Construct the given water column. */
  QVar(double h_, double hx_, double hy_)
      : h(h_), hx(hx_), hy(hy_) {
  }
  // More operators?
};

/** @class Mesh
 * @brief A template for 3D triangular meshes.
 *
 * Users can add triangles and retrieve nodes, edges, and triangles.
 */
template <typename N, typename E, typename T>
class Mesh {
  // HW3: YOUR CODE HERE
  // Write all typedefs, public member functions, private member functions,
  //   inner classes, private member data, etc.hello




 public:
  /** Type of indexes and sizes. Return type of Mesh::num_nodes(). */
  typedef unsigned size_type;
  typedef Graph<N, E> GraphType;
  typedef typename GraphType::node_type Node;
  typedef typename GraphType::edge_type Edge;
  typedef typename GraphType::node_type node_type;
  typedef typename GraphType::node_iterator node_iterator;
  typedef typename GraphType::edge_iterator edge_iterator;
  /** Predeclaration of Triangle type. */
  class Triangle;
  /** Type of triangle iterators, which iterate over all adjacent mesh triangles of a triangle. */
  class IncidentIterator_Triangle;  
  /** Synonym for IncidentIterator_Triangle */
  typedef IncidentIterator_Triangle incidentiterator_triangle;

  Mesh() {}

  /** Return the number of nodes in the mesh. */
  size_type num_nodes() const {
    return graph_.num_nodes();
  }

  /** Return the number of edges in the mesh. */
  size_type num_edges() const {
    return graph_.num_edges();
  }

  /** Return the number of triangles in the mesh. */
  size_type num_triangles() const {
    return tri_vec.size();
  }

  class Triangle
  {
  public:
    /**construc an invalid triangle*/
    Triangle(){
    }
    /**Access the i node of the triangle
     * @pre 0<=i<3
     * @post result.index() == mesh_->tri_vec[uid_].nodes[i]
     * Complexity = O(1)
     */
    Node node(size_type i) {
      size_type node_uid = mesh_->tri_vec[uid_].nodes[i];
      return mesh_->graph_.node(node_uid);
    }
    /**Access the first node of the triangle
     * @pre 0<=i<3
     * @post result.index() == mesh_->tri_vec[uid_].edges[i]
     * Complexity = O(1)
     */
    Edge edge(size_type i) {
      size_type edge_uid = mesh_->tri_vec[uid_].edges[i];
      return mesh_->graph_.edge(edge_uid);
    }

    /**return the area of this triangle
     *Comlexity = O(1)
     */
    double area() const{
      double x0 = mesh_->tri_vec[uid_].nodes[0].position().x;
      double y0 = mesh_->tri_vec[uid_].nodes[0].position().y;
      double x1 = mesh_->tri_vec[uid_].nodes[1].position().x;
      double y1 = mesh_->tri_vec[uid_].nodes[1].position().y;
      double x2 = mesh_->tri_vec[uid_].nodes[2].position().x;
      double y2 = mesh_->tri_vec[uid_].nodes[2].position().y;
      return (x0*y1+x1*y2+x2*y0-x1*y0-x2*y1-x0*y2)/2;
    }
    /**Access the Q of the triangle
     *Complexity = O(1)
     */
    QVar Q(){
      return mesh_->tri_vec[uid_].Q;
    }
    /**return the index of the triangle
     *Complexity = O(1)
     */
    size_type index(){
      return uid_;
    }

    /**return the unit normal vector of the edge
        *that is opposite to node i
      */
    Point normal(size_type i){
      return mesh_->tri_vec[uid_].n[i];
    }

    /**return the flux of the edge
        *that is opposite to node i
      */
    QVar F(size_type i){
      return mesh_->tri_vec[uid_].F[i];
    }



    /* Return an iterator points to the first adjacent triangle of this triangle.
     * @post result->incident_i ==0
     * @post result->mesh_ == this->mesh_
     * @post result->triangle_uid == this->uid_
     *
     *Complexity = O(1)
         */
    incidentiterator_triangle triangle_begin() const {
      return incidentiterator_triangle(this, uid_, 0);
    }
  
        /* Return an iterator points to the last adjacent triangle of this triangle.
     * @post result->incident_i == this->mesh_->edges.size()
     * @post result->mesh_ == this->mesh_
     * @post result->triangle_uid == this->uid_
     *
     *Complexity = O(1)
         */
    incidentiterator_triangle triangle_end() const {
      return incidentiterator_triangle(this, uid_, 2);
    }

  private:
    friend class Mesh;
    Triangle(Mesh* mesh, size_type uid):
      mesh_(mesh), uid_(uid) {}
    size_type uid_;
    Mesh* mesh_;
  };

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value_in The value associated with this node
   * @post new graph_.size() == old graph_.size() + 1
   * @post result_node.index() == old size()
   * Complexity: O(1)
   */
  Node add_node(const Point& position, const N& value_in = N()) {
    return graph_.add_node(position, value_in);
  }

  /** Add an triangle to the graph
   *@param[in] a, b, c are all valid nodes in the trigraph
   *@return The triangle that has been added
   *@post new size() = old size() +1
   * Complexity: O(d) d is the degree of a node
   */
  Triangle add_triangle(const Node& a, const Node& b, const Node& c){
    T new_triangle;
    //add edge in graph
    graph_.add_edge(a, b);
    new_triangle.nodes[0] = c.index();

    //calculate normal for new edges
    //update triangle id in edge data
    //pushback triange data
    (void) a, b, c;
    return Triangle();
  }
/*
Edge add_edge(const Node& a, const Node& b) {
    size_type node1_uid = a.uid_;
    size_type node2_uid = b.uid_;
    //check if edge exists
    if (has_edge(a, b)){
      if (node1_uid < node2_uid)
        return Edge(this, node1_uid, node2_uid);
      else
        return Edge(this, node2_uid, node1_uid);
    }
    //if not, add a new edge
    nodes[node1_uid].neighbors.push_back(node2_uid);
    nodes[node1_uid].edgevalues.push_back(edges.size());
    nodes[node2_uid].neighbors.push_back(node1_uid);
    nodes[node2_uid].edgevalues.push_back(edges.size());
    edges.push_back(internal_edge());
    //update edge size
    edgesize_++;
    if (node1_uid < node2_uid)
      return Edge(this, node1_uid, node2_uid);
    else
      return Edge(this, node2_uid, node1_uid);
  }
*/










  /** Determine if this Triangle belongs to this Graph
    * @return True if @a t is currently a Triangle of this Graph
    *
    * Complexity: O(1).
    */
  bool has_triangle(const Triangle& t){
    //compare the trigraph_ and uid_
    (void) t;
    return false;
  }

    /** Return the triangle with index @a i.
     * @pre 0 <= @a i < size()
     * @post result.index() == i
     * Complexity: O(1).
     */
  Triangle triange(size_type i){
    (void) i;
    return Triangle();
  }
  
  /** Return the first triangle adjacent to the input edge.
     * @param[in] e, an edge inside the mesh
     * @pre @a e.index() < graph_->num_edges()
     * @post result.uid_ == graph_->get_edge_value(edge_uid).F_triangle1.
     * Complexity: O(1).
     */
  Triangle triangle1(const Edge& e){
    (void) e;
    return Triangle();
  }

  /** Return the second triangle adjacent to the input edge.
     * @param[in] e, an edge inside the mesh
     * @pre @a e.index() < graph_->num_edges()
     * @post result.uid_ == graph_->get_edge_value(edge_uid).F_triangle2.
     * Complexity: O(1).
     */
  Triangle triangle2(const Edge& e){
    (void) e;
    return Triangle();
  }


  /* Return an iterator points to the first node in graph_.
     * Complexity: O(1).
     */
  node_iterator  node_begin(){
    return graph_.node_begin();
  }

    /* Return an iterator points to the last node in graph_.
     * Complexity: O(1).
     */
  node_iterator node_end(){
    return graph_.node_end();
  }
    
    /* Return an iterator points to the first edge in graph_.
     * Complexity: O(1).
     */
  edge_iterator edge_begin(){
    return graph_.edge_begin();
  }

    /* Return an iterator points to the last edge in graph_.
     * Complexity: O(1).
     */
  edge_iterator edge_end(){
    return graph_.edge_end();
  }


  //Iterator that iterates through the adjacent triangle of a node
  class IncidentIterator_Node
  {
  public:
    // These type definitions help us use STL's iterator_traits.
      /** Element type. */
      typedef Triangle value_type;
      /** Type of pointers to elements. */
      typedef Triangle* pointer;
      /** Type of references to elements. */
      typedef Triangle& reference;
      /** Iterator category. */
      typedef std::input_iterator_tag iterator_category;
      /** Difference between iterators */
      typedef std::ptrdiff_t difference_type;
    

    /** Construct an invalid IncidentIterator_Node*/
    IncidentIterator_Node(){
    }

    /* Return the Triangle pointed by the IncidentIterator_Node
     * @post return.uid_ = t_
     */
    Triangle operator*() const{
      return mesh_->triangle(t_);
    }

    /* return the IncidentIterator_Node that points to the
     * next incident triangle
     * @pre t_ old != -1
     * @post If t_ old has two adjacent triangles, init_ != last_ 
       * t_ new points to the adjacent triangle which does not 
       * contain node last_.
       * If t_ old has two adjacent triangles, init_ == last_, t_ new = -1
       * If t_ old has one adjacent triangles, t2_!=-1, t_ new = t2_
       * If t_ old has one adjacent triangles, t2_ == -1, t_ new = -1
     */
    IncidentIterator_Node& operator++(){
      //get an incident edge of the node
      //record the last_ = first_ = incident_edge.node().index()
      //t_ = incident_edge.value().triangle1 t2_ = incident_edge.value().triangle2
      //get the new t_ and update last_
      //if last_ == init_ set t_ to invalid to mark the end
      //if meets the border. set last_ = init_ and explore t2_
      //when t2_ meets the end, set t_ to invilid
      return IncidentIterator_Node();
      }

      //True if n_ and t_ are the same
      bool operator==(const IncidentIterator_Node& x) const{
        (void) x;
        return false;
      }
  private:
    friend class Mesh;
    Mesh* mesh_;
    size_type n_; //The uid of the node in the center
    size_type t_; //The uid of the triangle pointed to
    size_type t2_; //the uid of the triangle to be explored when t_ reaches boundary
    size_type last_; //the uid of the other node of the edge just discovered
    size_type init_; //the uid of the other node of the first edge

    /** Construct a valid IncidentIterator_Node*/
    IncidentIterator_Node( const Mesh* m, size_type n, size_type t, size_type t2, size_type last, size_type init): mesh_(m),n_(n),t_(t),t2_(t2),last_(last),init_(init){
    }
  };



  /* Get the iterator that points to the first adjacent triangle of a node
   * @param n The node in the center
   * @return the iterator that points to the first adjacent triangle of @a n
   * @post result.n_ == n.uid_
   * @post result.t_ result.t2_ are the adjacent triangles of edge(n_, result.last_)
   * @post result.last_ == result.first_
   */
  IncidentIterator_Node triangle_begin(Node n){
    (void) n;
    return IncidentIterator_Node();
  }

  /* Get the iterator that points to an invalid triangle of a node
   * @param n The node in the center
   * @return the iterator that points to an invalid triangle of @a n
   * @post result.n_ == n.uid_
   * @post result.t_ == -1
   */
  IncidentIterator_Node triangle_end(Node n){
    (void) n;
    return IncidentIterator_Node();
  }


  //Iterators:
  //all triangles
  //adjacent triangles of a node
  //adjacent triangles of a triangle

  //Iterator that iterates through the adjacent triangle of a triangle
  /** @class Mesh::IncidentIterator_Triangle
     * @brief Interator Class for triangles incident to a triangle. A forward iterator.
     * @RI graph_pointer != nullptr && incident_i <= 2
     */
  class IncidentIterator_Triangle
  {
  public:
    // These type definitions help us use STL's iterator_traits.
      /** Element type. */
      typedef Triangle value_type;
      /** Type of pointers to elements. */
      typedef Triangle* pointer;
      /** Type of references to elements. */
      typedef Triangle& reference;
      /** Iterator category. */
      typedef std::input_iterator_tag iterator_category;
      /** Difference between iterators */
      typedef std::ptrdiff_t difference_type;
    

    /** Construct an invalid IncidentIterator_Triangle*/
      IncidentIterator_Triangle(){
      }

      /**return the Triangle pointed by the IncidentIterator_Triangle
       *@pre this != nullptr
       *@post result is a triange adjacent to Triangle(this->mesh_,this->triangle_uid)
       *@post result is one of the two adjacent triangles to the ith edge of Triangle(this->mesh_,this->triangle_uid)  
       *@post result.uid_ != this->triangle_uid
       *
       *Complexity = O(1)
       */
      Triangle operator*() const{
        // size_type edge_uid = mesh_->edges[incident_i];
        // pair<size_type,size_type> pair = graph_->get_edge_value(edge_uid);
        // size_type t_uid;
        // if(F_triange1 == triangle_uid)
        //   t_uid = F_triange2;
        // else
        //   t_uid = F_triange1;
        // return Triangle(mesh_,t_uid);

        return Triangle();
      }

      /*return the IncidentIterator_Triangle that points to the next incident edge
       *@pre this->incident_i < 2
       *@post result.incident_i > this->incident_i
       *@post result.incident_i <=2
       *
       *Complexity = O(1)
       */
      IncidentIterator_Triangle& operator++(){
        return IncidentIterator_Triangle();
      }

        /*return true if the two triangles are same and false if not
         *@param[in] x, reference to a triangle 
       *@pre @a x.uid_ >=0
       *@pre triangle_uid >=0
       *@post result==true if @a x.uid_ == this->triangle_uid
       *@post result==false if @a x.uid_ != this->triangle_uid
       *
       *Complexity = O(1)
       */
      bool operator==(const IncidentIterator_Triangle& x) const{
        (void) x;
        return false;
      }
  private:
    
    Mesh* mesh_;  //poiter pointing to this Mesh
    size_type triangle_uid;  //the uid of this triangle
    /** incident_i represents the ith edge of this triangle*/
    size_type incident_i;

        /** Construct a valid IncidentIterator_Triangle*/
        IncidentIterator_Triangle(const Mesh* m, size_type uid, size_type i): mesh_(const_cast<Mesh*>(m)),triangle_uid(uid),incident_i(i){
      }
  };


private:
  GraphType graph_; //a Graph Class object storing all the nodes and edges
  std::vector<T> tri_vec; // a vector storing all the triangles' data
};
