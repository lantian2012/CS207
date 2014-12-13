/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"
#include "Viewer_Extension.hpp"
#include "collision.hpp"

#include "Mesh.hpp"
#include "Point.hpp"
#include "CollisionDetector.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point velocity;  //< Node velocity
  double mass;     //< Node mass
  QVar Q;
  NodeData():Q(0, 0, 0){}
};

/** Custom structure of data to store with Edges */
struct EdgeData{
  double L;
  double K;
};

/** Custom structure of data to store with Triangles */
struct TriData
{
  Point n; //the outward surface normal vector
  QVar Q;  //Qk: the average value of Q inside this triangle
  std::vector<QVar> F; //The transition over an edge
  TriData():Q(QVar()), F(3, QVar()){}
};


typedef Mesh<NodeData, EdgeData, TriData> MeshType;
typedef typename MeshType::node_type Node;
typedef typename MeshType::edge_type Edge;


// This is one of the example object that can be floating on the water
// It should contain weights, length, its location and speed available
// Again this is just a simple example of floating objects which can not be shown in viewer
// Main fuctionality is generalized using template in the latter part for genaric usage
// Small ship object used in the game provided
struct Ship {
  double Weight;
  double length;
  double center_x;
  double center_y;
  double speed_x;
  double speed_y;
  double Area;

  /** Default constructor for ship
   *  Create a generic ship that has weight, length, center_x, center_y, speed_x, speed_y set as default value
   *  Area can be changed and is calculated depending on the user need
   **/
  Ship(): Weight(1.5), length(1), center_x(0.0), center_y(0.0), speed_x(0), speed_y(0) {
    double e = length/2;
    Area = 3.1416 * e * e;
  }

  /** Constructor for ship
   * @param[in] double Weight
   * @param[in] double length
   * @param[in] double center_x
   * @param[in] double center_y
   * @param[in] double speed_x
   * @param[in] double speed_y
   *
   * Create a generic ship that has weight, length, center_x, center_y, speed_x, speed_y
   * Area can be changed and is calculated depending on the user need
   * 
   **/
  Ship(double w, double a, double c_x, double c_y, double s_x, double s_y)
      : Weight(w), length(a), center_x(c_x), center_y(c_y), speed_x(s_x), speed_y(s_y) {
    double e = length/2;
    Area = 3.1416 * e * e;
  }

  /** Check if a point is under the ship area
   * @param[in] Point a point position to check
   * @param[out] bool true if it is under the ship, false if it is not
   * 
   * This function check if a point is under the ship. 
   * This can be changed based on the shape of the ship and direction of the ship
   */
  bool checkCoverNode(Point p) {
    if (norm(p - Point(center_x, center_y, 0))
       < sqrt(Area/3.1416)) {
       return true;
    } else {
      return false;
    }
  }

  /** Move the ship location in time dt
   * @param[in] double dt small change in time period
   * 
   * This fucntion changes the ship location in short time dt based on 
   * current speed current location
   */
  void move(double dt) {
    center_x += speed_x * dt;
    center_y += speed_y * dt;
  }
};

/** Source calculcator of partial derivative of x by varying position x and y
   * @param[in] double x
   * @param[in] double y
   * @param[out] double partial derivative of x at position x and y
   *
   * This function can be changed to tailer different shapes of source
   */
double dx_value(double x, double y){
  if (x > 0) {
    x = -x;
  }
  return 0.2*0.12*2*x+0*y;
};

/** Source calculcator of partial derivative of y by varying position x and y
   * @param[in] double x
   * @param[in] double y
   * @param[out] double partial derivative of y at position x and y
   *
   * This function can be changed to tailer different shapes of source
   */
double dy_value(double x, double y){
  if (y > 0) {
    y = -y;
  }
  return 0.2*0.15*2*y+0*x;
};


template <typename OBJ>
struct EdgeFluxCalculator {
  QVar operator()(double nx, double ny, double dt,
                  const QVar& qk, const QVar& qm, Point p1, Point p2, std::vector<OBJ>& obj_vector) {
    // Normalize the (nx,ny) vector
    double n_length = std::sqrt(nx*nx + ny*ny);
    nx /= n_length;
    ny /= n_length;

    // The velocities normal to the edge
    double wm = (qm.hx*nx + qm.hy*ny) / qm.h;
    double wk = (qk.hx*nx + qk.hy*ny) / qk.h;

    // Lax-Wendroff local dissipation coefficient
    double vm = sqrt(grav*qm.h) + sqrt(qm.hx*qm.hx + qm.hy*qm.hy) / qm.h;
    double vk = sqrt(grav*qk.h) + sqrt(qk.hx*qk.hx + qk.hy*qk.hy) / qk.h;
    double a  = dt * std::max(vm*vm, vk*vk);

    // Helper values
    double scale = 0.5 * n_length;
    double gh2   = 0.5 * grav * (qm.h*qm.h + qk.h*qk.h);

    for (unsigned i = 0; i < obj_vector.size(); ++i) {
      if (obj_vector[i].checkCoverNode(p1)) {
        gh2 = gh2  +  obj_vector[i].Weight*qm.h/1/obj_vector[i].Area;
      }
      if (obj_vector[i].checkCoverNode(p2)) {
        gh2 = gh2 + obj_vector[i].Weight*qk.h/1/obj_vector[i].Area;
      }
    }
    // Simple flux with dissipation for stability
    return QVar(scale * (wm*qm.h  + wk*qk.h)           - a * (qm.h  - qk.h),
                scale * (wm*qm.hx + wk*qk.hx + gh2*nx) - a * (qm.hx - qk.hx),
                scale * (wm*qm.hy + wk*qk.hy + gh2*ny) - a * (qm.hy - qk.hy));
  }
};

/** Node position function object for use in the SDLViewer. */
struct NodePosition {
  template <typename NODE>
  Point operator()(const NODE& n) {
    // Change this to plot something other than the
    // positions of the nodes
    return Point(n.position().x, n.position().y, n.value().Q.h);
  }
};

/** Integrate a hyperbolic conservation law defined over the mesh m
 * with flux functor f by dt in time.
 */
template <typename MESH, typename FLUX, typename OBJ>
double hyperbolic_step2(MESH& m, FLUX& f, double t, double dt, std::vector<OBJ>& obj_vector) {
  // Step the finite volume model in time by dt.
  // Pseudocode:
  // Compute all fluxes. (before updating any triangle Q_bars)
  // For each triangle, update Q_bar using the fluxes as in Equation 8.
  //  NOTE: Much like symp_euler_step, this may require TWO for-loops
  for (auto i = m.edge_begin(); i != m.edge_end(); ++i) {

    if ((*i).triangle1().index() != (unsigned) -1 && (*i).triangle2().index() != (unsigned) -1) {
      MeshType::Triangle trik = (*i).triangle1();
      MeshType::Triangle trim = (*i).triangle2();
      unsigned int edge_k = 0;
      unsigned int edge_m = 0;
      //which edge (*i) is in trik and trim
      while(trik.node(edge_k).index()== (*i).node1().index() 
        || trik.node(edge_k).index()== (*i).node2().index() )
        ++edge_k;
      while(trim.node(edge_m).index()== (*i).node1().index() 
        || trim.node(edge_m).index()== (*i).node2().index() )
        ++edge_m;
      QVar flux = f(trik.normal(edge_k).x, trik.normal(edge_k).y, dt, trik.value().Q, trim.value().Q, (*i).node1().position(), (*i).node2().position(), obj_vector);
      trik.value().F[edge_k] = flux;
      trim.value().F[edge_m] = -flux;
    } else {
      MeshType::Triangle trik;
      if ((*i).triangle1().index() != (unsigned) -1)
        trik = (*i).triangle1();
      else
        trik = (*i).triangle2();
      unsigned int edge_k = 0;
      while(trik.node(edge_k).index()== (*i).node1().index() 
        || trik.node(edge_k).index()== (*i).node2().index() )
        ++edge_k;
      QVar flux = f(trik.normal(edge_k).x, trik.normal(edge_k).y, dt, trik.value().Q, QVar(trik.value().Q.h, 0, 0), (*i).node1().position(), (*i).node2().position(), obj_vector);
      trik.value().F[edge_k] = flux;
    }
  }

  for(auto i = m.triangle_begin(); i != m.triangle_end(); ++i){
    QVar sum = QVar(0, 0, 0);
    Point center = Point(0,0,0);
    for (int j = 0; j < 3; ++j){
      sum += (*i).value().F[j];
      center += (*i).node(j).position();
    }
    center = center/3;
    (*i).value().Q = (*i).value().Q-dt/(*i).area()*sum + 
      dt * QVar(0, -grav*(*i).value().Q.h*dx_value(center.x,center.y), 
                   -grav*(*i).value().Q.h*dy_value(center.x,center.y));
  }
  return t + dt;
}

/** Convert the triangle-averaged values to node-averaged values for viewing. */
template <typename MESH>
void post_process(MESH& m) {
  // Translate the triangle-averaged values to node-averaged values
  // Implement Equation 9 from your pseudocode here
  for (auto it = m.node_begin(); it != m.node_end(); ++it){
    double sumarea=0;
    QVar sumQ = QVar(0, 0, 0);
    for(auto j = m.triangle_begin(*it); j != m.triangle_end(*it); ++j){
      sumarea += (*j).area();
      sumQ += (*j).value().Q * (*j).area();
    }
    (*it).value().Q = sumQ/sumarea;
    (*it).position().z = (*it).value().Q.h-1;
  }
}

// Construct a Color functor and view with the SDLViewer
template <typename Ship>
struct ColorFunctor {
  std::vector<Ship>& obj_vector;
  
  template <typename NODE>
  CS207::Color operator()(const NODE& n) const {
    if (norm(n.position() - Point(0.75,1,0)) < 0.1) {
      return CS207::Color(0.9,0.2,0.8);
    }
    bool check = false;
    if (obj_vector.size() >= 1 && obj_vector[0].checkCoverNode(n.position())) {
      return CS207::Color(0.1,0.9,0.1);
    }
    if (obj_vector.size() >= 2) {
      for (unsigned i = 1; i < obj_vector.size(); ++i) {
        if (obj_vector[i].checkCoverNode(n.position())) {
          check = true;
        }
      }
    }
    if (check) {
      return CS207::Color(0.9,0.1,0.2);
    }
    double h = n.value().Q.h;
    if (h < 1.01) {
      return CS207::Color(0.0,0.0,1.0);
    } else if (h >= 1.2) {
      return CS207::Color(1.0,1.0,1.0);
    } else {
      return CS207::Color((h-1)*4.8, (h-1)*4.8, 1.0);
    }
  };
};

template <typename Pred, typename It, typename G>
class filter_iterator
    : private equality_comparable<filter_iterator<Pred,It,G>> {
 public:
  // Get all of the iterator traits and make them our own
  typedef typename std::iterator_traits<It>::value_type        value_type;
  typedef typename std::iterator_traits<It>::pointer           pointer;
  typedef typename std::iterator_traits<It>::reference         reference;
  typedef typename std::iterator_traits<It>::difference_type   difference_type;
  typedef typename std::input_iterator_tag                     iterator_category;

  typedef filter_iterator<Pred,It,G> self_type;

  // Constructor
  //filter_iterator<Pred,Iter>(p, it, end, mesh, list)
  filter_iterator(const Pred& p, const It& first, const It& last,const G& m, std::vector<unsigned> l)
      : p_(p), it_(first), end_(last),mesh(m),list(l) {
  }

  /** Return a new edge object filter_iterator itself pointing to
   *
   * @pre (*this) valid
   *
   * Complexity: O(1).
   */
  value_type operator*() const {
    return *it_;
  }
  /** Return a new edge object filter_iterator itself pointing to
   */
  self_type& operator++() {
    do {
      ++it_;
    } while (it_ != end_ && !p_(mesh,list,*it_));
    return *this;
  }
  /** Test whether this IncidentIterator and @a st are equal.
   *
   * Equal IncidentIterator have the same graph and the same iterator position and end position.
   */
  bool operator==(const self_type& st) const {
    return ((st.it_ == it_) && (st.end_ == end_));
  }

 private:
  Pred p_;
  It it_;
  It end_;
  G mesh; 
  std::vector<unsigned> list;
};

/** Helper function for constructing filter_iterators.
 *
 * Usage:
 * // Construct an iterator that filters odd values out and keeps even values.
 * std::vector<int> a = ...;
 * auto it = make_filtered(a.begin(), a.end(), [](int k) {return k % 2 == 0;});
 */
template <typename Pred, typename Iter,typename G>
filter_iterator<Pred,Iter,G> make_filtered(const Iter& it, const Iter& end,
                                         const Pred& p,const G& mesh, std::vector<unsigned> list) {
  return filter_iterator<Pred,Iter,G>(p, it, end, mesh, list);
}

//template<typename G,typename NODE>
struct MyPredicate{
  template<typename G,typename NODE>
  bool operator()(G& g1,std::vector<unsigned> list1,const NODE& n) {
    double z0=0;
   for (auto it1 = list1.begin(); it1 != list1.end(); ++it1){
      Node node = g1.node(*it1);
      z0 += node.position().z;
      //if(node.position().z>z0)
      //  z0=node.position().z;
    }
    if (list1.size()==0)
      return 1;
    z0/=list1.size();
    //std::cout<< (n.position().z >z0) << std::endl;
    return ( n.position().z >z0 );
  }
};

template<typename G>
struct PlaneConstraint
{
  PlaneConstraint(double h): height(h) {} 
  void operator()(G& g, double){
    for (auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      Node node = (*it);
      if (dot(node.position(), Point(0, 0, 1)) < height){
        node.position().elem[2] = height;
        node.value().velocity.elem[2] = 0;
      }
    }
  }
  double height;
};

template<typename G>
struct CollisionConstraint
{
  void operator()(G& g1,G& g2, std::vector<unsigned> list1,std::vector<unsigned> list2){
    Point center1 = Point(0, 0, 0);
      for (auto it=g1.node_begin(); it != g1.node_end(); ++it){
        center1 += (*it).position()/g1.num_nodes();
      }
    Point center2 = Point(0, 0, 0);
      for (auto it=g2.node_begin(); it != g2.node_end(); ++it){
        center2 += (*it).position()/g2.num_nodes();
      }
    Point center0 = (center1+center2)/2;
    Point n1 = (center0-center1)/norm(center0-center1);
    Point n2 = (center0-center2)/norm(center0-center2);
 
    for (auto it1 = list1.begin(); it1 != list1.end(); ++it1){
      Node node = g1.node(*it1);
      Point p = node.position();
      Point p1 = center0-p;
      Point v = node.value().velocity;
      Point v1 = center0-v;  
      node.position() = dot(n1,p1)*n1+p;
      node.value().velocity = v -dot(n1,v1)*(-n1);             
    }
    
    for (auto it2 = list2.begin(); it2 != list2.end(); ++it2){
      Node node = g2.node(*it2);
      Point p = node.position();
      Point p1 = center0-p;
      Point v = node.value().velocity;
      Point v1 = center0-v;  
      node.position() = dot(n2,p1)*n2+p;
      node.value().velocity = v -dot(n2,v1)*(-n2);        
    }
  }
  //std::vector<unsigned> list1;
  //std::vector<unsigned> list2;
};

template<typename G>
struct BoxConstraint
{
  BoxConstraint(double h1, double h2, double left,double right ): h_lower(h1),h_upper(h2),l(left),r(right) {} 
  void operator()(G& g, double){
    for (auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      Node node = (*it);
      if (dot(node.position(), Point(0, 0, 1)) < h_lower){
        node.position().elem[2] = h_lower;
        node.value().velocity.elem[2] = 0;
      }
      if (dot(node.position(), Point(0, 0, 1)) > h_upper){
        node.position().elem[2] = h_lower;
        node.value().velocity.elem[2] = 0;
      }
      if (dot(node.position(), Point(0, 1, 0)) < l){
        node.position().elem[1] = l;
        node.value().velocity.elem[1] = 0;
      }
      if (dot(node.position(), Point(0, 1, 0)) > r){
        node.position().elem[1] = r;
        node.value().velocity.elem[1] = 0;
      }
    }
  }
  double h_lower;
  double h_upper;
  double l;
  double r;

};

template<typename G>
struct SphereConstraint
{
  SphereConstraint(Point center, double radius):c(center), r(radius) {}
  void operator()(G& g, double){
    for (auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      Node node = (*it);
      if (norm(node.position()-c) < r){
        Point R = (node.position()-c)/norm(node.position()-c);
        node.position() = c + R*r;
        node.value().velocity = node.value().velocity - dot(node.value().velocity, R)*R;
      }
    }
  }
  Point c;
  double r;
};

template<typename G>
struct ConstantConstraint
{
  ConstantConstraint(Point P1, Point P2):p1(P1), p2(P2) {}
  void operator()(G& g, double){
    for (auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      auto n = *it;
      if (n.position() == p1 || n.position() == p2){
        n.value().velocity = Point(0, 0, 0);
      }
    }
  }
  Point p1;
  Point p2;
};

template<typename C1, typename C2, typename G>
struct CombinedConstraint
{
  C1 cons1;
  C2 cons2;
  CombinedConstraint(C1 c1=C1(), C2 c2=C2()):cons1(c1), cons2(c2){}
  void operator()(G& g, double){
    cons1(g, 0);
    cons2(g, 0);
  }
};

template<typename C1, typename C2, typename G>
CombinedConstraint<C1, C2, G> make_combined_constraint(C1 c1, C2 c2, G& g){
  (void) g;
  return CombinedConstraint<C1, C2, G>(c1, c2);
}



/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on Node
 *           at time @a t.
 */
template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the {n+1} node positions
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    n.position() += n.value().velocity * dt;
  }
  // Compute the {n+1} node velocities
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    n.value().velocity += force(n, t) * (dt / n.value().mass);
  }
  return t + dt;
}

//Force Function to calculate gravity
struct GravityForce
{
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return (Point(0, 0, -grav)*n.value().mass);
  }
};

//Force Function to calculate spring force
struct MassSpringForce
{
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point spring = Point(0, 0, 0);  //spring force
    //add up all spring forces
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
      auto incident = *it;

      if(incident.node1()==n)
        spring += ((incident.value().K)*(incident.node2().position()-incident.node1().position())/incident.length()*(incident.length()-incident.value().L));
      else
        spring += ((incident.value().K)*(incident.node1().position()-incident.node2().position())/incident.length()*(incident.length()-incident.value().L));
    }
    (void) t;
    return spring;
  }
};

//Force Function to calculate damp force
struct DampingForce
{
  DampingForce(double coef): c(coef) {}
  template <typename NODE>
  Point operator()(NODE n, double t){
    (void) t;
    return (-(c*n.value().velocity));
  }
  double c;
};


// The wind force
struct WindForce {
  WindForce(Point wind): w(wind) {}

  template <typename NODE>
  Point operator()(NODE n, double t) {
    double c = 0.00004;
    auto normal = Point(0,0,0);
    for (auto it=n.triangle_begin(); it!=n.triangle_end(); ++it){
      Point tnorm;
      tnorm = cross((*it).node(0).position()-(*it).node(1).position(), 
        (*it).node(0).position()-(*it).node(2).position());
      tnorm = tnorm/norm(tnorm);
      normal = normal + tnorm;
    }
    (void) t;
    return c*dot((w-n.value().velocity),normal)*normal;
  }
  Point w;
};


//the air pressure force
template <typename NODE, typename G>
struct PressureForce
{
  PressureForce(double p_out, double c, G* graph): P_out(p_out), C(c), g(graph) {}

  Point operator()(NODE n, double t) {

    //if n.index()==0, update the volume, center and normal vector,
    //P_diff, for each node
    if(n.index() == 0){
      //update the center
      center = Point(0, 0, 0);
      for (auto it=g->node_begin(); it != g->node_end(); ++it){
        center += (*it).position()/g->num_nodes();
      }

      //update the outward normal vector
      for (auto it=g->triangle_begin(); it != g->triangle_end(); ++it){
        Point tnorm;
        tnorm = cross((*it).node(0).position()-(*it).node(1).position(), 
          (*it).node(0).position()-(*it).node(2).position());
        if (dot(tnorm, center-(*it).node(0).position())>0)
          tnorm = -tnorm;
        tnorm = tnorm/norm(tnorm);
        (*it).value().n = tnorm;
      }
      //update the Volume
      V = 0;
      for (auto it=g->triangle_begin(); it != g->triangle_end(); ++it){
        V += (*it).value().n.z*(*it).area()*((*it).node(0).position().z +
          (*it).node(1).position().z + (*it).node(2).position().z)/3;
      }
      //P_diff
      P_diff = C/V - P_out;
    }

    //for any node, calculate the force
    Point force = Point(0, 0, 0);
    for (auto it=n.triangle_begin(); it != n.triangle_end(); ++it){
      force += P_diff*(*it).area()*(*it).value().n/3;
    }
    (void) t;
    //std::cout<<force<<std::endl;
    return force;
  }
private:
  double P_out; //the pressure outside the ball
  double C;  //nRT
  double V;  //the volumn of the ball
  double P_diff; //P_inside-P_out
  Point center; //The point in the center of the ball
  G* g;
};

//Force function which represents the combined effects of F1 and F2
template<typename F1,typename F2>
struct CombinedForce{
  F1 force1;
  F2 force2;
  CombinedForce(F1 f1=F1(), F2 f2=F2()):force1(f1), force2(f2){}
  template <typename NODE>
  Point operator() (NODE n, double t){
    return (force1(n, t)+force2(n, t));
  }
};


//Combine the effects of two forces
template<typename F1,typename F2>
CombinedForce<F1, F2> make_combined_force(F1 f1 = F1(), F2 f2 = F2()){
  return CombinedForce<F1, F2>(f1, f2);
}
//combine the effects of three forces
template<typename F1, typename F2, typename F3>
CombinedForce<CombinedForce<F1, F2>, F3> make_combined_force(F1 force1, F2 force2, F3 force3){
  return make_combined_force(make_combined_force(force1, force2), force3);
}






int main(int argc, char** argv) {
  // Check arguments
  if (argc < 5) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a mesh
  MeshType mesh;
  MeshType mesh2;

  std::vector<typename MeshType::node_type> mesh_node;
  // Read all Points and add them to the Mesh
  std::ifstream nodes_file(argv[1]);
  Point p;
  while (CS207::getline_parsed(nodes_file, p)) {
    mesh_node.push_back(mesh.add_node(p));
  }

  // Read all mesh triangles and add them to the Mesh
  std::ifstream tris_file(argv[2]);
  std::array<int,3> t;
  while (CS207::getline_parsed(tris_file, t)) {
    mesh.add_triangle(mesh_node[t[0]], mesh_node[t[1]], mesh_node[t[2]]);
  }

  std::vector<typename MeshType::node_type> mesh_node2;
  // Read all Points and add them to the Mesh
  std::ifstream nodes_file2(argv[3]);
  while (CS207::getline_parsed(nodes_file2, p)) {
    mesh_node2.push_back(mesh2.add_node(p));
  }

  // Read all mesh triangles and add them to the Mesh
  std::ifstream tris_file2(argv[4]);
  while (CS207::getline_parsed(tris_file2, t)) {
    mesh2.add_triangle(mesh_node2[t[0]], mesh_node2[t[1]], mesh_node2[t[2]]);
  }

  // Print out the stats
  std::cout << mesh.num_nodes() << " "
            << mesh.num_edges() << " "
            << mesh.num_triangles() << std::endl;
  std::cout << mesh2.num_nodes() << " "
            << mesh2.num_edges() << " "
            << mesh2.num_triangles() << std::endl;

  // set the position of the sphere
  for (auto it = mesh.node_begin(); it != mesh.node_end(); ++it){
    (*it).position() *= 1;
    (*it).position().z += 1.3;
  }

  for (auto it = mesh2.node_begin(); it != mesh2.node_end(); ++it){
    (*it).position() *= 5;
    //(*it).position().z += 1.2;
  }


  //set the mass and velocity of each Node for Sphere
  for (auto it = mesh.node_begin(); it != mesh.node_end(); ++it){
    (*it).value().mass = float(1)/mesh.num_nodes();
    (*it).value().velocity = Point(0, 0, 0);
  }


  //set K and L for each edge for Sphere
  for (auto it = mesh.node_begin(); it != mesh.node_end(); ++it)
  {
    for (auto j = (*it).edge_begin(); j != (*it).edge_end(); ++j){
       (*j).value().L = (*j).length();
       (*j).value().K = 400;
    }
  }

  for (auto it = mesh2.node_begin(); it != mesh2.node_end(); ++it) 
  {
    (*it).value().Q = QVar(1,0,0);
  }

  for (auto it=mesh2.triangle_begin(); it!=mesh2.triangle_end(); ++it) {
    (*it).value().Q = ((*it).node(0).value().Q + (*it).node(1).value().Q + (*it).node(2).value().Q)/3;
  }

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(mesh);
  auto node_map2 = viewer.empty_node_map(mesh2);
  viewer.launch();

  viewer.add_nodes(mesh.node_begin(), mesh.node_end(), node_map);
  viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);
  viewer.add_nodes(mesh2.node_begin(), mesh2.node_end(),node_map2);
  viewer.add_edges(mesh2.edge_begin(), mesh2.edge_end(), node_map2);

  viewer.center_view();

  double min_edge_length = (*mesh2.edge_begin()).length();
  
  for (auto iter = mesh2.edge_begin(); iter != mesh2.edge_end(); ++iter) {
    if (min_edge_length > (*iter).length()) {
      min_edge_length = (*iter).length();
    }
  }

  double max_height = mesh2.node(0).value().Q.h;

  for (auto iter = mesh2.node_begin(); iter != mesh2.node_end(); ++iter) {
    if (max_height < (*iter).value().Q.h) {

      max_height = (*iter).value().Q.h;
    }
  }
  double dt = 0.25 * min_edge_length / (sqrt(grav * max_height));
  dt = 0.0013;
  //dt = 0.0002;

  //Begin the mass-spring simulation
  double t_start = 0.0;
  double t_end   = 10.0;

  // Preconstruct a Flux functor
  EdgeFluxCalculator<Ship> f;
  Ship obj0 = Ship();
    
  std::vector<Ship> obj_vector;
  obj_vector.push_back(obj0);
  
  
  //three color parameter
  int color1 = 1; 
  int color2 = 1; 
  int color3 = 1; 
  
  //Create listener
  Pause_listener* pause = new Pause_listener(dt); 
  Speed_listener* speed = new Speed_listener(dt, dt); 
  XYZ_listener<MeshType>* xyz = new XYZ_listener<MeshType>(&mesh);
  //Color_listener* col = new Color_listener(&color1, &color2, &color3);
  
  //add listener
  viewer.add_listener(pause);
  viewer.add_listener(speed);
  viewer.add_listener(xyz);
  //viewer.add_listener(col);

  //Initialize forces
  WindForce wind_force(Point(0,0,0));
  PressureForce<typename MeshType::node_type, MeshType> pressure_force(0.2, 1, &mesh);
  DampingForce damp_force(float(1)/mesh.num_nodes());
  //auto force = GravityForce();
  //Initialize constriants
  auto constraint = PlaneConstraint<MeshType>(-4);
  //auto constraint = BoxConstraint<MeshType>(-2.2,2.0,-2.0,1.8);
  //auto constraint = make_combined_constraint(,)
  for (double t = t_start; t < t_end; t += dt) {

    constraint(mesh, 0);
    //auto collision_constrain = CollisionConstraint<MeshType>();
    CollisionDetector<MeshType> c;
    c.add_object(mesh);
    c.add_object(mesh2);
    c.check_collisions();
    std::vector<unsigned> collision;
    std::vector<unsigned> collision2;

    for (auto it=c.begin(); it!= c.end(); ++it){
      auto boom = *it;
      Node n = boom.n1;
      if (boom.mesh1 == &mesh)
        collision.push_back(n.index());
      if (boom.mesh1 == &mesh2)
        collision2.push_back(n.index());
    }

    //std::cout<<collision2.size()<<std::endl;

    wind_force.w.z = collision2.size() * 10;
    obj_vector[0].Weight = collision2.size() * 100;
    obj_vector[0].length = collision2.size() * 500;

    auto force = make_combined_force(MassSpringForce(), GravityForce(), make_combined_force(pressure_force, damp_force, wind_force));

    symp_euler_step(mesh, t, dt, force);
    // Step forward in time with forward Euler
    hyperbolic_step2(mesh2, f, t, dt, obj_vector);
    // Update node values with triangle-averaged values
    post_process(mesh2);

    viewer.set_label(t);
    
    //update with removed nodes
    //update viewer with new positions and new edges
    viewer.clear();
    node_map.clear();
    node_map2.clear();
    auto it_begin = make_filtered(mesh.node_begin(), mesh.node_end(),MyPredicate(),mesh2,collision2);
    auto it_end = make_filtered(mesh.node_end(), mesh.node_end(),MyPredicate(),mesh2,collision2);
    viewer.add_nodes(it_begin, it_end, color(color1, color2, color3), node_map);
    viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);
    //viewer.add_nodes(mesh.node_begin(), mesh.node_end(), color(color1, color2, color3), node_map);
    viewer.add_nodes(mesh2.node_begin(), mesh2.node_end(), color(color1, color2, color3), node_map2);
    viewer.add_edges(mesh2.edge_begin(), mesh2.edge_end(), node_map2);
    // These lines slow down the animation for small graphs
    if (mesh.num_nodes() < 100)
      CS207::sleep(0.001);
  }

  return 0;
}
