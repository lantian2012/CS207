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

#include "Mesh.hpp"
#include "Point.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** 
 * \struct NodeData
 *
 * \brief Custom structure of data to store with Nodes 
 * 
 * This struct stores data associated with each node.
 * Users can store their own data as members, and add
 * new member functions. However, users should NOT delete
 * the @a velocity and @a mass defined here. 
 */
struct NodeData {
  Point velocity;  ///< Node velocity
  double mass;     ///< Node mass
};

/** 
 * \struct EdgeData
 *
 * \brief Custom structure of data to store with Edges
 * 
 * This struct stores data associated with each edge.
 * Users can store their own data as members, and add
 * new member functions. However, users should NOT delete
 * the @a L and @a K defined here. 
 */
struct EdgeData{
  double L;  ///< Edge length
  double K;  ///< Edge spring constant (stiffness)
};

/** 
 * \struct TriData
 *
 * \brief Custom structure of data to store with Triangles
 * 
 * This struct stores data associated with each triangle.
 * Users can store their own data as members, and add
 * new member functions. However, users should NOT delete
 * the @a n defined here. 
 */
struct TriData
{
  Point n; ///<the outward surface normal vector of the triangle
};

typedef Mesh<NodeData, EdgeData, TriData> MeshType;
typedef typename MeshType::node_type Node;
typedef typename MeshType::edge_type Edge;

/** 
 * \struct PlaneConstraint
 *
 * \brief Constraint that acts as a horizontal plane
 * 
 * This struct implements a constraint that acts like
 * a plane. This struct can define a plane at any height.
 * 
 * \tparam G A class that satisfies the graph concept. 
 */
template<typename G>
struct PlaneConstraint
{
  /** \brief Constructor Function
   *  
   * \param[in] h the height of the plane
   */
  PlaneConstraint(double h): height(h) {}

  /** \brief Constrain nodes with plane constraints.
    *  
    * \param[in, out] g An object that satisfies the graph concept.
    * 
    * \param[in] t Time.(Not used here)
    *
    * This function sets the velocity and position of nodes, so 
    * nodes will act as they reach a plane and bounce back.
    */
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

  double height; ///< The height of the plane.
};


/** 
 * \struct BoxConstraint
 *
 * \brief Constraint that acts as a box
 * 
 * This struct implements a constraint that acts like
 * a box. This struct can define four planes at any posistion.
 * 
 * \tparam G A class that satisfies the graph concept. 
 */
template<typename G>
struct BoxConstraint
{
  /** \brief Constructor Function
   *  
   * \param[in] h1 The height of the plane at the bottom
   *
   * \param[in] h2 The height of the plane at the top 
   *
   * \param[in] left The posistion of the plane on the left
   *
   * \param[in] right The posistion of the plane on the right
   */
  BoxConstraint(double h1, double h2, double left,double right ): h_lower(h1),h_upper(h2),l(left),r(right) {} 
  
  /** \brief Constrain nodes with box constraints.
    *  
    * \param[in, out] g An object that satisfies the graph concept.
    * 
    * \param[in] t Time.(Not used here)
    *
    * This function sets the velocity and position of nodes, so 
    * nodes will act as they reach a box and bounce back.
    */
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
  double h_lower; ///< The height of the plane at the bottom
  double h_upper; ///< The height of the plane at the top
  double l; ///< The posistion of the plane on the left
  double r; ///< The posistion of the plane on the right

};

/** 
 * \struct ConstantConstraint
 *
 * \brief Constraint that pins two points of a graph
 * 
 * This struct implements a constraint that pins two nodes of 
 * a graph. This struct can define two points at any posistion.
 * 
 * \tparam G A class that satisfies the graph concept. 
 */
template<typename G>
struct ConstantConstraint
{
  /** \brief Constructor Function
   *  
   * \param[in] P1 The first point to be pinned
   *
   * \param[in] P2 The second point to be pinned
   */
  ConstantConstraint(Point P1, Point P2):p1(P1), p2(P2) {}

  /** \brief Pins two nodes
    *  
    * \param[in, out] g An object that satisfies the graph concept.
    * 
    * \param[in] t Time.(Not used here)
    *
    * This function sets the velocity and position of two nodes, so 
    * the two nodes will not move.
    */
  void operator()(G& g, double){
    for (auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      auto n = *it;
      if (n.position() == p1 || n.position() == p2){
        n.value().velocity = Point(0, 0, 0);
      }
    }
  }
  Point p1; ///< The first point to be pinned
  Point p2; ///< The second point to be pinned
};


/** 
 * \struct CombinedConstraint
 *
 * \brief Constraint that combines two constraints
 * 
 * This struct implements a constraint that pins two nodes of 
 * a graph. This struct can define two points at any posistion.
 * 
 * \tparam G A class that satisfies the graph concept. 
 */
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
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a mesh
  MeshType mesh;

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

  // Print out the stats
  std::cout << mesh.num_nodes() << " "
            << mesh.num_edges() << " "
            << mesh.num_triangles() << std::endl;


  //set the mass and velocity of each Node
  for (auto it = mesh.node_begin(); it != mesh.node_end(); ++it){
    (*it).value().mass = float(1)/mesh.num_nodes();
    (*it).value().velocity = Point(0, 0, 0);
  }


  //set K and L for each edge
  for (auto it = mesh.node_begin(); it != mesh.node_end(); ++it)
  {
    for (auto j = (*it).edge_begin(); j != (*it).edge_end(); ++j){
       (*j).value().L = (*j).length();
       (*j).value().K = 8000;
    }
  }

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(mesh);
  viewer.launch();

  // Add nodes and edges to the viewer
  viewer.add_nodes(mesh.node_begin(), mesh.node_end(), node_map);
  viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);
  viewer.center_view();

  //Begin the mass-spring simulation
  double dt = 0.0002;
  double t_start = 0.0;
  double t_end   = 10.0;

  //Initialize forces

  WindForce wind_force(Point(10,80,60));
  PressureForce<typename MeshType::node_type, MeshType> pressure_force(1, 600, &mesh);
  DampingForce damp_force(float(1)/mesh.num_nodes());
  auto force = make_combined_force(MassSpringForce(), GravityForce(), make_combined_force(pressure_force, damp_force, wind_force));
  //Initialize constriants
  //auto constraint = PlaneConstraint<MeshType>(-2.5);
  auto constraint = BoxConstraint<MeshType>(-2.2,2.0,-2.0,1.8);
  //auto constraint = make_combined_constraint(,)
  
  for (double t = t_start; t < t_end; t += dt) {

    // Constrain the nodes' velocity and position
    constraint(mesh, 0);

    // Update the position and velocity of nodes 
    symp_euler_step(mesh, t, dt, force);

    viewer.set_label(t);
    
    //update with removed nodes
    //clear teh viewer's node
    viewer.clear();
    node_map.clear();
    //update viewer with new positions and new edges
    viewer.add_nodes(mesh.node_begin(), mesh.node_end(), node_map);
    viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);

    // These lines slow down the animation for small graphs
    if (mesh.num_nodes() < 100)
      CS207::sleep(0.001);
  }

  return 0;
}
