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

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point velocity;  //< Node velocity
  double mass;     //< Node mass
};

/** Custom structure of data to store with Edges */
struct EdgeData{
  double L;
  double K;
};


typedef Mesh<NodeData, EdgeData, bool> MeshType;
typedef typename MeshType::node_type Node;
typedef typename MeshType::edge_type Edge;

template<typename G>
struct PlaneConstraint
{
  void operator()(G& g, double){
    for (auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      Node node = (*it);
      if (dot(node.position(), Point(0, 0, 1)) < -0.75){
        node.position().elem[2] = -0.75;
        node.value().velocity.elem[2] = 0;
      }
    }
  }
};

template<typename G>
struct SphereConstraint
{
  Point c = Point(0.5, 0.5, -0.5);
  double r = 0.15;
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
};

template<typename G>
struct ConstantConstraint
{
  void operator()(G& g, double){
    for (auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      auto n = *it;
      if (n.position() == Point(-1, 1, 0) || n.position() == Point(1, 1, 0)){
        n.value().velocity = Point(0, 0, 0);
      }
    }
  }
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
  auto constraint = make_combined_constraint(make_combined_constraint(ConstantConstraint<G>(), SphereConstraint<G>(), g), PlaneConstraint<G>(), g);
  //auto constraint = ConstantConstraint<G>();
  constraint(g, 0);
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


//Force Function
struct Problem1Force {
  /** Return the force being applied to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  Point operator()(Node n, double t) {
    //constrain the corners
    //if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
    //  return Point(0, 0, 0);
    Point spring = Point(0, 0, 0);  //spring force
    Point gravity;   //gravity force
    gravity = Point(0, 0, -grav)*n.value().mass;
    //add up all spring forces
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
      Edge incident = *it;
      spring += ((incident.value().K)*(incident.node2().position()-incident.node1().position())/incident.length()*(incident.length()-incident.value().L));
    }
    (void) t;
    return (spring+gravity);
  }
};

//Force Function to calculate gravity
struct GravityForce
{
  Point operator()(Node n, double t) {
    (void) t;
    return (Point(0, 0, -grav)*n.value().mass);
  }
};

//Force Function to calculate spring force
struct MassSpringForce
{
  Point operator()(Node n, double t) {
    Point spring = Point(0, 0, 0);  //spring force
    //add up all spring forces
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
      Edge incident = *it;

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
  Point operator()(Node n, double ){
    return (-(c*n.value().velocity));
  }
  static double c;
};
double DampingForce::c = 0;

//Force function which represents the combined effects of F1 and F2
template<typename F1,typename F2>
struct CombinedForce{
  F1 force1;
  F2 force2;
  CombinedForce(F1 f1=F1(), F2 f2=F2()):force1(f1), force2(f2){}
  Point operator() (Node n, double t){
    (void) t;
    return (force1(n, 0)+force2(n, 0));
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

 #if 1
  std::vector<typename MeshType::node_type> mesh_node;
 #endif

  // Read all Points and add them to the Mesh
  std::ifstream nodes_file(argv[1]);
  Point p;
  while (CS207::getline_parsed(nodes_file, p)) {
    // HW4B: Need to implement add_node before this can be used!
#if 1
    mesh_node.push_back(mesh.add_node(p));
#endif
  }

  // Read all mesh triangles and add them to the Mesh
  std::ifstream tris_file(argv[2]);
  std::array<int,3> t;
  while (CS207::getline_parsed(tris_file, t)) {
    // HW4B: Need to implement add_triangle before this can be used!
#if 1
    mesh.add_triangle(mesh_node[t[0]], mesh_node[t[1]], mesh_node[t[2]]);
#endif
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
       (*j).value().K = 100;
    }
  }

  //set the damping force constriant
  DampingForce::c = float(1)/mesh.num_nodes();
  
  // Print out the stats
  std::cout << mesh.num_nodes() << " " << mesh.num_edges() << std::endl;

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(mesh);
  viewer.launch();

  viewer.add_nodes(mesh.node_begin(), mesh.node_end(), node_map);
  viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);

  viewer.center_view();

  // Begin the mass-spring simulation
  double dt = 0.001;
  double t_start = 0.0;
  double t_end   = 5.0;

  for (double t = t_start; t < t_end; t += dt) {
    //std::cout << "t = " << t << std::endl;
    symp_euler_step(mesh, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce()));

    // Update viewer with nodes' new positions
    //viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
    viewer.set_label(t);
    
    //update with removed nodes
    //clear teh viewer's node
    viewer.clear();
    node_map.clear();
    //update viewer with new positions and new edges
    viewer.add_nodes(mesh.node_begin(), mesh.node_end(), node_map);
    viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);

    // These lines slow down the animation for small graphs, like grid0_*.
    // Feel free to remove them or tweak the constants.
    if (mesh.num_nodes() < 100)
      CS207::sleep(0.001);
  }

  return 0;
}
