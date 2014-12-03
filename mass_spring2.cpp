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

#include "Graph.hpp"
#include "Point.hpp"
#include "Mesh.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point velocity;  //< Node velocity
  double mass;     //< Node mass
  double c; //Damping coefficient
};
/**Custom structure of data to store with Edges*/
struct EdgeData{
	double L;
	double K;
};
// Define your Graph type
typedef Mesh<NodeData, EdgeData, bool> MeshType;
typedef typename MeshType::node_type Node;
typedef typename MeshType::edge_type Edge;


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
//template <typename G, typename F>
//double symp_euler_step(G& g, double t, double dt, F force) {

template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
   //Compute the {n+1} node positions
  constraint(g, t);
  	for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    if(n.position() == Point(1,0.5,0) || n.position() == Point(1,-0.5,0)) continue;
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position()=n.position()+(n.value().velocity) * dt;
  }
  // adding constraints to graph
  
  // Compute the {n+1} node velocities
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    if(n.position() == Point(1,0.5,0) || n.position() == Point(1,-0.5,0)) continue;
    n.value().velocity += force(n, t) * (dt / n.value().mass);
  }
  return t + dt;
}


/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force being applied to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
   /*Initial Position*/
  Point operator()(Node n, double t) {
    // HW2 #1: YOUR CODE HERE
    (void) t;
    Point sum;
    double K = 100;
    if(n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
    	return Point(0,0,0);
	for(auto it = n.edge_begin(); it != n.edge_end(); ++it){
	    Edge e = *it;
	    Point xi = e.node1().position();
	    Point xj = e.node2().position();
  	    sum += -K*((xi-xj)/(e.length()))*(e.length()-e.value().L);		
  	}
  	sum += Point(0,0,-grav) * n.value().mass;
    return sum;
  }
};

/** the gravity structure */
struct GravityForce{
    Point operator()(Node n, double t){
        (void) t;
        return Point(0,0,-grav) * n.value().mass;
    }
};
/** the MassSpringForce structure */

struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    auto f = Point(0,0,0);
    for (auto it=n.edge_begin(); it!=n.edge_end(); ++it){
        if ((*it).node1()==n){
          f = f + (-1.0)*(*it).value().K*(n.position()-(*it).node2().position())/(*it).length()*((*it).length()-(*it).value().L);
        }
        if ((*it).node2()==n){
          f = f + (-1.0)*(*it).value().K*(n.position()-(*it).node1().position())/(*it).length()*((*it).length()-(*it).value().L);
        }
    }
    (void) t;
    return f;
  }
};


/** the DampingForce structure */
struct DampingForce{
    Point operator()(Node n, double t){
	(void) t;
	return (-1.0) * n.value().c * n.value().velocity;
    }
};



/**Constraint for a plane z = -0.75
 *a node violates thie constraint if xi*(0,0,1) < -0.75
 *when a node violates this constraint
 *set the position to the nearest point on the plane
 *set the z-component of the Node velocity to zero
 */
struct PlaneConstraint {
  template<typename G>
  void operator()(G& g, double t) {
  (void) t;
    for (auto it=g.node_begin(); it!=g.node_end(); ++it){
      Node n = *it;
      Point p = Point(0,0,1);
      if (dot(n.position(),p) < -0.75){
        n.position()=Point(n.position().x, n.position().y, -0.75);
        n.value().velocity.z = 0;
      }
    }
  }
};

/**Constraint for a sphere centered at c = (0.5,0.5,-0.5) with radus r = 0.15
 * a node violates thie constraint if |xi - c| < r
 * when a node violates this constraint
 * set the position to the nearest point on the surface of the sphere
 * set the component of the velocity that is normal to the sphere's surface to zero
 */
struct SphereConstraint {
  template<typename G>
  void operator()(G& g, double t) {
  	(void) t;
    Point c = Point(0.5,0.5,-0.5);
    double r = 0.15;
    for (auto it=g.node_begin(); it!=g.node_end(); ++it){
      Node n = *it;
      double distance = norm(n.position()-c);
      Point R = (n.position()-c)/distance;
      if (distance < r){
        n.position() = c + ((n.position()-c)/distance)*r;
        n.value().velocity -= dot(n.value().velocity, R)*R;
      }
    }
  }
};
/**Constraint for a sphere centered at c = (0.5,0.5,-0.5) with radus r = 0.15
 * a node violates thie constraint if |xi - c| < r
 * when a node violates this constraint
 * remove the node and all of its edges
 */
/*struct SphereConstraintRemove {
  template <typename G>
    void operator()(G& g, double t) {
  	(void) t;
    Point c = Point(0.5,0.5,-0.5);
    double r = 0.15;
    for (auto it=g.node_begin(); it!=g.node_end(); ++it){
      Node n = *it;
      double distance = norm(n.position()-c);
      if (distance < r)   g.remove_node(n);
    }
  }
};*/

// Combine two templated functors and return the sum of the two forces in the functor
template <typename F1, typename F2>
struct force2 {
  F1 f1_;
  F2 f2_;
  force2(const F1& f1, const F2& f2): f1_(f1), f2_(f2) {
  }
  Point operator()(Node n, double t) {
    return f1_(n,t) + f2_(n,t);
  }
};  

// Combine two templated functors and return the sum of the three forces in the functor
template <typename F1, typename F2>
force2<F1, F2> make_combined_force(const F1& f1, const F2& f2){
  return force2<F1, F2>(f1, f2);
}

// Combine three templated functors and return the sum of the three forces in the functor
template <typename F1, typename F2, typename F3>
struct force3 {
  F1 f1_;
  F2 f2_;
  F3 f3_;
  force3(const F1& f1, const F2& f2, const F3& f3): f1_(f1), f2_(f2), f3_(f3) {
  }

  Point operator()(Node n, double t) {
    return f1_(n,t) + f2_(n,t) + f3_(n,t);
  }
};  


template <typename F1, typename F2, typename F3>
force3<F1, F2, F3> make_combined_force(const F1& f1, const F2& f2, const F3& f3){
  return force3<F1, F2, F3>(f1, f2, f3);
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

  // Set initial conditions for your nodes, if necessary.
 for(auto it = mesh.node_begin(); it!= mesh.node_end(); ++it){
 	(*it).value().mass = 1.0/mesh.num_nodes();
 	(*it).value().velocity = Point(0,0,0);
 	(*it).value().c = 1.0/mesh.num_nodes();
 }
 for(auto it = mesh.edge_begin(); it!= mesh.edge_end(); ++it){
 	(*it).value().L = (*it).length();
 	(*it).value().K = 100;
 }
  // Construct Forces/Constraints

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
    //symp_euler_step(mesh, t, dt, MassSpringForce());
    symp_euler_step(mesh, t, dt, make_combined_force(GravityForce(), MassSpringForce(),DampingForce()),PlaneConstraint());
    viewer.clear();
    node_map.clear();
    // Update viewer with nodes' new positions
   viewer.add_nodes(mesh.node_begin(), mesh.node_end(), node_map);
    viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);
viewer.set_label(t);

    // These lines slow down the animation for small graphs, like grid0_*.
    // Feel free to remove them or tweak the constants.
    //if (mesh.num_nodes() < 100)
      CS207::sleep(0.001);
  }

  return 0;
}
