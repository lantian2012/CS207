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
   * \param[in] h2 The height of the plane at the top 
   * \param[in] left The posistion of the plane on the left
   * \param[in] right The posistion of the plane on the right
   */
  BoxConstraint(double h1, double h2, double left,double right ): h_lower(h1),h_upper(h2),l(left),r(right) {} 
  
  /** \brief Constrain nodes with box constraints.
    *  
    * \param[in, out] g An object that satisfies the graph concept.
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
   * \param[in] P2 The second point to be pinned
   */
  ConstantConstraint(Point P1, Point P2):p1(P1), p2(P2) {}

  /** \brief Pins two nodes
    *  
    * \param[in, out] g An object that satisfies the graph concept.
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
 * This struct implements a constraint that combines the  
 * effects of two constraints.
 * 
 * \tparam G A class that satisfies the graph concept. 
 * \tparam C1 A class that satisfies the constraint concept.
 * \tparam C2 A class that satisfies the constraint concept.
 */
template<typename C1, typename C2, typename G>
struct CombinedConstraint
{
  C1 cons1; ///< The first constraint
  C2 cons2; ///< The second constraint

  /** \brief Constructor Function
   *  
   * \param[in] c1 The first constraint
   * \param[in] c2 The first constraint
   */
  CombinedConstraint(C1 c1=C1(), C2 c2=C2()):cons1(c1), cons2(c2){}
  /** \brief Constrain the node with the effect of two constraints
    *  
    * \param[in, out] g An object that satisfies the graph concept.
    * \param[in] t Time.(Not used here)
    *
    * This function sets the velocity and position nodes, so nodes 
    * will act according to the combined effect of two constraints.
    */
  void operator()(G& g, double){
    cons1(g, 0);
    cons2(g, 0);
  }
};

/** \brief Make a combined constraint from two constraints
  *  
  * \param[in] g An object that satisfies the graph concept.
  * \param[in] c1 The first constraint
  * \param[in] c2 The second constraint
  *
  * \return A constraint that acts as the combined effect of two constraints.
  */
template<typename C1, typename C2, typename G>
CombinedConstraint<C1, C2, G> make_combined_constraint(C1 c1, C2 c2, G& g){
  (void) g;
  return CombinedConstraint<C1, C2, G>(c1, c2);
}



/** \brief Change nodes according to the symplectic Euler method
 * \param[in,out] g      Graph
 * \param[in]     t      The current time (useful for time-dependent forces)
 * \param[in]     dt     The time step
 * \param[in]     force  Function object defining the force per node
 * \return the next time step (usually @a t + @a dt)
 *
 * \tparam G A class that satisfies the graph concept. 
 * \tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on Node
 *           at time @a t.
 *
 * Change a graph's nodes according to a step of the symplectic Euler
 * method with the given node force.
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


/** 
 * \struct GravityForce
 *
 * \brief Force Function to calculate gravity
 */
struct GravityForce
{
  /** \brief Calculate the gravity of a node 
    *  
    * \param[in] n The node to calculate the gravity
    * \param[in] t Time.(Not used here)
    *
    * \tparam NODE a class that satisfies the Node concept
    *
    * \return The gravity of the node
    */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return (Point(0, 0, -grav)*n.value().mass);
  }
};

/** 
 * \struct MassSpringForce
 *
 * \brief Force Function to the spring force on a node
 */
struct MassSpringForce
{
  /** \brief Calculate the spring force on a node
    *  
    * \param[in] n The node to calculate the gravity
    * \param[in] t Time.(Not used here)
    *
    * \tparam NODE a class that satisfies the Node concept
    *
    * \return the spring force of the node
    *
    * This function calculates the spring force of a node
    * The force is calcuated by Hook's law. The initial length
    * of the spring can be set by the user.
    */
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


/** 
 * \struct DampingForce
 *
 * \brief Force Function to calculate damp force
 */
struct DampingForce
{
  /** \brief Constructor Function
   *  
   * \param[in] coef The damping coefficient
   */
  DampingForce(double coef): c(coef) {}
  /** \brief Calculate the damping force on a node
    *  
    * \param[in] n The node to calculate the gravity
    * \param[in] t Time.(Not used here)
    *
    * \tparam NODE a class that satisfies the Node concept
    *
    * \return the damping force of the node
    */
  template <typename NODE>
  Point operator()(NODE n, double t){
    (void) t;
    return (-(c*n.value().velocity));
  }
  double c; ///< The damping coefficient
};


/** 
 * \struct WindForce
 *
 * \brief Force Function to calculate wind force
 */
struct WindForce {
  /** \brief Constructor Function
   *  
   * \param[in] wind The direction and strength of the wind force
   */
  WindForce(Point wind): w(wind) {}

  /** \brief Calculate the wind force on a node
    *  
    * \param[in] n The node to calculate the gravity
    * \param[in] t Time.(Not used here)
    *
    * \tparam NODE a class that satisfies the Node concept
    *
    * \return the wind force of the node
    * 
    * The wind force of a  node is calculated to simulate the aggregated
    * effect of the wind force on adjacent triangles of the node.
    */
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

  Point w; ///< The strength and direction of the wind
};


/** 
 * \struct PressureForce
 *
 * \brief Force Function to the air pressure force 
 */
template <typename NODE, typename G>
struct PressureForce
{
  /** \brief Constructor Function
   *  
   * \param[in] p_out The airpressure outside the sphere
   * \param[in] c The amout of air inside the sphere
   * \param[in] graph Pointer to the sphere mesh to calculate air pressure on
   */
  PressureForce(double p_out, double c, G* graph): P_out(p_out), C(c), g(graph) {}

  /** \brief Calculate the pressure force on a node
    *  
    * \param[in] n The node to calculate the gravity
    * \param[in] t Time.(Not used here)
    *
    * \tparam NODE a class that satisfies the Node concept
    *
    * \return the pressure force of the node
    * 
    * The pressure force of a node is calculated by aggreggating the 
    * pressure force of the adjacent triangles of the node. The force
    * on one triangle is: F = (c/V - p_out)*area. V is the volume of 
    * the sphere, and area is the area of the triangle.
    */
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
    return force;
  }
private:
  double P_out; ///< the pressure outside the ball
  double C;  ///< nRT
  double V;  ///< the volume of the ball
  double P_diff; ///< P_inside-P_out
  Point center; ///< The point in the center of the ball
  G* g; ///< Pointer to the sphere
};

//Force function which represents the combined effects of F1 and F2
/** 
 * \struct CombinedForce
 *
 * \brief Constraint that combines two forces
 * 
 * This struct implements a force that combines the  
 * effects of two forces.
 * 
 * \tparam F1 A class that satisfies the force concept.
 * \tparam F2 A class that satisfies the force concept.
 */
template<typename F1,typename F2>
struct CombinedForce{
  F1 force1; ///< The first force
  F2 force2; ///< The second force

  /** \brief Constructor Function
   *  
   * \param[in] c1 The first force
   * \param[in] c2 The first force
   */
  CombinedForce(F1 f1=F1(), F2 f2=F2()):force1(f1), force2(f2){}

  /** \brief Get the force on the node with the effect of two force
    *  
    * \param[in] n The node to calculate the gravity
    * \param[in] t Time.(Not used here)
    *
    * \tparam NODE a class that satisfies the Node concept
    *
    * \return the combined force of the node
    *
    * This function calculates the force on a node. The force is 
    * the same calculating two forces separately and adding them together
    */
  template <typename NODE>
  Point operator() (NODE n, double t){
    return (force1(n, t)+force2(n, t));
  }
};


/** \brief Make a combined force from two forces
  *  
  * \param[in] f1 The first force
  * \param[in] f2 The second force
  *
  * \return A force that acts as the combined effect of two forces.
  */
template<typename F1,typename F2>
CombinedForce<F1, F2> make_combined_force(F1 f1 = F1(), F2 f2 = F2()){
  return CombinedForce<F1, F2>(f1, f2);
}

/** \brief Make a combined force from three forces
  *  
  * \param[in] f1 The first force
  * \param[in] f2 The second force
  * \param[in] f3 The third force
  *
  * \return A force that acts as the combined effect of three forces.
  */
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
