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

/** Custom structure of data to store with Triangles */
struct TriData
{
  Point n; //the outward surface normal vector
};

typedef Mesh<NodeData, EdgeData, TriData> MeshType;
typedef typename MeshType::node_type Node;
typedef typename MeshType::edge_type Edge;

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
       (*j).value().K = 100;
    }
  }

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(mesh);
  viewer.launch();

  viewer.add_nodes(mesh.node_begin(), mesh.node_end(), node_map);
  viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);

  viewer.center_view();
  



  //Begin the mass-spring simulation
  double dt = 0.0002;
  double t_start = 0.0;
  double t_end   = 10.0;
  
  
  //three color parameter
  int color1 = 1; 
  int color2 = 1; 
  int color3 = 1; 
  
  //Create listener
  Pause_listener* pause = new Pause_listener(dt); 
  Speed_listener* speed = new Speed_listener(dt, dt); 
  XYZ_listener<MeshType>* xyz = new XYZ_listener<MeshType>(&mesh);
  Color_listener* col = new Color_listener(&color1, &color2, &color3);
  
  //add listener
  viewer.add_listener(pause);
  viewer.add_listener(speed);
  viewer.add_listener(xyz);
  viewer.add_listener(col);

  //Initialize forces
  WindForce wind_force(Point(-60,-60,0));
  PressureForce<typename MeshType::node_type, MeshType> pressure_force(1, 500, &mesh);
  DampingForce damp_force(float(1)/mesh.num_nodes());
  //auto force = make_combined_force(MassSpringForce(), GravityForce(), make_combined_force(pressure_force, damp_force, wind_force));
  auto force = make_combined_force(MassSpringForce(), GravityForce(), make_combined_force(damp_force, wind_force));
  //Initialize constriants
  //auto constraint = PlaneConstraint<MeshType>(-2);
  //auto constraint = BoxConstraint<MeshType>(-2.2,2.0,-2.0,1.8);
  
  //auto constraint = ConstantConstraint<MeshType>(Point(-1,1,0),Point(1,1,0));
  auto constraint = make_combined_constraint(PlaneConstraint<MeshType>(-1),ConstantConstraint<MeshType>(Point(-1,1,0),Point(1,1,0)), mesh);
  for (double t = t_start; t < t_end; t += dt) {

    constraint(mesh, 0);
    symp_euler_step(mesh, t, dt, force);

    viewer.set_label(t);
    
    //update with removed nodes
    //clear teh viewer's node
    viewer.clear();
    node_map.clear();
    //update viewer with new positions and new edges
    viewer.add_nodes(mesh.node_begin(), mesh.node_end(), color(color1, color2, color3), node_map);
    viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);

    // These lines slow down the animation for small graphs
    if (mesh.num_nodes() < 100)
      CS207::sleep(0.001);
  }

  return 0;
}
