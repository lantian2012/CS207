/**
 * @file shallow_water.cpp
 * Implementation of a shallow water system using Mesh
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D point list (one per line) defined by three doubles
 * Second file: Triangles (one per line) defined by 3 indices into the point list
 */


#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"
#include "Point.hpp"
#include <fstream>
#include <cmath>

#include "Mesh.hpp"


// Standard gravity (average gravity at Earth's surface) in meters/sec^2
static constexpr double grav = 9.80665;
typedef unsigned size_type;


// HW4B: Placeholder for Mesh Type!
// Define NodeData, EdgeData, TriData, etc
// or redefine for your particular Mesh

struct NodeData
{
  QVar Q;
  NodeData():Q(0.0,0.0,0.0){}
};
struct EdgeData
{
  size_type triangle1;
  size_type triangle2;
  EdgeData():triangle1(-1), triangle2(-1){}
};

/** @struct Mesh::TriData
   * information associated with triangles
   */    
struct TriData{
  std::vector<size_type> nodes; //a vector storing the uids of the three nodes of this triangle
  std::vector<size_type> edges; //a vector storing the uids of the three edges of this triangle
  QVar Q;  //Qk: the average value of Q inside this triangle
  double area;  //the area of this triangle
  std::vector<QVar> F; //The transition over an edge
  std::vector<Point> n; //The unit normal vector of 3 edges
  /**construct an invalid TriData*/
  TriData(): nodes(3,0),edges(3,0),Q(QVar()),area(-1), F(3,QVar(0,0,0)){
  }
};
typedef Mesh<NodeData,EdgeData,TriData> MeshType;


/** Function object for calculating shallow-water flux.
 *          |n
 *   T_k    |---> n = (nx,ny)   T_m
 *   QBar_k |                   QBar_m
 *          |
 * @param[in] nx, ny Defines the 2D outward normal vector n = (@a nx, @a ny)
 *            from triangle T_k to triangle T_m. The length of n is equal to the
 *            the length of the edge, |n| = |e|.
 * @param[in] dt The time step taken by the simulation. Used to compute the
 *               Lax-Wendroff dissipation term.
 * @param[in] qk The values of the conserved variables on the left of the edge.
 * @param[in] qm The values of the conserved variables on the right of the edge.
 * @return The flux of the conserved values across the edge e
 */
struct EdgeFluxCalculator {
  QVar operator()(double nx, double ny, double dt,
                  const QVar& qk, const QVar& qm) {
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
    // HW4B: You may change this to plot something other than the
    // positions of the nodes
    return Point(n.position().x, n.position().y, n.value().Q.h);
    //return n.position();
  }
};


/** Integrate a hyperbolic conservation law defined over the mesh m
 * with flux functor f by dt in time.
 */
template <typename MESH, typename FLUX>
double hyperbolic_step(MESH& m, FLUX& f, double t, double dt) {
  // HW4B: YOUR CODE HERE
  // Step the finite volume model in time by dt.

  // Pseudocode:
  // Compute all fluxes. (before updating any triangle Q_bars)
  // For each triangle, update Q_bar using the fluxes as in Equation 8.
  //  NOTE: Much like symp_euler_step, this may require TWO for-loops
  
  // Implement Equation 7 from your pseudocode here.
  for (auto it=m.triangle_begin(); it!=m.triangle_end(); ++it) {
    QVar temp_sum = QVar(0.0,0.0,0.0);
    int j=0;
    for (auto init = (*it).triangle_begin(); init!=(*it).triangle_end(); ++init) {  
        if((*init).index() !=unsigned(-1)){                         
          QVar temp = f((*it).normal(j).x, (*it).normal(j).y, dt, (*it).Q(), (*init).Q());
          temp *= dt;
          temp /= (*it).area();
          temp_sum += temp;
        }
        else{
          QVar temp = f((*it).normal(j).x, (*it).normal(j).y, dt, (*it).Q(), QVar((*it).Q().h,0.0,0.0));
          temp *= dt;
          temp /= (*it).area();
          temp_sum += temp;
        }
        ++j;
    }
    (*it).Q() -= temp_sum;
  }
  return t + dt;
  

  
  /*
  for(auto i = m.edge_begin(); i != m.edge_end(); ++i){
    if ((*i).value().triangle1 != (unsigned) -1 && (*i).value().triangle2 != (unsigned) -1 ){
      MeshType::Triangle trik = m.triangle((*i).value().triangle1);
      MeshType::Triangle trim = m.triangle((*i).value().triangle2);
      unsigned int edge_k = 0;
      unsigned int edge_m = 0;
      //which edge (*i) is in trik and trim
      while(trik.node(edge_k).index()== (*i).node1().index() 
        || trik.node(edge_k).index()== (*i).node2().index() )
        ++edge_k;
      while(trim.node(edge_m).index()== (*i).node1().index() 
        || trim.node(edge_m).index()== (*i).node2().index() )
        ++edge_m;
      QVar flux = f(trik.normal(edge_k).x, trik.normal(edge_k).y, dt, trik.Q(), trim.Q());
      trik.F(edge_k) = flux;
      trim.F(edge_m) = -flux;
    }
    else{
      MeshType::Triangle trik;
      if ((*i).value().triangle1 != (unsigned) -1)
        trik = m.triangle((*i).value().triangle1);
      else
        trik = m.triangle((*i).value().triangle2);
      unsigned int edge_k = 0;
      while(trik.node(edge_k).index()== (*i).node1().index() 
        || trik.node(edge_k).index()== (*i).node2().index() )
        ++edge_k;
      QVar flux = f(trik.normal(edge_k).x, trik.normal(edge_k).y, dt, trik.Q(), QVar(trik.Q().h, 0, 0));
      trik.F(edge_k) = flux;
    }
  }

  for(auto i = m.triangle_begin(); i != m.triangle_end(); ++i){
    QVar sum = QVar(0, 0, 0);
    for (int j = 0; j < 3; ++j){
      sum += (*i).F(j);
    }
    //std::cout<<"tbefore=="<<(*i).Q().h<<" "<<(*i).Q().hx<<" "<<(*i).Q().hy<<" "<<std::endl;
    (*i).Q() = (*i).Q()-dt/(*i).area()*sum;
    //std::cout<<"tafter=="<<(*i).Q().h<<" "<<(*i).Q().hx<<" "<<(*i).Q().hy<<" "<<std::endl;
  }
  
  return t + dt;*/
};

/** Convert the triangle-averaged values to node-averaged values for viewing. */
template <typename MESH>
void post_process(MESH& m) {
  // HW4B: Post-processing step
  // Translate the triangle-averaged values to node-averaged values
  // Implement Equation 9 from your pseudocode here

  
  for (auto it = m.node_begin(); it != m.node_end(); ++it){
    double sumarea=0;
    QVar sumQ = QVar(0, 0, 0);
    for(auto j = m.triangle_begin(*it); j != m.triangle_end(*it); ++j){
      sumarea += (*j).area();
      sumQ += (*j).Q() * (*j).area();
    }
    (*it).value().Q = sumQ/sumarea;
  }
};


int main(int argc, char* argv[])
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: shallow_water NODES_FILE TRIS_FILE\n";
    exit(1);
  }

  MeshType mesh;
  // HW4B: Need node_type before this can be used!
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

  // HW4B Initialization
  // Set the initial conditions
  // Perform any needed precomputation

  // Case 1: pebble being thrown in a pond
  /*for (auto it=mesh.node_begin(); it!=mesh.node_end(); ++it) {
    (*it).value().Q = QVar(0.0,0.0,0.0);
    (*it).value().Q.h = 1.0-0.75*exp(-80.0*(pow(((*it).position().x-0.75),2)+pow((*it).position().y,2)));
  }*/
  
  
  
  // Case 2: a large column of water
  for (auto it=mesh.node_begin(); it!=mesh.node_end(); ++it) {
    (*it).value().Q = QVar(0.0,0.0,0.0);
    double temp = ((pow(((*it).position().x-0.75),2)+pow((*it).position().y,2))-pow(0.15,2));
    if (temp < 0) {
      (*it).value().Q.h = 1.0+0.75;
    }
    else {
      (*it).value().Q.h = 1.0;
    }
  }


/*
  // Case 3: a dam break
  for (auto it=mesh.node_begin(); it!=mesh.node_end(); ++it) {
    (*it).value().Q = QVar(0.0,0.0,0.0);
    if ((*it).position().x < 0) {
      (*it).value().Q.h = 1.0+0.75;
    }
    else {
      (*it).value().Q.h = 1.0;
    }
  }*/



  // Set triangle values
  for (auto it=mesh.triangle_begin(); it!=mesh.triangle_end(); ++it) {
    (*it).Q() = QVar(0.0,0.0,0.0);
    (*it).Q() += (*it).node(0).value().Q;
    (*it).Q() += (*it).node(1).value().Q;
    (*it).Q() += (*it).node(2).value().Q;
    (*it).Q() /= 3.0;
  }

 

  
  
  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  viewer.launch();

  // HW4B: Need to define Mesh::node_type and node/edge iterator
  // before these can be used!
#if 1
  auto node_map = viewer.empty_node_map(mesh);
  viewer.add_nodes(mesh.node_begin(), mesh.node_end(),
                   CS207::DefaultColor(), NodePosition(), node_map);
  viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);
#endif
  viewer.center_view();

  // HW4B: Timestep
  // CFL stability condition requires dt <= dx / max|velocity|
  // For the shallow water equations with u = v = 0 initial conditions
  //   we can compute the minimum edge length and maximum original water height
  //   to set the time-step
  // Compute the minimum edge length and maximum water height for computing dt
  double min_edge_length =( *mesh.edge_begin()).length();
  for (auto it=mesh.edge_begin(); it!=mesh.edge_end(); ++it) {
    if ((*it).length() < min_edge_length) {
      min_edge_length = (*it).length();
    }
  }
  double max_height = 0.0;
  for (auto it=mesh.node_begin(); it!=mesh.node_end(); ++it) {
    if ((*it).value().Q.h > max_height) {
      max_height = (*it).value().Q.h;
    }
  }
  
#if 1
  double dt = 0.25 * min_edge_length / (sqrt(grav * max_height));
#else
  // Placeholder!! Delete me when min_edge_length and max_height can be computed!
  double dt = 0.1;
#endif
  double t_start = 0;
  double t_end = 10;

  // Preconstruct a Flux functor
  EdgeFluxCalculator f;

  // Begin the time stepping
  for (double t = t_start; t < t_end; t += dt) {
    // Step forward in time with forward Euler
    hyperbolic_step(mesh, f, t, dt);

    // Update node values with triangle-averaged values
    post_process(mesh);

    // Update the viewer with new node positions
    // HW4B: Need to define node_iterators before these can be used!
#if 1
    viewer.add_nodes(mesh.node_begin(), mesh.node_end(),
                     CS207::DefaultColor(), NodePosition(), node_map);
#endif
    viewer.set_label(t);

    // These lines slow down the animation for small meshes.
    // Feel free to remove them or tweak the constants.
    if (mesh.num_nodes() < 100)
      CS207::sleep(0.05);
  }

  return 0;
}
