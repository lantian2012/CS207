/**
 * @file meshed_mass_spring.cpp
 * Implementation of mass-spring system using mesh
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Triangles (one per line) defined by 3 indices into the point list
 * 
 */

#include <fstream>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Mesh.hpp"
#include "Point.hpp"
#include "Meshed_mass_spring.hpp"

typedef Mesh<NodeData, EdgeData, TriData> MeshType;
typedef typename MeshType::node_type Node;
typedef typename MeshType::edge_type Edge;





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
