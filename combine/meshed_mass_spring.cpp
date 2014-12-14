/**
 * @file meshed_mass_spring.cpp
 * Implementation of mass-spring system using Mesh
 *
 * @brief Reads in four files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Triangles (one per line) defined by 3 indices into the point list
 * Third file: 3D Points (one per line) defined by three doubles
 * Fourth file: Triangles (one per line) defined by 3 indices into the point list
 */

#include <fstream>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"
#include "Viewer_Extension.hpp"

#include "Mesh.hpp"
#include "Point.hpp"
#include "CollisionDetector.hpp"
#include "Meshed_mass_spring.hpp"


typedef Mesh<NodeData, EdgeData, TriData> MeshType;
typedef typename MeshType::node_type Node;
typedef typename MeshType::edge_type Edge;





int main(int argc, char** argv) {
  // Check arguments
  if (argc < 5) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct the first mesh
  MeshType mesh;
  //construct the second mesh
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
  
  //move the second mesh to the specified position
  for(auto it = mesh2.node_begin();it!=mesh2.node_end();++it){
    (*it).position().elem[1] +=4 ;
    (*it).position().elem[2] +=4 ;
  }

  // Print out the stats
  std::cout << mesh.num_nodes() << " "
            << mesh.num_edges() << " "
            << mesh.num_triangles() << std::endl;
  std::cout << mesh2.num_nodes() << " "
            << mesh2.num_edges() << " "
            << mesh2.num_triangles() << std::endl;



  //set the mass and velocity of each Node for the first mesh
  for (auto it = mesh.node_begin(); it != mesh.node_end(); ++it){
    (*it).value().mass = float(1)/mesh.num_nodes();
    (*it).value().velocity = Point(0, 10, 10);
  }
  //set the mass and velocity of each Node for the second mesh
  for (auto it = mesh2.node_begin(); it != mesh2.node_end(); ++it){
    (*it).value().mass = float(1)/mesh.num_nodes();
    (*it).value().velocity = Point(0, -10, -10);
  }

  //set K and L for each edge
  for (auto it = mesh.node_begin(); it != mesh.node_end(); ++it)
  {
    for (auto j = (*it).edge_begin(); j != (*it).edge_end(); ++j){
       (*j).value().L = (*j).length();
       (*j).value().K = 16000;
    }
  }

  for (auto it = mesh2.node_begin(); it != mesh2.node_end(); ++it)
  {
    for (auto j = (*it).edge_begin(); j != (*it).edge_end(); ++j){
       (*j).value().L = (*j).length();
       (*j).value().K = 16000;
    }
  }


  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(mesh);
  viewer.launch();

  viewer.add_nodes(mesh.node_begin(), mesh.node_end(), node_map);
  viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);
  viewer.add_nodes(mesh2.node_begin(), mesh2.node_end(), node_map);
  viewer.add_edges(mesh2.edge_begin(), mesh2.edge_end(), node_map);

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
  WindForce wind_force(Point(0,100,200));
  PressureForce<typename MeshType::node_type, MeshType> pressure_force(1, 2500, &mesh);
  DampingForce damp_force(float(1)/mesh.num_nodes());
  //auto force = make_combined_force(MassSpringForce(), GravityForce(), make_combined_force(pressure_force, damp_force, wind_force));
  auto force = make_combined_force(MassSpringForce(), make_combined_force(pressure_force, damp_force, wind_force));
  
  //Initialize constriants
  auto constraint = PlaneConstraint<MeshType>(-4);

  //simulation processing   
  for (double t = t_start; t < t_end; t += dt) {
    constraint(mesh, 0);
    constraint(mesh2, 0);
    //define a collision constrain
    auto collision_constrain = CollisionConstraint<MeshType>();
    //add the forces to the meshs at each dt
    symp_euler_step(mesh, t, dt, force);
    symp_euler_step(mesh2, t, dt, force);
    
    //detec the collision betweent the two meshes
    CollisionDetector<MeshType> c;
    c.add_object(mesh);
    c.add_object(mesh2);
    c.check_collisions();
    std::vector<unsigned> collision;
    std::vector<unsigned> collision2;

    //find the corresponding mesh for each node
    for (auto it=c.begin(); it!= c.end(); ++it){
      auto boom = *it;
      Node n = boom.n1;
      if (boom.mesh1 == &mesh)
        collision.push_back(n.index());
      if (boom.mesh1 == &mesh2)
        collision2.push_back(n.index());
    }

    //add the collision constrain to the meshes
    collision_constrain(mesh,mesh2,collision,collision2);
    viewer.set_label(t);
    
    //update with removed nodes
    //clear teh viewer's node
    viewer.clear();
    node_map.clear();
    //update viewer with new positions and new edges
    viewer.add_nodes(mesh.node_begin(), mesh.node_end(), color(color1, color2, color3), node_map);
    viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);
    viewer.add_nodes(mesh2.node_begin(), mesh2.node_end(), color(color1, color2, color3), node_map);
    viewer.add_edges(mesh2.edge_begin(), mesh2.edge_end(), node_map);

    // These lines slow down the animation for small graphs
    if (mesh.num_nodes() < 100)
      CS207::sleep(0.001);
  }

  return 0;
}
