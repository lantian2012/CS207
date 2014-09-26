/**
 * @file shortest_path.cpp
 * Test script for using our templated Graph to determine shortest paths.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <vector>
#include <fstream>
#include <queue>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Graph.hpp"


/** Comparator that compares the distance from a given point p.
 */
struct MyComparator {
   Point p_;
   MyComparator(const Point& p) : p_(p) {
   };

   template <typename NODE>
   bool operator()(const NODE& node1, const NODE& node2) const {
    double d1;
    double d2;
    d1 = (node1.position().x - p_.x)*(node1.position().x - p_.x)
           +(node1.position().y - p_.y)*(node1.position().y - p_.y)
           +(node1.position().z - p_.z)*(node1.position().z - p_.z);
    d2 = (node2.position().x - p_.x)*(node2.position().x - p_.x)
           +(node2.position().y - p_.y)*(node2.position().y - p_.y)
           +(node2.position().z - p_.z)*(node2.position().z - p_.z);
    return d1 < d2;
  }
};


/** Calculate shortest path lengths in @a g from the nearest node to @a point.
 * @param[in,out] g Input graph
 * @param[in] point Point to find the nearest node to.
 * @post Graph has modified node values indicating the minimum path length
 *           to the nearest node to @a point
 * @post Graph nodes that are unreachable to the nearest node to @a point have
 *           the value() -1.
 * @return The maximum path length found.
 *
 * Finds the nearest node to @a point and treats that as the root node for a
 * breadth first search.
 * This sets node's value() to the length of the shortest path to
 * the root node. The root's value() is 0. Nodes unreachable from
 * the root have value() -1.
 */
int shortest_path_lengths(Graph<int>& g, const Point& point) {
  // HW1 #4: YOUR CODE HERE
  std::cout<<"entered BFS"<<std::endl;
  std::queue<Point::size_type> BFSqueue;
  std::vector<bool> visited(g.size(), false);      //if visited[i] node i has been visited. default false
  unsigned int max = 0;  //keep track of the maxium distance 
  Graph<int>::node_iterator root_it = std::min_element(g.node_begin(), g.node_end(), MyComparator(point));
  //set all node values to be -1
  for(Graph<int>::node_iterator it=g.node_begin(); it != g.node_end(); ++it){
    (*it).value() = -1;
  }
  //set root node value to be 0
  (*root_it).value() = 0;
  visited[(*root_it).index()] = true;
  BFSqueue.push((*root_it).index());
  std::cout<<"to enter while"<<std::endl;
  while(!BFSqueue.empty()){
    Point::size_type current = BFSqueue.front();
//std::cout<<"current"<<current<<"  "<<g.node(current).degree()<<std::endl;
    BFSqueue.pop();
    for (Graph<int>::incident_iterator it = g.node(current).edge_begin(); it != g.node(current).edge_end(); ++it){
      Point::size_type neighbor = (*it).node2().index();
//std::cout<<"incident:"<<neighbor<<std::endl;
      if (!visited[neighbor]){
        g.node(neighbor).value() = g.node(current).value() + 1;
//std::cout<<"current:"<<current<<"  "<<g.node(current).value()<<"neighbor:"<<neighbor<<std::endl;
        max = g.node(neighbor).value();
        visited[neighbor] = true;
        BFSqueue.push(neighbor);
      }
    }
  }
  return max;
}



int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a Graph
  typedef Graph<int> GraphType;
  GraphType graph;
  std::vector<GraphType::node_type> nodes;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CS207::getline_parsed(tets_file, t))
    for (unsigned i = 1; i < t.size(); ++i)
      for (unsigned j = 0; j < i; ++j)
        graph.add_edge(nodes[t[i]], nodes[t[j]]);

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  viewer.launch();

  // HW1 #4: YOUR CODE HERE
  // Use shortest_path_lengths to set the node values to the path lengths
  // Construct a Color functor and view with the SDLViewer
  unsigned int max; 
  max = shortest_path_lengths(graph, Point(-1, 0, 1));
  std::cout<<"max="<<max<<std::endl;
  class Make_distance{
  public:
    Make_distance(unsigned int max):max_(max){
    }
    CS207::Color operator() (Graph<int>::Node n){
      if (n.value() != -1)
        return CS207::Color::make_heat(1- double(n.value())/max_);
      else
        return CS207::Color::make_heat(1);
    }
  private:
    unsigned int max_;
  };
  auto node_map = viewer.empty_node_map(graph);
  //viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_nodes(graph.node_begin(), graph.node_end(),Make_distance(max), node_map);
  viewer.center_view();
  return 0;
}
