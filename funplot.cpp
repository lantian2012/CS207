#include <vector>
#include <fstream>
#include <queue>
#include <iterator>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Graph.hpp"

/** An iterator that skips over elements of another iterator based on whether
 * those elements satisfy a predicate.
 *
 * Given an iterator range [@a first, @a last) and a predicate @a pred,
 * this iterator models a filtered range such that all i with
 * @a first <= i < @a last and @a pred(*i) appear in order of the original range.
 */
template <typename Pred, typename It>
class filter_iterator
    : private equality_comparable<filter_iterator<Pred,It>> {
 public:
  // Get all of the iterator traits and make them our own
  typedef typename std::iterator_traits<It>::value_type        value_type;
  typedef typename std::iterator_traits<It>::pointer           pointer;
  typedef typename std::iterator_traits<It>::reference         reference;
  typedef typename std::iterator_traits<It>::difference_type   difference_type;
  typedef typename std::input_iterator_tag                     iterator_category;

  typedef filter_iterator<Pred,It> self_type;

  // Constructor
  filter_iterator(const Pred& p, const It& first, const It& last)
      : p_(p), it_(first), end_(last) {
    // HW1 #4: YOUR CODE HERE
    while(it_ != end_ && !p_(*it_))
      ++it_;
  }

  // HW1 #4: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  //return the value_type object pointed by the filter_iterator
  value_type operator*() const{
    return *it_;
  }
  /*return the filter_Iterator that points to the next value_type object
   *and satisfies the predicate
   */
  self_type& operator++(){
    if (it_ == end_)
      return *this;
    do{
      ++it_;
    } while (it_ != end_ && !p_(*it_));
    return *this;
  }
  //True if (it_ == last): When filter_iterator
  //reaches its last position
  /*Test whether this filter_iterator is the same as @a x
   * two filter_iterators are the same if they point to the 
     * same value_type object 
   */
  bool operator==(const self_type& x) const{
    return (it_ == x.it_);
  }

 private:
  Pred p_;
  It it_;
  It end_;
};

/** Helper function for constructing filter_iterators.
 *
 * Usage:
 * // Construct an iterator that filters odd values out and keeps even values.
 * std::vector<int> a = ...;
 * auto it = make_filtered(a.begin(), a.end(), [](int k) {return k % 2 == 0;});
 */
template <typename Pred, typename Iter>
filter_iterator<Pred,Iter> make_filtered(const Iter& it, const Iter& end,
                                         const Pred& p) {
  return filter_iterator<Pred,Iter>(p, it, end);
}

// HW1 #4: YOUR CODE HERE
// Specify and write an interesting predicate on the nodes.
// Explain what your predicate is intended to do and test it.
// If you'd like you may create new nodes and tets files.
//Delete all isolated nodes
//Only show half of structure
struct MyPredicate{
  template <typename NODE>
  bool operator()(const NODE& n) {
    return (n.degree()!=0 && n.position().x<n.position().z);
  }
};

/** Test predicate for HW1 #4 */
struct SlicePredicate {
  template <typename NODE>
  bool operator()(const NODE& n) {
    return n.position().x < 0;
  }
};

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
  std::queue<Point::size_type> BFSqueue;        //Queue for BFS
  std::vector<bool> visited(g.size(), false);      //if (visited[i]==true) node i has been visited. default false
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
  while(!BFSqueue.empty()){
    //dequeue the first element and examine its neighbors
    Point::size_type current = BFSqueue.front();
    BFSqueue.pop();
    for (Graph<int>::incident_iterator it = g.node(current).edge_begin(); it != g.node(current).edge_end(); ++it){
      Point::size_type neighbor = (*it).node2().index();
      if (!visited[neighbor]){
        //get the value for unvisited neighbors
        g.node(neighbor).value() = g.node(current).value() + 1;
        max = g.node(neighbor).value();
        //mark the neighbor as visited
        visited[neighbor] = true;
        //enqueue the neighbor
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
  /* color functor that colors a graph's nodes using a heat map, based on 
   * the integer value stored in the nodes 
   */ 
  class Make_distance{
  public:
    Make_distance(unsigned int max):max_(max){
    }
    CS207::Color operator() (Graph<int>::Node n){
      if (n.value() != -1)
        return CS207::Color::make_heat(double(n.value())/max_);
      else
        return CS207::Color::make_heat(1);
    }
  private:
    unsigned int max_;
  };
  auto node_map = viewer.empty_node_map(graph);
  auto it_begin = make_filtered(graph.node_begin(), graph.node_end(),MyPredicate());
  auto it_end = make_filtered(graph.node_end(), graph.node_end(),MyPredicate());
  viewer.add_nodes(it_begin, it_end, Make_distance(max), node_map);
  viewer.center_view();
  return 0;
}



