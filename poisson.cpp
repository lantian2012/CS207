/**
 * @file poisson.cpp
 * Test script for treating the Graph as a MTL Matrix
 * and solving a Poisson equation.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles.
 * Second file: Eges (one per line) defined by 2 indices into the point list
 *              of the first file.
 *
 * Launches an SDLViewer to visualize the solution.
 */

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"

#include "Graph.hpp"
#include "Point.hpp"
#include "BoundingBox.hpp"
#include "math.h"
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>


struct NodeData
{
  double value;
  bool onBoundary;
};

typedef Graph<NodeData,double> GraphType;  //<  DUMMY Placeholder

// HW3: YOUR CODE HERE
// Define a GraphSymmetricMatrix that maps
// your Graph concept to MTL's Matrix concept. This shouldn't need to copy or
// modify Graph at all!
class GraphSymmetricMatrix
{
public:
  GraphSymmetricMatrix(GraphType* graph):g_(graph) {}

  template <typename VectorIn, typename VectorOut, typename Assign>
  void mult(const VectorIn& v, VectorOut& w, Assign) const{
    for (auto it = g_->node_begin(); it != g_->node_end(); ++it)
    {
      unsigned int i = (*it).index();
      double sum = 0;
      if ((*it).value().onBoundary){
        sum = v[i];
      }
      else{
        for (auto j = (*it).edge_begin(); j != (*it).edge_end(); ++j){
          if (!(*j).node2().value().onBoundary){
            sum += v[(*j).node2().index()];
          }
        }
        sum -= (*it).degree() * v[i];
      }
      Assign::apply(w[i], sum);
    }
  }

  template <typename Vector>
  mtl::vec::mat_cvec_multiplier<GraphSymmetricMatrix, Vector> operator*(const Vector& v) const {
    return mtl::vec::mat_cvec_multiplier<GraphSymmetricMatrix, Vector>(*this, v);
  }
  unsigned int size() const{
    return g_->size();
  }

private:
  GraphType* g_;
};
/** The number of elements in the matrix*/
inline std::size_t size(const GraphSymmetricMatrix& A){
  return A.size()*A.size();
}
/** The number of rows in the matrix*/
inline std::size_t num_rows(const GraphSymmetricMatrix& A){
  return A.size();
}
/** The number of columns in the matrix*/
inline std::size_t num_cols(const GraphSymmetricMatrix& A){
  return A.size();
}
namespace mtl{
namespace ashape{
  /**Define GraphSymmetricMatrix to be a non-svalar type*/
template<>
struct ashape_aux<GraphSymmetricMatrix>
{
  typedef nonscal type;
};
}

/**  GraphSymmetricMatrix implements the Collection concept
 *  with value_type and size_type*/
template<>
struct Collection<GraphSymmetricMatrix>
{
  typedef double value_type;
  typedef unsigned size_type;
};
}


/** Remove all the nodes in graph @a g whose posiiton is contained within
 * BoundingBox @a bb
 * @post For all i, 0 <= i < @a g.num_nodes(),
 *        not bb.contains(g.node(i).position())
 */
void remove_box(GraphType& g, const BoundingBox& bb) {
  // HW3: YOUR CODE HERE
  auto it = g.node_begin();
  while (it != g.node_end()){
    if (bb.contains((*it).position()))
    {
      g.remove_node(it);
    }
    else
      ++it;
  }
  return;
}

struct zPosition {
  template <typename NODE>
  Point operator()(const NODE& node) {
    Point pos = node.position();
    pos.elem[2] = node.value().value;
    return pos;
  }
};

struct HeatColor{
  CS207:: Color operator() (GraphType::node_type node){
    return CS207::Color::make_heat(1);
  }
};


   template <typename Real, typename PointFn, typename ColorFn>
   class visual_iteration : public itl::basic_iteration<Real> 
   {
       typedef itl::basic_iteration<Real> super;
       typedef visual_iteration self;
     public:
   
       template <class Vector>
       visual_iteration(GraphType* graph, CS207::SDLViewer* viewer, mtl::dense_vector<double>* x, const Vector& r0, int max_iter_, Real tol_, Real atol_ = Real(0), int cycle_ = 100)
         : super(r0, max_iter_, tol_, atol_), cycle(cycle_), last_print(-1), g_(graph), viewer_(viewer), x_(x)
       {
       }
       

       bool finished() { return super::finished(); }

       template <typename T>
       bool finished(const T& r) 
       {
           bool ret= super::finished(r);
           for (auto it = g_->node_begin(); it != g_->node_end(); ++it){
            (*it).value().value = (*x_)[(*it).index()];
           }
           auto node_map = viewer_->empty_node_map(*g_);
           viewer_->clear();
           viewer_->add_nodes(g_->node_begin(), g_->node_end(), ColorFn(), PointFn(), node_map);
           viewer_->add_edges(g_->edge_begin(), g_->edge_end(), node_map);
           return ret;
       }

       inline self& operator++() { ++this->i; return *this; }
       
       inline self& operator+=(int n) { this->i+= n; return *this; }

       operator int() const { return error_code(); }

       bool is_multi_print() const { return multi_print; }

       void set_multi_print(bool m) { multi_print= m; }

       int error_code() const 
       {
           if (!this->my_suppress)
               std::cout << "finished! error code = " << this->error << '\n'
                   << this->iterations() << " iterations\n"
                   << this->resid() << " is actual final residual. \n"
                   << this->relresid() << " is actual relative tolerance achieved. \n"
                   << "Relative tol: " << this->rtol_ << "  Absolute tol: " << this->atol_ << '\n'
                   << "Convergence:  " << pow(this->relresid(), 1.0 / double(this->iterations())) << std::endl;
           return this->error;
       }
     protected:
       int        cycle, last_print;
       GraphType* g_;
       bool       multi_print = false;
       CS207::SDLViewer* viewer_;
       mtl::dense_vector<double>* x_;
   };


int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE EDGES_FILE\n";
    exit(1);
  }

  // Define an empty Graph
  GraphType graph;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  std::vector<typename GraphType::node_type> node_vec;
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    node_vec.push_back(graph.add_node(2*p - Point(1,1,0)));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CS207::getline_parsed(tets_file, t)) {
    graph.add_edge(node_vec[t[0]], node_vec[t[1]]);
    graph.add_edge(node_vec[t[0]], node_vec[t[2]]);
    graph.add_edge(node_vec[t[1]], node_vec[t[3]]);
    graph.add_edge(node_vec[t[2]], node_vec[t[3]]);
  }

  // Get the edge length, should be the same for each edge
  double h = graph.edge(0).length();

  // Make holes in our Graph
  remove_box(graph, BoundingBox(Point(-0.8+h,-0.8+h,-1), Point(-0.4-h,-0.4-h,1)));
  remove_box(graph, BoundingBox(Point( 0.4+h,-0.8+h,-1), Point( 0.8-h,-0.4-h,1)));
  remove_box(graph, BoundingBox(Point(-0.8+h, 0.4+h,-1), Point(-0.4-h, 0.8-h,1)));
  remove_box(graph, BoundingBox(Point( 0.4+h, 0.4+h,-1), Point( 0.8-h, 0.8-h,1)));
  remove_box(graph, BoundingBox(Point(-0.6+h,-0.2+h,-1), Point( 0.6-h, 0.2-h,1)));

  // HW3: YOUR CODE HERE
  // Define b using the graph, f, and g.
  class Boundary{
  public:
    double operator() (GraphType::node_type n){
      if (norm_inf(n.position()) == 1)
        return 0;
      else if ((norm_inf(n.position()-Point(0.6, 0.6, 0)) < 0.2) || (norm_inf(n.position()-Point(-0.6, -0.6, 0)) < 0.2))
        return -0.2;
      else if (BoundingBox(Point(-0.6,-0.2,-1), Point(0.6,0.2,1)).contains(n.position()))
        return 1;
      else
        return -10;  //A marker to tell a node is NOT on the boundary
    }
  };

  class Force{
  public:
    double operator() (GraphType::node_type n){
      return 5*cos(norm_1(n.position()));
    }
  };
  Boundary g;
  Force f;
  mtl::dense_vector<double> b(graph.size());

  graph.setFlag<Boundary>(g);

  for (auto it = graph.node_begin(); it != graph.node_end(); ++it){
    if ((*it).value().onBoundary)
      b[(*it).index()] = g(*it);
    else{
      b[(*it).index()] = h*h*f(*it);
      for (auto j = (*it).edge_begin(); j != (*it).edge_end(); ++j){
        if ((*j).node2().value().onBoundary)
          b[(*it).index()] -= g((*j).node2());
      }
    }
  }



  // Construct the GraphSymmetricMatrix A using the graph
  GraphSymmetricMatrix A(&graph);
  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(graph);
  viewer.launch();

  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();
  // Solve Au = b using MTL.
  itl::pc::identity<GraphSymmetricMatrix> P(A);
  //itl::cyclic_iteration<double> iter(b, 1000, 1.e-10, 0.0, 50);
  mtl::dense_vector<double> x(graph.size(), 0.0);
  visual_iteration<double, zPosition, HeatColor> iter(&graph, &viewer, &x, b, 1000, 1.e-10, 0.0, 5);
  
  cg(A, x, b, P, iter);

  return 0;
}







