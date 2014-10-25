/**
 * @file mtl_test.cpp
 * Test script for interfacing with MTL4 and it's linear solvers.
 */

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>


// Define a IdentityMatrix that interfaces with MTL

struct IdentityMatrix
{
	IdentityMatrix(unsigned int col):num_rows_(col), num_cols_(col){
		size_ = num_rows_*num_cols_;
	}

	/** Helper function to perform multiplication. Allows for delayed
	 * evaluation of results.
	 * Assign::apply(a, b) resolves to an assignment operation such as * a += b, a -= b, or a = b.
	 * @pre @a size(v) == size(w) */
	template <typename VectorIn, typename VectorOut, typename Assign>
	void mult(const VectorIn& v, VectorOut& w, Assign) const{
		for (unsigned int i = 0; i < num_rows_; ++i)
		{
	 		Assign::apply(w[i], v[i]);
		}
	}

	/** Matrix - vector multiplication forwards to MTL ’s lazy mat_cvec_multiplier operator */
	template <typename Vector>
	mtl::vec::mat_cvec_multiplier<IdentityMatrix, Vector> operator*(const Vector& v) const {
		return mtl::vec::mat_cvec_multiplier<IdentityMatrix, Vector>(*this, v);
	}
	unsigned int size() const {
		return size_;
	}
	unsigned int num_rows() const{
		return num_rows_;
	}
private:
	unsigned int size_;
	unsigned int num_rows_;
	unsigned int num_cols_;
};
/** The number of elements in the matrix*/
inline std::size_t size(const IdentityMatrix& A){
	return A.size();
}
/** The number of rows in the matrix*/
inline std::size_t num_rows(const IdentityMatrix& A){
	return A.num_rows();
}
/** The number of columns in the matrix*/
inline std::size_t num_cols(const IdentityMatrix& A){
	return A.num_rows();
}

/** Traits that MTL uses to determine properties of our IdentityMatrix . */
namespace mtl{
namespace ashape{
	/**Define IdentityMatrix to be a non-svalar type*/
template<>
struct ashape_aux<IdentityMatrix>
{
	typedef nonscal type;
};
} //end namespace ashape

/**  IdentityMatrix implements the Collection concept
 *  with value_type and size_type*/
template<>
struct Collection<IdentityMatrix>
{
	typedef double value_type;
	typedef unsigned size_type;
};
} //end namesapce mtl


int main()
{
	// HW3: YOUR CODE HERE
	// Construct an IdentityMatrix and "solve" Ix = b
	// using MTL's conjugate gradient solver
	const unsigned int N=15;
	IdentityMatrix I(N);

	// Create a preconditioner
	itl::pc::identity<IdentityMatrix> P(I);

	// Set the b such that x==2.7 is the solution
    mtl::dense_vector<double> x(N, 2.7), b(N);
    b = I * x;

    // start with x=0
    x = 0;

    itl::cyclic_iteration<double> iter(b, 100, 1.e-11, 0.0, 5);
    
    cg(I, x, b, P, iter);
    //output the solve
    std::cout<<x<<std::endl;
    return 0;
}
