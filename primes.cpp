#include "CS207/Util.hpp"
#include "math.h"
#include <set>

/** Return true iff @a n is prime.
 * @pre @a n >= 0
 */
bool is_prime(int n)
{
  assert(n >= 0);
  int n_sqrt;
  static std::set<int> primes;
  n_sqrt = (int) (sqrt(n)+0.5);
  if (n < 2)
    return false;
  if (n > 2){
    //check if n is in the prime table
    std::set<int>::iterator it;
    it = primes.find(n);
    if (primes.size() != 0 && it != primes.end())
      return true;
    //if n is not in the prime table, check every number
    int i = 2;
    while (i <= n_sqrt){
      if (n % i == 0)
        return false;
      i++;
    }
    //store the prime 
    primes.insert(n);
    return true;
  }
  //a special case where n==2
  else{
    primes.insert(n);
    return true;
  }
}

int main()
{
  while (!std::cin.eof()) {
    // How many primes to test? And should we print them?
    std::cerr << "Input Number: ";
    int n = 0;
    CS207::getline_parsed(std::cin, n);
    if (n <= 0)
      break;

    std::cerr << "Print Primes (y/n): ";
    char confirm = 'n';
    CS207::getline_parsed(std::cin, confirm);
    bool print_primes = (confirm == 'y' || confirm == 'Y');

    CS207::Clock timer;

    // Loop and count primes from 2 up to n
    int num_primes = 0;
    for (int i = 2; i <= n; ++i) {
      if (is_prime(i)) {
        ++num_primes;
        if (print_primes)
          std::cout << i << std::endl;
      }
    }

    double elapsed_time = timer.seconds();

    std::cout << "There are " << num_primes
              << " primes less than or equal to " << n << ".\n"
              << "Found in " << (1000 * elapsed_time) << " milliseconds.\n\n";
  }

  return 0;
}
