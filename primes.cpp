#include "CS207/Util.hpp"
#include "math.h"

/** Return true iff @a n is prime.
 * @pre @a n >= 0
 * @pre For all i, 0 <= i < n, is_prime(i) has been previously called.
 */
bool is_prime(int n)
{
  assert(n >= 0);
  int n_sqrt;
  static std::vector<int> primes;
  n_sqrt = (int) (sqrt(n)+0.5);
  if (n < 2)
    return false;
  if (n > 2){
  int i = 0;
  while (i < primes.size() && primes[i] <= n_sqrt){
    if (n % primes[i] == 0)
      return false;
    i++;
    }
  //only store the prime if it is not previously stored
  if (n > primes.back())
    primes.push_back(n);
  return true;
  }
  //a special case where n==2
  else{
  if (primes.size() == 0)
    primes.push_back(2);
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
