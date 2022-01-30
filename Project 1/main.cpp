#include <iostream>
#include <functional>
#include <cmath>
#include <cassert>
#include "tlinalg.hpp"

// Function declarations
int main();

////////////////////////////////////////////
// PROJECT
// This is the function you need to implement
////////////////////////////////////////////
template <unsigned int n>
vec<n> markov_chain(
  matrix<n, n> A,
  vec<n> v0,
  double eps_step,
  unsigned int max_iterations
);

int main() {
 vec<5> v0{ 1.0, 0.0, 0.0, 0.0, 0.0 };

  matrix<5, 5> A{
    {0.3957, 0.1931, 0.0224, 0.8002, 0.4276},
    {0.8426, 0.4123, 0.9964, 0.3864, 0.6946},
    {0.7730, 0.7306, 0.1065, 0.3964, 0.9449},
    {0.2109, 0.7501, 0.4547, 0.7366, 0.3298},
    {0.6157, 0.8470, 0.4711, 0.3926, 0.8364}
  };

  // This should throw an exception
  try {
    std::cout << markov_chain<5>( A, v0, 1e-5, 100 )
              << std::endl;
  } catch ( std::invalid_argument &e ) {
    std::cout << "A is not stochastic" << std::endl;
  }

  // Make 'A' into a markov_chain matrix
  for ( unsigned int j{ 0 }; j < 5; ++j ) {
    double column_sum{ 0.0 };

    for ( unsigned int i{ 0 }; i < 5; ++i ) {
      column_sum += A( i, j );
    }

    for ( unsigned int i{ 0 }; i < 5; ++i ) {
      A( i, j ) /= column_sum;
    }
  }

  // This should print
  //  [0.139434 0.065835 0.010921 0.295037 0.132249;
  //   0.296910 0.140568 0.485788 0.142467 0.214827;
  //   0.272385 0.249088 0.051923 0.146154 0.292240;
  //   0.074316 0.255736 0.221686 0.271588 0.102001;
  //   0.216956 0.288773 0.229682 0.144753 0.258683]
  std::cout << A << std::endl;

  // This should print
  //     [0.123697 0.247392 0.202221 0.193653 0.233038]'
  std::cout << markov_chain<5>( A, v0, 1e-5, 100 )
            << std::endl;

  // Change 'A' so that the column sums are still 1.0,
  // but there is a negative entry in (0, 0).
  //  - Ethan Maeda noted that the second should be A( 1, 0 )
  A( 0, 0 ) -= 1.1;
  A( 1, 0 ) += 1.1;

  // This should throw an exception
  try {
    std::cout << markov_chain<5>( A, v0, 1e-5, 100 )
              << std::endl;
  } catch ( std::invalid_argument &e ) {
    std::cout << "A is not stochastic" << std::endl;
  }

  // PROJECT Question 5
  //
  matrix<3, 3> B{ {0.084886, 0.360784, 0.257636}, {0.428801, 0.514730, 0.529242}, {0.486313, 0.124486, 0.213122} };   // Stochastic matrix
  vec<3> u3{ 0.2, 0.3, 0.5 };    // Stochastic vector
  std::cout << markov_chain<3>( B, u3, 1e-5, 1000 ) << std::endl;

  return 0;
}

////////////////////////////////////////////
// PROJECT
// This is the function you need to
// implement
////////////////////////////////////////////

template <unsigned int n>
vec<n> markov_chain( matrix<n, n> A, vec<n> v0, double eps_step, unsigned int max_iterations) {
  //iterate through rows
  for (unsigned int i = 0; i < n; i++) {
    float sum = 0;
    //iterate through columns
    for (unsigned int j = 0; j < n; j++) {
      //Check if all entries are non-negative
      if (A(i, j) <= 0) {
        throw std::invalid_argument{"Not all entries are non-negative"};
      }
      else {
        if ( std::fabs( 1.0 - sum )  < eps_step / n ) {
          sum = sum + A(i, j);
        }
      }
    }
  }

  for (unsigned int k = 1; k <= max_iterations; k++) {
		//Vk+1 <--- A*Vk
		vec<n> v1 = A * v0;
    if ( norm(v1 - v0) < eps_step) {
      return v1;
    }
    else {
      v0 = v1;
    }
    if (k == max_iterations) {
      throw std::invalid_argument{"Max iterations reached"};
    }
  }  

  return vec<n>{};
}