#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

mat partial_trace(mat density_matrix, int first_spin, int second_spin) 				// Function to get partial trace over all but two spins
{
  mat output_matrix = zeros<mat>(2, 2);								// Initialize the output matrix

  for (int i = 0; i < density_matrix.n_rows; i++) 						// Loop over all spins
  {
    if (i != first_spin && i != second_spin) 							// If the current spin is not the first or second spin
    {
      output_matrix += density_matrix.row(i);							// Update the output matrix
    }
  }
  
  return output_matrix;
}

int main() 
{
 
  mat density_matrix = {0.5, 0.5, 0.5, 0.5};							// Initialize the density matrix

  mat output_matrix = partial_trace(density_matrix, 0, 1);					 // Get the partial trace over all but two spins

  cout << output_matrix << endl;

  return 0;
}

