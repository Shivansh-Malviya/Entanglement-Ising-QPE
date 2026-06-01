/* Derivative of the concurrence - To obtain the inflexion point, ergo the quantum phase transition */
#include<iostream>
#include<armadillo>
#include<cmath>
#include<complex>
#include<string>
#include<fstream>
#include<vector>

#define QICLIB_DONT_USE_NLOPT
#include <QIClib>

using namespace std;
using namespace arma;
using namespace qic;

#define cout std::cout


void hz(int n, sp_mat &Hz)
{
/* For Sigma_Z elements in H */

	int i, j;
	sp_mat I(2, 2), Z(2, 2), temp(2, 2);
	 	
	I(0, 0) = 1;								
	I(1, 1) = 1;
	Z(0, 0) = 1;
	Z(1, 1) = -1;
	temp = I;
	
	for(i = 1; i <= n; i++)					
	{	
		//cout<< "\n Hz, first loop " << i;
		temp = I;

		for(j = 1; j <= n; j++)										// For 1 entire tensor product in the summation
		{
			if(i == j && j == 1)									// 1st element of the 1st tensor product
				temp = Z;

		  	else if(j > 1)
		     {
				if(i == j)										// Position of the Z matrix
					temp = kron(Z, temp);

				else												// Filling the rest with I	
					temp = kron(I, temp);
			}

			else
					temp = I;
					
			//cout<< "\n j : " << j << endl;
		}

		Hz = Hz + temp;											// Adding each final tensor product
		
	}
     	//cout<<"\n Hz = \n"<<Hz;
}


void hx(int n, sp_mat &Hx)
{
/* For sigma_X */    
	
	int i, j;
	sp_mat I(2,2), X(2, 2), temp, S;
	
	I(0, 0) = 1;								
	I(1, 1) = 1;
	X(0, 1) = 1;
	X(1, 0) = 1;
	temp = I;
	
	S = kron(X, X);
	
	for(i = 1; i < n; i++)				
	{
		//cout<< "\n Hx, first loop " << i;
		temp = I;

		for(j = 1; j < n; j++)										// different tensor-product elements of the summation
		{
			if(i == j && j == 1)									// For the first 2 elements of the first tensor product
				temp = S;

			else if(j > 1)
			{
				if(i == j)
					temp = kron(S, temp);
			
				else
					temp = kron(I, temp);
			}

			else
				temp = I;
		}

		Hx = Hx + temp;
	}
	
		//cout<<"\n Hx = \n"<<Hx;


     //PBC
     temp = X;

     for(j = 2; j <= n - 1; j++)
		temp = kron(temp, I);

	temp = kron(temp, X);

	Hx += temp;

     //cout<<"\n Hx = \n"<<Hx;
}


cx_mat trace(int k, cx_vec psi)
{
	int r = 4;
	cx_mat rho12;
	rho12.zeros(4, 4);
	
	for(int i = 0; i < (k/r); i++)
	{
		for(int j = 0; j < r; j++)
		{                
			for(int l = 0; l < r; l++)
				rho12((i*r + j)%r , (i*r + l)%r) += psi(i*r + j, 0)*psi(i*r + l, 0);                               
		}
	}
	
	return rho12; 

}




int main()
{

	system("clear");
	
	unsigned int i, j, n, k, p;
	unsigned int min, num;
	float J = 1.0, h = 1.0;
	float concurrence, derivative;
	char choice;
	string filename;
	
	cout<<"\n\n\n Enter number of spins(>= 3) : ";
	cin>>n;
L :	k = pow(2, n);													// Dimension of Hilbert space	
	cout<<"\n\n k = " << k;
	
	ofstream fsp;
	
	vec eigenvalues;
     uvec subsystemInd, sortedIndices;
     vector<unsigned int> subsystemind;
     vector<float> conc, dconc, lambda;
	
	cx_mat rho12, rhot, con, eigenvectors, rho, Y;
	//sp_vec psi(k);
	cx_vec psi;
	sp_mat H(k, k), rho12_sp(4, 4), Hz(k, k), Hx(k, k);	
	mat eigen_vectors;
	
	filename = "sdt_con" + to_string(n) + ".dat";
	fsp.open(filename, ios::trunc);
	//cout<< "\n Filename : "<<filename;

	hx(n, Hx);
	hz(n, Hz);									
    
	auto start = std::chrono::high_resolution_clock::now();
	
/* Concurrence */
	
	Y = spm.S(2);	
	
	for (h = 0.1; h <= 3; h = h + 0.05)								// Running h from 0.1 to 3 with intervals of 0.05
	{
    		H = h*Hz - J*Hx;											// Defining H with J = 1 and h as parameter
			
    		eigs_sym(eigenvalues, eigen_vectors, H, 1, "sa");							// Getting Eigenvectors and Eigenvalues for the corresponding value of h
		
		//min = eigenvalues.index_min();								// to get the index of the ground state in Hamiltonian
		psi = conv_to<cx_vec>::from(eigen_vectors.col(0));				
		
		rho12 = trace(k, psi);							
		rhot = kron(Y, Y) * rho12 * kron(Y, Y);							// Used in the formula
		con = sqrtmat(sqrtmat(rho12) * rhot * sqrtmat(rho12));				// The formula for concurrence
		
		eig_sym(eigenvalues, eigenvectors, con);						// Calculating the eigenvectors and eigenvalues of "con" matrix
   			
  		concurrence = eigenvalues(3) - eigenvalues(2) - eigenvalues(1) - eigenvalues(0);	// because eigenvalues in ascending order			
   						
		if(concurrence < 0)
	    		concurrence = 0;										// By definition			

		fsp<<J/h<<"   "<<concurrence<<endl;							// Store values in the file
	
		//cout<<"\n Flag j"<<i<<endl;
		//cout<<"\n Flag k"<<i<<endl;
		
		conc.push_back(concurrence);
      	lambda.push_back(J/h); 
    			    			
	}
    
   	fsp.close();
   	
   	cout<<"\n Concurrence data for "<<n<<" spin chain written to : "<<filename;


/* Derivative Block */
	filename = "sdt_dcon" + to_string(n) + ".dat";
	fsp.open(filename, ios::trunc);
	
	for(i = 1; i < conc.size(); i++)									// 59 elements in conc which will yield 58 averages. But the index goes from 0 to 58
	{
		derivative = (conc[i] - conc[i-1]) / (lambda[i] - lambda[i-1]);
		dconc.push_back(derivative);
	}
	 
	i = 0; 
	for(h = 0.1; h < 2.95 ; h = h + 0.05)								
		fsp << (J/h + J/(h+0.05)) / 2 << "   "<<dconc[i++]<<endl;			// Taking average of lambda because we calculated the derivative at the midpoint
	  
	cout<< "\n Data for derivative of concurrence written to : "<<filename; 
	fsp.close();
	
	auto end = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> duration = end - start;
	
	cout << "\n Execution time: " << duration.count() << " seconds" << endl;

/* Choice Block 												// Append at the end of the program after every modification  	
M:   cout<<"\n\n Continue calculation? [y/n] : ";
   	cin>>choice;
   	
   	if(choice == 'y')
   		goto L;
   		
   	else if(choice != 'n' && choice != 'y')
   	{
   		cout<<"\n Invalid Input! Choose again. ";
   		goto M;
   	}
*/

	n = n + 2;
	goto L;
   	
}
