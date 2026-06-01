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

int main()
{

	system("clear");
	
	unsigned int i, j, n, k, p;
	unsigned int min, num;
	float J = 1.0, h = 1.0;
	float concurrence, derivative;
	char choice;
	string filename;
	
L:	cout<<"\n\n\n Enter number of spins(between 3 and 12) : ";
	cin>>n;
	k = pow(2, n);													// Dimension of Hilbert space	

	ofstream fcon, fdcon;
	
	vec eigenvalues;
     uvec subsystemInd, sortedIndices;
     vector<unsigned int> subsystemind;
     vector<float> conc, dconc, lambda;
	
	cx_mat I, X, Y, Z, temp, S;
	cx_mat rho12, rhot, con, eigenvectors, rho;
	cx_vec psi;
	
	filename = "con" + to_string(n) + ".dat";
	fcon.open(filename);

	// Paulie Matrices
	I = spm.S(0);											
	X = spm.S(1);
	Y = spm.S(2);
	Z = spm.S(3);
	
	// cout<<I<<"\n "<<X<<"\n "<<Y<<"\n "<<Z<<"\n ";
	
	temp = I;
	
	//cout<<"\n Flag 1";

	cx_mat H, Hz, Hx;												

	H.zeros(k, k);
	//cout<<"\n H = \n"<<H;
	Hz.zeros(k, k);
	//cout<<"\n Hz = \n"<<Hz;
	Hx.zeros(k, k);
	//cout<<"\n Hx = \n"<<Hx;
	

/* Block for defining the Hamiltonian */

	/* For Sigma_Z elements in H */
	for(i = 1; i <= n; i++)					
	{	
		temp = I;

		for(j = 1; j <= n; j++)										// For 1 entire tensor product in the summation
		{
			if(i == j && j == 1)									// 1st element of the 1st tensor product
				temp = Z;

		  	else if(j>1)
		     {
				if(i == j)										// Position of the Z matrix
					temp = arma::kron(Z, temp);

				else												// Filling the rest with I	
					temp = arma::kron(I, temp);
			}

			else
					temp = I;
		}

		Hz = Hz + temp;											// Adding each final tensor product
		
	}
     	//cout<<"\n Hz = \n"<<Hz;
     
     
	/* For sigma_X */    
	S = arma::kron(X, X);
	
	for(i = 1; i < n; i++)				
	{
		temp = I;

		for(j = 1; j < n; j++)										// different tensor-product elements of the summation
		{
			if(i == j && j == 1)									// For the first 2 elements of the first tensor product
				temp = S;

			else if(j > 1)
			{
				if(i == j)
					temp = arma::kron(S,temp);
			
				else
					temp = arma::kron(I,temp);
			}

			else
				temp = I;
		}

		Hx = Hx + temp;
	}
	
		//cout<<"\n Hx = \n"<<Hx;


     //PBC
     temp = X;

     for(j = 2; j < n; j++)
	{
		temp = arma::kron(temp, I);
	}

	temp = arma::kron(temp, X);

	Hx = Hx + temp;

     //cout<<"\n Hx = \n"<<Hx;
    
    
/* Concurrence */
	
	for(num = 3; num <= n; num++)
	{
		subsystemind.push_back(num);									// To append the indices of the subsystem to vec subsystemind; starting from 3
	}
	subsystemInd = conv_to<uvec>::from(subsystemind);						// ??

	for (h = 0.1; h <= 3; h = h + 0.05)								// Running h from 0.1 to 3 with intervals of 0.05
	{
    		H = h*Hz - J*Hx;											// Defining H with J = 1 and h as parameter
			
    		eig_sym(eigenvalues, eigenvectors, H);							// Getting Eigenvectors and Eigenvalues for the corresponding value of h
		
		min = eigenvalues.index_min();								// to get the index of the ground state in Hamiltonian
		psi = eigenvectors.col(min);
		rho = psi*psi.t();											// Computing the outer product
		
		rho12 = TrX(rho, subsystemInd);								// Traces rho over the subsystem subsystemInd. A function from QIClib.
		rhot = kron(Y, Y)*rho12*kron(Y, Y);							// Used in the formula
		con = sqrtmat(sqrtmat(rho12)*rhot*sqrtmat(rho12));				// The formula for concurrence
		
		eig_sym(eigenvalues, eigenvectors, con);						// Calculating the eigenvectors and eigenvalues of "con" matrix
		
		sortedIndices = arma::sort_index(eigenvalues, "descend");			// Sorting eigenvalues in descending order
		eigenvalues = eigenvalues(sortedIndices);
		eigenvectors = eigenvectors.cols(sortedIndices);
   			
  		concurrence = eigenvalues[0] - eigenvalues[2] - eigenvalues[3] - eigenvalues[1];	// because eigenvalues in desceding order			
   						
		if(concurrence < 0)
	    		concurrence = 0;										// By definition			

		fcon<<J/h<<"   "<<concurrence<<endl;							// Store values in the file
	
		//cout<<"\n Flag j"<<i<<endl;
		//cout<<"\n Flag k"<<i<<endl;
		
		conc.push_back(concurrence);
      	lambda.push_back(J/h); 
    			    			
	}
    
   	fcon.close();
   	
   	cout<<"\n Concurrence data for "<<n<<" spin chain written to : "<<filename;


/* Derivative Block */
	filename = "dcon" + to_string(n) + ".dat";
	fdcon.open(filename);
	
	for(i = 1; i < conc.size(); i++)									// 59 elements in conc which will yield 58 averages. But the index goes from 0 to 58
	{
		derivative = (conc[i] - conc[i-1]) / (lambda[i] - lambda[i-1]);
		dconc.push_back(derivative);
	}
	 
	i = 0; 
	for(h = 0.1; h < 2.95; h = h + 0.05)								
		fdcon << (J/h + J/(h+0.05)) / 2 << "   "<<dconc[i++]<<endl;			// Taking average of lambda because we calculated the derivative at the midpoint
	  
	cout<< "\n Data for derivative of concurrence written to : "<<filename; 
	fdcon.close();


/* Choice Block */													// Append at the end of the program after every modification  	
M:   cout<<"\n\n Continue calculation? [y/n] : ";
   	cin>>choice;
   	
   	if(choice == 'y')
   		goto L;
   		
   	else if(choice != 'n' && choice != 'y')
   	{
   		cout<<"\n Invalid Input! Choose again. ";
   		goto M;
   	}
   	
}
