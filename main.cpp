/*
 * main.cpp
 *
 *  Created on: 27-Jan-2016
 *      Author: arnab
 */
#include <iostream>
#include <armadillo>
#include <string>

using namespace std;
using namespace arma;

void ReadInput(string filename)
{
	cout << "Not implemented yet!\n";
}

int main()
{
	cout << "Starting program...\n";

	//------Setting parameters--------------------
	float b = 0.95; // for b% singular values
	int r_cap = 0; // Initialized to 0 currently
	int d = 0; // Initialized to 0 currently
	//------Reading input data--------------------
	ReadInput("InputFile.mat");

    //------Initialization part-------------------
	mat X = randu<mat>(5,5);  // This to be replaced by input data matrix M;

	mat U,V,P0_cap,Ptrain_cap,Tt_cap;
	vec s;

	svd(U,s,V,X);

	float sigvalsum = accu(s); // Calculate the sum of the singular values
	float bsigvalsum = 0;
	int bsigvalidx = 0;

	// Computing the index of the b% left singular vectors
	for(unsigned int i=0; i < s.size();++i)
	{
		bsigvalsum +=s.at(i);
		if(bsigvalsum/sigvalsum > b)
		{
			bsigvalidx = i;
			break;
		}
	}

	// Computing initial approximate basis of M
	P0_cap = U.cols(0,bsigvalidx-1);

	// Initializing other parameters and variables
	r_cap = rank(P0_cap);
	d = 10*r_cap;
	Ptrain_cap = P0_cap;

	cout << Ptrain_cap << endl << d << endl;

	return 0;
}



