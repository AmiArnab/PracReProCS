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

	ReadInput("InputFile.mat");

	mat X = randu<mat>(5,5);

	mat U;
	vec s;
	mat V;

	svd(U,s,V,X);
	cout << U << endl;
	cout << s << endl;
	cout << V << endl;
	return 0;
}



