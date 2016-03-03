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
	// Needed library to read mat files
	cout << "Data reading function not implemented yet!\n";
}

bool NoSupportChange(mat t_cap_t_2,mat t_cap_t_1)
{
	float noofcomelem = 0, frac = 0;
	int res;
	for(unsigned int i=0;i<t_cap_t_2.n_cols;++i)
	{
		res = 0;
		for(unsigned int j=0;j<t_cap_t_1.n_cols;++j)
		{
			vec t2 = t_cap_t_2.col(i);
			vec t1 = t_cap_t_1.col(j);
			uvec comp = (t2 == t1);
			res = sum(comp);
			if(res != 0) // CHANGE, if equal, then sum = 1, else 0
			{
				noofcomelem = noofcomelem + 1;
			}
		}
	}

	frac = noofcomelem/t_cap_t_1.n_cols;

	if(frac < 0.5 )
	{
		return true;
	}
	else
	{
		return false;
	}
}

int l1Minimization(mat &St_cs, mat y_t, mat Phit_cap,double epsilon)
{
	cout << "Not implemented yet!\n";
	return 0;
}

int Thresh(mat Tt,mat St_cs, double omega)
{
	double max = St_cs.max();
	Tt = clamp(St_cs,omega,max);
	return 0;
}

double SupCardDiff(mat Tt_capold,mat Tt_cap)
{

	int res;
	float cardratio = 0;
	int Ttoldcard = Tt_capold.n_cols;
	int Ttcard = Tt_cap.n_cols;
	int numercard = Ttoldcard;
	for(unsigned int i=0;i<Tt_capold.n_cols;++i)
	{
		res = 0;
		vec t2 = Tt_capold.col(i);
		for(unsigned int j=0;j<Tt_cap.n_cols;++j)
		{
			vec t1 = Tt_cap.col(j);
			uvec comp = (t2 == t1);
			res = sum(comp);
			if(res == comp.size())
			{
				--numercard;
			}
		}
	}
	return numercard/Ttcard;
}

int Wl1Minimization(mat St_cs,mat y_t,mat Phit_cap,double epsilon,double lambda)
{
	cout << "Not implemented yet!\n";
	return 0;
}

int Cardinality(mat Tt_capold)
{
	return Tt_capold.n_cols;
}

int Prune(mat Tadd_cap,mat St_cs,long card)
{
	cout << "Not implemented yet!\n";
	return 0;
}

int LS(mat Stadd_cap,mat y_t, mat Phit_cap,mat Tadd_cap)
{
	cout << "Not implemented yet!\n";
	return 0;
}

int main()
{
	cout << "Starting program...\n";

	//------Setting parameters--------------------
	float b = 0.95; // for b% singular values
	int r_cap = 0; // Initialized to 0 currently
	int d = 0; // Initialized to 0 currently
	int alpha = 20; // Update interval
	long t = 100,t_train=40;

	//------Reading input data--------------------
	ReadInput("InputFile.mat");

    //------Initialization part-------------------
	mat Mtrain = randu<mat>(5,5);  // This to be replaced by input data matrix M;
	mat Mt = randu<mat>(5,5);

	mat U,V,P0_cap,Ptrain_cap,Tt_cap,Tt_capold;
	mat Lt_cap,Lt_capold, St_cap;
	vec s;

	svd(U,s,V,Mtrain);

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
	Tt_cap = zeros<mat>(1,1);
	Tt_capold = zeros<mat>(1,1);

	cout << "Initialization done!\nr_cap = "<< r_cap << endl;
	cout << "Starting recovery and update phase!\n";

	// When t > t_train
	//while(true)
	//{
		//Perpendicular projection
		mat I = eye(P0_cap.n_rows,P0_cap.n_rows);
		mat Phit_cap = I - (P0_cap * P0_cap.t());
		mat y_t = Phit_cap * Mt;
		mat St_cs,Tadd_cap,tempmat, Stadd_cap;
		bool supchange;
		double epsilon, omega, sqomega;
		double lambda;
		long card;
		supchange = NoSupportChange(Tt_capold,Tt_cap);
		if(supchange)
		{
			tempmat = y_t * Lt_capold;
			epsilon = norm(tempmat,2);
			l1Minimization(St_cs,y_t,Phit_cap,epsilon);
			sqomega = (norm(Mt) * norm(Mt))/Mt.n_rows;
			omega = sqrt(sqomega);
			Tt_capold = Tt_cap;
			Thresh(Tt_cap,St_cs,omega);
		}
		else
		{
			lambda = SupCardDiff(Tt_capold,Tt_cap);
			tempmat = y_t * Lt_capold;
			epsilon = norm(tempmat,2);
			Wl1Minimization(St_cs,y_t,Phit_cap,epsilon,lambda);
			card = Cardinality(Tt_cap);
			Prune(Tadd_cap,St_cs,1.4*card);
			LS(Stadd_cap,y_t,Phit_cap,Tadd_cap);
			sqomega = (norm(Mt) * norm(Mt))/Mt.n_rows;
			omega = sqrt(sqomega);
			Tt_capold = Tt_cap;
			Thresh(Tt_cap,Stadd_cap,omega);
		}
		LS(St_cap,y_t,Phit_cap,Tt_cap);

		//Estimate Lt
		Lt_cap = Mt - St_cap;

		//Update Pt

		if((t-t_train)/alpha == 0)
		{
			mat Utemp,Vtemp;
			vec stemp;
			svd(Utemp,stemp,Vtemp,Mt);

			double tempsigvalsum = accu(s); // Calculate the sum of the singular values
			double tempbsigvalsum = 0;
			int tempbsigvalidx = 0;

			for(unsigned int i=0; i < stemp.size();++i)
				{
					tempbsigvalsum +=stemp.at(i);
					if(tempbsigvalsum/tempsigvalsum > r_cap)
					{
						tempbsigvalidx = i;
						break;
					}
				}
			P0_cap = Utemp.cols(0,tempbsigvalidx-1);
		}

	//}
	return 0;
}



