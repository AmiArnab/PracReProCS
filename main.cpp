#include <iostream>
#include <armadillo>
#include <string>
#include <complex>
#include <csv_parser/csv_parser.hpp>

using namespace std;
using namespace arma;

void ReadInput(string filename)
{
    // Needed library to read mat files
    cout << "Can not open " << filename << "\nData reading function not implemented yet!\n";

    //io::CSVReader<1000> csvfile(filename.data());
    csv_parser p;

}

bool NoSupportChange(uvec t_cap_t_2,uvec t_cap_t_1)
{
    int count = 0;
    float frac = 0.0;

    sort(t_cap_t_1.begin(),t_cap_t_1.end());
    sort(t_cap_t_2.begin(),t_cap_t_2.end());

    int t2 = *(max_element(t_cap_t_2.begin(),t_cap_t_2.end()));
    int t1 = *(max_element(t_cap_t_1.begin(),t_cap_t_1.end()));
    int maxitem = (t1>t2)?t1:t2;

    for(int i = 0;i<maxitem;++i)
    {
        bool in_t1 = binary_search(t_cap_t_1.begin(),t_cap_t_1.end(),i);
        bool in_t2 = binary_search(t_cap_t_2.begin(),t_cap_t_2.end(),i);
        if(in_t1 && in_t2)
        {
            ++count;
        }
    }

    frac = (float)count/t_cap_t_2.size();

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
    // Augmented Lagrange Multiplier Method
    double t = 0,t_old = 0,mu=0,rho = 1.5,tau=0, error=0;
    colvec z = Col<double>(y_t.n_cols);
    colvec x = randu<colvec>(y_t.n_cols);//Col<double>(y_t.n_cols);
    colvec z_old = Col<double>(y_t.n_cols);
    colvec x_old = Col<double>(y_t.n_cols);
    colvec y = Col<double>(y_t.n_cols);
    colvec y_old = Col<double>(y_t.n_cols);
    colvec u = Col<double>(y_t.n_cols);
    colvec u_old = Col<double>(y_t.n_cols);
    bool converged1 = false,converged2 = false;

    mat A_temp = Phit_cap.t()*Phit_cap;

    cx_vec eigval = eig_gen(A_temp);
    tau = abs(eigval.max());

    t = 1;
    z = x;
    u = x;
    y = y_t.col(0);// CHANGE!!!!!!
    while(!converged1)
    {
        while(!converged2)
        {
            colvec th_temp = z - (1/tau)*(Phit_cap.t()*((Phit_cap*z)-y_t-((1/mu)*y)));
            double lambda = 1/(mu*tau) * 0.5;
            for(size_t i = 0;i<th_temp.n_rows;++i)
            {
                double tmpval = th_temp.at(i);
                if(tmpval > lambda)
                {
                    th_temp.at(i) = tmpval - lambda;
                }
                else if(tmpval < -lambda)
                {
                    th_temp.at(i) = tmpval + lambda;
                }
            }
            u_old = u;
            u = th_temp;
            t_old = t;
            t = 0.5*(1+sqrt(1+4*t*t));
            z = u + ((t_old-1)/t)*(u-u_old);

            colvec evec2 = y_t - u;
            double error2 = norm(evec2,2);
            converged2 = (error2 < epsilon);

        }
        x = u;
        y = y_old + mu*(y_t-(Phit_cap*x));
        mu = rho*mu;

        colvec evec1 = y_t - u;
        double error1 = norm(evec1,2);
        converged1 = (error1 < epsilon);
    }
    return 0;
}

int Thresh(uvec &Tt,mat St_cs, double omega)
{
    uvec Tt_temp = Col<uword>(St_cs.n_cols);
    uword count = 0;
    for(uword i=0;i<St_cs.n_cols;++i)
    {
        Col<double> tempcol = St_cs.col(i);
        double magofcol = norm(tempcol,2);
        if(magofcol > omega)
        {
            Tt_temp.at(count) = i;
            ++count;
        }
    }
    Tt=Tt_temp.head(count);
    return 0;
}

float SupCardDiff(uvec Tt_capold,uvec Tt_cap)
{
    int count = 0;

    sort(Tt_cap.begin(),Tt_cap.end());
    sort(Tt_capold.begin(),Tt_capold.end());

    int t2 = *(max_element(Tt_capold.begin(),Tt_capold.end()));
    int t1 = *(max_element(Tt_cap.begin(),Tt_cap.end()));
    int maxitem = (t1>t2)?t1:t2;

    for(int i = 0;i<maxitem;++i)
    {
        bool in_t1 = binary_search(Tt_cap.begin(),Tt_cap.end(),i);
        bool in_t2 = binary_search(Tt_capold.begin(),Tt_capold.end(),i);
        if(!in_t1 && in_t2)
        {
            ++count;
        }
    }
    return ((float)count)/Tt_capold.size();
}

int Wl1Minimization(mat St_cs,mat y_t,mat Phit_cap,double epsilon,float lambda)
{
    cout << "Not implemented yet!\n";
    return 0;
}

uword Cardinality(uvec Tt_capold)
{
    return Tt_capold.size();
}

int Prune(uvec &Tadd_cap,mat St_cs,long card)
{
    Row<double> sumofcols = sum(St_cs);
    uvec srtidx = sort_index(sumofcols,1);
    Tadd_cap.clear();

    srtidx.resize(card);
    sort(srtidx.begin(),srtidx.end());

    Tadd_cap = srtidx;

    return 0;
}

int LS(mat Stadd_cap,mat y_t, mat Phit_cap,uvec Tadd_cap)
{
    mat At,sqMat,sqMatInv,psinv,xt,xcap;
    for(uword i=0;i<Tadd_cap.size();++i)
    {
        At.insert_cols(i,Phit_cap.col(Tadd_cap.at(i)));
    }
    sqMat = trans(At)*At;
    sqMatInv = inv(sqMat);
    psinv = sqMatInv*trans(At);
    xt = psinv*y_t;
    xcap = Mat<double>(Phit_cap);
    xcap.zeros();
    for(uword i=0;i<Tadd_cap.size();++i)
    {
        xcap.col(Tadd_cap.at(i)) = xt.at(i);
    }
    Stadd_cap = xcap;
    return 0;
}

int main()
{
    cout << "Starting program...\n";

    //------Setting parameters--------------------
    float b = 0.95; // for b% singular values
    uword r_cap = 0; // Initialized to 0 currently
    int alpha = 20; // Update interval
    long t = 100,t_train=40;

    //------Reading input data--------------------
    ReadInput("InputFile.mat"); //Not Implemented Yet!!

    //------Initialization part-------------------
    mat Mtrain = randu<mat>(5,5);  // This to be replaced by input data matrix M;
    //mat Mt = randu<mat>(5,5);
    colvec Mt;

    mat U,V,P0_cap,Ptrain_cap;
    //mat Lt_cap,Lt_capold, St_cap;
    colvec Lt_cap,Lt_capold, St_cap;
    uvec Tt_cap,Tt_capold;
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
    r_cap = arma::rank(P0_cap);
    Ptrain_cap = P0_cap;
    Tt_cap = zeros<uvec>(1);
    Tt_capold = zeros<uvec>(1);
    Lt_capold = Mtrain.col(Mtrain.n_cols-1);

    cout << "Initialization done!\n";
    cout << "Starting recovery and update phase!\n";

    // When t > t_train
    //while(true)
    //{
        //Perpendicular projection
        mat I = eye(P0_cap.n_rows,P0_cap.n_rows);
        mat Phit_cap = I - (P0_cap * P0_cap.t());
        //mat y_t = Phit_cap * Mt;
        colvec y_t = Phit_cap*Mt;
        //mat St_cs,tempmat, Stadd_cap;
        colvec St_cs,tempmat,Stadd_cap;
        uvec Tadd_cap;
        bool supchange;
        double epsilon, omega, sqomega;
        float lambda;
        uword card;
        supchange = NoSupportChange(Tt_capold,Tt_cap);
        if(supchange)
        {
            tempmat = Phit_cap * Lt_capold;
            epsilon = norm(tempmat,2);
            /*l1Minimization(St_cs,y_t,Phit_cap,epsilon);
            sqomega = (norm(Mt) * norm(Mt))/Mt.n_rows;
            omega = sqrt(sqomega);
            Tt_capold = Tt_cap;
            Thresh(Tt_cap,St_cs,omega);*/
        }
        /*else
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

    //}*/
    return 0;
}
