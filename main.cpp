#include <iostream>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/calib3d.hpp>
#include <opencv2/imgproc.hpp>
#include <armadillo>
#include <complex>
#include <string>
#include <vector>
#include "csvparser.h"

using namespace std;
using namespace arma;

int ReadInput(string filename, mat &data)
{
    double **mydata;
    int noofrows = 0;
    int noofcols = 0;

    CsvParser *csvparser = CsvParser_new(filename.data(), ",", 0);
    CsvRow *row;
    vector<double *> myrow;
    while ((row = CsvParser_getRow(csvparser)) )
    {
        const char **rowFields = CsvParser_getFields(row);
        noofcols = CsvParser_getNumFields(row);
        double *temp = new double[noofcols];
        for (int i = 0 ; i < noofcols ; i++)
        {
            temp[i] = atof(rowFields[i]);
        }
        myrow.push_back(temp);
        CsvParser_destroy_row(row);
    }

    CsvParser_destroy(csvparser);

    noofrows = myrow.size();
    mydata = new double* [noofrows];

    for(int i = 0;i<noofrows;++i)
    {
        *(mydata+i)  = new double[noofcols];
        for(int j = 0;j<noofcols;++j)
        {
            *((*(mydata+i))+j) = (myrow.at(i))[j];
        }
    }

    mat Mdat = Mat<double>(noofrows,noofcols);

    for(int i = 0;i<noofrows;++i)
    {
        for(int j = 0;j<noofcols;++j)
        {
            Mdat(i,j) = *((*(mydata+i))+j);
        }
    }

    data = Mdat;

    for(int i = 0;i<noofrows;++i)
    {
        delete [] *(mydata+i);
    }

    delete [] mydata;

    return 0;
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

float SupportCardinalityDifference(uvec Tt_capold,uvec Tt_cap)
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

int Thresh(uvec &Tt,colvec St_cs, double omega)
{
    uvec Tt_temp = Col<uword>(St_cs.n_rows);
    uword count = 0;
    for(uword i=0;i<St_cs.n_elem;++i)
    {
        if(St_cs.at(i) > omega)
        {
            Tt_temp.at(count) = i;
            ++count;
        }
    }
    Tt=Tt_temp.head(count);
    return 0;
}

uword Cardinality(uvec Tt_capold)
{
    return Tt_capold.n_elem;
}

int Prune(uvec &Tadd_cap,colvec St_cs,uword card) //Change
{
    uvec T = Col<uword>(St_cs.n_elem);
    T = sort_index(St_cs,"descend");
    T.resize(card);

    Tadd_cap.clear();
    Tadd_cap = T;

    return 0;
}

int LeastSquareEst(colvec &Stadd_cap,colvec y_t, mat Phit_cap,uvec Tadd_cap)
{
    mat At = zeros<mat>(Phit_cap.n_rows,Tadd_cap.n_elem);
    mat sqMat,sqMatInv,psinv;
    colvec xt,xcap;
    for(uword i=0;i<Tadd_cap.n_elem;++i)
    {
        At.col(i) = Phit_cap.col(Tadd_cap.at(i));
    }
    sqMat = trans(At)*At;
    sqMatInv = inv(sqMat);
    psinv = sqMatInv*trans(At);
    xt = psinv*y_t;
    xcap = Col<double>(y_t.n_rows);
    xcap.zeros();
    for(uword i=0;i<Tadd_cap.n_rows;++i)
    {
        xcap.at(Tadd_cap.at(i),0) = xt.at(i,0);
    }
    Stadd_cap = xcap;
    return 0;
}

colvec Shrink(colvec v,double mu)
{
    colvec temp = Col<double>(v.n_rows);
    temp.zeros();
    for(size_t i =0;i<v.n_rows;++i)
    {
        if(v.at(i,0) > mu)
        {
            temp.at(i,0) = v.at(i,0) - mu;
        }
        else if(v.at(i,0) < -mu)
        {
            temp.at(i,0) = v.at(i,0) + mu;
        }
        else
        {
            temp.at(i,0) = 0;
        }
    }
    return temp;
}

int L1Minimization(colvec &St_cs,colvec y_t,mat A,double epsilon)
{
    double error = 100, delta = 0.05, mu = 0.3;
    colvec x = randu<colvec>(y_t.n_rows);
    colvec errmat = Col<double>(y_t.n_rows);
    colvec u = Col<double>(y_t.n_rows);
    colvec v = Col<double>(y_t.n_rows);
    colvec v_old = Col<double>(y_t.n_rows);
    u.zeros();
    v.zeros();
    v_old.zeros();
    while(error > epsilon)
    {
        v = v_old + trans(A)*(y_t - A*u);
        u = delta*Shrink(v,mu);
        errmat = y_t - A*u;
        error = norm(errmat,2);
        v_old = v;
    }
    St_cs = u;
    return 0;
}

int WtdL1Minimization(colvec &St_cs, colvec y_t, mat A, uvec T,double epsilon, float lambda)
{
    double error = 100, delta = 0.05, mu = 0.3;
    colvec x = randu<colvec>(y_t.n_rows);
    colvec errmat = Col<double>(y_t.n_rows);
    colvec u = Col<double>(y_t.n_rows);
    colvec v = Col<double>(y_t.n_rows);
    colvec v_old = Col<double>(y_t.n_rows);
    u.zeros();
    v.zeros();
    v_old.zeros();

    while(error > epsilon)
    {
        v = v_old+trans(A)*(y_t - A*x);
        x = delta*Shrink(v,mu);
        u = x;
        for(uword i=0;i<T.n_elem;++i)
        {
            x.at(T.at(i)) = lambda * x.at(T.at(i));
        }
        errmat = y_t - A*u;
        error = norm(errmat,2);
        v_old = v;
    }
    St_cs = u;

    return 0;
}

void showImage(colvec rawdat,int cols,int rows,string name)
{
    cv::Mat outimg;
    outimg.create(rows,cols,CV_8UC1);
    uword count = 0;
    for(int i=0;i<rows;++i)
    {
        for(int j=0;j<cols;++j)
        {
            outimg.at<uchar>(i,j) = (uchar)((unsigned int)rawdat.at(count));
            count++;
        }
    }
    cv::imshow(name.data(),outimg);
    cv::waitKey(0);
}

int main()
{
    cout << "Starting program...\n";

    //------Setting parameters--------------------
    float b = 0.95; // for b% singular values
    uword r_cap = 0; // Initialized to 0 currently
    int alpha = 20; // Update interval
    uword total = 100;
    uword t=0,d=0;
    mat Mtrain, M, Ltd,Ltdold;

    //------Reading input data--------------------
    ReadInput("curtaintraindata.csv",Mtrain); //Reading training data
    ReadInput("curtainimagedata.csv",M); //Reading image data

    //showImage(Mtrain.col(0),64,80,"background");
    colvec Mt;

    mat U,V,P0_cap,Ptrain_cap;
    colvec Lt_cap,Lt_capold, St_cap, Mtrainmean,Lfinal,Sfinal;
    uvec Tt_cap,Tt_capold;
    vec s;

    //Mtrainmean = mean(Mtrain,1);
    //Mtrain.each_col() -= Mtrainmean;
    //M.each_col() -= Mtrainmean;

    bool decom = svd(U,s,V,Mtrain); //Perform SVD
    cout << "Initial SVD done!\n";

    if(!decom)
    {
        cout << "SVD failed!\nClosing program!\n";
        return 1;
    }

    float sigvalsum = accu(s); // Calculate the sum of the singular values
    float bsigvalsum = 0;
    int bsigvalidx = 0;

    // Computing the index of the b% left singular vectors
    for(unsigned int i=0; i < s.n_elem;++i)
    {
        bsigvalsum +=s.at(i);
        if(bsigvalsum/sigvalsum > b)
        {
            bsigvalidx = i;
            break;
        }
    }
    cout << "Computation of b% left singular values is complete!\n";

    // Computing initial approximate basis of M
    P0_cap = U.cols(0,bsigvalidx-1);
    cout << "Computation of initial approximate basis is complete!\n";

    // Initializing other parameters and variables
    r_cap = arma::rank(P0_cap);
    Ptrain_cap = P0_cap;
    Tt_cap = zeros<uvec>(1);
    Tt_capold = zeros<uvec>(1);
    Lt_capold = Mtrain.col(Mtrain.n_cols-1);
    d = 100 + r_cap; //Should be 10*r_cap

    Ltd = Mtrain.cols(Mtrain.n_cols-d,Mtrain.n_cols-1); //Initialize for updating P in update phase
    Ltdold = Ltd;

    cout << "Initialization done!\nStarting recovery and update phase!\n";

    // When t > t_train
    t = 0; //Keep it in mind!!
    while(t<total)
    {
        //Perpendicular projection
        Mt = M.col(t);
        mat I = eye(P0_cap.n_rows,P0_cap.n_rows);
        mat Phit_cap = I - (P0_cap * P0_cap.t());

        colvec y_t = Phit_cap*Mt;
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
            L1Minimization(St_cs,y_t,Phit_cap,epsilon);
            sqomega = (norm(Mt) * norm(Mt))/Mt.n_rows;
            omega = sqrt(sqomega);
            Tt_capold = Tt_cap;
            Thresh(Tt_cap,St_cs,omega);
        }
        else
        {
            lambda = SupportCardinalityDifference(Tt_capold,Tt_cap);
            tempmat = Phit_cap * Lt_capold;
            epsilon = norm(tempmat,2);
            WtdL1Minimization(St_cs,y_t,Phit_cap,Tt_cap,epsilon,lambda);
            card = Cardinality(Tt_cap);
            Prune(Tadd_cap,St_cs,1.4*card);
            LeastSquareEst(Stadd_cap,y_t,Phit_cap,Tadd_cap);
            sqomega = (norm(Mt) * norm(Mt))/Mt.n_rows;
            omega = sqrt(sqomega);
            Tt_capold = Tt_cap;
            Thresh(Tt_cap,Stadd_cap,omega);
        }

        //Least square estimate
        LeastSquareEst(St_cap,y_t,Phit_cap,Tt_cap);

        //Estimate Lt
        Lt_cap = Mt - St_cap;

        //Update Pt
        if((t/alpha) == 0)
        {

            Ltd.cols(0,Ltd.n_cols-2) = Ltdold.cols(1,Ltdold.n_cols-1);
            Ltd.col(Ltd.n_cols-1) = Lt_cap;

            mat Utemp,Vtemp;
            vec stemp;
            svd(Utemp,stemp,Vtemp,Ltd);

            P0_cap = Utemp.cols(0,r_cap-1);
            Ltdold = Ltd;
        }
        cout << t << endl;
        ++t; // Next iteration
    }
    //showImage(Mtrain.col(0),64,80,"Training");
    Sfinal = St_cap;// + Mtrainmean;
    showImage(Sfinal,64,80,"Sparse");
    //showImage(Lt_cap,64,80,"Dense");
    return 0;
}
