#include <iostream>
#include <cmath>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TObject.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <TSpline.h>
#include <vector>
#include <random>
#include <unistd.h>
#include <TMatrixD.h>
#include <chrono>
#include <TROOT.h>
#include <ROOT/TThreadedObject.hxx>
#include <ROOT/TTreeProcessorMT.hxx>
#include <TTreeReaderValue.h>
#include <TDecompLU.h>


std::vector<double> Legendre_coefficients(int m) {
    if (m == 0)
    {
        return {1};
    }
    if (m == 1)
    {
        return {0, 1};
    }

    // Initialize with zero, only at most (half + 1) of the terms will be changed later
    std::vector<double> coeffs(m + 1, 0);

    // Consider some form of memoization instead of this recursion
    std::vector<double> v = Legendre_coefficients(m - 1);
    std::vector<double> u = Legendre_coefficients(m - 2);

    // using literals of floating point type, 'm' is promoted by the compiler
    double a = (2.0 * m - 1.0) / m;
    double b = (m - 1.0) / m;

    int first = 1;
    // If 'm' is even, so is (m - 2) and its first element is zero. It can be skipped.
    // It also avoids annoying '-0' in the output
    if ( m % 2 == 0 )
    {
        coeffs[0] = -b * u[0];
        first = 2;
    }
    for (int i = first; i < m - 1; i += 2)
    {
        coeffs[i] = (a * v[i - 1] - b * u[i]);
    }
    coeffs[m] = a * v[m - 1];

    return coeffs;
}

double ErfInv(double x){
    double tt1, tt2, lnx, sgn;
    sgn = (x < 0) ? -1.0f : 1.0f;

    x = (1 - x)*(1 + x);        // x = 1 - x*x;
    lnx = log(x);

    tt1 = 2/(M_PI*0.147) + 0.5f * lnx;
    tt2 = 1/(0.147) * lnx;

    return(sgn*sqrt(-tt1 + sqrt(tt1*tt1 - tt2)));
}

void invert_matrix(int ncols, TMatrixD *inv_Hankel, std::vector<double> &vector,
                   double &hankel_o, Double_t tol, double *det,
                   int r, int order){
    try {
        auto temp = inv_Hankel->GetSub(0, ncols-1, 0, ncols-1);
        double hankel[ncols*ncols];

        if (!TDecompLU::InvertLU(temp, tol, det)) {
            throw temp.GetNcols()-1;
        }

        temp.GetMatrix2Array(hankel);

        for (int i = 0; i < ncols; ++i) {
            if (std::abs(hankel[i * order + i]) < tol || std::abs(hankel[i * order + i]) > 1.e14) {
                throw temp.GetNcols()-1;
            }
        }

        if (hankel_o != 0.) {
            hankel[0] = -1. / hankel_o;
        }

        for (int i = 0; i < ncols; ++i) {
            for (int j = 0; j < ncols; ++j) {
                vector.at(r * order * order + i * order + j) = hankel[i * order + j];
            }
        }
    }
    catch (int i){
        if (i > 1){
            invert_matrix(i, inv_Hankel, vector,
                          hankel_o, tol, det,
                          r, order);
        }
        else { 
            vector.at(r * order * order) = -1. / hankel_o;
        }
    }
}

void set_inv_hankel(int order, int resolution, std::vector<double> & vector) {
    vector.clear();
    vector.resize(resolution*order*order, 0.);  vector.at(0) = 1.;

    auto *inv_Hankel = new TMatrixD(order, order);
    Double_t inv_hankel[order * order];
    int nr_angles;
    double hankel_o;
    double weight = 2.;
    double *det = nullptr;
    Double_t tol = 2.2204e-16;
    int ncols = order;

    double min_area = M_PI * 1/2.;

    for (int r = 1; r < resolution; ++r) {

        if (r==resolution-1) weight = 1.; else weight = 2.;
        nr_angles = static_cast<int>(M_PI * (r + 1 / 2.) / min_area );

        for (int i = 0; i < order; ++i) {
            for (int j = 0; j < order; ++j) {
                inv_hankel[i * order + j] = 0.;
                for (int ang = 0; ang < nr_angles; ++ang) {
                    inv_hankel[i * order + j] -= weight * pow(std::cos(ang * M_PI / (nr_angles - 1)), i + j);
                }
            }
        }

        hankel_o = inv_hankel[0];
        inv_Hankel->ResizeTo(order,order);
        inv_Hankel->SetMatrixArray(inv_hankel);

        ncols = inv_Hankel->GetNcols();

        invert_matrix(ncols, inv_Hankel, vector, hankel_o, tol, det, r, order);
    }

    inv_Hankel->Clear();
    delete inv_Hankel;
}

void set_inv_hankel_pyabel(int order, int resolution, std::vector<double> & vector){
    vector.clear();
    vector.resize(resolution*order*order,0.);

    std::vector<int> rad_bins(resolution+3);
    for (int i = 0; i < resolution + 3; ++i) {
        rad_bins.at(i) = i;
    }

    int grid_size = resolution+1;
    int index;
    std::vector<double> pc((2*order-1)*(resolution+2),0.);
    double hold_r;

    // Create a grid of grid_sizeXgrid_size of powers of cosine to the angles in the grid
    for (int i = -grid_size; i < grid_size; ++i) {
        for (int j = 0; j < grid_size+1; ++j) {
            hold_r = pow(i*i + j*j,1./2.);

            // Find the bin where this radius lies in
            auto lower_r = std::lower_bound(rad_bins.begin(), rad_bins.end(), hold_r);
            index = (int) std::distance(rad_bins.begin(), lower_r);
            if (index >= resolution+2){
                index = resolution+1;
            }

            // Integrate to the radial bins for each order n
            for (int n = 0; n < 2*order-1; ++n) {
                if (i!=0 && j!=0){
                    pc.at(n * resolution + index - 1) -= pow(i / hold_r, n) * ((double)index-hold_r) * 2;
                    pc.at(n * resolution + index) -= pow(i / hold_r, n) * (1.-((double)index-hold_r)) * 2;
                }
                else if ((j==0 && i!=0) || (j==grid_size && i!=0)) {
                    pc.at(n*resolution + index-1) -= pow(i / hold_r,n) * (1.-((double)index-hold_r));
                    pc.at(n*resolution + index) -= pow(i / hold_r,n) * ((double)index-hold_r);
                }
            }
        }
    }

    for (int r = 0; r < 1; ++r) {
        for (int i = 0; i < order; ++i) {
            for (int j = 0; j < order; ++j) {
                vector.at(r * order * order + i * order + j) = 0.;
            }
        }
        vector.at(r * order * order) = 1.;
    }

    // Now make a hankel matrix for each radius and invert it
    double hankel[order*order];
    auto *inv_Hankel = new TMatrixD(order, order);
    double hankel_o;
    for (int r = 1; r < resolution; ++r) {


        // Make the hankel matrix
        for (int i = 0; i < order; ++i) {
            for (int j = 0; j < order; ++j) {
                if (i == 0 && j == 0 && pc.at((i + j) * resolution + r) == 0.) {
                    hankel[i * order + j] = 1.;
                }
                else {
                    hankel[i * order + j] = pc.at((i + j) * resolution + r);
                }
            }
        }
        hankel_o = hankel[0];

        inv_Hankel->SetMatrixArray(hankel);
        inv_Hankel->Invert();
        inv_Hankel->GetMatrix2Array(hankel);

        if (hankel_o != 0.) {
            hankel[0] = - 1. / hankel_o;
        }

        for (int i = 0; i < order; ++i) {
            for (int j = 0; j < order; ++j) {
                vector.at(r * order * order + i * order + j) = hankel[i * order + j];
            }
        }
    }

    inv_Hankel->Clear();
    delete inv_Hankel;
}

void set_inv_legendre(int order, std::vector<double> & vector) {
    //Initialise to zero
    Double_t  inv_legendre[order*order];
    for (auto & ent : inv_legendre) {
        ent = 0;
    }
    std::vector<double> legendre_coef;
    int count;
    int index;
    for (int col = 0; col < order; ++col) {
        count = 0;
        index = col;
        legendre_coef = Legendre_coefficients(index);
        for (auto & row : legendre_coef){
            inv_legendre[order*count+col] = row;
            count++;
        }
    }
    auto * inv_Legendre = new TMatrixT<double>(order,order);
    inv_Legendre->SetMatrixArray(inv_legendre);
    inv_Legendre->Invert();
    inv_Legendre->GetMatrix2Array(inv_legendre);
    copy(&inv_legendre[0], &inv_legendre[order*order], back_inserter(vector));

    inv_Legendre->Clear();
    delete inv_Legendre;
}

void set_basis_rbasex(int order, int resolution, std::vector<double> & proj_basis) {
    // p_{\rho r;n} basis
    // Initialise array
    proj_basis.resize(order*resolution*resolution);
    for (int n = 0; n < order; ++n) {
        for (int rho = 0; rho < resolution; ++rho){
            for (int r = 0; r < resolution; ++r) {
                if (n==0 && r==0 && rho > 0){
                    proj_basis.at(rho*resolution + r + n*resolution*resolution) = 2;
                }
                else if (r == rho){
                    proj_basis.at(rho*resolution + r + n*resolution*resolution) = 1;
                }
                else {
                    proj_basis.at(rho*resolution + r + n*resolution*resolution) = 0;
                }
            }
        }
    }


    std::vector<double> R;
    std::vector<double> z;
    std::vector<double> frac;
    std::vector<double> frac_n;
    std::vector<double> f;
    std::vector<std::vector<double>> F;
    double r2;

    for (int r = 1; r < resolution + 1; ++r){
        r2 = r*r;

        // Creating all R, rho, z as well as initial F values needed in the sum
        for (int i = r-1; i < resolution + 2; ++i){
            R.emplace_back(i);
            if (i == r-1){
                frac.emplace_back(1);
                z.emplace_back(0);
                f.emplace_back((r * log(r)) / 2.);
            }
            else {
                frac.emplace_back(r / static_cast<double>(i));
                z.emplace_back(pow(i*i-r2,1./2.));
                f.emplace_back((z.at(z.size()-1) / r * i + r * log(z.at(z.size()-1)+i)) / 2.);
            }
        }
        F.emplace_back(f);
        F.emplace_back(z);
        f.clear();

        if (order >= 1){
            // Loop over radial coordinate
            for (int i = r-1; i < resolution + 2; ++i) {
                if (i == r-1){
                    f.emplace_back(r * log(r));
                }
                else {
                    f.emplace_back(r * log(pow(i*i-r2,1./2.)+i));
                }
            }
            F.emplace_back(f);
            f.clear();
        }

        if (order >= 2){
            // Loop over radial coordinate
            for (int i = r-1; i < resolution + 2 ; ++i) {
                if (i == r-1){
                    f.emplace_back(0);
                }
                else {
                    f.emplace_back(r * acos(r / static_cast<double>(i)));
                }
            }
            F.emplace_back(f);
            f.clear();
        }

        if (order >= 3) {
            // Loop over radial coordinate
            for (int i = 0 ; i < resolution + 3 - r; ++i) {
                f.emplace_back(z.at(i) * frac.at(i));
            }
            F.emplace_back(f);
            f.clear();
        }

        if (order >=4){
            frac_n = frac;
            // Loop over all orders >=4, from 4-2 to order-1
            for (int n = 2; n < order-1; ++n) {
                // frac_n = (r / \rho)^n
                for (int i = 0; i < frac_n.size(); ++i) {
                    frac_n.at(i) = frac_n.at(i) * frac.at(i);
                }
                //Calculating F[n+2]
                // Loop over radial coordinate
                for (int i = 0; i < resolution + 3 - r; ++i) {
                    f.emplace_back((z.at(i) * frac_n.at(i) + (n - 1) * F.at(n+1).at(i)) / n);
                }
                F.emplace_back(f);
                f.clear();
            }
        }

        //Compute the projected basis set p_{\rho r;n} basis for all R and n
        for (int n = 1; n < order+1; ++n){
            // Calculating all needed projections
            for (int i = 0; i < R.size(); ++i) {
                f.emplace_back(r * F.at(n-1).at(i) - R.at(i) * F.at(n).at(i));
            }

            // Summing projections to basis
            for (int rho = r; rho < resolution; rho++) {
                proj_basis.at(rho*resolution + r + (n-1)*resolution*resolution) = 2. * (2. * f.at(rho - r + 1) - f.at(rho - r + 2) - f.at(rho - r));
                                                                //rho >= r,              at R                    at R-1              at R+1
            }
            f.clear();
        }
        R.clear();
        z.clear();
        F.clear();
        frac_n.clear();
        frac.clear();
    }
}

void set_inv_abel_transform(int order, int resolution, std::vector<double> & abel_trans_mat) {
    //Calculate the basis set
    abel_trans_mat.resize(order*resolution*resolution);
    std::vector<double> proj_basis;
    set_basis_rbasex(order, resolution,proj_basis);

    auto * inv_pn = new TMatrixT<double>(resolution,resolution);
    Double_t pn[resolution*resolution];

    for (int n = 0; n < order; ++n) {
        // Setting pn
        for (int rho = 0; rho < resolution; ++rho) {
            for (int r = 0; r < resolution; ++r) {
                pn[rho*resolution + r] = proj_basis.at(rho*resolution + r + n*resolution*resolution);
            }
        }
        //Inverting pn
        inv_pn->SetMatrixArray(pn);
        inv_pn->Invert();
        inv_pn->GetMatrix2Array(pn);
        //Setting inverse abel transformation matrix of order n
        for (int rho = 0; rho < resolution; ++rho) {
            for (int r = 0; r < resolution; ++r) {
                abel_trans_mat.at(rho*resolution + r + n*resolution*resolution) = pn[r*resolution + rho];
            }
        }
    }
    inv_pn->Clear();
    delete inv_pn;
}

void get_angular_intensity(int order, int resolution, std::vector<double> & vector,
                           std::vector<double> & hankel, std::vector<double> & abel_trans_mat,
                           std::vector<double> & rad_bin_val, std::vector<double> & data) {
    vector.clear();
    vector.resize(order*resolution);
    double sum;
    double sum2;
    double bin_size = (rad_bin_val.at(resolution-2)-rad_bin_val.at(0))/resolution;
    double bin_size_sq = pow(bin_size,2.);
    for (int n = 0; n < order; ++n) {
        for (int r = 0; r < resolution; ++r) {
            sum = 0;
            for (int j = 0; j < resolution; ++j) {
                for (int i = 0; i < order; ++i) {
                    // sum_ji p_rjn * H_jni * d_ij
                    sum += abel_trans_mat.at(r*resolution + j + n*resolution*resolution) *
                           hankel.at(j*order*order + n*order + i) * data.at(i*resolution + j);
                }
            }
            vector.at(n*resolution + r) = sum;
        }
    }
}

void get_legendre_coef(int order, int resolution, std::vector<double> & vector,
                       std::vector<double> & angular_i, std::vector<double> & legendre) {
    vector.clear();
    vector.resize(order*resolution);
    double sum;
    for (int n = 0; n < order; ++n) {
        for (int r = 0; r < resolution; ++r) {
            sum = 0;
            for (int k = 0; k < order; ++k) {
                for (int j = 0; j < resolution; ++j) {
                    for (int i = 0; i < order; ++i) {
                        // sum_kji C_nk * p_rjk * H_ki * d_ij
                        sum += legendre.at(n*order + k) * angular_i.at(k*resolution + r);
                    }
                }
            }
            vector.at(n*resolution + r) = sum;
        }
    }
}

void get_intensity_err(int order, int resolution, std::vector<double> & vector, std::vector<double> & hankel,
                       std::vector<double> & legendre, std::vector<double> & abel_trans_mat,
                       std::vector<double> & data_err, std::vector<double> & bin_val, std::vector<double> &data) {
    vector.clear();
    vector.resize(resolution);
    double sum;
    double jsum;
    for (int i = 0; i < resolution; ++i) {
        sum = 0;
        for (int n = 0; n < order; ++n) {
            //Sum over the ik'th elemnt of Jacobian
            for (int k = 0; k < resolution; ++k) {
                // calculating Jacobian
                jsum = 0;
                for (int z = 0; z < order; ++z) {
                    jsum += legendre.at(z) * abel_trans_mat.at(i*resolution + k + z*resolution*resolution) *
                            hankel.at(k*order*order + z*order + n);
                }
                sum += jsum*jsum*data_err.at(n*resolution + k);
            }
        }
        vector.at(i) = pow(sum*4*M_PI*pow(bin_val.at(i), 2.),1./2.);
    }
}

void get_center_average(TTree * tree, int & entries,
                        double & x_center, double & y_center, double & center_error, double & x_y_variance, 
                        double & trig_time_upper, double & trig_time_lower,
                        double & event_time_upper, double & event_time_lower,
                        double & peak_sum_upper, double & peak_sum_lower,
                        std::vector<Double_t> * X, std::vector<Double_t> * Y,
                        std::vector<Double_t> *trig_time,
                        std::vector<Double_t> *event_time,
                        std::vector<Double_t> *peak_sum,
                        TBranch *X_branch,
                        TBranch *Y_branch,
                        TBranch *trig_time_branch,
                        TBranch *event_time_branch,
                        TBranch *peak_sum_branch
                        )  {
    double sum_x = 0;
    double count_x = 0;
    double sum_y = 0;
    //Binning x and y-axis
    for (int i = 0; i < entries; i++) {
        X_branch->GetEntry(tree->LoadTree(i));
        Y_branch->GetEntry(tree->LoadTree(i));
        trig_time_branch->GetEntry(tree->LoadTree(i));
        event_time_branch->GetEntry(tree->LoadTree(i));
        peak_sum_branch->GetEntry(tree->LoadTree(i));

        // Bin the data to the resolution
        for (int j = 0; j < X->size(); j++) {
            // Define time window
            if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                event_time->at(j) < event_time_upper && event_time->at(j) > event_time_lower &&
                peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {
                sum_x += X->at(j);
                count_x += 1.;

                sum_y += Y->at(j);
            }
        }
    }

    //Locating center
    x_center = sum_x / count_x; // average x
    y_center = sum_y / count_x; // average y
    center_error = x_y_variance/count_x+x_y_variance; // Variance in center + Variance in data points
}

void get_center_convolution(TTree * tree, int & entries,
                        double & x_center, double & y_center, double & center_error, double & x_y_variance, 
                        double & trig_time_upper, double & trig_time_lower,
                        double & event_time_upper, double & event_time_lower,
                        double & peak_sum_upper, double & peak_sum_lower,
                        std::vector<Double_t> * X, std::vector<Double_t> * Y,
                        std::vector<Double_t> *trig_time,
                        std::vector<Double_t> *event_time,
                        std::vector<Double_t> *peak_sum,
                        TBranch *X_branch,
                        TBranch *Y_branch,
                        TBranch *trig_time_branch,
                        TBranch *event_time_branch,
                        TBranch *peak_sum_branch,
                        Double_t & sq_min, Double_t & sq_max
                        ) {
    std::vector<int> x_axis = {0};
    std::vector<int> y_axis = {0};
    std::vector<Double_t> bins = {10};

    // Create x_axis and y_axis
    for (int i = 0; i < entries; i++) {
        X_branch->GetEntry(tree->LoadTree(i));
        Y_branch->GetEntry(tree->LoadTree(i));
        trig_time_branch->GetEntry(tree->LoadTree(i));
        event_time_branch->GetEntry(tree->LoadTree(i));
        peak_sum_branch->GetEntry(tree->LoadTree(i));

        // Bin the data to the resolution
        for (int j = 0; j < X->size(); j++) {
            // Define time window
            if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                event_time->at(j) < event_time_upper && event_time->at(j) > event_time_lower &&
                peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {
                // Bin the data
                while (X->at(j) >= bins.back() || Y->at(j) >= bins.back()){
                    bins.emplace_back(bins.back() + 0.2);
                    x_axis.emplace_back(0);
                    y_axis.emplace_back(0);
                }
                while (X->at(j) <= bins.at(0) || Y->at(j) <= bins.at(0)) {
                    bins.insert(bins.begin(),bins.at(0) - 0.2);
                    x_axis.insert(x_axis.begin(),0);
                    y_axis.insert(y_axis.begin(),0);
                }
                if ((X->at(j) <= bins.back()) && (X->at(j) >= bins.at(0))) {
                    auto lower = std::lower_bound(bins.begin(), bins.end(), X->at(j));
                    x_axis.at(std::distance(bins.begin(), lower))++;
                }
                if ((Y->at(j) <= bins.back()) && (Y->at(j) >= bins.at(0))) {
                    auto lower = std::lower_bound(bins.begin(), bins.end(), Y->at(j));
                    y_axis.at(std::distance(bins.begin(), lower))++;
                }
            }
        }
    }
    
    int sum;
    int x_max = 0;
    int y_max = 0;

    // Convolution of x_axis with x_axis inverted at point i
    for (int i = 0; i < x_axis.size()*2 ; ++i) {
        sum = 0;
        if (i % 2 == 0){
            for (int j = 0; j <= i; ++j) {
                
                if (j < x_axis.size() && ((i-j) < x_axis.size())){
                    sum += x_axis.at(j) * x_axis.at((i - j));
                }
            }
            if (sum > x_max){
                x_max = sum;
                x_center = bins.at(int(i/2));
            }
        }
    }

    // Convolution of y_axis with y_axis inverted at point i
    for (int i = 0; i < y_axis.size()*2 ; ++i) {
        sum = 0;
        if (i % 2 == 0){
            for (int j = 0; j <= i; ++j) {
                
                if (j < y_axis.size() && ((i-j) < y_axis.size())){
                    sum += y_axis.at(j) * y_axis.at((i - j));
                }
            }
            if (sum > y_max){
                y_max = sum;
                y_center = bins.at(int(i/2));
            }
        }
    }

    center_error = x_y_variance;
}

void set_profile(int resolution, int nr_angles, int ang, std::vector<double> & circ_data, std::vector<double> &profile){
    profile.clear();
    profile.resize(resolution);
    for (int i = 0; i < resolution; ++i) {
        profile.at(i) = circ_data.at(ang*resolution + i);
    }
}

struct least_sq_function {
    double operator()(std::vector<double> &profile_c, std::vector<double> &rad_bin_val_c,
            std::vector<double> &previous_c,
            int &resolution_c,const double *x) const {
        const double rad_corr = x[0];
        const double amp_corr = x[1];

        double radial_f[resolution_c];
        double profile_f[resolution_c];
        for (int i = 0; i < resolution_c; ++i) {
            radial_f[i] = rad_bin_val_c.at(i)*rad_corr;
            profile_f[i] = profile_c.at(i);
            if (std::isnan(radial_f[i])) radial_f[i] = 0.;
            if (std::isnan(profile_f[i])) profile_f[i] = 0;
        }

        auto *spline = new TSpline3("test", radial_f, profile_f, resolution_c,"b1e1");

        //Summing squares of difference
        double sum = 0;
        for (int i = 0; i < resolution_c; ++i) {
            if (std::isnan(rad_bin_val_c.at(i))) {
                sum += pow(previous_c.at(i), 2.);
            }
            else {
                sum += pow(spline->Eval(rad_bin_val_c.at(i))*amp_corr - previous_c.at(i),2.);
            }
        }
        spline->Clear();
        delete spline;
        return sum;
    }
};

void least_sq_min(int resolution, double & rad_corr, double & amp_corr,
                  std::vector<double> &previous, std::vector<double> &profile,
                  std::vector<double> &rad_bin_val) {

    ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit");
    min->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2
    min->SetMaxIterations(100000);  // for GSL
    min->SetTolerance(0.00001);
    min->SetPrintLevel(-1);


    std::function<double(std::vector<double>&, std::vector<double>&, std::vector<double>&,
                int&, const double*)> unboundFct = least_sq_function();
    std::function<double(const double*)> boundFct =
            std::bind(unboundFct,profile,rad_bin_val,previous,resolution, std::placeholders::_1);

    //boundFct(xx);
    ROOT::Math::Functor fct(boundFct,2);

    min->SetFunction(fct);

    double step[2] = {.01, .01};

    // Set the free variables to be minimized !
    min->SetLimitedVariable(0, "rad_corr", rad_corr, step[0], 0.5, 1.5);
    min->SetLimitedVariable(1, "amp_corr", amp_corr, step[1], 0.1, 10);

    // do the minimization
    min->Minimize();

    rad_corr = min->X()[0];
    amp_corr = min->X()[1];

    min->Clear();
}

void get_radial_corrections(int resolution, int nr_angles,
                            std::vector<double> & circ_data, std::vector<double> & rad_bin_val,
                            std::vector<double> & corrections){
    corrections.clear();
    corrections.resize(2*nr_angles); //Here the answer will be appended

    double rad_corr = 1;
    double amp_corr = 1;
    corrections.at(0) = rad_corr;
    corrections.at(nr_angles) = amp_corr;

    std::vector<double> previous;
    std::vector<double> profile;
    set_profile(resolution, nr_angles, 1, circ_data, previous);

    double radial_sol[resolution];
    double profile_sol[resolution];

    //Running over all angles
    for (int ang = 1; ang < nr_angles; ++ang) {
        set_profile(resolution, nr_angles, ang, circ_data, profile);

        //Do a least square fit on profile and the previous profile'
        least_sq_min(resolution, rad_corr, amp_corr, previous, profile, rad_bin_val);

        //Append result
        corrections.at(ang) = rad_corr;
        corrections.at(nr_angles + ang) = amp_corr;

        //Update previous
        for (int i = 0; i < resolution; ++i) {
            radial_sol[i] = rad_bin_val.at(i) * rad_corr;
            profile_sol[i] = profile.at(i);
            if (std::isnan(radial_sol[i])) radial_sol[i] = 0.;
            if (std::isnan(profile_sol[i])) profile_sol[i] = 0.;
        }

        auto *spline = new TSpline3("test", radial_sol, profile_sol, resolution, "b1e1");
        for (int i = 0; i < resolution; ++i) {
            if (std::isnan(rad_bin_val.at(i))) {
                previous.at(i) = spline->Eval(0) * amp_corr;
            }
            else {
                previous.at(i) = spline->Eval(rad_bin_val.at(i)) * amp_corr;
            }
        }
    }
}

void import_radial_corrections(int nr_angles, std::vector<double> & rad_correction, std::string import_file) {
    rad_correction.clear();
    rad_correction.resize(2*nr_angles); //Here the answer will be appended;

    std::ifstream file(import_file);
    std::string line;
    std::string delimiter_char = ",";
    size_t pos;

    int line_count = -1;
    int col_count;
    char * ptr;

    while (getline (file, line)) {
        if (line_count == -1)line_count++;
        else {
            pos = 0;
            col_count = 0;
            while ((pos = line.find(delimiter_char)) != std::string::npos) {
                rad_correction.at(col_count*nr_angles + line_count) = std::strtod(line.substr(0, pos).c_str(), & ptr);
                line.erase(0, pos + delimiter_char.length());
                col_count++;
            }
            line_count++;
        }
    }

    // Close the file
    file.close();

}

int import_nr_angular_bins(std::string import_file) {
    std::ifstream file(import_file);
    std::string line;

    getline(file, line);
    char * ptr;
    // Close the file
    file.close();

    return std::strtol(line.c_str(), & ptr, 10);
}

int main(int argc, char *argv[]){

    auto start_time = std::chrono::steady_clock::now();
    auto end_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_time-start_time;
    std::chrono::duration<double> running_sum = elapsed_seconds;
    auto frame_start_time = std::chrono::steady_clock::now();
    auto frame_end_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> frame_elapsed_seconds = frame_end_time-frame_start_time;
    std::chrono::duration<double> frame_running_sum = frame_elapsed_seconds;

    /*--------            Reading input              --------*/
    std::string root_file;
    std::vector<int> settings_val_vec(5,0);
    std::vector<double> settings_time_vec(7,0);

    //Init values of settings
    settings_val_vec.at(0) = 200; //Resolution
    settings_val_vec.at(1) = 2; // Order
    settings_val_vec.at(2) = 250; // VMI Resolution
    settings_val_vec.at(3) = 0; // Circularize
    settings_val_vec.at(4) = 100; // Bins for Circularize

    settings_time_vec.at(0) = 0; // Trigger time lower
    settings_time_vec.at(1) = 1.e-4; // Trigger time upper
    settings_time_vec.at(2) = 0; // Event time lower
    settings_time_vec.at(3) = 1; // Event time upper
    settings_time_vec.at(4) = 0; // Max radius
    settings_time_vec.at(5) = 0; // Peak sum lower
    settings_time_vec.at(6) = 10000; // Peak sum upper

    // Background mode
    std::vector<double> mode(4,0);
    bool normal = true; mode.at(0) = 1;
    bool no_laser = false;
    bool time_subtraction = false;
    Double_t time_sub_lower = 0;
    Double_t time_sub_upper = 0;

    // Centering mode
    std::string center_mode = "convolution"; // average, convolution
    bool center_mode_changed = false;

    // Settings change
    bool mode_changed = false;
    bool event_time_set = false;
    bool peak_sum_set = false;
    bool resolution_changed = false;
    bool order_changed = false;
    bool vmi_resolution_changed = false;
    bool trigger_time_changed = false;
    bool event_time_changed = false;
    bool peak_sum_changed = false;
    bool circularise_changed = false;
    bool circularise = false;
    bool circularise_imported = false;
    bool bootstrap = false;
    int bootstrap_nr_samp = 0;
    int nr_threads = 1;
    double x_y_variance = pow(0.02*20,2.); // 2 % of the radius of the MCP?
    int gif_nr_frames = 10;
    std::string circularise_import;
    char * ptr;

    if (argc > 1) {
        root_file = argv[1];
    }
    else {
        std::cout << "ERROR" << std::endl;
        std::cout << "Check you specified an argument" << std::endl;
        std::cout << "./Main <data_file_location>" << std::endl;
        assert(true);
    }

    auto *gif_file = new TFile((root_file.substr(0,root_file.size()-5)+std::string("_gif.root")).c_str(), "recreate");
    auto *data_file = new TFile(root_file.c_str(), "read");

    // If settings exist, load them.
    bool setting_exist = data_file->GetListOfKeys()->Contains("SETTINGS_VALUES");
    if (setting_exist){
        auto *settings_val_tree = (TTree *) data_file->Get("SETTINGS_VALUES");
        auto *settings_time_tree = (TTree *) data_file->Get("SETTINGS_TIMES");
        auto *settings_mode_tree = (TTree *) data_file->Get("SETTINGS_MODE");

        int settings_val;
        double settings_time;
        double settings_mode;

        settings_val_tree->SetBranchAddress("settings_val", &settings_val);
        settings_time_tree->SetBranchAddress("settings_time", &settings_time);
        settings_mode_tree->SetBranchAddress("settings_mode", &settings_mode);

        // Setting resolution, order, and VMI resolution
        for (int i = 0; i < settings_val_tree->GetEntries(); ++i) {
            settings_val_tree->GetEntry(i);
            settings_val_vec.at(i) = settings_val;
        }
        // Setting trigger time, and event time
        for (int i = 0; i < settings_time_tree->GetEntries(); ++i) {
            settings_time_tree->GetEntry(i);
            settings_time_vec.at(i) = settings_time;
        }
        // Change back to normal if other mode is not given
        for (int i = 0; i < settings_mode_tree->GetEntries(); ++i){
            settings_mode_tree->GetEntry(i);
            mode.at(i) = settings_mode;
        }
        if (mode.at(0) == 1) normal = true; else no_laser = true;
        if (mode.at(1) == 1) {
            time_subtraction = true;
            time_sub_lower = mode.at(2);
            time_sub_upper = mode.at(3);
        }
        event_time_set = true;
        peak_sum_set = true;
    }

    // Checking the arguments given and checking if they changed from previous settings.
    argc--;
    argv++;
    int c;
    while ((c = getopt(argc, argv, "r:o:v:t:T:e:E:b:c:C:n:m:g:p:P:lLs:S:")) != -1) {
        switch (c) {
            case 'l':
                mode_changed = true;
                no_laser = true;
                mode.at(0) = 0;
                normal = false;
                break;
            case 'L':
                mode_changed = true;
                no_laser = false;
                mode.at(0) = 1;
                normal = true;
                break;
            case 's':
                if (strtod(optarg, &ptr) == -1) {
                    mode_changed = true;
                    time_subtraction = false; mode.at(1) = 0;
                }
                else {
                    time_sub_lower = strtod(optarg, &ptr);
                    mode.at(2) = time_sub_lower;
                    time_subtraction = true;
                    mode.at(1) = 1;
                }
                break;
            case 'S':
                mode_changed = true;
                time_sub_upper = strtod(optarg, & ptr); mode.at(3) = time_sub_upper;
                time_subtraction = true;  mode.at(1) = 1;
                break;
            case 'r':
                settings_val_vec.at(0) = std::strtol(optarg, & ptr, 10);
                resolution_changed = true;
                break;
            case 'o':
                settings_val_vec.at(1) = std::strtol(optarg, & ptr, 10)+1;
                order_changed = true;
                break;
            case 'v':
                settings_val_vec.at(2) = std::strtol(optarg, & ptr, 10);
                vmi_resolution_changed = true;
                break;
            case 't':
                settings_time_vec.at(0) = strtod(optarg, & ptr);
                trigger_time_changed = true;
                break;
            case 'T':
                settings_time_vec.at(1) = strtod(optarg, & ptr);
                trigger_time_changed = true;
                break;
            case 'e':
                settings_time_vec.at(2) = strtod(optarg, & ptr);
                event_time_changed = true;
                break;
            case 'E':
                settings_time_vec.at(3) = strtod(optarg, & ptr);
                event_time_changed = true;
                break;
            case 'p':
                settings_time_vec.at(5) = strtod(optarg, & ptr);
                peak_sum_changed = true;
                break;
            case 'P':
                settings_time_vec.at(6) = strtod(optarg, & ptr);
                peak_sum_changed = true;
                break;
            case 'b':
                bootstrap_nr_samp = std::strtol(optarg, & ptr, 10);
                bootstrap = true;
                break;
            case 'c':
                if (circularise_imported) {
                    break;
                }
                else {
                    settings_val_vec.at(4) = std::strtol(optarg, &ptr, 10);
                    if (settings_val_vec.at(4) > 0) {
                        circularise_changed = true;
                        settings_val_vec.at(3) = true;
                    }
                    if (settings_val_vec.at(4) == -1) {
                        circularise_changed = true;
                        settings_val_vec.at(3) = false;
                    }
                    break;
                }
            case 'C':
                circularise_import = std::string(optarg);
                settings_val_vec.at(4) = import_nr_angular_bins(circularise_import);
                circularise_changed = true;
                circularise_imported = true;
                settings_val_vec.at(3) = true;
                break;
            case 'n':
                nr_threads = std::strtol(optarg, & ptr, 10);
                break;
            case 'm':
                settings_time_vec.at(4) = strtod(optarg, & ptr);
                break;
            case 'g':
                gif_nr_frames = std::strtol(optarg, & ptr, 10);
                break;
            default:
                abort();
        }
    }

    int resolution = settings_val_vec.at(0);
    int norder = settings_val_vec.at(1);
    int VMI_res = settings_val_vec.at(2);
    circularise = settings_val_vec.at(3);
    int nr_angular_bins = settings_val_vec.at(4);

    double trig_time_lower = settings_time_vec.at(0);
    double trig_time_upper = settings_time_vec.at(1);
    double event_time_lower = settings_time_vec.at(2);
    double event_time_upper = settings_time_vec.at(3);
    double max_rad = settings_time_vec.at(4);
    double peak_sum_lower = settings_time_vec.at(5);
    double peak_sum_upper = settings_time_vec.at(6);

    /*--------            Setting up data structures              --------*/

    Int_t injection;
    Bool_t shutter;
    std::vector<Double_t> *trig_time = nullptr;
    std::vector<Double_t> *event_time = nullptr;
    std::vector<Double_t> *peak_sum = nullptr;
    std::vector<Double_t> *X = nullptr;
    std::vector<Double_t> *Y = nullptr;
    Double_t max_val;

    TBranch *injection_branch = nullptr;
    TBranch *shutter_branch = nullptr;
    TBranch *trig_time_branch = nullptr;
    TBranch *event_time_branch = nullptr;
    TBranch *peak_sum_branch = nullptr;
    TBranch *X_branch = nullptr;
    TBranch *Y_branch = nullptr;
    TBranch *max_branch = nullptr;

    // Setting up root file
    auto *tree = (TTree *) data_file->Get("DATA");
    auto *max_tree = (TTree *) data_file->Get("MAXVALUES");
    int entries = static_cast<int>(tree->GetEntries());

    // Set up the branches of the root tree
    tree->SetBranchAddress("injection", &injection, &injection_branch);
    tree->SetBranchAddress("shutter", &shutter, &shutter_branch);
    tree->SetBranchAddress("trig_time", &trig_time, &trig_time_branch);
    tree->SetBranchAddress("event_time", &event_time, &event_time_branch);
    tree->SetBranchAddress("peak_sum", &peak_sum, &peak_sum_branch);
    tree->SetBranchAddress("X", &X, &X_branch);
    tree->SetBranchAddress("Y", &Y, &Y_branch);
    max_tree->SetBranchAddress("max_values", &max_val, &max_branch); //x, y, r, event time


    // If event_time is not set, look for the maximum time
    if (!event_time_set) {
        max_branch->GetEntry(3);
        event_time_upper = max_val;
    }
    if (!peak_sum_set) {
        max_branch->GetEntry(6);
        peak_sum_upper = max_val;
    }

    /*--------          Frame dependencies           --------*/

    std::string frame_name("FRAME");
    double frame_size = (event_time_upper-event_time_lower) / gif_nr_frames;
    double frame_lower;
    double frame_upper;

    /*--------            Radial binning dependencies             --------*/

    //Making raw VMI matrix
    std::vector<int> VMI_on(VMI_res * VMI_res, 0);
    std::vector<int> VMI_off(VMI_res * VMI_res, 0);
    std::vector<int> VMI_fold(VMI_res * VMI_res, 0);
    std::vector<double> VMI_circ(VMI_res * VMI_res, 0);
    std::vector<double> VMI_bin_x(VMI_res);
    std::vector<double> VMI_bin_y(VMI_res);
    int center_col_row = (int)VMI_res/2;

    // Making our data binning n x bins matrix
    std::vector<double> rad_bin_shutter_on;
    std::vector<double> rad_bin_shutter_off;
    std::vector<double> rad_bin_shutter_on_err;
    std::vector<double> rad_bin_shutter_off_err;
    std::vector<double> rad_bin_shutter_on_first;
    std::vector<double> rad_bin_shutter_on_last;
    std::vector<double> rad_bin_shutter_on_err_first;
    std::vector<double> rad_bin_shutter_on_err_last;
    std::vector<double> rad_bin_shutter_off_first;
    std::vector<double> rad_bin_shutter_off_last;
    std::vector<double> rad_bin_shutter_off_err_first;
    std::vector<double> rad_bin_shutter_off_err_last;
    if (!time_subtraction) {
        rad_bin_shutter_on.resize(resolution * norder, 0);
        rad_bin_shutter_off.resize(resolution * norder, 0);
        rad_bin_shutter_on_err.resize(resolution * norder, 0.);
        rad_bin_shutter_off_err.resize(resolution * norder, 0.);
    }
    else {
        if (normal) {
            rad_bin_shutter_on_first.resize(resolution * norder, 0);
            rad_bin_shutter_on_last.resize(resolution * norder, 0);
            rad_bin_shutter_on_err_first.resize(resolution * norder, 0.);
            rad_bin_shutter_on_err_last.resize(resolution * norder, 0.);
        }
        rad_bin_shutter_off_first.resize(resolution * norder, 0);
        rad_bin_shutter_off_last.resize(resolution * norder, 0);
        rad_bin_shutter_off_err_first.resize(resolution * norder, 0.);
        rad_bin_shutter_off_err_last.resize(resolution * norder, 0.);
    }
    std::vector<double> rad_bin_tot(resolution * norder, 0);
    std::vector<double> rad_bin_err(resolution * norder);
    std::vector<double> rad_bin_val(resolution + 1);
    double rad_bin_size;
    Double_t hold_r;
    Double_t hold_y;
    Double_t hold_x;
    Double_t hold_cor;
    Double_t hold_cos;
    Double_t hold_weight;
    double injections_laser_on;
    double injections_laser_off;
    double injections_time_first;
    double injections_time_last;
    double injections_off_time_first;
    double injections_off_time_last;

    // Find the maximal radius in the data, maximum x,y and store data on axis.
    max_branch->GetEntry(1);
    Double_t sq_max = max_val;
    max_branch->GetEntry(0);
    if (max_val > sq_max) {
        sq_max = max_val;
    }
    max_branch->GetEntry(5);
    Double_t sq_min = max_val;
    max_branch->GetEntry(4);
    if (max_val > sq_max) {
        sq_min = max_val;
    }

    /*--------          Finding center dependencies           --------*/
    int index;

    double center_error = 0;

    double x_center;
    double y_center;
    double biggest_diff;

    Double_t rad_max;

    /*--------            Circularize image dependencies        --------*/
    std::vector<double> circ_rad_bin_on;
    std::vector<double> circ_rad_bin_off;
    std::vector<double> circ_rad_bin_on_first;
    std::vector<double> circ_rad_bin_on_last;
    std::vector<double> circ_rad_bin_off_first;
    std::vector<double> circ_rad_bin_off_last;
    std::vector<double> circ_rad_bin_tot;
    std::vector<double> angular_bin_val;
    if (nr_angular_bins == -1)
        angular_bin_val.resize(100);
    else
        angular_bin_val.resize(nr_angular_bins);
    std::vector<double> rad_corrections;
    int ang_index;
    double max_rad_corr;
    if (circularise) {
        if (!time_subtraction) {
            circ_rad_bin_on.resize(nr_angular_bins * resolution, 0);
            circ_rad_bin_off.resize(nr_angular_bins * resolution, 0);
        }
        else {
            if (normal) {
                circ_rad_bin_on_first.resize(nr_angular_bins * resolution, 0);
                circ_rad_bin_on_last.resize(nr_angular_bins * resolution, 0);
            }
            circ_rad_bin_off_first.resize(nr_angular_bins * resolution, 0);
            circ_rad_bin_off_last.resize(nr_angular_bins * resolution, 0);
        }
        circ_rad_bin_tot.resize(nr_angular_bins * resolution, 0);
    }

    /*--------          Abel inversion basis dependencies           --------*/
    std::vector<double> abel_inv_mat; //Getting the inverse abel matrix
    set_inv_abel_transform(norder, resolution, abel_inv_mat);
    std::vector<double> inv_hankel; //Getting the inverse hankel matrix
    std::vector<double> inv_legendre; //Getting the inverse legendre matrix
    set_inv_legendre(norder, inv_legendre);

    std::vector<double> I_n;
    std::vector<double> legendre_coefficients;
    std::vector<double> I(resolution);
    std::vector<double> I_err(resolution);

    std::vector<double> legendre_polynomials(norder*norder,0.);
    std::vector<double> place_holder_vec;

    std::vector<double> inv_vmi(VMI_res * VMI_res, 0.);
    double vmi_relative_res;

    double inv_vmi_rad_max;
    set_inv_hankel(norder,resolution,inv_hankel);

    /*--------          Sampling dependencies           --------*/
    int nr_samples = bootstrap_nr_samp;
    int count;
    int nr_obs = 0;
    double weight_val;
    int lower_index;
    std::vector<double> intensity(resolution);
    std::vector<double> ref_radial;
    std::vector<double> I_iter;
    std::vector<double> z; //Bias correction
    std::vector<double> iteration;
    std::vector<double> I_lower;
    std::vector<double> I_upper;
    if (bootstrap){
        intensity.resize(resolution);
        ref_radial.resize(resolution);
        I_iter.resize(nr_samples*resolution);
        z.resize(resolution); //Bias correction
        iteration.resize(nr_samples);
        I_lower.resize(resolution,0);
        I_upper.resize(resolution,0);
    }

    std::vector<double> weights;
    double new_radial[resolution];
    double new_intensity[resolution];

    // Random generator
    auto gen = std::mt19937{std::random_device{}()};
    std::poisson_distribution<> weight(1);

    double norm_percent;
    TSpline3 *spline;

    // The processor on the dataframe
    ROOT::EnableImplicitMT(nr_threads);
    ROOT::TTreeProcessorMT tp(*tree); //For multi-threading

    /*--------          Iteration on raw data           --------*/
    {
        std::cout << "Finding center and Circularization\n " << std::endl;
        /*--------            Finding center             --------*/

        // Edits x_cent, y_center and center_error
        if (center_mode == "average") {
        get_center_average(tree, entries,
                            x_center, y_center, center_error, x_y_variance, 
                            trig_time_upper, trig_time_lower,
                            event_time_upper, event_time_lower,
                            peak_sum_upper, peak_sum_lower,
                            X, Y, trig_time, event_time, peak_sum,
                            X_branch, Y_branch, trig_time_branch, event_time_branch, peak_sum_branch); 
        }
        if (center_mode == "convolution") {
        get_center_convolution(tree, entries,
                            x_center, y_center, center_error, x_y_variance, 
                            trig_time_upper, trig_time_lower,
                            event_time_upper, event_time_lower,
                            peak_sum_upper, peak_sum_lower,
                            X, Y, trig_time, event_time, peak_sum,
                            X_branch, Y_branch, trig_time_branch, event_time_branch, peak_sum_branch, sq_min, sq_max);
        }
        
        biggest_diff = std::abs(x_center - sq_max);
        if (std::abs(x_center - sq_min) > biggest_diff) {
            biggest_diff = std::abs(x_center - sq_min);
        }
        if (std::abs(y_center - sq_min) > biggest_diff) {
            biggest_diff = std::abs(y_center - sq_min);
        }
        if (std::abs(y_center - sq_max) > biggest_diff) {
            biggest_diff = std::abs(y_center - sq_max);
        }
        
        //Setting maximum radius
        rad_max = biggest_diff * pow(2., 1. / 2.);

        //remaking VMI bins
        for (int i = 0; i < VMI_res; ++i) {
            VMI_bin_x.at(i) = i * (biggest_diff * 2) / (VMI_res - 1) - (biggest_diff);
            VMI_bin_y.at(i) = i * (biggest_diff * 2) / (VMI_res - 1) - (biggest_diff);
        }

        /*--------            Reading from file             --------*/
        rad_bin_size = rad_max / resolution;
        // Make a vector of the radial values
        for (int i = 0; i < resolution+1; i++) {
            rad_bin_val.at(i) = i * rad_bin_size;
        }
        // Make a vector of the angular values
        for (int i = 0; i < nr_angular_bins; i++) {
            angular_bin_val.at(i) = i * 2*M_PI / (nr_angular_bins-1) - M_PI;
        }

        if (circularise) {
            if (circularise_imported) {
                import_radial_corrections(nr_angular_bins, rad_corrections, circularise_import);

                // Getting the maximal radial correction
                max_rad_corr = 1;
                for (int i = 0; i < nr_angular_bins; ++i) {
                    if (max_rad_corr < rad_corrections.at(i))
                        max_rad_corr = rad_corrections.at(i);
                }


                if (max_rad < rad_max && max_rad != 0) rad_max = max_rad;
                // ReMake a vector of the radial values
                rad_bin_size = rad_max * max_rad_corr / resolution;
                for (int i = 0; i < resolution + 1; i++) {
                    rad_bin_val.at(i) = i * rad_bin_size;
                }

                //remaking VMI bins
                for (int i = 0; i < VMI_res; ++i) {
                    VMI_bin_x.at(i) =
                            i * (biggest_diff * max_rad_corr * 2) / (VMI_res - 1) - (biggest_diff * max_rad_corr);
                    VMI_bin_y.at(i) =
                            i * (biggest_diff * max_rad_corr * 2) / (VMI_res - 1) - (biggest_diff * max_rad_corr);
                }
            }
            else {
                /*--------            Circularize image             --------*/
                injections_laser_on = 0;
                injections_laser_off = 0;
                injections_time_first = 0;
                injections_time_last = 0;

                // Do different things for different modes
                if (normal) {
                    if (!time_subtraction) {
                        //Read the data for circularization
                        for (int i = 0; i < entries; i++) {
                            X_branch->GetEntry(tree->LoadTree(i));
                            Y_branch->GetEntry(tree->LoadTree(i));
                            trig_time_branch->GetEntry(tree->LoadTree(i));
                            event_time_branch->GetEntry(tree->LoadTree(i));
                            peak_sum_branch->GetEntry(tree->LoadTree(i));
                            shutter_branch->GetEntry(tree->LoadTree(i));

                            // Bin the data to the resolution
                            if (shutter) {
                                injections_laser_on += 1.;
                                for (int j = 0; j < X->size(); j++) {
                                    // Define time window
                                    if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                        event_time->at(j) < event_time_upper && event_time->at(j) > event_time_lower &&
                                        peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {

                                        hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                        hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                        // Calculate the radius
                                        hold_r = pow(hold_x * hold_x + hold_y * hold_y, 1. / 2.);

                                        // Find the bin where this radius lies in
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r);
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                        if (index >= resolution + 1) {
                                            index = resolution;
                                        }

                                        //Find the bin where the theta lies in
                                        if (hold_x >= 0) { //Upper plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                            angular_bin_val.end(),
                                                                            std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        } else { //Lower plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                            angular_bin_val.end(),
                                                                            -std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        }

                                        // Add to bin
                                        circ_rad_bin_on.at(ang_index * resolution + index - 1) +=
                                                (hold_r - rad_bin_val.at(index - 1)) / rad_bin_size;
                                        circ_rad_bin_on.at(ang_index * resolution + index) +=
                                                (1. + rad_bin_val.at(index - 1) - hold_r) / rad_bin_size;

                                    }
                                }
                            } else {
                                injections_laser_off += 1.;
                                for (int j = 0; j < X->size(); j++) {
                                    // Define time window
                                    if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                        event_time->at(j) < event_time_upper && event_time->at(j) > event_time_lower &&
                                        peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {

                                        hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                        hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                        // Calculate the radius
                                        hold_r = pow(hold_x * hold_x + hold_y * hold_y, 1. / 2.);

                                        // Find the bin where this radius lies in
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r);
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                        if (index >= resolution + 1) {
                                            index = resolution;
                                        }

                                        //Find the bin where the theta lies in
                                        if (hold_x >= 0) { //Upper plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                            angular_bin_val.end(),
                                                                            std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        } else { //Lower plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                            angular_bin_val.end(),
                                                                            -std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        }

                                        // Add to bin
                                        circ_rad_bin_off.at(ang_index * resolution + index - 1) -=
                                                (hold_r - rad_bin_val.at(index - 1)) / rad_bin_size;
                                        circ_rad_bin_off.at(ang_index * resolution + index) -=
                                                (1. + rad_bin_val.at(index - 1) - hold_r) / rad_bin_size;

                                    }
                                }
                            }
                        }
                        // Removing background
                        for (int i = 0; i < resolution * nr_angular_bins; ++i) {
                            circ_rad_bin_tot.at(i) = circ_rad_bin_on.at(i) - injections_laser_on / injections_laser_off *
                                                                             circ_rad_bin_off.at(i);
                        }
                    }
                    else {
                        //Read the data for circularization
                        for (int i = 0; i < entries; i++) {
                            X_branch->GetEntry(tree->LoadTree(i));
                            Y_branch->GetEntry(tree->LoadTree(i));
                            trig_time_branch->GetEntry(tree->LoadTree(i));
                            event_time_branch->GetEntry(tree->LoadTree(i));
                            peak_sum_branch->GetEntry(tree->LoadTree(i));
                            shutter_branch->GetEntry(tree->LoadTree(i));

                            // Bin the data to the resolution
                            if (shutter) {
                                injections_laser_on += 1.;
                                for (int j = 0; j < X->size(); j++) {
                                    // Define time window
                                    if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                        event_time->at(j) < event_time_upper && event_time->at(j) > event_time_lower &&
                                        peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {
                                        injections_time_first += 1;

                                        hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                        hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                        // Calculate the radius
                                        hold_r = pow(hold_x * hold_x + hold_y * hold_y, 1. / 2.);

                                        // Find the bin where this radius lies in
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r);
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                        if (index >= resolution + 1) {
                                            index = resolution;
                                        }

                                        //Find the bin where the theta lies in
                                        if (hold_x >= 0) { //Upper plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                            angular_bin_val.end(),
                                                                            std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        } else { //Lower plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                            angular_bin_val.end(),
                                                                            -std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        }

                                        // Add to bin
                                        circ_rad_bin_on_first.at(ang_index * resolution + index - 1) +=
                                                (hold_r - rad_bin_val.at(index - 1)) / rad_bin_size;
                                        circ_rad_bin_on_first.at(ang_index * resolution + index) +=
                                                (1. + rad_bin_val.at(index - 1) - hold_r) / rad_bin_size;

                                    }
                                    if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                        event_time->at(j) < time_sub_upper && event_time->at(j) > time_sub_lower &&
                                        peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {
                                        injections_time_last += 1;

                                        hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                        hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                        // Calculate the radius
                                        hold_r = pow(hold_x * hold_x + hold_y * hold_y, 1. / 2.);

                                        // Find the bin where this radius lies in
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r);
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                        if (index >= resolution + 1) {
                                            index = resolution;
                                        }

                                        //Find the bin where the theta lies in
                                        if (hold_x >= 0) { //Upper plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                            angular_bin_val.end(),
                                                                            std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        } else { //Lower plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                            angular_bin_val.end(),
                                                                            -std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        }

                                        // Add to bin
                                        circ_rad_bin_on_last.at(ang_index * resolution + index - 1) +=
                                                (hold_r - rad_bin_val.at(index - 1)) / rad_bin_size;
                                        circ_rad_bin_on_last.at(ang_index * resolution + index) +=
                                                (1. + rad_bin_val.at(index - 1) - hold_r) / rad_bin_size;

                                    }
                                }
                            } else {
                                injections_laser_off += 1.;
                                for (int j = 0; j < X->size(); j++) {
                                    // Define time window
                                    if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                        event_time->at(j) < event_time_upper && event_time->at(j) > event_time_lower &&
                                        peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {

                                        hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                        hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                        // Calculate the radius
                                        hold_r = pow(hold_x * hold_x + hold_y * hold_y, 1. / 2.);

                                        // Find the bin where this radius lies in
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r);
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                        if (index >= resolution + 1) {
                                            index = resolution;
                                        }

                                        //Find the bin where the theta lies in
                                        if (hold_x >= 0) { //Upper plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                            angular_bin_val.end(),
                                                                            std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        } else { //Lower plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                            angular_bin_val.end(),
                                                                            -std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        }

                                        // Add to bin
                                        circ_rad_bin_off_first.at(ang_index * resolution + index - 1) -=
                                                (hold_r - rad_bin_val.at(index - 1)) / rad_bin_size;
                                        circ_rad_bin_off_first.at(ang_index * resolution + index) -=
                                                (1. + rad_bin_val.at(index - 1) - hold_r) / rad_bin_size;

                                    }
                                    if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                        event_time->at(j) < time_sub_upper && event_time->at(j) > time_sub_lower &&
                                        peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {

                                        hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                        hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                        // Calculate the radius
                                        hold_r = pow(hold_x * hold_x + hold_y * hold_y, 1. / 2.);

                                        // Find the bin where this radius lies in
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r);
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                        if (index >= resolution + 1) {
                                            index = resolution;
                                        }

                                        //Find the bin where the theta lies in
                                        if (hold_x >= 0) { //Upper plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                            angular_bin_val.end(),
                                                                            std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        } else { //Lower plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                            angular_bin_val.end(),
                                                                            -std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        }

                                        // Add to bin
                                        circ_rad_bin_off_last.at(ang_index * resolution + index - 1) -=
                                                (hold_r - rad_bin_val.at(index - 1)) / rad_bin_size;
                                        circ_rad_bin_off_last.at(ang_index * resolution + index) -=
                                                (1. + rad_bin_val.at(index - 1) - hold_r) / rad_bin_size;

                                    }
                                }
                            }
                        }
                        // Removing background
                        for (int i = 0; i < resolution * nr_angular_bins; ++i) {
                            circ_rad_bin_tot.at(i) = (circ_rad_bin_on_first.at(i) - injections_laser_on / injections_laser_off *
                                                                                    circ_rad_bin_off_first.at(i)) -
                                                     injections_time_first / injections_time_last * (circ_rad_bin_on_last.at(i) - injections_laser_on / injections_laser_off *
                                                                                                                                  circ_rad_bin_off_last.at(i));
                        }
                    }
                }
                else if (no_laser) {
                    if (!time_subtraction) {
                        //Read the data for circularization
                        for (int i = 0; i < entries; i++) {
                            X_branch->GetEntry(tree->LoadTree(i));
                            Y_branch->GetEntry(tree->LoadTree(i));
                            trig_time_branch->GetEntry(tree->LoadTree(i));
                            event_time_branch->GetEntry(tree->LoadTree(i));
                            peak_sum_branch->GetEntry(tree->LoadTree(i));
                            shutter_branch->GetEntry(tree->LoadTree(i));

                            // Bin the data to the resolution
                            if (!shutter) {
                                injections_laser_off += 1.;
                                for (int j = 0; j < X->size(); j++) {
                                    // Define time window
                                    if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                        event_time->at(j) < event_time_upper && event_time->at(j) > event_time_lower &&
                                        peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {

                                        hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                        hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                        // Calculate the radius
                                        hold_r = pow(hold_x * hold_x + hold_y * hold_y, 1. / 2.);

                                        // Find the bin where this radius lies in
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r);
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                        if (index >= resolution + 1) {
                                            index = resolution;
                                        }

                                        //Find the bin where the theta lies in
                                        if (hold_x >= 0) { //Upper plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                            angular_bin_val.end(),
                                                                            std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        } else { //Lower plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                            angular_bin_val.end(),
                                                                            -std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        }

                                        // Add to bin
                                        circ_rad_bin_off.at(ang_index * resolution + index - 1) -=
                                                (hold_r - rad_bin_val.at(index - 1)) / rad_bin_size;
                                        circ_rad_bin_off.at(ang_index * resolution + index) -=
                                                (1. + rad_bin_val.at(index - 1) - hold_r) / rad_bin_size;

                                    }
                                }
                            }
                        }
                        // Removing background
                        for (int i = 0; i < resolution * nr_angular_bins; ++i) {
                            circ_rad_bin_tot.at(i) = circ_rad_bin_off.at(i);
                        }
                    }
                    else {
                        //Read the data for circularization
                        for (int i = 0; i < entries; i++) {
                            X_branch->GetEntry(tree->LoadTree(i));
                            Y_branch->GetEntry(tree->LoadTree(i));
                            trig_time_branch->GetEntry(tree->LoadTree(i));
                            event_time_branch->GetEntry(tree->LoadTree(i));
                            peak_sum_branch->GetEntry(tree->LoadTree(i));
                            shutter_branch->GetEntry(tree->LoadTree(i));

                            // Bin the data to the resolution
                            if (!shutter) {
                                injections_laser_off += 1.;
                                for (int j = 0; j < X->size(); j++) {
                                    // Define time window
                                    if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                        event_time->at(j) < event_time_upper && event_time->at(j) > event_time_lower &&
                                        peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {
                                        injections_time_first += 1;

                                        hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                        hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                        // Calculate the radius
                                        hold_r = pow(hold_x * hold_x + hold_y * hold_y, 1. / 2.);

                                        // Find the bin where this radius lies in
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r);
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                        if (index >= resolution + 1) {
                                            index = resolution;
                                        }

                                        //Find the bin where the theta lies in
                                        if (hold_x >= 0) { //Upper plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                            angular_bin_val.end(),
                                                                            std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        } else { //Lower plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                            angular_bin_val.end(),
                                                                            -std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        }

                                        // Add to bin
                                        circ_rad_bin_off_first.at(ang_index * resolution + index - 1) -=
                                                (hold_r - rad_bin_val.at(index - 1)) / rad_bin_size;
                                        circ_rad_bin_off_first.at(ang_index * resolution + index) -=
                                                (1. + rad_bin_val.at(index - 1) - hold_r) / rad_bin_size;

                                    }
                                    if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                        event_time->at(j) < time_sub_upper && event_time->at(j) > time_sub_lower &&
                                        peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {
                                        injections_time_last += 1;

                                        hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                        hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                        // Calculate the radius
                                        hold_r = pow(hold_x * hold_x + hold_y * hold_y, 1. / 2.);

                                        // Find the bin where this radius lies in
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r);
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                        if (index >= resolution + 1) {
                                            index = resolution;
                                        }

                                        //Find the bin where the theta lies in
                                        if (hold_x >= 0) { //Upper plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                            angular_bin_val.end(),
                                                                            std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        } else { //Lower plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                            angular_bin_val.end(),
                                                                            -std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        }

                                        // Add to bin
                                        circ_rad_bin_off_last.at(ang_index * resolution + index - 1) -=
                                                (hold_r - rad_bin_val.at(index - 1)) / rad_bin_size;
                                        circ_rad_bin_off_last.at(ang_index * resolution + index) -=
                                                (1. + rad_bin_val.at(index - 1) - hold_r) / rad_bin_size;

                                    }
                                }
                            }
                        }
                        // Removing background
                        for (int i = 0; i < resolution * nr_angular_bins; ++i) {
                            circ_rad_bin_tot.at(i) = circ_rad_bin_off_first.at(i) - injections_time_first / injections_time_last * circ_rad_bin_off_last.at(i);
                        }
                    }
                }

                // Get corrections
                get_radial_corrections(resolution, nr_angular_bins, circ_rad_bin_tot, rad_bin_val, rad_corrections);

                // Getting the maximal radial correction
                max_rad_corr = 1;
                for (int i = 0; i < nr_angular_bins; ++i) {
                    if (max_rad_corr < rad_corrections.at(i))
                        max_rad_corr = rad_corrections.at(i);
                }


                if (max_rad < rad_max && max_rad != 0) rad_max = max_rad;
                // ReMake a vector of the radial values
                rad_bin_size = rad_max * max_rad_corr / resolution;
                for (int i = 0; i < resolution + 1; i++) {
                    rad_bin_val.at(i) = i * rad_bin_size;
                }


                //remaking VMI bins
                for (int i = 0; i < VMI_res; ++i) {
                    VMI_bin_x.at(i) =
                            i * (biggest_diff * max_rad_corr * 2) / (VMI_res - 1) - (biggest_diff * max_rad_corr);
                    VMI_bin_y.at(i) =
                            i * (biggest_diff * max_rad_corr * 2) / (VMI_res - 1) - (biggest_diff * max_rad_corr);
                }
            }
        }
        else {
            if (max_rad < rad_max && max_rad!=0) rad_max = max_rad;
            // ReMake a vector of the radial values
            rad_bin_size = rad_max / resolution;
            for (int i = 0; i < resolution+1; i++) {
                rad_bin_val.at(i) = i * rad_bin_size;
            }
        }
    }
    /*--------          Frames           --------*/
    double max_intensity = 0;
    double min_intensity = -1e-16;
    double max_radial_tot = 0;
    double min_radial_tot = -1e-16;
    double max_radius = 0;
    double min_radius = 0;
    std::cout << "Starting frame calculations\n " << std::endl;
    for (int fr = 0; fr < gif_nr_frames; fr++){
        frame_start_time = std::chrono::steady_clock::now();
        data_file->cd();
        frame_lower = fr * frame_size;
        frame_upper = (fr+1) * frame_size;

        /*--------            Reset dependencies             --------*/
        if (!time_subtraction) {
            rad_bin_shutter_on.clear();
            rad_bin_shutter_off.clear();
            rad_bin_shutter_on.resize(resolution * norder, 0);
            rad_bin_shutter_off.resize(resolution * norder, 0);

            rad_bin_shutter_on_err.clear();
            rad_bin_shutter_off_err.clear();
            rad_bin_shutter_on_err.resize(resolution * norder, 0);
            rad_bin_shutter_off_err.resize(resolution * norder, 0);
        }
        else {
            if (normal) {
                rad_bin_shutter_on_first.clear();
                rad_bin_shutter_on_last.clear();
                rad_bin_shutter_on_first.resize(resolution * norder, 0);
                rad_bin_shutter_on_last.resize(resolution * norder, 0);
                rad_bin_shutter_on_err_first.clear();
                rad_bin_shutter_on_err_last.clear();
                rad_bin_shutter_on_err_first.resize(resolution * norder, 0);
                rad_bin_shutter_on_err_last.resize(resolution * norder, 0);
            }
            rad_bin_shutter_off_first.clear();
            rad_bin_shutter_off_last.clear();
            rad_bin_shutter_off_first.resize(resolution * norder, 0);
            rad_bin_shutter_off_last.resize(resolution * norder, 0);
            rad_bin_shutter_off_err_first.clear();
            rad_bin_shutter_off_err_last.clear();
            rad_bin_shutter_off_err_first.resize(resolution * norder, 0);
            rad_bin_shutter_off_err_last.resize(resolution * norder, 0);
        }
        rad_bin_tot.clear();
        rad_bin_val.clear();
        rad_bin_err.clear();
        I_n.clear();
        legendre_coefficients.clear();
        VMI_on.clear();
        VMI_off.clear();
        VMI_fold.clear();
        VMI_circ.clear();
        rad_bin_val.resize(resolution+1);
        rad_bin_tot.resize(resolution * norder, 0);
        rad_bin_err.resize(resolution * norder);
        VMI_on.resize(VMI_res * VMI_res, 0);
        VMI_off.resize(VMI_res * VMI_res, 0);
        VMI_fold.resize(VMI_res * VMI_res, 0);
        VMI_circ.resize(VMI_res * VMI_res, 0);

        injections_laser_on = 0;
        injections_laser_off = 0;
        injections_time_last = 0;
        injections_time_first = 0;
        injections_off_time_last = 0;
        injections_off_time_first = 0;

        /*--------            Reading from file             --------*/
        if (!circularise) {
            // ReMake a vector of the radial values
            rad_bin_size = rad_max / resolution;
            for (int i = 0; i < resolution+1; i++) {
                rad_bin_val.at(i) = i * rad_bin_size;
            }
        }
        else {
            // ReMake a vector of the radial values
            rad_bin_size = rad_max * max_rad_corr / resolution;
            for (int i = 0; i < resolution+1; i++) {
                rad_bin_val.at(i) = i * rad_bin_size;
            }
        }

        // Reading data depending on mode
        if (normal) {
            if (!time_subtraction) {
                //Reading the data
                for (int i = 0; i < entries; i++) {
                    X_branch->GetEntry(tree->LoadTree(i));
                    Y_branch->GetEntry(tree->LoadTree(i));
                    trig_time_branch->GetEntry(tree->LoadTree(i));
                    event_time_branch->GetEntry(tree->LoadTree(i));
                    peak_sum_branch->GetEntry(tree->LoadTree(i));
                    shutter_branch->GetEntry(tree->LoadTree(i));

                    // Bin the data to the resolution
                    if (shutter) {
                        injections_laser_on += 1.;
                        for (int j = 0; j < X->size(); j++) {
                            // Define time window
                            if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                event_time->at(j) < frame_upper && event_time->at(j) > frame_lower &&
                                peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {

                                hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                if (circularise) {
                                    // Calculate the radius and radially correct it
                                    hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                    if (hold_r < rad_max) {

                                        //Find the bin where the theta lies in
                                        if (hold_x >= 0) { //Upper plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                            std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        } else { //Lower plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                            -std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        }

                                        //Add to the VMI image
                                        auto lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                        hold_x);
                                        auto lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                        hold_y);
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_on.at(index) += 1;

                                        //Add to circularised VMI image
                                        lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                   hold_x * rad_corrections.at(ang_index));
                                        lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                   hold_y * rad_corrections.at(ang_index));
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_circ.at(index) += 1;

                                        //Find index correction for integral
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r * rad_corrections.at(ang_index));
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                        if (index < resolution){
                                            hold_cor = (hold_r * rad_corrections.at(ang_index) - rad_bin_val.at(index - 1));
                                            hold_cos = hold_y / hold_r;
                                            // Add sin(theta)cos n(theta) to bin n
                                            for (int n = 0; n < norder; ++n) {
                                                rad_bin_shutter_on.at(n * resolution + index - 1) +=
                                                        hold_cor / rad_bin_size * pow(hold_cos, n);
                                                rad_bin_shutter_on.at(n * resolution + index) +=
                                                        (1. - hold_cor / rad_bin_size) * pow(hold_cos, n);

                                                rad_bin_shutter_on_err.at(n * resolution + index - 1) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * hold_cor+
                                                                                                      pow(hold_cos, 2.) * hold_r * rad_corrections.at(ang_index) -
                                                                                                      n * hold_cor * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                rad_corrections.at(ang_index) * (n-1) * hold_r -
                                                                                n * rad_bin_val.at(index - 1)) ,2.));

                                                rad_bin_shutter_on_err.at(n * resolution + index) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * (hold_cor-rad_bin_size) -
                                                                                                      pow(hold_cos, 2.) * hold_r * rad_corrections.at(ang_index) -
                                                                                                      n * (hold_cor-rad_bin_size) * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                rad_corrections.at(ang_index) * (n-1) * hold_r -
                                                                                n * (rad_bin_val.at(index - 1)+rad_bin_size)) ,2.));
                                            }
                                        }
                                    }

                                }
                                else {

                                    // Calculate the radius
                                    hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                    if (hold_r < rad_max) {

                                        // Find the bin where this radius lies in
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r);
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                        if (index < resolution){
                                            hold_cor = (hold_r - rad_bin_val.at(index - 1));
                                            hold_cos = hold_y / hold_r;
                                            // Add sin(theta)cos n(theta) to bin n
                                            for (int n = 0; n < norder; ++n) {
                                                rad_bin_shutter_on.at(n * resolution + index - 1) +=
                                                        hold_cor / rad_bin_size * pow(hold_cos, n);
                                                rad_bin_shutter_on.at(n * resolution + index) +=
                                                        (1. - hold_cor / rad_bin_size) * pow(hold_cos, n);

                                                rad_bin_shutter_on_err.at(n * resolution + index - 1) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * hold_cor+
                                                                                                      pow(hold_cos, 2.) * hold_r -
                                                                                                      n * hold_cor * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                (n-1) * hold_r -
                                                                                n * rad_bin_val.at(index - 1)) ,2.));

                                                rad_bin_shutter_on_err.at(n * resolution + index) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * (hold_cor-rad_bin_size) -
                                                                                                      pow(hold_cos, 2.) * hold_r -
                                                                                                      n * (hold_cor-rad_bin_size) * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                (n-1) * hold_r -
                                                                                n * (rad_bin_val.at(index - 1)+rad_bin_size)) ,2.));
                                            }
                                        }

                                        //Add to the VMI image
                                        auto lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                        hold_x);
                                        auto lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                        hold_y);
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_on.at(index) += 1;
                                    }
                                }
                            }
                        }
                    } else {
                        injections_laser_off += 1.;
                        for (int j = 0; j < X->size(); j++) {
                            // Define time window
                            if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                event_time->at(j) < frame_upper && event_time->at(j) > frame_lower &&
                                peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {

                                hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                if (circularise) {

                                    // Calculate the radius and radially correct it
                                    hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                    if (hold_r < rad_max*max_rad_corr) {
                                        //Find the bin where the theta lies in
                                        if (hold_x >= 0) { //Upper plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                            std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        } else { //Lower plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                            -std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        }

                                        //Add to the VMI image
                                        auto lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                        hold_x);
                                        auto lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                        hold_y);
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_off.at(index) += 1;

                                        //Add to circularised VMI image
                                        lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                   hold_x * rad_corrections.at(ang_index));
                                        lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                   hold_y * rad_corrections.at(ang_index));
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_circ.at(index) -= 1;

                                        //Find index correction for integral
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r * rad_corrections.at(ang_index));
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                        if (index < resolution){
                                            hold_cor = (hold_r * rad_corrections.at(ang_index) - rad_bin_val.at(index - 1));
                                            hold_cos = hold_y / hold_r;
                                            // Add sin(theta)cos n(theta) to bin n
                                            for (int n = 0; n < norder; ++n) {
                                                rad_bin_shutter_off.at(n * resolution + index - 1) +=
                                                        hold_cor / rad_bin_size * pow(hold_cos, n);
                                                rad_bin_shutter_off.at(n * resolution + index) +=
                                                        (1. - hold_cor / rad_bin_size) * pow(hold_cos, n);

                                                rad_bin_shutter_off_err.at(n * resolution + index - 1) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * hold_cor+
                                                                                                      pow(hold_cos, 2.) * hold_r * rad_corrections.at(ang_index) -
                                                                                                      n * hold_cor * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                rad_corrections.at(ang_index) * (n-1) * hold_r -
                                                                                n * rad_bin_val.at(index - 1)) ,2.));

                                                rad_bin_shutter_off_err.at(n * resolution + index) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * (hold_cor-rad_bin_size) -
                                                                                                      pow(hold_cos, 2.) * hold_r * rad_corrections.at(ang_index) -
                                                                                                      n * (hold_cor-rad_bin_size) * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                rad_corrections.at(ang_index) * (n-1) * hold_r -
                                                                                n * (rad_bin_val.at(index - 1)+rad_bin_size)) ,2.));
                                            }
                                        }
                                    }

                                }
                                else {
                                    // Calculate the radius
                                    hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                    if (hold_r < rad_max) {

                                        // Find the bin where this radius lies ini
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r);
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                        if (index < resolution){
                                            hold_cor = (hold_r - rad_bin_val.at(index - 1));
                                            hold_cos = hold_y / hold_r;
                                            // Add sin(theta)cos n(theta) to bin n
                                            for (int n = 0; n < norder; ++n) {
                                                rad_bin_shutter_off.at(n * resolution + index - 1) +=
                                                        hold_cor / rad_bin_size * pow(hold_cos, n);
                                                rad_bin_shutter_off.at(n * resolution + index) +=
                                                        (1. - hold_cor / rad_bin_size) * pow(hold_cos, n);

                                                rad_bin_shutter_off_err.at(n * resolution + index - 1) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * hold_cor+
                                                                                                      pow(hold_cos, 2.) * hold_r -
                                                                                                      n * hold_cor * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                (n-1) * hold_r -
                                                                                n * rad_bin_val.at(index - 1)) ,2.));

                                                rad_bin_shutter_off_err.at(n * resolution + index) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * (hold_cor-rad_bin_size) -
                                                                                                      pow(hold_cos, 2.) * hold_r -
                                                                                                      n * (hold_cor-rad_bin_size) * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                (n-1) * hold_r -
                                                                                n * (rad_bin_val.at(index - 1)+rad_bin_size)) ,2.));
                                            }
                                        }

                                        //Add to the VMI image
                                        auto lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                        hold_x);
                                        auto lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                        hold_y);
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_off.at(index) += 1;
                                    }
                                }
                            }
                        }
                    }
                }

                // Removing background and saving the data
                for (int i = 0; i < resolution * norder; ++i) {
                    // Sum, r_on - inject_on / inject_off * r_off
                    rad_bin_tot.at(i) = rad_bin_shutter_on.at(i) - injections_laser_on / injections_laser_off *
                                                                   rad_bin_shutter_off.at(i);

                    // Error in this context means the variance*
                    // Error propagation of the sum above, r_on_err + (inject_on / inject_off)^2 * r_off_err
                    // Assume no error in injection count
                    rad_bin_err.at(i) = rad_bin_shutter_on_err.at(i) + pow(injections_laser_on / injections_laser_off,2.) *
                                                                       rad_bin_shutter_off_err.at(i);
                }
            }
            else {
                //Reading the data
                for (int i = 0; i < entries; i++) {
                    X_branch->GetEntry(tree->LoadTree(i));
                    Y_branch->GetEntry(tree->LoadTree(i));
                    trig_time_branch->GetEntry(tree->LoadTree(i));
                    event_time_branch->GetEntry(tree->LoadTree(i));
                    peak_sum_branch->GetEntry(tree->LoadTree(i));
                    shutter_branch->GetEntry(tree->LoadTree(i));

                    // Bin the data to the resolution
                    if (shutter) {
                        injections_laser_on += 1.;
                        for (int j = 0; j < X->size(); j++) {
                            // Define time window
                            if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                event_time->at(j) < frame_upper && event_time->at(j) > frame_lower &&
                                peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {
                                injections_time_first += 1;

                                hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                if (circularise) {
                                    // Calculate the radius and radially correct it
                                    hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                    if (hold_r < rad_max) {

                                        //Find the bin where the theta lies in
                                        if (hold_x >= 0) { //Upper plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                            std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        } else { //Lower plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                            -std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        }

                                        //Add to the VMI image
                                        auto lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                        hold_x);
                                        auto lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                        hold_y);
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_on.at(index) += 1;

                                        //Add to circularised VMI image
                                        lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                   hold_x * rad_corrections.at(ang_index));
                                        lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                   hold_y * rad_corrections.at(ang_index));
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_circ.at(index) += 1;

                                        //Find index correction for integral
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r * rad_corrections.at(ang_index));
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                        if (index < resolution){
                                            hold_cor = (hold_r * rad_corrections.at(ang_index) - rad_bin_val.at(index - 1));
                                            hold_cos = hold_y / hold_r;
                                            // Add sin(theta)cos n(theta) to bin n
                                            for (int n = 0; n < norder; ++n) {
                                                rad_bin_shutter_on_first.at(n * resolution + index - 1) +=
                                                        hold_cor / rad_bin_size * pow(hold_cos, n);
                                                rad_bin_shutter_on_first.at(n * resolution + index) +=
                                                        (1. - hold_cor / rad_bin_size) * pow(hold_cos, n);

                                                rad_bin_shutter_on_err_first.at(n * resolution + index - 1) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * hold_cor+
                                                                                                      pow(hold_cos, 2.) * hold_r * rad_corrections.at(ang_index) -
                                                                                                      n * hold_cor * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                rad_corrections.at(ang_index) * (n-1) * hold_r -
                                                                                n * rad_bin_val.at(index - 1)) ,2.));

                                                rad_bin_shutter_on_err_first.at(n * resolution + index) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * (hold_cor-rad_bin_size) -
                                                                                                      pow(hold_cos, 2.) * hold_r * rad_corrections.at(ang_index) -
                                                                                                      n * (hold_cor-rad_bin_size) * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                rad_corrections.at(ang_index) * (n-1) * hold_r -
                                                                                n * (rad_bin_val.at(index - 1)+rad_bin_size)) ,2.));
                                            }
                                        }
                                    }

                                }
                                else {

                                    // Calculate the radius
                                    hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                    if (hold_r < rad_max) {

                                        // Find the bin where this radius lies in
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r);
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                        if (index < resolution){
                                            hold_cor = (hold_r - rad_bin_val.at(index - 1));
                                            hold_cos = hold_y / hold_r;
                                            // Add sin(theta)cos n(theta) to bin n
                                            for (int n = 0; n < norder; ++n) {
                                                rad_bin_shutter_on_first.at(n * resolution + index - 1) +=
                                                        hold_cor / rad_bin_size * pow(hold_cos, n);
                                                rad_bin_shutter_on_first.at(n * resolution + index) +=
                                                        (1. - hold_cor / rad_bin_size) * pow(hold_cos, n);

                                                rad_bin_shutter_on_err_first.at(n * resolution + index - 1) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * hold_cor+
                                                                                                      pow(hold_cos, 2.) * hold_r -
                                                                                                      n * hold_cor * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                (n-1) * hold_r -
                                                                                n * rad_bin_val.at(index - 1)) ,2.));

                                                rad_bin_shutter_on_err_first.at(n * resolution + index) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * (hold_cor-rad_bin_size) -
                                                                                                      pow(hold_cos, 2.) * hold_r -
                                                                                                      n * (hold_cor-rad_bin_size) * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                (n-1) * hold_r -
                                                                                n * (rad_bin_val.at(index - 1)+rad_bin_size)) ,2.));
                                            }
                                        }

                                        //Add to the VMI image
                                        auto lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                        hold_x);
                                        auto lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                        hold_y);
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_on.at(index) += 1;
                                    }
                                }
                            }
                            if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                event_time->at(j) < time_sub_upper && event_time->at(j) > time_sub_lower &&
                                peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {
                                injections_time_last += 1;

                                hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                if (circularise) {
                                    // Calculate the radius and radially correct it
                                    hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                    if (hold_r < rad_max) {

                                        //Find the bin where the theta lies in
                                        if (hold_x >= 0) { //Upper plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                            std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        } else { //Lower plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                            -std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        }

                                        //Add to the VMI image
                                        auto lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                        hold_x);
                                        auto lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                        hold_y);
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_on.at(index) -= 1;

                                        //Add to circularised VMI image
                                        lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                   hold_x * rad_corrections.at(ang_index));
                                        lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                   hold_y * rad_corrections.at(ang_index));
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_circ.at(index) -= 1;

                                        //Find index correction for integral
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r * rad_corrections.at(ang_index));
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                        if (index < resolution){
                                            hold_cor = (hold_r * rad_corrections.at(ang_index) - rad_bin_val.at(index - 1));
                                            hold_cos = hold_y / hold_r;
                                            // Add sin(theta)cos n(theta) to bin n
                                            for (int n = 0; n < norder; ++n) {
                                                rad_bin_shutter_on_last.at(n * resolution + index - 1) +=
                                                        hold_cor / rad_bin_size * pow(hold_cos, n);
                                                rad_bin_shutter_on_last.at(n * resolution + index) +=
                                                        (1. - hold_cor / rad_bin_size) * pow(hold_cos, n);

                                                rad_bin_shutter_on_err_last.at(n * resolution + index - 1) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * hold_cor+
                                                                                                      pow(hold_cos, 2.) * hold_r * rad_corrections.at(ang_index) -
                                                                                                      n * hold_cor * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                rad_corrections.at(ang_index) * (n-1) * hold_r -
                                                                                n * rad_bin_val.at(index - 1)) ,2.));

                                                rad_bin_shutter_on_err_last.at(n * resolution + index) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * (hold_cor-rad_bin_size) -
                                                                                                      pow(hold_cos, 2.) * hold_r * rad_corrections.at(ang_index) -
                                                                                                      n * (hold_cor-rad_bin_size) * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                rad_corrections.at(ang_index) * (n-1) * hold_r -
                                                                                n * (rad_bin_val.at(index - 1)+rad_bin_size)) ,2.));
                                            }
                                        }
                                    }

                                }
                                else {

                                    // Calculate the radius
                                    hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                    if (hold_r < rad_max) {

                                        // Find the bin where this radius lies in
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r);
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                        if (index < resolution){
                                            hold_cor = (hold_r - rad_bin_val.at(index - 1));
                                            hold_cos = hold_y / hold_r;
                                            // Add sin(theta)cos n(theta) to bin n
                                            for (int n = 0; n < norder; ++n) {
                                                rad_bin_shutter_on_last.at(n * resolution + index - 1) +=
                                                        hold_cor / rad_bin_size * pow(hold_cos, n);
                                                rad_bin_shutter_on_last.at(n * resolution + index) +=
                                                        (1. - hold_cor / rad_bin_size) * pow(hold_cos, n);

                                                rad_bin_shutter_on_err_last.at(n * resolution + index - 1) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * hold_cor+
                                                                                                      pow(hold_cos, 2.) * hold_r -
                                                                                                      n * hold_cor * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                (n-1) * hold_r -
                                                                                n * rad_bin_val.at(index - 1)) ,2.));

                                                rad_bin_shutter_on_err_last.at(n * resolution + index) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * (hold_cor-rad_bin_size) -
                                                                                                      pow(hold_cos, 2.) * hold_r -
                                                                                                      n * (hold_cor-rad_bin_size) * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                (n-1) * hold_r -
                                                                                n * (rad_bin_val.at(index - 1)+rad_bin_size)) ,2.));
                                            }
                                        }

                                        //Add to the VMI image
                                        auto lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                        hold_x);
                                        auto lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                        hold_y);
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_on.at(index) -= 1;
                                    }
                                }
                            }
                        }
                    } else {
                        injections_laser_off += 1.;
                        for (int j = 0; j < X->size(); j++) {
                            // Define time window
                            if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                event_time->at(j) < frame_upper && event_time->at(j) > frame_lower &&
                                peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {
                                injections_off_time_first += 1;

                                hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                if (circularise) {

                                    // Calculate the radius and radially correct it
                                    hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                    if (hold_r < rad_max*max_rad_corr) {
                                        //Find the bin where the theta lies in
                                        if (hold_x >= 0) { //Upper plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                            std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        } else { //Lower plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                            -std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        }

                                        //Add to the VMI image
                                        auto lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                        hold_x);
                                        auto lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                        hold_y);
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_off.at(index) += 1;

                                        //Add to circularised VMI image
                                        lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                   hold_x * rad_corrections.at(ang_index));
                                        lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                   hold_y * rad_corrections.at(ang_index));
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_circ.at(index) -= 1;

                                        //Find index correction for integral
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r * rad_corrections.at(ang_index));
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                        if (index < resolution){
                                            hold_cor = (hold_r * rad_corrections.at(ang_index) - rad_bin_val.at(index - 1));
                                            hold_cos = hold_y / hold_r;
                                            // Add sin(theta)cos n(theta) to bin n
                                            for (int n = 0; n < norder; ++n) {
                                                rad_bin_shutter_off_first.at(n * resolution + index - 1) +=
                                                        hold_cor / rad_bin_size * pow(hold_cos, n);
                                                rad_bin_shutter_off_first.at(n * resolution + index) +=
                                                        (1. - hold_cor / rad_bin_size) * pow(hold_cos, n);

                                                rad_bin_shutter_off_err_first.at(n * resolution + index - 1) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * hold_cor+
                                                                                                      pow(hold_cos, 2.) * hold_r * rad_corrections.at(ang_index) -
                                                                                                      n * hold_cor * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                rad_corrections.at(ang_index) * (n-1) * hold_r -
                                                                                n * rad_bin_val.at(index - 1)) ,2.));

                                                rad_bin_shutter_off_err_first.at(n * resolution + index) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * (hold_cor-rad_bin_size) -
                                                                                                      pow(hold_cos, 2.) * hold_r * rad_corrections.at(ang_index) -
                                                                                                      n * (hold_cor-rad_bin_size) * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                rad_corrections.at(ang_index) * (n-1) * hold_r -
                                                                                n * (rad_bin_val.at(index - 1)+rad_bin_size)) ,2.));
                                            }
                                        }
                                    }

                                }
                                else {
                                    // Calculate the radius
                                    hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                    if (hold_r < rad_max) {

                                        // Find the bin where this radius lies ini
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r);
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                        if (index < resolution){
                                            hold_cor = (hold_r - rad_bin_val.at(index - 1));
                                            hold_cos = hold_y / hold_r;
                                            // Add sin(theta)cos n(theta) to bin n
                                            for (int n = 0; n < norder; ++n) {
                                                rad_bin_shutter_off_first.at(n * resolution + index - 1) +=
                                                        hold_cor / rad_bin_size * pow(hold_cos, n);
                                                rad_bin_shutter_off_first.at(n * resolution + index) +=
                                                        (1. - hold_cor / rad_bin_size) * pow(hold_cos, n);

                                                rad_bin_shutter_off_err_first.at(n * resolution + index - 1) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * hold_cor+
                                                                                                      pow(hold_cos, 2.) * hold_r -
                                                                                                      n * hold_cor * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                (n-1) * hold_r -
                                                                                n * rad_bin_val.at(index - 1)) ,2.));

                                                rad_bin_shutter_off_err_first.at(n * resolution + index) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * (hold_cor-rad_bin_size) -
                                                                                                      pow(hold_cos, 2.) * hold_r -
                                                                                                      n * (hold_cor-rad_bin_size) * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                (n-1) * hold_r -
                                                                                n * (rad_bin_val.at(index - 1)+rad_bin_size)) ,2.));
                                            }
                                        }

                                        //Add to the VMI image
                                        auto lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                        hold_x);
                                        auto lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                        hold_y);
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_off.at(index) += 1;
                                    }
                                }
                            }
                            if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                event_time->at(j) < time_sub_upper && event_time->at(j) > time_sub_lower &&
                                peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {
                                injections_off_time_last += 1;

                                hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                if (circularise) {

                                    // Calculate the radius and radially correct it
                                    hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                    if (hold_r < rad_max*max_rad_corr) {
                                        //Find the bin where the theta lies in
                                        if (hold_x >= 0) { //Upper plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                            std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        } else { //Lower plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                            -std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        }

                                        //Add to the VMI image
                                        auto lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                        hold_x);
                                        auto lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                        hold_y);
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_off.at(index) -= 1;

                                        //Add to circularised VMI image
                                        lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                   hold_x * rad_corrections.at(ang_index));
                                        lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                   hold_y * rad_corrections.at(ang_index));
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_circ.at(index) += 1;

                                        //Find index correction for integral
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r * rad_corrections.at(ang_index));
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                        if (index < resolution){
                                            hold_cor = (hold_r * rad_corrections.at(ang_index) - rad_bin_val.at(index - 1));
                                            hold_cos = hold_y / hold_r;
                                            // Add sin(theta)cos n(theta) to bin n
                                            for (int n = 0; n < norder; ++n) {
                                                rad_bin_shutter_off_last.at(n * resolution + index - 1) +=
                                                        hold_cor / rad_bin_size * pow(hold_cos, n);
                                                rad_bin_shutter_off_last.at(n * resolution + index) +=
                                                        (1. - hold_cor / rad_bin_size) * pow(hold_cos, n);

                                                rad_bin_shutter_off_err_last.at(n * resolution + index - 1) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * hold_cor+
                                                                                                      pow(hold_cos, 2.) * hold_r * rad_corrections.at(ang_index) -
                                                                                                      n * hold_cor * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                rad_corrections.at(ang_index) * (n-1) * hold_r -
                                                                                n * rad_bin_val.at(index - 1)) ,2.));

                                                rad_bin_shutter_off_err_last.at(n * resolution + index) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * (hold_cor-rad_bin_size) -
                                                                                                      pow(hold_cos, 2.) * hold_r * rad_corrections.at(ang_index) -
                                                                                                      n * (hold_cor-rad_bin_size) * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                rad_corrections.at(ang_index) * (n-1) * hold_r -
                                                                                n * (rad_bin_val.at(index - 1)+rad_bin_size)) ,2.));
                                            }
                                        }
                                    }

                                }
                                else {
                                    // Calculate the radius
                                    hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                    if (hold_r < rad_max) {

                                        // Find the bin where this radius lies ini
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r);
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                        if (index < resolution){
                                            hold_cor = (hold_r - rad_bin_val.at(index - 1));
                                            hold_cos = hold_y / hold_r;
                                            // Add sin(theta)cos n(theta) to bin n
                                            for (int n = 0; n < norder; ++n) {
                                                rad_bin_shutter_off_last.at(n * resolution + index - 1) +=
                                                        hold_cor / rad_bin_size * pow(hold_cos, n);
                                                rad_bin_shutter_off_last.at(n * resolution + index) +=
                                                        (1. - hold_cor / rad_bin_size) * pow(hold_cos, n);

                                                rad_bin_shutter_off_err_last.at(n * resolution + index - 1) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * hold_cor+
                                                                                                      pow(hold_cos, 2.) * hold_r -
                                                                                                      n * hold_cor * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                (n-1) * hold_r -
                                                                                n * rad_bin_val.at(index - 1)) ,2.));

                                                rad_bin_shutter_off_err_last.at(n * resolution + index) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * (hold_cor-rad_bin_size) -
                                                                                                      pow(hold_cos, 2.) * hold_r -
                                                                                                      n * (hold_cor-rad_bin_size) * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                (n-1) * hold_r -
                                                                                n * (rad_bin_val.at(index - 1)+rad_bin_size)) ,2.));
                                            }
                                        }

                                        //Add to the VMI image
                                        auto lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                        hold_x);
                                        auto lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                        hold_y);
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_off.at(index) -= 1;
                                    }
                                }
                            }
                        }
                    }
                }

                // Removing background and saving the data
                for (int i = 0; i < resolution * norder; ++i) {
                    // Sum, r_on - inject_on / inject_off * r_off
                    rad_bin_tot.at(i) = (rad_bin_shutter_on_first.at(i) - injections_laser_on / injections_laser_off *
                                                                          rad_bin_shutter_off_first.at(i)) - (injections_time_first / injections_time_last *
                                                                                                              rad_bin_shutter_on_last.at(i) - injections_off_time_first / injections_off_time_last * injections_laser_on / injections_laser_off *
                                                                                                                                              rad_bin_shutter_off_last.at(i));

                    // Error in this context means the variance*
                    // Error propagation of the sum above, r_on_err + (inject_on / inject_off)^2 * r_off_err
                    // Assume no error in injection count
                    rad_bin_err.at(i) = rad_bin_shutter_on_err_first.at(i) + pow(injections_laser_on / injections_laser_off,2.) *
                                                                             rad_bin_shutter_off_err_first.at(i) +
                                        pow(injections_time_first / injections_time_last, 2.) *
                                        rad_bin_shutter_on_err_last.at(i) + pow(injections_off_time_first / injections_off_time_last * injections_laser_on / injections_laser_off,2.) *
                                                                            rad_bin_shutter_off_err_last.at(i);
                }
            }
        }
        else if (no_laser) {
            if (!time_subtraction) {
                //Reading the data
                for (int i = 0; i < entries; i++) {
                    X_branch->GetEntry(tree->LoadTree(i));
                    Y_branch->GetEntry(tree->LoadTree(i));
                    trig_time_branch->GetEntry(tree->LoadTree(i));
                    event_time_branch->GetEntry(tree->LoadTree(i));
                    peak_sum_branch->GetEntry(tree->LoadTree(i));
                    shutter_branch->GetEntry(tree->LoadTree(i));

                    // Bin the data to the resolution
                    if (!shutter) {
                        injections_laser_off += 1.;
                        for (int j = 0; j < X->size(); j++) {
                            // Define time window
                            if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                event_time->at(j) < frame_upper && event_time->at(j) > frame_lower &&
                                peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {

                                hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                if (circularise) {

                                    // Calculate the radius and radially correct it
                                    hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                    if (hold_r < rad_max*max_rad_corr) {
                                        //Find the bin where the theta lies in
                                        if (hold_x >= 0) { //Upper plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                            std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        } else { //Lower plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                            -std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        }

                                        //Add to the VMI image
                                        auto lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                        hold_x);
                                        auto lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                        hold_y);
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_off.at(index) += 1;

                                        //Add to circularised VMI image
                                        lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                   hold_x * rad_corrections.at(ang_index));
                                        lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                   hold_y * rad_corrections.at(ang_index));
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_circ.at(index) += 1;

                                        //Find index correction for integral
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r * rad_corrections.at(ang_index));
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                        if (index < resolution){
                                            hold_cor = (hold_r * rad_corrections.at(ang_index) - rad_bin_val.at(index - 1));
                                            hold_cos = hold_y / hold_r;
                                            // Add sin(theta)cos n(theta) to bin n
                                            for (int n = 0; n < norder; ++n) {
                                                rad_bin_shutter_off.at(n * resolution + index - 1) +=
                                                        hold_cor / rad_bin_size * pow(hold_cos, n);
                                                rad_bin_shutter_off.at(n * resolution + index) +=
                                                        (1. - hold_cor / rad_bin_size) * pow(hold_cos, n);

                                                rad_bin_shutter_off_err.at(n * resolution + index - 1) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * hold_cor+
                                                                                                      pow(hold_cos, 2.) * hold_r * rad_corrections.at(ang_index) -
                                                                                                      n * hold_cor * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                rad_corrections.at(ang_index) * (n-1) * hold_r -
                                                                                n * rad_bin_val.at(index - 1)) ,2.));

                                                rad_bin_shutter_off_err.at(n * resolution + index) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * (hold_cor-rad_bin_size) -
                                                                                                      pow(hold_cos, 2.) * hold_r * rad_corrections.at(ang_index) -
                                                                                                      n * (hold_cor-rad_bin_size) * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                rad_corrections.at(ang_index) * (n-1) * hold_r -
                                                                                n * (rad_bin_val.at(index - 1)+rad_bin_size)) ,2.));
                                            }
                                        }
                                    }

                                }
                                else {
                                    // Calculate the radius
                                    hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                    if (hold_r < rad_max) {

                                        // Find the bin where this radius lies ini
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r);
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                        if (index < resolution){
                                            hold_cor = (hold_r - rad_bin_val.at(index - 1));
                                            hold_cos = hold_y / hold_r;
                                            // Add sin(theta)cos n(theta) to bin n
                                            for (int n = 0; n < norder; ++n) {
                                                rad_bin_shutter_off.at(n * resolution + index - 1) +=
                                                        hold_cor / rad_bin_size * pow(hold_cos, n);
                                                rad_bin_shutter_off.at(n * resolution + index) +=
                                                        (1. - hold_cor / rad_bin_size) * pow(hold_cos, n);

                                                rad_bin_shutter_off_err.at(n * resolution + index - 1) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * hold_cor+
                                                                                                      pow(hold_cos, 2.) * hold_r -
                                                                                                      n * hold_cor * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                (n-1) * hold_r -
                                                                                n * rad_bin_val.at(index - 1)) ,2.));

                                                rad_bin_shutter_off_err.at(n * resolution + index) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * (hold_cor-rad_bin_size) -
                                                                                                      pow(hold_cos, 2.) * hold_r -
                                                                                                      n * (hold_cor-rad_bin_size) * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                (n-1) * hold_r -
                                                                                n * (rad_bin_val.at(index - 1)+rad_bin_size)) ,2.));
                                            }
                                        }

                                        //Add to the VMI image
                                        auto lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                        hold_x);
                                        auto lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                        hold_y);
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_off.at(index) += 1;
                                    }
                                }
                            }
                        }
                    }
                }

                // Removing background and saving the data
                for (int i = 0; i < resolution * norder; ++i) {
                    // Sum, r_on - inject_on / inject_off * r_off
                    rad_bin_tot.at(i) = rad_bin_shutter_off.at(i);

                    // Error in this context means the variance*
                    // Error propagation of the sum above, r_on_err + (inject_on / inject_off)^2 * r_off_err
                    // Assume no error in injection count
                    rad_bin_err.at(i) = rad_bin_shutter_off_err.at(i);
                }
            }
            else {
                //Reading the data
                for (int i = 0; i < entries; i++) {
                    X_branch->GetEntry(tree->LoadTree(i));
                    Y_branch->GetEntry(tree->LoadTree(i));
                    trig_time_branch->GetEntry(tree->LoadTree(i));
                    event_time_branch->GetEntry(tree->LoadTree(i));
                    peak_sum_branch->GetEntry(tree->LoadTree(i));
                    shutter_branch->GetEntry(tree->LoadTree(i));

                    // Bin the data to the resolution
                    if (!shutter) {
                        injections_laser_off += 1.;
                        for (int j = 0; j < X->size(); j++) {
                            // Define time window
                            if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                event_time->at(j) < frame_upper && event_time->at(j) > frame_lower &&
                                peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {

                                hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                if (circularise) {

                                    // Calculate the radius and radially correct it
                                    hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                    if (hold_r < rad_max*max_rad_corr) {
                                        //Find the bin where the theta lies in
                                        if (hold_x >= 0) { //Upper plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                            std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        } else { //Lower plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                            -std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        }

                                        //Add to the VMI image
                                        auto lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                        hold_x);
                                        auto lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                        hold_y);
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_off.at(index) += 1;

                                        //Add to circularised VMI image
                                        lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                   hold_x * rad_corrections.at(ang_index));
                                        lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                   hold_y * rad_corrections.at(ang_index));
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_circ.at(index) += 1;

                                        //Find index correction for integral
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r * rad_corrections.at(ang_index));
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                        if (index < resolution){
                                            hold_cor = (hold_r * rad_corrections.at(ang_index) - rad_bin_val.at(index - 1));
                                            hold_cos = hold_y / hold_r;
                                            // Add sin(theta)cos n(theta) to bin n
                                            for (int n = 0; n < norder; ++n) {
                                                rad_bin_shutter_off_first.at(n * resolution + index - 1) +=
                                                        hold_cor / rad_bin_size * pow(hold_cos, n);
                                                rad_bin_shutter_off_first.at(n * resolution + index) +=
                                                        (1. - hold_cor / rad_bin_size) * pow(hold_cos, n);

                                                rad_bin_shutter_off_err_first.at(n * resolution + index - 1) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * hold_cor+
                                                                                                      pow(hold_cos, 2.) * hold_r * rad_corrections.at(ang_index) -
                                                                                                      n * hold_cor * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                rad_corrections.at(ang_index) * (n-1) * hold_r -
                                                                                n * rad_bin_val.at(index - 1)) ,2.));

                                                rad_bin_shutter_off_err_first.at(n * resolution + index) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * (hold_cor-rad_bin_size) -
                                                                                                      pow(hold_cos, 2.) * hold_r * rad_corrections.at(ang_index) -
                                                                                                      n * (hold_cor-rad_bin_size) * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                rad_corrections.at(ang_index) * (n-1) * hold_r -
                                                                                n * (rad_bin_val.at(index - 1)+rad_bin_size)) ,2.));
                                            }
                                        }
                                    }

                                }
                                else {
                                    // Calculate the radius
                                    hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                    if (hold_r < rad_max) {

                                        // Find the bin where this radius lies ini
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r);
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                        if (index < resolution){
                                            hold_cor = (hold_r - rad_bin_val.at(index - 1));
                                            hold_cos = hold_y / hold_r;
                                            // Add sin(theta)cos n(theta) to bin n
                                            for (int n = 0; n < norder; ++n) {
                                                rad_bin_shutter_off_first.at(n * resolution + index - 1) +=
                                                        hold_cor / rad_bin_size * pow(hold_cos, n);
                                                rad_bin_shutter_off_first.at(n * resolution + index) +=
                                                        (1. - hold_cor / rad_bin_size) * pow(hold_cos, n);

                                                rad_bin_shutter_off_err_first.at(n * resolution + index - 1) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * hold_cor+
                                                                                                      pow(hold_cos, 2.) * hold_r -
                                                                                                      n * hold_cor * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                (n-1) * hold_r -
                                                                                n * rad_bin_val.at(index - 1)) ,2.));

                                                rad_bin_shutter_off_err_first.at(n * resolution + index) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * (hold_cor-rad_bin_size) -
                                                                                                      pow(hold_cos, 2.) * hold_r -
                                                                                                      n * (hold_cor-rad_bin_size) * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                (n-1) * hold_r -
                                                                                n * (rad_bin_val.at(index - 1)+rad_bin_size)) ,2.));
                                            }
                                        }

                                        //Add to the VMI image
                                        auto lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                        hold_x);
                                        auto lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                        hold_y);
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_off.at(index) += 1;
                                    }
                                }
                            }
                            if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                event_time->at(j) < time_sub_upper && event_time->at(j) > time_sub_lower &&
                                peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {

                                hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                if (circularise) {

                                    // Calculate the radius and radially correct it
                                    hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                    if (hold_r < rad_max*max_rad_corr) {
                                        //Find the bin where the theta lies in
                                        if (hold_x >= 0) { //Upper plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                            std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        } else { //Lower plane
                                            auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                            -std::acos(hold_y / hold_r));
                                            ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                        }

                                        //Add to the VMI image
                                        auto lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                        hold_x);
                                        auto lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                        hold_y);
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_off.at(index) -= 1;

                                        //Add to circularised VMI image
                                        lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                   hold_x * rad_corrections.at(ang_index));
                                        lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                   hold_y * rad_corrections.at(ang_index));
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_circ.at(index) -= 1;

                                        //Find index correction for integral
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r * rad_corrections.at(ang_index));
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                        if (index < resolution){
                                            hold_cor = (hold_r * rad_corrections.at(ang_index) - rad_bin_val.at(index - 1));
                                            hold_cos = hold_y / hold_r;
                                            // Add sin(theta)cos n(theta) to bin n
                                            for (int n = 0; n < norder; ++n) {
                                                rad_bin_shutter_off_last.at(n * resolution + index - 1) +=
                                                        hold_cor / rad_bin_size * pow(hold_cos, n);
                                                rad_bin_shutter_off_last.at(n * resolution + index) +=
                                                        (1. - hold_cor / rad_bin_size) * pow(hold_cos, n);

                                                rad_bin_shutter_off_err_last.at(n * resolution + index - 1) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * hold_cor+
                                                                                                      pow(hold_cos, 2.) * hold_r * rad_corrections.at(ang_index) -
                                                                                                      n * hold_cor * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                rad_corrections.at(ang_index) * (n-1) * hold_r -
                                                                                n * rad_bin_val.at(index - 1)) ,2.));

                                                rad_bin_shutter_off_err_last.at(n * resolution + index) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * (hold_cor-rad_bin_size) -
                                                                                                      pow(hold_cos, 2.) * hold_r * rad_corrections.at(ang_index) -
                                                                                                      n * (hold_cor-rad_bin_size) * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                rad_corrections.at(ang_index) * (n-1) * hold_r -
                                                                                n * (rad_bin_val.at(index - 1)+rad_bin_size)) ,2.));
                                            }
                                        }
                                    }

                                }
                                else {
                                    // Calculate the radius
                                    hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                    if (hold_r < rad_max) {

                                        // Find the bin where this radius lies ini
                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                        hold_r);
                                        index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                        if (index < resolution){
                                            hold_cor = (hold_r - rad_bin_val.at(index - 1));
                                            hold_cos = hold_y / hold_r;
                                            // Add sin(theta)cos n(theta) to bin n
                                            for (int n = 0; n < norder; ++n) {
                                                rad_bin_shutter_off_last.at(n * resolution + index - 1) +=
                                                        hold_cor / rad_bin_size * pow(hold_cos, n);
                                                rad_bin_shutter_off_last.at(n * resolution + index) +=
                                                        (1. - hold_cor / rad_bin_size) * pow(hold_cos, n);

                                                rad_bin_shutter_off_err_last.at(n * resolution + index - 1) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * hold_cor+
                                                                                                      pow(hold_cos, 2.) * hold_r -
                                                                                                      n * hold_cor * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                (n-1) * hold_r -
                                                                                n * rad_bin_val.at(index - 1)) ,2.));

                                                rad_bin_shutter_off_err_last.at(n * resolution + index) +=
                                                        center_error / pow(rad_bin_size,2.) * ( //Variance
                                                                pow(// d/dy term
                                                                        pow(hold_cos, n) * hold_r * ( n * (hold_cor-rad_bin_size) -
                                                                                                      pow(hold_cos, 2.) * hold_r -
                                                                                                      n * (hold_cor-rad_bin_size) * pow(hold_cos,2)) , 2.) +
                                                                pow(// d/dx term
                                                                        hold_x * pow(hold_cos, n) * pow(hold_r,2.) * (
                                                                                (n-1) * hold_r -
                                                                                n * (rad_bin_val.at(index - 1)+rad_bin_size)) ,2.));
                                            }
                                        }

                                        //Add to the VMI image
                                        auto lower_x = std::lower_bound(VMI_bin_x.begin(), VMI_bin_x.end(),
                                                                        hold_x);
                                        auto lower_y = std::lower_bound(VMI_bin_y.begin(), VMI_bin_y.end(),
                                                                        hold_y);
                                        index = (int) std::distance(VMI_bin_y.begin(), lower_y) * VMI_res +
                                                (int) std::distance(VMI_bin_x.begin(), lower_x);
                                        VMI_off.at(index) -= 1;
                                    }
                                }
                            }
                        }
                    }
                }

                // Removing background and saving the data
                for (int i = 0; i < resolution * norder; ++i) {
                    // Sum, r_on - inject_on / inject_off * r_off
                    rad_bin_tot.at(i) = rad_bin_shutter_off_first.at(i) - injections_time_first / injections_time_last * rad_bin_shutter_off_last.at(i);

                    // Error in this context means the variance*
                    // Error propagation of the sum above, r_on_err + (inject_on / inject_off)^2 * r_off_err
                    // Assume no error in injection count
                    rad_bin_err.at(i) = rad_bin_shutter_off_err_first.at(i)  +
                                        pow(injections_time_first / injections_time_last, 2.) * rad_bin_shutter_off_err_last.at(i);
                }
            }
        }

        /*--------            Folding image             --------*/
        if (circularise) {
            for (int col = 0; col < VMI_res; ++col) {
                for (int row = 0; row < VMI_res; ++row) {
                    //check for which side (col) is in
                    if (col >= center_col_row) { //First
                        VMI_fold.at(row * VMI_res + col) += (int)(VMI_circ.at(row * VMI_res + col)) +
                                                            (int)(VMI_circ.at(row * VMI_res - col + 2*center_col_row));
                    }
                    else {
                        VMI_fold.at(row * VMI_res + col) += (int)(VMI_circ.at(row * VMI_res + col)) +
                                                            (int)(VMI_circ.at(row * VMI_res + (VMI_res-1) - col));
                    }
                }
            }
        }
        else {
            if (normal) {
                for (int col = 0; col < VMI_res; ++col) {
                    for (int row = 0; row < VMI_res; ++row) {
                        //check for which side (col) is in
                        if (col >= center_col_row) { //First
                            VMI_fold.at(row * VMI_res + col) += (int) (VMI_on.at(row * VMI_res + col)) +
                                                                (int) (VMI_on.at(
                                                                        row * VMI_res - col + 2 * center_col_row));
                        } else {
                            VMI_fold.at(row * VMI_res + col) += (int) (VMI_on.at(row * VMI_res + col)) +
                                                                (int) (VMI_on.at(
                                                                        row * VMI_res + (VMI_res - 1) - col));
                        }
                    }
                }
            }
            else if (no_laser){
                for (int col = 0; col < VMI_res; ++col) {
                    for (int row = 0; row < VMI_res; ++row) {
                        //check for which side (col) is in
                        if (col >= center_col_row) { //First
                            VMI_fold.at(row * VMI_res + col) += (int) (VMI_off.at(row * VMI_res + col)) +
                                                                (int) (VMI_off.at(
                                                                        row * VMI_res - col + 2 * center_col_row));
                        } else {
                            VMI_fold.at(row * VMI_res + col) += (int) (VMI_off.at(row * VMI_res + col)) +
                                                                (int) (VMI_off.at(
                                                                        row * VMI_res + (VMI_res - 1) - col));
                        }
                    }
                }
            }
        }

        /*--------            Calculating Abel inversion             --------*/
        //Calculating inverse angular intensity distribution
        get_angular_intensity(norder, resolution, I_n,
                              inv_hankel, abel_inv_mat, rad_bin_val, rad_bin_tot);

        // Calculating all coefficients of legendre polynomials
        get_legendre_coef(norder, resolution, legendre_coefficients,
                          I_n, inv_legendre);

        // Taking out the intensity
        for (int i = 0; i < resolution; ++i) {
            I.at(i) = 4 * M_PI * pow(rad_bin_val.at(i), 2.) * legendre_coefficients.at(i);
            intensity.at(i) = I.at(i);
        }

        // Calculate the systematic error in intensity
        get_intensity_err(norder, resolution, I_err,
                          inv_hankel, inv_legendre, abel_inv_mat, rad_bin_err, rad_bin_val, rad_bin_tot);

        //Calculating inverse VMI image
        vmi_relative_res = biggest_diff / (VMI_res-1) * 2.;
        for (int x = 0; x < VMI_res; ++x) {
            for (int y = 0; y < VMI_res; ++y) {
                // Calculate the radius
                hold_r = pow((x - VMI_res / 2.) * (x - VMI_res / 2.) + (y - VMI_res / 2.) * (y - VMI_res / 2.),
                             1. / 2.) * vmi_relative_res;

                // Find the bin where this radius lies in
                auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                hold_r);
                index = (int) std::distance(rad_bin_val.begin(), lower_r);
                if (index >= resolution) {
                    index = resolution-1;
                }
                if (hold_r!=0. && hold_r < rad_max) {
                    for (int n = 0; n < norder; ++n) {
                        inv_vmi.at(y * VMI_res + x) += pow(rad_bin_val.at(index), 2.) * I_n.at(n * resolution + index) * pow((x - VMI_res / 2.) / hold_r*vmi_relative_res, n);
                    }
                }
            }
        }

        /*--------            Saving data in data structures             --------*/
        // Setting up new trees to save the binned data
        gif_file->cd();

        auto *gif_tree = new TTree((frame_name+std::to_string(fr)).c_str(),(frame_name+std::to_string(fr)).c_str());

        double raw_image_val;
        double raw_off_image_val;
        double inverse_image_val;
        double circular_image_val;
        double intensity_val;
        double intensity_err_val;
        double time_val;
        double settings_val;
        double rad_bin_tot_val;
        double rad_bin_val_val;
        double rad_bin_err_val;
        double I_val_err_l;
        double I_val_err_u;
        double I_iter_val;

        auto *raw_image_branch = gif_tree->Branch("raw_image", &raw_image_val);
        auto *raw_off_image_branch = gif_tree->Branch("raw_off_image", &raw_off_image_val);
        auto *inverse_image_branch = gif_tree->Branch("inv_image", &inverse_image_val);
        auto *circular_image_branch = gif_tree->Branch("circ_image", &circular_image_val);
        auto *intensity_branch = gif_tree->Branch("intensity", &intensity_val);
        auto *intensity_err_branch = gif_tree->Branch("intensity_err", &intensity_err_val);
        auto *time_branch = gif_tree->Branch("time", &time_val);
        auto *settings_branch = gif_tree->Branch("settings", &settings_val);
        auto *rad_tot_branch = gif_tree->Branch("radial_bin_TOTAL", &rad_bin_tot_val);
        auto *rad_val_branch = gif_tree->Branch("radial_bin_values", &rad_bin_val_val);
        auto *rad_err_branch = gif_tree->Branch("radial_bin_error", &rad_bin_err_val);

        // Saving the data
        for (auto &ent: VMI_on) {
            raw_image_val = ent;
            raw_image_branch->Fill();
        }
        for (auto &ent: VMI_off){
            raw_off_image_val = ent;
            raw_off_image_branch->Fill();
        }
        for (auto &ent: VMI_fold) {
            circular_image_val = ent;
            circular_image_branch->Fill();
        }
        for (auto &ent:inv_vmi) {
            if (ent < 0.){
                inverse_image_val = 0;
            } else {
                inverse_image_val = ent;
            }
            inverse_image_branch->Fill();
        }
        for (auto &ent: I) {
            if (ent > max_intensity) max_intensity = ent;
            if (ent < min_intensity) min_intensity = ent;
            intensity_val = ent;
            intensity_branch->Fill();
        }
        for (auto &ent: I_err) {
            intensity_err_val = ent;
            intensity_err_branch->Fill();
        }
        for (auto &ent: rad_bin_tot) {
            if (ent > max_radial_tot) max_radial_tot = ent;
            if (ent < min_radial_tot) min_radial_tot = ent;
            rad_bin_tot_val = ent;
            rad_tot_branch->Fill();
        }
        for (int i=0; i < resolution; i++) {
            rad_bin_val_val = rad_bin_val.at(i);
            rad_val_branch->Fill();
        }
        for (auto &ent: rad_bin_err) {
            rad_bin_err_val = ent;
            rad_err_branch->Fill();
        }

        time_val = frame_lower; time_branch->Fill();
        time_val = frame_upper; time_branch->Fill();

        settings_val = resolution; settings_branch->Fill();
        settings_val = rad_bin_val.at(resolution-1); settings_branch->Fill();
        settings_val = VMI_res; settings_branch->Fill();

        if (rad_bin_val.at(resolution-1) > max_radius) max_radius = rad_bin_val.at(resolution-1);

        gif_tree->Write(nullptr, TObject::kOverwrite);

        /*--------            Bootstrap statistical errors             --------*/
        if (bootstrap){
            data_file->cd();

            ref_radial.clear();
            ref_radial.resize(resolution);
            for (int i = 0; i < resolution; i++){
                ref_radial.at(i) = rad_bin_val.at(i);
            }

            if (!time_subtraction) {
                rad_bin_shutter_on.clear();
                rad_bin_shutter_off.clear();
                rad_bin_shutter_on.resize(resolution * norder, 0);
                rad_bin_shutter_off.resize(resolution * norder, 0);

                rad_bin_shutter_on_err.clear();
                rad_bin_shutter_off_err.clear();
                rad_bin_shutter_on_err.resize(resolution * norder, 0);
                rad_bin_shutter_off_err.resize(resolution * norder, 0);
            }
            else {
                if (normal) {
                    rad_bin_shutter_on_first.clear();
                    rad_bin_shutter_on_last.clear();
                    rad_bin_shutter_on_first.resize(resolution * norder, 0);
                    rad_bin_shutter_on_last.resize(resolution * norder, 0);
                    rad_bin_shutter_on_err_first.clear();
                    rad_bin_shutter_on_err_last.clear();
                    rad_bin_shutter_on_err_first.resize(resolution * norder, 0);
                    rad_bin_shutter_on_err_last.resize(resolution * norder, 0);
                }
                rad_bin_shutter_off_first.clear();
                rad_bin_shutter_off_last.clear();
                rad_bin_shutter_off_first.resize(resolution * norder, 0);
                rad_bin_shutter_off_last.resize(resolution * norder, 0);
                rad_bin_shutter_off_err_first.clear();
                rad_bin_shutter_off_err_last.clear();
                rad_bin_shutter_off_err_first.resize(resolution * norder, 0);
                rad_bin_shutter_off_err_last.resize(resolution * norder, 0);
            }
            rad_bin_tot.clear();
            rad_bin_val.clear();
            rad_bin_err.clear();
            I_n.clear();
            legendre_coefficients.clear();
            rad_bin_val.resize(resolution+1);
            rad_bin_tot.resize(resolution * norder, 0);
            rad_bin_err.resize(resolution * norder);

            injections_laser_on = 0;
            injections_laser_off = 0;
            injections_time_last = 0;
            injections_time_first = 0;
            injections_off_time_last = 0;
            injections_off_time_first = 0;

            I_iter.clear();
            z.clear(); //Bias correction
            iteration.clear();
            I_lower.clear();
            I_upper.clear();

            I_iter.resize(nr_samples*resolution);
            z.resize(resolution); //Bias correction
            iteration.resize(nr_samples);
            I_lower.resize(resolution,0);
            I_upper.resize(resolution,0);

            std::cout << "Bootstrap sampling frame " << fr+1 << "\n" << std::endl;
            // Multithread or not
            if (nr_threads > 1) {
                for (int sample = 0; sample < nr_samples; ++sample) {
                    start_time = std::chrono::steady_clock::now();
                    /*--------            Reset dependencies             --------*/
                    if (!time_subtraction) {
                        rad_bin_shutter_on.clear();
                        rad_bin_shutter_off.clear();
                        rad_bin_shutter_on.resize(resolution * norder, 0);
                        rad_bin_shutter_off.resize(resolution * norder, 0);
                    }
                    else {
                        if (normal) {
                            rad_bin_shutter_on_first.clear();
                            rad_bin_shutter_on_last.clear();
                            rad_bin_shutter_on_first.resize(resolution * norder, 0);
                            rad_bin_shutter_on_last.resize(resolution * norder, 0);
                        }
                        rad_bin_shutter_off_first.clear();
                        rad_bin_shutter_off_last.clear();
                        rad_bin_shutter_off_first.resize(resolution * norder, 0);
                        rad_bin_shutter_off_last.resize(resolution * norder, 0);
                    }
                    rad_bin_tot.clear();
                    rad_bin_val.clear();
                    I_n.clear();
                    legendre_coefficients.clear();
                    rad_bin_val.resize(resolution + 1);
                    rad_bin_tot.resize(resolution * norder, 0);

                    injections_laser_on = 0;
                    injections_laser_off = 0;
                    injections_time_last = 0;
                    injections_time_first = 0;
                    injections_off_time_last = 0;
                    injections_off_time_first = 0;

                    /*--------            Reading from file             --------*/
                    if (!circularise) {
                        // ReMake a vector of the radial values
                        rad_bin_size = rad_max / resolution;
                        for (int i = 0; i < resolution+1; i++) {
                            rad_bin_val.at(i) = i * rad_bin_size;
                            if (i!=resolution)
                                new_radial[i] = i * rad_bin_size;
                        }
                    }
                    else {
                        // ReMake a vector of the radial values
                        rad_bin_size = rad_max * max_rad_corr / resolution;
                        for (int i = 0; i < resolution+1; i++) {
                            rad_bin_val.at(i) = i * rad_bin_size;
                            if (i!=resolution)
                                new_radial[i] = i * rad_bin_size;
                        }
                    }

                    if (normal) {
                        if (!time_subtraction){
                            // Creating thread safe objects
                            ROOT::TThreadedObject<std::vector<double>> th_rad_bin_on_obj(resolution * norder, 0.);
                            ROOT::TThreadedObject<std::vector<double>> th_rad_bin_off_obj(resolution * norder, 0.);
                            ROOT::TThreadedObject<double> th_injection_on_obj(0.);
                            ROOT::TThreadedObject<double> th_injection_off_obj(0.);

                            // The multi-threaded reader of the TTree
                            auto thread_data_read = [&](TTreeReader &reader) {
                                auto th_rad_bin_on = th_rad_bin_on_obj.Get();
                                auto th_rad_bin_off = th_rad_bin_off_obj.Get();
                                auto th_injection_on = th_injection_on_obj.Get();
                                auto th_injection_off = th_injection_off_obj.Get();

                                TTreeReaderValue<int> injection_reader(reader, "injection");
                                TTreeReaderValue<bool> shutter_reader(reader, "shutter");
                                TTreeReaderValue<std::vector<double>> trig_time_reader(reader, "trig_time");
                                TTreeReaderValue<std::vector<double>> event_time_reader(reader, "event_time");
                                TTreeReaderValue<std::vector<double>> peak_sum_reader(reader, "event_time");
                                TTreeReaderValue<std::vector<double>> X_reader(reader, "X");
                                TTreeReaderValue<std::vector<double>> Y_reader(reader, "Y");

                                while (reader.Next()) {
                                    auto X = *X_reader;
                                    auto Y = *Y_reader;
                                    auto trig_time = *trig_time_reader;
                                    auto event_time = *event_time_reader;
                                    auto peak_sum = *peak_sum_reader;
                                    auto shutter = *shutter_reader;

                                    double th_hold_r;
                                    double th_hold_x;
                                    double th_hold_y;
                                    double th_hold_cos;
                                    double th_hold_cor;
                                    double th_hold_weight;
                                    int th_index;
                                    int th_ang_index;

                                    // Bin the data to the resolution
                                    if (shutter) {
                                        th_hold_weight = weight(gen);
                                        *th_injection_on += th_hold_weight;
                                        for (int j = 0; j < X.size(); j++) {
                                            // Define time window
                                            if (trig_time.at(j) < trig_time_upper && trig_time.at(j) > trig_time_lower &&
                                                event_time.at(j) < frame_upper && event_time.at(j) > frame_lower &&
                                                peak_sum.at(j) < peak_sum_upper && peak_sum.at(j) > peak_sum_lower) {

                                                th_hold_x = ((X.at(j) - x_center) + (Y.at(j) - y_center)) / M_SQRT2;
                                                th_hold_y = (-(X.at(j) - x_center) + (Y.at(j) - y_center)) / M_SQRT2;

                                                if (circularise) {
                                                    // Calculate the radius and radially correct it
                                                    th_hold_r = pow(th_hold_x*th_hold_x + th_hold_y*th_hold_y, 1. / 2.);

                                                    if (th_hold_r < rad_max) {

                                                        th_hold_cos = th_hold_y /th_hold_r;

                                                        // Find the bin where this radius lies in
                                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                        th_hold_r);
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                                        //Find the bin where the theta lies in
                                                        if (th_hold_x >= 0) { //Upper plane
                                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                                            angular_bin_val.end(),
                                                                                            std::acos(th_hold_cos));
                                                            th_ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                        } else { //Lower plane
                                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                                            angular_bin_val.end(),
                                                                                            -std::acos(
                                                                                                    th_hold_cos));
                                                            th_ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                        }

                                                        //Find index correction for integral
                                                        lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                   th_hold_r * rad_corrections.at(th_ang_index));
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                                        if (th_index < resolution) {

                                                            th_hold_cor = (th_hold_r *rad_corrections.at(th_ang_index) -
                                                                           rad_bin_val.at(th_index - 1))/
                                                                          rad_bin_size;

                                                            // Add sin(theta)cos n(theta) to bin n
                                                            for (int n = 0; n < norder; ++n) {
                                                                th_rad_bin_on->at(n * resolution + th_index - 1) +=
                                                                        th_hold_weight *
                                                                        th_hold_cor *
                                                                        pow(th_hold_cos, n);
                                                                th_rad_bin_on->at(n * resolution + th_index) +=
                                                                        th_hold_weight *
                                                                        (1. - th_hold_cor) *
                                                                        pow(th_hold_cos, n);;
                                                            }
                                                        }
                                                    }
                                                } else {

                                                    // Calculate the radius
                                                    th_hold_r = pow(th_hold_x*th_hold_x + th_hold_y*th_hold_y, 1. / 2.);

                                                    if (th_hold_r < rad_max) {

                                                        // Find the bin where this radius lies in
                                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                        th_hold_r);
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                                        if (th_index < resolution) {

                                                            th_hold_cos = th_hold_y /th_hold_r;
                                                            th_hold_cor = (th_hold_r - rad_bin_val.at(th_index -1)) / rad_bin_size;

                                                            // Add sin(theta)cos n(theta) to bin n
                                                            for (int n = 0; n < norder; ++n) {
                                                                th_rad_bin_on->at(n * resolution + th_index - 1) +=
                                                                        th_hold_weight *
                                                                        th_hold_cor *
                                                                        pow(th_hold_cos, n);
                                                                th_rad_bin_on->at(n * resolution + th_index) +=
                                                                        th_hold_weight *
                                                                        (1. - th_hold_cor) *
                                                                        pow(th_hold_cos, n);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    } else {
                                        th_hold_weight = weight(gen);
                                        *th_injection_off += th_hold_weight;
                                        for (int j = 0; j < X.size(); j++) {
                                            // Define time window
                                            if (trig_time.at(j) < trig_time_upper && trig_time.at(j) > trig_time_lower &&
                                                event_time.at(j) < frame_upper && event_time.at(j) > frame_lower &&
                                                peak_sum.at(j) < peak_sum_upper && peak_sum.at(j) > peak_sum_lower) {

                                                th_hold_x = ((X.at(j) - x_center) + (Y.at(j) - y_center)) / M_SQRT2;
                                                th_hold_y = (-(X.at(j) - x_center) + (Y.at(j) - y_center)) / M_SQRT2;

                                                if (circularise) {

                                                    // Calculate the radius and radially correct it
                                                    th_hold_r = pow(th_hold_x*th_hold_x + th_hold_y*th_hold_y, 1. / 2.);

                                                    if (th_hold_r < rad_max) {

                                                        th_hold_cos = th_hold_y /th_hold_r;

                                                        // Find the bin where this radius lies in
                                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                        th_hold_r);
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                                        //Find the bin where the theta lies in
                                                        if (th_hold_x >= 0) { //Upper plane
                                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                                            angular_bin_val.end(),
                                                                                            std::acos(
                                                                                                    th_hold_cos));
                                                            th_ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                        } else { //Lower plane
                                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                                            angular_bin_val.end(),
                                                                                            -std::acos(th_hold_cos));
                                                            th_ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                        }

                                                        //Find th_index correction for integral
                                                        lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                   th_hold_r * rad_corrections.at(th_ang_index));
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                                        if (th_index < resolution) {

                                                            th_hold_cor = (th_hold_r * rad_corrections.at(th_ang_index) -
                                                                           rad_bin_val.at(th_index -1))/
                                                                          rad_bin_size;

                                                            // Add sin(theta)cos n(theta) to bin n
                                                            for (int n = 0; n < norder; ++n) {
                                                                th_rad_bin_off->at(n * resolution + th_index) +=
                                                                        th_hold_weight *
                                                                        th_hold_cor *
                                                                        pow(th_hold_cos, n);
                                                                th_rad_bin_off->at(n * resolution + th_index - 1) +=
                                                                        th_hold_weight *
                                                                        (1. - th_hold_cor) *
                                                                        pow(th_hold_cos, n);
                                                            }
                                                        }
                                                    }

                                                } else {
                                                    // Calculate the radius
                                                    th_hold_r = pow(th_hold_x*th_hold_x + th_hold_y*th_hold_y, 1. / 2.);

                                                    if (th_hold_r < rad_max) {

                                                        // Find the bin where this radius lies ini
                                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                        th_hold_r);
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                                        if (th_index < resolution) {

                                                            th_hold_cos = th_hold_y /th_hold_r;
                                                            th_hold_cor = (th_hold_r - rad_bin_val.at(th_index -1))/
                                                                          rad_bin_size;

                                                            // Add sin(theta)cos n(theta) to bin n
                                                            for (int n = 0; n < norder; ++n) {
                                                                th_rad_bin_off->at(n * resolution + th_index) +=
                                                                        th_hold_weight *
                                                                        th_hold_cor *
                                                                        pow(th_hold_cos, n);
                                                                th_rad_bin_off->at(n * resolution + th_index - 1) +=
                                                                        th_hold_weight *
                                                                        (1. - th_hold_cor) *
                                                                        pow(th_hold_cos, n);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            };

                            tp.Process(thread_data_read);

                            // Merge threads
                            for (int i = 0; i < nr_threads; ++i) {
                                // Merge vectors
                                auto st = *th_rad_bin_on_obj.GetAtSlot(i);
                                for (int j = 0; j < resolution * norder; ++j) {
                                    rad_bin_shutter_on.at(j) += st.at(j);
                                }
                                st = *th_rad_bin_off_obj.GetAtSlot(i);
                                for (int j = 0; j < resolution * norder; ++j) {
                                    rad_bin_shutter_off.at(j) += st.at(j);
                                }

                                //Merge injection counts
                                auto ent = th_injection_on_obj.GetAtSlot(i);
                                injections_laser_on += *ent;
                                ent = th_injection_off_obj.GetAtSlot(i);
                                injections_laser_off += *ent;
                            }

                            // Removing background and saving the data
                            for (int i = 0; i < resolution * norder; ++i) {
                                rad_bin_tot.at(i) = rad_bin_shutter_on.at(i) - injections_laser_on / injections_laser_off *
                                                                               rad_bin_shutter_off.at(i);
                            }
                        }
                        else {
                            // Creating thread safe objects
                            ROOT::TThreadedObject<std::vector<double>> th_rad_bin_on_first_obj(resolution * norder, 0.);
                            ROOT::TThreadedObject<std::vector<double>> th_rad_bin_off_first_obj(resolution * norder, 0.);
                            ROOT::TThreadedObject<std::vector<double>> th_rad_bin_on_last_obj(resolution * norder, 0.);
                            ROOT::TThreadedObject<std::vector<double>> th_rad_bin_off_last_obj(resolution * norder, 0.);
                            ROOT::TThreadedObject<double> th_injection_on_obj(0.);
                            ROOT::TThreadedObject<double> th_injection_off_obj(0.);
                            ROOT::TThreadedObject<double> injections_time_last_obj(0.);
                            ROOT::TThreadedObject<double> injections_time_first_obj(0.);
                            ROOT::TThreadedObject<double> injections_off_time_last_obj(0.);
                            ROOT::TThreadedObject<double> injections_off_time_first_obj(0.);

                            // The multi-threaded reader of the TTree
                            auto thread_data_read = [&](TTreeReader &reader) {
                                auto th_rad_bin_on_first = th_rad_bin_on_first_obj.Get();
                                auto th_rad_bin_off_first = th_rad_bin_off_first_obj.Get();
                                auto th_rad_bin_on_last = th_rad_bin_on_last_obj.Get();
                                auto th_rad_bin_off_last = th_rad_bin_off_last_obj.Get();
                                auto th_injection_on = th_injection_on_obj.Get();
                                auto th_injection_off = th_injection_off_obj.Get();
                                auto th_injections_time_last = injections_time_last_obj.Get();
                                auto th_injections_time_first = injections_time_first_obj.Get();
                                auto th_injections_off_time_last = injections_off_time_last_obj.Get();
                                auto th_injections_off_time_first = injections_off_time_first_obj.Get();

                                TTreeReaderValue<int> injection_reader(reader, "injection");
                                TTreeReaderValue<bool> shutter_reader(reader, "shutter");
                                TTreeReaderValue<std::vector<double>> trig_time_reader(reader, "trig_time");
                                TTreeReaderValue<std::vector<double>> event_time_reader(reader, "event_time");
                                TTreeReaderValue<std::vector<double>> peak_sum_reader(reader, "event_time");
                                TTreeReaderValue<std::vector<double>> X_reader(reader, "X");
                                TTreeReaderValue<std::vector<double>> Y_reader(reader, "Y");

                                while (reader.Next()) {
                                    auto X = *X_reader;
                                    auto Y = *Y_reader;
                                    auto trig_time = *trig_time_reader;
                                    auto event_time = *event_time_reader;
                                    auto peak_sum = *peak_sum_reader;
                                    auto shutter = *shutter_reader;

                                    double th_hold_r;
                                    double th_hold_x;
                                    double th_hold_y;
                                    double th_hold_cos;
                                    double th_hold_cor;
                                    double th_hold_weight;
                                    int th_index;
                                    int th_ang_index;

                                    // Bin the data to the resolution
                                    if (shutter) {
                                        th_hold_weight = weight(gen);
                                        *th_injection_on += th_hold_weight;
                                        for (int j = 0; j < X.size(); j++) {
                                            // Define time window
                                            if (trig_time.at(j) < trig_time_upper && trig_time.at(j) > trig_time_lower &&
                                                event_time.at(j) < frame_upper && event_time.at(j) > frame_lower &&
                                                peak_sum.at(j) < peak_sum_upper && peak_sum.at(j) > peak_sum_lower) {
                                                *th_injections_time_first += th_hold_weight;

                                                th_hold_x = ((X.at(j) - x_center) + (Y.at(j) - y_center)) / M_SQRT2;
                                                th_hold_y = (-(X.at(j) - x_center) + (Y.at(j) - y_center)) / M_SQRT2;

                                                if (circularise) {
                                                    // Calculate the radius and radially correct it
                                                    th_hold_r = pow(th_hold_x*th_hold_x + th_hold_y*th_hold_y, 1. / 2.);

                                                    if (th_hold_r < rad_max) {

                                                        th_hold_cos = th_hold_y /th_hold_r;

                                                        // Find the bin where this radius lies in
                                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                        th_hold_r);
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                                        //Find the bin where the theta lies in
                                                        if (th_hold_x >= 0) { //Upper plane
                                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                                            angular_bin_val.end(),
                                                                                            std::acos(th_hold_cos));
                                                            th_ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                        } else { //Lower plane
                                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                                            angular_bin_val.end(),
                                                                                            -std::acos(
                                                                                                    th_hold_cos));
                                                            th_ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                        }

                                                        //Find index correction for integral
                                                        lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                   th_hold_r * rad_corrections.at(th_ang_index));
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                                        if (th_index < resolution) {

                                                            th_hold_cor = (th_hold_r *rad_corrections.at(th_ang_index) -
                                                                           rad_bin_val.at(th_index - 1))/
                                                                          rad_bin_size;

                                                            // Add sin(theta)cos n(theta) to bin n
                                                            for (int n = 0; n < norder; ++n) {
                                                                th_rad_bin_on_first->at(n * resolution + th_index - 1) +=
                                                                        th_hold_weight *
                                                                        th_hold_cor *
                                                                        pow(th_hold_cos, n);
                                                                th_rad_bin_on_first->at(n * resolution + th_index) +=
                                                                        th_hold_weight *
                                                                        (1. - th_hold_cor) *
                                                                        pow(th_hold_cos, n);;
                                                            }
                                                        }
                                                    }
                                                } else {

                                                    // Calculate the radius
                                                    th_hold_r = pow(th_hold_x*th_hold_x + th_hold_y*th_hold_y, 1. / 2.);

                                                    if (th_hold_r < rad_max) {

                                                        // Find the bin where this radius lies in
                                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                        th_hold_r);
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                                        if (th_index < resolution) {

                                                            th_hold_cos = th_hold_y /th_hold_r;
                                                            th_hold_cor = (th_hold_r - rad_bin_val.at(th_index -1)) / rad_bin_size;

                                                            // Add sin(theta)cos n(theta) to bin n
                                                            for (int n = 0; n < norder; ++n) {
                                                                th_rad_bin_on_first->at(n * resolution + th_index - 1) +=
                                                                        th_hold_weight *
                                                                        th_hold_cor *
                                                                        pow(th_hold_cos, n);
                                                                th_rad_bin_on_first->at(n * resolution + th_index) +=
                                                                        th_hold_weight *
                                                                        (1. - th_hold_cor) *
                                                                        pow(th_hold_cos, n);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            if (trig_time.at(j) < trig_time_upper && trig_time.at(j) > trig_time_lower &&
                                                event_time.at(j) < time_sub_upper && event_time.at(j) > time_sub_lower &&
                                                peak_sum.at(j) < peak_sum_upper && peak_sum.at(j) > peak_sum_lower) {
                                                *th_injections_time_last += th_hold_weight;

                                                th_hold_x = ((X.at(j) - x_center) + (Y.at(j) - y_center)) / M_SQRT2;
                                                th_hold_y = (-(X.at(j) - x_center) + (Y.at(j) - y_center)) / M_SQRT2;

                                                if (circularise) {
                                                    // Calculate the radius and radially correct it
                                                    th_hold_r = pow(th_hold_x*th_hold_x + th_hold_y*th_hold_y, 1. / 2.);

                                                    if (th_hold_r < rad_max) {

                                                        th_hold_cos = th_hold_y /th_hold_r;

                                                        // Find the bin where this radius lies in
                                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                        th_hold_r);
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                                        //Find the bin where the theta lies in
                                                        if (th_hold_x >= 0) { //Upper plane
                                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                                            angular_bin_val.end(),
                                                                                            std::acos(th_hold_cos));
                                                            th_ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                        } else { //Lower plane
                                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                                            angular_bin_val.end(),
                                                                                            -std::acos(
                                                                                                    th_hold_cos));
                                                            th_ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                        }

                                                        //Find index correction for integral
                                                        lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                   th_hold_r * rad_corrections.at(th_ang_index));
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                                        if (th_index < resolution) {

                                                            th_hold_cor = (th_hold_r *rad_corrections.at(th_ang_index) -
                                                                           rad_bin_val.at(th_index - 1))/
                                                                          rad_bin_size;

                                                            // Add sin(theta)cos n(theta) to bin n
                                                            for (int n = 0; n < norder; ++n) {
                                                                th_rad_bin_on_last->at(n * resolution + th_index - 1) +=
                                                                        th_hold_weight *
                                                                        th_hold_cor *
                                                                        pow(th_hold_cos, n);
                                                                th_rad_bin_on_last->at(n * resolution + th_index) +=
                                                                        th_hold_weight *
                                                                        (1. - th_hold_cor) *
                                                                        pow(th_hold_cos, n);;
                                                            }
                                                        }
                                                    }
                                                } else {

                                                    // Calculate the radius
                                                    th_hold_r = pow(th_hold_x*th_hold_x + th_hold_y*th_hold_y, 1. / 2.);

                                                    if (th_hold_r < rad_max) {

                                                        // Find the bin where this radius lies in
                                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                        th_hold_r);
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                                        if (th_index < resolution) {

                                                            th_hold_cos = th_hold_y /th_hold_r;
                                                            th_hold_cor = (th_hold_r - rad_bin_val.at(th_index -1)) / rad_bin_size;

                                                            // Add sin(theta)cos n(theta) to bin n
                                                            for (int n = 0; n < norder; ++n) {
                                                                th_rad_bin_on_last->at(n * resolution + th_index - 1) +=
                                                                        th_hold_weight *
                                                                        th_hold_cor *
                                                                        pow(th_hold_cos, n);
                                                                th_rad_bin_on_last->at(n * resolution + th_index) +=
                                                                        th_hold_weight *
                                                                        (1. - th_hold_cor) *
                                                                        pow(th_hold_cos, n);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    } else {
                                        th_hold_weight = weight(gen);
                                        *th_injection_off += th_hold_weight;
                                        for (int j = 0; j < X.size(); j++) {
                                            // Define time window
                                            if (trig_time.at(j) < trig_time_upper && trig_time.at(j) > trig_time_lower &&
                                                event_time.at(j) < frame_upper && event_time.at(j) > frame_lower &&
                                                peak_sum.at(j) < peak_sum_upper && peak_sum.at(j) > peak_sum_lower) {
                                                *th_injections_off_time_first += th_hold_weight;

                                                th_hold_x = ((X.at(j) - x_center) + (Y.at(j) - y_center)) / M_SQRT2;
                                                th_hold_y = (-(X.at(j) - x_center) + (Y.at(j) - y_center)) / M_SQRT2;

                                                if (circularise) {

                                                    // Calculate the radius and radially correct it
                                                    th_hold_r = pow(th_hold_x*th_hold_x + th_hold_y*th_hold_y, 1. / 2.);

                                                    if (th_hold_r < rad_max) {

                                                        th_hold_cos = th_hold_y /th_hold_r;

                                                        // Find the bin where this radius lies in
                                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                        th_hold_r);
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                                        //Find the bin where the theta lies in
                                                        if (th_hold_x >= 0) { //Upper plane
                                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                                            angular_bin_val.end(),
                                                                                            std::acos(
                                                                                                    th_hold_cos));
                                                            th_ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                        } else { //Lower plane
                                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                                            angular_bin_val.end(),
                                                                                            -std::acos(th_hold_cos));
                                                            th_ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                        }

                                                        //Find th_index correction for integral
                                                        lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                   th_hold_r * rad_corrections.at(th_ang_index));
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                                        if (th_index < resolution) {

                                                            th_hold_cor = (th_hold_r * rad_corrections.at(th_ang_index) -
                                                                           rad_bin_val.at(th_index -1))/
                                                                          rad_bin_size;

                                                            // Add sin(theta)cos n(theta) to bin n
                                                            for (int n = 0; n < norder; ++n) {
                                                                th_rad_bin_off_first->at(n * resolution + th_index) +=
                                                                        th_hold_weight *
                                                                        th_hold_cor *
                                                                        pow(th_hold_cos, n);
                                                                th_rad_bin_off_first->at(n * resolution + th_index - 1) +=
                                                                        th_hold_weight *
                                                                        (1. - th_hold_cor) *
                                                                        pow(th_hold_cos, n);
                                                            }
                                                        }
                                                    }

                                                } else {
                                                    // Calculate the radius
                                                    th_hold_r = pow(th_hold_x*th_hold_x + th_hold_y*th_hold_y, 1. / 2.);

                                                    if (th_hold_r < rad_max) {

                                                        // Find the bin where this radius lies ini
                                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                        th_hold_r);
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                                        if (th_index < resolution) {

                                                            th_hold_cos = th_hold_y /th_hold_r;
                                                            th_hold_cor = (th_hold_r - rad_bin_val.at(th_index -1))/
                                                                          rad_bin_size;

                                                            // Add sin(theta)cos n(theta) to bin n
                                                            for (int n = 0; n < norder; ++n) {
                                                                th_rad_bin_off_first->at(n * resolution + th_index) +=
                                                                        th_hold_weight *
                                                                        th_hold_cor *
                                                                        pow(th_hold_cos, n);
                                                                th_rad_bin_off_first->at(n * resolution + th_index - 1) +=
                                                                        th_hold_weight *
                                                                        (1. - th_hold_cor) *
                                                                        pow(th_hold_cos, n);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            if (trig_time.at(j) < trig_time_upper && trig_time.at(j) > trig_time_lower &&
                                                event_time.at(j) < time_sub_upper && event_time.at(j) > time_sub_lower &&
                                                peak_sum.at(j) < peak_sum_upper && peak_sum.at(j) > peak_sum_lower) {
                                                *th_injections_off_time_last += th_hold_weight;

                                                th_hold_x = ((X.at(j) - x_center) + (Y.at(j) - y_center)) / M_SQRT2;
                                                th_hold_y = (-(X.at(j) - x_center) + (Y.at(j) - y_center)) / M_SQRT2;

                                                if (circularise) {

                                                    // Calculate the radius and radially correct it
                                                    th_hold_r = pow(th_hold_x*th_hold_x + th_hold_y*th_hold_y, 1. / 2.);

                                                    if (th_hold_r < rad_max) {

                                                        th_hold_cos = th_hold_y /th_hold_r;

                                                        // Find the bin where this radius lies in
                                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                        th_hold_r);
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                                        //Find the bin where the theta lies in
                                                        if (th_hold_x >= 0) { //Upper plane
                                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                                            angular_bin_val.end(),
                                                                                            std::acos(
                                                                                                    th_hold_cos));
                                                            th_ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                        } else { //Lower plane
                                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                                            angular_bin_val.end(),
                                                                                            -std::acos(th_hold_cos));
                                                            th_ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                        }

                                                        //Find th_index correction for integral
                                                        lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                   th_hold_r * rad_corrections.at(th_ang_index));
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                                        if (th_index < resolution) {

                                                            th_hold_cor = (th_hold_r * rad_corrections.at(th_ang_index) -
                                                                           rad_bin_val.at(th_index -1))/
                                                                          rad_bin_size;

                                                            // Add sin(theta)cos n(theta) to bin n
                                                            for (int n = 0; n < norder; ++n) {
                                                                th_rad_bin_off_last->at(n * resolution + th_index) +=
                                                                        th_hold_weight *
                                                                        th_hold_cor *
                                                                        pow(th_hold_cos, n);
                                                                th_rad_bin_off_last->at(n * resolution + th_index - 1) +=
                                                                        th_hold_weight *
                                                                        (1. - th_hold_cor) *
                                                                        pow(th_hold_cos, n);
                                                            }
                                                        }
                                                    }

                                                } else {
                                                    // Calculate the radius
                                                    th_hold_r = pow(th_hold_x*th_hold_x + th_hold_y*th_hold_y, 1. / 2.);

                                                    if (th_hold_r < rad_max) {

                                                        // Find the bin where this radius lies ini
                                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                        th_hold_r);
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                                        if (th_index < resolution) {

                                                            th_hold_cos = th_hold_y /th_hold_r;
                                                            th_hold_cor = (th_hold_r - rad_bin_val.at(th_index -1))/
                                                                          rad_bin_size;

                                                            // Add sin(theta)cos n(theta) to bin n
                                                            for (int n = 0; n < norder; ++n) {
                                                                th_rad_bin_off_last->at(n * resolution + th_index) +=
                                                                        th_hold_weight *
                                                                        th_hold_cor *
                                                                        pow(th_hold_cos, n);
                                                                th_rad_bin_off_last->at(n * resolution + th_index - 1) +=
                                                                        th_hold_weight *
                                                                        (1. - th_hold_cor) *
                                                                        pow(th_hold_cos, n);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            };

                            tp.Process(thread_data_read);

                            // Merge threads
                            for (int i = 0; i < nr_threads; ++i) {
                                // Merge vectors
                                auto st = *th_rad_bin_on_first_obj.GetAtSlot(i);
                                for (int j = 0; j < resolution * norder; ++j) {
                                    rad_bin_shutter_on_first.at(j) += st.at(j);
                                }
                                st = *th_rad_bin_off_first_obj.GetAtSlot(i);
                                for (int j = 0; j < resolution * norder; ++j) {
                                    rad_bin_shutter_off_first.at(j) += st.at(j);
                                }
                                st = *th_rad_bin_on_last_obj.GetAtSlot(i);
                                for (int j = 0; j < resolution * norder; ++j) {
                                    rad_bin_shutter_on_last.at(j) += st.at(j);
                                }
                                st = *th_rad_bin_off_last_obj.GetAtSlot(i);
                                for (int j = 0; j < resolution * norder; ++j) {
                                    rad_bin_shutter_off_last.at(j) += st.at(j);
                                }

                                //Merge injection counts
                                auto ent = th_injection_on_obj.GetAtSlot(i);
                                injections_laser_on += *ent;
                                ent = th_injection_off_obj.GetAtSlot(i);
                                injections_laser_off += *ent;
                                ent = injections_time_last_obj.GetAtSlot(i);
                                injections_time_last += *ent;
                                ent = injections_time_first_obj.GetAtSlot(i);
                                injections_time_first += *ent;
                                ent = injections_off_time_last_obj.GetAtSlot(i);
                                injections_off_time_last += *ent;
                                ent = injections_off_time_first_obj.GetAtSlot(i);
                                injections_off_time_first += *ent;
                            }

                            // Removing background and saving the data
                            for (int i = 0; i < resolution * norder; ++i) {
                                // Sum, r_on - inject_on / inject_off * r_off
                                rad_bin_tot.at(i) = (rad_bin_shutter_on_first.at(i) - injections_laser_on / injections_laser_off *
                                                                                      rad_bin_shutter_off_first.at(i)) - (injections_time_first / injections_time_last * rad_bin_shutter_on_last.at(i) -
                                                                                                                          injections_off_time_first / injections_off_time_last * injections_laser_on / injections_laser_off *
                                                                                                                          rad_bin_shutter_off_last.at(i));
                            }
                        }
                    }
                    else if (no_laser) {
                        if (!time_subtraction){
                            // Creating thread safe objects
                            ROOT::TThreadedObject<std::vector<double>> th_rad_bin_off_obj(resolution * norder, 0.);

                            // The multi-threaded reader of the TTree
                            auto thread_data_read = [&](TTreeReader &reader) {
                                auto th_rad_bin_off = th_rad_bin_off_obj.Get();

                                TTreeReaderValue<int> injection_reader(reader, "injection");
                                TTreeReaderValue<bool> shutter_reader(reader, "shutter");
                                TTreeReaderValue<std::vector<double>> trig_time_reader(reader, "trig_time");
                                TTreeReaderValue<std::vector<double>> event_time_reader(reader, "event_time");
                                TTreeReaderValue<std::vector<double>> peak_sum_reader(reader, "event_time");
                                TTreeReaderValue<std::vector<double>> X_reader(reader, "X");
                                TTreeReaderValue<std::vector<double>> Y_reader(reader, "Y");

                                while (reader.Next()) {
                                    auto X = *X_reader;
                                    auto Y = *Y_reader;
                                    auto trig_time = *trig_time_reader;
                                    auto event_time = *event_time_reader;
                                    auto peak_sum = *peak_sum_reader;
                                    auto shutter = *shutter_reader;

                                    double th_hold_r;
                                    double th_hold_x;
                                    double th_hold_y;
                                    double th_hold_cos;
                                    double th_hold_cor;
                                    double th_hold_weight;
                                    int th_index;
                                    int th_ang_index;

                                    // Bin the data to the resolution
                                    if (!shutter) {
                                        th_hold_weight = weight(gen);
                                        for (int j = 0; j < X.size(); j++) {
                                            // Define time window
                                            if (trig_time.at(j) < trig_time_upper && trig_time.at(j) > trig_time_lower &&
                                                event_time.at(j) < frame_upper && event_time.at(j) > frame_lower &&
                                                peak_sum.at(j) < peak_sum_upper && peak_sum.at(j) > peak_sum_lower) {

                                                th_hold_x = ((X.at(j) - x_center) + (Y.at(j) - y_center)) / M_SQRT2;
                                                th_hold_y = (-(X.at(j) - x_center) + (Y.at(j) - y_center)) / M_SQRT2;

                                                if (circularise) {
                                                    // Calculate the radius and radially correct it
                                                    th_hold_r = pow(th_hold_x*th_hold_x + th_hold_y*th_hold_y, 1. / 2.);

                                                    if (th_hold_r < rad_max) {

                                                        th_hold_cos = th_hold_y /th_hold_r;

                                                        // Find the bin where this radius lies in
                                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                        th_hold_r);
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                                        //Find the bin where the theta lies in
                                                        if (th_hold_x >= 0) { //Upper plane
                                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                                            angular_bin_val.end(),
                                                                                            std::acos(th_hold_cos));
                                                            th_ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                        } else { //Lower plane
                                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                                            angular_bin_val.end(),
                                                                                            -std::acos(
                                                                                                    th_hold_cos));
                                                            th_ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                        }

                                                        //Find index correction for integral
                                                        lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                   th_hold_r * rad_corrections.at(th_ang_index));
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                                        if (th_index < resolution) {

                                                            th_hold_cor = (th_hold_r *rad_corrections.at(th_ang_index) -
                                                                           rad_bin_val.at(th_index - 1))/
                                                                          rad_bin_size;

                                                            // Add sin(theta)cos n(theta) to bin n
                                                            for (int n = 0; n < norder; ++n) {
                                                                th_rad_bin_off->at(n * resolution + th_index - 1) +=
                                                                        th_hold_weight *
                                                                        th_hold_cor *
                                                                        pow(th_hold_cos, n);
                                                                th_rad_bin_off->at(n * resolution + th_index) +=
                                                                        th_hold_weight *
                                                                        (1. - th_hold_cor) *
                                                                        pow(th_hold_cos, n);;
                                                            }
                                                        }
                                                    }
                                                } else {

                                                    // Calculate the radius
                                                    th_hold_r = pow(th_hold_x*th_hold_x + th_hold_y*th_hold_y, 1. / 2.);

                                                    if (th_hold_r < rad_max) {

                                                        // Find the bin where this radius lies in
                                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                        th_hold_r);
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                                        if (th_index < resolution) {

                                                            th_hold_cos = th_hold_y /th_hold_r;
                                                            th_hold_cor = (th_hold_r - rad_bin_val.at(th_index -1)) / rad_bin_size;

                                                            // Add sin(theta)cos n(theta) to bin n
                                                            for (int n = 0; n < norder; ++n) {
                                                                th_rad_bin_off->at(n * resolution + th_index - 1) +=
                                                                        th_hold_weight *
                                                                        th_hold_cor *
                                                                        pow(th_hold_cos, n);
                                                                th_rad_bin_off->at(n * resolution + th_index) +=
                                                                        th_hold_weight *
                                                                        (1. - th_hold_cor) *
                                                                        pow(th_hold_cos, n);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            };

                            tp.Process(thread_data_read);

                            // Merge threads
                            for (int i = 0; i < nr_threads; ++i) {
                                // Merge vectors
                                auto st = *th_rad_bin_off_obj.GetAtSlot(i);
                                for (int j = 0; j < resolution * norder; ++j) {
                                    rad_bin_shutter_off.at(j) += st.at(j);
                                }
                            }

                            // Removing background and saving the data
                            for (int i = 0; i < resolution * norder; ++i) {
                                rad_bin_tot.at(i) = rad_bin_shutter_off.at(i);
                            }
                        }
                        else {
                            // Creating thread safe objects
                            ROOT::TThreadedObject<std::vector<double>> th_rad_bin_off_first_obj(resolution * norder, 0.);
                            ROOT::TThreadedObject<std::vector<double>> th_rad_bin_off_last_obj(resolution * norder, 0.);
                            ROOT::TThreadedObject<double> injections_time_last_obj(0.);
                            ROOT::TThreadedObject<double> injections_time_first_obj(0.);

                            // The multi-threaded reader of the TTree
                            auto thread_data_read = [&](TTreeReader &reader) {
                                auto th_rad_bin_off_first = th_rad_bin_off_first_obj.Get();
                                auto th_rad_bin_off_last = th_rad_bin_off_last_obj.Get();
                                auto th_injections_time_last = injections_time_last_obj.Get();
                                auto th_injections_time_first = injections_time_first_obj.Get();

                                TTreeReaderValue<int> injection_reader(reader, "injection");
                                TTreeReaderValue<bool> shutter_reader(reader, "shutter");
                                TTreeReaderValue<std::vector<double>> trig_time_reader(reader, "trig_time");
                                TTreeReaderValue<std::vector<double>> event_time_reader(reader, "event_time");
                                TTreeReaderValue<std::vector<double>> peak_sum_reader(reader, "event_time");
                                TTreeReaderValue<std::vector<double>> X_reader(reader, "X");
                                TTreeReaderValue<std::vector<double>> Y_reader(reader, "Y");

                                while (reader.Next()) {
                                    auto X = *X_reader;
                                    auto Y = *Y_reader;
                                    auto trig_time = *trig_time_reader;
                                    auto event_time = *event_time_reader;
                                    auto peak_sum = *peak_sum_reader;
                                    auto shutter = *shutter_reader;

                                    double th_hold_r;
                                    double th_hold_x;
                                    double th_hold_y;
                                    double th_hold_cos;
                                    double th_hold_cor;
                                    double th_hold_weight;
                                    int th_index;
                                    int th_ang_index;

                                    // Bin the data to the resolution
                                    if (!shutter) {
                                        th_hold_weight = weight(gen);
                                        for (int j = 0; j < X.size(); j++) {
                                            // Define time window
                                            if (trig_time.at(j) < trig_time_upper && trig_time.at(j) > trig_time_lower &&
                                                event_time.at(j) < frame_upper && event_time.at(j) > frame_lower &&
                                                peak_sum.at(j) < peak_sum_upper && peak_sum.at(j) > peak_sum_lower) {
                                                *th_injections_time_first += th_hold_weight;

                                                th_hold_x = ((X.at(j) - x_center) + (Y.at(j) - y_center)) / M_SQRT2;
                                                th_hold_y = (-(X.at(j) - x_center) + (Y.at(j) - y_center)) / M_SQRT2;

                                                if (circularise) {
                                                    // Calculate the radius and radially correct it
                                                    th_hold_r = pow(th_hold_x*th_hold_x + th_hold_y*th_hold_y, 1. / 2.);

                                                    if (th_hold_r < rad_max) {

                                                        th_hold_cos = th_hold_y /th_hold_r;

                                                        // Find the bin where this radius lies in
                                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                        th_hold_r);
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                                        //Find the bin where the theta lies in
                                                        if (th_hold_x >= 0) { //Upper plane
                                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                                            angular_bin_val.end(),
                                                                                            std::acos(th_hold_cos));
                                                            th_ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                        } else { //Lower plane
                                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                                            angular_bin_val.end(),
                                                                                            -std::acos(
                                                                                                    th_hold_cos));
                                                            th_ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                        }

                                                        //Find index correction for integral
                                                        lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                   th_hold_r * rad_corrections.at(th_ang_index));
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                                        if (th_index < resolution) {

                                                            th_hold_cor = (th_hold_r *rad_corrections.at(th_ang_index) -
                                                                           rad_bin_val.at(th_index - 1))/
                                                                          rad_bin_size;

                                                            // Add sin(theta)cos n(theta) to bin n
                                                            for (int n = 0; n < norder; ++n) {
                                                                th_rad_bin_off_first->at(n * resolution + th_index - 1) +=
                                                                        th_hold_weight *
                                                                        th_hold_cor *
                                                                        pow(th_hold_cos, n);
                                                                th_rad_bin_off_first->at(n * resolution + th_index) +=
                                                                        th_hold_weight *
                                                                        (1. - th_hold_cor) *
                                                                        pow(th_hold_cos, n);;
                                                            }
                                                        }
                                                    }
                                                } else {

                                                    // Calculate the radius
                                                    th_hold_r = pow(th_hold_x*th_hold_x + th_hold_y*th_hold_y, 1. / 2.);

                                                    if (th_hold_r < rad_max) {

                                                        // Find the bin where this radius lies in
                                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                        th_hold_r);
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                                        if (th_index < resolution) {

                                                            th_hold_cos = th_hold_y /th_hold_r;
                                                            th_hold_cor = (th_hold_r - rad_bin_val.at(th_index -1)) / rad_bin_size;

                                                            // Add sin(theta)cos n(theta) to bin n
                                                            for (int n = 0; n < norder; ++n) {
                                                                th_rad_bin_off_first->at(n * resolution + th_index - 1) +=
                                                                        th_hold_weight *
                                                                        th_hold_cor *
                                                                        pow(th_hold_cos, n);
                                                                th_rad_bin_off_first->at(n * resolution + th_index) +=
                                                                        th_hold_weight *
                                                                        (1. - th_hold_cor) *
                                                                        pow(th_hold_cos, n);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            if (trig_time.at(j) < trig_time_upper && trig_time.at(j) > trig_time_lower &&
                                                event_time.at(j) < time_sub_upper && event_time.at(j) > time_sub_lower &&
                                                peak_sum.at(j) < peak_sum_upper && peak_sum.at(j) > peak_sum_lower) {
                                                *th_injections_time_last += th_hold_weight;

                                                th_hold_x = ((X.at(j) - x_center) + (Y.at(j) - y_center)) / M_SQRT2;
                                                th_hold_y = (-(X.at(j) - x_center) + (Y.at(j) - y_center)) / M_SQRT2;

                                                if (circularise) {
                                                    // Calculate the radius and radially correct it
                                                    th_hold_r = pow(th_hold_x*th_hold_x + th_hold_y*th_hold_y, 1. / 2.);

                                                    if (th_hold_r < rad_max) {

                                                        th_hold_cos = th_hold_y /th_hold_r;

                                                        // Find the bin where this radius lies in
                                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                        th_hold_r);
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                                        //Find the bin where the theta lies in
                                                        if (th_hold_x >= 0) { //Upper plane
                                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                                            angular_bin_val.end(),
                                                                                            std::acos(th_hold_cos));
                                                            th_ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                        } else { //Lower plane
                                                            auto lower_t = std::lower_bound(angular_bin_val.begin(),
                                                                                            angular_bin_val.end(),
                                                                                            -std::acos(
                                                                                                    th_hold_cos));
                                                            th_ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                        }

                                                        //Find index correction for integral
                                                        lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                   th_hold_r * rad_corrections.at(th_ang_index));
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                                        if (th_index < resolution) {

                                                            th_hold_cor = (th_hold_r *rad_corrections.at(th_ang_index) -
                                                                           rad_bin_val.at(th_index - 1))/
                                                                          rad_bin_size;

                                                            // Add sin(theta)cos n(theta) to bin n
                                                            for (int n = 0; n < norder; ++n) {
                                                                th_rad_bin_off_last->at(n * resolution + th_index - 1) +=
                                                                        th_hold_weight *
                                                                        th_hold_cor *
                                                                        pow(th_hold_cos, n);
                                                                th_rad_bin_off_last->at(n * resolution + th_index) +=
                                                                        th_hold_weight *
                                                                        (1. - th_hold_cor) *
                                                                        pow(th_hold_cos, n);;
                                                            }
                                                        }
                                                    }
                                                } else {

                                                    // Calculate the radius
                                                    th_hold_r = pow(th_hold_x*th_hold_x + th_hold_y*th_hold_y, 1. / 2.);

                                                    if (th_hold_r < rad_max) {

                                                        // Find the bin where this radius lies in
                                                        auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                        th_hold_r);
                                                        th_index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                                        if (th_index < resolution) {

                                                            th_hold_cos = th_hold_y /th_hold_r;
                                                            th_hold_cor = (th_hold_r - rad_bin_val.at(th_index -1)) / rad_bin_size;

                                                            // Add sin(theta)cos n(theta) to bin n
                                                            for (int n = 0; n < norder; ++n) {
                                                                th_rad_bin_off_last->at(n * resolution + th_index - 1) +=
                                                                        th_hold_weight *
                                                                        th_hold_cor *
                                                                        pow(th_hold_cos, n);
                                                                th_rad_bin_off_last->at(n * resolution + th_index) +=
                                                                        th_hold_weight *
                                                                        (1. - th_hold_cor) *
                                                                        pow(th_hold_cos, n);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            };

                            tp.Process(thread_data_read);

                            // Merge threads
                            for (int i = 0; i < nr_threads; ++i) {
                                // Merge vectors
                                auto st = *th_rad_bin_off_first_obj.GetAtSlot(i);
                                for (int j = 0; j < resolution * norder; ++j) {
                                    rad_bin_shutter_off_first.at(j) += st.at(j);
                                }
                                st = *th_rad_bin_off_last_obj.GetAtSlot(i);
                                for (int j = 0; j < resolution * norder; ++j) {
                                    rad_bin_shutter_off_last.at(j) += st.at(j);
                                }

                                //Merge injection counts
                                auto ent = injections_time_last_obj.GetAtSlot(i);
                                injections_time_last += *ent;
                                ent = injections_time_first_obj.GetAtSlot(i);
                                injections_time_first += *ent;
                            }

                            // Removing background and saving the data
                            for (int i = 0; i < resolution * norder; ++i) {
                                // Sum, r_on - inject_on / inject_off * r_off
                                rad_bin_tot.at(i) = rad_bin_shutter_off_first.at(i) - injections_time_first / injections_time_last * rad_bin_shutter_off_last.at(i);
                            }
                        }
                    }

                    /*--------            Calculating Abel inversion             --------*/
                    //Calculating inverse angular intensity distribution
                    get_angular_intensity(norder, resolution, I_n,
                                          inv_hankel, abel_inv_mat, rad_bin_val, rad_bin_tot);

                    // Calculating all coefficients of legendre polynomials
                    get_legendre_coef(norder, resolution, legendre_coefficients,
                                      I_n, inv_legendre);

                    // Taking out the intensity
                    for (int i = 0; i < resolution; ++i) {
                        new_intensity[i] = 4 * M_PI * pow(rad_bin_val.at(i), 2.) * legendre_coefficients.at(i);
                    }

                    //Making a spline to evaluate the intensity at the old radial values
                    spline = new TSpline3("test", new_radial, new_intensity, resolution, "b1e1");

                    //Save the sample statistics.
                    for (int i = 0; i < resolution; ++i) {
                        I_iter.at(sample * resolution + i) = spline->Eval(ref_radial.at(i));
                    }

                    spline->Clear();

                    end_time = std::chrono::steady_clock::now();
                    elapsed_seconds = end_time - start_time;
                    running_sum += elapsed_seconds;

                    std::cout << "Sample " << sample + 1 << " done in " << elapsed_seconds.count() << " s\n" <<
                              "Expected done in " << (nr_samples - sample) * running_sum.count() / 60. / (sample+1) << " min."
                              << std::endl;
                }
            }
            else {
                for (int sample = 0; sample < nr_samples; ++sample) {
                    start_time = std::chrono::steady_clock::now();
                    /*--------            Reset dependencies             --------*/

                    if (!time_subtraction) {
                        rad_bin_shutter_on.clear();
                        rad_bin_shutter_off.clear();
                        rad_bin_shutter_on.resize(resolution * norder, 0);
                        rad_bin_shutter_off.resize(resolution * norder, 0);
                    }
                    else {
                        if (normal) {
                            rad_bin_shutter_on_first.clear();
                            rad_bin_shutter_on_last.clear();
                            rad_bin_shutter_on_first.resize(resolution * norder, 0);
                            rad_bin_shutter_on_last.resize(resolution * norder, 0);
                        }
                        rad_bin_shutter_off_first.clear();
                        rad_bin_shutter_off_last.clear();
                        rad_bin_shutter_off_first.resize(resolution * norder, 0);
                        rad_bin_shutter_off_last.resize(resolution * norder, 0);
                    }
                    rad_bin_tot.clear();
                    rad_bin_val.clear();
                    I_n.clear();
                    legendre_coefficients.clear();
                    rad_bin_val.resize(resolution+1);
                    rad_bin_tot.resize(resolution * norder, 0);

                    injections_laser_on = 0;
                    injections_laser_off = 0;
                    injections_time_last = 0;
                    injections_time_first = 0;
                    injections_off_time_last = 0;
                    injections_off_time_first = 0;

                    /*--------            Reading from file             --------*/
                    if (!circularise) {
                        // ReMake a vector of the radial values
                        rad_bin_size = rad_max / resolution;
                        for (int i = 0; i < resolution+1; i++) {
                            rad_bin_val.at(i) = i * rad_bin_size;
                            if (i!=resolution)
                                new_radial[i] = i * rad_bin_size;
                        }
                    }
                    else {
                        // ReMake a vector of the radial values
                        rad_bin_size = rad_max * max_rad_corr / resolution;
                        for (int i = 0; i < resolution+1; i++) {
                            rad_bin_val.at(i) = i * rad_bin_size;
                            if (i!=resolution)
                                new_radial[i] = i * rad_bin_size;
                        }
                    }

                    // Reading data depending on mode
                    if (normal) {
                        if (!time_subtraction) {
                            //Reading the data
                            for (int i = 0; i < entries; i++) {
                                X_branch->GetEntry(tree->LoadTree(i));
                                Y_branch->GetEntry(tree->LoadTree(i));
                                trig_time_branch->GetEntry(tree->LoadTree(i));
                                event_time_branch->GetEntry(tree->LoadTree(i));
                                peak_sum_branch->GetEntry(tree->LoadTree(i));
                                shutter_branch->GetEntry(tree->LoadTree(i));

                                // Bin the data to the resolution
                                if (shutter) {
                                    hold_weight = weight(gen);
                                    injections_laser_on += hold_weight;
                                    for (int j = 0; j < X->size(); j++) {
                                        // Define time window
                                        if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                            event_time->at(j) < frame_upper && event_time->at(j) > frame_lower &&
                                            peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {

                                            hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                            hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                            if (circularise) {
                                                // Calculate the radius and radially correct it
                                                hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                                if (hold_r < rad_max) {

                                                    //Find the bin where the theta lies in
                                                    if (hold_x >= 0) { //Upper plane
                                                        auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                                        std::acos(hold_y / hold_r));
                                                        ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                    } else { //Lower plane
                                                        auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                                        -std::acos(hold_y / hold_r));
                                                        ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                    }

                                                    //Find index correction for integral
                                                    auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                    hold_r * rad_corrections.at(ang_index));
                                                    index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                                    if (index < resolution){
                                                        hold_cor = (hold_r * rad_corrections.at(ang_index) -
                                                                    rad_bin_val.at(index - 1)) / rad_bin_size;
                                                        hold_cos = hold_y / hold_r;
                                                        // Add sin(theta)cos n(theta) to bin n
                                                        for (int n = 0; n < norder; ++n) {
                                                            rad_bin_shutter_on.at(n * resolution + index - 1) +=
                                                                    hold_weight *
                                                                    hold_cor *
                                                                    pow(hold_cos, n);
                                                            rad_bin_shutter_on.at(n * resolution + index) +=
                                                                    hold_weight *
                                                                    (1. - hold_cor) *
                                                                    pow(hold_cos, n);
                                                        }
                                                    }
                                                }

                                            }
                                            else {

                                                // Calculate the radius
                                                hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                                if (hold_r < rad_max) {

                                                    // Find the bin where this radius lies in
                                                    auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                    hold_r);
                                                    index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                                    if (index < resolution){
                                                        hold_cor = (hold_r - rad_bin_val.at(index - 1))/
                                                                   rad_bin_size;
                                                        hold_cos = hold_y / hold_r;
                                                        // Add sin(theta)cos n(theta) to bin n
                                                        for (int n = 0; n < norder; ++n) {
                                                            rad_bin_shutter_on.at(n * resolution + index - 1) +=
                                                                    hold_weight *
                                                                    hold_cor *
                                                                    pow(hold_cos, n);
                                                            rad_bin_shutter_on.at(n * resolution + index) +=
                                                                    hold_weight *
                                                                    (1. - hold_cor)  *
                                                                    pow(hold_cos, n);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                } else {
                                    hold_weight = weight(gen);
                                    injections_laser_off += hold_weight;
                                    for (int j = 0; j < X->size(); j++) {
                                        // Define time window
                                        if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                            event_time->at(j) < frame_upper && event_time->at(j) > frame_lower &&
                                            peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {

                                            hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                            hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                            if (circularise) {

                                                // Calculate the radius and radially correct it
                                                hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                                if (hold_r < rad_max*max_rad_corr) {
                                                    //Find the bin where the theta lies in
                                                    if (hold_x >= 0) { //Upper plane
                                                        auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                                        std::acos(hold_y / hold_r));
                                                        ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                    } else { //Lower plane
                                                        auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                                        -std::acos(hold_y / hold_r));
                                                        ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                    }

                                                    //Find index correction for integral
                                                    auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                    hold_r * rad_corrections.at(ang_index));
                                                    index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                                    if (index < resolution){
                                                        hold_cor = (hold_r*rad_corrections.at(ang_index) - rad_bin_val.at(index - 1)) /
                                                                   rad_bin_size;
                                                        hold_cos = hold_y / hold_r;
                                                        // Add sin(theta)cos n(theta) to bin n
                                                        for (int n = 0; n < norder; ++n) {
                                                            rad_bin_shutter_off.at(n * resolution + index) +=
                                                                    hold_weight *
                                                                    hold_cor *
                                                                    pow(hold_cos, n);
                                                            rad_bin_shutter_off.at(n * resolution + index - 1) +=
                                                                    hold_weight *
                                                                    (1. - hold_cor) *
                                                                    pow(hold_cos, n);
                                                        }
                                                    }
                                                }

                                            }
                                            else {
                                                // Calculate the radius
                                                hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                                if (hold_r < rad_max) {

                                                    // Find the bin where this radius lies ini
                                                    auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                    hold_r);
                                                    index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                                    if (index < resolution){
                                                        hold_cor = (hold_r - rad_bin_val.at(index - 1)) /
                                                                   rad_bin_size;
                                                        hold_cos = hold_y / hold_r;
                                                        // Add sin(theta)cos n(theta) to bin n
                                                        for (int n = 0; n < norder; ++n) {
                                                            rad_bin_shutter_off.at(n * resolution + index) +=
                                                                    hold_weight *
                                                                    hold_cor *
                                                                    pow(hold_cos, n);
                                                            rad_bin_shutter_off.at(n * resolution + index - 1) +=
                                                                    hold_weight *
                                                                    (1. - hold_cor) *
                                                                    pow(hold_cos, n);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                            // Removing background and saving the data
                            for (int i = 0; i < resolution * norder; ++i) {
                                // Sum, r_on - inject_on / inject_off * r_off
                                rad_bin_tot.at(i) = rad_bin_shutter_on.at(i) - injections_laser_on / injections_laser_off *
                                                                               rad_bin_shutter_off.at(i);
                            }
                        }
                        else {
                            //Reading the data
                            for (int i = 0; i < entries; i++) {
                                X_branch->GetEntry(tree->LoadTree(i));
                                Y_branch->GetEntry(tree->LoadTree(i));
                                trig_time_branch->GetEntry(tree->LoadTree(i));
                                event_time_branch->GetEntry(tree->LoadTree(i));
                                peak_sum_branch->GetEntry(tree->LoadTree(i));
                                shutter_branch->GetEntry(tree->LoadTree(i));

                                // Bin the data to the resolution
                                if (shutter) {
                                    hold_weight = weight(gen);
                                    injections_laser_on += hold_weight;
                                    for (int j = 0; j < X->size(); j++) {
                                        // Define time window
                                        if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                            event_time->at(j) < frame_upper && event_time->at(j) > frame_lower &&
                                            peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {
                                            injections_time_first += hold_weight;

                                            hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                            hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                            if (circularise) {
                                                // Calculate the radius and radially correct it
                                                hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                                if (hold_r < rad_max) {

                                                    //Find the bin where the theta lies in
                                                    if (hold_x >= 0) { //Upper plane
                                                        auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                                        std::acos(hold_y / hold_r));
                                                        ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                    } else { //Lower plane
                                                        auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                                        -std::acos(hold_y / hold_r));
                                                        ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                    }

                                                    //Find index correction for integral
                                                    auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                    hold_r * rad_corrections.at(ang_index));
                                                    index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                                    if (index < resolution){
                                                        hold_cor = (hold_r * rad_corrections.at(ang_index) -
                                                                    rad_bin_val.at(index - 1)) / rad_bin_size;
                                                        hold_cos = hold_y / hold_r;
                                                        // Add sin(theta)cos n(theta) to bin n
                                                        for (int n = 0; n < norder; ++n) {
                                                            rad_bin_shutter_on_first.at(n * resolution + index - 1) +=
                                                                    hold_weight *
                                                                    hold_cor *
                                                                    pow(hold_cos, n);
                                                            rad_bin_shutter_on_first.at(n * resolution + index) +=
                                                                    hold_weight *
                                                                    (1. - hold_cor) *
                                                                    pow(hold_cos, n);
                                                        }
                                                    }
                                                }

                                            }
                                            else {

                                                // Calculate the radius
                                                hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                                if (hold_r < rad_max) {

                                                    // Find the bin where this radius lies in
                                                    auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                    hold_r);
                                                    index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                                    if (index < resolution){
                                                        hold_cor = (hold_r - rad_bin_val.at(index - 1))/ rad_bin_size;
                                                        hold_cos = hold_y / hold_r;
                                                        // Add sin(theta)cos n(theta) to bin n
                                                        for (int n = 0; n < norder; ++n) {
                                                            rad_bin_shutter_on_first.at(n * resolution + index - 1) +=
                                                                    hold_weight *
                                                                    hold_cor *
                                                                    pow(hold_cos, n);
                                                            rad_bin_shutter_on_first.at(n * resolution + index) +=
                                                                    hold_weight *
                                                                    (1. - hold_cor) *
                                                                    pow(hold_cos, n);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                            event_time->at(j) < time_sub_upper && event_time->at(j) > time_sub_lower &&
                                            peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {
                                            injections_time_last += hold_weight;

                                            hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                            hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                            if (circularise) {
                                                // Calculate the radius and radially correct it
                                                hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                                if (hold_r < rad_max) {

                                                    //Find the bin where the theta lies in
                                                    if (hold_x >= 0) { //Upper plane
                                                        auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                                        std::acos(hold_y / hold_r));
                                                        ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                    } else { //Lower plane
                                                        auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                                        -std::acos(hold_y / hold_r));
                                                        ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                    }

                                                    //Find index correction for integral
                                                    auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                    hold_r * rad_corrections.at(ang_index));
                                                    index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                                    if (index < resolution){
                                                        hold_cor = (hold_r * rad_corrections.at(ang_index) - rad_bin_val.at(index - 1)) / rad_bin_size;
                                                        hold_cos = hold_y / hold_r;
                                                        // Add sin(theta)cos n(theta) to bin n
                                                        for (int n = 0; n < norder; ++n) {
                                                            rad_bin_shutter_on_last.at(n * resolution + index - 1) +=
                                                                    hold_weight *
                                                                    hold_cor *
                                                                    pow(hold_cos, n);
                                                            rad_bin_shutter_on_last.at(n * resolution + index) +=
                                                                    hold_weight *
                                                                    (1. - hold_cor) *
                                                                    pow(hold_cos, n);
                                                        }
                                                    }
                                                }

                                            }
                                            else {

                                                // Calculate the radius
                                                hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                                if (hold_r < rad_max) {

                                                    // Find the bin where this radius lies in
                                                    auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                    hold_r);
                                                    index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                                    if (index < resolution){
                                                        hold_cor = (hold_r - rad_bin_val.at(index - 1)) / rad_bin_size;
                                                        hold_cos = hold_y / hold_r;
                                                        // Add sin(theta)cos n(theta) to bin n
                                                        for (int n = 0; n < norder; ++n) {
                                                            rad_bin_shutter_on_last.at(n * resolution + index - 1) +=
                                                                    hold_weight *
                                                                    hold_cor *
                                                                    pow(hold_cos, n);
                                                            rad_bin_shutter_on_last.at(n * resolution + index) +=
                                                                    hold_weight *
                                                                    (1. - hold_cor) *
                                                                    pow(hold_cos, n);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                } else {
                                    hold_weight = weight(gen);
                                    injections_laser_off += hold_weight;
                                    for (int j = 0; j < X->size(); j++) {
                                        // Define time window
                                        if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                            event_time->at(j) < frame_upper && event_time->at(j) > frame_lower &&
                                            peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {
                                            injections_off_time_first += hold_weight;

                                            hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                            hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                            if (circularise) {

                                                // Calculate the radius and radially correct it
                                                hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                                if (hold_r < rad_max*max_rad_corr) {
                                                    //Find the bin where the theta lies in
                                                    if (hold_x >= 0) { //Upper plane
                                                        auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                                        std::acos(hold_y / hold_r));
                                                        ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                    } else { //Lower plane
                                                        auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                                        -std::acos(hold_y / hold_r));
                                                        ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                    }

                                                    //Find index correction for integral
                                                    auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                    hold_r * rad_corrections.at(ang_index));
                                                    index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                                    if (index < resolution){
                                                        hold_cor = (hold_r*rad_corrections.at(ang_index) - rad_bin_val.at(index - 1)) /
                                                                   rad_bin_size;
                                                        hold_cos = hold_y / hold_r;
                                                        // Add sin(theta)cos n(theta) to bin n
                                                        for (int n = 0; n < norder; ++n) {
                                                            rad_bin_shutter_off_first.at(n * resolution + index) +=
                                                                    hold_weight *
                                                                    hold_cor *
                                                                    pow(hold_cos, n);
                                                            rad_bin_shutter_off_first.at(n * resolution + index - 1) +=
                                                                    hold_weight *
                                                                    (1. - hold_cor) *
                                                                    pow(hold_cos, n);
                                                        }
                                                    }
                                                }

                                            }
                                            else {
                                                // Calculate the radius
                                                hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                                if (hold_r < rad_max) {

                                                    // Find the bin where this radius lies ini
                                                    auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                    hold_r);
                                                    index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                                    if (index < resolution){
                                                        hold_cor = (hold_r - rad_bin_val.at(index - 1)) /
                                                                   rad_bin_size;
                                                        hold_cos = hold_y / hold_r;
                                                        // Add sin(theta)cos n(theta) to bin n
                                                        for (int n = 0; n < norder; ++n) {
                                                            rad_bin_shutter_off_first.at(n * resolution + index) +=
                                                                    hold_weight *
                                                                    hold_cor *
                                                                    pow(hold_cos, n);
                                                            rad_bin_shutter_off_first.at(n * resolution + index - 1) +=
                                                                    hold_weight *
                                                                    (1. - hold_cor) *
                                                                    pow(hold_cos, n);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                            event_time->at(j) < time_sub_upper && event_time->at(j) > time_sub_lower &&
                                            peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {
                                            injections_off_time_last += hold_weight;

                                            hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                            hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                            if (circularise) {

                                                // Calculate the radius and radially correct it
                                                hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                                if (hold_r < rad_max*max_rad_corr) {
                                                    //Find the bin where the theta lies in
                                                    if (hold_x >= 0) { //Upper plane
                                                        auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                                        std::acos(hold_y / hold_r));
                                                        ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                    } else { //Lower plane
                                                        auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                                        -std::acos(hold_y / hold_r));
                                                        ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                    }

                                                    //Find index correction for integral
                                                    auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                    hold_r * rad_corrections.at(ang_index));
                                                    index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                                    if (index < resolution){
                                                        hold_cor = (hold_r * rad_corrections.at(ang_index) - rad_bin_val.at(index - 1)) /
                                                                   rad_bin_size;
                                                        hold_cos = hold_y / hold_r;
                                                        // Add sin(theta)cos n(theta) to bin n
                                                        for (int n = 0; n < norder; ++n) {
                                                            rad_bin_shutter_off_last.at(n * resolution + index) +=
                                                                    hold_weight *
                                                                    hold_cor *
                                                                    pow(hold_cos, n);
                                                            rad_bin_shutter_off_last.at(n * resolution + index - 1) +=
                                                                    hold_weight *
                                                                    (1. - hold_cor) *
                                                                    pow(hold_cos, n);
                                                        }
                                                    }
                                                }

                                            }
                                            else {
                                                // Calculate the radius
                                                hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                                if (hold_r < rad_max) {

                                                    // Find the bin where this radius lies ini
                                                    auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                    hold_r);
                                                    index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                                    if (index < resolution){
                                                        hold_cor = (hold_r - rad_bin_val.at(index - 1)) /
                                                                   rad_bin_size;
                                                        hold_cos = hold_y / hold_r;
                                                        // Add sin(theta)cos n(theta) to bin n
                                                        for (int n = 0; n < norder; ++n) {
                                                            rad_bin_shutter_off_last.at(n * resolution + index) +=
                                                                    hold_weight *
                                                                    hold_cor *
                                                                    pow(hold_cos, n);
                                                            rad_bin_shutter_off_last.at(n * resolution + index - 1) +=
                                                                    hold_weight *
                                                                    (1. - hold_cor) *
                                                                    pow(hold_cos, n);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                            // Removing background and saving the data
                            for (int i = 0; i < resolution * norder; ++i) {
                                // Sum, r_on - inject_on / inject_off * r_off
                                rad_bin_tot.at(i) = (rad_bin_shutter_on_first.at(i) - injections_laser_on / injections_laser_off *
                                                                                      rad_bin_shutter_off_first.at(i)) - (injections_time_first / injections_time_last * rad_bin_shutter_on_last.at(i) -
                                                                                                                          injections_off_time_first / injections_off_time_last * injections_laser_on / injections_laser_off *
                                                                                                                          rad_bin_shutter_off_last.at(i));
                            }
                        }
                    }
                    else if (no_laser) {
                        if (!time_subtraction) {
                            //Reading the data
                            for (int i = 0; i < entries; i++) {
                                X_branch->GetEntry(tree->LoadTree(i));
                                Y_branch->GetEntry(tree->LoadTree(i));
                                trig_time_branch->GetEntry(tree->LoadTree(i));
                                event_time_branch->GetEntry(tree->LoadTree(i));
                                peak_sum_branch->GetEntry(tree->LoadTree(i));
                                shutter_branch->GetEntry(tree->LoadTree(i));

                                // Bin the data to the resolution
                                if (!shutter) {
                                    hold_weight = weight(gen);
                                    for (int j = 0; j < X->size(); j++) {
                                        // Define time window
                                        if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                            event_time->at(j) < frame_upper && event_time->at(j) > frame_lower &&
                                            peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {

                                            hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                            hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                            if (circularise) {
                                                // Calculate the radius and radially correct it
                                                hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                                if (hold_r < rad_max) {

                                                    //Find the bin where the theta lies in
                                                    if (hold_x >= 0) { //Upper plane
                                                        auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                                        std::acos(hold_y / hold_r));
                                                        ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                    } else { //Lower plane
                                                        auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                                        -std::acos(hold_y / hold_r));
                                                        ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                    }

                                                    //Find index correction for integral
                                                    auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                    hold_r * rad_corrections.at(ang_index));
                                                    index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                                    if (index < resolution){
                                                        hold_cor = (hold_r * rad_corrections.at(ang_index) -
                                                                    rad_bin_val.at(index - 1)) / rad_bin_size;
                                                        hold_cos = hold_y / hold_r;
                                                        // Add sin(theta)cos n(theta) to bin n
                                                        for (int n = 0; n < norder; ++n) {
                                                            rad_bin_shutter_off.at(n * resolution + index - 1) +=
                                                                    hold_weight *
                                                                    hold_cor *
                                                                    pow(hold_cos, n);
                                                            rad_bin_shutter_off.at(n * resolution + index) +=
                                                                    hold_weight *
                                                                    (1. - hold_cor) *
                                                                    pow(hold_cos, n);
                                                        }
                                                    }
                                                }

                                            }
                                            else {

                                                // Calculate the radius
                                                hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                                if (hold_r < rad_max) {

                                                    // Find the bin where this radius lies in
                                                    auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                    hold_r);
                                                    index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                                    if (index < resolution){
                                                        hold_cor = (hold_r - rad_bin_val.at(index - 1))/
                                                                   rad_bin_size;
                                                        hold_cos = hold_y / hold_r;
                                                        // Add sin(theta)cos n(theta) to bin n
                                                        for (int n = 0; n < norder; ++n) {
                                                            rad_bin_shutter_off.at(n * resolution + index - 1) +=
                                                                    hold_weight *
                                                                    hold_cor *
                                                                    pow(hold_cos, n);
                                                            rad_bin_shutter_off.at(n * resolution + index) +=
                                                                    hold_weight *
                                                                    (1. - hold_cor)  *
                                                                    pow(hold_cos, n);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                            // Removing background and saving the data
                            for (int i = 0; i < resolution * norder; ++i) {
                                // Sum, r_on - inject_on / inject_off * r_off
                                rad_bin_tot.at(i) = rad_bin_shutter_off.at(i);
                            }
                        }
                        else {
                            //Reading the data
                            for (int i = 0; i < entries; i++) {
                                X_branch->GetEntry(tree->LoadTree(i));
                                Y_branch->GetEntry(tree->LoadTree(i));
                                trig_time_branch->GetEntry(tree->LoadTree(i));
                                event_time_branch->GetEntry(tree->LoadTree(i));
                                peak_sum_branch->GetEntry(tree->LoadTree(i));
                                shutter_branch->GetEntry(tree->LoadTree(i));

                                // Bin the data to the resolution
                                if (!shutter) {
                                    hold_weight = weight(gen);
                                    for (int j = 0; j < X->size(); j++) {
                                        // Define time window
                                        if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                            event_time->at(j) < frame_upper && event_time->at(j) > frame_lower &&
                                            peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {
                                            injections_time_first += hold_weight;

                                            hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                            hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                            if (circularise) {
                                                // Calculate the radius and radially correct it
                                                hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                                if (hold_r < rad_max) {

                                                    //Find the bin where the theta lies in
                                                    if (hold_x >= 0) { //Upper plane
                                                        auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                                        std::acos(hold_y / hold_r));
                                                        ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                    } else { //Lower plane
                                                        auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                                        -std::acos(hold_y / hold_r));
                                                        ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                    }

                                                    //Find index correction for integral
                                                    auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                    hold_r * rad_corrections.at(ang_index));
                                                    index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                                    if (index < resolution){
                                                        hold_cor = (hold_r * rad_corrections.at(ang_index) -
                                                                    rad_bin_val.at(index - 1)) / rad_bin_size;
                                                        hold_cos = hold_y / hold_r;
                                                        // Add sin(theta)cos n(theta) to bin n
                                                        for (int n = 0; n < norder; ++n) {
                                                            rad_bin_shutter_off_first.at(n * resolution + index - 1) +=
                                                                    hold_weight *
                                                                    hold_cor *
                                                                    pow(hold_cos, n);
                                                            rad_bin_shutter_off_first.at(n * resolution + index) +=
                                                                    hold_weight *
                                                                    (1. - hold_cor) *
                                                                    pow(hold_cos, n);
                                                        }
                                                    }
                                                }

                                            }
                                            else {

                                                // Calculate the radius
                                                hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                                if (hold_r < rad_max) {

                                                    // Find the bin where this radius lies in
                                                    auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                    hold_r);
                                                    index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                                    if (index < resolution){
                                                        hold_cor = (hold_r - rad_bin_val.at(index - 1))/ rad_bin_size;
                                                        hold_cos = hold_y / hold_r;
                                                        // Add sin(theta)cos n(theta) to bin n
                                                        for (int n = 0; n < norder; ++n) {
                                                            rad_bin_shutter_off_first.at(n * resolution + index - 1) +=
                                                                    hold_weight *
                                                                    hold_cor *
                                                                    pow(hold_cos, n);
                                                            rad_bin_shutter_off_first.at(n * resolution + index) +=
                                                                    hold_weight *
                                                                    (1. - hold_cor) *
                                                                    pow(hold_cos, n);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                                            event_time->at(j) < time_sub_upper && event_time->at(j) > time_sub_lower &&
                                            peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {
                                            injections_time_last += hold_weight;

                                            hold_x = ((X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;
                                            hold_y = (-(X->at(j) - x_center) + (Y->at(j) - y_center)) / M_SQRT2;

                                            if (circularise) {
                                                // Calculate the radius and radially correct it
                                                hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                                if (hold_r < rad_max) {

                                                    //Find the bin where the theta lies in
                                                    if (hold_x >= 0) { //Upper plane
                                                        auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                                        std::acos(hold_y / hold_r));
                                                        ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                    } else { //Lower plane
                                                        auto lower_t = std::lower_bound(angular_bin_val.begin(), angular_bin_val.end(),
                                                                                        -std::acos(hold_y / hold_r));
                                                        ang_index = (int) std::distance(angular_bin_val.begin(), lower_t);
                                                    }

                                                    //Find index correction for integral
                                                    auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                    hold_r * rad_corrections.at(ang_index));
                                                    index = (int) std::distance(rad_bin_val.begin(), lower_r);
                                                    if (index < resolution){
                                                        hold_cor = (hold_r * rad_corrections.at(ang_index) - rad_bin_val.at(index - 1)) / rad_bin_size;
                                                        hold_cos = hold_y / hold_r;
                                                        // Add sin(theta)cos n(theta) to bin n
                                                        for (int n = 0; n < norder; ++n) {
                                                            rad_bin_shutter_off_last.at(n * resolution + index - 1) +=
                                                                    hold_weight *
                                                                    hold_cor *
                                                                    pow(hold_cos, n);
                                                            rad_bin_shutter_off_last.at(n * resolution + index) +=
                                                                    hold_weight *
                                                                    (1. - hold_cor) *
                                                                    pow(hold_cos, n);
                                                        }
                                                    }
                                                }

                                            }
                                            else {

                                                // Calculate the radius
                                                hold_r = pow(hold_x*hold_x + hold_y*hold_y, 1. / 2.);

                                                if (hold_r < rad_max) {

                                                    // Find the bin where this radius lies in
                                                    auto lower_r = std::lower_bound(rad_bin_val.begin(), rad_bin_val.end(),
                                                                                    hold_r);
                                                    index = (int) std::distance(rad_bin_val.begin(), lower_r);

                                                    if (index < resolution){
                                                        hold_cor = (hold_r - rad_bin_val.at(index - 1)) / rad_bin_size;
                                                        hold_cos = hold_y / hold_r;
                                                        // Add sin(theta)cos n(theta) to bin n
                                                        for (int n = 0; n < norder; ++n) {
                                                            rad_bin_shutter_off_last.at(n * resolution + index - 1) +=
                                                                    hold_weight *
                                                                    hold_cor *
                                                                    pow(hold_cos, n);
                                                            rad_bin_shutter_off_last.at(n * resolution + index) +=
                                                                    hold_weight *
                                                                    (1. - hold_cor) *
                                                                    pow(hold_cos, n);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                            // Removing background and saving the data
                            for (int i = 0; i < resolution * norder; ++i) {
                                // Sum, r_on - inject_on / inject_off * r_off
                                rad_bin_tot.at(i) = rad_bin_shutter_off_first.at(i) - injections_time_first / injections_time_last * rad_bin_shutter_off_last.at(i);
                            }
                        }
                    }

                    /*--------            Calculating Abel inversion             --------*/
                    //Calculating inverse angular intensity distribution
                    get_angular_intensity(norder, resolution, I_n,
                                          inv_hankel, abel_inv_mat, rad_bin_val, rad_bin_tot);

                    // Calculating all coefficients of legendre polynomials
                    get_legendre_coef(norder, resolution, legendre_coefficients,
                                      I_n, inv_legendre);

                    // Taking out the intensity
                    for (int i = 0; i < resolution; ++i) {
                        new_intensity[i] = 4 * M_PI * pow(rad_bin_val.at(i), 2.) * legendre_coefficients.at(i);
                    }

                    //Making a spline to evaluate the intensity at the old radial values
                    spline = new TSpline3("test", new_radial, new_intensity, resolution, "b1e1");

                    //Save the sample statistics.
                    for (int i = 0; i < resolution; ++i) {
                        I_iter.at(sample*resolution + i) = spline->Eval(ref_radial.at(i));
                    }

                    spline->Clear();

                    end_time = std::chrono::steady_clock::now();
                    elapsed_seconds = end_time - start_time;
                    running_sum += elapsed_seconds;

                    std::cout << "Sample " << sample + 1 << " done in " << elapsed_seconds.count() << " s\n" <<
                              "Expected done in " << (nr_samples - sample) * running_sum.count() / 60. / (sample+1) << " min."
                              << std::endl;
                }
            }

            /*--------            Making CI based on samples             --------*/
            std::cout << "Sampling done" << std::endl;
            std::cout << "Calculating CI frame " << fr+1 << std::endl;
            count = 0;
            //Sorting the sample values (making the cdf)
            double place_holder;
            for (int i = 0; i < resolution; ++i) {
                for (int k = nr_samples-1; k >= 0; --k) {
                    for (int j = 0; j < k; ++j) {
                        if (I_iter.at(j * resolution + i) > I_iter.at((j + 1) * resolution + i)) {
                            place_holder = I_iter.at(j * resolution + i);
                            I_iter.at(j * resolution + i) = I_iter.at((j + 1) * resolution + i);
                            I_iter.at((j + 1) * resolution + i) = place_holder;
                        }
                    }
                }
            }

            //Calculating the bias correction z_0
            for (int i = 0; i < resolution; ++i) {
                //Take out the iteration needed
                for (int j = 0; j < nr_samples; ++j) {
                    iteration.at(j) = I_iter.at(j*resolution + i);
                }
                // Find the bin where the original I lie in
                auto lower_bin = std::lower_bound(iteration.begin(), iteration.end(),intensity.at(i));
                // Get the cdf(original I)
                place_holder = static_cast<double>(std::distance(iteration.begin(), lower_bin))/(nr_samples-1);
                z.at(i) = M_SQRT2 * ErfInv(2 * place_holder - 1); // Inverse Norm(0,1) CDF of cdf(original I)
            }

            //Calculating empirical influence function of I **Not possible with nearest neighbor integration **Not viable anyways

            //Calculating the acceleration constant a **Not possible with nearest neighbor integration **Not viable anyways

            //Calculating upper and lower bounds on 2*std % CI
            double upper_percent = 0.841345;
            double lower_percent = 0.158655;
            norm_percent = M_SQRT2* ErfInv(2*upper_percent-1);  //z^alpha, alpha = 0.841345, mean+1std

            for (int i = 0; i < resolution; ++i) {
                //Catch if statistic is outside of sample cdf
                if (std::isinf(z.at(i)) || std::isnan(z.at(i))){
                    if (std::signbit(z.at(i))) {
                        I_lower.at(i) = I_iter.at(0 * resolution + i);

                        lower_index = static_cast<int>((nr_samples - 1) * upper_percent);
                        I_upper.at(i) = I_iter.at(lower_index * resolution + i) +
                                        // + however far along a straight line to the next point
                                        (I_iter.at((lower_index + 1) * resolution + i) -
                                         I_iter.at(lower_index * resolution + i)) *
                                        (upper_percent * (nr_samples - 1) - lower_index);
                    }
                    else {
                        lower_index = static_cast<int>((nr_samples - 1) * lower_percent);
                        I_lower.at(i) = I_iter.at(lower_index * resolution + i) +
                                        // + however far along a straight line to the next point
                                        (I_iter.at((lower_index + 1) * resolution + i) -
                                         I_iter.at(lower_index * resolution + i)) *
                                        (lower_percent * (nr_samples - 1) - lower_index);;

                        I_upper.at(i) = I_iter.at((nr_samples-1) * resolution + i);
                    }
                }
                    //Calculate the CI
                else {
                    //      Norm(z_alpha) = 1/2 * [1 + erf (z_alpha/sqrt(2))],   z_alpha = z_0 + (z_0 + z^alpha) / (1 - a * (z_0 + z^alpha))
                    //       z^alpha is the alpha percentile of the normal distribution
                    place_holder = (nr_samples - 1) * 1. / 2. *
                                   (1. + std::erf((2 * z.at(i) - norm_percent) / M_SQRT2));
                    lower_index = static_cast<int>(place_holder);
                    if (lower_index >= nr_samples-1){
                        lower_index = nr_samples-2;
                    }

                    I_lower.at(i) = I_iter.at(lower_index * resolution + i) +
                                    // + however far along a straight line to the next point
                                    (I_iter.at((lower_index + 1) * resolution + i) -
                                     I_iter.at(lower_index * resolution + i)) *
                                    (place_holder - lower_index);


                    place_holder = (nr_samples - 1) * 1. / 2. *
                                   (1. + std::erf((2 * z.at(i) + norm_percent) / M_SQRT2));
                    lower_index = static_cast<int>(place_holder);
                    if (lower_index >= nr_samples-1){
                        lower_index = nr_samples-2;
                    }

                    I_upper.at(i) = I_iter.at(lower_index * resolution + i) +
                                    // + however far along a straight line to the next point
                                    (I_iter.at((lower_index + 1) * resolution + i) -
                                     I_iter.at(lower_index * resolution + i)) *
                                    (place_holder - lower_index);
                }
            }

            gif_file->cd();
            auto *I_err_l_branch = gif_tree->Branch("intensity_err_lower", &I_val_err_l);
            auto *I_err_u_branch = gif_tree->Branch("intensity_err_upper", &I_val_err_u);
            auto *I_iter_branch = gif_tree->Branch("iteration", &I_iter_val);
            for (auto &ent: I_lower) {
                I_val_err_l = ent;
                I_err_l_branch->Fill();
            }
            for (auto &ent: I_upper) {
                I_val_err_u = ent;
                I_err_u_branch->Fill();
            }
            for (auto &ent: I_iter) {
                I_iter_val = ent;
                I_iter_branch->Fill();
            }

            gif_tree->Write(nullptr, TObject::kOverwrite);
        }

        frame_end_time = std::chrono::steady_clock::now();
        frame_elapsed_seconds = frame_end_time - frame_start_time;
        frame_running_sum += frame_elapsed_seconds;

        std::cout << "Frame " << fr + 1 << " done in " << frame_elapsed_seconds.count() << " s\n" <<
                  "Expected done in " << (gif_nr_frames - fr + 1) * frame_running_sum.count() / 60. / (fr + 1)<< " min."
                  << std::endl;
    }


    double settings_val;
    auto *gif_tree = new TTree("SETTINGS","SETTINGS");
    auto *settings_branch = gif_tree->Branch("settings", &settings_val);
    settings_val = min_intensity; settings_branch->Fill();
    settings_val = max_intensity; settings_branch->Fill();
    settings_val = max_radius; settings_branch->Fill();
    settings_val = gif_nr_frames; settings_branch->Fill();
    settings_val = min_radial_tot; settings_branch->Fill();
    settings_val = max_radial_tot; settings_branch->Fill();
    gif_tree->Write(nullptr, TObject::kOverwrite);

    // Closing up
    gif_file->Close();
    data_file->Close();

    std::cout << "All done\n" << std::endl;

    return 0;
}