#include"PowerSpectrum.h"

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert,
    double A_s,
    double n_s,
    double kpivot_mpc) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert),
  A_s(A_s),
  n_s(n_s),
  kpivot_mpc(kpivot_mpc)
{}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(){

  generate_bessel_function_splines();

  line_of_sight_integration();

  auto cell_TT = solve_for_cell(thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
  
  if(Constants.polarization){
    auto cell_EE = solve_for_cell(thetaE_ell_of_k_spline, thetaE_ell_of_k_spline);
    cell_EE_spline.create(ells, cell_EE, "Cell_EE_of_ell");
    
    auto cell_TE = solve_for_cell(thetaT_ell_of_k_spline, thetaE_ell_of_k_spline);
    cell_TE_spline.create(ells, cell_TE, "Cell_TE_of_ell");
  }
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());
   
  double eta_0 = cosmo->eta_of_x(0.0);
  Vector z_array = Utils::linspace(0.0, k_max*eta_0, n_k);
  for(size_t i = 0; i < ells.size(); i++){
    Vector j_ell_array(z_array.size());
    const int ell = ells[i];
    #pragma omp parallel for num_threads(omp_get_max_threads())
    for(int j = 0; j < z_array.size(); j++){
           j_ell_array[j] = Utils::j_ell(ell, z_array[j]);
    }

    // Make the j_ell_splines[i] spline
    Spline ell_spline{"ell"};
    ell_spline.create(z_array, j_ell_array, "ell");
    j_ell_splines[i] = ell_spline;
  }

  Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array,
    std::function<double(double,double)> &source_function){
  Utils::StartTiming("lineofsight");
    
  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));

  //Pre computing stuff to avoid reduntant calculations inside the loop
  Vector x_array = Utils::linspace(Constants.x_start, Constants.x_end, n_x);
  const double eta0 = cosmo->eta_of_x(0.0);
  Vector eta_array(n_x);
  for(size_t ix = 0; ix < n_x; ix++) {
      eta_array[ix] = cosmo->eta_of_x(x_array[ix]);
  }

  Vector2D source_array(k_array.size(), Vector(n_x));
  for(size_t ix = 0; ix < n_x; ix++) {
      for(size_t ik = 0; ik < k_array.size(); ik++) {
          source_array[ik][ix] = source_function(x_array[ix], k_array[ik]);
      }
  }

  #pragma omp parallel for num_threads(omp_get_max_threads())
  for(size_t il = 0; il < ells.size(); il++){
    for(size_t ik = 0; ik < k_array.size(); ik++){
    
    // Trapezoidal rule states that int_a^b f(x) dx ~ (b-a)/2 * (f(a) + f(b)), so we use this approximation at every x-step
      double Theta_l = 0;
      double x0;
      double x1 = x_array[0];
      double f0;
      double f1 = source_array[ik][0] * j_ell_splines[il](k_array[ik] * (eta0 - eta_array[0]));
      
      for(size_t ix=0; ix < n_x-1; ix++){
      // We avoid computing both f(a) and f(b) at every step by letting the new f(a) equal the old f(b) and recomputing f(b)
      // This should roughly halve the computation time that this loop takes
        x0 = x1;
        x1 = x_array[ix+1];
        f0 = f1;
        f1 = source_array[ik][ix+1] * j_ell_splines[il](k_array[ik] * (eta0 - eta_array[ix+1]));
        Theta_l += (f1 + f0) * (x1 - x0);
      }
      result[il][ik] = Theta_l * 0.5;
    }
  }

  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(){
  const int nells      = ells.size();
  // Here we need an N >~ 6
  int N = 10;
  const double dk = 2.0 * pi / (cosmo->eta_of_x(0.0) * N);
  const int nk = int( (k_max - k_min) / dk) + 1;
  Vector k_array = Utils::linspace(k_min, k_max, nk);

  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);
  thetaE_ell_of_k_spline = std::vector<Spline>(nells);

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };
  
  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);
  
  for(int il = 0; il < nells; il++){
    thetaT_ell_of_k_spline[il].create(k_array, thetaT_ell_of_k[il]);
  }

  if(Constants.polarization){
    // Make a function returning the source function for polarization
    std::function<double(double,double)> source_function_E = [&](double x, double k){
      return pert->get_Source_E(x,k);
    };

    // Do the line of sight integration
    Vector2D thetaE_ell_of_k = line_of_sight_integration_single(k_array, source_function_E);

    for(int il = 0; il < nells; il++){
      double ell = ells[il];
      double factor = sqrt( (ell+2.0) * (ell-1.0) * ell * (ell+1.0) );
      for(int ik = 0; ik < nk; ik++){
        thetaE_ell_of_k[il][ik] *= factor;
      }
      thetaE_ell_of_k_spline[il].create(k_array, thetaE_ell_of_k[il]);
    }
  }
}

//====================================================
// Compute Cell (could be TT or TE or EE) 
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================
Vector PowerSpectrum::solve_for_cell(
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
  
  const int nells = ells.size();
  int N = 64;
  const double dk = 2.0 * pi / (cosmo->eta_of_x(0.0) * N);
  const int nk = int( (k_max - k_min) / dk) + 1;
  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), nk);
  
  Vector result(nells);
  // Use 2 pi instead of 4 pi because the trapezoid rule tells us to divide by 2
  const double Cl_prefactor = 2.0 * pi * A_s * pow( kpivot_mpc / Constants.Mpc , 1.0 - n_s);

  for(int il = 0; il < nells; il++){
    double Cell = 0.0;
    
    // We use the log of the k values for integration
    double log_k0;
    double log_k1 = log_k_array[0];
    
    double f0;
    double f1 = pow(k_min, n_s - 1.0) * f_ell_spline[il](k_min) * g_ell_spline[il](k_min);
    double k1 = exp(log_k1);
    
    // Integrand: P(k) * Theta^2. P(k) = A_s * (k/kpivot)^(ns-1)
    // Note: pow(k1, n_s - 1.0) is the k part of the primordial spectrum
    f0 = pow(k1, n_s - 1.0) * f_ell_spline[il](k1) * g_ell_spline[il](k1);

    for(int ik = 0; ik < nk-1; ik++){
      log_k0 = log_k1;
      log_k1 = log_k_array[ik+1];
      
      k1 = exp(log_k1);
      f0 = f1;
      f1 = pow(k1, n_s - 1.0) * f_ell_spline[il](k1) * g_ell_spline[il](k1);
      
      // Integrate over d(ln k)
      Cell += (f1 + f0) * (log_k1 - log_k0);
      
    }
    result[il] = Cell * Cl_prefactor;
  }
  
  return result;
}

//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow( Constants.Mpc * k / kpivot_mpc , n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
//====================================================


double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const{
  double c = Constants.c;
  double k = k_mpc/Constants.Mpc;
  double a = exp(x);
  double P_primordial = 2.0*pi*pi*primordial_power_spectrum(k)/(k*k*k);
  double Phi = pert->get_Phi(x,k);
  double Omega_M_0 = cosmo->get_OmegaCDM(0.0) + cosmo->get_OmegaB(0.0);
  double H_0 = cosmo->get_H0();
  double Delta_M = (c*c*k*k*Phi*(2.0/3.0)*a)/(Omega_M_0*H_0*H_0);
  double pofk = Delta_M*Delta_M*P_primordial;
  return pofk;
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}
double PowerSpectrum::get_cell_TE(const double ell) const{
  return cell_TE_spline(ell);
}
double PowerSpectrum::get_cell_EE(const double ell) const{
  return cell_EE_spline(ell);
}

//====================================================
// Output the cells to file
//====================================================

void PowerSpectrum::output_ell(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    double normfactorN = (ell * (ell+1)) / (2.0 * M_PI) 
      * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);
    fp << ell                                 << " ";
    fp << cell_TT_spline( ell ) * normfactor  << " ";
    if(Constants.polarization){
      fp << cell_EE_spline( ell ) * normfactor  << " ";
      fp << cell_TE_spline( ell ) * normfactor  << " ";
    }
    else{
      fp << 0 << " ";
      fp << 0 << " ";
    }
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}

void PowerSpectrum::output_k(std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int k_number = 1000;
  const double h = cosmo->get_h();
  auto kvalues = Utils::linspace(Constants.k_min*Constants.Mpc, Constants.k_max*Constants.Mpc, k_number);
  auto print_data = [&] (const double k_mpc) {
    
    fp << k_mpc/h                                         << " ";
    fp << get_matter_power_spectrum(0.0, k_mpc)*h*h*h/(Constants.Mpc*Constants.Mpc*Constants.Mpc) << " ";
    fp << thetaT_ell_of_k_spline[9](k_mpc/Constants.Mpc)  << " ";   //ell=15
    fp << thetaT_ell_of_k_spline[32](k_mpc/Constants.Mpc) << " ";   //ell=500
    fp << thetaT_ell_of_k_spline[42](k_mpc/Constants.Mpc) << " ";   //ell=1000
    fp << thetaE_ell_of_k_spline[9](k_mpc/Constants.Mpc)  << " ";   //ell=15
    fp << thetaE_ell_of_k_spline[32](k_mpc/Constants.Mpc) << " ";   //ell=500
    fp << thetaE_ell_of_k_spline[42](k_mpc/Constants.Mpc) << " ";   //ell=1000 
    fp << "\n";
  };
  std::for_each(kvalues.begin(), kvalues.end(), print_data);
}