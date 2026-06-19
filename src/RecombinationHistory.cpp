#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp = 0) :
  cosmo(cosmo),
  Yp(Yp)
{
    //Physical constants
    const double k_b         = Constants.k_b;
    const double G           = Constants.G;
    const double m_e         = Constants.m_e;
    const double hbar        = Constants.hbar;
    const double c           = Constants.c;
    const double epsilon_0   = Constants.epsilon_0;

    const double H0 = cosmo->get_H0();
    const double alpha = 1/137.0359992; // Fine structure constant

    // const double phi_2 = 0.448 * log(epsilon_0/(k_b*T_b));
    // const double alpha = 1/137.0359992; // Fine structure constant
    // const double alpha_2 = 64 * sqrt(pi/27) * alpha*alpha*hbar*hbar/(m_e*m_e*c) * sqrt(epsilon_0/(k_b*T_b)) * phi_2;
    // const double beta = alpha_2 * pow(m_e*k_b*T_b/(2*pi*hbar*hbar), 1.5) * exp(-epsilon_0/(k_b*T_b));
    // const double beta_2 = beta * exp(3*epsilon_0/(4*k_b*T_b));
    // const double H = cosmo->H_of_x_over_H0(x) * H0;
    // const double n_1s = n_b * (1 - X_e);
    // const double lambda_alpha = H * pow(3*epsilon_0/(hbar*c), 3) / (64*pi*pi*n_1s);
    // const double C_r = (lambda_2s1s + lambda_alpha) / (lambda_2s1s + lambda_alpha + beta_2);

    // Pre compute constants we will need in the Peebles equation
    rho_c0 = 3*H0*H0 / (8*pi*G);
    Lambda_alpha_coefficient = pow(3*epsilon_0/(hbar*c), 3) / (64*pi*pi);
    beta_coefficient =  pow(m_e*k_b/(2*pi*hbar*hbar), 1.5);
    alpha_coefficient = 64 * sqrt(pi/27) * pow(alpha*hbar/m_e, 2) / c * 0.448 * sqrt(epsilon_0/k_b);
}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  
  //=============================================================================
  // TODO: Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector Xe_arr(npts_rec_arrays);
  Vector Xe_reionization(npts_rec_arrays);
  Vector Xe_saha_arr(npts_rec_arrays);
  Vector ne_arr(npts_rec_arrays);
  Vector ne_reionization(npts_rec_arrays);

  Xe_arr[0] = 1.0; // Initial condition for Xe at early times (before recombination) is 1.0
  Xe_reionization[0] = 1.0;
  Xe_saha_arr[0] = 1.0;
  double n_H0 = cosmo->get_OmegaB(0.0) * 3*cosmo->get_H0()*cosmo->get_H0() / (8*pi*Constants.G * Constants.m_H);
  ne_arr[0] = n_H0; // Initial condition for ne at early times (before recombination) is the baryon number density today

  // Calculate recombination history
  bool saha_regime = true;

  double Xe_current;
  double ne_current;

  for(int i = 0; i < npts_rec_arrays; i++){

    //==============================================================
    // TODO: Get X_e from solving the Saha equation so
    // implement the function electron_fraction_from_saha_equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    Xe_current = Xe_ne_data.first;
    ne_current = Xe_ne_data.second;
    Xe_saha_arr[i] = Xe_current;
    ne_arr[i] = ne_current;
    
    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit)
      saha_regime = false;

    if(saha_regime){
      
      Xe_arr[i] = Xe_current;
      ne_arr[i] = ne_current;
      Xe_reionization[i] = Xe_arr[i]; // No reionization in the Saha regime since we are at early times
      ne_reionization[i] = ne_arr[i];

    } else {
      // Compute X_e from current time til today by solving the Peebles equation
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){

        return rhs_peebles_ode(x, Xe, dXedx);
      };
      // We need a reduced x-array that starts at the current x-value since we are solving from the current x-value to today
      Vector reduced_x_array;
      for(int j = i; j < npts_rec_arrays; j++){
        reduced_x_array.push_back(x_array[j]);
      }
      Vector initial_condition{Xe_current};
      peebles_Xe_ode.solve(dXedx, reduced_x_array, initial_condition);
      Vector Xe_peebles_array = peebles_Xe_ode.get_data_by_component(0);
      
      double z_reion = 8.0;
      double delta_z_reion = 0.5;
      double f_He = Yp / (4*(1 - Yp));
      double y_reion = pow((1 + z_reion),1.5);
      double delta_y_reion = 1.5 * sqrt(1 + z_reion) * delta_z_reion;
      
      // Store results and exit loop
      for(int j = i; j < npts_rec_arrays; j++){
      auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[j]);
  
      // Electron fraction and number density
      Xe_current = Xe_ne_data.first;
      ne_current = Xe_ne_data.second;
      Xe_saha_arr[j] = Xe_current;
      Xe_arr[j] = Xe_peebles_array[j - i];
      ne_arr[j] = Xe_arr[j] * n_H0 * exp(-3*x_array[j]);
      
      double y = exp(-3*x_array[j]/2);
      Xe_reionization[j] = Xe_arr[j] + 0.5 * (1 + f_He) * (1 + tanh((y_reion - y)/delta_y_reion));
      ne_reionization[j] = Xe_reionization[j] * n_H0 * exp(-3*x_array[j]);
      }
      break; // Exit the for-loop since we've solved the rest
    }  
  }
  
  // Find x_decoupling where Xe = 0.1
  int idx_decoupling = 0;
  double min_diff = std::abs(Xe_arr[0] - 0.1);
  for(int i = 1; i < npts_rec_arrays; i++){
    double diff = std::abs(Xe_arr[i] - 0.1);
    if(diff < min_diff){
      min_diff = diff;
      idx_decoupling = i;
    }
  }
  x_decoupling_recombination = x_array[idx_decoupling];

  int idx_decoupling_Saha = 0;
  min_diff = std::abs(Xe_saha_arr[0] - 0.1);
  for(int i = 1; i < npts_rec_arrays; i++){
    double diff = std::abs(Xe_saha_arr[i] - 0.1);
    if(diff < min_diff){
      min_diff = diff;
      idx_decoupling_Saha = i;
    }
  }
  x_decoupling_Saha_recombination = x_array[idx_decoupling_Saha];

  Xe_of_x_spline.create(x_array, Xe_arr, "Xe_of_x");
  Xe_reionization_of_x_spline.create(x_array, Xe_reionization, "Xe_reionization_of_x");
  Xe_saha_of_x_spline.create(x_array, Xe_saha_arr, "Xe_saha_of_x");

  Utils::EndTiming("Xe");
}

std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a           = exp(x);
 
  // Physical constants
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0;
  const double H0_over_h   = Constants.H0_over_h;

  const double T_b = cosmo->get_TCMB(x);
  const double OmegaB_0 = cosmo->get_OmegaB(0.0);
  const double H0 = cosmo->get_H0();

  
  double n_b = OmegaB_0 * rho_c0 / (m_H*a*a*a);
  
  // Electron fraction and number density
  double Xe = 2/(1 + sqrt(1 + 4*n_b*pow(2*pi*hbar*hbar/(m_e*k_b*T_b), 1.5)*exp(epsilon_0/(k_b*T_b))));  
  double ne = Xe * n_b;

  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double sigma_T     = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;

  

  // Cosmological parameters
  const double T_b = cosmo->get_TCMB(x);
  const double n_b = cosmo->get_OmegaB(0.0) * rho_c0 / (m_H*a*a*a);

  const double alpha_2 = alpha_coefficient * log(epsilon_0/(k_b*T_b)) / sqrt(T_b);
  const double beta = beta_coefficient * alpha_2 * T_b * sqrt(T_b) * exp(-epsilon_0/(k_b*T_b));
  const double beta_2 = beta * exp(3*epsilon_0/(4*k_b*T_b));
  const double H = cosmo->H_of_x_over_H0(x) * cosmo->get_H0();
  const double n_1s = std::max(n_b * (1.0 - X_e), 1e-10);
  const double lambda_alpha = Lambda_alpha_coefficient * H / n_1s;
  // const double C_r = (lambda_2s1s + lambda_alpha) / (lambda_2s1s + lambda_alpha + beta_2);

  double num = lambda_2s1s + lambda_alpha;
  double den = num + beta_2;
  double C_r = 1.0; 
  if (den > 0) C_r = num / den;
  if (std::isinf(num) || std::isinf(beta_2)) {
      // If things are exploding, C_r is effectively the ratio of the escape rates
      C_r = lambda_alpha / (lambda_alpha + beta_2); 
      if (std::isnan(C_r)) C_r = 1.0; // Fallback for pure Saha regime
}
  
  dXedx[0] = C_r * (beta * (1 - X_e) - n_b * alpha_2 * X_e*X_e) / H;

  return GSL_SUCCESS;
}

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over 
  // Since the IC is at x=0 (tau(0) = 0) the array should go from 0.0 -> x_start 
  const int npts = 1000;
  Vector x_array = Utils::linspace(0.0, x_start, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    double n_e = ne_of_x(x);
    double H = cosmo->H_of_x_over_H0(x) * cosmo->get_H0();

    // Set the derivative for photon optical depth
    // We flip the sign of the derivative because we integrate backwards
    dtaudx[0] = -Constants.c * n_e * Constants.sigma_T / H;

    return GSL_SUCCESS;
  };

  ODEFunction dtaudx_reion = [&](double x, const double *tau, double *dtaudx){

    double n_e = ne_reionization_of_x(x);
    double H = cosmo->H_of_x_over_H0(x) * cosmo->get_H0();

    // Set the derivative for photon optical depth
    // We flip the sign of the derivative because we integrate backwards
    dtaudx[0] = -Constants.c * n_e * Constants.sigma_T / H;

    return GSL_SUCCESS;
  };

  ODEFunction dtaudx_Saha = [&](double x, const double *tau, double *dtaudx){

    double n_e = electron_fraction_from_saha_equation(x).second;
    double H = cosmo->H_of_x_over_H0(x) * cosmo->get_H0();

    // Set the derivative for photon optical depth
    // We flip the sign of the derivative because we integrate backwards
    dtaudx[0] = -Constants.c * n_e * Constants.sigma_T / H;

    return GSL_SUCCESS;
  };

  ODESolver tau_ode;
  Vector initial_condition{0.0}; // tau(0) = 0
  tau_ode.solve(dtaudx, x_array, initial_condition);
  Vector tau_array = tau_ode.get_data_by_component(0);

  ODESolver tau_reion_ode;
  tau_reion_ode.solve(dtaudx_reion, x_array, initial_condition);
  Vector tau_reion_array = tau_reion_ode.get_data_by_component(0);

  ODESolver tau_Saha_ode;
  tau_Saha_ode.solve(dtaudx_Saha, x_array, initial_condition);
  Vector tau_Saha_array = tau_Saha_ode.get_data_by_component(0);

  tau_of_x_spline.create(x_array, tau_array, "tau_of_x");
  tau_reionization_of_x_spline.create(x_array, tau_reion_array, "tau_reionization");
  tau_saha_of_x_spline.create(x_array, tau_Saha_array, "tau_saha");

  Vector g_tilde_array(npts);
  Vector g_tilde_reionization_array(npts);
  Vector g_tilde_Saha_array(npts);

  for(int i = 0; i < npts; i++){
    g_tilde_array[i] = -dtaudx_of_x(x_array[i]) * exp(-tau_array[i]);
    g_tilde_reionization_array[i] = -dtaudx_reionization_of_x(x_array[i]) * exp(-tau_reion_array[i]);
    g_tilde_Saha_array[i] = -tau_saha_of_x_spline.deriv_x(x_array[i]) * exp(-tau_Saha_array[i]);
  }

  // Find x_last_scattering where g_tilde is maximum
  int idx_last_scattering = 0;
  double max_g_tilde = g_tilde_array[0];
  for(int i = 1; i < npts; i++){
    if(g_tilde_array[i] > max_g_tilde){
      max_g_tilde = g_tilde_array[i];
      idx_last_scattering = i;
    }
  }
  x_decoupling_scattering = x_array[idx_last_scattering];

  int idx_last_scattering_Saha = 0;
  double max_g_tilde_Saha = g_tilde_Saha_array[0];
  for(int i = 1; i < npts; i++){
    if(g_tilde_Saha_array[i] > max_g_tilde_Saha){
      max_g_tilde_Saha = g_tilde_Saha_array[i];
      idx_last_scattering_Saha = i;
    }
  }
  x_decoupling_Saha_scattering = x_array[idx_last_scattering_Saha];

  g_tilde_of_x_spline.create(x_array, g_tilde_array, "g_tilde_of_x");
  g_tilde_reionization_of_x_spline.create(x_array, g_tilde_reionization_array, "g_tilde_reionization_of_x");

  Utils::EndTiming("opticaldepth");

  // Compute the sound horizon at decoupling and spline it as well

  x_array = Utils::linspace(x_start, x_end, npts); // Here we want to integrate forward, so we redefine x_array

  const double R_minus1 = (3*cosmo->get_OmegaB(0.0)) / (4*cosmo->get_OmegaR(0.0));
  ODESolver sound_horizon_ode;

  ODEFunction dsdx = [&](double x, const double *s, double *dsdx){

    double c_s = Constants.c / sqrt(3*(1 + R_minus1*exp(x)));

    dsdx[0] = c_s / cosmo->Hp_of_x(x);

    return GSL_SUCCESS;
  };

  Vector initial_condition_sound_horizon{Constants.c * sqrt(1 / (3 * (1 + R_minus1))) / cosmo->Hp_of_x(x_start)};
  sound_horizon_ode.solve(dsdx, x_array, initial_condition_sound_horizon);

  Vector s_array = sound_horizon_ode.get_data_by_component(0);
  s_of_x_spline.create(x_array, s_array, "s_of_x");

}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{

  return tau_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{

  return tau_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::tau_reionization_of_x(double x) const{
  return tau_reionization_of_x_spline(x);
}

double RecombinationHistory::dtaudx_reionization_of_x(double x) const{

  return tau_reionization_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddtauddx_reionization_of_x(double x) const{

  return tau_reionization_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{

  return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{

  return g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::g_tilde_reionization_of_x(double x) const{
  return g_tilde_reionization_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_reionization_of_x(double x) const{

  return g_tilde_reionization_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_reionization_of_x(double x) const{

  return g_tilde_reionization_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::Xe_of_x(double x) const{

  return Xe_of_x_spline(x);
}

double RecombinationHistory::Xe_reionization_of_x(double x) const{

  return Xe_reionization_of_x_spline(x);
}

double RecombinationHistory::ne_of_x(double x) const{

  double n_H = cosmo->get_OmegaB(0.0) * 3*cosmo->get_H0()*cosmo->get_H0() / (8*pi*Constants.G * Constants.m_H) * exp(-3*x);

  return Xe_of_x_spline(x) * n_H;
}

double RecombinationHistory::ne_reionization_of_x(double x) const{

  double n_H = cosmo->get_OmegaB(0.0) * 3*cosmo->get_H0()*cosmo->get_H0() / (8*pi*Constants.G * Constants.m_H) * exp(-3*x);

  return Xe_reionization_of_x_spline(x) * n_H;
}

double RecombinationHistory::s_of_x(double x) const{
  return s_of_x_spline(x);
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  double Mpc = Constants.Mpc;
  double Myr = Constants.Gyr * 1000.0;
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << "\n";
  std::cout << "x_decoupling_recombination: " << x_decoupling_recombination << "\n";
  std::cout << "x_decoupling_recombination_Saha: " << x_decoupling_Saha_recombination << "\n";
  std::cout << "z_decoupling_recombination: " << exp(-x_decoupling_recombination) - 1 << "\n";
  std::cout << "z_decoupling_recombination_Saha: " << exp(-x_decoupling_Saha_recombination) - 1 << "\n";
  std::cout << "t_decoupling_recombination: " << cosmo->t_of_x(x_decoupling_recombination) * Myr << " Myr\n";
  std::cout << "t_decoupling_recombination_Saha: " << cosmo->t_of_x(x_decoupling_Saha_recombination) * Myr << " Myr\n";
  std::cout << "s_decoupling_recombination: " << s_of_x(x_decoupling_recombination) / Mpc << " Mpc\n";
  std::cout << "s_decoupling_recombination_Saha: " << s_of_x(x_decoupling_Saha_recombination) / Mpc << " Mpc\n";
  std::cout << "\n";
  std::cout << "x_decoupling_scattering: " << x_decoupling_scattering << "\n";
  std::cout << "x_decoupling_Saha_scattering: " << x_decoupling_Saha_scattering << "\n";
  std::cout << "z_decoupling_scattering: " << exp(-x_decoupling_scattering) - 1 << "\n";
  std::cout << "z_decoupling_Saha_scattering: " << exp(-x_decoupling_Saha_scattering) - 1 << "\n";
  std::cout << "t_decoupling_scattering: " << cosmo->t_of_x(x_decoupling_scattering) * Myr << " Myr\n";
  std::cout << "t_decoupling_Saha_scattering: " << cosmo->t_of_x(x_decoupling_Saha_scattering) * Myr << " Myr\n";
  std::cout << "s_decoupling_scattering: " << s_of_x(x_decoupling_scattering) / Mpc << " Mpc\n";
  std::cout << "s_decoupling_Saha_scattering: " << s_of_x(x_decoupling_Saha_scattering) / Mpc << " Mpc\n";
  std::cout << "\n";
  std::cout << "Freeze_out X_e today: " << Xe_of_x(0.0) << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 5000;
  const double x_min   = x_start;
  const double x_max   = x_end;
  const double Myr = Constants.Gyr * 1000.0;
  const double Mpc = Constants.Mpc;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << ne_of_x(x)           << " ";
    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)      << " ";
    fp << ddtauddx_of_x(x)     << " ";
    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x) << " ";
    fp << Xe_reionization_of_x(x)           << " ";
    fp << ne_reionization_of_x(x)           << " ";
    fp << tau_reionization_of_x(x)          << " ";
    fp << dtaudx_reionization_of_x(x)       << " ";
    fp << ddtauddx_reionization_of_x(x)     << " ";
    fp << g_tilde_reionization_of_x(x)      << " ";
    fp << dgdx_tilde_reionization_of_x(x)   << " ";
    fp << ddgddx_tilde_reionization_of_x(x) << " ";
    fp << Xe_saha_of_x_spline(x) << " ";
    fp << s_of_x(x) / Mpc << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

