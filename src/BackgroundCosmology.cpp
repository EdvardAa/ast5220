#include "BackgroundCosmology.h"
#include <numbers>

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM, 
    double OmegaK,
    double Neff, 
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaK(OmegaK),
  Neff(Neff), 
  TCMB(TCMB)
{
  double pi = std::numbers::pi;
  // Hubble parameter today
  H0 = h * Constants.H0_over_h;
  // Photon density today
  OmegaR = 16/90. * pow(pi,3) * pow(TCMB*Constants.k_b, 4) / pow(Constants.hbar, 3) / pow(Constants.c, 5) * Constants.G / pow(H0, 2);
  // Neutrino density today
  OmegaNu = Neff * 7/8. * pow(4/11., 4/3.) * OmegaR;
  // Dark energy density today
  OmegaLambda = 1.0 - OmegaK - OmegaB - OmegaCDM - OmegaR - OmegaNu;
}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
  // Utils::StartTiming("Eta");
    
  Vector x_array = Utils::linspace(x_start, x_end, 1000);

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){

    detadx[0] = Constants.c / Hp_of_x(x);

    return GSL_SUCCESS;
  };

  ODEFunction dtdx = [&](double x, const double *t, double *dtdx){

    dtdx[0] = 1.0 / (H0*H_of_x_over_H0(x));

    return GSL_SUCCESS;
  };

  double initial_eta = 0.0;
  Vector initial_condion{initial_eta};

  ODESolver ode;
  ode.solve(detadx, x_array, initial_condion);

  Vector eta_array = ode.get_data_by_component(0);

  eta_of_x_spline.create(x_array, eta_array, "eta_of_x");

  ODESolver ode_t;
  ode_t.solve(dtdx, x_array, initial_condion);

  Vector t_array = ode_t.get_data_by_component(0);

  t_of_x_spline.create(x_array, t_array, "t_of_x");

  // Utils::EndTiming("Eta");
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x_over_H0(double x) const{
  double a_minus1 = exp(-x);
  double a_minus2 = a_minus1 * a_minus1;
  double H_over_H0 = sqrt((OmegaB + OmegaCDM) * a_minus2*a_minus1 + (OmegaR + OmegaNu) * a_minus2*a_minus2 + OmegaK * a_minus2 + OmegaLambda);
  return H_over_H0;
}

double BackgroundCosmology::Hp_of_x(double x) const{
  double a = exp(x);
  double Hp = H0 * H_of_x_over_H0(x) * a;
  return Hp;
}

double BackgroundCosmology::dHpdx_of_x(double x) const{
  double a_minus1 = exp(-x);
  double a_minus2 = a_minus1 * a_minus1;
  double dHpdx = -H0*H0/(2*Hp_of_x(x)) * ((OmegaB + OmegaCDM)*a_minus1 + 2*(OmegaR + OmegaNu) * a_minus2 - 2*OmegaLambda / a_minus2);
  return dHpdx;
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{
  double a = exp(x);
  double Hp = Hp_of_x(x);
  double ddHpddx = H0*H0/(2*Hp) * (dHpdx_of_x(x)/Hp * ((OmegaB + OmegaCDM)/a + 2*(OmegaR + OmegaNu) / (a*a) - 2*OmegaLambda *a*a)
                                                    + ((OmegaB + OmegaCDM)/a + 4*(OmegaR + OmegaNu) / (a*a) + 4*OmegaLambda *a*a) );

  return ddHpddx;
}

double BackgroundCosmology::get_OmegaB(double x) const{ 
  double a_minus3 = exp(-3*x);
  double H_over_H0 = H_of_x_over_H0(x);
  return OmegaB * a_minus3 / (H_over_H0 * H_over_H0);
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  double a_minus4 = exp(-4*x);
  double H_over_H0 = H_of_x_over_H0(x);
  return OmegaR * a_minus4 / (H_over_H0 * H_over_H0);
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  double a_minus4 = exp(-4*x);
  double H_over_H0 = H_of_x_over_H0(x);
  return OmegaNu * a_minus4 / (H_over_H0 * H_over_H0);
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  double a_minus3 = exp(-3*x);
  double H_over_H0 = H_of_x_over_H0(x);
  return OmegaCDM * a_minus3 / (H_over_H0 * H_over_H0);
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  double H_over_H0 = H_of_x_over_H0(x);
  return OmegaLambda / (H_over_H0 * H_over_H0);
}

double BackgroundCosmology::get_OmegaK(double x) const{ 
  double aH_over_H0 = exp(x)*H_of_x_over_H0(x);
  return OmegaK / (aH_over_H0 * aH_over_H0);
}
    
double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  double r = get_comoving_distance_of_x(x);
  double a = exp(x);
  return r / a;
}
double BackgroundCosmology::get_comoving_distance_of_x(double x) const{
  double eta0 = eta_of_x(0.0);
  double eta = eta_of_x(x);
  double chi = eta0 - eta;
  double H0chiOverC = H0*chi/Constants.c;
  double sqrt_OK = sqrt(abs(OmegaK));

  if (OmegaK < 0.0){
    return chi * sin(sqrt_OK*H0chiOverC) / (sqrt_OK * H0chiOverC);
  }
  else if (OmegaK > 0.0){
    return chi * sinh(sqrt_OK*H0chiOverC) / (sqrt_OK * H0chiOverC);
  }
  else{
    return chi;
  } 
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::t_of_x(double x) const{
  return t_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB(double x) const{ 
  if(x == 0.0) return TCMB;
  return TCMB * exp(-x); 
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK      << "\n";
  std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << "Age of the universe: " << t_of_x(0.0) * Constants.Gyr << " Gyr\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -15.0;
  const double x_max =  0.0;
  const int    n_pts =  100;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)        << " ";
    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << get_OmegaB(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x)      << " ";
    fp << get_OmegaNu(x)     << " ";
    fp << get_OmegaK(x)      << " ";
    fp << t_of_x(x)          << " ";
    fp << ddHpddx_of_x(x)     << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

