#include"Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}

//====================================================
// Do all the solving
//====================================================

void Perturbations::solve(){

  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
  // compute_source_functions();
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations(){
  Utils::StartTiming("integrateperturbation");

  Vector k_array(n_k);
  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  Vector delta_CDM_array(n_x * n_k,0.0);
  Vector delta_b_array(n_x * n_k,0.0);
  Vector v_CDM_array(n_x * n_k,0.0);
  Vector v_b_array(n_x * n_k,0.0);
  Vector Phi_array(n_x * n_k,0.0);
  Vector Psi_array(n_x * n_k,0.0);
  Vector Pi_array(n_x * n_k,0.0);
  Vector2D Theta_array(Constants.n_ell_theta, Vector(n_x * n_k, 0.0));
  Vector2D Theta_p_array(Constants.n_ell_thetap, Vector(n_x * n_k, 0.0));
  Vector2D Nu_array(Constants.n_ell_neutrinos, Vector(n_x * n_k, 0.0));
  // We set the polarization multipoles and the l>1 photon multipoles to zero in the tight coupling regime, so we can just fill the array with zeros from the start
  
  // We use logarithmic spacing for the k-array
  Vector k_array_log = Utils::linspace(log(k_min), log(k_max), n_k);
  for(int ik = 0; ik < n_k; ik++){
    k_array[ik] = exp(k_array_log[ik]);
  }

  // Loop over all wavenumbers
  for(int ik = 0; ik < n_k; ik++){

    // Progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << (100*ik+100)/n_k << "% " << std::flush;
      if(ik == n_k-1) std::cout << std::endl;
    }

    // Current value of k
    double k = k_array[ik];

    // Find value to integrate to
    double x_end_tight = get_tight_coupling_time(k);

    // Find the index in the x-array where we switch from tight coupling to the full system
    auto end_tight = std::upper_bound(x_array.begin(), x_array.end(), x_end_tight);
    Vector x_array_tc(x_array.begin(), end_tight);
    Vector x_array_full(end_tight, x_array.end());
    x_array_full.insert(x_array_full.begin(), x_array_tc.back());
    

    auto y_tight_coupling_ini = set_ic(x_start, k);
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };

    // Integrating from x_start -> x_end_tight
    ODESolver ode_tc;
    ode_tc.solve(dydx_tight_coupling, x_array_tc, y_tight_coupling_ini);

    auto y_full_ini = set_ic_after_tight_coupling(ode_tc.get_final_data(), x_array_tc.back(), k);
    
    
    // The full ODE system
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };
    ODESolver ode_full;
    ode_full.solve(dydx_full, x_array_full, y_full_ini);


    Vector2D y_tc = ode_tc.get_data();
    Vector2D y_full = ode_full.get_data();

    const bool polarization = Constants.polarization;
    const bool neutrinos = Constants.neutrinos;
    const double c = Constants.c;
    const double H0 = cosmo->get_H0();
    const double OmegaR = cosmo->OmegaR;
    const double OmegaNu = cosmo->OmegaNu;

    for(int ix=0; ix < x_array_tc.size(); ix++){
      double Hp = cosmo->Hp_of_x(x_array_tc[ix]);
      double a = exp(x_array_tc[ix]);
      double ck_over_HpT = c*k/(Hp*rec->dtaudx_of_x(x_array_tc[ix]));
      delta_CDM_array[ix + n_x * ik] = y_tc[ix][Constants.ind_deltacdm_tc];
      delta_b_array[ix + n_x * ik] = y_tc[ix][Constants.ind_deltab_tc];
      v_CDM_array[ix + n_x * ik] = y_tc[ix][Constants.ind_vcdm_tc];
      v_b_array[ix + n_x * ik] = y_tc[ix][Constants.ind_vb_tc];
      Phi_array[ix + n_x * ik] = y_tc[ix][Constants.ind_Phi_tc];
      Theta_array[0][ix + n_x * ik] = y_tc[ix][Constants.ind_start_theta_tc];
      Theta_array[1][ix + n_x * ik] = y_tc[ix][Constants.ind_start_theta_tc + 1];
      if(polarization){
        Theta_array[2][ix + n_x * ik] = - 8. * ck_over_HpT * Theta_array[1][ix + n_x * ik] / 15.;
        Theta_p_array[0][ix + n_x * ik] = 5./4. * Theta_array[2][ix + n_x * ik];
        Theta_p_array[1][ix + n_x * ik] = - ck_over_HpT * Theta_array[2][ix + n_x * ik] / 4.;
        Theta_p_array[2][ix + n_x * ik] = Theta_array[2][ix + n_x * ik] / 4.;
        Pi_array[ix + n_x * ik] += Theta_p_array[0][ix + n_x * ik] + Theta_p_array[2][ix + n_x * ik];
        for(int l=3; l<Constants.n_ell_thetap_tc; l++){
          Theta_p_array[l][ix + n_x * ik] = -l/(2.*l+1.) * ck_over_HpT * Theta_p_array[l-1][ix + n_x * ik];
        }
      }
      else{
        Theta_array[2][ix + n_x * ik] = - 20. * ck_over_HpT * Theta_array[1][ix + n_x * ik] / 45.;
      }
      for(int l=3; l<Constants.n_ell_theta; l++){
        Theta_array[l][ix + n_x * ik] = -l/(2.*l+1.) * ck_over_HpT * Theta_array[l-1][ix + n_x * ik];
      }
      Psi_array[ix + n_x * ik] = -Phi_array[ix + n_x * ik] - 12.0 * pow(H0/(c*k*a), 2) * OmegaR * Theta_array[2][ix + n_x * ik];
      Pi_array[ix + n_x * ik] += Theta_array[2][ix + n_x * ik];

      if(neutrinos){
        for(int l=0; l<Constants.n_ell_neutrinos_tc; l++){
          Nu_array[l][ix + n_x * ik] = y_tc[ix][Constants.ind_start_nu_tc + l];
        }
        Psi_array[ix + n_x * ik] -= 12.0 * pow(H0/(c*k*a), 2) * OmegaNu * Nu_array[2][ix + n_x * ik];
      }
    }

    int start_idx = (int)x_array_tc.size() - 1; // index of overlap point
    for(int ix = 1; ix < (int)x_array_full.size(); ix++){
        // ix=0 is the overlap point already written by TC loop, skip it
        int global_ix = start_idx + ix;  // goes from start_idx+1 to n_x-1
        delta_CDM_array[global_ix + n_x * ik] = y_full[ix][Constants.ind_deltacdm];
        delta_b_array[global_ix + n_x * ik]   = y_full[ix][Constants.ind_deltab];
        v_CDM_array[global_ix + n_x * ik]     = y_full[ix][Constants.ind_vcdm];
        v_b_array[global_ix + n_x * ik]       = y_full[ix][Constants.ind_vb];
        Phi_array[global_ix + n_x * ik]       = y_full[ix][Constants.ind_Phi];
        for(int l=0; l<Constants.n_ell_theta; l++)
            Theta_array[l][global_ix + n_x * ik] = y_full[ix][Constants.ind_start_theta + l];
        double a = exp(x_array_full[ix]);
        Psi_array[global_ix + n_x * ik] = -Phi_array[global_ix + n_x * ik]
            - 12.0 * pow(H0/(c*k*a), 2) * OmegaR * Theta_array[2][global_ix + n_x * ik];
        Pi_array[global_ix + n_x * ik] = Theta_array[2][global_ix + n_x * ik];
        if(polarization){
            for(int l=0; l<Constants.n_ell_thetap; l++)
                Theta_p_array[l][global_ix + n_x * ik] = y_full[ix][Constants.ind_start_thetap + l];
            Pi_array[global_ix + n_x * ik] +=
                Theta_p_array[0][global_ix + n_x * ik] + Theta_p_array[2][global_ix + n_x * ik];
        }
        if(neutrinos){
            for(int l=0; l<Constants.n_ell_neutrinos; l++)
                Nu_array[l][global_ix + n_x * ik] = y_full[ix][Constants.ind_start_nu + l];
            Psi_array[global_ix + n_x * ik] -=
                12.0 * pow(H0/(c*k*a), 2) * OmegaNu * Nu_array[2][global_ix + n_x * ik];
        }
    }

    
  }
  Utils::EndTiming("integrateperturbation");

  delta_cdm_spline.create(x_array, k_array, delta_CDM_array, "delta_cdm_spline");
  delta_b_spline.create(x_array, k_array, delta_b_array, "delta_b_spline");
  v_cdm_spline.create(x_array, k_array, v_CDM_array, "v_cdm_spline");
  v_b_spline.create(x_array, k_array, v_b_array, "v_b_spline");
  Phi_spline.create(x_array, k_array, Phi_array, "Phi_spline");
  Psi_spline.create(x_array, k_array, Psi_array, "Psi_spline");
  Pi_spline.create(x_array, k_array, Pi_array, "Pi_spline");

  Theta_spline = std::vector<Spline2D>(Constants.n_ell_theta);
  Theta_p_spline = std::vector<Spline2D>(Constants.n_ell_thetap);
  Nu_spline = std::vector<Spline2D>(Constants.n_ell_neutrinos);
  for(int l=0; l<Constants.n_ell_theta; l++){
    Theta_spline[l].create(x_array, k_array, Theta_array[l], "Theta_spline_l"+std::to_string(l));
  }
  if(Constants.polarization){
    for(int l=0; l<Constants.n_ell_thetap; l++){
      Theta_p_spline[l].create(x_array, k_array, Theta_p_array[l], "Theta_p_spline_l"+std::to_string(l));
    }
  }
  if(Constants.neutrinos){
    for(int l=0; l<Constants.n_ell_neutrinos; l++){
      Nu_spline[l].create(x_array, k_array, Nu_array[l], "Nu_spline_l"+std::to_string(l));
    }
  }


}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b          =  y_tc[Constants.ind_vb_tc];
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  double *Theta        = &y_tc[Constants.ind_start_theta_tc];
  double *Nu           = &y_tc[Constants.ind_start_nu_tc];

  double c = Constants.c;
  double OmegaNu = cosmo->OmegaNu;
  double OmegaR = cosmo->OmegaR;
  double f_nu = OmegaNu / (OmegaR + OmegaNu);
  double Psi = -1/(1.5 + 0.4*f_nu);

  double a = exp(x);
  double H0 = cosmo->get_H0();
  double Hp = cosmo->Hp_of_x(x);
  double dtaudx = rec->dtaudx_of_x(x);

  Phi = -(1 + 0.4*f_nu) * Psi;
  delta_cdm = -1.5*Psi;
  delta_b = delta_cdm;
  
  v_cdm = -c*k / (2.0 * Hp) * Psi;
  v_b = v_cdm;

  Theta[0] = -0.5*Psi;
  Theta[1] = c*k/(6.0*Hp)*Psi;

  if(neutrinos){
    Nu[0] = -0.5*Psi;
    Nu[1] = c*k/(6.0*Hp)*Psi;
    Nu[2] = -c*c * k*k * a*a * (Psi + Phi) / (12.0 * H0*H0 * OmegaNu);

    for(int l=3; l<n_ell_neutrinos_tc; l++){
      Nu[l] = c*k / (2.*l+1) / Hp * Nu[l-1];
    }
  }

  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);

  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double *Theta_tc        = &y_tc[Constants.ind_start_theta_tc];
  const double *Nu_tc           = &y_tc[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm];
  double &delta_b         =  y[Constants.ind_deltab];
  double &v_cdm           =  y[Constants.ind_vcdm];
  double &v_b             =  y[Constants.ind_vb];
  double &Phi             =  y[Constants.ind_Phi];
  double *Theta           = &y[Constants.ind_start_theta];
  double *Theta_p         = &y[Constants.ind_start_thetap];
  double *Nu              = &y[Constants.ind_start_nu];

  double c = Constants.c;
  double Hp = cosmo->Hp_of_x(x);
  double dtaudx = rec->dtaudx_of_x(x);

  Phi = Phi_tc;
  delta_cdm = delta_cdm_tc;
  delta_b = delta_b_tc;
  v_cdm = v_cdm_tc;
  v_b = v_b_tc;

  Theta[0] = Theta_tc[0];
  Theta[1] = Theta_tc[1];
  
  if(!polarization){
    Theta[2] = -20.0 * c*k /(45.0 * Hp * dtaudx) * Theta_tc[1];
  }
  if(polarization){
    Theta[2] = -8.0 * c*k /(15.0 * Hp * dtaudx) * Theta_tc[1];
  }
  for(int l=3; l<n_ell_theta; l++){
    Theta[l] = -l/(2.*l+1.) * c*k / (Hp*dtaudx) * Theta[l-1];
  }


  // SET: Photon polarization perturbations (Theta_p_ell)
  if(polarization){
    Theta_p[0] = 5./4. * Theta[2];
    Theta_p[1] = -c*k / (4. * Hp * dtaudx) * Theta[2];
    Theta_p[2] = Theta[2] / 4.;
    for(int l=3; l<n_ell_thetap; l++){
      Theta_p[l] = -l/(2.*l+1.) * c*k / (Hp*dtaudx) * Theta_p[l-1];
    }
  }

  // SET: Neutrino perturbations (N_ell)
  if(neutrinos){
    for(int l=0; l<n_ell_neutrinos_tc; l++){
      Nu[l] = Nu_tc[l];
    }
    for(int l=n_ell_neutrinos_tc; l<n_ell_neutrinos; l++){
        Nu[l] = 0.0; // Just in case we use more neutrinos in the full regime than in the tight coupling regime, we set the higher multipoles to zero
    }
  }
  // std::cout << "Initial conditions for k = " << k << ", x = " << x << ", H0 " << cosmo->get_H0() << ":\n";
  // std::cout << "  delta_cdm: " << delta_cdm << "\n";
  // std::cout << "  delta_b: " << delta_b << "\n";
  // std::cout << "  v_cdm: " << v_cdm << "\n";
  // std::cout << "  v_b: " << v_b << "\n";
  // std::cout << "  Phi: " << Phi << "\n";
  // for(int l = 0; l < n_ell_theta; l++){
  //   std::cout << "  Theta[" << l << "]: " << Theta[l] << "\n";
  // }
  // if(polarization){
  //   for(int l = 0; l < n_ell_thetap; l++){
  //     std::cout << "  Theta_p[" << l << "]: " << Theta_p[l] << "\n";
  //   }
  // }
  // if(neutrinos){
  //   for(int l = 0; l < n_ell_neutrinos; l++){
  //     std::cout << "  Nu[" << l << "]: " << Nu[l] << "\n";
  //   }
  // }

  return y;
}

//====================================================
// The time when tight coupling end
//====================================================

double Perturbations::get_tight_coupling_time(const double k) const{
  double ck = Constants.c * k;
  double dx = (x_end - x_start) / (n_x - 1.0);  // match linspace grid exactly
  double x = x_start;
  
  while(x < x_end){
    double dtaudx = rec->dtaudx_of_x(x);
    double Hp = cosmo->Hp_of_x(x);
    if(std::abs(dtaudx) <= 10.0 * std::max(1.0, ck / Hp))
      break;
    x += dx;
  }

  if(x > -7.0) x = -7.0;
  // Snap x to the nearest grid point
  int idx = (int)std::round((x - x_start) / dx);
  return x_start + idx * dx;
}

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");

  //=============================================================================
  // TODO: Make the x and k arrays to evaluate over and use to make the splines
  //=============================================================================
  // ...
  // ...
  Vector k_array;
  Vector x_array;

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  Vector SE_array(k_array.size() * x_array.size());

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;

      //=============================================================================
      // TODO: Compute the source functions
      //=============================================================================
      // Fetch all the things we need...
      // const double Hp       = cosmo->Hp_of_x(x);
      // const double tau      = rec->tau_of_x(x);
      // ...
      // ...

      // Temperatur source
      ST_array[index] = 0.0;

      // Polarization source
      if(Constants.polarization){
        SE_array[index] = 0.0;
      }
    }
  }

  // Spline the source functions
  ST_spline.create (x_array, k_array, ST_array, "Source_Temp_x_k");
  if(Constants.polarization){
    SE_spline.create (x_array, k_array, SE_array, "Source_Pol_x_k");
  }

  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double *Theta           = &y[Constants.ind_start_theta_tc];
  const double *Nu              = &y[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double *dThetadx        = &dydx[Constants.ind_start_theta_tc];
  double *dNudx           = &dydx[Constants.ind_start_nu_tc];

  // cosmological parameters and variables
  double ck = Constants.c * k;
  double H0 = cosmo->get_H0();
  double H0_square = H0*H0;
  double Hp = cosmo->Hp_of_x(x);
  double dHpdx = cosmo->dHpdx_of_x(x);
  double keta = k*cosmo->eta_of_x(x);
  double keta_inv = (keta > 1e-8) ? 1.0 / keta : 0.0; // To avoid division by zero when k*eta is very small
  double ck_over_Hp = ck / Hp;
  double OmegaNu = cosmo->OmegaNu;
  double OmegaR = cosmo->OmegaR;
  double OmegaB = cosmo->OmegaB;
  double OmegaCDM = cosmo->OmegaCDM;
  double a_inv = exp(-x);
  double a2_inv = a_inv * a_inv;
  double R = 4.*OmegaR*a_inv / (3.*OmegaB);
  double H0_over_cka = H0 * a_inv / (Constants.c*k);
  double H0_over_cka_square = H0_over_cka * H0_over_cka;
  
  // Recombination variables
  double dtaudx = rec->dtaudx_of_x(x);
  double ddtauddx = rec->ddtauddx_of_x(x);

  // Neutrino multipoles
  
    double Theta2;
    if(polarization){
      Theta2 = -8.0 * ck_over_Hp /(15.0 * dtaudx) * Theta[1];
    }
    else{
      Theta2 = -20.0 * ck_over_Hp /(45.0 * dtaudx) * Theta[1];
    }

  double Psi;
  if(neutrinos){
    Psi = -Phi - 12.* H0_over_cka_square * (OmegaR*Theta2 + OmegaNu*Nu[2]);
    dPhidx = Psi - ck_over_Hp * ck_over_Hp / 3. * Phi + H0_square / (2*Hp*Hp) * (OmegaCDM*a_inv*delta_cdm + OmegaB*a_inv*delta_b + 4.*OmegaR*a2_inv*Theta[0] + 4.*OmegaNu*a2_inv*Nu[0]);
    dNudx[0] = -ck_over_Hp * Nu[1] - dPhidx;
    dNudx[1] = ck_over_Hp/3. * (Nu[0] - 2*Nu[2] + Psi);
    int l_max = n_ell_neutrinos_tc-1;
    for(int l=2; l<l_max; l++){
      dNudx[l] = ck_over_Hp / (2.*l+1) * (l*Nu[l-1] - (l+1)*Nu[l+1]);
    }
    dNudx[l_max] = ck_over_Hp * (Nu[l_max-1] - (l_max+1)*keta_inv*Nu[l_max]);
  }
  else{
    Psi = -Phi - 12.* H0_over_cka_square * OmegaR*Theta2;
    dPhidx = Psi - ck_over_Hp * ck_over_Hp / 3. * Phi + H0_square / (2*Hp*Hp) * (OmegaCDM*a_inv*delta_cdm + OmegaB*a_inv*delta_b + 4.*OmegaR*a2_inv*Theta[0]);
  }

  // Scalar quantities
  ddelta_cdmdx = ck_over_Hp * v_cdm - 3.0 * dPhidx;
  dv_cdmdx = -v_cdm - ck_over_Hp * Psi;
  ddelta_bdx = ck_over_Hp * v_b - 3.0 * dPhidx;
  dThetadx[0] = -ck_over_Hp * Theta[1] - dPhidx;
  double q = (-((1-R)*dtaudx + (1+R)*ddtauddx)*(3*Theta[1] + v_b) - ck_over_Hp*Psi + (1 - dHpdx/Hp) * ck_over_Hp * (2*Theta2 - Theta[0]) - ck_over_Hp * (-ck_over_Hp * Theta[1] - dPhidx)) / ((1+R)*dtaudx + dHpdx/Hp - 1);
  dv_bdx = (-v_b - ck_over_Hp * Psi + R*(q + ck_over_Hp*(2*Theta2 - Theta[0]) - ck_over_Hp * Psi)) / (1+R);
  dThetadx[1] = (q - dv_bdx) / 3.;

  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){

  // Index and number of the different quantities
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double *Theta           = &y[Constants.ind_start_theta];
  const double *Theta_p         = &y[Constants.ind_start_thetap];
  const double *Nu              = &y[Constants.ind_start_nu];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double *dThetadx        = &dydx[Constants.ind_start_theta];
  double *dTheta_pdx      = &dydx[Constants.ind_start_thetap];
  double *dNudx           = &dydx[Constants.ind_start_nu];

  // cosmological parameters and variables
  double ck = Constants.c * k;
  double H0 = cosmo->get_H0();
  double H0_square = H0*H0;
  double Hp = cosmo->Hp_of_x(x);
  double dHpdx = cosmo->dHpdx_of_x(x);
  double keta = k*cosmo->eta_of_x(x);
  double keta_inv = (keta > 1e-8) ? 1.0 / keta : 0.0; // To avoid division by zero when k*eta is very small
  double ck_over_Hp = ck / Hp;
  double OmegaNu = cosmo->OmegaNu;
  double OmegaR = cosmo->OmegaR;
  double OmegaB = cosmo->OmegaB;
  double OmegaCDM = cosmo->OmegaCDM;
  double a_inv = exp(-x);
  double a2_inv = a_inv * a_inv;
  double R = 4.*OmegaR*a_inv / (3.*OmegaB);
  double H0_over_cka = H0 * a_inv / (Constants.c*k);
  double H0_over_cka_square = H0_over_cka * H0_over_cka;
  
  // Recombination variables
  double dtaudx = rec->dtaudx_of_x(x);

  double Psi;
  if(neutrinos){
    Psi = -Phi - 12.* H0_over_cka_square * (OmegaR*Theta[2] + OmegaNu*Nu[2]);
    dPhidx = Psi - ck_over_Hp * ck_over_Hp / 3. * Phi + H0_square / (2*Hp*Hp) * (OmegaCDM*a_inv*delta_cdm + OmegaB*a_inv*delta_b + 4.*OmegaR*a2_inv*Theta[0] + 4.*OmegaNu*a2_inv*Nu[0]);
    dNudx[0] = -ck_over_Hp * Nu[1] - dPhidx;
    dNudx[1] = ck_over_Hp/3. * (Nu[0] - 2*Nu[2] + Psi);
    int l_max = n_ell_neutrinos-1;
    for(int l=2; l<l_max; l++){
      dNudx[l] = ck_over_Hp / (2.*l+1) * (l*Nu[l-1] - (l+1)*Nu[l+1]);
    }
    dNudx[l_max] = ck_over_Hp * (Nu[l_max-1] - (l_max+1)*keta_inv*Nu[l_max]);
  }
  else{
    Psi = -Phi - 12.* H0_over_cka_square * (OmegaR*Theta[2]);
    dPhidx = Psi - ck_over_Hp * ck_over_Hp / 3. * Phi + H0_square / (2*Hp*Hp) * (OmegaCDM*a_inv*delta_cdm + OmegaB*a_inv*delta_b + 4.*OmegaR*a2_inv*Theta[0]);
  }
  // if((k <= k_min + 1e-3)){
  // std::cout << "k: " << k << ", x: " << x << ", Phi: " << Phi << ", Psi: " << Psi << ", dPhidx: " << dPhidx << ", term 1: " << ck_over_Hp * ck_over_Hp / 3. * Phi << ", term 2: " << H0_square / (2*Hp*Hp) * (OmegaCDM*a_inv*delta_cdm + OmegaB*a_inv*delta_b + 4.*OmegaR*a2_inv*Theta[0]) << std::endl;
  // }
  // Scalar quantities
  ddelta_cdmdx = ck_over_Hp * v_cdm - 3.0 * dPhidx;
  dv_cdmdx = -v_cdm - ck_over_Hp * Psi;
  ddelta_bdx = ck_over_Hp * v_b - 3.0 * dPhidx;
  dv_bdx = - v_b - ck_over_Hp * Psi + dtaudx*R*(3*Theta[1] + v_b);
  double slip = 3.0 * Theta[1] + v_b;
  
  int l_max_theta = n_ell_theta-1;
  // Photon multipoles
  double Pi;
  if(polarization){
    Pi = Theta[2] + Theta_p[0] + Theta_p[2];
  }
  else{
    Pi = Theta[2];
  }
  dThetadx[0] = -ck_over_Hp * Theta[1] - dPhidx;
  dThetadx[1] = ck_over_Hp/3. * (Theta[0] - 2*Theta[2] + Psi) + dtaudx * (Theta[1] + v_b/3.);
  dThetadx[2] = ck_over_Hp / 5. * (2*Theta[1] - 3*Theta[3]) + dtaudx * (Theta[2] - Pi/10.);
  for(int l=3; l<l_max_theta; l++){
    dThetadx[l] = ck_over_Hp / (2.*l+1.) * (l*Theta[l-1] - (l+1)*Theta[l+1]) + dtaudx * Theta[l];
  }
  dThetadx[l_max_theta] = ck_over_Hp * (Theta[l_max_theta-1] - (l_max_theta+1)*keta_inv*Theta[l_max_theta]) + dtaudx * Theta[l_max_theta];
  
  int l_max_thetap = n_ell_thetap-1;
  // Photon polarization multipoles
  if(polarization){
    dTheta_pdx[0] = - ck_over_Hp*Theta_p[1] + dtaudx * (Theta_p[0] - Pi/2.);
    dTheta_pdx[1] = ck_over_Hp/3. * (Theta_p[0] - 2*Theta_p[2]) + dtaudx * Theta_p[1];
    dTheta_pdx[2] = ck_over_Hp/5. * (2*Theta_p[1] - 3*Theta_p[3]) + dtaudx * (Theta_p[2] - Pi/10.);
    for(int l=3; l<l_max_thetap; l++){
      dTheta_pdx[l] = ck_over_Hp / (2.*l+1.) * (l*Theta_p[l-1] - (l+1)*Theta_p[l+1]) + dtaudx * Theta_p[l];
    }
    dTheta_pdx[l_max_thetap] = ck_over_Hp * (Theta_p[l_max_thetap-1] - (l_max_thetap+1)*keta_inv*Theta_p[l_max_thetap]) + dtaudx * Theta_p[l_max_thetap];
  }
  
  // std::cout << "k: " << k << ", Psi: " << Psi << ", Theta[2]: " << Theta[2] << ", Nu[2]: " << Nu[2] << ", delta_cdm: " << delta_cdm << ", v_cdm: " << v_cdm << std::endl;
  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================

double Perturbations::get_delta_cdm(const double x, const double k) const{
  return delta_cdm_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_v_cdm(const double x, const double k) const{
  return v_cdm_spline(x,k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Pi_spline(x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Source_E(const double x, const double k) const{
  return SE_spline(x,k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}
double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
  return Theta_p_spline[ell](x,k);
}
double Perturbations::get_Nu(const double x, const double k, const int ell) const{
  return Nu_spline[ell](x,k);
}

//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start, x_end, npts);
  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " ";
    fp << get_Theta(x,k,0)   << " ";
    fp << get_Theta(x,k,1)   << " ";
    fp << get_Theta(x,k,2)   << " ";
    fp << get_Phi(x,k)       << " ";
    fp << get_Psi(x,k)       << " ";
    fp << get_Pi(x,k)        << " ";
    fp << get_delta_cdm(x,k)   << " ";
    fp << get_delta_b(x,k)     << " ";
    fp << get_v_cdm(x,k)       << " ";
    fp << get_v_b(x,k)         << " ";
    fp << get_Nu(x,k,0)        << " ";
    fp << get_Nu(x,k,1)        << " ";
    fp << get_Nu(x,k,2)        << " ";
    fp << get_Theta_p(x,k,0)   << " ";
    fp << get_Theta_p(x,k,1)   << " ";
    fp << get_Theta_p(x,k,2)   << " ";
    // fp << get_Source_T(x,k)  << " ";
    // fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    // fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    // fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

