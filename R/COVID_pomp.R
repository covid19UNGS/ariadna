# for fitting E0 added E_init parameter that would need to be discarded
# define the changes for a time step forward
sir_step <- Csnippet("
                     // adjust betat for social distancing interventions
                     double betat;
                     betat = beta0; // everything else is no intervention
                     
                     // adjust contact rates if isolation of symptomatic cases is in place
                     double iso_m = 0;
                     double iso_s = 0;

                    
                     // if import rate is above zero, draw importations, assuming they are perfectly balanced with departures of susceptible individuals
                     double import = 0;
                     if(import_rate > 0){
                      import = fmin(rpois(import_rate*dt), S);
                     }
                     // tracking of total imported, removing them them from susceptibles
                     import_total += import;
                     S -= import;
                    
                     // calculate transition numbers
                     double dSE = rbinom(S, 1-exp(-betat*(Ca*Ia/N + Cp*Ip/N + iso_m*Cm*Im/N + iso_s*Cs*Is/N)*dt)); 
                     
                     // gamma = 1/periodo de latencia
                     // alpha = proporcion de asyntomaticos
                     double rateE[2];
                     double dE_all[2];
                     rateE[0] = alpha*gamma; // going to asymtomatic
                     rateE[1] = (1-alpha)*gamma; // going to presymptomatic
                     reulermultinom(2, E, rateE, dt, &dE_all);
                     double dEIa = dE_all[0];
                     double dEIp = dE_all[1];
                     double dIaR = rbinom(Ia, 1 - exp(-lambda_a*dt));
                     
                     // Desde Ip va a Ih=se van al hospital, Ic=Se van a la casa
                     //
                     double rateIp[2];
                     double dIp_all[2];
                     rateIp[0] = mu*lambda_p;      // van a la casa
                     rateIp[1] = (1-mu)*lambda_p;  // van al hospital
                     reulermultinom(2, Ip, rateIp, dt, &dIp_all);
                     double dIpIc = dIp_all[0];
                     double dIpIh = dIp_all[1];
                     
                     // De Ic (casa) a recuperados, 1/lambda_m periodo de recuperacion de Ic 
                     //
                     double dIcR = rbinom(Ic, 1 - exp(-lambda_m*dt));

                     // De Los hospitalizados el parametro vu define si va ITU
                     //    variables = It (terapia intensiva) Iu va a cama comun
                     // definimos parametro lambda_h 
                     //double rateIp[2];
                     //double dIp_all[2];
                     rateIp[0] = vu*lambda_h;     // va a terapia intensiva
                     rateIp[1] = (1-vu)*lambda_h; // va a cama comun
                     reulermultinom(2, Ih, rateIp, dt, &dIp_all);
                     double dIpIt = dIp_all[0];
                     double dIpIu = dIp_all[1];

                     // De terapia intensiva a Death It --> D
                     // delta proporcion de It que mueren
                     // lambda_t tasa de salida de Iu 1/It es el tiempo de internacion en UTI
                     //
                     double rateIs[2];
                     double dIs_all[2];
                     rateIs[0] = delta*lambda_t;     // hospitalized en ICU ultimately going to death
                     rateIs[1] = (1-delta)*lambda_t; // hospitalized en ICU ultimately going to recovered
                     reulermultinom(2, It, rateIs, dt, &dIs_all);
                     double dItHd = dIs_all[0];
                     double dItHr = dIs_all[1];
                     
                     // De cama comun a Death Iu --> D o no
                     // deltu proporcion de Iu que mueren
                     // lambda_u tasa de salida de Iu 1/Iu es el tiempo de internacion en UTI
                     //
                     double rateIs[2];
                     double dIs_all[2];
                     rateIs[0] = deltu*lambda_u;     // hospitalized en ICU ultimately going to death
                     rateIs[1] = (1-deltu)*lambda_u; // hospitalized en ICU ultimately going to recovered
                     reulermultinom(2, Iu, rateIs, dt, &dIs_all);
                     double dIuHd = dIs_all[0];
                     double dIuHr = dIs_all[1];


                     // update the compartments
                     S  -= dSE; // susceptible 
                     E  += dSE - dEIa - dEIp + import; // exposed
                     Ia += dEIa - dIaR; // infectious and asymptomatic
                     Ip += dEIp - dIpIs - dIpIm; // infectious and pre-symptomatic
                     Is += dIpIs - dIsHd - dIsHr; // infectious and severe symptoms (that will be hospitalized)
                     Im += dIpIm - dImR; // infectious and minor symptoms
                     I   = Ia + Ip + Im + Is; // total number of infected
                     I_new_sympt += dIpIs + dIpIm; // total number of newly symptomatic
                     Hr += dIsHr - dHrR; // hospitalized that will recover
                     Hd += dIsHd - dHdD; // hospitalizations that will die
                     H   = Hr + Hd; // total hospitalizations
                     R  += dHrR + dImR + dIaR; // recovered
                     D  += dHdD; // fatalities
                     D_new += dHdD; // daily fatalities
                     H_new += dIsHr + dIsHd; // daily new hospitalizations
                     if(intervention == 3 & H >= thresh_H_start) thresh_crossed = 1;
                     else if(intervention == 3 & thresh_crossed == 1 & H < thresh_H_end) thresh_crossed = 0;
                     ")

# define the initial set up, currently, every is susceptible except the exposed people
sir_init <- Csnippet("
                     S = N-E0;
                     E = E0;
                     Ia = 0;
                     Ip = 0;
                     Is = 0;
                     Im = 0;
                     I = 0;
                     I_new_sympt = 0;
                     Hr = 0;
                     Hd = 0;
                     H = Hd + Hr;
                     R = 0;
                     D = 0;
                     D_new = 0;
                     H_new = 0;
                     thresh_crossed = 0;
                     import_total = 0;
                     ")

# define random simulator of measurement
rmeas_deaths <- Csnippet("double tol = 1e-16;
                   deaths = rpois(D_new + tol);
                  ")
# define evaluation of model prob density function
dmeas_deaths <- Csnippet("double tol = 1e-16;
                   lik = dpois(deaths, D_new + tol, give_log);
                  ")
# define random simulator of measurement
rmeas_hosp <- Csnippet("double tol = 1e-16;
                   hosp = rpois(H + tol);
                  ")
# define evaluation of model prob density function
dmeas_hosp <- Csnippet("double tol = 1e-16;
                   lik = dpois(hosp, H + tol, give_log);
                  ")

# parameters to transform
if (fitting) {
  if (fit_to_sip) {
par_trans <- parameter_trans(log = c("beta0", "import_rate"),
                            logit = c("soc_dist_level_sip"))

param_names <- c(
   "beta0"
  , "Ca", "Cp", "Cs", "Cm"
  , "alpha"
  , "mu"
  , "delta"
  , "gamma"
  , "lambda_a", "lambda_s", "lambda_m", "lambda_p"
  , "rho_r"
  , "rho_d"
  , "N"
  , "E0"
  , "soc_dist_level_sip"
  , "import_rate"
)
  } else {

par_trans <- parameter_trans(log = c("beta0", "import_rate"))

param_names <- c(
   "beta0"
  , "Ca", "Cp", "Cs", "Cm"
  , "alpha"
  , "mu"
  , "delta"
  , "gamma"
  , "lambda_a", "lambda_s", "lambda_m", "lambda_p"
  , "rho_r"
  , "rho_d"
  , "N"
  , "E0"
  , "import_rate"
)
      
}

} else {
par_trans <- parameter_trans(log = c("beta0"))  

param_names <- c(
   "beta0"
  , "Ca", "Cp", "Cs", "Cm"
  , "alpha"
  , "mu"
  , "delta"
  , "gamma"
  , "lambda_a", "lambda_s", "lambda_m", "lambda_p"
  , "rho_r"
  , "rho_d"
  , "N"
  , "E0"
  , "import_rate"
)

}

# variables that should be zeroed after each obs
accum_names = c("D_new", "H_new", "I_new_sympt")

# state variables
state_names = c(
    "S" , "E" , "Ia"
  , "Ip", "Is", "Im"
  , "I" , "I_new_sympt"
  , "H" , "Hr", "Hd"
  , "R" , "D" 
  , "D_new", "H_new" 
  , "thresh_crossed"
  , "import_total"
)


## R0 here just based on the simple transmission rate / recovery rate (weighted by the probability of going into different classes)
covid_R0 <- function (beta0est, fixed_params, sd_strength, prop_S) {
## transmission rate
 R <-   beta0est * prop_S * sd_strength * 
    (                
## proportion * time in asymptomatic
      fixed_params["alpha"] * fixed_params["Ca"] * (1/fixed_params["lambda_a"]) +                  
## proportion * time in mildly symptomatic
      (1 - fixed_params["alpha"]) * fixed_params["mu"] * ((1/fixed_params["lambda_p"]) + (1/fixed_params["lambda_m"])) +    
## proportion * time in severely symptomatic
      (1 - fixed_params["alpha"]) * (1 - fixed_params["mu"]) * ((1/fixed_params["lambda_p"]) + (1/fixed_params["lambda_s"]))      
      )
 
 unlist(R)
}

