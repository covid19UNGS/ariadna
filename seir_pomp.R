#
#
# 
require(tidyverse)
require(pomp)

csv_fname <- "https://raw.githubusercontent.com/covid19UNGS/proyecciones/master/Data/coronavirus_ar.csv"
cor <- read_csv(csv_fname) %>% dplyr::select(fecha:TDFdia)  %>% mutate(fecha=ymd(fecha), dias =as.numeric( fecha - min(fecha))) %>%
        mutate(importadosdia=importados-lag(importados))
cor$importadosdia[1] <- 1
cor <- cor %>% mutate(localesdia=casosdia - importadosdia, CABAdia=ifelse(is.na(CABAdia),0,CABAdia))


seir_step <- Csnippet("
  double dN_SE = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_EI = rbinom(E,1-exp(-mu_EI*dt));
  double dN_IR = rbinom(I,1-exp(-mu_IR*dt));
  S -= dN_SE;
  E += dN_SE - dN_EI;
  I += dN_EI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")


seir_init <- Csnippet("
  S = nearbyint(eta*N);
  E = 0;
  I = 1;
  R = nearbyint((1-eta)*N);
  H = 0;
")

measSIR %>%
  pomp(
    rprocess=euler(seir_step,delta.t=1/7),
    rinit=seir_init,
    paramnames=c("N","Beta","mu_EI","mu_IR","rho","eta"),
    statenames=c("S","E","I","R","H")
  ) -> measSEIR
