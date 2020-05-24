#
#
# 
library(tidyverse)
library(lubridate)
library(pomp)
theme_set(theme_bw())
options(stringsAsFactors=FALSE)
stopifnot(packageVersion("pomp")>="2.1")
set.seed(594709947L)

#
# Lee datos de repositorio covid mapache
#
corRaw<-read_csv('https://docs.google.com/spreadsheets/d/16-bnsDdmmgtSxdWbVMboIHo5FRuz76DBxsz_BbsEVWA/export?format=csv&id=16-bnsDdmmgtSxdWbVMboIHo5FRuz76DBxsz_BbsEVWA&gid=0')

cor <- corRaw %>% filter( osm_admin_level_4 =="CABA") %>% transmute(CABAdia=nue_casosconf_diff,fecha=dmy(fecha),dias =as.numeric( fecha - min(fecha)))

#
# Lee datos de repositorio Datos Abiertos Ministerio (1ra version) 
#
csv_fname <- "https://github.com/covid19UNGS/proyecciones/raw/master/Data/Covid19Casos.csv"
corRaw <- read_csv(csv_fname) 
cor <- corRaw %>% filter(clasificacion_resumen=="Confirmado",provincia_carga=="CABA",) %>% 
  transmute(fecha=ymd(fecha_fis),dias =as.numeric( fecha - min(fecha))) %>% group_by(dias) %>% summarise(CABADia=n())



#
# Lee archivo de datos de repositiorio covid19UNGS
#
csv_fname <- "https://raw.githubusercontent.com/covid19UNGS/proyecciones/master/Data/coronavirus_ar.csv"
cor <- read_csv(csv_fname) %>% dplyr::select(fecha:TDFdia)  %>% mutate(fecha=ymd(fecha), dias =as.numeric( fecha - min(fecha))) %>%
        mutate(importadosdia=importados-lag(importados))
cor$importadosdia[1] <- 1
cor <- cor %>% mutate(localesdia=casosdia - importadosdia, CABAdia=ifelse(is.na(CABAdia),0,CABAdia))

cor %>% select(dias,CABAdia) -> meas

meas %>%  ggplot(aes(x=dias,y=CABAdia))+
  geom_line()+
  geom_point()

#
# Formula SIR
#


sir_step <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-mu_IR*dt));
  S -= dN_SI;
  I += dN_SI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")

sir_init <- Csnippet("
  S = nearbyint(eta*N);
  I = 1;
  R = nearbyint((1-eta)*N);
  H = 0;
")

dmeas <- Csnippet("
  lik = dbinom(CABAdia,H,rho,give_log);
")

rmeas <- Csnippet("
  CABAdia = rbinom(H,rho);
")

meas %>%
  pomp(times = "dias",
       t0=0,
       rprocess=euler(sir_step,delta.t=1/7),
       rinit=sir_init,
       rmeasure=rmeas,
       dmeasure=dmeas,
       accumvars="H",
       statenames=c("S","I","R","H"),
       paramnames=c("Beta","mu_IR","N","eta","rho")
  ) -> measSIR

measSIR %>%
  simulate(params=c(Beta=0.48,mu_IR=0.0714,rho=0.5,eta=0.99,N=15000000),
           nsim=20,format="data.frame",include.data=TRUE) -> sims
# R0 = 0.48/0.074 -> 6.5

sims %>%
  ggplot(aes(x=dias,y=CABAdia,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color=FALSE) + scale_y_log10()

#' 


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
