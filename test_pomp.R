
#
# Simulamos el modelo de Ricker con variabilidad lognormal 
#
library(pomp)

simulate(t0=0, times=1:20,
         params=c(r=1.2,K=200,sigma=0.1,N_0=50),
         rinit=function (N_0, ...) {
           c(N=N_0)
         },
         rprocess=discrete_time(
           function (N, r, K, sigma, ...) {
             eps <- rnorm(n=1,mean=0,sd=sigma)
             c(N=r*N*exp(1-N/K+eps))
           },
           delta.t=1
         )
) -> sim1


ggplot(data=as.data.frame(sim1),aes(x=time,y=N))+
  geom_line()


simulate(t0=0, times=1:20,
         params=c(r=1.2,K=200,sigma=0.1,N_0=50,b=0.05),
         rinit=function (N_0, ...) {
           c(N=N_0)
         },
         rprocess=discrete_time(
           function (N, r, K, sigma, ...) {
             eps <- rnorm(n=1,mean=0,sd=sigma)
             c(N=r*N*exp(1-N/K+eps))
           },
           delta.t=1
         ),
         rmeasure=function (N, b, ...) {
           c(Y=rpois(n=1,lambda=b*N))
         }
) -> sim2

plot(sim2)
ggplot(data=gather(
  as(sim2,"data.frame"),
  variable,value,-time),
  aes(x=time,y=value,color=variable))+
  geom_line()

simulate(sim2,nsim=20) -> sims

ggplot(data=gather(
  as.data.frame(sims),
  variable,value,Y,N
),
aes(x=time,y=value,color=variable,
    group=interaction(.id,variable)))+
  geom_line()+
  facet_grid(variable~.,scales="free_y")+
  labs(y="",color="")

parus %>%
  ggplot(aes(x=year,y=pop))+
  geom_line()+geom_point()+
  expand_limits(y=0)

#
# Parus Ricker and Logistic continuous
#

parus %>%
  pomp(
    times="year", t0=1960,
    rinit=function (N_0, ...) {
      c(N=N_0)
    },
    rprocess=discrete_time(
      function(N, r, K, sigma, ...) {
        eps <- rnorm(n=1,mean=0,sd=sigma)
        c(N=r*N*exp(1-N/K+eps))
      },
      delta.t=1
    ),
    rmeasure=function (N, b, ...) {
      c(pop=rpois(n=1,lambda=b*N))
    }
  ) -> rick

#
#
#
vpstep <- function (N, r, K, sigma, delta.t, ...) {
  dW <- rnorm(n=1,mean=0,sd=sqrt(delta.t))
  c(N = N + r*N*(1-N/K)*delta.t + sigma*N*dW)
}

rick %>% pomp(rprocess=euler(vpstep,delta.t=1/365)) -> vp

vp %>%
  simulate(
    params=c(r=0.5,K=2000,sigma=0.1,b=0.1,N_0=2000),
    format="data.frame", include.data=TRUE, nsim=5) %>%
  mutate(ds=case_when(.id=="data"~"data",TRUE~"simulation")) %>%
  ggplot(aes(x=year,y=pop,group=.id,color=ds))+
  geom_line()+
  labs(color="")


Csnippet("
  pop = rpois(b*N);  
  ") -> rmeas

Csnippet("
  N = N_0;
  ") -> rinit

Csnippet("
  double eps = rnorm(0,sigma);
  N = r*N*exp(1-N/K+eps);
") -> rickstepC

Csnippet("
  double dW = rnorm(0,sqrt(dt));
  N += r*N*(1-N/K)*dt+sigma*N*dW;
") -> vpstepC

parus %>%
  pomp(
    times="year", t0=1960,
    rinit=rinit,
    rmeasure=rmeas,
    rprocess=discrete_time(rickstepC,delta.t=1),
    statenames="N",
    paramnames=c("r","K","sigma","b","N_0")
  ) -> rickC

parus %>%
  pomp(
    times="year", t0=1960,
    rinit=rinit,
    rmeasure=rmeas,
    rprocess=euler(vpstepC,delta.t=1/365),
    statenames="N",
    paramnames=c("r","K","sigma","b","N_0")
  ) -> vpC

Csnippet("
  lik = dpois(pop,b*N,give_log);
") -> dmeas

rickC %>%
  pfilter(Np=1000,
          params=c(r=1.2,K=2000,sigma=0.3,N_0=1600,b=0.1),
          dmeasure=dmeas,
          paramnames="b",statenames="N") -> pfrick

logLik(pfrick)


sobolDesign(
  lower=c(r=0,K=100,sigma=0,N_0=150,b=1),
  upper=c(r=5,K=600,sigma=2,N_0=150,b=1),
  nseq=100
) -> guesses
guesses

vpC %>%
  mif2(
    params=guesses[1,],
    Np=1000,
    Nmif=20,
    dmeasure=dmeas,
    partrans=parameter_trans(log=c("r","K","sigma","N_0")),
    rw.sd=rw.sd(r=0.02,K=0.02,sigma=0.02,N_0=ivp(0.02)),
    cooling.fraction.50=0.5,
    paramnames=c("r","K","sigma","N_0","b"),
    statenames=c("N")
  ) -> mf1

plot(mf1)
