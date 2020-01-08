library(nlrx)

library(arm)
library(tidyverse)


library(modelr)
library(RColorBrewer)
library(cowplot)

###inference
library(coda)
library(rstan)
library(bayesplot)
library(brms)
rstan_options(auto_write = TRUE)
options(mc.cores = 2)
library(tidybayes)

###add criterion to decide whether landscape is 2sided or 1 sided and
### set up to one sided for this paper
### indeed if not core patch may receive immigrants from 2 patches rather than 1
### may contribute to lower decay of core genetic diversity in simu compared to expe

####
#proof of principle that reduced connectivity restricts dispersal at low N, due to binomial properties
expand.grid(dispersal.rate=c(0.05,0.1,0.2,0.4,0.8),N=1:450) %>% 
ggplot()+
  geom_line(aes(x=N,y=1-pbinom(0,N,dispersal.rate),col=factor(dispersal.rate)))+
  scale_y_continuous(name="Probability of at least 1 individual dispersing")+
  scale_x_log10(name="Population size before dispersal phase")+
  scale_color_discrete(name="Nominal emigration probability") +
  theme_classic()

expand.grid(dispersal.rate=c(0.05,0.1,0.2,0.4,0.8),N=1:450) %>% 
  ggplot()+
  geom_line(aes(x=N,y=pbinom(2*N*dispersal.rate,N,dispersal.rate,lower.tail = FALSE),col=factor(dispersal.rate)))+
  scale_y_continuous(name="Probability of dispersing > twice the nominal rate")+
  scale_x_log10(name="Population size before dispersal phase")+
  scale_color_discrete(name="Emigration probability") +
  theme_classic()

expand.grid(dispersal.rate=c(0.05,0.1,0.2,0.4,0.8),N=1:450) %>% 
  ggplot()+
  geom_line(aes(x=N,y=qbinom(0.5,N,dispersal.rate)/N,col=factor(dispersal.rate)))+
  scale_y_continuous(name="Median dispersal rate")+
  scale_x_log10(name="Population size before dispersal phase")+
  scale_color_discrete(name="Nominal emigration probability") +
  theme_classic()
### the effect vanishes around N=200 and become strong at N<100 so probably too rare to be detected in previous studies with K in tens of thousands
### so even with low dispersal founding units rarely below the 100s

###SIDEBAR EMPIRICAL DATA
###use model predicted rate of diversity loss and N at front
### to estimate diversity exponent birzu 2018 PNAS fig 5
### should be lower the more pushsed, and = -1 for truly pushed

### TO DO: really think about the reduced connec mechanism (addd DDD on eggs?)
### do the empirical stuff
###count more than 1 patch for the diversity?

nl_object <- nl(nlversion = "6.1.0", nlpath = "C:/Program Files/NetLogo 6.1.0",
               modelpath = "D:/Maxime/Documents/POSTDOC INRA SOPHIA/IBM/range-expansion-ibm-pp-2019/range_expansion.nlogo", 
               jvmmem = 1024)

###designing the experiment
trt_grid=tibble(treatment=c("control","weak Allee effect (a = 0.95)","strong Allee effect (a = 20)",
                            "density-dependent dispersal",
                            "reduced connectedness"),
                start_allee_thres=c(0,0.95,20,0,0),
                slope_disp_mean=c(0,0,0,1,0), ##leads to roughly doubling of disp between 0 and K, realistic
                disp0_mean=c(0.2,0.2,0.2,0.2,0.1)
) %>% 
  mutate(K=450,fecundity=5) %>% 
  mutate(velocity_fisher=2*sqrt(
    log(fecundity)*(1 - 1 / K) * (1 - start_allee_thres / 1)*  ##fecundity term at N = 1 individual
                                  0.5*invlogit(logit(disp0_mean)+(1/K)*slope_disp_mean)) ##dispersal term
    )

#### dispersal formula is birzu prob because not bounded
#### our own is bounded and close to what empiricists would do in DDD situations

###when 0<a<1, weak allee effect (cf moreljournel paper). Do one treatment with a =0.5
### stop studying multiple K. Do only K =450 as in experiment

##strong allee effect always pushed even if velocity fisher not defined see birzu 2018 PNAS
nreplicates=50

duration=100
###############################

nl_object@experiment <- experiment(
  expname="test",
  outpath="D:/Maxime/Documents/POSTDOC INRA SOPHIA/IBM/range-expansion-ibm-pp-2019/output",
  repetition = 1,
  tickmetrics = "true",
  idsetup = "setup",
  idgo="go",
  stopcond="not any? turtles",
  runtime=duration,
  evalticks= 1:duration,
  metrics=c("ticks"),
  metrics.patches = c("pxcor","N_predispersal","N_postdispersal","N_allele0","N_allele1"),
  constants=list(
    "trait_variation"= "\"reshuffled\"",   ###careful with the way string variable/csts must be entered
    #"K"=450,
    #"start_allee_thres"=5,
    "duration"=duration,
    "logit_disp0_sd"=0,
    #"disp0_mean"=0.2, ###disp rate of 0.2 at low densities, matches empirical in general and tricho
    "slope_disp_sd"=0
    #"slope_disp_mean"=0,
    #"fecundity"=1.2,
    ),
  
  variables=list( ###needs to be a nested list
    "start_allee_thres"=list(values=trt_grid$start_allee_thres), ##20 matches Marjorie's PCI paper
    "slope_disp_mean"=list(values=trt_grid$slope_disp_mean),
    "disp0_mean"=list(values=trt_grid$disp0_mean),
    "fecundity"=list(values=trt_grid$fecundity),
    "K"=list(values=trt_grid$K)
  )
  
)


###NB in morel journel model is default dspersal rate is fixed, not even binomially distributed, fixed
### can be coded by ranking indiv in a patch by their ID and saying first N disperse
### but weird
### better imho (and more novel) to stick with binomial and more stochastic than binomial
### in addition allow to discuss discrepancy

##i'm using the growth function of morel journal oikos, but not teh density function because unsuitable
## it is the altwegg one, which posits that the DI dispersal rate is the rate at K, rather than at low densities
## (NB leads to >1 dispersal rate when negative DDD if unchecked)

set.seed(42) ##we set seed here to guarantee the seeds selected below are the same everytime
nl_object@simdesign<-simdesign_distinct(nl=nl_object,nseeds=nreplicates)

test=run_nl_all(nl_object)



setsim(nl_object, "simoutput")<-test

##test2=unnest_simoutput(nl_object) ##some conflicts when extricting patch and turtle variables?

output_patches=unnest(test,cols=c(metrics.patches)) %>% 
  left_join(trt_grid)%>% 
  mutate(replicateID=paste(K,treatment,`random-seed`))%>% 
  mutate(Location_full=paste(replicateID,pxcor,ticks))

tab=output_patches %>% filter(N_postdispersal>0) %>% 
  group_by(treatment,replicateID,ticks) %>% 
  filter(pxcor>(max(pxcor)-1)|pxcor==0) %>%
  mutate(is.front=pxcor!=0) %>% 
  group_by(is.front,add=TRUE) %>% 
  summarise(N_allele0=sum(N_allele0),N_allele1=sum(N_allele1),N_postdispersal=sum(N_postdispersal)) %>%
  mutate(Hexp=2*(N_allele0/N_postdispersal)*(N_allele1/N_postdispersal),
         P0=(N_allele0/N_postdispersal),P1=(N_allele1/N_postdispersal)) %>% 
  #group_by(treatment,replicateID,ticks,is.front) %>%
  #summarise(Hexp=mean(Hexp)) %>% 
ungroup()

ggplot(tab) +geom_line(aes(x=ticks,y=variance,group=replicateID,col=treatment))+facet_wrap(~treatment+is.front)
tab %>% group_by(ticks,treatment,is.front) %>% summarise(Hexp=mean(Hexp)) %>% ggplot() +geom_line(aes(x=ticks,y=Hexp,col=treatment))+facet_wrap(~is.front)


mod=brm(bf(Hexp|trunc(lb=0,ub=0.5)~(0.5*exp(-(ticks)*decay)),
           #nlf(H0~exp(logH0)),
           nlf(decay~exp(logdecay)),
           logdecay~0+treatment+treatment:is.front+(is.front|1|replicateID),
           #logH0~1+(1|1|replicateID),
           nl=TRUE,family=student),
        data=tab,
        iter=200,chains=2,
        prior=c(set_prior("normal(0,5)",nlpar=c("logdecay"),class="b"),
                set_prior("normal(0,5)",nlpar=c("logdecay"),class="sd"),
                set_prior("normal(0,1)",class="sigma"),
                set_prior("gamma(2,0.1)",class="nu"),
                set_prior("lkj(2)",class="cor")
        ))




####approche variance in fraction (gandhi 2019 pnas). works by "averaging across" replicates, no zeros in dataset
#### simplifies problems caused in Hexp when there are zeros (beta, loglinear can't do anything)

tab2=tab %>% group_by(ticks,treatment,is.front) %>% summarise(Hexp=mean(Hexp),varP1=var(P1),varP2=var(P2)) 
mod=brm(bf(varP1~log(0.25*(1-exp(-(ticks)*decay))),
           #nlf(H0~exp(logH0)),
           nlf(decay~exp(logdecay)),
           logdecay~0+treatment:is.front,
           #logH0~1+(1|1|replicateID),
           nl=TRUE,family=student(log)),
        data=tab2,
        iter=2000,chains=4,
        prior=c(set_prior("normal(0,5)",nlpar=c("logdecay"),class="b"),
                set_prior("gamma(2,0.1)",class="nu")
        ))


fit=posterior_samples(mod) %>% select(contains("logdecay")) %>%
  pivot_longer(everything(),values_to="value",names_to="name") %>%
  mutate(is.front=str_detect(name,"is.frontTRUE"),treatment=str_remove(name,pattern="b_logdecay_treatment")) %>% 
  mutate(treatment=str_remove(treatment,":is.frontFALSE")) %>%
  mutate(treatment=str_remove(treatment,":is.frontTRUE")) %>% select(-name)
ggplot(fit)+geom_violin(aes(x=factor(is.front),y=exp(value)))+
  facet_grid(cols=vars(treatment))+
  scale_y_log10(name=expression(paste("genetic diversity decay rate  ",lambda)))+
  scale_x_discrete("Location",labels=c("core","front"))


####extract front position info from both front
output_patches %>% filter(N_postdispersal>0) %>% 
  group_by(replicateID,ticks) %>% 
  mutate(maxfront=max(pxcor)) %>% 
  ungroup() %>% 
  select(maxfront,ticks,replicateID,treatment,K,velocity_fisher) %>% 
  distinct() %>% 
  ggplot()+
  geom_line(aes(x=ticks,y=maxfront/ticks,group=replicateID))+
  facet_wrap(~treatment+K)+
  geom_hline(aes(yintercept=(3/(2*sqrt(2)))*velocity_fisher),col="red")+
  geom_hline(aes(yintercept=velocity_fisher),col="blue")



tabspeed=output_patches %>% filter(N_postdispersal>0) %>% 
  group_by(treatment,replicateID,ticks) %>% 
  select(replicateID,pxcor,ticks,treatment,K,velocity_fisher) %>% 
  filter(pxcor==max(pxcor)) %>% 
  group_by(replicateID) %>% 
  arrange(ticks) %>% 
  mutate(front=pxcor) %>% 
  mutate(front_prev=lag(front)) %>% 
  mutate(front_prev=replace_na(front_prev,0)) %>% 
  mutate(speed=abs(front)/ticks) %>% 
  mutate(stalled=front_prev>=front) %>% 
  ungroup()

ggplot(tabspeed)+
  geom_line(aes(x=ticks,y=speed,group=replicateID))+
  facet_wrap(~treatment+K)+
  geom_hline(aes(yintercept=(3/(2*sqrt(2)))*velocity_fisher),col="red")+
  geom_hline(aes(yintercept=velocity_fisher),col="blue")


mod=brm(bf(speed~log(speedasym+(speedstart-speedasym)*exp(-(ticks)*rate)),
           nlf(rate~exp(lograte)),nlf(speedasym~exp(logspeedasym)),nlf(speedstart~exp(logspeedstart)),
           logspeedasym~0+treatment+(1|1|replicateID),#+(1|2|replicateID:is.leftfront),
           logspeedstart~0+treatment+(1|1|replicateID),#+(1|2|replicateID:is.leftfront),
           lograte~0+treatment+(1|1|replicateID),#+(1|2|replicateID:is.leftfront),
           #sigma~0+treatment,
           family=lognormal,nl=TRUE),
        data=subset(tabspeed,ticks<=100),chains=4,iter=2000,
        prior=c(
          set_prior("normal(0,1)",nlpar=c("logspeedstart"),class="b"),
          set_prior("normal(0,1)",nlpar=c("logspeedasym"),class="b"),
          set_prior("normal(0,2)",nlpar=c("lograte"),class="b"),
          set_prior("normal(0,1)",nlpar=c("logspeedasym","lograte","logspeedstart"),class="sd"),
         set_prior("normal(0,1)",class="sigma"),
         #set_prior("normal(0,5)",dpar="sigma",class="b"),
          set_prior("lkj(2)",class="cor")
        ),control=list(adapt_delta=0.8,max_treedepth=10))

fronts$speed=fronts$front/fronts$Generation
###replace the empiric front by this one:
###model 1/rate instead (farther from 0, easier to converge)
### 3/rate: 0.05 of delta remainin, 6/rate, 0.01
mod=brm(bf(speed~log(speedasym+(speedstart-speedasym)*exp(-Generation*rate)),nlf(rate~exp(lograte)),nlf(speedasym~exp(logspeedasym)),nlf(speedstart~exp(logspeedstart)),
           logspeedasym~(Mix2+Mix3)*is.pushed+(1|1|landscapeID),
           logspeedstart~(Mix2+Mix3)*is.pushed+(1|1|landscapeID),lograte~(Mix2+Mix3)*is.pushed+(1|1|landscapeID),
           family=lognormal,nl=TRUE),
        data=subset(fronts,Generation>0),chains=4,iter=1000,
        prior=c(
          set_prior("normal(0,1)",nlpar=c("logspeedasym","lograte","logspeedstart"),class="b"),
          set_prior("exponential(1)",nlpar=c("logspeedasym","lograte","logspeedstart"),class="sd"),
          set_prior("lkj(2)",class="cor"),
          set_prior("exponential(1)",class="sigma")
        ),control=list(adapt_delta=0.95,max_treedepth=15))


newdata=tabspeed %>% ungroup() %>% select(treatment) %>% distinct() %>% 
  mutate(replicateID=tabspeed$replicateID[1],ticks=1) %>% 
  add_fitted_draws(mod,re_formula=NA,nlpar="logspeedasym") %>% 
  mutate(.value=exp(.value)) %>% 
  median_hdi(.width=0.95)


ggplot(tabspeed)+
  geom_line(aes(x=ticks,y=speed,group=replicateID),col="grey")+
  facet_grid(cols=vars(treatment))+
  geom_hline(aes(yintercept=(3/(2*sqrt(2)))*velocity_fisher),col="red")+
  geom_hline(aes(yintercept=velocity_fisher),col="blue")+
  geom_pointinterval(data=newdata,aes(x=100,y=.value))


prior_genetics=c(
  set_prior("normal(0,1)",nlpar=c("logdecay","logH0"),class="b"),
  set_prior("exponential(1)",nlpar=c("logdecay","logH0"),class="sd"),
  set_prior("lkj(2)",class="cor"),
  set_prior("exponential(1)",class="sigma")
)
bf_genetics=bf((Hexp)~log(exp(logH0)*exp(-exp(logdecay)*Generation)),
               logH0~(Mix2+Mix3)+(1|1|landscapeID),
               logdecay~is.front*is.pushed*(Mix2+Mix3)+(is.front|1|landscapeID),
               family=student("log"),nl=TRUE)

mod_genetics=brm(bf_genetics,
                 data=data_genetics,
                 chains=4,iter=2000,warmup=1000,prior=c(prior_genetics,set_prior("gamma(2,0.1)",class="nu")),
                 seed=42,control=list(adapt_delta=0.8,max_treedepth=10)
)

                     




##beta version
prior_genetics=c(
  set_prior("normal(0,1)",nlpar=c("logdecay","logH0"),class="b"),
  set_prior("exponential(1)",nlpar=c("logdecay","logH0"),class="sd"),
  set_prior("lkj(2)",class="cor")
)

bf_genetics=bf((Hexp)~(inv_logit(logH0)*exp(-exp(logdecay)*Generation)),
               logH0~(Mix2+Mix3)+(1|1|landscapeID),
               logdecay~is.front*is.pushed*(Mix2+Mix3)+(is.front|1|landscapeID),phi~1,
               family=Beta(identity),nl=TRUE)

mod_gen=brm(bf_genetics,
                 data=data_genetics,
                 chains=4,iter=2000,warmup=1000,prior=c(prior_genetics),
                 seed=42,control=list(adapt_delta=0.8,max_treedepth=10)
)


compare_ic(loo(mod_gen),loo(mod_genetics))