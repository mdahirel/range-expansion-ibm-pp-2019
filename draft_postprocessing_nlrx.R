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
##test2=unnest_simoutput(nl_object) ##some conflicts when extricting patch and turtle variables?

output_patches=unnest(test,cols=c(metrics.patches)) %>% 
  left_join(trt_grid)%>% 
  mutate(replicateID=paste(K,treatment,`random-seed`))%>% 
  mutate(Location_full=paste(replicateID,pxcor,ticks))

tab=output_patches %>% filter(N_postdispersal>0) %>% 
  group_by(treatment,replicateID,ticks,K) %>% 
  filter(pxcor>(max(pxcor)-1)|pxcor==0) %>%
  mutate(is.front=pxcor!=0) %>% 
  group_by(is.front,add=TRUE) %>% 
  summarise(N_allele0=sum(N_allele0),N_allele1=sum(N_allele1),N_postdispersal=sum(N_postdispersal)) %>%
  mutate(Hexp=2*(N_allele0/N_postdispersal)*(N_allele1/N_postdispersal),
         P0=(N_allele0/N_postdispersal),P1=(N_allele1/N_postdispersal)) %>% 
  mutate(denom_Hexp=(N_postdispersal^2),nume_Hexp=2*N_allele0*N_allele1)
#group_by(treatment,replicateID,ticks,is.front) %>%
#summarise(Hexp=mean(Hexp)) %>% 
ungroup()

ggplot(tab) +geom_line(aes(x=ticks,y=Hexp,group=replicateID,col=treatment))+facet_wrap(~treatment+is.front)
tab %>% group_by(ticks,treatment,is.front,K) %>% summarise(Hexp=mean(Hexp)) %>% ggplot() +geom_line(aes(x=ticks,y=Hexp,col=treatment))+facet_wrap(~is.front+K)


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

tab2=tab %>% group_by(ticks,treatment,is.front,K) %>% summarise(Hexp=mean(Hexp),varP1=var(P1),varP0=var(P0)) 
ggplot(tab2)+geom_line(aes(x=ticks,y=varP1,col=treatment))+facet_wrap(~is.front+K)


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
  group_by(replicateID,ticks,treatment,K,velocity_fisher) %>% 
  summarize(maxfront=max(pxcor)) %>% 
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

