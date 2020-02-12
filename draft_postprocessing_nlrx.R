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
library(patchwork)
##test2=unnest_simoutput(nl_object) ##some conflicts when extricting patch and turtle variables?



output_patches<-exp1_2020 %>% 
  dplyr::select(ticks,start_allee_thres,slope_disp_mean,disp0_mean,fecundity,K,`random-seed`,metrics.patches) %>% 
  unnest(cols=c(metrics.patches)) %>% 
  left_join(trt_grid)%>% 
  mutate(replicateID=paste(K,treatment,`random-seed`))%>% 
  mutate(Location_full=paste(replicateID,pxcor,ticks))


tab=output_patches %>% filter(N_postdispersal>0) %>% 
  group_by(treatment,replicateID,ticks,K) %>% 
  filter(pxcor>(max(pxcor)-1)|pxcor==0) %>%
  mutate(is.front=pxcor!=0,seedID=`random-seed`) %>% 
  group_by(is.front,seedID,add=TRUE) %>% 
  summarise(N_allele0=sum(N_allele0),N_allele1=sum(N_allele1),N_postdispersal=sum(N_postdispersal)) %>%
  mutate(Hexp=2*(N_allele0/N_postdispersal)*(N_allele1/N_postdispersal),
         P0=(N_allele0/N_postdispersal),P1=(N_allele1/N_postdispersal)) %>% 
  mutate(denom_Hexp=(N_postdispersal^2),nume_Hexp=2*N_allele0*N_allele1) %>% 
#group_by(treatment,replicateID,ticks,is.front) %>%
#summarise(Hexp=mean(Hexp)) %>% 
ungroup()%>% 
  filter(treatment != "strong Allee effect (a = 20)") 

ggplot(tab) +geom_line(aes(x=ticks,y=Hexp,group=replicateID,col=treatment))+facet_wrap(~treatment+is.front+K)
tab %>% group_by(ticks,treatment,is.front,K) %>% summarise(Hexp=mean(Hexp)) %>% 
  ggplot() +geom_line(aes(x=ticks,y=Hexp,col=treatment))+facet_wrap(~is.front+K)

####approche variance in fraction (gandhi 2019 pnas). works by "averaging across" replicates, no zeros in dataset
#### simplifies problems caused in Hexp when there are zeros (beta, loglinear can't do anything)

tab_genetics=tab %>% group_by(ticks,treatment,is.front,K) %>% summarise(medHexp=median(Hexp),Hexp=mean(Hexp),varP1=var(P1),varP0=var(P0)) %>% 
  mutate(K=factor(K))
ggplot(tab_genetics)+geom_line(aes(x=ticks,y=varP1,col=treatment))+facet_wrap(~is.front+K)


prior_genetics <- c(
  set_prior("normal(0,1)",nlpar=c("logdecay","logitHasym"),class="b"),
  set_prior("normal(0,1)",class="sigma")
)

bf_genetics <- bf(varP1~log(Hasym*(1-exp(-(ticks)*decay))),
                  nlf(Hasym~inv_logit(logitHasym)),
                  nlf(decay~10^(logdecay)),
                  logitHasym~1,
                  #logitHasym~0+treatment,
                  logdecay~0+treatment:is.front:K,
                  nl=TRUE,family=student(log))

mod_genetics=brm(bf_genetics,
        data=tab_genetics,
        iter=2000,chains=4,
        prior=c(prior_genetics,set_prior("gamma(2,0.1)",class="nu"))
        )


tab_genetics %>% ungroup() %>% 
  select(treatment,is.front,K) %>% 
  distinct() %>% 
  expand_grid(ticks=(1:100)) %>% 
  add_fitted_draws(mod_genetics,re_formula = NA) %>% 
  group_by(ticks, treatment) %>% 
  ggplot()+
  geom_line(data=tab_genetics,aes(x=ticks, y=varP1,col=treatment))+
  stat_lineribbon(aes(x=ticks,y=.value),.width=0.95,alpha=0.25)+
  scale_y_continuous("Variance in fraction")+
  facet_grid(cols=vars(treatment),rows=vars(is.front,K))+
  theme_bw()



tab_genetics %>% ungroup() %>% 
  select(treatment,is.front,K) %>% 
  distinct() %>% 
  mutate(ticks=1) %>% 
  add_fitted_draws(mod_genetics,re_formula = NA,nlpar="logdecay") %>% 
  mutate(.value=exp(.value)) %>% 
  ggplot()+
  geom_violin(aes(x=is.front,fill=K,col=K,y=.value))+
  facet_grid(cols=vars(treatment))+
  scale_y_continuous(name=expression(paste("genetic diversity decay rate  ",lambda)))+
  scale_x_discrete("Location",labels=c("core","front"))+
  theme_bw()


tab_genetics %>% ungroup() %>% 
  select(treatment,is.front,K) %>% 
  distinct() %>% 
  mutate(ticks=1) %>% 
  add_fitted_draws(mod_genetics,re_formula = NA,nlpar="logdecay") %>% 
  mutate(.value=exp(.value)) %>% group_by(treatment,K,is.front) %>% compare_levels(variable=.value,by=K) %>% median_hdi()

tab_speed=output_patches %>% filter(N_postdispersal>0) %>%
  group_by(treatment,replicateID,ticks) %>% 
  select(replicateID,pxcor,ticks,treatment,K,velocity_fisher, seedID = `random-seed`) %>% 
  filter(pxcor==max(pxcor)) %>% 
  group_by(replicateID) %>% 
  #arrange(ticks) %>% 
  mutate(front=pxcor) %>% 
  #mutate(front_prev=lag(front)) %>% 
  #mutate(front_prev=replace_na(front_prev,0)) %>% 
  mutate(speed=abs(front)/ticks) %>% 
  #mutate(stalled=front_prev>=front) %>% 
  ungroup()%>% 
  filter(treatment != "strong Allee effect (a = 20)") 

ggplot(tab_speed)+
  geom_line(aes(x=ticks,y=front,group=replicateID))+
  facet_wrap(~treatment+K)+
  geom_abline(aes(intercept=0,slope=(3/(2*sqrt(2)))*velocity_fisher),col="red")+
  geom_abline(aes(intercept=0,slope=velocity_fisher),col="blue")



prior_front <- c(
  #set_prior("normal(0,1)",nlpar=c("logspeedstart"),class="b"),
  set_prior("normal(0,1)",nlpar=c("logspeedasym"),class="b"),
  set_prior("normal(0,1)",nlpar=c("lograte"),class="b"),
  set_prior("normal(0,1)",nlpar=c("logspeedasym","lograte"),class="sd"),
  set_prior("normal(0,1)",class="sigma"),
  set_prior("lkj(2)",class="cor")
)



bf_front <- bf(front~log(ticks*(speedasym+(1-speedasym)*exp(-(ticks-1)/rate))),
   nlf(rate~10*exp(lograte)),
   nlf(speedasym~inv_logit(logspeedasym)),
   #nlf(speedstart~exp(logspeedstart)),
   logspeedasym~0+treatment:factor(K)+(1|1|replicateID),
   #logspeedstart~0+treatment+(1|1|replicateID),
   lograte~0+treatment:factor(K)+(1|1|replicateID),
   family=lognormal,nl=TRUE)

mod_front=brm(bf_front,
        data=subset(tab_speed,ticks>0),chains=4,iter=2000,
        prior=prior_front,control=list(adapt_delta=0.8,max_treedepth=10))


p1 <- tab_speed %>% 
  ungroup() %>% 
  mutate(replicateID= unique(replicateID)[1]) %>% 
  select(treatment,replicateID,K) %>% 
  distinct() %>% 
  expand_grid(ticks=(1:100)) %>% 
  add_fitted_draws(mod_front,re_formula=NA) %>% 
  mutate(.value_front=.value) %>% 
  group_by(ticks, treatment) %>% 
  ggplot()+
  geom_line(data=tab_speed,aes(x=ticks, y=front,group=replicateID),col="grey",alpha=0.2)+
  stat_lineribbon(aes(x=ticks,y=.value_front),.width=0.95,alpha=0.5)+
  scale_y_continuous("Front location (patches from release)")+
  scale_x_continuous("time (generations)")+
  facet_grid(cols=vars(treatment),rows=vars(paste("K = ",K,sep="")),space= "free_x",scale="free_x")+
  geom_abline(data=trt_grid,aes(intercept=0,slope=(3/(2*sqrt(2)))*velocity_fisher),col="red",lty=2)+
  geom_abline(data=trt_grid,aes(intercept=0,slope=velocity_fisher),col="blue",lty=2)

p2<- tab_speed %>% 
  ungroup() %>% 
  mutate(replicateID= replicateID[1],ticks=1) %>% 
  select(treatment,ticks,replicateID,K) %>% 
  distinct() %>%
  add_fitted_draws(mod_front,nlpar="logspeedasym",re_formula=NA) %>% 
  mutate(.value=invlogit(.value)) %>% 
  ggplot()+
  geom_ribbon(data=trt_grid,aes(x=3*(-1+K/225),ymin = velocity_fisher,ymax=(3/(sqrt(2)*2))*velocity_fisher),col="grey80")+
  geom_jitter(data=tab_speed %>% filter(ticks==100),aes(x=factor(K),y=front/ticks),col="grey",alpha=0.3)+
  #geom_hline(data=trt_grid,aes(yintercept = (3/(sqrt(2)*2))*velocity_fisher),lty=2)+
  geom_violin(aes(x=factor(K),group=K,y=.value),alpha=0.25)+
  scale_y_continuous("Asymptotic velocity (posterior)",limits = c(0,1)) +
  scale_x_discrete("equilibrium population size K")+
  facet_grid(cols=vars(treatment),space= "free_x",scale="free_x")

tab_speed %>% 
  ungroup() %>% 
  mutate(replicateID= replicateID[1],ticks=1) %>% 
  select(treatment,ticks,replicateID,K) %>% 
  distinct() %>%
  add_fitted_draws(mod_front,nlpar="logspeedasym",re_formula=NA) %>%
  mutate(.value=invlogit(.value))%>% 
  left_join(trt_grid) %>%
  mutate(.value=.value/velocity_fisher) %>% 
  mean_hdi()
 

p1/p2 & theme_bw() & theme(legend.position="none")

tab_speed %>% 
  ungroup() %>% 
  mutate(replicateID= replicateID[1],ticks=1) %>% 
  select(treatment,ticks,replicateID,K) %>% 
  distinct() %>%
  add_fitted_draws(mod_front,nlpar="logspeedasym",re_formula=NA) %>% 
  mutate(.value=invlogit(.value)) %>% group_by(treatment,K) %>% compare_levels(variable=.value,by=K) %>% median_hdi()








#### for analysing size of new propagules:
#### works only if all ticks are recorded, not one out of two

tab_propa <- output_patches %>% filter(N_postdispersal>0) %>% 
  group_by(treatment,replicateID,ticks,K) %>% 
  mutate(front=max(pxcor),seedID=`random-seed`) %>%
  filter(pxcor >= (front-1) ) %>% 
  mutate(is.front=pxcor==front) %>% 
  group_by(treatment,replicateID,is.front,K) %>% 
  arrange(ticks) %>% 
  mutate(front_prev=lag(front))%>%
  mutate(front_prev=replace_na(front_prev,0)) %>% 
  filter(front>front_prev & ticks >1) %>%  ##we only keep case where the front advanced (new pops) 
  ###remove G1 because dispersal in only one direction so messier for comparison to deterministic
  select(replicateID,front,is.front,front_prev,ticks,N_predispersal,N_postdispersal) %>% 
  pivot_wider(names_from=is.front,values_from=c(N_predispersal,N_postdispersal)) %>% 
  left_join(trt_grid) %>% 
  filter(treatment == "control" | treatment == "reduced connectedness") %>% 
  mutate(N_founding_actual=N_postdispersal_TRUE,N_founding_deterministic=N_predispersal_FALSE*disp0_mean/2) %>% 
  ungroup() %>% 
  select(N_founding_actual,N_founding_deterministic,treatment,replicateID,ticks,K) %>% 
  mutate(ratio=N_founding_actual/N_founding_deterministic) %>% 
  mutate(K=factor(K))

mod_propa1=brm(N_founding_actual~0+treatment:K+(1|replicateID),
              family=negbinomial,data=tab_propa,
              prior=c(set_prior("normal(0,1)",class="b"),
                      set_prior("normal(0,1)",class="sd")),chains=4,iter=100)

 
mod_propa2=brm(ratio~0+treatment:K+(1|replicateID),
               family=lognormal,data=tab_propa,
               prior=c(set_prior("normal(0,1)",class="b"),
                       set_prior("normal(0,1)",class="sd"),
                       set_prior("normal(0,1)",class="sigma")),chains=2,iter=100)

 
  
tab_propa %>% 
    mutate(replicateID=replicateID[1]) %>% 
    select(treatment,K,replicateID) %>% distinct() %>% add_fitted_draws(mod_propa2,re_form=NA) %>% 
    ggplot()+
  geom_boxplot(data=tab_propa,aes(x=treatment,y=ratio),alpha=0.1)+
  geom_violin(aes(x=treatment,y=.value))+
    facet_wrap(~K)
###double check weird behavior of add_fitted_draws here (AND ELSEWHERE)

tab_propa %>% 
  mutate(replicateID=replicateID[1]) %>% 
  select(treatment,K,replicateID) %>% distinct() %>% add_fitted_draws(mod_propa2,re_formula=NA,scale="linear") %>% 
  ggplot()+
  #geom_boxplot(data=tab_propa,aes(x=treatment,y=ratio),alpha=0.1)+
  geom_eye(aes(x=treatment,y=exp(.value)))+ 
  scale_y_continuous("new population size[observed] / new population size [deterministic] (posterior mean)")+
  scale_x_discrete("landscape type")+
  facet_wrap(~paste ("K = ",K))+theme_bw()



tab_propa %>% 
  mutate(replicateID=replicateID[1]) %>% 
  select(treatment,K,replicateID) %>% distinct() %>% add_fitted_draws(mod_propa1,re_formula=NA,scale="linear") %>% 
  ggplot()+
  #geom_boxplot(data=tab_propa,aes(x=treatment,y=ratio),alpha=0.1)+
  geom_eye(aes(x=treatment,y=exp(.value)))+ 
  scale_y_continuous("new population size[observed] (posterior)")+
  scale_x_discrete("landscape type")+
  facet_wrap(~paste ("K = ",K))+theme_bw()
  