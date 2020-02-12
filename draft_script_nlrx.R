library(nlrx)
library(arm)
library(tidyverse)

###add criterion to decide whether landscape is 2sided or 1 sided and
### set up to one sided for this paper
### indeed if not core patch may receive immigrants from 2 patches rather than 1
### may contribute to lower decay of core genetic diversity in simu compared to expe

####

nl_object <- nl(nlversion = "6.1.0", nlpath = "C:/Program Files/NetLogo 6.1.0",
               modelpath = "D:/Maxime/Documents/POSTDOC INRA SOPHIA/IBM/range-expansion-ibm-pp-2019/range_expansion.nlogo", 
               jvmmem = 1024)

###designing the experiment
trt_grid=tibble(treatment=c("control","weak Allee effect (a = 0.95)",
                            "density-dependent dispersal",
                            "reduced connectedness"),
                start_allee_thres=c(0,0.95,0,0),
                slope_disp_mean=c(0,0,1,0), ##leads to roughly doubling of disp between 0 and K, realistic
                disp0_mean=c(0.2,0.2,0.2,0.1)
) %>% 
  mutate(fecundity = 5) %>% 
  expand_grid(K=c(225,450)) %>% 
  ###fecundity = 5 ; 
  ## close to Trichogramma conditions given their sex ratio of about 50% and assuming 24 or 48 h of egg laying, as in experiments see papers
  mutate(velocity_fisher=2*sqrt(
    log(fecundity)*(1 - 1 / K) * (1 - start_allee_thres / 1)*  ##fecundity term at N = 1 individual
                                  0.5*invlogit(logit(disp0_mean)+(1/K)*slope_disp_mean)) ##dispersal term
    )

#### dispersal formula is birzu prob because not bounded
#### our own is bounded and close to what empiricists would do in DDD situations

###when 0<a<1, weak allee effect (cf moreljournel paper). Do one treatment with a =0.5
### stop studying multiple K. Do only K =450 as in experiment

##strong allee effect always pushed even if velocity fisher not defined see birzu 2018 PNAS
nreplicates=100

duration=100
###############################

nl_object@experiment <- experiment(
  expname="experiment-art1-2020",
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
##i'm using the growth function of morel journal oikos, but not teh density function because unsuitable
## it is the altwegg one, which posits that the DI dispersal rate is the rate at K, rather than at low densities
## (NB leads to >1 dispersal rate when negative DDD if unchecked)

set.seed(42) ##we set seed here to guarantee the seeds selected below are the same everytime
nl_object@simdesign<-simdesign_distinct(nl=nl_object,nseeds=nreplicates)

exp1_2020=run_nl_all(nl_object)

setsim(nl_object, "simoutput")<-exp1_2020

