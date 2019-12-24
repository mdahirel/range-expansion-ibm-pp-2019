

library(nlrx)

library(arm)
library(tidyverse)

nl_object <- nl(nlversion = "6.1.0", nlpath = "C:/Program Files/NetLogo 6.1.0",
               modelpath = "D:/Maxime/Documents/POSTDOC INRA SOPHIA/IBM/range-expansion-ibm-pp-2019/range_expansion.nlogo", 
               jvmmem = 1024)

nl_object@experiment <- experiment(
  expname="test",
  outpath="D:/Maxime/Documents/POSTDOC INRA SOPHIA/IBM/range-expansion-ibm-pp-2019/output",
  repetition = 1,
  tickmetrics = "true",
  idsetup = "setup",
  idgo="go",
  stopcond="not any? turtles",
  runtime=300,
  evalticks= c(1:10,seq(25,300,25)),
  metrics=c("ticks","count turtles"),
  metrics.turtles = list("turtles"=c("pxcor","pxcor_birth","neutral_allele","logit_disp0","disp_slope","adult")),
  #variables=list("trait_variation"= c("\"evolutionary\"", "\"reshuffled\"" )),
  constants=list(
    "trait_variation"= "\"reshuffled\"",   ###carfeul with the way string variable/csts must be entered
    "competition_type"="\"beverton_holt_like\"",
    "allee_effects_repro"="\"no\"",
    "density_dependent_dispersal"= "false",
    "K"=100,
    "allele_number"=20,
    "start_allee_thres"=5,
    "logit_disp0_sd"=1,
    "logit_disp0_mean"=-1,
    "fecundity"=2,
    "spatial_burnin" = -1   ###change it to a dirft-only/dispersal allowed switch ###seed first patches, not only the first one
  )  ##so that would mean allowed to move in original range
  
)

nl_object@simdesign<-simdesign_simple(nl=nl_object,nseeds=5)

test=run_nl_all(nl_object)



setsim(nl_object, "simoutput")<-test
#test2=unnest_simoutput(nl_object) ##buggy cause the function has not integrated yet updates to main unnest

test2=unnest(test,cols=c(metrics.turtles)) 


test2 %>% 
  group_by(`random-seed`,ticks,pxcor,adult) %>% 
  summarise(Nturtles = length(pxcor),
         meandisp=invlogit(mean(`logit_disp0`)),
         Nallele=length(unique(neutral_allele))) %>% 
  ungroup() %>% 
  filter(adult==1) %>% 
  ggplot()+geom_line(aes(x=pxcor_birth,y=Nturtles,group=`random-seed`))+facet_wrap(~ticks)

                     