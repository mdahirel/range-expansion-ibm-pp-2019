

library(nlrx)

library(arm)
library(tidyverse)


library(modelr)
library(RColorBrewer)
library(cowplot)

library(vegan)
library(RVAideMemoire)
library(matrixStats)
library(readr)
library(adegenet)
library(PopGenReport)
library(pegas)
library(poppr)
library(ggthemes)

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
  runtime=200,
  evalticks= c(1:10,seq(25,200,25)),
  metrics=c("ticks"),
  metrics.patches = c("pxcor","N_predispersal","N_postdispersal"),
  metrics.turtles = list("turtles"=c("neutral_allele","neutral_allele_shuffled","parentID","birth_patch","xcor")),
  constants=list(
    "trait_variation"= "\"reshuffled\"",   ###careful with the way string variable/csts must be entered
    #"K"=500,
    "allele_number"=20,
    #"start_allee_thres"=5,
    "duration"=200,
    "logit_disp0_sd"=0,
    "logit_disp0_mean"=0,
    "slope_disp_sd"=0,
    "slope_disp_mean"=0,
    "fecundity"=2
  ),
  
  variables=list( ###needs to be a nested list
    "K"=list(values=c(50,100,250,500,1000,50,100,250,500,1000)),
    "start_allee_thres"=list(values=c(0,0,0,0,0,20,20,20,20,20))
  )
  
)

set.seed(42) ##we set seed here to guarantee the seeds selected below are the same everytime
nl_object@simdesign<-simdesign_distinct(nl=nl_object,nseeds=10)

test=run_nl_all(nl_object)



setsim(nl_object, "simoutput")<-test

##test2=unnest_simoutput(nl_object) ##some conflicts when extricting patch and turtle variables?

output_turtles=unnest(test,cols=c(metrics.turtles)) %>% 
  select(-metrics.patches) %>% 
  mutate(replicateID=paste("K",K,"a",start_allee_thres,`random-seed`))%>% 
  mutate(Location_full=paste(replicateID,xcor,ticks))%>% 
  mutate(neutral_allele=str_pad(neutral_allele,3,pad="0"),
         neutral_allele_shuffled=str_pad(neutral_allele_shuffled,3,pad="0"))


output_patches=unnest(test,cols=c(metrics.patches)) %>% 
  select(-metrics.turtles) %>% 
  mutate(replicateID=paste("K",K,"a",start_allee_thres,`random-seed`)) %>% 
  mutate(Location_full=paste(replicateID,pxcor,ticks)) 


ttt=output_turtles %>% 
  group_by(K,start_allee_thres,replicateID,birth_patch,ticks) %>% 
  summarise(Neff=n_distinct(parentID)) %>% 
  ungroup() %>% 
  mutate(pxcor=birth_patch) %>% 
  select(-birth_patch)


output_turtles_corefront=output_turtles %>% 
  ungroup() %>% 
  group_by(ticks,replicateID) %>% 
  mutate(minfront=min(xcor),maxfront=max(xcor)) %>%
  filter(xcor==minfront |xcor==0 |xcor==maxfront) %>% 
  ungroup()
  
###some patches colonised by only 1 individual
###problematic for diversity; ultra noisy, use the two frontmost??


###isolate the genetics markers ####NB: needs to decide how null alleles are treated, esp individuals with only null alleles
genmat<- df2genind(select(output_turtles_corefront,neutral_allele),####all allelic colunns starts with P
                   ncode=3,pop=output_turtles_corefront$Location_full)
###the genmat is way too big, preprocess bu using only very front and core patches??

genmat=genind2genpop(genmat)

genmat_shuffled<- df2genind(select(output_turtles_corefront,neutral_allele_shuffled),####all allelic colunns starts with P
                   ncode=3,pop=output_turtles_corefront$Location_full)
###the genmat is way too big, preprocess bu using only very front and core patches??

genmat_shuffled=genind2genpop(genmat_shuffled)


data_genetics <- tibble(Hexp=Hs(genmat),Hexp_shuffled=Hs(genmat_shuffled),Location_full=names(Hs(genmat))) %>% 
  left_join(output_patches)







                     