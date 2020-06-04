# An individual-based model of range expansions

A Netlogo model written for my projects on the ecological and evolutionary dynamics of range expansions.

At the moment (v0.5, March 2020), it allows for both neutral genetic diversity and dispersal to evolve, but the latter is not correctly implemented yet and can (and should) be turned off for most applications. A better implementation of phenotypic evolution will come later

For details on how the model is set up: (excerpts of the Info pane content in-model)

## WHAT IS IT?

This is a general spatial range expansion model, designed to study what leads to pushed versus pushed range expansions, and the possible evolutionary consequences.

This model is designed to operate at "low" densities (equilibrium population size < 1000), at least several orders of magnitude lower than the ones classically used in the theoretical literature, which are adapted to microbial species, but may not be applicable to macroscopic species with lower population sizes.

From an initial population, haploid individuals reproduce, compete and disperse and as a result the species colonizes the landscape.

The model allows for density-dependency in dispersal and growth (including Allee effects), as well as individual variation in dispersal and growth.

When there is individual variation in traits, evolution can be activated (individuals inherit their parental allele) or turned down (trait values are drawn at random every generation).

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

### Set-up phase
The model is initiated by placing K adult individuals in the patch of coordinates pxcor 0 and pycor 0, and setting their phenotypic traits from the distribution described by parameters, and by drawing the allele at the neutral locus from a Bernoulli distribution with p = 0.5.

### Go phase
Individuals then live the following life cycle:

-once adults, they disperse or not based on their dispersal-density reaction norm and current patch population size

-they reproduce clonally and transmit their neutral allele (and in the evolutionary setting trait alleles to offspring
the fecundity formula directly gives the number of offspring post-competition and accounting for Allee effects in one step, to avoid creating individuals that would then be killed

-they die


## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

### Global parameters

The following global parameters can be set up in the Interface tab:

- *trait_variation*: two settings; whether or not phenotypic traits are re-drawn at random at each generation, blocking evolutionary change
- *K*: (average) carrying capacity/ equilibrium population density
- *duration*: duration of a run, in generations
- *disp0_mean*: average dispersal rate at the start of the run
- *slope_disp_mean*: slope of the relationship between relative density (population size/K) and dispersal rate **on the logit scale**
- *logit_disp0_sd*: standard deviation of the distribution of dispersal rates at the start of the run **on the logit scale**
- *slope_disp_sd*: standard deviation of the distribution of dispersal reaction norm slopes at the start of the run **on the logit scale**
- *fecundity*: average fecundity at the start of the run
- *start_allee_thres*: Allee threshold at the start of the run (always >= 0; 0: no Allee effect, >1: strong Allee effect, between 0 and 1: weak Allee effect)

### Landscape window

The landscape window shows the progression of the current wave through the landscape. The whiter the patch, the larger the population, empty patches are black.

### Report graphs

Because the model is designed to be analysed through R (package nlrx), there is only one graph in the Interface tab, showing the total metapopulation size through time. But more can easily be added
(suggestions: front position (or maximal x coordinate of patch with at least 1 individual) through time
Genetic diversity (formula for expected heterozygosity when 2 alleles = 2pq where p and q are allelic proportions) in front patch through time)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)

written by Maxime Dahirel, from an initial Matlab model by Marjorie Haond
