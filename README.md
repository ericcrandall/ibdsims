---
title: "README.md"
output: html_document
---
# Simulating and Analyzing Marine Isolation By Distance Models
This is an open notebook of my endeavors to simulate gene flow and demography in a marine system, and then analyze the results with various methods to see how they perform. I am using [IBDsim](http://www1.montpellier.inra.fr/CBGP/software/ibdsim/) for the simulations.

#####9/5/16
To start with, the idea is to simulate a "coastline" of 100 demes, with sampling occurring every 10 demes of 10,000 effective individuals. 10,000 generations ago, the species existed in a restricted range of 10 demes each with 1000 individuals before expanding to present size with logistic expansion (r = 0.3; Jensen et al 2012).  Dispersal for each generation occurs with a zeta distribution (the discrete version of a pareto distribution) where the probability of moving k steps is $$ f(x) = M / k^n $$ where M is the proportion of emigrants (I have a 2 year generation time, with adults reproducing, so 0.5 individuals creating gametes = M) and n is the zeta shape parameter (higher values means tigher distribution and less kurtosis). 

I am also going to try a variety of markers, including mtDNA, msats and SNPS.  Tonight I simulated two datasets with these parameters, one with msats and one for mtDNA. The SNP simulation didn't work, possibly because my computer choked when trying to simulate 5000 SNPS? It works for 10 SNPS, but not 100 or 1000.

Anyway, here are the simulation parameters, with some notes about choices.

```
%%%%%%%%%% SIMULATION PARAMETERS %%%%%%%
Data_File_Name = pareto2
GenePop_File_Extension = .gen  #needed to be read by adegenet
Run_Number = 10
Random_Seeds = 12345
Pause = Never
Genepop = True
Migrate = True
Migraine_Settings = True
Nexus_File_Format = Haplotypes_and_Individuals


%%%%%%%%%% MARKERS PARAMETERS %%%%%%%%%
Locus_Number = 10
Mutation_Rate = 0.0005  #standard literature value for msat mutation rate
Mutation_Model = GSM
Allelic_Lower_Bound = 1
Allelic_Upper_Bound = 200
Allelic_State_MRCA = 90
Repeated_Motif_Size = 3
Geometric_Variance_In_GSM = 0.36 #the default for generalized stepwise model (GSM)

%%%%%%%%%% VARIOUS COMPUTATION OPTIONS %%%
DiagnosticTables = Iterative_Identity_Probability, Hexp, Fis, Seq_stats, Prob_Id_Matrix, Effective_Dispersal, Iterative_Statistics, Allelic_Variance,

%%%%%%%%%% DEMOGRAPHIC OPTIONS %%%%%%%%
Ploidy = Diploid
Lattice_Boundaries = Absorbing
Total_Range_Dispersal = True
Lattice_SizeX = 100
Lattice_SizeY = 1
Ind_Per_Pop = 10000
Void_Nodes = 1
Sample_SizeX = 10
Sample_SizeY = 1
Min_Sample_CoordinateX = 5  #so sampling doesn't happen at the edge, 
Min_Sample_CoordinateY = 1
Ind_Per_Pop_Sampled = 20
Void_Sample_Node = 10 # it will happen at 5, 15, 25...
Dispersal_Distribution = Pareto  #this is actually a zeta distribution (discrete analog of Pareto)
Immigration_Control = Simple_1D_Product
Total_Emigration_Rate = 0.5
Dist_max = 100
Pareto_Shape = 2 #total emigration rate of 0.5 and pareto shape of 2 means 1% of individuals from a given deme will make it 10 demes downstream.
Continuous_Deme_Size_Variation = Logistic
Logistic_Growth_Rate = 0.3
Continuous_Lattice_Size_Variation = Logistic
Lattice_Logistic_Growth_Rate = 0.3


% 1 DEMOGRAPHIC CHANGE
New_Demographic_Phase_At = 10000
Lattice_SizeX = 10
Ind_Per_Pop = 1000
```
## 9/9/16 
Submitted three rounds of three replicate migrate runs for the msat dataset. Each round contains seven models, a stepping stone model, a single parameter stepping stone and an island model all 10 sampled demes as well as demes lumped into twos.

```
[sgulamhu@campusrocks3 ~]$ cd migrate_ibdsim4
[sgulamhu@campusrocks3 ~/migrate_ibdsim4]$ ../submitter_all_ibd.sh
10island/
Your job 1345 ("migrate") has been submitted
10stepping.stone/
Your job 3217 ("migrate") has been submitted
10stepping.stone1/
Your job 1347 ("migrate") has been submitted
5island/
Your job 1348 ("migrate") has been submitted
5stepping.stone/
Your job 1349 ("migrate") has been submitted
5stepping.stone1/
Your job 1350 ("migrate") has been submitted
panmixia/
Your job 1351 ("migrate") has been submitted
done submitting 7 runs
[sgulamhu@campusrocks3 ~/migrate_ibdsim4]$ cd ../migrate_ibdsim5
[sgulamhu@campusrocks3 ~/migrate_ibdsim5]$ ../submitter_all_ibd.sh
10island/
Your job 3219 ("migrate") has been submitted
10stepping.stone/
Your job 1353 ("migrate") has been submitted
10stepping.stone1/
Your job 2553 ("migrate") has been submitted
5island/
Your job 3218 ("migrate") has been submitted
5stepping.stone/
Your job 1356 ("migrate") has been submitted
5stepping.stone1/
Your job 1853 ("migrate") has been submitted
panmixia/
Your job 1358 ("migrate") has been submitted
done submitting 7 runs

[sgulamhu@campusrocks3 ~/migrate_ibdsim6]$ ../submitter_all_ibd.sh
10island/
Your job 1857 ("migrate") has been submitted
10stepping.stone/
Your job 2983 ("migrate") has been submitted
10stepping.stone1/
Your job 1361 ("migrate") has been submitted
5island/
Your job 1362 ("migrate") has been submitted
5stepping.stone/
Your job 1363 ("migrate") has been submitted
5stepping.stone1/
Your job 1364 ("migrate") has been submitted
panmixia/
Your job 1365 ("migrate") has been submitted
done submitting 7 runs

[sgulamhu@campusrocks3 ~/migrate_singlepop1]$ ../submitter_all_ibd.sh
10island/
Your job 2556 ("migrate") has been submitted
10stepping.stone/
Your job 3035 ("migrate") has been submitted
10stepping.stone1/
Your job 2558 ("migrate") has been submitted
5island/
Your job 2559 ("migrate") has been submitted
5stepping.stone/
Your job 2560 ("migrate") has been submitted
5stepping.stone1/
Your job 2561 ("migrate") has been submitted
panmixia/
Your job 2562 ("migrate") has been submitted
done submitting 7 runs

[sgulamhu@campusrocks3 ~/migrate_mtdna_ibdsim1]$ ../submitter_all_ibd.sh
10island/
Your job 3010 ("migrate") has been submitted
10stepping.stone/
Your job 3011 ("migrate") has been submitted
10stepping.stone1/
Your job 3012 ("migrate") has been submitted
5island/
Your job 3013 ("migrate") has been submitted
5stepping.stone/
Your job 3014 ("migrate") has been submitted
5stepping.stone1/
Your job 3015 ("migrate") has been submitted
panmixia/
Your job 3016 ("migrate") has been submitted
done submitting 7 runs

[sgulamhu@campusrocks3 ~/migrate_mtdna_ibdsim1]$ ../submitter_all_ibd.sh
10island/
Your job 3010 ("migrate") has been submitted
10stepping.stone/
Your job 3011 ("migrate") has been submitted
10stepping.stone1/
Your job 3012 ("migrate") has been submitted
5island/
Your job 3013 ("migrate") has been submitted
5stepping.stone/
Your job 3014 ("migrate") has been submitted
5stepping.stone1/
Your job 3015 ("migrate") has been submitted
panmixia/
Your job 3220 ("migrate") has been submitted
done submitting 7 runs

[sgulamhu@campusrocks3 ~/migrate_mtdna_ibdsim2]$ ../submitter_all_ibd.sh
10island/
Your job 3018 ("migrate") has been submitted
10stepping.stone/
Your job 3019 ("migrate") has been submitted
10stepping.stone1/
Your job 3020 ("migrate") has been submitted
5island/
Your job 3021 ("migrate") has been submitted
5stepping.stone/
Your job 3022 ("migrate") has been submitted
5stepping.stone1/
Your job 3023 ("migrate") has been submitted
panmixia/
Your job 3221 ("migrate") has been submitted
done submitting 7 runs

[sgulamhu@campusrocks3 ~/migrate_mtdna_ibdsim3]$ ../submitter_all_ibd.sh
10island/
Your job 3026 ("migrate") has been submitted
10stepping.stone/
Your job 3027 ("migrate") has been submitted
10stepping.stone1/
Your job 3028 ("migrate") has been submitted
5island/
Your job 3029 ("migrate") has been submitted
5stepping.stone/
Your job 3030 ("migrate") has been submitted
5stepping.stone1/
Your job 3031 ("migrate") has been submitted
panmixia/
Your job 3222 ("migrate") has been submitted
done submitting 7 runs

Single Panmictic population run for mtDNA
sgulamhu@campusrocks3 ~/migrate_mtdna_singlepop]$ ../submitter_all_ibd.sh
10island/
Your job 3230 ("migrate") has been submitted
10stepping.stone/
Your job 3231 ("migrate") has been submitted
10stepping.stone1/
Your job 3232 ("migrate") has been submitted
5island/
Your job 3233 ("migrate") has been submitted
5stepping.stone/
Your job 3234 ("migrate") has been submitted
5stepping.stone1/
Your job 3235 ("migrate") has been submitted
panmixia/
Your job 3236 ("migrate") has been submitted
done submitting 7 runs




```
