# Data-Repository-for-EBMforDisTransOnRandGraphs-An-Application-to-Mitigate-a-Syphilis-Outbreak
Data Repository for the submitted article:

Edge-based Modeling for Disease Transmission on Random Graphs: An Application to Mitigate a Syphilis Outbreak

arXiv: https://arxiv.org/abs/2410.13024

Submitted and accepted by Royal Society Open Science (RSOS): https://royalsocietypublishing.org/journal/rsos

## Authorship: 
Sicheng Zhao (McMaster University, Queen's University, KFL&A Public Health)

Sahar Saeed (Queen's University)

Megan Carter (KFL&A Public Health, Queen's University)

Bradley Stoner (Centers for Disease Control and Prevention, Queen's University) 

Maggie Hoover (Queen's University, KFL&A Public Health)

Hugh Guan (KFL&A Public Health)

Felicia M.G. Magpantay (Queen's University)

## Corresponding Authors
Article: Felicia Magpantay (felicia.magpantay@queensu.ca)

Repository: Sicheng Zhao (20sz11@queensu.ca, zhaos126@mcmaster.ca)

## File Description
1. [Degree Seq 19-23 Init.csv](Degree_Seq_19-23_Init.csv): Aggregated degree sequences initially collected from data and modified/adjusted as described to article.

2. [Time Distribution 19-23_ReInfection.xlsx](Time_Distribution_19-23_ReInfection.xlsx): Aggregated and anonymized data for syphilis stages of reported cases and encounter date, together with the "recovery date" after adjustment described in the article. Used to generate the incident data and the recovery time distribution

3. [Time Test Adj.csv](Time_Test_Adj.csv): Biweekly aggregated reported casecount after adjustment for treament duration used for fitting the parameters of model.

4. [R-EpiNecPerco](R-EpiNetPerco): functions from R-EpiNetPerco used to percolation process calculation in the article

5. [Heatmap_INITFIT.R](Heatmap_INITFIT.R): .R code for using MLE to find optimized initial condition $\hat{I}_0$ for models.
    - Please contact with [Repository Author](zhaos126@mcmaster.ca) for corresponding datasets (too large for github repo).

6. [Heatmap_HR500.R](Heatmap_HR500.R): .R code for generating heatmap of MLE to find optimized $(\hat{\beta},\hat{p(\gamma})$ for models based on different initial conditions $I_0$.
    - [3QT_500_Init117_MASIR_combined.csv](3QT_500_Init117_MASIR_combined.csv): Log-likelihood results for varying parameter set from simulation of $I_0=117$ ma-SIR model.
    - [3QT_500_Init27_MASIR_combined.csv](3QT_500_Init27_MASIR_combined.csv): Log-likelihood results for varying parameter set from simulation of $I_0=27$ ma-SIR model.
    - [3QT_500_Init27_combined.csv](3QT_500_Init27_MASIR_combined.csv): Log-likelihood results for varying parameter set from simulation of $I_0=27$ network model.

7. [Heatmap_HR100.R](Heatmap_HR100.R): .R code for generating heatmap of MLE focused on neighbourhood of optimized $(\hat{\beta},\hat{p(\gamma})$ for models.
    - [3QT_HighRes100_Init117_MASIR_combined.csv](3QT_HighRes100_Init117_MASIR_combined.csv): Log-likelihood results for varying parameter set from simulation of $I_0=117$ ma-SIR model.
    - [3QT_HighRes100_Init27_MASIR_combined.csv](3QT_HighRes100_Init27_MASIR_combined.csv): Log-likelihood results for varying parameter set from simulation of $I_0=27$ ma-SIR model.
    - [3QT_HighRes100_Init27_combined.csv](3QT_HighRes100_Init27_MASIR_combined.csv): Log-likelihood results for varying parameter set from simulation of $I_0=27$ network model.

8. .R code for sensitive analysis:
    - [Heatmap_NPlus.R](Heatmap_NPlus.R): analysis of population size $N \times (1+0.05) $.
    - [Heatmap_NMinus.R](Heatmap_NMinus.R): analysis of population size $N \times (1-0.05) $.
    - [Sensitive_HR100.R](Sensitive_HR100.R): analysis of other model parameters $\beta, p, \alpha$.

9. Data for validation:
    - [Valid_data_fit.csv](Valid_data_fit.csv): 75% portion of randomly selected bi-weekly reported case counts as fitting set (n=87) for validation.
    - [Valid_data_test.csv](Valid_data_test.csv): the rest 25% portion reported case counts as testing set (n=29) for validation.

10. [Heatmap_INITFIT_valid.R](Heatmap_INITFIT_valid.R): .R code for using MLE to find optimized initial condition $\hat{I}_0$ for split network models on fitting data set (n=87).
    - Please contact with [Repository Author](zhaos126@mcmaster.ca) for corresponding datasets (too large for github repo).

11. [Heatmap_HR100_valid.R](Heatmap_HR100_valid.R): .R code for comparing split network model (n=87) with full network model (n=116) with validation on testing data set (n=29)
    - [Valid_HighRes100_Init27_combined.csv](Valid_HighRes100_Init27_combined.csv): Log-likelihood results for varying parameter set from simulation of $I_0=27$ split network model.
