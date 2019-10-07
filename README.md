# cm-dag
Code for the compartmental model/DAG project: Havumaki J and Eisenberg MC, 2019. [Mathematical modeling of directed acyclic graphs to explore competing causal mechanisms underlying epidemiological study data](https://www.medrxiv.org/content/10.1101/19007922v1), medRxiv preprint.

The age-structured obesity paradox model can be ran in R using the model2_main.R file.

The code will sample parameter values using Latin Hypercube Sampling and simulate a yearlong cohort study for each parameter set. Once each study is completed, a simulated dataset is created (by caculating person-time and number of deaths for each disease state) and mortality rate ratios (comparing normal weight to obese) for ever-smokers and never-smokers are calculated.

The obesity paradox occurs when, for a given simualated study (sampled parameter set), normal weight ever-smokers have higher mortality rates than their obese counterparts (i.e., mortality rate ratio >1) AND normal weight never-smokers have lower mortality rates than their obese counterparts (i.e., mortality rate ratio <1).
