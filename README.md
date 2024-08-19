# Bayesian Sample Size Re-estimation for Response Adaptive Randomized Trials

This code is to accompany "Sample Sizes for Randomized Controlled Trials Utilizing Bayesian Response Adaptive Randomization for Continuous Outcomes" by Aslanyan et. al. 

Change directories as needed, but keep all the files in the same directory.

## Getting Started

Clone this repo into your  desired directory

```
git close git@github.com:vahanaslanyan/bayesian_sse_rar.git
```
## Figures 1 and 2
1. Execute `simulations_delta.R`. This will give create 8 csv files with names starting with `changing_delta*`. 

2. Execute `simulations_eta.R`, `simulations_zeta.R`, and `simulations_xi.R`. These will create the csv files with names starting with `changing_eta*`, `changing_zeta*`, and `changing_xi*`, respectively. 

3. Once all the code is executed, execute `plots.R` and it will produce create and export Figures 1 and 2 in svg format. 

Data for NHST scenarios can be provided upon request. There is no code for these because they were generated using [PASS software](https://www.ncss.com/software/pass/). 

## Table 2

4. To get the values for Scenarios 1-5, execute `simulations_markdown.qmd`. This will produce a pdf with the results.

5. To get the values for Scenarios 6-10, execute `simulations_markdown610.qmd`. This will produce a pdf with the results. 

6. To get the values for Scenarios 11-15, execute `simulations_markdown1115.qmd`. This will produce a pdf with the results. 

7. To get the values for Scenarios 16-20, execute `simulations_markdown_lecanemab.qmd`. This will produce a pdf with the results. 

8. To get the values for Scenarios 21-25, execute `simulations_markdown_Tideglusib.qmd`. This will produce a pdf with the results. 

