# ComboPred: prioritizing synergistic, selective, and clinially actionable mutli-targeting drug combinations for individual patients.

## Introduction
ComboPred is part of drug combination prediction and testing (DCPT) platform developed by computational systems medicine group at FIMM, University of Helsinki (<https://www.fimm.fi/en/research/groups/aittokallio>). ComboPred combines the ex vivo sensitivity profiles of the patient’s response to the single compounds, with exome-seq and RNA-seq profiles of the same patient, through a systems-wide compound-target interaction network from our crowdsourcing bioactivity data platform, DrugTargetCommons (DTC) (<https://drugtargetcommons.fimm.fi/>) to predict combination response. With the predicted combination response, ComboPred calculates combination synergy score with highest single agent (HSA) model.  

![](man/figures/combopred.png)

## Usage
Data preparation: drug-target interaction matrix, gene expression data for each individual sample, mutation data for each individual sample, drug responses (DSS).  Example data can be found in <https://github.com/hly89/ComboPred/tree/master/data>. And the results are saved csv files.

Run the following command to predict combination response and obtain synergy score and results are save as csv files:
```
ComboPred (gex, mut, dtm,
                      model.iteration = 6, response, 
                      patient.index = 1, 
                      control.index = 1)
```
## Support 
If you have any problems in using ComboPred, please contace **Liye He** (liye.he@helsinki.fi)

## Copyright & License

Code released under the MIT license.

