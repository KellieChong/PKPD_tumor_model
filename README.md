# PKPD_tumor_model
## Abstract
Non Small Cell Lung Cancer is a highly malignant form of lung cancer that takes the lives of over a million people annually, worldwide. It is caused by mutations in multiple genes, including the epidermal growth factor receptor and tumor suppressor p53 gene. Individuals who smoke are at a high risk of developing Non Small Cell Lung Cancer. Current treatment strategies include surgery to remove the tumor mass in early stages, coupled with radiotherapy, immunotherapy, and/or chemotherapy therapy to shrink and kill the tumor. Due to the heterogeneity of the disease, patient response to treatment varies, making it extremely difficult to formulate a generalized treatment plan. The use of mathematical models that can predict tumor growth dynamics as a function of patient response can help alleviate this issue and improve patient prognosis. In this paper we aim to reconstruct the minimal pharmacodynamic and pharmacokinetic compartment model described in Ionescu et al.’s paper. Followed by the extension of introducing a time-dependent necrosis and level of effect decay to accurately account for drug trapping and tumor mutation, leading to decreased drug efficacy over time. Further, preliminary results including immune responses as a function of tumor size and radiotherapy dose from the extended model will be discussed in the scope of planning for optimal radio-immunotherapy timing and treatment regimens.

## PK/PD Model Simulation Usage 

The file project.m can be used to produce the plots and simulations from Ghita et al.'s (2022) paper. The proposed system is shown in the figure below.
The first two blocks describe the immunotherapy, and chemotherapy (anti-angiogenesis) dose response curves over time in days. An equivalent figure is produced from the second block for radiotherapy.
![image](https://github.com/KellieChong/PKPD_tumor_model/assets/56052705/a631b67c-faf2-4c2f-a4fe-3abcafd974ad)


The following blocks are used to solve the system of ODEs for the original model and the extended model.

## Immune Response Model
The linear quadratic immune repsonse model was originally proposed by Cho et al. (2023) The solutions to the system of ODEs from Ghita et al's model is used as an input to the immune response model and can be found in the results folder. The immune response over the simulation time span is plotted in the resultant figure along with the model's total tumor volume.

