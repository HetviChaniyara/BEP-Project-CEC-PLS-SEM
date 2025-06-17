# BEP-Project-CEC-PLS-SEM

## Description

This repository includes the implementation of Cardinality and Equality Constrained Partial Least Squares Structural Equation Modelling (CEC-PLS-SEM), a new approach aiming to tackle the issue of weight instability of PLS-SEM under high-dimensional low sample size (HDLSS) settings. The simulation study compares the performance of the new approach to Sparse Generalised Canonical Correlation Analysis (SGCCA). The repository contains the functions developed in R to implement the method and the code for the simulation study.

## Repository Structure

The repository is organised into two main folders:

-   Scripts/: Contains the R scripts for data generation and the functions required to run both methods.

    -   Psparse_Data.R : Script to obtain data with sparse factor loadings.

    -   Wsparse_Data.R : Script to obtain data with sparse weights.

    -   CEC_PLS_SEM_Functions.R : Functions for the CEC-PLS-SEM method.

    -   SGCCA_Function.R : Functions for the SGCCA method

-   Demo/ : Provides the scripts to run the method for the datasets and summarises the results for the methods.

    -   CEC_PLS_SEM.R : Script to run CEC-PLS-SEM and summarises results by averaging across 5 iterations per condition. Change the path to the data folder to obtain results under both data generation schemes.

    -   SGCCA.R : Script to run SGCCA and summarise results by averaging across 5 iterations per condition.

## How to Run This Project

-   **Step 1:** Clone the repository

-   **Step 2:** Generate the data by running Scripts/Psparse_Data.R and Scripts/Wsparse_Data.R

-   **Step 3:** Move to the demo folder and run the scripts CEC_PLS_SEM.R and SGCCA.R. Change the data directory in the files to run the scripts for both data generation schemes.

## **Acknowledgements**

This project is part of the Bachelor End Project and is submitted in partial fulfillment of the requirements of the degree of Bachelor of Science at Eindhoven University of Technology and Tilburg University under the supervision of Prof. Dr. Katrijn Van Deun. This publication is part of the project SEM2.0 (with project numbers 406.22.GO.022 and VI.C.231.092) of the Open Competition and Talent research programs financed by the Dutch Research Council (NWO), awarded to Prof. Dr. Katrijn Van Deun.

## Authors

Bachelor Project by: Hetvi Chaniyara

Under the supervision of Prof. Dr. Katrijn Van Deun
