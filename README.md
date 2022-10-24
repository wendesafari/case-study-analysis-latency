# Source code and data for the paper: "Latency function estimation under the mixture cure model when the cure status is available"


**Authors:**

Wende Clarence Safari, University of A Coruña, E-mail: wende.safari@udc.es

Ignacio López-de-Ullibarri, University of A Coruña, E-mail: ilu@udc.es

M. Amalia Jácome, University of A Coruña, E-mail: majacome@udc.es


## Content 

An outline of the files structure, which details the R scripts and data, is:

1.- case_study_analysis.R

The code in this script is used to analyze the simulated data that mimics the COVID-19 data described and analyzed in results of the proposed estimator in Figure 3
in the paper. The script, 
    
     + Loads all necessary libraries.
     + Sources function file. Check out the functions_definitions.R script to see all the functions called. 
     + Loads the .csv simulated datasets.
    	     

2.- functions_definitions.R

Within brc_data_analysis.R script, we call for a script functions_definitions.R that loads functions needed in the computations of 

	+ the proposed estimator for the conditional survival function, the proposed cure probability estimator and the proposed latency estimator.
  +  the bootstrap bandwidths computations for the proposed estimator


   
  
3.- simulated_data.csv 

This file contains the simulated dataset.
