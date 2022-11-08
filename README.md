# Source code and data to mimic results in Section 4 of the paper: "Latency function estimation under the mixture cure model when the cure status is available"


**Authors:**

Wende Clarence Safari, University of A Coruña, E-mail: wende.safari@udc.es

Ignacio López-de-Ullibarri, University of A Coruña, E-mail: ilu@udc.es

M. Amalia Jácome, University of A Coruña, E-mail: majacome@udc.es


## Content 

An outline of the file structure, which details the R scripts and data, is:

1.- case_study_analysis.R

The code in this script describe and analyze Section 4 of the paper. The script, 
    
     + Loads all necessary libraries.
     + Sources function file. Check out the functions_definitions.R script to see all the functions called. 
     + Loads the .csv simulated datasets.
     + Describe the data and perform the analysis.
    	     

2.- functions_definitions.R

Within case_study_analysis.R script, we call for a script functions_definitions.R that loads functions needed in the computations of 

  + the proposed estimator for the conditional survival function, the cure probability estimator and the latency estimator.
  +  the bootstrap bandwidths computations for the proposed estimator.
   
  
3.- simulated_data.csv 

This file contains the simulated dataset.
