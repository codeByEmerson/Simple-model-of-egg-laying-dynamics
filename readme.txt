Simple model of egg-laying dynamics

Published as part of: Short-lived neurotransmitter and long-lived neuropeptide signaling combine to pattern activity in the C. elegans egg-laying circuit
Publication Authors: Emerson R. Santiago, Michael R. Koelle*
*Corresponding author. Email: michael.koelle@yale.edu  

Code written by: Emerson R. Santiago

Goal: 
   Demonstrate that an exponentially decaying inhibitory signal +
    stochastic short-lived activating signal can combine to yield
    similar on-off dynamics as exhibited by the actual egg-laying
    circuit.
    
++++++++++++++++++++++++    
To run Model 1 wildtype
++++++++++++++++++++++++

The code is already setup to run a simulation of a wild-type egg-laying system with Model 1. This is set for publication-quality scale, so the code will take some time to resolve. Run runner_new.R from the command line as:

Rscript Model_1.R

++++++++++++++++++++++++    
To run Model 2 wildtype
++++++++++++++++++++++++   

The code is already setup to run a simulation of a wild-type egg-laying system with Model 2. This is set for publication-quality scale, so the code will take some time to resolve. Run runner_new.R from the command line as:

Rscript Model_2_wt.R 
    
++++++++++++++++++++++++++++++++++++++++    
To confirm functionality for custom runs
++++++++++++++++++++++++++++++++++++++++

The code is already setup to run a simulation of a wild-type egg-laying system. This is set for publication-quality scale, so the code will take some time to resolve. Run runner_new.R from the command line as:

Rscript runner_new.R

+++++++++++++++++++++++++++
To customize run parameters
+++++++++++++++++++++++++++

Edit the file runner_new.R with your preferred output configuration by adding/removing a line such as the following:

output<-runModel(STOPTIME=300000, Is=5.55, As=1, Ihl=3, Ad=1, Ap=0.0764, capA=TRUE, randIs=FALSE, randAs=FALSE)
analyzeRun(output)

Ensure that tMin has been set appropriately for your selected genotype. Ensure that STOPTIME matches in your output-call and when set in the file (for basic testing, STOPTIME=60000 will be sufficient and much faster).

Run runner_new.R with R from the command line:

Rscript runner_new.R

+++++++++++++++++++++++++++
Parameter names and meaning
+++++++++++++++++++++++++++

STOPTIME=300000      | How long to simulate for, 1 time step = 1 minute simulated
Is=5.55              | How much inhibitor to release at a time
As=1                 | How much activator to release at a time. Doubles as Amax in Model 1
Ihl=3                | Inhibitor half-life
Ad=1                 | Decay parameter for activator. Rate of linear decay in Model 1. Half-life decay in Model 2.
Ap=0.0764            | Odds of releasing A at a given time-step.
capA=TRUE            | Whether to cap A at As. TRUE in Model 1, FALSE in Model 2.
randIs=FALSE         | Whether to randomize level of Inhibitor released, FALSE in Model 1, TRUE in Model 2.
randAs=FALSE         | Whether to randomize level of Activator released, FALSE in Model 1, TRUE in Model 2.
expDecay=TRUE        | Whether to use exponential decay for Activator. FALSE in Model 1, TRUE in Model 2.
randN=3              | Range of variation in activator release. Unset for Model 1, Default is 6 in Model 2.



+++++++++++++++++++++++++++++++++++++++++++++++++
Scanning performance of various parameter choices
+++++++++++++++++++++++++++++++++++++++++++++++++

In addition to running single instances of the model, you can scan a range of parameter choices. The current implementation only supports scans over combinations of EITHER Activator Step Size + Activator decay (As, Ad) OR Inhibitor step size + Inhibitor half life (Is, Ihl).

By default, calling 'Rscript gridscanner.R' will run a scan of As and Ad in the wild-type model.

To run a predefined scan for a given model genotype (akin to Fig. 5 in the publication), simply open the gridscanner.R file, and uncomment your desired genotype as the base of the file. (Make sure to comment out the wild-type As/Ad scan!).

To run a custom scan, add a new call to gridscanner_AsAd (scan combinations of As/Ad) or gridscanner_IsIhl (scan combinations of Is/Ihl).

Upon conclusion of a scan, two outputs will be produced:
 - A .csv file containing a record of the scan. Further scans under the same name will re-use and extend this file, allowing a small initial scan region to be expanded, or an initial coarse scan to be made fine-grained. 
 - A .pdf file with various graphs related to the scan results. Contents of the .pdf file can be modified by altering the file 'plotParamScan.R'
 
By default, each parameter combination is run for 25000 time steps (simulated minutes). This setting can be changed in the file 'runParamScan.R' by altering the variable: STOPTIME=25000 in all places that it occurs.

The folder titled 'real_data' contains the history-dependence slopes from actual experimental data for each of the genotypes listed, corresponding to Fig. 4 in the publication.

