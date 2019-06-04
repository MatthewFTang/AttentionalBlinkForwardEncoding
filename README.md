# AttentionalBlinkForwardEncoding
Underlying code for manuscript: Tang, M. F., Ford, L., Arabzadeh, E., Enns, J. T., Visser, T. A. W., &amp; Mattingley, J. B. (2019). Neural dynamics of the attentional blink revealed by encoding orientation selectivity during rapid visual presentation. bioRxiv, 595355. http://doi.org/10.1101/595355


Readme

MATLAB code for analysing ‘Neural dynamics of the attentional blink revealed by encoding orientation selectivity during rapid visual presentation’.

System requirements 
MATLAB 2018a or 2018b
Code has been tested on Mac OSX 10.13 and 10.14

Installation guide
For each script, change the variable ‘dataFolder’ to the location of the folder containing the code. The script should add the Functions file, which contains code from other group that is needed for the analysis. All other functions should be contained directly in the scripts themselves. The scripts are standard MATLAB code so have no installation time. 

Demo.
To run each script, after changing the variable dataFolder, press F5 or run and it should produce the figures from the manuscript but with a smaller number of participants. 

The scripts for the behavioural data (Figures 2 and 3) should run very quickly (less than 2 mins). The scripts for EEG data will take around 20 mins for each participant on a MacBook Pro. Will be considerably faster when run on a cluster. 
