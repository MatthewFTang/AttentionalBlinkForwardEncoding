# AttentionalBlinkForwardEncoding
Underlying code for manuscript: Tang, M. F., Ford, L., Arabzadeh, E., Enns, J. T., Visser, T. A. W., &amp; Mattingley, J. B. (2019). Neural dynamics of the attentional blink revealed by encoding orientation selectivity during rapid visual presentation. bioRxiv, 595355. http://doi.org/10.1101/595355

The data are avaiable at : https://osf.io/f9g6h/

Readme

System requirements 
MATLAB 
Code has been tested on Mac OSX and Linux

Installation guide
For each script, change the variable ‘dataFolder’ to the location of the folder containing the code. The script should add the Functions file, which contains code from other group that is needed for the analysis. All other functions should be contained directly in the scripts themselves. The scripts are standard MATLAB code so have no installation time. 

Figure 8B requires EEGlab and Feildtrip to be on the MATLAB path


The scripts for the behavioural data (Figures 2 and 3) should run very quickly (less than 1 mins). The scripts for EEG data will take around 20 mins for each participant on an older MacBook Pro. Will be considerably faster when run on a cluster. 
