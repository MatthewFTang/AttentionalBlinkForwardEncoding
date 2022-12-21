# AttentionalBlinkForwardEncoding
Underlying code for : Tang, M. F., Ford, L., Arabzadeh, E., Enns, J. T., Visser, T. A. W., & Mattingley, J. B. (2020). Neural dynamics of the attentional blink revealed by encoding orientation selectivity during rapid visual presentation. **Nature Communications**, 11(1), 434. http://doi.org/10.1038/s41467-019-14107-z

# Abstract 

The human brain is inherently limited in the information it can make consciously accessible. When people monitor a rapid stream of visual items for two targets, they typically fail to see the second target if it occurs within 200–500 ms of the first, a phenomenon called the attentional blink (AB). The neural basis for the AB is poorly understood, partly because conventional neuroimaging techniques cannot resolve visual events displayed close together in time. Here we introduce an approach that characterises the precise effect of the AB on behaviour and neural activity. We employ multivariate encoding analyses to extract feature-selective information carried by randomly-oriented gratings. We show that feature selectivity is enhanced for correctly reported targets and suppressed when the same items are missed, whereas irrelevant distractor items are unaffected. The findings suggest that the AB involves both short- and long-range neural interactions between visual representations competing for access to consciousness.


## Usage guide 

The data are available [here](https://osf.io/f9g6h/) <br>

### System requirements 
MATLAB <br>
Code has been tested on Mac OSX and Linux

### Getting started
Modify the dataFolder variable in the Setup.m script to give the location of the code file. By default, the data should be located in the same directory as the code, but you can change this by modifying the **SetUp.m** script. The script should add the Functions file, which contains code from other group that is needed for the analysis. All other functions should be contained directly in the scripts themselves.

**Figure 8B** requires EEGlab and Fieldtrip to be on the MATLAB path


The scripts for the behavioural data (Figures 2 and 3) should run very quickly (less than 1 min). The scripts for EEG data will take around 20 mins for each participant on an older MacBook Pro. Will be considerably faster when run on a cluster. 
