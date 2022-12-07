# Pain Metric Analysis
MATLAB functions and scripts to extract and visualize pain and behavioral REDCAP survey data collected alongside Summit RC+S neural data. Initial processing includes
extraction of REDCAP survey data and transformation of data into a data format for analysis, plotting, and visualization.

**Background**: UCSF teams investigating adaptive neurostimulation administered through Summit RC+S devices need an analysis and visualization framework for pain and psychological data.

**Aim**: To synthesize MATLAB functions and scripts for analyzing and visualizing REDCAP surveys that collect pain metrics and behavioral data.

**Collaborators**: Ryan.Leriche@ucsf.edu, Joanna.Lin@ucsf.edu, Ana.Shaughnessy@ucsf.edu, Prasad.Shirvalkar@ucsf.edu (open for more colleagues to join...)

Note: This research is funded by UH# HEAL Grant; IRB Study: G190160.


Currently:
| Functions to edit/create  | Person |
| ------------- | ------------- |
| RCS07 plot_timeline() | Ryan|
| stim groups by contacts-freq-amp-PW-cyc to motivate aDBS settings| Ryan |
| incrp RCS02's textlogs into stim group changes via 'RCS_logs()'| Ryan|
| ensuring quality of RCS neural data via rcsPlotter| Joanna |
| statistical framework-- Bayer factor for null comparisons | Joanna & Ryan|

Up Next:
| Functions to edit/create  | Person |
| ------------- | ------------- |
| merge RCS pt visit dates as metadata call(from PNLprasad) into PNL (this) branch  |  |
| have plotting fxns use makedatebaseRCS (from openmind-consortium/prasad) rather than makedatabaseRCS_Ryan |   |
| 'plot_timeline()' with stim parameters as background colors (wait for 'wrt_stim_REDcap' struc (how we format pain reports to given stim settings) to be finalized| |
| 'pain_versus()' residual plots for correlations| |
| 'pain_versus()' linear mixed model w/ pts as random variables (group-analysis of relationships btwn pain metrics)| |
