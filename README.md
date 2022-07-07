# Analysis-rcs-data
Matlab functions and scripts to extract and visualize pain and behavioral REDCAP survey data collected alongside Summit RC+S neural data. Initial processing includes
extraction of REDCAP survey data and transformation of data into a data format for analysis, plotting, and visualization.

**Background**: UCSF teams investigating adaptive neurostimulation administered through Summit RC+S devices need an analysis and visualization framework for pain and psychological data.

**Aim**: To synthesize matlab functions and scripts for analyzing and visualizing REDCAP surveys that collect pain metrics and behavioral data.

**Collaborators**: Ryan.Leriche@ucsf.edu, Joanna.Lin@ucsf.edu, Ana.Shaughnessy@ucsf.edu, Prasad.Shirvalkar@ucsf.edu (open for more colleagues to join...)

Note: This research is funded by UH# HEAL Grant; IRB Study: G190160.

# Overview of Variables

**Independent Variables** 
- **Stimulation programs** 
   - A, B, C, vs D
   - Note that EVEN within a patient, the stim program letters are not always consistent
   - utilize metadata from 'ProcessRCS()' and/or'RCS_pain_preprocess()' to define N of unique stim program
   - stim programs can be grouped (useful to  reduce the number of comparisons made)
- **Side** 
   - L caudate, R thalamus, vs both
- **Stim parameters**        
   - polarity, pulse width, frequency, cycling, etc.
- **Days on fentanyl patch** (Dependent on patient)
   - day 1, day 2, day 3

**IJ suggestion** (see Slack)
- Left Caudate Bipolar Stim group A: Contacts 9+11- (variable amplitudes - worth looking at 0.5, 1, 1.5, 2 mA respectively) 125 hz, 300 microsecond pulse width
- Right Thalamus Bipolar Stim group A: Contacts 9+11- (variable amplitudes - worth looking at 0.5, 1, 1.5, 2 mA respectively) 125 hz, 300 microsecond pulse width
- best to compare these two regions for now and possibly expand into including some further Left caudate stim at other contacts

**Dependent Variables**
- **Daily**                   
   - NRS-intensity, NRS-worst VAS-intensity, VAS-unpleasantness, MPQ-affective, MPQ-somatic and MPQ-total)
- **Weekly**
   - AE reporting, C-SSRS, Hamd6 Self-rating Scale for Depression
- **Monthly**
   - AE reporting, C-SSRS, Beck Depression Inventory (BDI), Beck Anxiety Inventory (BAI), Rand 36 Item SF Health Survey Instrument, IMS-25, Promis SfV11 Global Health, Clinical Global Impressions Scale

Note: Assuming sufficient (>10) number of pain metrics (see Shirvalkar et al., 2020 where min of 10 trials is used to justify the duration of Stage 0 which allows generalizability of the streaming sessions indexed in this script to 
subsequent neurophy analysis.
- If this assumption fails, reconsider grouping.

# Approaches 
**Visualization Approach**
- Timeline 
- Histogram
- Select days prior


**Statistical Approach**
- **Individual**
   - RCS04: (N-way ind ANOVA, regression, etc.)
   - RCSXX: (2-way ind ANOVA - stim progam BY side)
- **Group-level**
   - TBD

# Progress Notes

| Functions to edit/create  | Person Working on it |
| ------------- | ------------- |
| Color plot background according to stage| |
| Select a time window for plot|  |
| Correlation analyses||
| | |

