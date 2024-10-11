## NDNF Interneurons and Sleep

___________________
### Author:
Christine Muheim ((christine.muheim@gmail.de)
___________________

#### General Info
Aim of the code is to analyze calcium activity of NDNF interneurons in layer 1 of V1 of mice over the 24h day 
and in different sleep stages. The data used to run this script contains timestamped and sleep-scored calcium recordings, 
aquired with nVoke or nVista systems from Inscopix/Bruker. 

#### Experimental Design:
NDNF-cre mice were injected with a AAV-GECI, equipped with a cranial window over the injection site and EEG/EMG electrodes. 
For details, see the Manuscript Methods section.

After adequate recovery time, animals were mounted with a baseplate to securely but reversibly attache the miniaturized microscope.

Animals were recorded over 48h, whereof the first 24hrs were undisturbed (baseline, BL), followed by a 6h sleep deprivation (SD) 
and 18hrs recovery period. SD started at ZT0 when the lights were turned on. 

___________________

#### Analysis:

Data was preprocessed based on the recommendations of Inscopix and is not part of this analysis code. 
Each ROI contains dF/Fo normalized data. Normalization of the calcium signal was within each recording session. 
The EEG data was synchronized to the timestamped calcium data to align the sleep scoring. The synchronization is not part of this code.
 
The calcium data is first separated into hours and then into sleep stages.
For each ROI, the normalized calcium signal was summed per hour and state and normalized by the number of frames in this state per hour. 

Preferred state evaluation:

Mean across the summed hourly values is calculated per vigilance state.
The vigilance state with the maximal value is defined as the "preferred" state (or -ON state).

#### Visualization
The code provided generates the same quantitave figures as the manuscript. Mot inluded is Fig 1a.


ROIs are then grouped based on this preferred state and statistically compared to the mean activity in the other vigilance states.

 
ROIs are kept in these groups. Mean activity at ZT6-12 at baseline is compared to mean activity at ZT6-12 after SD.
