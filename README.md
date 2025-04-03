This repository contains the scripts used to analyze grid cell properties from the paper "Tiling of large-scaled environments by grid cells requires experience",
available at  https://www.biorxiv.org/content/10.1101/2025.02.16.638536v1.

Data on neural spiking were collected from 4 mice (mouse1 = fuego, mouse2 = huevo, mouse3 = samwise, mouse4 = lima).  
Data are organized by sessions as CMBHOME objects in MATLAB. Each CMBHOME object contains spikes from 4 tetrodes sorted as single units. 
For complete information, see the wiki at https://github.com/wchapman/CMBHOME/wiki

File name organization: 
Each file name contains the name of the mouse, the date of the recording, the duration of the recording (in minutes), 
the arena dimension, the height of the walls, and the session number if the session was repeated.

Example:
20230831_25min_45x45cmHighWalls_ses2 indicates that the data were acquired on Aug 31, 2023, in a recording that lasted 25 min, 
where mice explored and arena with a floor size of 45 x 45 cm2 and 30 cm high walls. This was the second recording in this arena on that day.
 
Code was custom-written and adapted from toolboxes from The Hasselmo Lab: https://github.com/hasselmonians/CMBHOME and the Moser group: https://bitbucket.org/cnc-ntnu/bnt/src/master/. 
Custom-written code is available on https://github.com/dannenberglab/Tiling-of-large-scaled-environments-by-grid-cells.

Installation:

Download this repository.

Download dataset https://doi.org/10.13021/ORC2020/OJT2JL

Download BNT from Moser group 

Download CMBHOME 

Main script: "callGrid" analyzes data from two recording sessions (45x45 small maze) and (60x60 medium maze) 
create example figures of all analysis performed in the paper and computing grid cell properties and spatial crosscorrelation. 

Please reference the following citation when using this dataset in your work:
Gutiérrez-Guzmán BE, Hernández-Pérez JJ, Dannenberg H (2025) Tiling of large-scaled environments by grid cells requires experience.



