# eye_movement_analysis
Code for eye movement analysis of human firefly task data.
It has been tested and can run in any version of MATLAB between 2018 and 2022.
No specific hardware requirements needed.



**In order to run:**

1) To run the full dataset: download full dataset EyeMovementDataset.mat from https://gin.g-650 node.org/akis_stavropoulos/belief_embodiment_through_eye_movements_facilitates_memory-651 guided_navigation.
Otherwise, run the demo dataset under /demo_dataset\EyeMovementDataset_demo.mat, that includes data from one subject.

3) Download all content of the repo.

4) Add to path both the main_repo and helper_functions folders.

5) You can find the main script under /main_repo with the name EyesEmbodyInternalBeliefs_code.m.

6) In line 9 of that script, change the data folder path to correspond to the folder you saved the dataset in (demo or full dataset), e.g., data_folder = 'C:\Users\me\Data\';

7) Run this script and the data analysis will start and the figures shown in the manuscript will be produced and should look identical.

* It should take less than 5 minutes for the whole code to run (including loading of dataset).

