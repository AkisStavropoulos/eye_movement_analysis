# eye_movement_analysis
Code for eye movement analysis of human firefly task data.
It has been tested and can run in any version of MATLAB between 2018 and 2022.
No specific hardware requirements needed.



**In order to run:**

1) Full and demo datasets are uploaded on https://gin.g-650 node.org/akis_stavropoulos/belief_embodiment_through_eye_movements_facilitates_memory-651 guided_navigation.

- To run the **full dataset**: download EyeMovementDataset.mat.

- To run the **demo dataset**: download EyeMovementDataset_demo.mat

2) Download all content of the repo.

3) Add to path both the main_repo and helper_functions folders.

4) You can find the main script under /main_repo with the name EyesEmbodyInternalBeliefs_code.m.

5) In line 9 of that script, change the data folder path to correspond to the folder you saved the dataset in (demo or full dataset), e.g., data_folder = 'C:\Users\me\Data\';

6) Run this script and the data analysis will start and the figures shown in the manuscript will be produced and should look identical.

* It should take less than 5 minutes for the whole code to run (including loading of dataset).

