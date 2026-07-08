# To calculate the corrected binding free energies
## Step 1 prepare the PPE inputs: PPE = post-process energy correction
 - modify path in ```PPEprep30_px_631.sh```
 	* MM_path --> to obtain `TS` value from MM run
	* PPE_path --> to extract FMO energy
	* hope --> Harmonic oscilator (H.O.) constant
 - output, for example 2i0d_ac_67_PPE
## Step 2 determine the corrected binding FE
 - using calc_BFE.sh
## Verification: compare with binding_FE_fullprotein.csv
 the binding FE from this procedure should match qmvm2_feprocess column
