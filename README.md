PLEASE READ THIS FILE CAREFULLY FIRST BEFORE RUNNING ANY OF THE CODES.

PLEASE NOTE THAT FOR NGC3031(M81), WHICH IS USED AS A TRAINING GALAXY IN THE PAPER, THE NOTATION RULE IS 'Halo Fieldxx = Field(xx+1)'. IN THIS PAGE, WE PROVIDE THE LATTER FORMAT OF FILENAMES. (e.g. Halo Field04 used in Figure 3, is actully Field05 on GHOSTS page and this GitHub/fits folder).

#description of folder

'fits' -- *.fits files grabbed from GHOSTS survey

'csv' -- *.csv files before spatial clipping

'clipped_csv' -- *.csv files after spatial clipping (we only provide results of version 4, but anyone can have results of other version or DIY his/hers via code '2_spatial_clipping.py')

‘cmd' -- *.png files showing CMDs for each galaxy before clipping (we only provide results of version 4, but anyone can have results of other version or DIY his/hers via code '2_spatial_clipping.py')

'clipped_cmd' -- *.png files showing CMDs for each galaxy after clipping (we only provide results of version 4, but anyone can have results of other version or DIY his/hers via code '2_spatial_clipping.py')

'clipping_map' -- *.png files showing clipping process for each field (we only provide results of version 4, but anyone can have results of other version or DIY his/hers via code '2_spatial_clipping.py')

'info' -- IMPORTANT! *.csv files showing basic information of each field we use, containing 'red_lim' and 'faint_lim' columns representing the red and green lines used in Figure 2, and 'slope_bc/ac' and 'inter_bc/ac' columns representing the slope and intercept of color band before clipping (bc) and after clipping (ac)

'detection' -- IMPORTANT! *.csv files showing our tip detection results of each field we use, including our defined tip properties.

PLEASE NOTE: there are version numbers in 'info' and 'detection' folders, which correspond to descriptions in 'version_code2.csv' and 'version_code3.csv'.

'plots' -- the plots in our paper and codes to generate those.

#description of data files

'coordinate.txt' -- ra, dec info, output of '1_fundamental_analysis.py'

'ghosts_extinction.txt' -- foreground extinction EBV data obtained using 'coordinate.txt', and the website for obtaining this is listed in '1_fundamental_analysis.py'

'ghosts_analysis.csv' -- output of '1_fundamental_analysis.py', including extinction, edd tip, and ra dec information of each field

'version_code2.csv' -- descriptions of parameters used in '2_spatial_clipping.py'. one can add lines in this file to DIY his/hers spatial clipping rules. the meaning of each column is stated as follows: 1)version--version number; 2)spatial_clip--'no' represents no spatial clipping, 'xxp' means clipping of xx%, which retains regions with blue star number density <= xx%, '(a simple float)' (e.g. 0.3) means the clipping retains regions with blue star number density <= that float value (0.3 for the example); 3)slope_min--the minimum of searching range for the slope of color band; 4)slope_max--the maximum of searching range for the slope of color band.
(e.g. one DIY example is '6,1.0,-7.5,-5', another is '7,15p,-6.5,-4' and enter version number 6 or 7 in '2_spatial_clipping.py')

'version_code3.csv' -- descriptions of parameters used in '3_tip_detection.py'. one can add lines in this file to DIY his/hers detection rules. the meaning of each column is stated as follows: 1)version--version number (NOTE THAT version x.y uses the spatial clipping results of version x in 'version_code2.csv'); 2)band--'yes' or 'no' stand for apply or not apply a color band; 3)smoothing--the smoothing scale; 4)weighting--choose from 'simple', 'poisson', or 'hatt'; 5)ratio_cut--set restriction for minimum R; 6)nbt_cut--set restriction for minimum N_{-, 1.0}.
(e.g. one DIY example is '5.2,no,0.08,poisson,3,200')

#description of code

'1_fundamental_analysis.py' -- turn 'FITS' into 'CSV', generate 'ghosts_analysis.csv' and 'csv'

'2_spatial_clipping.py' -- perform spatial clipping, generate 'ghosts_info_vx.csv'
			-- also drawing CMDs, clipped CMDs, and clipping maps
			-- please change version number at LINE 83 (make sure there is a corresponding line in 'version_code2.csv')

'3_tip_detection.py' -- perform tip detection, generate 'ghosts_detection_vx.y.csv'
		     -- please change version number at LINE 117 (make sure there is a corresponding line in 'version_code3.csv')

other codes in the 'plots' file --  PLEASE READ THROUGH THE CODE BEFORE RUNNING. THERE ARE INDICATORS OF WHERE TO CHANGE VERSION NUMBER. DIFFERENT VERSION NUMBER USES DIFFERENT 'info' AND 'detection' FILE, AND MIGHT GIVE COMPLETELY DIFFERENT RESULTS.
