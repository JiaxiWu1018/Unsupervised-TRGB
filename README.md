#description of folder

'fits' -- *.fits files grabbed from GHOSTS survey

'csv' -- *.csv files before spatial clipping

'clipped_csv' -- *.csv files after spatial clipping

‘cmd' -- *.png files showing CMDs for each galaxy before clipping

'clipped_cmd' -- *.png files showing CMDs for each galaxy after clipping

'clipping_mag' -- *.png files showing clipping process for each field


#description of data files

'coordinate.txt' -- ra, dec info, output of '1_fundamental_analysis.py'

'ghosts_extinction.txt' -- foreground extinction EBV data obtained using 'coordinate.txt'

'ghosts_analysis.csv' -- output of '1_fundamental_analysis.py', including extinction, edd tip, and ra dec information of each field


#description of code

'1_fundamental_analysis.py' -- turn 'FITS' into 'CSV', generate 'ghosts_analysis.csv' and 'csv'

'2_spatial_clipping.py' -- perform spatial clipping, generate 'ghosts_info_vx.csv'
			-- also drawing CMDs, clipped CMDs, and clipping maps
			-- change parameters at line 47 and 58

'3_tip_detection.py' -- perform tip detection, generate 'ghosts_detection_vx.y.csv'
		     -- change parameters at line 120, 138 and 145


#description of ghosts_info version -- output of '2_spatial_clipping.py'

'ghosts_info_v1.csv'
	a copy from ghosts_analysis, no spatial clipping and band information

'ghosts_info_v2.csv'
	only slope_bc, inter_bc columns (bc means before clipping)
	10% spatial clipping and [-7,0] band searching version
	generate clipped_csv/{filename}_10p.csv

'ghosts_info_v3.csv'
	contain slope_ac, inter_ac columns (ac means after clipping)
	same clipping and searching criteria as v2, but another band searching after spatial clipping
	apply the second band searching result as color cut


#description of ghosts_detection version -- output of '3_tip_detection.py', vx.y using info_vx data

'ghosts_detection_v1.1.csv'
	ratio_min, nbt_min, peak_min = 3, 0, 0.6
	sp_clip, band = False, False
	smooth, weighting = 0.1, 'hatt'
	ratio_cut, nbt_cut = 0, 0

'ghosts_detection_v2.1.csv'
	ratio_min, nbt_min, peak_min = 3, 0, 0.6
	sp_clip, band = True, False
	smooth, weighting = 0.1, 'hatt'
	ratio_cut, nbt_cut = 0, 0

'ghosts_detection_v3.1.csv'
	ratio_min, nbt_min, peak_min = 3, 0, 0.6
	sp_clip, band = True, True
	smooth, weighting = 0.1, 'hatt'
	ratio_cut, nbt_cut = 0, 0

'ghosts_detection_v3.2.csv'
	ratio_min, nbt_min, peak_min = 3, 0, 0.6
	sp_clip, band = True, True
	smooth, weighting = 0.1, 'hatt'
	ratio_cut, nbt_cut = 1.5, 50

'ghosts_detection_v3.3.csv'
	ratio_min, nbt_min, peak_min = 3, 0, 0.6
	sp_clip, band = True, True
	smooth, weighting = 0.1, 'hatt'
	ratio_cut, nbt_cut = 4, 0

'ghosts_detection_v3.4.csv'
	ratio_min, nbt_min, peak_min = 3, 0, 0.6
	sp_clip, band = True, True
	smooth, weighting = 0.1, 'hatt'
	ratio_cut, nbt_cut = 4, 100
