!module edp_dat

	character*32 acc_dat
	data acc_dat/'acc.txt'/

	character*32 pdbstd_dat
	data pdbstd_dat/'pdbstd.txt'/

	character*(108) edp_data0, edp_data
	character*1 delimiter, wildcard, wildcard1
	integer window_size,window_size0 ,verbose, status
	data window_size0/20/

	common /cmm_data/ in_unit, verbose, lverbose, window_size
	1	,edp_data ,edp_data0, status 
	1	,delimiter ,wildcard ,wildcard1

!	common /cmm_data/ in_unit, verbose, lverbose, window_size &
!		,edp_data ,edp_data0 &
!		,delimiter ,wildcard ,wildcard1, status 

!end module edp_dat
!copyright by x. cai zhang

