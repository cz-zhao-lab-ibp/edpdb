!module edp_file

	character*(108) pdb_in, pdb_in2, pdb_in3, pdb_out
	1 	, rtn_in, rtn_out 
	1	, edp_in, edp_out, scr_out, acc_in, std_in
	1 	, kin_out, seq_out 
	1	, update_in,log_out

	character*28 cm1, cm2, cm3

	common /file_in_out/ pdb_in, cm1, ncm1 
	1	, pdb_in2, cm2, ncm2 
	1	, pdb_in3, cm3, ncm3 
	1	, pdb_out, rtn_in, rtn_out 
	1	, edp_in, edp_out, scr_out, acc_in, std_in
	1	, kin_out, seq_out 
	1	, update_in, log_out
	
	character*(64) guess(2)
	data guess/'vWkwy -- (9999), XLKBIRTSG YB c. xZR aSZMT.' 
	1 ,' zGGVMGRLM: gSRH EVIHRLM LU vWkwy RH uivv ULI ZXZWVNRXH.'/
                   
!end module edp_file
