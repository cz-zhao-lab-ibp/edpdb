! If you are working on a Lunix or Unix system, after modifying this file, 
! type "touch *.f" before using "make" to recompile the program. 

!module edp_dim

	parameter (max_atom=40000)
!	parameter (max_atom=100000)
	parameter (max_res=8000)	
!	parameter (max_res=15000)	
	parameter (max_symm=48)
!	parameter (max_symm=24)
	parameter (max_gr=40)
!	parameter (max_gr=20)
	parameter (max_vector=20)
	parameter (max_num_zones=32)
	parameter (max_num_chars=132)

	parameter (m_main_para=31)

	parameter (max_param=64, max_reserved_param=8)
	parameter (max_udks=10)

	parameter (max_map=256)		! volume calculation
	parameter (max_l=2400)		! mcs clique search
	parameter (max_lv=96)		! match3d clique search
	
	parameter (max_m=20)		! max_num_main_chain_atoms

	parameter (max_rt=40, max_r_a= 300) !vdw radii

!end module edp_dim
