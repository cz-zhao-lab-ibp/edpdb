!module edp_main
!	use edp_dim
	include 'edp_dim.inc'

	common /coor/ n_atom,lf(max_atom),atom(max_atom)
	1 , n_res,res(max_res),ires(max_res),ijk(max_res) 
	1 , x(max_atom),y(max_atom),z(max_atom),w(max_atom),b(max_atom) 
	1 , iorder(max_atom), aa_seq(max_atom)
	integer n_atom, n_res, ijk, iorder, aa_seq
	real x,y,z,w,b
	logical*1 lf
	character ires*5,atom*4,res*3 

!see subroutinte readf()
	parameter (max_jnk= 256*72/4)
	common /cmm_text_string/ text, jnk_text
	character*(72) text(max_atom)
	integer jnk_text(max_jnk)

	common /error_msg/ errmsg, max_err, n_err, syntax, n_of_syn
	character*(max_num_chars) errmsg
	character*(max_num_chars) syntax(32)
	integer max_err, n_err, n_of_syn

	common /initial/ zone_untouch, atom_untouch, residue_untouch 
	1	,incl, inter, iadd
	logical zone_untouch, atom_untouch, residue_untouch
	integer iadd
	logical incl,inter
	
	common /cmm_group/ 
	1	n_groupa(0:max_gr), igroupa(max_atom,0:max_gr) 
	1	, cgroup, i_pipe
	integer n_groupa, igroupa, i_pipe
	character*4 cgroup(max_gr)
	

	common /cmm_vector/ vectors(7,max_vector) ,cvector(max_vector)
	character*4 cvector
	real vectors

	character*4 version
	data version/'v13a'/

!end module edp_main
