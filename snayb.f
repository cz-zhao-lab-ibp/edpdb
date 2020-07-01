

	logical function do_they_match(i,j,mch)
!  	==============================
	include 'edp_main.inc'
	integer i,j,mch

!	if(text(i)(1:4).eq.text(j)(1:4)) then	!???? what is it for????
	if(i.eq.j) then
	  do_they_match=.false.
	else if(mch.eq.0) then	
	  if(w(i).eq.w(j)) then
	    do_they_match=.true.
	  else
	    do_they_match=.false.
	    endif
	else if(mch.gt.0) then
	  if(atom(i)(:mch).eq.atom(j)(:mch)) then
	    do_they_match=.true.
	  else
	    do_they_match=.false.
	    endif
	else
	  if(atom(i)(:-mch).ne.atom(j)(:-mch)) then
	    do_they_match=.true.
	  else
	    do_they_match=.false.
	    endif
	  endif
	end

	subroutine mcs_atm(*)
!  	==================
	include 'edp_main.inc'
!	use edp_main
	include 'edp_dat.inc'
!	use edp_dat

	parameter 	(max_v=max_l/3)
!	integer, parameter :: max_v=max_l/3
	real xa(3,max_v), vva(max_v,max_v)
	real xb(3,max_v), vvb(max_v,max_v)
	integer	   name_a(2,max_v),  name_b(2,max_v),  label(4,max_l)
	integer mch

	common /clique_vector/ label, rms_min

	real eps(3)

	logical t(max_l,max_l)
	common /cmm_mcs/ is_angle,num_total,t, min_clique, max_cliques
	1 ,jnk(max_l*2-2)

	logical do_they_match

	n_of_syn=2					!000515
	syntax(1)='syntax:' 
	syntax(2)=
	1'clique group_id.s cmin.i rms_cutoff.r eps.r '//
	1'[max_#_cliques.i [mch_char.i [mch_txt.i]]'

	igr= match_l( max_gr, cgroup)
	if( igr .le. 0) then
	  errmsg=' errmsg: Define the group first.'
	  return 1
	  endif

	call read_ar(3,eps,*901,*901)
	min_clique=eps(1)
	rms_min=eps(2)

	max_cliques=10
	call read_ai(1,max_cliques,*901,*40)

40	mch=1
	call read_ai(1,mch,*901,*45)
	
45	no_dup=0
	call read_ai(1,no_dup,*901,*50)

50	if(mch.lt.-4.or.mch.gt.4) return 1
	if(verbose .ge. 2) then
	  write(6,1005) 
	1     cgroup(igr), min_clique, rms_min, eps(3), max_cliques
	  if(mch.gt.0) then
	    write(6,1006) mch
	  else if(mch.lt.0) then
	    write(6,1007) -mch
	  else 
	    write(6,1008) 
	   endif
	  if(no_dup.eq.1) then
	    write(6,1019)
	  else
	    write(6,1018) 
	   endif
	  endif
1005	format(
	1' clique-I2> match against group ',a/
	1' clique-I2> min# of atoms in a clique=',i5/
	1' clique-I2> rms cutoff to the output list =',f8.3/
	1' clique-I2> rms cutoff to a vector pair selection=',f8.3/
	1' clique-I2> max# of cliques to be listed =',i5)

1006	format(
	1' clique-I2> match is based on the first',i2,
	1' character(s) in the atom names')
1007	format(
	1' clique-I2> match is based on different first',i2,
	1' character(s) in the atom names')
1008	format(
	1' clique-I2> match is based on identical weight')
1018	format(
	1' clique-I2> duplications in cliques are allowed')
1019	format(
	1' clique-I2> only output cliques that do not contain duplications')

	num_a=0
	do ir=1,n_res
	do i=ijk(ir),ijk(ir+1)-1
	  if(lf(i)) then
	    num_a=num_a+1
	    if(num_a.gt.max_v) then
	      write(6,*) 'clique-W> too many ON atoms.'
	      if(verbose .ge. 3) write(6,*)
	1 '%EdPDB-I3- increase max_l in edp_dim.inc (currently'
	1 ,max_l,')'
	      return
	      endif 
	    xa(1,num_a)=x(i)
	    xa(2,num_a)=y(i)
	    xa(3,num_a)=z(i)
	    name_a(1,num_a)=i
	    name_a(2,num_a)=ir
	    endif
	  enddo
	  enddo
	num_b=n_groupa(igr)
	if(num_b.gt.max_v) then
	  write(6,*) 
	1 'clique-W> too many atoms in group [',cgroup(igr),'].'
	  if(verbose .ge. 3) write(6,*) 
	1 '%EdPDB-I3- increase max_l in edp_dim.inc (currently'
	1 ,max_l,')'
	  return
	  endif 

	do ii=1,num_b
	  i=igroupa( ii,igr)
	  xb(1,ii)=x(i)
	  xb(2,ii)=y(i)
	  xb(3,ii)=z(i)
	  name_b(1,ii)=i
	  name_b(2,ii)=aa_seq(igroupa(ii,igr))
	  enddo

	write(6,1009)   num_a,cgroup(igr),num_b
1009	format(' clique> #of on atoms =',i5,
	1 ', #of atoms in group ',a,' =',i5)

	if(min_clique.lt.3.or.min_clique.gt.min(num_a,num_b)) then
	  write(6,*) 
	1 'clique-W> UNDONE: inproper min#_clique in the command line.'
	  return 
	  endif

	call make_vv(num_a,xa,vva)
	call make_vv(num_b,xb,vvb)

	n0=0
	do i=1,num_a
	do j=1,num_b
	  if( do_they_match(name_a(1,i),name_b(1,j),mch) ) then
	    n0=n0+1
	    if(n0.gt.max_l) then
	      write(6,*)'clique-W> UNDONE: too many potential matches.'
	      if(verbose .ge. 3 ) write(6,*) 
	1 '%EdPDB-I3- increase max_l in edp_dim.inc (currently'
	1 ,max_l,')'
	      return
	      endif 
	    label(1,n0)=name_a(1,i)
	    label(2,n0)=name_a(2,i)
	    label(3,n0)=name_b(1,j)
	    label(4,n0)=name_b(2,j)
	n1=0
	do i1=1,num_a
	do j1=1,num_b
	  if( do_they_match(name_a(1,i1),name_b(1,j1),mch) ) then
	    n1=n1+1
	    if(n1.gt.max_l) then
	      write(6,*) 
	1 'clique-W> UNDONE: too many potential matches.'
	      if(verbose .ge. 3 ) write(6,*) 
	1 '%EdPDB-I3- increase max_l in edp_dim.inc (currently'
	1 ,max_l,')'
	      return
	      endif 
	    t(n0,n1)=.true.
	    if(abs(vva(i,i1) - vvb(j,j1)).gt.eps(3)) t(n0,n1)=.false.
!	    if(i.eq.i1 .or. j.eq.j1 ) t(n0,n1)=.false.
	    endif
	  enddo
	  enddo
	    t(n0,n0)=.true.
	    endif
	  enddo
	  enddo
	num_total=n0
	call clique(no_dup)
	call typelist(48)
	return
901	errmsg=' errmsg: input min_clique, rms_cutoff, eps & max_#_cliques.'
	return 1
	end

	subroutine make_vv(num_a,xa,vva)
!  	==================
	include 'edp_dim.inc'
!	use edp_dim

	common /cmm_mcs/ is_angle
	parameter 	(max_v=max_l/3)
!	integer, parameter :: max_v=max_l/3
	real xa(3,max_v), vva(max_v,max_v)

	if(is_angle.eq.1) then
	 do i=  1,num_a
	 do j=i+1,num_a
	  dr=polar_diff_func(xa(1,i),xa(1,j))
	  vva(j,i)=dr	  
	  vva(i,j)=dr
	  enddo
	  enddo
	else ! is_angle .eq. 0
	 do i=  1,num_a
	 do j=i+1,num_a
	  dx=xa(1,i)-xa(1,j)
	  dy=xa(2,i)-xa(2,j)
	  dz=xa(3,i)-xa(3,j)
	  dr=sqrt(dx*dx+dy*dy+dz*dz)
	
	  vva(j,i)=dr	  
	  vva(i,j)=dr
	  enddo
	  enddo
	endif
	end

	subroutine after_clique(num_node,node_name,num_cliques,rms, no_dup)
!  	=======================
	include 'edp_main.inc'
	include 'edp_dat.inc'

	parameter 	(max_v=max_l/3)
!	integer, parameter :: max_v=max_l/3
	integer  node_name(num_node), k(max_l)
	real s(3,max_v), p(3,max_v), weight(max_v)
	data weight/max_v*1./

	integer label(4,max_l)

	common /cmm_mcs/ is_angle
	real tmp1(3,3), tmp2(3,3)
	common /clique_vector/ label, rms_min 
	character*20 text_c(max_l)
	real trn(3,4), a_polar(3)

	sum=0.
	kk=0
	do i=1,num_node
	  ni=node_name(i)
	  n1=label(1,ni)
	  n2=label(3,ni)
	  if(is_angle.eq.1) then
	    call polar_matrix(x(n1),y(n1), z(n1),tmp1)
	    call polar_matrix(x(n2),y(n2), z(n2),tmp2)
	    do k1=1,3
	     kk=kk+1
	     do k2=1,3
	      s(k2,kk)=tmp1(k1,k2)
	      p(k2,kk)=tmp2(k1,k2)
	      sum=sum+(s(k2,kk)-p(k2,kk))**2
	     enddo
	    enddo
	  else
	  kk=kk+1
	  s(1,i)= x(n1)
	  s(2,i)= y(n1)
	  s(3,i)= z(n1)
	  p(1,i)= x(n2)
	  p(2,i)= y(n2)
	  p(3,i)= z(n2)
	  sum=sum+ (x(n1)-x(n2))**2+ (y(n1)-y(n2))**2+ (z(n1)-z(n2))**2
	  endif 
	  enddo
	sum=sqrt(sum/kk)
	
	rms=0.0 ! this is necessary for getting rid of det(M) < 0.0
	call edp_rms_doit(rms,kk,s,p,weight,trn,0
	1  ,.false.,.true.)
!	if(is_angle.eq.1) rms=2.0*asind(0.5*rms)
	if(rms.lt.0.) return

	if(rms.lt.rms_min) then
	    call sort_in2(num_node,node_name,k)
	    do ii=1,num_node
	      i=node_name(k(ii)) 
	      n1=label(1,i)
	      n2=label(3,i)
	      m1=label(2,i)
	      m2=label(4,i)
	      text_c(ii)=atom(n1)//ires(m1)//atom(n2)//ires(m2)
	if(no_dup.eq.1) then
	 do jj=ii+1,num_node
	      j=node_name(k(jj)) 
	      nj1=label(1,j)
	      nj2=label(3,j)
	      mj1=label(2,j)
	      mj2=label(4,j)
	      if(atom(n1).eq.atom(nj1).and.ires(m1).eq.ires(mj1)) return
	      if(atom(n2).eq.atom(nj2).and.ires(m2).eq.ires(mj2)) return
	  enddo
	 endif      
	      enddo
	    status= num_node	      
	    write(48,1003) rms,(i,text_c(i),i=1,num_node)
1003	    format(/' rms <=',g12.4/(2x,i4,2X,': ',a))
	    if(verbose .ge. 2 ) then
	      call get_polar(trn,a_polar)
	      write(48, 1001)
	1  a_polar(1), a_polar(2), a_polar(3), 
	1  trn(1,4), trn(2,4), trn(3,4),
	1  sum
c	1  nint(a_polar(1)), nint(a_polar(2)), nint(a_polar(3)),
c	1  nint(trn(1,4)), nint(trn(2,4)), nint(trn(3,4)),      
	      endif
1001	format(/' rtn polar ',6f7.1,' ! curr.rms=',g12.4)
	    num_cliques=num_cliques + 1
	    endif
	end

	subroutine clique(no_dup)
!  	=================
	include 'edp_dim.inc'
!	use edp_dim
	include 'edp_dat.inc'
!	use edp_dat

	integer d, num_c(max_l), name_c(max_l,max_l)
	integer  name_d(max_l)
	logical  t(max_l,max_l)
	common /cmm_mcs/ is_angle,num_total,t, min_clique, max_cliques
	1 ,jnk(max_l*2-2)

	if(verbose .ge. 6) write(6,1069) 	!000504
1069	format(' clique-I6> references:'
	1/'  (1) Grindley HM, Artymiuk PJ, Rice DW, Willett P. (1993)'
	1/'  Identification of tertiary structure resemblance in proteins '
	1/'  using a maximal common subgraph isomorphism algorithm.'
	1/'  J Mol Biol., 229(3):707-21.'
	1/'  (2) Bron, C & Kerbosch J (1973).'
	1/'  Algorithm 457, finding all cliques of an undirected graph.'
	1/'  Commum. A.C.M 16, 575-577.')

	n_true=0
	do i=1,num_total
	do j=1,num_total
	  if(t(j,i)) n_true=n_true+1
c	  name_c(j,i)=0
	  enddo
c	  name_d(i)=0
c	  num_c(i)=0
	  enddo

	if(verbose .ge. 2 ) then
	  write(6,1004) num_total,num_total,n_true
1004	format(/' clique-I2> The matrix is ',i6,' x ',i6,
	1 ' in dimension, with',i8,'T'/)
	  endif

!  	initialize:
	num_cliques=0
	d=1
	do i=1,num_total
	  name_c(i,d)=i
	  enddo
	num_c(d)=num_total

!  	find the node of the maximum # of connections
200	num_ce=num_c(d)
	m1=0
	do i=1,num_ce
	  m0=0
	  do j= 1,num_ce
	    if(t(name_c(j,d),name_c(i,d))) m0=m0+1
	    enddo
	  if(m0.gt.m1) then
	    m1=m0
	    i1=i
	    endif
	  enddo
	if(m1.lt.min_clique-d ) goto 600

!  	goto d+1 level
400	name_c_i=name_c(i1,d)
	name_d(d)=name_c_i
	num_ce=num_c(d)
	m1=0
	do i= 1,num_ce
	  ii=abs(name_c(i,d))
	  if(t(name_c_i,ii)) then
	    if(ii.ne.name_c_i) then
	      m1=m1+1
	      name_c(m1,d+1)=ii
	      endif
	    name_c(i,d)=-ii
	    endif
	  enddo
	d=d+1
	num_c(d)=m1
	if(m1.gt.0) goto 200

!  	Great! found a clique
	if(d-1.ge.min_clique) then
	  call after_clique(d-1,name_d,num_cliques,rms, no_dup)
	  if(num_cliques.ge.max_cliques) then
	    write(6,*) 
	1 'clique-W> there may be more cliques than those listed.'
	    return
	    endif
	  endif

!  	Come back to d-1 level
600	d=d-1
	if(d.le.0) return		! finished
	num_ce=num_c(d)
	do i=1,num_ce
	  if(name_c(i,d).gt.0) then
	    i1=i
	    goto 400
	    endif
	   enddo
	goto 600
	end

copyright by X. Cai Zhang

	subroutine dfcell_b(cell0,irl0)
!  	===================
	include 'edp_main.inc'
!	use edp_main

	real cell0(6), cell, trn0, trn1
	common /cmm_cell/ cell(6), irl, trn0(3,3), trn1(3,3)
	common /cmm_symm/ num_symm
	1 ,jnk1(3,4,max_symm) ,jnk2(3,4,max_symm)

	character*(max_num_chars) txt
	logical nword0
	common /cmm_txt/n_len,txt,ib,ie
	save	!030422

	do i=1,6
	  cell(i)=cell0(i)
	  enddo
	irl=irl0

3100	if( (cell(1)*cell(2)*cell(3)*cell(4)*cell(5)*cell(6).le.0.) .or.
	1 (cell(4)+cell(5)+cell(6) .gt. 360.0) )then
	  write(6,*) '%EdPDB-E- undefined or wrong cell parameters'
	  return
	  endif
	call trnslnb0(cell,trn0,irl)
	call trnslnb1(cell,trn1,irl)
	num_symm=0	
	return

	entry dfcell(*)
c	===============
	n_of_syn=12					!000515
	syntax(1)='syntax:'
	syntax(2)='1) cell' 
	syntax(3)='2) cell a.r b.r c.r alpha.r beta.r gamma.r convention.i'
	syntax(4)='   where convention is 1-8:'
	syntax(5)='     1: x//a*, y//b, z//(a* X b) '
	syntax(6)='     2: x//(b X c*), y//b, z//c* '
	syntax(7)='     3: x//(b* X c), y//b*, z//c '
	syntax(8)='     4: x//a*, y//(c X a*), z//c '
	syntax(9)='     5: x//a, y//(c* X a), z//c* '
	syntax(10)='     6: x//a, y//b*, z//(a X b*) '
	syntax(11)=
	1'     7: x//(a-b), y//(a+b-2c), z//(a+b+c) -- only for R+ lattice'
	syntax(12)=
	1'     8: x//(a-c), y//(2b-a-c), z//(a+b+c) -- only for R- lattice'
 
	ie0=ie
	if(nword0(n_len,txt,ib,ie)) then
	  call shcell
	  return
	  endif

	ie=ie0
	call read_ar(6, cell, *9001, *3900)
3900	if(cell(1).le.0..or.
	1  cell(2).le.0..or.
	1  cell(3).le.0..or.
	1  cell(4).le.0..or.
	1  cell(5).le.0..or.
	1  cell(6).le.0.) return 1
	if(cell(4).ne.90..or.cell(5).ne.90..or.cell(6).ne.90.) then
	  call read_ai(1,irl,*9001,*3901)
3901	  if(irl.lt.1.or.irl.gt.8) goto 9001
	else
	  irl=1
	  endif
	if( (cell(1)*cell(2)*cell(3)*cell(4)*cell(5)*cell(6).le.0.) .or.
	1 (cell(4)+cell(5)+cell(6) .gt. 360.0) )then
	  write(6,*) '%EdPDB-E- undefined or wrong cell parameters'
	  return
	  endif
	call trnslnb0(cell,trn0,irl)
	call trnslnb1(cell,trn1,irl)
	num_symm=0	
	return
9001	return 1
	end

	subroutine shcell
!  	=================
	include 'edp_dat.inc'
!	use edp_dat
	include 'edp_dim.inc'
!	use edp_dim

	character*25 rule(8)
	real cell
	common /cmm_cell/ cell(6), irl
	1 ,jnk3(3,3) ,jnk4(3,3)

	common /cmm_symm/ num_symm
	1 ,jnk1(3,4,max_symm) ,jnk2(3,4,max_symm)

	data rule/
	1 ' x//a*, y//b, z//(a* x b)'
	1,' x//(b x c*), y//b, z//c*'
	1,' x//(b* x c), y//b*, z//c'
	1,' x//a*, y//(c x a*), z//c'
	1,' x//a, y//(c* x a), z//c*'
	1,' x//a, y//b*, z//(a x b*)'
	1,'x//a-b,y//a+b-2c,z//a+b+c'
	1,'x//a-c,y//2b-a-c,z//a+b+c'/

	logical no_cell

901	if(num_symm.lt.0) then
	  if(no_cell()) then
	    write(6,*) 'cell-W> undefined'
	    return
	  else
	    goto 901
	    endif
	else
	  write(6,1012) cell, irl, rule(irl)
1012	  format(' cell ',6f7.2,i3, ' ! ',a25)
	  if(verbose .ge. 2 ) call check_cell(cell)
	  endif
	end

	subroutine get_dist(xc,yc,zc,xp,yp,zp,d,r,*)
c	===================
	x=abs(xp-xc)
	if(x.gt.d) return 1
	y=abs(yp-yc)
	if(y.gt.d) return 1
	z=abs(zp-zc)
	if(z.gt.d) return 1

	r= x*x+ y*y+ z*z
	if(r.gt.d*d) return 1
	r=sqrt(r)
	end

	subroutine snayb(*)
!  	==================
c	snayb selects the atoms in a crystal sphere centered at a given atom of a given
c	residue with a specified radias.
c 
	include 'edp_main.inc'
!	use edp_main
	include 'edp_dat.inc'
!	use edp_dat
	include 'edp_file.inc'
!	use edp_file

!	character*(*) note
	character*(108) file_name

	character*(max_num_chars) txt
	logical nword 
	common /cmm_txt/n_len,txt,ib,ie

	logical ierr, ls, ln,unitary
	common /cmm_cell/ cell(6), irl, trn0(3,3), trn1(3,3)
	real sym(3,4), 
	1 arr(3,3), vec(3),  ps(3)
	1 ,ox(3), px(3), qx(3), rx(3)
	character*32 symm_txt(max_symm)
	
	dimension trna(3,4),trnb(3,3)
	character*32 symm_txta
	common /cmm_symm/ num_symm, symmtx(3,4,max_symm), 
	1 trn2(3,4,max_symm)
	
	character*(10) string_init, string_punch
	data string_init/'initialize'/, string_punch/'punch'/

	data ox/0.5,0.5,0.5/
	logical no_cell

	logical no_symm
	external no_symm

	logical  open_file
	external open_file
!	save symm_txt
	save	!030422

	n_of_syn=3					!000515
	syntax(1)='syntax:' 
	syntax(2)='1) snayb radius.r [res_id.s [atom_name.s]]' 
	syntax(3)='2) snayb radius.r center x.r y.r z.r '

	if(num_symm.le.0) then
	  if(no_symm()) then
	    write(6,*) 
	1 '%EdPDB-W- UNDONE: symmetry is undefined'
	    return
	    endif
	  end if

!  	igr:    id of the buffer
        call dfgroup(igr,*9001)

!  	input search_radius
	call read_ar(1, rad, *9001, *9001)
		
!  	input atom_id
	call get_one_atom( ia, cx, cy, cz)
        if(ia.lt.0 ) return 1
	
1002	format(' snayb-I> the occ will be changed to the dist.',
	1' symm. op. on the given atom.')

100	ln=incl
	r=rad*rad
	
	j1=n_groupa(igr)
	if(j1.le.0) return
	xc0=0.
	yc0=0.
	zc0=0.
	rmax0=0.
	do j=1,j1
          i=igroupa( j,igr)
	  xc0=xc0+x(i)
	  yc0=yc0+y(i)
	  zc0=zc0+z(i)
	  rmax0=max(rmax0,x(i)*x(i)+y(i)*y(i)+z(i)*z(i))
	  enddo
	xc0=xc0/j1
	yc0=yc0/j1
	zc0=zc0/j1
	rmax0=sqrt(rmax0)

	do is=1,num_symm
	  do i=1,3
	    vec(i)= trn2(i,1,is)*cx+
	1           trn2(i,2,is)*cy+
	1           trn2(i,3,is)*cz+
	1           trn2(i,4,is)
	    enddo 

	  do 800 ja=-2,2
	  do 800 jb=-2,2
	  do 800 jc=-2,2
	    ls=.false.
	    do i=1,3
	      ps(i)=vec(i)+ trn1(i,1)*ja+trn1(i,2)*jb+trn1(i,3)*jc
	      end do
c
c	rps is the distance between the geometry centre and the new position 
c 	 of the atom
c
	    rps=sqrt((ps(1)-xc0)**2+(ps(2)-yc0)**2+(ps(3)-zc0)**2)
	    if(rps.gt.rmax0+rad) goto 800 

	    if( abs(ps(1)-cx).lt.1.e-3.and.
	1       abs(ps(2)-cy).lt.1.e-3.and.
	1       abs(ps(3)-cz).lt.1.e-3) goto 800

	    do 200 j=1,j1
              i=igroupa( j,igr)
	      dx=abs(ps(1)-x(i))
	      dy=abs(ps(2)-y(i))
	      dz=abs(ps(3)-z(i))
	      if(dx.gt.rad) goto 200
	      if(dy.gt.rad) goto 200
	      if(dz.gt.rad) goto 200
	      d=dx*dx+dy*dy+dz*dz
	      if(d.gt.r) goto 200
	      lf(i)=incl
	      if(incl) w(i)=sqrt(d)
	      ls=.true.
200	      continue
	    if(ls) then		! new symmetry
	      if(ln) write(6,1002) 	! MMI information 
	      ln=.false.
	      write(6,1030) is, symm_txt(is), ja,jb,jc
	      if(verbose.ge.3) then
 	        write(6,1230) is,ja,jb,jc
1230	format(' snayb-I3> To apply the transformation, type:'
	1 /'   rtn symmetry',4i4)
	      endif
	      endif
800	    continue
	  enddo
	return
1030	format(' snayb> symmetry #',i2,': ',a32,
	1 ' plus [',i2,',',i2,',',i2,']')
9001	return 1

	entry snaybr(*)
c	==================
	n_of_syn=2					!000515
	syntax(1)='syntax:'
	syntax(2)='snaybr distance.r res_id.s '

	if(num_symm.le.0) then
	  if(no_symm()) then
	    write(6,*)
	1 '%EdPDB-W- UNDONE: symmetry is undefined'
	    return
	    endif
	  endif

!  	input search_radius
	call read_ar(1, rad, *9001, *9001)
c-		
	call find(n_res,ires,jq)
	if(jq.le.0) return 1

2100	r=rad*rad
	ijkjq=ijk(jq)
	ijkjq1=ijk(jq+1)-1

	j1=n_groupa(igr)
	if(j1.le.0) return
	xc0=0.
	yc0=0.
	zc0=0.
	rmax0=0.
	do j=1,j1
          i=igroupa( j,igr)
	  xc0=xc0+x(i)
	  yc0=yc0+y(i)
	  zc0=zc0+z(i)
	  rmax0=max(rmax0,x(i)*x(i)+y(i)*y(i)+z(i)*z(i))
	  enddo
	xc0=xc0/j1
	yc0=yc0/j1
	zc0=zc0/j1
	rmax0=sqrt(rmax0)

	do 2900 is=1,num_symm
	if(abs(trn2(1,1,is)-1.).lt.1.e-4.and.
	1  abs(trn2(2,2,is)-1.).lt.1.e-4.and.
	1  abs(trn2(3,3,is)-1.).lt.1.e-4) then
	 unitary=.true.
	else
	 unitary=.false.
	end if
 
	do 2800 ja=-2,2
	do 2800 jb=-2,2
	do 2800 jc=-2,2
	do i=1,3
	 vec(i)= trn1(i,1)*ja+trn1(i,2)*jb+trn1(i,3)*jc+ trn2(i,4,is)
	end do

	if(abs(vec(1)).lt.1.e-4.and.
	1  abs(vec(2)).lt.1.e-4.and.
	1  abs(vec(3)).lt.1.e-4.and.unitary) goto 2800

	ls=.false.
	do 2700 iq=ijkjq, ijkjq1
	 cx=x(iq)
	 cy=y(iq)
	 cz=z(iq)
	 do i=1,3
	  ps(i)=vec(i)+
	1 trn2(i,1,is)*cx+ trn2(i,2,is)*cy+ trn2(i,3,is)*cz
	 end do
c
c	rps is the distance between the geometry centre and the new position 
c 	 of the atom
c
	 rps=sqrt((ps(1)-xc0)**2+(ps(2)-yc0)**2+(ps(3)-zc0)**2)
	 if(rps.gt.rmax0+rad) goto 2700
	do 2200 i=1,n_atom
c	if(lf(i).ne.incl) then
	dx=ps(1)-x(i)
	dy=ps(2)-y(i)
	dz=ps(3)-z(i)
	if(dx.gt.rad) goto 2200
	if(dy.gt.rad) goto 2200
	if(dz.gt.rad) goto 2200
	d=dx*dx+dy*dy+dz*dz
	if(d.gt.r) goto 2200
	lf(i)=incl
	ls=.true.
c	end if
2200	continue
2700	continue
	if(ls)	then
	  write(6,1130) symm_txt(is),ja,jb,jc
	  if(verbose.ge.3) write(6,1230) is,ja,jb,jc
	endif
1130	format(' snaybr> symmetry #',i2,': ',a32,
	1' plus [',i2,',',i2,',',i2,'], on the res.')
2800	continue
2900	continue
	return


	entry dfsymm(*)
c	entry dfsymm(*,note)
c	====================
c	n_of_syn=4					!000515
c	syntax(1)='syntax:'
c	syntax(2)='symmetry '
c	syntax(3)='symmetry operator (e.g. -x,-y,z)'
c	syntax(4)='symmetry initialize'

	  if(num_symm.lt.0) then
	    if(no_cell()) then
	      write(6,*) 
	1 '%EdPDB-W- UNDONE: cell parameters must be defined first'
	      return
	    endif
	  endif

	if(nword(n_len,txt,ib,ie)) then 
	  goto 902
	endif
	
	i=ltrim(txt(ib:ie))
	if(i.le.0) goto 902

	if (txt(ib:ie).eq.string_punch(1:min(i,len(string_punch)))) then
	  if(nword(n_len,txt,ib,ie)) return 1
	  if(.not.open_file(22,txt(ib:ie),'unknown','.edp')) return 1
	  rtn_out=txt(ib:ie) 
	  if(nword(n_len,txt,ib,ie)) then 
	    call chksymm(22,num_symm,symmtx,trn2)
	  else if(txt(ib:ie).eq.'o') then
	    call punch_sym_o(22,num_symm,symmtx)
	  else
	    return 1
	    endif	    
	  close(22)
	  return
	else if(txt(ib:ie).eq.string_init(1:min(i,len(string_init))))then 
	  if(verbose .ge. 3) then
	    write(6,*) 'symmetry-W3> ',num_symm,' operators are erased.'
	  endif
 	  num_symm=0
	  return
	  endif
	  
	if(num_symm.eq.max_symm) then
	  write(6,*) 
	1 'symmetry-W> UNDONE: too many symmetry operators.'
	  if(verbose .ge. 3 ) write(6,*)
	1 '%EdPDB-I3- increase max_symm in edp_dim.inc (currently'
	1 ,max_symm,')'
 	  return
	  endif

	num_symm=num_symm+1
	symm_txt(num_symm)=txt(ib:n_len)

	call symm(symm_txt(num_symm),32,sym,ierr)
	if(ierr) then
	  num_symm=num_symm-1
	  return 1
	  endif

	trn2(1,4,num_symm)= 
	1 trn1(1,1)*sym(1,4)+trn1(1,2)*sym(2,4)+trn1(1,3)*sym(3,4)
	trn2(2,4,num_symm)=
	1 trn1(2,1)*sym(1,4)+trn1(2,2)*sym(2,4)+trn1(2,3)*sym(3,4)
	trn2(3,4,num_symm)= 
	1 trn1(3,1)*sym(1,4)+trn1(3,2)*sym(2,4)+trn1(3,3)*sym(3,4)
	call axbeqc(trn1,sym,arr)
	call axbeqc(arr,trn0,trn2(1,1,num_symm))

	do j=1,4
	do i=1,3
	  symmtx(i,j,num_symm)=sym(i,j)
	  enddo
	  enddo
	return

	entry shsymm
!  	============
902	if(num_symm.le.0) then
	  if(no_symm()) then 
	    write(6,*) 'symmetry-W> undefined'
	    return
	    endif
	  endif
	
	call chksymm(0,num_symm,symmtx,trn2)
	return

	entry get_trn(isymm,trna,trnb,symm_txta)
c	=============
	if(isymm.ge.0) then
	  if(num_symm.le.0) then
	    if(no_symm()) then
	      write(6,*) 
	1 '%EdPDB-W- UNDONE: symmetry is undefined'
	      return
	      endif
	  endif
	else 
	  if(num_symm.lt.0) then
	    if(no_cell()) then
	      write(6,*) 
	1 '%EdPDB-W- get_cell UNDONE: cell parameters are undefined'
	      return
	      endif
	    endif
	  end if

	if(isymm.le.-3) then  ! see writef in readf.for
	 trna(1,1)=cell(1)
	 trna(2,1)=cell(2)
	 trna(3,1)=cell(3)
	 trna(1,2)=cell(4)
	 trna(2,2)=cell(5)
	 trna(3,2)=cell(6)
	endif

	isymm=isymm+1
	if(isymm.lt.0) then
	  do i=1,3
	  do j=1,3
	    trnb(j,i)=trn0(j,i)
	  end do
	  end do
	else
	  do i=1,3
	  do j=1,3
	    trnb(j,i)=trn1(j,i)
	  end do
	  end do
	end if
	
	if(isymm.le.0) return

	if(isymm.gt.num_symm) then
	  if(isymm.le.0)  write(6,*) 
	1 '%EdPDB-W- UNDONE: symmetry is undefined'
	  isymm=0
	  return
	end if

	symm_txta=symm_txt(isymm)

	do i=1,4
	do j=1,3
	  trna(j,i)=trn2(j,i,isymm)
	end do
	end do
	return

	entry move2o(*)
c	============
	n_of_syn=2					!000515
	syntax(1)='syntax:'
	syntax(2)=
	1'movecenter [filename.s] [fx1.r fy1.r fz1.r [fx2.r fy2.r fz2.r]]' 

	if(num_symm.le.0) then
	  if(no_symm()) then
	    write(6,*)
	1 '%EdPDB-W- UNDONE: symmetry is undefined'
	    return
	    endif
	  endif

	if(nword(n_len,txt,ib,ie) ) then
	  file_name='rtn_.txt'
	else
	  file_name=txt(ib:ie)
	  endif
	if(.not.open_file(29,file_name,'unknown','.txt')) then
	  errmsg=' errmsg: failure in opening the file '
	  return 1
	  endif
	rtn_out=file_name

	px(1)=0.
	px(2)=0.
	px(3)=0.
	call read_ar(3, ox(1), *5008,*5001)
5001	call read_ar(3, px(1), *5008,*5004)

c	define centriod of the on atoms
5004	vec(1)=0.
	vec(2)=0.
	vec(3)=0.
	na=0
	do i=1,n_atom
	  if(lf(i)) then
	    na=na+1
	    vec(1)=vec(1)+ x(i)
	    vec(2)=vec(2)+ y(i)
	    vec(3)=vec(3)+ z(i)
	    endif
	  enddo
	if(na.le.0) then
	  write(6,*) 'move-W>  UNDONE: there is no on atom.'
	  return
	  endif
	vec(1)=vec(1)/na	
	vec(2)=vec(2)/na	
	vec(3)=vec(3)/na

c	px -- the 1st center (ox)
c	ps -- the 2nd center
	do i=1,3
	 ps(i)= trn1(i,1)*px(1)+trn1(i,2)*px(2)+trn1(i,3)*px(3)
	enddo
	do i=1,3
	 px(i)= trn1(i,1)*ox(1)+trn1(i,2)*ox(2)+trn1(i,3)*ox(3)
	enddo

	d_min1=1.e8
	d_min2=1.e8
	do is=1,num_symm
	  do i=1,3
	    ox(i)= trn2(i,1,is)*vec(1)
	1         +trn2(i,2,is)*vec(2)
	1         +trn2(i,3,is)*vec(3)
	1         +trn2(i,4,is)
	    enddo 
	  do ja=-2,2
	  do jb=-2,2
	  do jc=-2,2
	    do i=1,3
	      qx(i)= ox(i)+ trn1(i,1)*ja+trn1(i,2)*jb+trn1(i,3)*jc
	    enddo
	    dx=sqrt( (px(1)-qx(1))**2+(px(2)-qx(2))**2+(px(3)-qx(3))**2)
	    ds=sqrt( (ps(1)-qx(1))**2+(ps(2)-qx(2))**2+(ps(3)-qx(3))**2)

	    if(abs(d_min1-dx).lt.1.e-3) then
	      if(d_min2-ds.gt.1.e-3) then
	        d_min1=dx
	        d_min2=ds
	        min_is=is
	        min_ja=ja
	        min_jb=jb
	        min_jc=jc
	        rx(1)=qx(1)
	        rx(2)=qx(2)
	        rx(3)=qx(3)
	      endif 
	    else if(dx.lt.d_min1) then
	      d_min1=dx
	      d_min2=ds
	      min_is=is
	      min_ja=ja
	      min_jb=jb
	      min_jc=jc
	      rx(1)=qx(1)
	      rx(2)=qx(2)
	      rx(3)=qx(3)
	    endif
	  enddo
	  enddo
	  enddo	
	enddo

	write(6,1050) min_is,symm_txt(min_is), min_ja,min_jb,min_jc
	if(verbose.ge.3)  write(6,1240) min_is, min_ja,min_jb,min_jc
	1 ,rtn_out(:ltrim(rtn_out))
1240	format(' move-I3> To apply the transformation, type:'
	1 /'   rtn symmetry',4i4,'       ! or '
	1 /'   rtn file ',a)
	
1050	format(' move> symmetry #',i2,': ',a32,
	1 ' plus [',i2,',',i2,',',i2,']')

	write(29,1051) ((trn2(i,j,min_is),j=1,3),i=1,3)
	1, (trn1(i,1)*min_ja +trn1(i,2)*min_jb +trn1(i,3)*min_jc
	1  +trn2(i,4,min_is) , i=1,3)
1051	format(3e16.8)
 
	write(29,1050)min_is,symm_txt(min_is), min_ja,min_jb,min_jc
	close (29)

c	save the center (fractional) for the next run.
	do i=1,3
	  ox(i)= trn0(i,1)*rx(1)+trn0(i,2)*rx(2)+trn0(i,3)*rx(3)
	  enddo
	return
	
5008	return 1

	end 

c	this subroutine read cell dimensions from the input pdb file, 
c	  if possible.
	logical function no_cell()
!  	=========++++++++=======
	include 'edp_dim.inc'
!	use edp_dim

	character*(max_num_chars) txt_save, txt
	common /cmm_txt/n_len,txt,ib,ie

	logical first
	data first/.true./

	real trn(3,4), trn0(3,4), cell(6)
!	save txt_save
!	save n_len_save, ib_save, ie_save

	n_len_save=n_len
	txt_save=txt
	ib_save=ib
	ie_save=ie

	no_cell=.true.
	if(.not.first) return
	first=.false.

	irl=0
	rewind (1)
100	read(1,'(a)',end=900) txt
	n_len=ltrim(txt)
	ie=6
	
	if(txt(1:6).eq.'CRYST1') then 
	  call read_ar(6, cell, *900, *900)
	  if(cell(1).le.0..or.
	1    cell(2).le.0..or.
	1    cell(3).le.0..or.
	1    cell(4).le.0..or.
	1    cell(5).le.0..or.
	1    cell(6).le.0.) goto 900
	  irl=6
	  if(cell(4).eq.90..and.cell(5).eq.90..and.cell(6).eq.90.) goto 300
	else if(txt(1:6).eq.'SCALE1') then
	  do i=1,4
	    call read_ar(1, trn(1,i), *900, *900)
	    enddo
	else if(txt(1:6).eq.'SCALE2') then
	  do i=1,4
	    call read_ar(1, trn(2,i), *900, *900)
	    enddo
	else if(txt(1:6).eq.'SCALE3') then
	  do i=1,4
	    call read_ar(1, trn(3,i), *900, *900)
	    enddo
	  goto 250
	else if(txt(1:4).eq.'ATOM') then
	  if(irl.ne.0) then 
	    write(6,*) 
	1'%EdPDB-W- the alignment convention is guessed to be #',irl
	    goto 300
	    endif
	  goto 900
	  endif
	goto 100

250	if(irl.le.0) goto 900
	irl_min=0
	sum_min=1.0E+06
	do i=1,8
	  call trnslnb0(cell,trn0,i)
	  sum=0.
	  cell_max=0.
	  do j=1,3
	  cell_max=max(cell(j),cell_max)
!	  do k=1,4
	  do k=1,3
	    sum=sum+abs(trn0(j,k)-trn(j,k))
	    enddo
	    enddo
	  if(sum.lt.sum_min) then
	    irl_min=i
	    sum_min=sum
	    endif
	  if(sum*cell_max.le.5.e-3) then
	    irl=i
	    goto 300
	    endif
	  enddo
	write(6,'(a)') 
	1 ' %EdPDB-W- double check the deorthoganization matrix!'
	write(6,'(a)') 
	1 '          the guessed alignment could be wrong.'
	irl=irl_min
c	goto 900

300	call dfcell_b(cell, irl)
	write(6,'(a)') 
	1' cell-I> cell parameters are read from the input pdb file.'
	no_cell=.false.

900	n_len=n_len_save
	txt=txt_save
	ib=ib_save
	ie=ie_save
	end

	logical function no_symm()
!  	================
	include 'edp_main.inc'
!	use edp_main	
	include 'edp_file.inc'
!	use edp_file	
	include 'edp_dat.inc'
!	use edp_dat

	integer np, npp
	character*8  pp
	character*72 pn
	common /ud_param/ npp(max_param), pp(max_param),
	1 np(max_param), pn(max_param)	

	character*(max_num_chars) txt_save
	integer n_len_save,ib_save,ie_save
!	save txt_save
!	save n_len_save, ib_save, ie_save

	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie
	logical nword

	logical  open_file
	external open_file
	logical  no_cell

	common /cmm_symm/ num_symm
	1 ,jnk1(3,4,max_symm) ,jnk2(3,4,max_symm)


	no_symm=.true.

	if(num_symm.lt.0) then
	  if(no_cell()) then
	    write(6,*) 
	1 '%EdPDB-W- UNDONE: cell parameters must be defined first'
	    return
	    endif
	  endif

	txt_save=txt
	n_len_save=n_len
	ib_save=ib
	ie_save=ie

	rewind (1)
	do while (.true.)
	 call read_txt(1,*200)
	 if(txt(1:6).eq.'cryst1') then 
	  j=1
	  txt(j:j)='.'
	  do i=56,68
	    if(txt(i:i).ne.' ') then
	      j=j+1
	      txt(j:j)=txt(i:i)
	      endif
	    enddo
	  if(j.gt.2) then
	    n_len=j
	    ib=0
	    ie=0

	    write(6,*) 
	1'symmetry-I> symmetry inf. is read from the input pdb file.'
	    write(6,*) 
	1'symmetry-I> symmetry inf. is read from '//txt(2:j)//'.edp'
	    call set_parameter(2,'sg',2,txt(2:j),j-1)
	    jou=19
	    edp_in= txt(2:j)
	    if(.not.open_file(jou,edp_in, 'old','.edp')) then
	      edp_in=edp_data(:ltrim(edp_data))//edp_in
	      if(.not.open_file(jou,edp_in, 'old','.edp')) then
	edp_in=edp_data(:ltrim(edp_data))//'symmetry/'//
	1 edp_in(ltrim(edp_data)+1:) 
!for vms, it may take something like the folowing
!	edp_in=trim(edp_data)//'.symmetry]'//edp_in(ltrim(edp_data)+1:) 
	       if(.not.open_file(jou,edp_in, 'old','.edp')) then
	        write(6,*) 
	1 'symmetry-W> failure in opening the file ['//
	1 edp_in(:ltrim(edp_in))//']'
	        return 
	        endif
	       endif
	      endif

150	    do while(.true.)
	      call read_txt(19,*200)
	      if(nword(n_len,txt,ib,ie)) goto 150
	      if(txt(ib:ib+3).ne.'symm') goto 900
	      call dfsymm(*900)
	      enddo
	    endif
	 else if(txt(1:4).eq.'ATOM') then 
	  goto 200
	  endif
	 enddo
900	write(6,*) 'symmetry-W> wrong parameters'
200	close (jou)
	n_len=n_len_save
	ib=ib_save
	ie=ie_save
	txt=txt_save
	if(num_symm.ge.1) no_symm=.false.
	return
	end
c****	end of snayb.for

copyright by X. Cai Zhang
