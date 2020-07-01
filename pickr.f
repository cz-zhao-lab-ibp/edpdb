	subroutine clique_sp(name_a,lnk)
!  	================================
	include 'edp_main.inc'
	include 'edp_dat.inc'

	parameter (io=48)
	parameter (max_v=50)
	integer	   name_a(max_v), jorder(max_v)
	integer    lnk(max_v,max_v)
	character*5 tmp(max_v)

	integer d, num_c(max_l), name_c(max_l,max_l)
	integer  name_d(max_l)
	logical  t(max_l,max_l)
	common /cmm_mcs/ is_angle,num_total,t, min_clique, max_cliques
	1 ,jnk(max_l*2-2)
	

	if(verbose .ge. 6) write(6,1069) 	!000504
1069	format(' sclique-I6> references: Bron, C & Kerbosch J (1973).'
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
1004	format(/' sclique-I2> The matrix is ',i6,' x ',i6,
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
!!	  call after_clique(d-1,name_d,num_cliques,rms)
	  num_node=d-1
	  num_cliques=num_cliques+1
	  call sort_in2(num_node, name_d, jorder)
	  write(io,1012) num_cliques
	1, (ires(name_a(name_d(jorder(i)))),i=1,num_node)
	  do j=1,num_node
	    do i=1,num_node
	      ii=lnk(name_d(jorder(j)), name_d(jorder(i)))
	      if(ii.le.0) then
	        tmp(i)='?'
	      else if(ii.eq.999) then
	        tmp(i)='U'
	      else
	        tmp(i)=ires(ii)
	      endif
	    enddo
	    write(io,1010) ires(name_a(name_d(jorder(j))))
	1,(tmp(i), i=1,num_node)
	  enddo
1010	format(1x,20a5)
1012	format(/' clique#',i4,': O(i)= O(j) x O(table)'
	1/1x,'i\\j=',t7,19a5)
	
	  if(num_cliques.ge.max_cliques) then
	    write(6,*) 
	1 'sclique-W> there may be more cliques than those listed.'
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
	
	integer function pp_pair(i,j,igr,epsl)
!	========================
	include 'edp_main.inc'
	include 'edp_dat.inc'
	real pv(3), pu(3), tmp1(3,3),tmp2(3,3),tmp3(3,3)
	
	pp_pair=0
	call polar_matrix(x(i),y(i),z(i),tmp1)
	call polar_matrix(x(j),y(j),-z(j),tmp2)
	call axbeqc(tmp2,tmp1,tmp3)
	call get_polar(tmp3,pv)

	if(igr.gt.0) then
	 n_group2= n_groupa(igr)
	 do kk=1,n_group2
	  k= igroupa(kk,igr)
	  pu(1)=x(k)
	  pu(2)=y(k)
	  pu(3)=z(k)
	  curr=polar_diff_func(pv,pu)
	  if(curr.lt.epsl) then
	    pp_pair=k
	    return
	  endif
	 enddo
	else
	  do k=1,n_atom
	    if(lf(k)) then
 	      pu(1)=x(k)
	      pu(2)=y(k)
	      pu(3)=z(k)
	      curr=polar_diff_func(pv,pu)
	      if(curr.lt.epsl) then
	        pp_pair=k
	        return
	      endif
	    endif
	  enddo
	endif
	if(i.eq.j) pp_pair=999
	end

	subroutine mcs_atm_sp(*)
!  	==================
	include 'edp_main.inc'
	include 'edp_dat.inc'
	include 'edp_file.inc'

	parameter (io=48)
	parameter (max_v=50)
	integer	   name_a(max_v)
	integer    lnk(max_v,max_v)

	real epsl(2)

	logical t(max_l,max_l)
	common /cmm_mcs/ is_angle,num_total,t, min_clique, max_cliques
	1 ,jnk(max_l*2-2)
	
	integer  pp_pair
	external pp_pair

	n_of_syn=2					!000515
	syntax(1)='syntax:' 
	syntax(2)='sclique group_id.s cmin.i  eps.r [max_#_cliques.i]'

	status=-2
	epsl(1)=4.0
	epsl(2)=3.0
	max_cliques=10
	igr= match_l( max_gr, cgroup)
	if( igr .lt. 0) then
	  errmsg=' errmsg: Define the group first.'
	  return 1
	else if(igr.eq.0) then
	  goto 50
	else
	  n_group2= n_groupa(igr)
	  if( n_group2.le.0) then
	    write(6,*) 
	1 'sclique-W> UNDONE: define group ['//cgroup(igr)//'] first.'
	    return
	    end if
	  endif

	call read_ar(2,epsl,*901,*50)

	call read_ai(1,max_cliques,*901,*50)

50	min_clique=epsl(1)
	if(verbose .ge. 2) then
	  write(6,1005)   cgroup(igr)
	1, min_clique, epsl(2), max_cliques
	  endif
1005	format(
	1' sclique-I2> match against group [',a,']'/
	1' sclique-I2> min# of atoms in a clique=',i5/
	1' sclique-I2> cutoff for the angle pair selection=',f8.3/
	1' sclique-I2> max# of cliques to be listed =',i5)

	status=-3
	num_a=0
	do i=1,n_atom
	  if(lf(i)) then
	    num_a=num_a+1
	    name_a(num_a)=i
	    if(num_a.gt. max_v) then
	      write(6,1010) max_v
1010	format(' sclique-W> UNDONE: too many on atoms.'
	1,' The max# allowed is ',i3)
	      return	
	    endif
	  endif
	enddo
	if(num_a.le.min_clique) then
	  write(6,1012) min_clique
1012	format(' sclique-W> UNDONE: too few ON atoms.'
	1,' The number should be not less than',i3)
	return
	endif	    
	
	status=0    
	ii=0
	do i=1,n_atom
	  if(lf(i)) then
	    ii=ii+1
	    jj=0
	    do j=1,n_atom
	      if(lf(j)) then
	        jj=jj+1
		if(ii.ge.jj) then
		  lnk(ii,jj)=pp_pair(i,j,igr,epsl(2))
	if(i.ne.j)lnk(jj,ii)=pp_pair(j,i,igr,epsl(2))
		  t(ii,jj)=lnk(ii,jj).gt.0.or.lnk(jj,ii).gt.0
		  t(jj,ii)=t(ii,jj)	! matrix "t" must be symmatrical 
		endif
	      endif
	    enddo ! j=1,n_atom
	  endif
	enddo ! i=1,n_atom

	if(verbose.ge.12) then
	 write(6,1022) num_a	 
	 do i=1,num_a
	  write(6,1020) (lnk(j,i),j=1,num_a)
	 enddo
	endif
1020	format(40i3)
1022	format(' sclique-I12> number of nodes=',i3,'. Connection Map:')

	num_total=num_a
	call clique_sp(name_a,lnk)
	call typelist(io)
	return
901	errmsg=' errmsg: input min_clique, eps & max_#_cliques.'
	return 1
	end

	subroutine sdist_p(*)
!	===============
! status= -3, no symmetry
!         -2, input error
!       >= 0, number of the last out put

! this subroutine is for analysing self rotation function solutions. 

	include 'edp_main.inc'
	include 'edp_dat.inc'
	include 'edp_file.inc'

	parameter (io=48, max_as= 2048)

	character*(max_num_chars) txt
	common /cmm_txt/ n_len,txt,ib,ie
	logical nword0

	character*6 cnum

	real xa(3), xb(3), trna(3,4), tmp1(3,3), tmp2(3,3)
	common /cmm_symm/ num_symm

	character*32 symm_txt

	n_of_syn=2					!000515
	syntax(1)='syntax:' 
	syntax(2)='psdist group_id.s dmin.r dmax.r '

	status=-3
	if(num_symm.le.0) then
	  write(6,1003) 
1003	format(' %EdPDB-W- UNDONE: symmetry is undefined')
	  return
	endif
	
	status=-2
	igr= match_l( max_gr, cgroup)
	if( igr .le. 0) then
	  errmsg=' errmsg: a group name needed'
	  return 1
	  end if
	n_group2= n_groupa(igr)
	if( n_group2.le.0) then
	  write(6,*) 
	1 'pdist-W> UNDONE: define group ['//cgroup(igr)//'] first.'
	  return
	  end if

	n_group1=0
	do i= 1, n_atom
	  if( lf(i)) then
	    n_group1=n_group1+1
	    endif
	  enddo
	if( n_group1 .le. 0) then
	  write(6,'(a)') ' dist-W>  UNDONE: no on atom exist.'
	  return
	  end if

chk	input dmin, dmax
	dmax=0.0
	call read_ar(1, dmin, *990, *49)
	call read_ar(1, dmax, *990, *49)
49	if(dmax.lt.dmin) then
	  rsmax=dmax
	  dmax=dmin
	  dmin=rsmax
	  endif

chk	input iskip
	iskip=1	! 990422
	imax=max_atom
!	storea= .false.
!	igo=0

	if(verbose.ge.4 ) then
	  write(6,'(a,f6.2)') ' psdist-I4> dmin=', dmin
	  write(6,'(a,f6.2)') ' psdist-I4> dmax=', dmax
	endif

	write(6,1009) 'pdist', n_group1, cgroup(igr), n_group2

	ksum=0
!	av=0.
!	sgm=0.
	isymm=0
	do ii=1,num_symm
	  call get_trn(isymm,trna,tmp1,symm_txt)
	  write(io,1030) isymm,symm_txt
1030	format(/' psdist> symmetry #',i2,': ',a)
	  	
	do ir= 1, n_res
	  do 110 i= ijk(ir), ijk(ir+1)-1
	    if( lf(i)) then
	      call polar_matrix(x(i), y(i), z(i), tmp1)
	      call axbeqc(trna,tmp1,tmp2)
	      call get_polar(tmp2,xa)

	      do 100 jj=1,n_group2
		j= igroupa(jj,igr)
		jr=aa_seq(j)
		if( abs(ir-jr) .lt. iskip) goto 100
	      xb(1)=x(j)
	      xb(2)=y(j)
	      xb(3)=z(j)
		curr=polar_diff_func(xa,xb)
		if(curr.lt.dmin.or.curr.gt.dmax) goto 100
!		av=av+curr
!		sgm=sgm+curr*curr
		write(io,1001) text(i)(14:27),text(j)(14:27), curr
!		if(igo.eq.4) then	! the option 'COPY' !010126
!		  text(i)(31:54)=text(j)(31:54)
!		  b(i)=b(j)
!		  w(i)=w(j)
!		else if(igo.eq.5) then	! the option 'MARK' !990422
!		  w(i)=w(j)
!		else 
!		  w(i)=w(i)+1.
!		  endif
	        ksum=ksum+1

!		if( storea) then
!	          igroupa( ksum,1)= j		  
!	          endif

		if( ksum .ge. imax) then
	          write(cnum,'(i6)') imax
		  errmsg=
	1' errmsg: the number of atom pairs is larger than '//cnum
		  return 1
		  endif
!	        if(.not.incl) then 
!		  lf(i)=.false.		! 021218
!		  goto 110
!		endif
100	      end do 
	    end if
110	  end do
	end do
	enddo ! i=1,num_symm

200	status=ksum
!	if( ksum .eq. 0)then
	  write(6,1008) ksum
!	else
!	  sumk=1./ksum
!	  av=av*sumk
!	  sgm=sgm*sumk
!	  sgm=sqrt( max(0.,sgm-av*av))
!	  write(6,1008)  ksum, av, sgm
!	  write(io,1008) ksum, av, sgm
!	end if

1001	format(1x,2(a14,','),'d= ',f8.3)
1008	format('!psdist>',i6,' matches are found':
	1', avd=',f8.3,', sgmd=',f8.3)
1009	format(' ',a,'> # of on atoms =',i5,
	1 ', # of atoms in group ',a,' =',i5)

	call typelist(io)
	return
990	return 1
	end
	
	subroutine mcs_atm_p(*)
!  	====================
!example:
! init; zone s1-s40 ; gr s	! from self  rotation function 
! init; zone r1-r15 ; gr r	! from cross rotation function
! init; load r 
! pclique s  4 3 3

	common /cmm_mcs/ is_angle
	is_angle=1
	call mcs_atm(*900)
	return
900	return 1
	end
	
	subroutine mcs_atm_n(*)
!  	====================
	common /cmm_mcs/ is_angle
	is_angle=0
	call mcs_atm(*900)
	return
900	return 1
	end
	
	subroutine edp_rms_p(*)
chk	====================
chk	assume the PDB file contains polar angles.
chk	this subroutine calculates best rotation to minimise difference 
chk	of two sets of polar angles. 

	include 'edp_main.inc'
	include 'edp_dat.inc'
	include 'edp_file.inc'
	
	logical identical, even_weight
	character*8	s_weight
	data s_weight/'weight'/

	integer match_l
	external match_l

	logical nword
	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie

 	real s(3,max_atom) ,p(3,max_atom), weight(max_atom), shift(12)
	1 ,sumx(3) ,sumxx(3)

	common /cmm_mcs/ is_angle
	real tmp1(3,3), tmp2(3,3)

!	logical  open_file
!	external open_file
	logical  open_file1
	external open_file1

	is_angle=1

	n_of_syn=2					!000515
	syntax(1)='syntax:' 
	syntax(2)='poverlay group_id.s [filename.s] [weight] '

	igr=match_l( max_gr, cgroup )
	if( igr .le. 0) return 1

	rtn_out='rtn_.txt'
	if(.not. open_file1(29, rtn_out, 'unknown','.txt')) return 1

	if(nword(n_len,txt,ib,ie) ) then
	  even_weight=.true. 		! equally weighted
	else if(txt(ib:ie).eq.s_weight
	1   (1:min(len(s_weight),ie-ib+1))) then
	  even_weight=.false.		! weighted according to w(i)
	else
	    return 1
	    endif

chk	initial rms=1. : weighted according to w(i)
chk	initial rms=0. : equally weighted

c	rms=-2. :  edp_rms aborted
	rms=-2.

	n_group0=0
	do i=1,n_atom
	 if(lf(i)) n_group0=n_group0+1
	end do
	n_group2=n_groupa(igr)

	write(6,1009)   n_group0,cgroup(igr),n_group2
1009	format(' rms>  #of on atoms (na)=',i5,
	1 ', #of atoms in group ',a,' (ng)=',i5)

	if(n_group0.ne.n_group2.or.n_group0.eq.0) then
	  errmsg= ' errmsg: na =/= ng or na=0'
	  status= -3
	  return 1
	end if

	jj=0
	identical=.true.
	sgm=0.0
	sumw=0.0
	do k=1,3
	  sumx(k)=0.0
	  sumxx(k)=0.0
	  enddo

	kk=0
	do i=1,n_atom
	  if(lf(i)) then
	    jj=jj+1
	    j=igroupa( jj,igr)
	    if(i.ne.j) identical=.false.
	    call polar_matrix(x(i),y(i), z(i),tmp1)
	    call polar_matrix(x(j),y(j), z(j),tmp2)
	    do k1=1,3
	     kk=kk+1
	     do k2=1,3
	      s(k2,kk)=tmp2(k1,k2)
	      p(k2,kk)=tmp1(k1,k2)
	     enddo
	     if(even_weight) then 
	      weight(kk)=1.
	     else
	      weight(kk)=w(i)
	     endif

	     sumw=sumw+ weight(kk)
	     do k=1,3
	      sumx(k) =sumx(k) + (p(k,kk)-s(k,kk))    *weight(kk)
	      sumxx(k)=sumxx(k)+ (p(k,kk)-s(k,kk))**2 *weight(kk)
	     enddo
	    enddo
	    
	    endif
	  enddo

	if(sumw.le.eps()) then
	  write(6,*) '%EdPDB-W- the weight is smaller than zero.'
	  return 1
	  endif

	do k=1,3
	  sumx(k) =sumx(k) /sumw
	  sumxx(k)=sumxx(k)/sumw
	  sgm= sgm + sumxx(k) -sumx(k)**2
	  enddo

	if(identical.or.sgm .lt. eps()*10.0) then
chk	  write(6,*) 'rms> the two (weighted) structures are identical.'
	  write(29,1018) 1.,0.,0.,0.,1.,0.,0.,0.,1.
	1  ,0.0, 0.0, 0.0
!	1  ,sumx(1),sumx(2),sumx(3) 
1018	  format(3(3f12.7/),3f12.5)
	  if(verbose.ge.3) write(6,1230) rtn_out(:ltrim(rtn_out))
1230	format(' poverlay-I3> To analyze the rotation, type:'
	1 /'   axis ',a)

!	  rms=sqrt(max(0.0,sgm))
!	  shift(1)=sqrt(sumxx(1)+sumxx(2)+sumxx(3))
!	  shift(2)=0.
	  rms=0.0
	  shift(1)=0.0
	  shift(2)=0.0
	else
	  rms=-1.0
	  call edp_rms_doit(rms, n_group0,s,p,weight, shift, 29
	1  ,.false., .false.)
!	  rms=2.0*asind(0.5*rms)
	  end if

	close (29)
	if(rms.ge.0.) write(6,1010) rms, shift(2)
1010	format(
	1' poverlay> rms=',e11.4,
	1'(deg.), shifted  by 0.0, and rotated by',e11.4)

	if(rtn_out.ne.'rtn_.txt') 
	1 write(6,1011) rtn_out(:ltrim(rtn_out))
1011	format(' poverlay>  the file ',a,' is created or overwritten.')
	end
	
	subroutine dist_p(*)
chk	===============
	include 'edp_main.inc'
	include 'edp_dat.inc'
	include 'edp_file.inc'

	parameter (io=48, max_as= 2048)

	character*(max_num_chars) txt
	common /cmm_txt/ n_len,txt,ib,ie
	logical nword0

	logical storea

	character*12 s_move(5)	!, s_load, s_copy, s_prekin
	character*6 cnum

	real xa(3), xb(3)

	data s_move/'load','move','punch_all','copy','mark'/

	n_of_syn=4					!000515
	syntax(1)='syntax:' 
	syntax(2)='pdist group_id.s dmin.r dmax.r '
	syntax(3)='  [skip.i max_output.i] '
	syntax(4)='  [(load, copy, move, mark, punch_all)]' 

	igr= match_l( max_gr, cgroup)
	if( igr .le. 0) then
	  errmsg=' errmsg: a group name needed'
	  return 1
	  end if
	n_group2= n_groupa(igr)
	if( n_group2.le.0) then
	  write(6,*) 
	1 'pdist-W> UNDONE: define group ['//cgroup(igr)//'] first.'
	  return
	  end if

	n_group1=0
	do i= 1, n_atom
	  if( lf(i)) then
	    n_group1=n_group1+1
	    endif
	  enddo
	if( n_group1 .le. 0) then
	  write(6,'(a)') ' dist-W>  UNDONE: no on atom exist.'
	  return
	  end if

chk	input dmin, dmax
	dmax=0.0
	call read_ar(1, dmin, *990, *49)
	call read_ar(1, dmax, *990, *49)
49	if(dmax.lt.dmin) then
	  rsmax=dmax
	  dmax=dmin
	  dmin=rsmax
	  endif

chk	input iskip
	iskip=1	! 990422
	imax=max_atom
	storea= .false.
	igo=0

	call read_ai(1, iskip, *990, *50)
	iskip= abs( iskip)
50	if( .not. nword0(n_len,txt,ib,ie) ) then
	  if( ib  .gt.ie) then
	    imax=max_atom
	  else
	    read(txt(ib:ie),*,err=990) imax
	    if( imax .lt. 0 .or. imax .gt. max_atom) then
	      write(cnum,1062) max_atom
1062	      format(i6)
	      errmsg=
	1' errmsg: the maximum # of listout should be less than'//cnum
	      return 1
	      endif
	    endif
	  igo= match_l1(5,s_move)
	  if( igo.ne.1 .and. igo.ne.4 .and. igo.ne.5) return 1
	  storea=igo.eq.1
	  endif

	if(verbose.ge.4 ) then
	  write(6,'(a,f6.2)') ' pdist-I4> dmin=', dmin
	  write(6,'(a,f6.2)') ' pdist-I4> dmax=', dmax
	  write(6,'(a,i6)'  ) ' pdist-I4> skip=', iskip
	  write(6,'(a,i6)')   ' pdist-I4> max_output=', imax
	
	  if(igo.eq.1) then
	    write(6,*) 
	1'pdist-I> the group SCR is stuffed with selected atoms'//
	1' from group ['//cgroup( igr)//']'
	  else  if(igo.eq.5) then
	    write(6,'(a)') 
	1' pdist-I> the  occ(w) of the ON atom is replaced with '
	    write(6,'(a)') 
	1'              that of the last matched atom from the group'
	    endif
	  endif

	write(6,1009) 'pdist', n_group1, cgroup(igr), n_group2

	ksum=0
	av=0.
	sgm=0.
	do ir= 1, n_res
	  do 110 i= ijk(ir), ijk(ir+1)-1
	    if( lf(i)) then
	      xa(1)=x(i)
	      xa(2)=y(i)
	      xa(3)=z(i)
	      w(i)=0.
	      do 100 jj=1,n_group2
		j= igroupa(jj,igr)
		jr=aa_seq(j)
		if( abs(ir-jr) .lt. iskip) goto 100
	      xb(1)=x(j)
	      xb(2)=y(j)
	      xb(3)=z(j)
		curr=polar_diff_func(xa,xb)
		if(curr.lt.dmin.or.curr.gt.dmax) goto 100
		av=av+curr
		sgm=sgm+curr*curr
		write(io,1001) text(i)(14:27),text(j)(14:27), curr
		if(igo.eq.4) then	! the option 'COPY' !010126
		  text(i)(31:54)=text(j)(31:54)
		  b(i)=b(j)
		  w(i)=w(j)
		else if(igo.eq.5) then	! the option 'MARK' !990422
		  w(i)=w(j)
		else 
		  w(i)=w(i)+1.
		  endif
	        ksum=ksum+1
		if( storea) then
	          igroupa( ksum,1)= j
		  
	          endif
		if( ksum .ge. imax) then
	          write(cnum,1062) imax
		  errmsg=
	1' errmsg: the number of atom pairs is larger than '//cnum
		  return 1
		  endif
	        if(.not.incl) then 
		  lf(i)=.false.		! 021218
		  goto 110
		endif
100	      end do 
	    end if
110	  end do
	end do

200	if( storea) n_groupa(1)=ksum
	if( ksum .eq. 0)then
	  write(6,1008) ksum
!	  write(6,1007) distance
	else
	  sumk=1./ksum
	  av=av*sumk
	  sgm=sgm*sumk
	  sgm=sqrt( max(0.,sgm-av*av))
	  write(6,1008)  ksum, av, sgm
!	  write(6,1007)  distance
	  write(io,1008) ksum, av, sgm
!	  write(io,1007) distance
	end if

1001	format(1x,2(a14,','),'d= ',f8.3)
!1007	format(
!	1'                      dist. between the mass centres =',f8.3)
1008	format('!pdist>',i6,' records are calculated':
	1', avd=',f8.3,', sgmd=',f8.3)
1009	format(' ',a,'> # of on atoms =',i5,
	1 ', # of atoms in group ',a,' =',i5)

	call typelist(io)
	return
990	return 1
	end

	subroutine pickr(*)
chk	================
	include 'edp_main.inc'
!	use edp_main
	include 'edp_dat.inc'
!	use edp_dat
	logical  nword ,except
	logical  lgo
	logical lf1(max_atom)

	character*(3) resi ,key(max_num_zones) ,sexcept*6
	character*(max_num_chars) txt, txt1
	common /cmm_txt/n_len,txt,ib,ie

	character*3 res_l(max_res)
	integer  kres(max_res)
	data sexcept/'except'/
	save	!030422

	n_of_syn=4					!000515
	syntax(1)='syntax:' 
	syntax(2)='1) residue' 
	syntax(3)='2) residue res_type1.s [res_type2.s ... ] '
	syntax(4)='3) residue except res_type1.s [res_type2.s ... ]' 

	if(residue_untouch) then
	  k_res=0
	  residue_untouch=.false.
	  endif

	ie0=ie
	if(k_res.le.0) then
	  txt1=txt
          n_len1=n_len
	  lgo=.false.
	  goto 800
	  end if

400	ie=ie0
	if(nword(n_len,txt,ib,ie)) goto 850	! list res. information
	ie=ie0

	do i=1, n_atom
	  lf1(i)=.false.
	  enddo

c	igr:	id of the basckit group, 0-4
	call dfgroup(igr,*901)

	j1=n_groupa(igr)
	do j=1,j1
	  i=igroupa( j,igr)
	  lf1(i)=.true.
	  enddo

600	jj=0
	except=.false.


601	call find(k_res,res_l,kk)
	if(kk) 6011,6014,6013
6011	if( jj .eq. 0 .and. sexcept(:ie-ib+1) .eq. txt(ib:ie) ) then
	  except=.true.
	  goto 601
	  end if
	if(verbose.ge.2 ) then
	  write(6,'(a)') 
	1' residue-W2> residue ['//txt(ib:ie)//'] does not exist'
	  endif
	goto 601

6013	jj=jj+1
	if(jj.gt.max_num_zones) goto 700
	key(jj)=txt(ib:ie)
	goto 601

6014	if(jj.eq.0 ) goto 900

	do i=1,n_res
	  resi=res(i)
	  do ii=1,jj
	    if(except) then
	      if(resi.eq.key(ii)) goto 602
	    else 
	      if(resi.eq.key(ii)) goto 6012
	      end if
	    end do
	  if(.not.except) goto 602
6012	  do kk=ijk(i),ijk(i+1)-1
	    if(lf1(kk)) then
	      lf(kk)=incl
	      endif
	    end do
602	  enddo
	return
c---
700	errmsg=' errmsg:  too many res_name in the input list'
	return 1

900	errmsg=' errmsg: unrecoganizible residue name(s)'
	return 1

901	errmsg=' errmsg: wrong zone information'
	return 1

	entry beforer
c	=============
	lgo=.true.
	txt1=txt
	ie1=ie
	n_len1=n_len
800	if(k_res.le.0) then
	  res_l(1)=res(1)
	  k_res=1
	  n_len=4
	  do i=2,n_res
	    txt=res(i)
	    ie=0
	    call find(k_res,res_l,j)
	    if(j.le.0) then
	      k_res=k_res+1
	      res_l(k_res)=res(i)
	      end if
	    end do
	  end if
	txt=txt1
	n_len=n_len1
	if(.not.lgo) goto 400
	ie=ie1
	return

850	do j=1,k_res
	  kres(j)=0
	  end do

	k=0
	do 750 ii=1,n_res
	  do 740 i=ijk(ii),ijk(ii+1)-1
	    if(lf(i)) then
	      resi=res(ii)
	      do j=1,k_res
	        if(resi.eq.res_l(j)) then
		  kres(j)=kres(j)+1
		  k=k+1	
	          goto 750
	          end if
		end do
		write(6,*) '%EdPDB-F- unexpected error in residue'
		write(6,*)
		call exit( 4)
		end if
740	continue
750	continue
	write(6,1004)  n_res,k, (res_l(i),min(kres(i),999),i=1,k_res)
1004	format(' residue>  total=',i5,', selected=',i5/
	1 ,(t10,6(a4,' [',i3,']')))
	end
chk***	the end of residue.for

copyright by X. Cai Zhang
