

	subroutine get_one_atom(ia, cx, cy, cz)
chk	=======================
chk	output:
chk	ia<0	: wrong residue/atom name
chk	ia=0	: center x y z
chk	ia>	: choose atom

	include 'edp_main.inc'
!	use edp_main
	character*4 katom
	real cxyz(3)

	character*(max_num_chars) txt
	logical nword
	common /cmm_txt/n_len,txt,ib,ie

	ia=-1
	call find(n_res,ires,j)
	if(j.lt.0) then
	  if(txt(ib:ie).eq.'center') then
	    call read_ar(3,cxyz,*900,*900)
	    cx=cxyz(1)
	    cy=cxyz(2)
	    cz=cxyz(3)
	    ia=0
	    return
	  else
	    return 
	    endif
	  endif

	if(j.eq.0) then 
	  do i=1,n_atom
	    if(lf(i)) then
	      ia=i
	      goto 100
	      endif
	    enddo
	  return 
	  endif

	if(nword(n_len,txt,ib,ie) ) then
	  ia=ijk(j)
	else
	  katom=txt(ib:ie)
	  if( .not. nword(n_len,txt,ib,ie)) return 
	  do i=ijk(j),ijk(j+1)-1
	    if(katom.eq.atom(i)) then
	      ia=i
	      goto 100
	      endif
	    enddo
	  errmsg=' errmsg: no atom with the given name found '
	  return 
	  endif
100	cx=x(ia)
	cy=y(ia)
	cz=z(ia)
900	end

	subroutine nayb(*)
chk	===============
	
c	nayb selects the atoms in a sphere centered at a given atom of a given
c	residue with a specified radias.
c 	the command for nayb is
c	nayb radias, residue id, atom_name
 
	include 'edp_main.inc'
!	use edp_main

	n_of_syn=3					!000515
	syntax(1)='syntax:' 
	syntax(2)='1) nayb radius.r [res_id.s [atom_name.s]]' 
	syntax(3)='2) nayb radius.r center x.r y.r z.r'
 
c	igr:	id of the basckit group, 0-4
	call dfgroup(igr,*901)

chk	input the search_radius
	call read_ar(1, rad, *902, *902)
 	
	call get_one_atom( ia, cx,cy,cz)
	if( ia.lt.0 ) return 1
1002	format(' nayb-I> the occ will be changed to the dist.')

	if(incl) write(6,1002) 
	r=rad*rad

	j1=n_groupa(igr)
	do 200 j=1,j1
	  i=igroupa( j,igr)
	dx=cx-x(i)
	dy=cy-y(i)
	dz=cz-z(i)
	if(dx.gt.rad) goto 200
	if(dy.gt.rad) goto 200
	if(dz.gt.rad) goto 200
	d=dx*dx+dy*dy+dz*dz
	if(d.gt.r) goto 200
	lf(i)=incl
	if(incl) w(i)=sqrt(d)
200	continue
	return

901	errmsg=' errmsg: wrong group/zone information'
902	return 1
	end

copyright by X. Cai Zhang
	subroutine pattern(*)
chk	=================
	include 'edp_main.inc'
!	use edp_main
	include 'edp_dat.inc'
!	use edp_dat

	character*(max_num_chars) tmp
	character*(1) template(32,max_num_chars), ic
	integer i_template(max_num_chars)

	character*(max_num_chars) txt
	common /cmm_txt/ n_len,txt,ib,ie
	logical nword 

	integer get_pattern
	external get_pattern

	n_of_syn=2					!000515
	syntax(1)='syntax:' 
	syntax(2)='pattern pattern.s [position.i]' 

	call dfgroup(igr0,*901)

	if(nword(n_len,txt,ib,ie)) return 1
	tmp=txt(ib:ie)
	
	if(get_pattern(tmp,template, i_template, n) .ne. 0 ) return 1
	n_g=n_groupa(igr0)-n

	k=22
	call read_ai(1, k, *902, *50)

50	do 100 i=0,n_g
	  do j=1,n
	    ic=text(igroupa(i+j,igr0))(k:k)
	    m1=i_template(j)
		if( m1.ge.0 ) then
	     do m=1,m1
		   if( template(m,j).eq.wildcard1) goto 60
		   if( template(m,j).eq.ic) goto 60
	     enddo
		 goto 100
		else 
	     do m=1,-m1
		   if( template(m,j).eq.ic) goto 100
	     enddo
		endif
60 	    enddo
	  do j=1,n
	    lf(igroupa(i+j,igr0))=incl	    
	    enddo
100	  enddo

	return
901	errmsg= 'errmsg: wrong group/zone information'
902	return 1
	end

	subroutine align3d(*)
chk	=====================
	include 'edp_main.inc'
!	use edp_main
	include 'edp_dat.inc'
! 	use edp_dat

	parameter (max_n=max_l)
	parameter (max_2n=max_n*2)
!	integer, parameter :: max_n=max_l
!	integer, parameter :: max_2n=max_n*2

	integer   ic1(max_2n), ic2(max_2n)
	common /cmm_align3d/ name1, name2, dist_cutoff, d2

	character*(14) text_i, text_j
	character*(max_num_chars) txt
	common /cmm_txt/ n_len,txt,ib,ie

!	save igr1, igr2, n_group0, n_group1, penalty_for_break	
	save	!030422
	
	n_of_syn=3					!000515
	syntax(1)='syntax:' 
	syntax(2)='load group2.s | align3d group1.s sub_group1.s '
	syntax(3)='  [distance_cutoff.r [penalty_to_break.r]] '

	call dfgroup(igr0,*901)
	n_group0=n_groupa(igr0)
	if(n_group0.gt.max_n) then
	  errmsg= ' errmsg: too many atoms to select'
	  return 1
	  endif

	igr1= match_l( max_gr, cgroup)
	if( igr1 .le. 0) then
	  errmsg=' errmsg: a group name is needed'
	  return 1
	  end if
	n_group1= n_groupa(igr1)
	if( n_group1.le.0) then
	  write(6,*) 
	1 'align3d-W> UNDONE: define group ['//cgroup(igr1)//'] first.'
	  return
	else if(n_group1.gt.max_n) then
	  errmsg=' errmsg: too many atoms in group ['//cgroup(igr1)//']'
	  return 1
	  end if

	call find_a_group(igr2)
	if( igr2 .le. 0) return 1

	dist_cutoff=2.0
	call read_ar(1, dist_cutoff, *902, *50)
50	d2 = dist_cutoff**2
	penalty_for_break = -d2
	call read_ar(1, penalty_for_break,  *902, *60)

60	name1=igr0	! for sc_align3d
	name2=igr1
	call edp_needleman(n_group0,n_group1,penalty_for_break)
	return

901	errmsg= ' errmsg: wrong group/zone information'
	return 1
902	errmsg= ' errmsg: wrong distance cutoff'
	return 1

	entry output_ndlmn(ic1,ic2,k,smax,kgap1,kgap2)
chk	==================
	n_g=0

	call e_typelist()

	do k0=k,1,-1
	  i = ic1(k0)
	  j = ic2(k0)

	  if(   i.gt.0
	1 .and. j.gt.0) then
	    d = sc_align3d(i,j)
	    i=igroupa(i,name1)
	    j=igroupa(j,name2)
 	    if( d .gt. 0.0) then
	      n_g=n_g+1
	      lf(i) = incl
	      igroupa(n_g,igr2) = j
	      if(verbose.ge.4 ) then
        write(48,1047) text(i)(14:27), text(j)(14:27), sqrt(d2-d)
1047    format(' align3d-I4> ',2(a14,','),'d= ',f8.3)
	        endif
	    else if(verbose.ge.5 ) then
        write(48,1057) text(i)(14:27), text(j)(14:27), sqrt(d2-d)
1057    format(' align3d-I5> ',2(a14,','),'d=>',f8.3)
	      endif	! d > 0 ?
	  else if(verbose.ge.6 ) then
	    text_i=' '
	    text_j=' '
	    if(i.gt.0) then
	      i=igroupa(i,name1)
	      text_i=text(i)(14:27)
	    else if(j.gt.0) then
	      j=igroupa(j,name2)
	      text_j=text(j)(14:27)
	      endif    
        if(text_i.eq.' ') then 
	  write(48,1067) text_i, text_j,-999.999
	else
	  write(48,1067) text_i, text_j, 999.999
	  endif
1067    format(' align3d-I6> ',2(a14,','),'d= ',f8.3)
	    endif ! i < 0 or j < 0
	  enddo
	  
	n_groupa(igr2)= n_g

500	if(verbose.ge.2 ) write(6,1008)    dist_cutoff, penalty_for_break
1008	format(' align3d-I2> distance_cutoff =',f5.1
	1,', penalty_for_a_gap =',f8.1)
	write(6,1009)   
	1               n_group0, kgap1
	1 ,cgroup(igr1),n_group1, kgap2
	1 ,n_g,smax
1009	format(
	1' align3d> # of atoms selected from =', i4,', # of gaps =',i4/
	1' align3d> # of atoms in group ',a,' =',i4,', # of gaps =',i4/
	1' align3d> # of atoms matched       =', i4,', with score=',f8.1)
	if(n_g.le.0 ) then
	  status= -3
	else
	  status= n_g
	endif

	end

	function sc_align3d(i0,j0)
chk	=================
	include 'edp_main.inc'
! 	use edp_main
	common /cmm_align3d/ name1, name2, d, d2

	i = igroupa(i0,name1)
	j = igroupa(j0,name2)

	sc_align3d=0.0

	dx = x(i)-x(j)
	if(abs(dx).gt.d) return
	dy = y(i)-y(j)
	if(abs(dy).gt.d) return
	dz = z(i)-z(j)
	if(abs(dz).gt.d) return
	rr = dx**2 + dy**2 + dz**2
	if(rr.gt.d2) return
	sc_align3d = d2 - rr
	end

	subroutine EDP_needleman(n1,n2,penalty_for_break)
chk	========================
	include 'edp_dim.inc'
!	use edp_dim
	include 'edp_dat.inc'
!	use edp_dat

	parameter (max_n=max_l)
	parameter (max_2n=max_n*2)
!	integer, parameter :: max_n=max_l
!	integer, parameter :: max_2n=max_n*2
c	integer ia0(max_n) ,ib0(max_n) 

c	real ad(max_n), bd(max_n)

c	! define function sc_align3d(i,j)

	DIMENSION SM(MAX_N,MAX_N), TK_A(MAX_N)
	integer   IC1(MAX_2N), IC2(MAX_2N)
	integer   IGO(MAX_N,MAX_N),  KGO_A(MAX_N)

	if(verbose.ge.6 ) write(6,1069) 	!

1069	format(' needleman-I6> reference:'
	1/'  Needleman SB, Wunsch CD.'
	1/'  A general method applicable to the search for similarities '
	1/'    in the amino acid sequence of two proteins.'
	1/'  J Mol Biol. 1970 Mar;48(3):443-53.')

	pe=0.1		! penalty_for_extrend_a_break
c	do i=1,max_n
c	  ia0(i)=i
c	  ib0(i)=i
c	  ad(i)=penalty_for_break
c	  bd(i)=penalty_for_break
c	  enddo

	do i=1,max_2n
	  ic1(i)=0
	  ic2(i)=0
	  enddo

c	IBJ=ib0(1)
	DO I= 1, N1
c	  SM(I,1)= sc_align3d(ia0(I),IBJ)
	  SM(I,1)= sc_align3d(i,1)
	  IGO(I,1)=0
	  END DO

c	IAI=ia0(1)
	DO J= 2, N2
c	  SM(1,J)= sc_align3d(IAI, ib0(J))
	  SM(1,J)= sc_align3d(1,j)
	  IGO(1,J)=0
	  END DO

c	IAI=ia0(2)
c	IBJ=ib0(2)
c	SM(2,2)=SM(1,1)+sc_align3d(IAI,IBJ)
	SM(2,2)=SM(1,1)+sc_align3d(2,2)
	IGO(2,2)=0

c	BDJ=BD(1)
	BDJ=penalty_for_break
	PEBDJ=PE*BDJ
	TL=-1024.
	DO I=3,N1
		JGO=0
		SMAX=SM(I-1,1)

		STMP=SM(I-2,1)+BDJ
		IF(SMAX.LT.STMP) THEN
		  SMAX=STMP
		  JGO=2-I
		END IF
		
		IF(SMAX.LT.TL) THEN
		  JGO=LGO
		  SMAX=TL
		  TL=TL+PEBDJ
		ELSE IF(TL.LT.STMP) THEN
		  TL=STMP+PEBDJ
		  LGO=2-I
		ELSE
		  TL=TL+PEBDJ
		END IF
	
c		SM(I,2)=SMAX+sc_align3d(ia0(I),IBJ)
		SM(I,2)=SMAX+sc_align3d(i,2)
		IGO(I,2)=JGO
	END DO

c	BDJ=AD(1)
	BDJ=penalty_for_break
	PEBDJ=PE*BDJ
	TL=-1024.
	DO I=3,N2
		JGO=0
		SMAX=SM(1,I-1)

		STMP=SM(1,I-2)+BDJ
		IF(SMAX.LT.STMP) THEN
		  SMAX=STMP
		  JGO= I-2
		END IF
		
		IF(SMAX.LT.TL) THEN
		  JGO=LGO
		  SMAX=TL
		  TL=TL+PEBDJ
		ELSE IF(TL.LT.STMP) THEN
		  TL=STMP+PEBDJ
		  LGO=I-2
		ELSE
		  TL=TL+PEBDJ
		END IF
	
c		SM(2,I)=SMAX+sc_align3d(IAI,ib0(I))
		SM(2,I)=SMAX+sc_align3d(2,i)
		IGO(2,I)=JGO
	END DO
C
	DO I=3,N1
		TK_A(I)=-1024.
	END DO
	DO J=3,N2
		J1=J-1
		J2=J-2
		J3=J-3
c		IBJ=ib0(J)
c		BDJ=BD(J)
		BDJ=penalty_for_break
		PEBDJ=PE*BDJ
		
		TL=-1024.
		DO I=3,N1
			I1=I-1
			I2=I-2
			I3=I-3
c			ADI=AD(I)
			ADI=penalty_for_break
			PEADI=PE*ADI

			JGO=0
			SMAX=SM(I1, J1)

			STMP=SM(I2,J1)+BDJ
			IF(SMAX.LT.STMP) THEN
			  SMAX=STMP
			  JGO=-I2
			END IF
		
			IF(SMAX.LT.TL) THEN
			  JGO=LGO
			  SMAX=TL
			  TL=TL+PEBDJ
			ELSE IF(TL.LT.STMP) THEN
			  TL=STMP+PEBDJ
			  LGO=-I2
			ELSE
			  TL=TL+PEBDJ
			END IF

			STMP=SM(I1,J2)+ADI
			IF(SMAX.LT.STMP) THEN
			  SMAX=STMP
			  JGO=J2	! why not -j2
			END IF
		
			TK=TK_A(I)
			KGO=KGO_A(I)
			IF(SMAX.LT.TK) THEN
			  JGO=KGO
			  SMAX=TK
			  TK=TK+PEADI
			ELSE IF(TK.LT.STMP) THEN
			  TK=STMP+PEADI
			  KGO=J2	! why not -j2
			ELSE
			  TK=TK+PEADI
			END IF
			TK_A(I)=TK
			KGO_A(I)=KGO

c			SM(I,J)=SMAX+sc_align3d(ia0(I),IBJ)
			SM(I,J)=SMAX+sc_align3d(i,j)
			IGO(I,J)=JGO
		END DO
	END DO

C
	SMAX=SM(N1,N2)
	IMAX=N1
	JMAX=N2
	DO I=1,N1
		IF(SMAX.LT.SM(I,N2)) THEN
			IMAX=I
			SMAX=SM(I,N2)
		END IF
	END DO

 	DO J= 1, N2
		IF(SMAX.LT.SM(N1,J)) THEN
			JMAX=J
			SMAX=SM(N1,J)
		END IF
	END DO

c	ratio=SMAX/MIN(N1,N2)

	K=0
	IF(JMAX.LT.N2) THEN
		IMAX=N1
		DO J=N2, JMAX+1, -1
			K=K+1
c			IC2(K)=ib0(J)
			IC2(K)=j
		END DO
	ELSE 
		DO I=N1, IMAX+1, -1
			K=K+1
c			IC1(K)=ia0(I)
			IC1(K)=i
		END DO
	END IF

C	'Search maximun score'
	kgap1=0
	kgap2=0

	DO WHILE(IMAX.GE.1.AND.JMAX.GE.1)
		K=K+1
c		IC1(K)=ia0(IMAX)
c		IC2(K)=ib0(JMAX)
		IC1(K)=imax
		IC2(K)=jmax
		JGO=IGO(IMAX,JMAX)
		IF(JGO) 100, 200, 300

100	KGAP2=KGAP2+1
	IS=IMAX-1
	IMAX=1-JGO
	DO I=IS,IMAX,-1
	K=K+1
c	IC1(K)=ia0(I)
	IC1(K)=i
	END DO
	GOTO 200

300	KGAP1=KGAP1+1
	JS=JMAX-1
	JMAX=JGO+1
	DO J= JS,JMAX,-1
	K=K+1
c	IC2(K)=ib0(J)
	IC2(K)=j
	END DO

200	IMAX=IMAX-1
	JMAX=JMAX-1
	END DO

	IF(IMAX.GT.1) THEN
		DO I=IMAX,1,-1
			K=K+1
c			IC1(K)=ia0(I)
			IC1(K)=i
		END DO
	ELSE IF(JMAX.GT.1) THEN
		DO J=JMAX,1,-1
			K=K+1
c			IC2(K)=ib0(J)
			IC2(K)=j
		END DO
	END IF
	call output_ndlmn(ic1,ic2,k,smax,kgap1,kgap2)
	end

	subroutine naybr(*)
chk	================
	include 'edp_main.inc'
!	use edp_main

	n_of_syn=2					!000515
	syntax(1)='syntax:' 
	syntax(2)='naybr distance.r res_id.s '

c	igr:	id of the basckit group, 0-4
	call dfgroup(igr,*901)

chk	input the search radius
	call read_ar(1, rad, *902, *902)

chk	input the residue_id
	call find(n_res,ires,nbr)
	if(nbr.le.0) return 1

	rad=rad*rad
	do i=ijk(nbr),ijk(nbr+1)-1
	  cx=x(i)
	  cy=y(i)
	  cz=z(i)

	  j1=n_groupa(igr)
	  do jr=1,j1
	    j=igroupa( jr,igr)
	    if(lf(j).neqv.incl) then
	      if((cx-x(j))**2+(cy-y(j))**2+(cz-z(j))**2.le.rad)
	1       lf(j)=incl
	      endif
	    enddo
	  enddo
	return

901	errmsg=' errmsg: wrong group/zone information'
902	return 1
	end

copyright by X. Cai Zhang

	subroutine ncac(trn,*)
chk	===============
	include 'edp_main.inc'
!	use edp_main
	include 'edp_dat.inc'
!	use edp_dat

	character*4 atom_name(3) 
	dimension trn(3,4), v1(3), v2(3), v3(3), u1(3), u2(3),
	1 u3(3), j_a(3,2), nr(3)
	logical no_2nd_group

	character*(max_num_chars) txt
	logical nword 
	common /cmm_txt/n_len,txt,ib,ie
	save	!030422

	iatom=1
	atom_name(1)='ca'
	atom_name(2)='n'
	atom_name(3)='c'
	nr(1)=0
	nr(2)=0
	nr(3)=0
	no_2nd_group=.false.

40	call find(n_res,ires,i)
	if(i.lt.0) then
	  write(6,*) 'rtn-W> wrong res._id'
	  goto 999
	else if( i.eq.0) then 
	  if(iatom.eq.1) then
	    write(6,*) 'rtn-W> res._id not exist'
	    goto 999
	  else
	    no_2nd_group=.true.
	    do ii=1,6
	      if(.not.nword(n_len,txt,ib,ie)) goto 999
	      enddo
	    goto 44
	    endif
	  endif

	do ii=1,3
	  if(.not.nword(n_len,txt,ib,ie)) atom_name(ii)=txt(ib:ie)
	  end do

	do ii=1,3
	  if(.not.nword(n_len,txt,ib,ie)) 
	1   read( txt(ib:ie),*,err=901) nr(ii)
!1061	  format(i<ie-ib+1>)
	  if(i+nr(ii).le.0.or.i+nr(ii).gt.n_res) goto 901
	  enddo

	do ii=1,3
	  j0=ijk(i+nr(ii))
	  j1=ijk(i+nr(ii)+1)-1
	  do j=j0, j1
	    if(atom(j).eq.atom_name(ii)) goto 43
	    enddo
	  goto 901
43	  j_a(ii,iatom)=j
	  enddo

	if(iatom.eq.1) then
	  iatom=2
	  goto 40
	  endif
c-
44	if(no_2nd_group) then
	  v1(1)=1.
	  v1(2)=0.
	  v1(3)=0.
	  v2(1)=0.
	  v2(2)=1.
	  v2(3)=0.
	  v3(1)=0.
	  v3(2)=0.
	  v3(3)=1.
	else
	  tmp= axbxc(j_a(1,2),v1,v2,v3,ierr)
	  endif
	if(ierr.ne.0) goto 999
	  tmp= axbxc(j_a(1,1),u1,u2,u3,ierr)
	if(ierr.ne.0) goto 999

1004	format(3f8.3)
	trn(1,1)= v1(1)*u1(1)+ v2(1)*u2(1)+ v3(1)*u3(1)
	trn(1,2)= v1(1)*u1(2)+ v2(1)*u2(2)+ v3(1)*u3(2)
	trn(1,3)= v1(1)*u1(3)+ v2(1)*u2(3)+ v3(1)*u3(3)
	trn(2,1)= v1(2)*u1(1)+ v2(2)*u2(1)+ v3(2)*u3(1)
	trn(2,2)= v1(2)*u1(2)+ v2(2)*u2(2)+ v3(2)*u3(2)
	trn(2,3)= v1(2)*u1(3)+ v2(2)*u2(3)+ v3(2)*u3(3)
	trn(3,1)= v1(3)*u1(1)+ v2(3)*u2(1)+ v3(3)*u3(1)
	trn(3,2)= v1(3)*u1(2)+ v2(3)*u2(2)+ v3(3)*u3(2)
	trn(3,3)= v1(3)*u1(3)+ v2(3)*u2(3)+ v3(3)*u3(3)

	if(no_2nd_group) then
	  x0=0.
	  y0=0.
	  z0=0.
	else
	  x0=x(j_a(1,2))
	  y0=y(j_a(1,2))
	  z0=z(j_a(1,2))
	  endif
	j1=j_a(1,1)
	trn(1,4)= x0 -(trn(1,1)*x(j1)+trn(1,2)*y(j1) +trn(1,3)*z(j1))
	trn(2,4)= y0 -(trn(2,1)*x(j1)+trn(2,2)*y(j1) +trn(2,3)*z(j1))
	trn(3,4)= z0 -(trn(3,1)*x(j1)+trn(3,2)*y(j1) +trn(3,3)*z(j1))
	return
901	write(6,*) 'rtn-W> UNDONE: wrong parameter'
	return 1
999	write(6,*) 'rtn-W> UNDONE: error during calc. the matrix.'
	return 1
	end

copyright by X. Cai Zhang

	logical function add(j_a,cj,rj,aj,tj)
chk	====================
	include 'edp_main.inc'
!	use edp_main
	dimension j_a(3), cj(3), v1(3),v2(3),v3(3)
	1 ,v10(3),v20(3),v30(3)
	real axbxc	!010318

	add=.false.
	goto 10

	entry axbxc(j_a,v10,v20,v30, ierr) 
c	===========
	ierr=1
	axbxc=0.0	!010318
	add=.true.

10	j1=j_a(1)
	j2=j_a(2)
	j3=j_a(3)

	v3(1)=x(j2)-x(j1)
	v3(2)=y(j2)-y(j1)
	v3(3)=z(j2)-z(j1)
	r= v3(1)*v3(1)+v3(2)*v3(2)+v3(3)*v3(3)
	if(r.lt.1.e-6) return
	r=1./sqrt(r)
	v3(1)=v3(1)*r
	v3(2)=v3(2)*r
	v3(3)=v3(3)*r

	v2(1)=x(j3)-x(j1)
	v2(2)=y(j3)-y(j1)
	v2(3)=z(j3)-z(j1)
	
	v1(1)=v2(2)*v3(3)-v3(2)*v2(3)
	v1(2)=v2(3)*v3(1)-v3(3)*v2(1)
	v1(3)=v2(1)*v3(2)-v3(1)*v2(2)
	r= v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3)
	if(r.lt.1.e-6) return
	r=1./sqrt(r)
	v1(1)=v1(1)*r
	v1(2)=v1(2)*r
	v1(3)=v1(3)*r

	v2(1)=v3(2)*v1(3)-v1(2)*v3(3)
	v2(2)=v3(3)*v1(1)-v1(3)*v3(1)
	v2(3)=v3(1)*v1(2)-v1(1)*v3(2)

	if(add) then
	do i=1,3
	v10(i)=v1(i)
	v20(i)=v2(i)
	v30(i)=v3(i)
	end do
	ierr=0
	return	
	end if

	a1=rj*sind(aj)*sind(tj)
	a2=rj*sind(aj)*cosd(tj)
	a3=rj*cosd(aj)

	cj(1)= v1(1)*a1+v2(1)*a2+v3(1)*a3 +x(j1)
	cj(2)= v1(2)*a1+v2(2)*a2+v3(2)*a3 +y(j1)
	cj(3)= v1(3)*a1+v2(3)*a2+v3(3)*a3 +z(j1)
	add=.true.
	end 

	subroutine newx(*)
chk	===============
	include 'edp_main.inc'
!	use edp_main
	include 'edp_file.inc'
!	use edp_file

	character*4 atoma(3), sxyz*3
	data atoma/3*' '/, sxyz/'xyz'/
	dimension ka(3), j_a(3), cj(3)
	data ka/3*0/
	logical la(3), add, lxyz
	data la /3*.true./
c-
	logical nword, nword0
	character*(max_num_chars) txt
	common /cmm_txt/ n_len,txt,ib,ie
	save	!030422

	n_of_syn=2					!000515
	syntax(1)='syntax:' 
	syntax(2)='newxyz [(a, b, c)]' 

	do i=1,3
	 if(atoma(i).eq.' ') then
	  write(6,*) 
	1'newx-W> UNDONE: define atom1 -  atom3 first, using dfnewx.'
	  return
	 end if
	end do

	if(nword(n_len,txt,ib,ie)) then
	 ist=0
	else if(txt(ib:ie).eq.'a') then
	 ist=1
	else if(txt(ib:ie).eq.'b') then
	 ist=2
	else if(txt(ib:ie).eq.'c') then
	 ist=3
	else
	 return 1
	end if

	if(nword(n_len,txt,ib,ie)) then
	 lxyz=.false.
	else if(txt(ib:ie).eq.sxyz(:ie-ib+1) ) then
	 lxyz=.true.
	else 
	 return 1
	endif
c-
c	call openlist(48)
c-
10	ksum=0
	ii0=min(ka(1),ka(2),ka(3))
	ii1=max(ka(1),ka(2),ka(3))	
	i0=1-ii0
	i1=n_res-ii1
c
	do 800 i=i0,i1
	 do ii=1,3
c
c-	  atoma
c
	  j0=ijk(i+ka(ii))
	  j1=ijk(i+ka(ii)+1)-1
	  do j=j0,j1
	   if(.not.la(ii).or.lf(j)) then
	    if(atom(j).eq.atoma(ii)) goto 100
	   end if
	  end do
	  goto 800
100	  j_a(ii)=j
	 end do
c
c-	calculate angle
c
	if( add(j_a,cj,rj,aj,tj)) then
	  ksum=ksum+1
	  if(ist.le.0) then
chk	    channel 4 is an open PDB file.
	    if(pdb_out.ne.'?') 
	1     write(4,1071)  text(j_a(1))(1:30), cj, 0.0, 99.99
!1071	    format('ATOM   ',a23,3f8.3,2f6.2)
1071	    format(a30,3f8.3,2f6.2)
	  else
	    igroupa( ksum,1)=j_a(ist)
	    if(lxyz) then 
	      x(j_a(ist))= cj(1)
	      y(j_a(ist))= cj(2)
	      z(j_a(ist))= cj(3)
	      endif
	    write( text(j_a(ist))(31:54),1061) cj
1061	    format(3f8.3)
	    endif
	  endif

800	continue
	if(ist.gt.0) n_groupa(1)=ksum
	write(6,1009) ksum
1009	format(' newx> ',i6,' records are calculated')
	return

	entry dfnewx(*)
c	=============
	n_of_syn=5					!000515
	syntax(1)='syntax:'
	syntax(2)='1) dfnewxyz' 
	syntax(3)='2) dfnewxyz atom_a.s atom_b.s atom_c.s '
	syntax(4)=
	1'    [reg_a.i reg_b.i reg_c.i] [status_a.l status_b.l status_c.l]' 
	syntax(5)=
	1'    [distance.r] [angle.r] [torsion_angle.r]' 

	rj=0.0
	aj=0.0
	tj=0.0
c
c	atoma
c
	do ii=1,3
	  if(nword(n_len,txt,ib,ie)) then
	    if(ii.gt.1) return 1
	    goto 901
	    end if
	  atoma(ii)=txt(ib:ie)
	  ka(ii)=0
	  la(ii)=.true.
	  end do
c
c	ka
c
	do ii=1,3
	  if(nword0(n_len,txt,ib,ie)) return
	  if(ib .le. ie) read( txt(ib:ie),*,err=900) ka(ii)
	  end do
c
c	la
c
	do ii=1,3
	  if(nword0(n_len,txt,ib,ie)) return
	  if(ib .le. ie) read( txt(ib:ib),1063,err=900) la(ii)
1063	  format(l1)
	  end do
c
c	rj, aj, tj
c
	call read_ar(1, rj, *900, *851)
851	call read_ar(1, aj, *900, *852)
852	call read_ar(1, tj, *900, *853)
853	return
900	return 1
c
	entry shnewx
c	===========
901	write(6,1002) atoma, ka, la, rj, aj, tj
1002	format(' newx>',t10,3a5,t30,3i5,t55,3(l1,1x),t65,f5.1,2f5.0)
	end

copyright by X. Cai Zhang


	subroutine order(*)
chk	================
	include 'edp_main.inc'
!	use edp_main
	real rl(max_atom)
	integer addess(max_atom), jorder(max_atom)

	logical fo
	data fo/.true./

	parameter  (num_lab=8)
!	integer, parameter  :: num_lab=7
	character*6 label(num_lab)
	data label/'w','-w','b','-b','dfres','swap','load','id'/

	n_of_syn=8					!000515
	syntax(1)='syntax:' 
	syntax(2)='1) sort'
	syntax(3)='2) sort [-]b'
	syntax(4)='3) sort [-]w'
	syntax(5)='4) sort dfres'
	syntax(6)='5) sort swap'
	syntax(7)='6) sort load [group_id(s).s]'
	syntax(8)='7) sort id'

	call find1(num_lab,label,i)
	if(i.lt.0) return 1
	if(i.eq.0) then
	  do j=1,n_atom
	    iorder(j)=j
	    end do
	  return
	  end if
	
	goto (100,100,100,100,500,600,700,100) i

100	k=0
	do jj=1,n_atom
	  j=iorder(jj)
	  if(lf(j)) then
	    k=k+1
	    if(i.eq.1) then
	      rl(k)=w(j)
	    else if(i.eq.2) then
	      rl(k)=-w(j)
	    else if(i.eq.3) then
	      rl(k)=b(j)
	    else if(i.eq.4) then	
	      rl(k)=-b(j)
	    else if(i.eq.8) then
	      read(text(j)(23:26),1001,err=990) rl(k)
	    endif
	    addess(k)=j
	  endif
	enddo
1001	format(f4.0)
	
	call sort_rl(k,rl,jorder)

	k=0
	do jj=1,n_atom
	  j=iorder(jj)
	  if(lf(j)) then
	    k=k+1
	    iorder(jj)=addess(jorder(k))
	  endif
	enddo
	return

500	if(fo) call std(fo)	! the first read stdpdb.txt only.
	call std(fo)
	return

600	call swap_on_gr(*990)
	return

700	call sort_by_loading_order(*990)
	return

990	return 1
	end

	subroutine sort_by_id(*)
	include 'edp_main.inc'
!	use edp_main
	logical lf1(max_atom)
	integer  order_tmp(max_atom)

	end 

	subroutine sort_by_loading_order(*)
chk	================================
	include 'edp_main.inc'
!	use edp_main
	logical lf1(max_atom)
	integer  order_tmp(max_atom)

	do i=1,n_atom
	  lf1(i)=.false.
	  order_tmp(i)=iorder(i)
	  enddo

	ist= match_l( max_gr, cgroup)
	if( ist .le. 0) return 1

	k=0
	do while (ist.gt.0)
	  jje=n_groupa(ist)
	  do  jj=1,jje
	    j=igroupa( jj,ist)
            if(.not.lf1(j)) then
	      lf1(j)=.true.
	      k=k+1
	      iorder(k)=j
	      endif
            end do
          ist= match_l( max_gr, cgroup)
          enddo
	
	do j=1,n_atom
	  i=order_tmp(j)
	  if(.not.lf1(i)) then
	    k=k+1
	    iorder(k)=i
	    endif
	  enddo

	if(k.ne.n_atom) then
	  errmsg=' errmsg: unexpect error in sort_by_loading_order()'
	  return 1
	  endif
        return
	end

	subroutine swap_on_gr(*)
chk	====================
chk	Swap the output/display order of the ON records with a given group
chk	  of records.

chk	suggestions:
chk	the ON      records should be in a block in the input sequence.
chk	the grouped records should be in a block in the input sequence.
	include 'edp_main.inc'
!	use edp_main

	logical a_on(max_atom), a_gr(max_atom) 
	integer  order_tmp(max_atom)
 
	k0= -1
	k1= 0
	do j=1,n_atom			! get the range of ON atoms
	  order_tmp(j)=iorder(j)
	  i=iorder(j)
	  if(Lf(i)) then
	    if(k0.le.-1) k0=j
	    k1=j
	    a_on(i)=.true.
	  else
	    a_on(i)=.false.
	    endif
	  enddo
	if(k1.le.0) then
	  errmsg= ' errmsg: no ON atom exists.'
	  return 1
	  endif

	call swap(*990)

	l0=-1
	l1= 0			
	do j=1,n_atom			! get the range of grouped atoms
	  i=iorder(j)
	  if(Lf(i)) then
	    a_gr(i)=.true.
	    if(l0.le.-1) l0=j
	    l1=j
	  else
	    a_gr(i)=.false.
	    endif
	  enddo

	if(l1.le.0) then
	  errmsg= ' errmsg: the swapped group is empty.'
	  return 1
	  endif

	j0=min(k0,l0)
	j1=max(k1,l1)
	k=j0-1

	if(j0.eq.k0) then
	  do j=l0,l1			! output grouped atoms(old) first
	    i=iorder(j)
	    if(a_gr(i)) then
	      k=k+1
	      order_tmp(k)=i
	      endif
	    enddo
	
	  do j=j0,j1				
	    i=iorder(j)
	    if(.not.a_on(i).and..not.a_gr(i)) then
	      k=k+1
	      order_tmp(k)=i
	    endif
	  enddo

	  do j=k0,k1
	    i=iorder(j)
	    if(a_on(i).and..not.a_gr(i)) then
	      k=k+1
	      order_tmp(k)=i
	      endif
	    enddo
	else
	  do j=k0,k1			! output ON(old) atoms first
	    i=iorder(j)
	    if(a_on(i).and..not.a_gr(i)) then
	      k=k+1
	      order_tmp(k)=i
	      endif
	    enddo

	  do j=j0,j1
	    i=iorder(j)
	    if(.not.a_on(i).and..not.a_gr(i)) then
	      k=k+1
	      order_tmp(k)=i
	    endif
	  enddo

	  do j=l0,l1
	    i=iorder(j)
	    if(a_gr(i)) then
	      k=k+1
	      order_tmp(k)=i
	      endif
	    enddo
	  endif

	do j=1,n_atom
	  iorder(j)=order_tmp(j)
	  enddo
	return
990	return 1
	end

	subroutine std(fo)
chk	=================

	include 'edp_main.inc'
!	use edp_main
	include 'edp_dat.inc'
!	use edp_dat
	include 'edp_file.inc'
!	use edp_file

	parameter (max_latom=50,max_lres=50)
!	integer, parameter :: max_latom=50,max_lres=50

	dimension kfk(6,max_lres), ifk(max_lres), cor(max_lres), j_a(6)
	1, t_rng(2,max_lres)
	logical cor, L_fk
	character*12 corr(2), atomj*4
	data corr/', corrected',' '/

	parameter (num_lab=2)
!	integer, parameter :: num_lab=2
	character*5 label(num_lab)
	data label/'dfres','dfork'/

	integer	iaa(max_lres), jorder(max_lres)
	real  	maa(max_lres)

	character*1 sequ2(max_res), si, laa1(max_lres)
	character*3 laa(max_lres) ,kres, sequ1(max_res)
	character*4 kaa(max_latom,max_lres), katom
	character*5 res_id

	data iaal /0/
	data iaa/max_lres*-1/
	data maa/max_lres*0.0/
	data laa/max_lres*' '/
	data laa1/max_lres*'u'/

	character*(108) file ,file1 
	character*(32) form 

	Logical check_fk
	external check_fk

	logical nword, nword0
	character*(max_num_chars) txt, txt0
	common /cmm_txt/ n_len,txt,ib,ie

	logical fo
	logical  open_file
	external open_file
	
	logical tolower
	common /main_para/ tolower, jnk(3)

	data  k_res, k_fk/1,1/

	real mw	! for molecular weight
	real amw(max_res)
	integer imw(max_res)
!	save l
	save	!030422
	
	if(.not.fo) goto 107

	txt0=txt 
	n_len0=n_len
	ib0=ib
	ie0=ie

	fo=.false.
	std_in=pdbstd_dat

80	if(.not.open_file(8,std_in,'old','.txt')) then
	  if(std_in .ne. pdbstd_dat) goto 106
	  std_in=edp_data(:ltrim(edp_data))//pdbstd_dat
	  goto 80
	  endif

c	if(i_pipe.gt.0) write(6,*)	!000427
c	1 '%EdPDB-W- the current command should not ',
c	1 'be followed by a pipe structure.'
	ij=0
102	call read_txt(8,*1061)
	call find1( num_lab,label,i)
	ij=ij+1

	if( i .eq. 1) then
	  goto 1021
	else if( i .eq. 2) then
	  goto 2110
	  end if
	goto 102

1021	call find(k_res,laa,i)
	if(i.gt.0)  goto 102
	i=k_res
	k_res=k_res+1
	laa(i)=txt(ib:ie)
	iaa(i)=0
103	if(nword(n_len,txt,ib,ie)) goto 102
	if(txt(ib:ib).eq.':') then
	laa1(i)=txt(ib+1:ib+1)
	else if(txt(ib:ib+2).eq.'mw=') then
	  read(txt(ib+3:ie),*) maa(i)
	else 
	iaa(i)=iaa(i)+1
	kaa(iaa(i),i)=txt(ib:ie)
	end if
	goto 103

106	write(6,*) 
	1'dfres-W> ERROR in opening'//
	1 edp_data(:ltrim(edp_data))//pdbstd_dat
1061	close (8)
	txt=txt0
	n_len=n_len0
	ib=ib0
	ie=ie0
	return

2110	if(txt(ib:ie).ne.'dfork') goto 102
	call find(k_res,laa,l)
	if(l) 2111,2111,2112
2111	write(6,*) 'dfres-W> undefined residue '//txt(ib:ie)
	goto 102
2112	ifk(k_fk)=l
	iaal=iaa(l)
	if(iaal.le.0) goto 2111
	j_fk=0
2113	call find(iaal, kaa(1,l), ii)
	if(ii) 2114, 2115, 2116
2114	write(6,*)  
	1 'dfres-W> undefined atom '//txt(ib:ie)
	goto 102
2115	write(6,*)  
	1 'dfres-W> error or insufficient dfork information'
	goto 102
2116	j_fk=j_fk+1
	kfk(j_fk,k_fk)=ii
	if(j_fk .lt. 6) goto 2113

	call read_ar(2,t_rng(1,k_fk), *2115, *2115)
	cor(k_fk)=.not.nword(n_len,txt,ib,ie)
	if(cor(k_fk)) cor(k_fk)=txt(ib:ib).eq.'y'
	k_fk=k_fk+1
	goto 102

	entry dfres(*)
c	==============
	call find(k_res,laa,i)
	if(i.eq.0) goto 901
	if(i.eq.-1) then
	i=k_res
	if(i.gt.max_lres) then
	  write(6,*) 
	1 'dfres-W> UNDONE: too many residue types'
	  if(verbose.ge.3 ) write(6,*)
	1 '%EdPDB-I3- increase max_lres in edp_dim.inc (currently'
	1 ,max_lres,')'
	  return 1
	end if	
	k_res=k_res+1
	laa(i)=txt(ib:ie)
	end if
	iaa(i)=0
203	if(nword(n_len,txt,ib,ie)) return
	if(txt(ib:ib).eq.':') then
	laa1(i)=txt(ib+1:ib+1)
	else
	iaa(i)=iaa(i)+1
	if(iaa(i).gt.max_latom) then
	  write( errmsg(24:28),1062) max_latom
1062	  format(i5)
	  errmsg(:23)=' errmsg: max_atom # is'
	  return 1
	end if	
	kaa(iaa(i),i)=txt(ib:ie)
	end if
	goto 203

107	l2=k_res-1

	do 200 i=1,n_res
	  kres=res(i)
	  ijk1=ijk(i)
	  ijk2=ijk(i+1)-1

	  do l=1,l2
	    if(kres.eq.laa(l)) goto 108
	    end do
	  write(6,1104) 
	1 ' dfres-W> undefined residue-name '//kres//ires(i)
1104	  format(a)
	  goto 109

108	  iaal=iaa(l)
	  m_fk=0			! to prevent infinit loop in check_fk
	  if(iaal.gt.0) goto 110
109	  do j=ijk1,ijk2
	    iorder(j)=j
	    end do
	    goto 200

110	  do k=1,iaal
	    jorder(k)=0
	    end do

	  do 111 j=ijk1,ijk2
	    katom=atom(j)
	    do k=1,iaal
	      if(kaa(k,l).eq.katom) then
	        if(jorder(k).ne.0) then
	          write(6,'(a)') 
	1 ' dfres-W> double records '//
	1 katom//' '//res(i)//' '//ires(i)
		  goto 109
	          end if
	        jorder(k)=j
	        lf(j)=.true.
	        goto 111 
		end if
	      end do
	      write(6,'(a)') 
	1       ' dfres-W> undefined atom-name '//
	1       katom//' '//res(i)//' '//ires(i)
	      goto 109
111	    enddo

	  kk=0
	  do 112 k=1,iaal
	    if(jorder(k).eq.0) goto 112
	    kk=kk+1
	    iorder(ijk1+kk-1)=jorder(k)
112	    enddo
	  if(kk.eq.iaal) then
	    m_fk=m_fk+1
	    if(m_fk.gt.20) then
	      write(6,1104) 
	1      'dfres-W> unexpected ERROR: infinit loop in check_fk at '//
	1      res(i)//' '//ires(i)
	      goto 200
	      endif
	    if(check_fk()) goto 110
	  else if(verbose.ge.1 ) then
	    write(6,1104)
	1       ' dfres-W1> missing atom(s); no fork is checked: '//
	1      res(i)//' '//ires(i)
	    endif
	
200	  enddo
	return

	entry	sequ(file,*)
c	============
	n_of_syn=2					!000515
	syntax(1)='syntax:'
	syntax(2)=
	1'sequence [output_seq_filename.s] [fortran_format.s] [p0.i p1.i] '
 
	file1=file
	if(nword0(n_len,txt,ib,ie)) goto 300
	if(ib  .gt.ie) goto 300
	file1=txt(ib:ie)
300	if(.not.open_file(7,file1,'unknown','.seq')) return 1
	seq_out=file1

	form='(5(1x,10a1))'
	if(nword(n_len,txt,ib,ie)) goto 305
	ic=ichar(delimiter)
	if(ichar(txt(ib:ib)).eq. ic) then
	  do j=ib+1,n_len
	    if(ichar(txt(j:j)).eq.ic ) then
	      ib=ib+1
	      ie=j-1
	      txt(j:j)=' '	!000303, it may screw the history file
	      goto 303
	      endif
	    enddo
	  endif
303	form=txt(ib:ie)
	write(6,*) 'sequence-I> format=',form

	if(nword0(n_len,txt,ib,ie)) goto 305
	if(ib  .gt.ie) goto 305

	if(txt(ib:ib).eq.'c') then
	  i0=22
	  i1=22 
	else if(txt(ib:ib).eq.'r') then
	  i0=18
	  i1=20 
	else 
	  return 1
	  endif 
	ii=0
	do 304 i=1,n_res
	  do j=ijk(i),ijk(i+1)-1
	    if(lf(j)) goto 3031
	    enddo
	   goto 304
3031	   ii=ii+1
	   sequ1(ii)= text(j)(i0:i1)
304	   enddo
	goto 330	

305	call std0		! read edp_data_pdbstd_dat
	ii=0
	do 320 i=1,n_res
	do j=ijk(i),ijk(i+1)-1
	if(lf(j)) goto 310
	end do
	goto 320
310	ii=ii+1
	do j=1, k_res-1
	if(laa(j).eq.res(i)) goto 312
	end do
	sequ1(ii)= 'u'
	goto 320
312	sequ1(ii)= laa1(j)	
320	continue
330	write(7,form,err=390) (sequ1(j),j=1,ii)
1051	format(5(1x,10a1))
	close (7)
	return

390	return 1

        entry get_sequence(igr,sequ2)        !991117
chk     ==================
        call std0               ! read edp_data_pdbstd_dat

	jje=n_groupa(igr)
	do  jj=1,jje
          i_curr= aa_seq(igroupa( jj,igr))
         do jt=1, k_res-1
            if(laa(jt).eq.res(i_curr)) goto 612
            enddo
          sequ2(jj)= 'u'
          goto 620
612       sequ2(jj)= laa1(jt)  
620       enddo
        return

	entry seq2pdb_output(sequ2,ns,ns1, *)	!991228
chk	====================
	call std0
	na1=0
	do i=1,ns
	  si=sequ2(i)
	  jc=ichar(si)
	  if(tolower) then
	    jc=ichar(si)
	    if(65.le.jc.and.jc.le.90) si=char(jc+32)
	  endif

	  if(si.ne.' ' .and. ichar(si).ne.0) then
	    do jt=1, k_res-1
	      if(si.eq.laa1(jt)) goto 702
	      enddo
	      if(verbose.ge.2) write(6,*) ns1,' unrecoganizible residue symbol ['//si//']'  !121004
	    goto 701
702	    do j=1,iaa(jt)
	      na1=na1+1
	      write(4,1071) na1, kaa(j,jt), laa(jt), ns1,99.,99.,99.,0.,99.
	      enddo  
	    ns1=ns1+1
	    endif
701       enddo 
1071	format('ATOM',i7,2x,a4,a3,2x,i4,4x,3f8.3,2f6.2)
	return
!701	errmsg=' errmsg: unrecoganizible residue symbol ['//si//']'
	if(verbose.ge.4 ) then
          write(6,1047) 
1047    format(
	1 ' seq2pdb-I4> the input sequence is converted to lower case,'
	1/'   unless the tolower option is turned off')
	endif
	return 1

	entry	get_mw(mw)
c	==============
	n_of_syn=2					!000515
	syntax(1)='syntax:'
	syntax(2)='mw [molecular_weight_in_Da]'

	call std0		! read edp_data_pdbstd_dat
	ii=0
	mw=0.0
	do 420 i=1,n_res
	do j=ijk(i),ijk(i+1)-1
	if(lf(j)) goto 410
	end do
	goto 420
410	do j=1, k_res-1
	if(laa(j).eq.res(i)) goto 412
	end do
	goto 420
412	mw=mw+maa(j)
420	continue
430	write(6,1052,err=390) mw
1052	format(' get_mw> mw=',f10.0)
	return
	
	entry	match_mw(mw, delta)
c	================
	n_of_syn=2					!000515
	syntax(1)='syntax:'
	syntax(2)='mw [molecular_weight_in_Da, delta]'

	call std0		! read edp_data/pdbstd.txt
	ii=0
	do 520 i=1,n_res
	do j=ijk(i),ijk(i+1)-1
	if(lf(j)) goto 510
	end do
	goto 520
510	do j=1, k_res-1
	if(laa(j).eq.res(i)) goto 512
	end do
	goto 520
512	ii=ii+1
	amw(ii)=maa(j)
	imw(ii)=i
520	continue

	do i=1,ii-1
	bmw=0
	do j=i,ii
	  bmw=bmw+amw(j)
	  if(abs(bmw-mw).le.delta) 
	1   write(6,1520) ires(imw(i)),ires(imw(j)),bmw
	enddo
	enddo
	return
1520	format(' match_mw> ',5x,a5,'- ',a5,f10.1)

	entry shdfres
c	=============
901	if(k_res.lt.2) write(6,*) 'res-W> undefined'
	do i=1,k_res-1
	write(6,1002) laa(i), laa1(i), maa(i),(kaa(j,i), j=1,iaa(i))
1002	format(' res> ',a3,' :',a1,' mw=',f5.1,(t26,15a4))
	end do
	return

	entry check_fk0(L_fk)	! check the chirality and labeling
c	===============
	L_fk=.false.
	do j= 1,k_fk-1  ! k_fk-1 is the number of dfork definitions.
	if(ifk(j).eq.l) then
		j_a(1)=jorder(kfk(1,j))
		j_a(2)=jorder(kfk(2,j))
		j_a(3)=jorder(kfk(3,j))
		j_a(4)=jorder(kfk(4,j))
		j_a(5)=jorder(kfk(5,j))
		j_a(6)=jorder(kfk(6,j))
		curr=f_tor(j_a)
		if(curr.lt.t_rng(1,j)) curr=curr+360.
		if(curr.lt.t_rng(1,j).or.curr.gt.t_rng(2,j)) then
		 ic=2
		 if(cor(j)) then
			ic=1
			atomj=atom(j_a(5))		! flip the atom names
			atom(j_a(5))=atom(j_a(6))	!  in the internal 
			atom(j_a(6))=atomj		!  memory.
			ii=iorder(j_a(5))
			iorder(j_a(5))=iorder(j_a(6))
			iorder(j_a(6))=ii
			text(j_a(5))(14:17)=atom(j_a(5)) ! modify the displayed 
			text(j_a(6))(14:17)=atom(j_a(6)) !  atom name.
			L_fk=.true.
		   end if
		 write(6,*) 
	1 'dfres> ERROR fork '//text(j_a(1))(18:27),corr(ic)
		end if
	end if
	end do
	end

	subroutine switchwb()
chk	===================
	include 'edp_main.inc'
!	use edp_main

	n_of_syn=2					!000515
	syntax(1)='syntax:' 
	syntax(2)='switchwb' 

	do i=1,n_atom
	r=w(i)
	w(i)=b(i)
	b(i)=r
	end do
	end 

	logical function check_fk()
chk	=========================
	logical L_fk
	call check_fk0(L_fk)
	check_fk=L_fk
	end

	subroutine std0
chk	===============
	logical log
	data log/.true./
	if(log) call std(log)
	end 

	subroutine sort_in2(n,m,k)
chk	========================
	integer  m,mi
	dimension m(n)
	integer  k(n)
	do 4 i=1,n
4	k(i)=i
	do 70 i=2,n
	if(m(k(i)).ge.m(k(i-1)))goto 70
	mi=m(k(i))
	if(mi.lt.m(k(1)))then
	i2=1
	goto 20
	end if
	i1=1
	i2=i-1
10	if(i2-i1.le.1)goto 20
	i3=(i1+i2)/2
	if(mi.ge.m(k(i3)))then
	i1=i3
	else
	i2=i3
	end if
	goto 10
20	k1=k(i)
	do 30 l=i2,i
	k2=k(l)
	k(l)=k1
30	k1=k2
70	continue
	end

	subroutine sort_rl(n,m,k)
chk	==================
	real m,mi
	dimension m(n)
	integer  k(n)
	do 4 i=1,n
4	k(i)=i
	do 70 i=2,n
	if(m(k(i)).ge.m(k(i-1)))goto 70
	mi=m(k(i))
	if(mi.lt.m(k(1)))then
	i2=1
	goto 20
	end if
	i1=1
	i2=i-1
10	if(i2-i1.le.1)goto 20
	i3=(i1+i2)/2
	if(mi.ge.m(k(i3)))then
	i1=i3
	else
	i2=i3
	end if
	goto 10
20	k1=k(i)
	do 30 l=i2,i
	k2=k(l)
	k(l)=k1
30	k1=k2
70	continue
	end
chk***	the end of order.for

copyright by X. Cai Zhang

