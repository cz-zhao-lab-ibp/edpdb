	subroutine main_chain(igo,*,*)
c	===========================
	include 'edp_main.inc'
!	use edp_main
	character*1 string1
	data string1/'/'/

	character*4 atomi

	character*4 dfm(max_m),dca, symb_ca
	data (dfm(i),i=1,4),dca/'n','ca','c','o','ca'/	
	data num_m/4/

	logical  nword
	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie

	common /main_para/ tolower, jou, echo_L, inter0, incl0, incl1,
	1 jjou, nest, n_length			!9
	1 ,jnk_main_para(m_main_para-9)
	save

	ib0=ib
	ie0=ie
	if(.not.nword(n_len,txt,ib,ie)) then
	  if(i_pipe.ge.0) then
	    errmsg=' errmsg: illegal piping syntax' 
	    return 1
	    endif
	    
	  if(txt(ib:ie).ne.'from'.and.txt(ib:ib).ne.string1) then
	    if(n_len+12.gt.max_num_chars) then
	      errmsg=
	1' errmsg: input line is too long for CA,MAIN and SIDE commands.'
	      return 1
	      endif
!	    txt(ib:n_len+12)='from {zone '//txt(ib:n_len)//'}'
!	    n_len=n_len+12
	    txt(ib:n_len+7)='| zone '//txt(ib:n_len)
	    n_len=n_len+7
	    n_length=n_length+7
		backspace(2)
		backspace(48)
		return 2
	    endif
	  endif
	ib=ib0
	ie=ie0

c	igr:	id of the basckit group, 0-4
	call dfgroup(igr,*900)
	j1=n_groupa(igr)

600	goto (601,602,603) igo
601	n_of_syn=2					!000515
	syntax(1)='syntax:'
	syntax(2)='main'

	do j=1,j1
	  i=igroupa( j,igr)

	  if(lf(i).neqv.incl.and.text(i)(1:4).eq.'atom') then
	   atomi=atom(i)
	   do ii=1,num_m
	    if(atomi.eq.dfm(ii)) then
	     lf(i)=incl
	     goto 6011
	    end if
	   end do
	  end if
6011      continue
	  end do
	return

602	n_of_syn=2					!000515
	syntax(1)='syntax:'
	syntax(2)='side'

	do j=1,j1
	  i=igroupa( j,igr)
	  if(lf(i).neqv.incl) then
	    atomi=atom(i)
	    do ii=1,num_m
	      if(atomi.eq.dfm(ii)) goto 6021
	      enddo
	    lf(i)=incl
	    endif
6021       continue
	   enddo
	return

603	n_of_syn=2					!000515
	syntax(1)='syntax:'
	syntax(2)='ca'

	do j=1,j1
	  i=igroupa( j,igr)
	  if(atom(i).eq.dca) then
	    lf(i)=incl
	    endif
	  end do
	return

900	return 1

	entry dfmain(*)
c	===============
	n_of_syn=2					!000515
	syntax(1)='syntax:'
	syntax(2)='dfmain [atom_1.s [atom_2.s ... ]]'

	if(nword(n_len,txt,ib,ie)) goto 901
	num_m=1
	dfm(num_m)=txt(ib:ie)
	do while (.not.nword(n_len,txt,ib,ie))
	 if(num_m.ge.max_m) then
	  errmsg=' errmsg: the max. num. of atom-types in main is 10.'
	  return 1
	 end if
        num_m=num_m+1
	dfm(num_m)=txt(ib:ie)
	end do
	return

	entry shmain
c	============
901	write(6,1005) ' main>',( dfm(i), i=1,num_m)
1005	format(a,t10,10a5)
	return

	entry dfca(*)
c	=============
	n_of_syn=2					!000515
	syntax(1)='syntax:'
	syntax(2)='dfca [atom_name.s] ' 

	if(nword(n_len,txt,ib,ie)) goto 902
	dca = txt(ib:ie)
	if(.not.nword(n_len,txt,ib,ie)) return 1
	return

	entry shca
c	==========
902	write(6,1005) ' ca> CA is defined as ['//dca//'].'
	return

	entry get_ca(symb_ca)
	symb_ca=dca
	return

	entry avb(*,ks)
c	===========
	n_of_syn=4					!000515
	syntax(1)='syntax:' 
	syntax(2)='avb  [(x, y, z or w)]' 
	syntax(3)='rmsw [(x, y, z or b)]' 
	syntax(4)='sumw [(x, y, z or b)]' 

	if(nword(n_len, txt, ib,ie)) then
	  jst=-1
	else if(txt(ib:ie).eq.'x') then
		jst=31
	else if(txt(ib:ie).eq.'y') then
		jst=39
	else if(txt(ib:ie).eq.'z') then
		jst=47
	else if(txt(ib:ie).eq.'w') then
		if(ks.ne.0) then
		 errmsg=' errmsg: only {x,y,z,b} allowed'
		 return 1
		end if
		jst=0
	else if(txt(ib:ie).eq.'b') then
		if(ks.eq.0) then
		 errmsg=' errmsg: only {x,y,z,w} allowed'
		 return 1
		end if
		jst=0
	else 
		return 1
	end if

	if(ks.eq.0) then 		!run avb
		total=0.0
		itt=0
		do 5001 i=1,n_res
			k=0
			kk=0
			av=0.
			do j=ijk(i), ijk(i+1)-1
				if(atom(j).eq.dca) k=j
				if(lf(j)) then
					av = av+b(j)
					kk=kk+1
				end if
			end do
			total=total+av
			itt=itt+kk
			if(kk.eq.0 .or. jst .eq. -1) goto 5001
			if(k.eq.0) then
	write(6,*)'avb-W> ABORTED: undefined ca in residue ['//res(i)//']'
				return
			end if
			if(jst.gt.1) then
			  write(text(k)(jst:jst+7),1061) av/kk
1061	format(f8.3)
			else
				w(k)=av/kk
			end if
5001		continue
		write(6,1003)  total/max(itt,1)
1003		format(' avb>     average of b= ',e12.5)

	else if(ks.eq.1) then	! run sumw
		total=0.0
		do 5002 i=1,n_res
			k=0
			kk=0
			av=0.
			do j=ijk(i), ijk(i+1)-1
				if(atom(j).eq.dca) k=j
				if(lf(j)) then
					av = av+w(j)
					kk=kk+1
				end if
			end do
			total=total + av
			if(kk.eq.0 .or. jst .eq. -1) goto 5002
			if(k.eq.0) then
	write(6,*) 'sumw-W> ABORTED: undefined ca in residue ['//res(i)//']'
 				return
			end if
			if(jst.gt.1) then
			  write(text(k)(jst:jst+7),1061) av
			else
				b(k)=av
			end if
5002		continue
		write(6,1001)  total
1001		format(' sumw> total summation= ',e12.5)
	else		! run rmsw
		total=0.0
		itt=0
		do 5003 i=1,n_res
			k=0
			kk=0
			av=0.
			do j=ijk(i), ijk(i+1)-1
				if(atom(j).eq.dca) k=j
				if(lf(j)) then
					av = av+ w(j)*w(j)
					kk=kk+1
				end if
			end do
			total=total+ av
			itt=itt+kk
			if(kk.eq.0 .or. jst .eq. -1) goto 5003
			if(k.eq.0) then
	write(6,*) 'rmsw-W> ABORTED: undefined ca in residue ['//res(i)//']'
 				return
			end if
			if(jst.gt.1) then
			  write(text(k)(jst:jst+7),1061) sqrt(av/kk)
			else
				b(k)=sqrt(av/kk)
			end if
5003		continue
		write(6,1007)  sqrt(total/max(itt,1))
1007		format(' rmsw>        rms of w= ',e12.5)
	end if
	end

copyright by X. Cai Zhang

	subroutine search_symm(a)
chk	======================
	include 'edp_main.inc'
!	use edp_main

	logical match_symm
	external match_symm

	common /cmm_symm/ num_symm, symmtx(3,4,max_symm),
	1 trn2(3,4,max_symm)
	real a(3,4), ai(3,4)

	do i=1,3
	do j=1,3
	  ai(i,j)=a(i,j)
	  enddo
	  enddo

	call mxinv(3, ai,ierr)

	ai(1,4)=-ai(1,1)*a(1,4)-ai(1,2)*a(2,4)-ai(1,3)*a(3,4)
	ai(2,4)=-ai(2,1)*a(1,4)-ai(2,2)*a(2,4)-ai(2,3)*a(3,4)
	ai(3,4)=-ai(3,1)*a(1,4)-ai(3,2)*a(2,4)-ai(3,3)*a(3,4)

	do i0=1,num_symm
	  if(match_symm(a,trn2(1,1,i0),i1,i2,i3) ) then
	    write(6,1002) i0,i1,i2,i3
	    goto 100
	    endif
	  enddo
	write(6,1003)

100	if(ierr.eq.0) then
	 do i0=1,num_symm
	  if(match_symm(ai,trn2(1,1,i0),i1,i2,i3) ) then
	    write(6,1004) i0,i1,i2,i3
	    return
	    endif
	  enddo
	 endif
	write(6,1005)
1002	format(' the equivalent               symm ',4i4)
1003	format(' the equivalent               symm undetermined')
1004	format(' the inverse                  symm ',4i4)
1005	format(' the inverse                  symm undetermined')
	end

	logical function match_symm(a,trn,i1,i2,i3) 
chk	===========================
	real a(3,4),trn(3,4), dt(3), da(3)
	common /cmm_cell/ cell(6), irl, trn0(3,3), trn1(3,3)
	data eps/1.e-4/

	match_symm=.false.

	do i=1,3
	do j=1,3
	  if(abs(a(i,j)-trn(i,j)).gt. eps) return
	  enddo
	  dt(i) = a(i,4)-trn(i,4)
	  enddo

	do i=1,3
	  da(i)=trn0(i,1)*dt(1)+trn0(i,2)*dt(2)+trn0(i,3)*dt(3)
	  if(abs(da(i)-float(nint(da(i)))).gt.eps) return
	  enddo
	i1=nint(da(1))
	i2=nint(da(2))
	i3=nint(da(3))
	match_symm=.true.
	end

	subroutine axisdist (*)
chk	===================
	include 'edp_main.inc'
!	use edp_main
	include 'edp_file.inc'
!	use edp_file

	logical  open_file1
	external open_file1
	character*2 s_xyz(3)
	data s_xyz/'x','y','z'/

	dimension a(3,4), a_dummy(3,4), f(3), euler(3)
	logical l_axis
	save

	n_of_syn=2					!000515
	syntax(1)='syntax:'
	syntax(2)='axis file_name.s [vector_id.s [(x,y,z)]]' 

	l_axis=.true.
	rtn_in='rtn_.txt'
        if(.not. open_file1(29, rtn_in, 'old','.txt')) return 1
    	read(29,*,end=900,err=900) ((a(i,j),j=1,3),i=1,3), (a(i,4),i=1,3)
        close (29)

	call get_vector_id(iv ,*10)
	i_axis=match_l(3,s_xyz)
	goto 10

	entry rtnout(a_dummy)
chk	============
	do i=1,3
	do j=1,4
	  a(i,j)=a_dummy(i,j)
	  enddo
	  enddo
	l_axis=.false.
	i_axis=0

10	call search_symm(a)
	call polar( a,ierr)
	if(ierr.gt.0) return
	call get_euler(a, euler)
	write(6,1001) euler

	call vector(1., a, f, *50)
	phi= (a(1,1) +a(2,2) +a(3,3))*.5 -.5
	if(abs(phi)-1..gt.1.e-5) then
	  write(6,1116) 
	  return
	else if(phi.ge.0.) then
	  phi=acosd(min(phi,1.)) 
	else 
	  phi=acosd(max(phi,-1.)) 
	  endif
	write(6,1106)  f,sign(phi,(a(1,3)-a(3,1))*f(2))

	call axis03 (a, f,tt1,tt2,tt3, screw_length, i_axis, *60)
	write(6,1119)  tt1,tt2,tt3,screw_length,a(1,4),a(2,4),a(3,4)
	1 , a(1,1)*a(1,4)+a(2,1)*a(2,4)+a(3,1)*a(3,4)
	1 , a(1,2)*a(1,4)+a(2,2)*a(2,4)+a(3,2)*a(3,4)
	1 , a(1,3)*a(1,4)+a(2,3)*a(2,4)+a(3,3)*a(3,4)
	if(l_axis) goto 100
	return

50	write(6,1116)
1116	format(
	1 ' the direction of rotation axis is undetermined'
	1/' the rotation angle             is undetermined'
	1/' the screw axis passes through   ( undetermined )'
	1/' the screw length               is undetermined')
	return
60	write(6,1117)
1117	format(
	1/' the screw axis passes through   ( undetermined )'
	1/' the screw length               is undetermined')
	return 

100	do i=1,n_atom
	  if(lf(i)) then
	    rx=x(i)-tt1
	    ry=y(i)-tt2
	    rz=z(i)-tt3
	    r=rx*f(1)+ry*f(2)+rz*f(3)
	    w(i)=sqrt(max(0.,rx*rx+ry*ry+rz*rz-r*r))
	    end if
	  end do

	if(iv.gt.0) then
	  vectors(1,iv)= tt1
	  vectors(2,iv)= tt2
	  vectors(3,iv)= tt3
	  vectors(4,iv)= f(1)
	  vectors(5,iv)= f(2)
	  vectors(6,iv)= f(3)
	  vectors(7,iv)= 1.0
	  endif

	write(6,*) 
	1'axis-I> the occs of on atoms have been changed to distance.'
1106	format(  ' the direction of rotation axis is ',3f10.3
	1       /' the rotation angle             is ',f10.3)
1119	format(  ' the screw axis passes through   ( ',3e10.2,')'
	1       /' the screw length               is ',f10.3,
	1       /' the translation (aft. rot.)   are ',3f10.3,
	1       /' the translation (bef. rot.)   are ',3f10.3)
1001	format(  ' the eulerian angle (z,y*,z**) are ',3f10.3) 
900	end

	subroutine get_euler(a, euler)
c	====================
	dimension euler(3), a(3,3)
	if(abs(a(3,3)-1.).lt.1.e-5) then
	  euler(1)=0.
	  euler(2)=0.
	  euler(3)=atan2d(a(2,1),a(1,1))
	else if(abs(a(3,3)+1.).lt.1.e-5) then
	  euler(1)=0.
	  euler(2)=180.
	  euler(3)=atan2d(a(1,2),a(2,2))
	else
	  if(abs(a(3,3)) .gt. 1.0 ) then
	    write(6,*) '%EdPDB-W- wrong matrix, m(3,3)=',a(3,3)
	    a(3,3)=sign(1.0,a(3,3))
	    endif
	  euler(1)=atan2d(a(2,3),a(1,3))
	  euler(2)=acosd(a(3,3)) 
	  euler(3)=atan2d(a(3,2),-a(3,1))
	  endif
	end  	   

	subroutine atoms_to_vector(*)
chk	=========================
	include 'edp_main.inc'
!	use edp_main

	integer atom_id(2)
	equivalence (i1, atom_id(1)), (i2, atom_id(2))

	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie

	call get_vector_id(iv,*900)

	call get_atoms(2, atom_id, ierr)
	if(ierr.ne.0) return 1

	vectors(1,iv)=x(i1)
	vectors(2,iv)=y(i1)
	vectors(3,iv)=z(i1)

	xc=x(i2)-x(i1)
	yc=y(i2)-y(i1)
	zc=z(i2)-z(i1)
	r=sqrt(xc**2+ yc**2+ zc**2)
	if(r.lt.eps()) then 
	  vectors(4,iv)=0.0
	  vectors(5,iv)=0.0
	  vectors(6,iv)=1.0
	  vectors(7,iv)=0.0
	else
	  vectors(4,iv)=xc/r
	  vectors(5,iv)=yc/r
	  vectors(6,iv)=zc/r
	  vectors(7,iv)=r
	  endif
	return 
900	return 1
	end

	subroutine get_vector_id(iv, *)
chk	========================
chk	create or overwrite a vector
chk	no name: 	return 1

	include 'edp_main.inc'
!	use edp_main

	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie

	iv=0
	jv=match_l(max_vector,cvector)

	if(ie-ib+1.gt.4) then
	  write(6,*) 
	1 '%EdPDB-W- only will the first four charactors be used as'
	1 ,' vector name.'
	  endif

	if(jv.eq.0) return 1
	if(jv.lt.0) then
	  do iv=1,max_vector
	    if(cvector(iv).eq.' ') then
	      cvector(iv)=txt(ib:ie)
	      return
	      endif
	    enddo
	  errmsg=' errmsg: too many vectors defined.'
	  return 1
	else
	  iv=jv
	  endif
	end

	subroutine find_vector_id(iv, *,*)
chk	========================
chk	find an existing vector
chk	wrong name: 	return 1
chk	no name:	return 2

	include 'edp_main.inc'
!	use edp_main

	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie

	iv=0
	jv=match_l(max_vector,cvector)

	if(ie-ib+1.gt.4) then
	  write(6,*) 
	1 '%EdPDB-W- only will the first four charactors be used as'
	1 ,' vector name.'
	  endif

	if(jv.lt.0) return 1
	if(jv.eq.0) return 2
	iv=jv
	end

	subroutine list_vector(*)
chk	========================
chk	wrong name: 	return 1
chk	no name:	return 2

	include 'edp_main.inc'
!	use edp_main

	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie

	real t(7)

	parameter (n_lab=8)
!	integer, parameter :: n_lab=8
	character*8 label(n_lab)
        data label/'list','by_xyz','vv','pv','vp'	!5
	1 ,'delete','av_cos','by_atom'/			!8

	n_of_syn=7					!000515
	syntax(1)='syntax:'
	syntax(2)='1) vector by_atom res_1.s atom_1.s res_2.s atom_2.s' 
	syntax(3)=
	1'2) vector by_xyz vector_id.s p1.r p2.r p3.r r1.r r2.r r3.r '//
     1'[length.r]' 
	syntax(4)='3) vector delete vector_id.s '
	syntax(5)='4) vector list [vector_id.s] '
	syntax(6)='5) vector pv vector_id.s [res_id.s [atom_name.s]]'
	syntax(7)='6) vector vp vector_id.s [res_id.s [atom_name.s]] '//
	1'[length.r]' 

	call find1(n_lab,label,iq)
        if(iq.lt.0) return 1

	if(iq.le.1) then
	  iv1=0
	  if(iq.eq.1) call find_vector_id(iv1,*900,*600)

600	  write(6,1002)
	  do i=1,max_vector
	    if(iv1.le.0 .or. iv1.eq.i) then
	      if(cvector(i).ne.' ') then
	        write(6,1001) cvector(i), (vectors(j,i),j=1,7),
	1       vectors(1,i)+vectors(4,i)*vectors(7,i),
	1       vectors(2,i)+vectors(5,i)*vectors(7,i),
	1       vectors(3,i)+vectors(6,i)*vectors(7,i)
	        endif
	      endif
	    enddo
1001	  format(1x,a4,3f7.2,',',3f7.3,',',f7.1,',',3f7.2)
1002	  format(' vector>     origin',16x,'vector',9x,'length',12x,'end')
	else if(iq.eq.2) then
	  call get_vector_id(iv,*900)
	  do j=1,6
	    t(j)=vectors(j,iv)
	    enddo
	  t(7)=1.
	  call read_ar(7,t,*900,*800)

800	  r=sqrt(t(4)**2+ t(5)**2+ t(6)**2)
	  if(abs(r-1.0).lt.4.0e-3.and.abs(t(7)).gt.eps()) then
	    do j=1,7
	      vectors(j,iv)=t(j)
	      enddo
	    if(t(7).lt.0.0) then
	      do j=4,7
	        vectors(j,iv)=-vectors(j,iv)
	        enddo
	      endif
	  else 
	    if(abs(t(7)).gt.eps().and.
	1      abs(t(7)-vectors(7,iv)).gt.eps()) then
	      write(6,*) 
	1 '%EdPDB-W- inconsistence was found in the input parameters.'
      write(6,*)
	1 '          the seventh parameter is discarded.'
	      endif

	    do j=1,3
	      vectors(j,iv)=t(j)
	      t(j+3)=t(j+3)-t(j)
	      enddo
	    r=sqrt(t(4)**2+ t(5)**2+ t(6)**2)
	    if(abs(r).lt.eps()) then
	      vectors(4,iv)=0.0
	      vectors(5,iv)=0.0
	      vectors(6,iv)=1.0
	      vectors(7,iv)=0.0
	    else
	      vectors(4,iv)=t(4)/r
	      vectors(5,iv)=t(5)/r
	      vectors(6,iv)=t(6)/r
	      vectors(7,iv)=r
	      endif
	    endif   	  
	else if(iq.eq.3) then
	  call vector_vector(*900)
	else if(iq.eq.4) then
	  call point_vector(*900)
	else if(iq.eq.5) then
	  call vector_point(*900)
	else if(iq.eq.6) then
	  call find_vector_id(iv,*900,*900)
	  cvector(iv)=' ' 
	else if(iq.eq.7) then
	  call av_cos(*900)
	else if(iq.eq.8) then
	  call atoms_to_vector(*900)
	else !if(i.eq.3) then
	  return 1
	  endif
	return
900	return 1
	end

c----	end of axisdist.for

copyright by X. Cai Zhang
