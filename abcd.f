
copyright by X. Cai Zhang
	subroutine abc(*)
chk	==============
	include 'edp_main.inc'
!	use edp_main

	character*4 atoma(3), store*1
	data atoma/3*' '/
	dimension ka(3), j_a(3)
	data ka/3*0/
	logical la(3), nword, nword0
	data la, amin, amax/3*.true., 0., 180./
c-
	character*(max_num_chars) txt
	common /cmm_txt/ n_len,txt,ib,ie
	save

	n_of_syn=2					!000515
	syntax(1)='syntax:'
	syntax(2)='abc [(a, b, c) [(x, y, z)]]'

	do i=1,3
	 if(atoma(i).eq.' ') then
	  write(6,*) 'abc-W> define atom1 -  atom3 first, using dfabc.'
	  return
	 end if
	end do

	if(nword(n_len,txt,ib,ie)) then
	ist=0
	goto 10
	else
	store=txt(ib:ie)
	end if

	if(store.eq.'a') then
	ist=1
	else if(store.eq.'b') then
	ist=2
	else if(store.eq.'c') then
	ist=3
	else
	return 1
	end if

	if(nword(n_len,txt,ib,ie)) then
	jst=0
	goto 10
	else
	store=txt(ib:ie)
	end if

	if(store.eq.'x') then
	jst=31
	else if(store.eq.'y') then
	jst=39
	else if(store.eq.'z') then
	jst=47
	else
	return 1
	end if
c-
c	call openlist(48)
c-
10	ksum=0
	av=0.
	sgm=0.
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
	curr=f_angl(j_a)
	if(curr.gt.amax.or.curr.lt.amin) goto 800
	av=av+curr
	sgm=sgm+curr*curr
	write(48,1001) (text(j_a(ii))(14:27), ii=1,3),
	1 curr
	ksum=ksum+1
	if(ist.gt.0) then
	  igroupa( ksum,1)=j_a(ist)
	  if(jst.gt.0) write(text(j_a(ist))(jst:jst+7),1061) curr
1061	  format(f8.3)
	  end if

1001	format(1x,3(a14,','),'a= ',f8.3)
800	continue
	if(ist.gt.0) n_groupa(1)=ksum
	if(ksum.eq.0)then
	write(6,1009) ksum
	else
	sumk=1./ksum
	av=av*sumk
	sgm=sgm*sumk
	write(6,1009) ksum, av, sqrt(sgm-av*av)
	write(48,1009)ksum, av, sqrt(sgm-av*av)
	end if
1009	format('^abc> ',i6,' records are calculated':
	1', ava=',f8.3,', sgma=',f8.3)
	call typelist(48)
	return

	entry dfabc(*)
chk	=============
	n_of_syn=4					!000515
	syntax(1)='syntax:'
	syntax(2)='1) dfabc '
	syntax(3)='2) dfabc atom_a.s atom_b.s atom_c.s '
	syntax(4)='    [reg_a reg_b reg_c] '//
	1'[status_a status_b status_c] [amin amax]' 
c
c	define atoma
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
	amin=0.
	amax=180.
c
c	ka
c
	do ii=1,3
	  if(nword0(n_len,txt,ib,ie)) return
	  if(ib .le. ie) read(txt(ib:ie),*,err=900) ka(ii)
	  end do
c
c	la
c
	do ii=1,3
	  if(nword0(n_len,txt,ib,ie)) return
	  if(ib .le. ie) read(txt(ib:ib),1063, err=900) la(ii)
1063	  format(l1)
	  enddo
c
c	amin, amax
c
	call read_ar(1, amin, *900, *201)
201	call read_ar(1, amax, *900, *202)
202	return
900	return 1
c
	entry shabc
chk	===========
901	write(6,1002) atoma, ka, la, amin, amax
1002	format(' abc>',t10,3a5,t30,3i5,t55,3(l1,1x),t65,2f7.1)
	end

	function f_angl(j_a)
chk	===============
	include 'edp_main.inc'
!	use edp_main
	dimension j_a(3),ax(3),ay(3),az(3)

	f_angl=999.
	do ii=1,3
	 ax(ii)=x(j_a(ii))
	 ay(ii)=y(j_a(ii))
	 az(ii)=z(j_a(ii))
	 end do

	x1=ax(1)-ax(2)
	y1=ay(1)-ay(2)
	z1=az(1)-az(2)
	x2=ax(3)-ax(2)
	y2=ay(3)-ay(2)
	z2=az(3)-az(2)
	r=sqrt((x1*x1+y1*y1+z1*z1)*(x2*x2+y2*y2+z2*z2))
	if(r.lt. eps()*10.0) return
	f_angl= acosd((x1*x2+y1*y2+z1*z2)/r)
	end
c
	function f_tor(j_a)
chk	==============
	include 'edp_main.inc'
!	use edp_main
	dimension j_a(4),ax(4),ay(4),az(4)
c
	do ii=1,4
	 ax(ii)=x(j_a(ii))
	 ay(ii)=y(j_a(ii))
	 az(ii)=z(j_a(ii))
	end do
10	a1=ax(2)-ax(1)
	b1=ax(3)-ax(2)
	c1=ax(4)-ax(3)
	a2=ay(2)-ay(1)
	b2=ay(3)-ay(2)
	c2=ay(4)-ay(3)
	a3=az(2)-az(1)
	b3=az(3)-az(2)
	c3=az(4)-az(3)
	axb1=a2*b3-a3*b2
	axb2=a3*b1-a1*b3
	axb3=a1*b2-a2*b1
	bxc1=b2*c3-b3*c2
	bxc2=b3*c1-b1*c3
	bxc3=b1*c2-b2*c1
	ab=axb1*axb1+axb2*axb2+axb3*axb3
	if(ab.lt.1.e-6) goto 900
	bc=bxc1*bxc1+bxc2*bxc2+bxc3*bxc3
	if(bc.lt.1.e-6) goto 900
	vol=b1*(axb2*bxc3-axb3*bxc2)+	
	1   b2*(axb3*bxc1-axb1*bxc3)+	
	1   b3*(axb1*bxc2-axb2*bxc1)
	f_tor= (axb1*bxc1+axb2*bxc2+axb3*bxc3)/sqrt(ab*bc)
	if(abs(f_tor).ge.1.) then
	ax(1)=ax(1)+0.001			! force the result away from single point
	goto 10
	end if
	f_tor= sign(acosd(f_tor),vol)
	return
900	f_tor=999.
	end

copyright by X. Cai Zhang

	subroutine abcd(*)
chk	===============
	include 'edp_main.inc'
!	use edp_main

	character*4 atoma(4), store*1
	data atoma/4*' '/
	dimension ka(4), j_a(4)
	data ka/4*0/
	logical la(4), nword, nword0
	data la, amin, amax/4*.true., -180.,180./

	character*(max_num_chars) txt
	common /cmm_txt/ n_len,txt,ib,ie
	save

	n_of_syn=2					!000515
	syntax(1)='syntax:'	
	syntax(2)='abcd [(a, b, c, d}, [(x, y, z)]]' 

	do i=1,4
	 if(atoma(i).eq.' ') then
	  write(6,*) 'abcd-W> define atom1 -  atom4 first, using dfabcd.'
	  return
	 end if
	end do

	if(nword0(n_len,txt,ib,ie)) then
	ist=0
	goto 10
	else
	store=txt(ib:ie)
	end if

	if(store.eq.'a') then
	ist=1
	else if(store.eq.'b') then
	ist=2
	else if(store.eq.'c') then
	ist=3
	else if(store.eq.'d') then
	ist=4
	else
	return 1
	end if

	if(nword0(n_len,txt,ib,ie)) then
	jst=0
	goto 10
	else
	store=txt(ib:ie)
	end if

	if(store.eq.'x') then
	jst=31
	else if(store.eq.'y') then
	jst=39
	else if(store.eq.'z') then
	jst=47
	else
	return 1
	end if

c-
c	call openlist(48)
c-
10	ksum=0
	av=0.
	sgm=0.
	ii0=min(ka(1),ka(2),ka(3),ka(4))
	ii1=max(ka(1),ka(2),ka(3),ka(4))	
	i0=1-ii0
	i1=n_res-ii1
c
	do 800 i=i0,i1
	 do ii=1,4
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
c-	calculate tor
c
	curr=f_tor(j_a)
	if(curr.lt.amin) curr=curr+360.
	if(curr.gt.amax.or.curr.lt.amin) goto 800

	ksum=ksum+1
	av=av+curr
	sgm=sgm+curr*curr
	write(48,1001) (text(j_a(ii))(14:27), ii=1,4),
	1 curr
1001	format(1x,4(a14,','),'t= ',f8.3)
	if(ist.gt.0) then
	  igroupa( ksum,1)=j_a(ist)
	  if(jst.gt.0) write(text(j_a(ist))(jst:jst+7),1061) curr
1061	  format(f8.3)
	  endif

800	continue
	if(ist.gt.0) n_groupa(1)=ksum
	if(ksum.eq.0)then
	write(6,1009) ksum
	else
	sumk=1./ksum
	av=av*sumk
	sgm= sqrt( max(sgm*sumk-av*av,0.))
	write(6,1009) ksum, av, sgm
	write(48,1009)ksum, av, sgm
	end if
1009	format('^abcd>',i6,' records are calculated':
	1', avt=',f8.3,', sgmt=',f8.3)
	call typelist(48)
	return

	entry dfabcd(*)
c	=============
	n_of_syn=5					!000515
	syntax(1)='syntax:'
	syntax(2)='1) dfabcd '
	syntax(3)='2) dfabcd atom_a.s atom_b.s atom_c.s atom_d.s'
	syntax(4)=
	1'    [reg_a.i reg_b.i reg_c.i reg_d.i]'
	syntax(5)=
	1'     [status_a.l status_b.l status_c.l status_d.l] [tmin.r tmax.r]' 

c
c	atoma
c
	do ii=1,4
	if(nword(n_len,txt,ib,ie)) then
	if(ii.gt.1) return 1
	goto 901
	end if
	atoma(ii)=txt(ib:ie)
	ka(ii)=0
	la(ii)=.true.
	end do
	amin=-180.
	amax=180.
c
c	ka
c
	do ii=1,4
	  if(nword0(n_len,txt,ib,ie)) return
	  if(ib .le. ie) read(txt(ib:ie),*,err=900) ka(ii)
c-u1062	  format(i<ie-ib+1>)
	  enddo
c
c	la
c
	do ii=1,4
	  if(nword0(n_len,txt,ib,ie)) return
	  if(ib .le. ie) read(txt(ib:ib),*,err=900) la(ii)
	  end do
c
c	amin, amax
c
        call read_ar(1, amin, *900, *201)
201     call read_ar(1, amax, *900, *202)
202     return
900     return 1

	entry shabcd
chk	============
901	write(6,1002) atoma, ka, la, amin, amax
1002	format(' abcd>',t10,4a5,t30,4i5,t55,4(l1,1x),t65,2f7.1)
	end

	subroutine ab(*)
chk	=============
	include 'edp_main.inc'
!	use edp_main

	character*4 atoma(2), store*1
	data atoma/2*' '/
	dimension ka(2), j_a(2)
	data ka/2*0/
	logical la(2) 
	data la, amin, amax/2*.true., 0., 99./

	logical nword, nword0
	character*(max_num_chars) txt
	common /cmm_txt/ n_len,txt,ib,ie
	save

	n_of_syn=2					!000515
	syntax(1)='syntax:'
	syntax(2)='ab [(a, b) [(x, y, z)]]'

	do i=1,2
	 if(atoma(i).eq.' ') then
	  write(6,*) 'ab-W> define atom1 -  atom2 first, using dfab.'
	  return
	 end if
	end do

	if(nword(n_len,txt,ib,ie)) then
	ist=0
	goto 10
	else
	store=txt(ib:ie)
	end if

	if(store.eq.'a') then
	ist=1
	else if(store.eq.'b') then
	ist=2
	else
	return 1
	end if

	if(nword(n_len,txt,ib,ie)) then
	jst=0
	goto 10
	else
	store=txt(ib:ie)
	end if

	if(store.eq.'x') then
	jst=31
	else if(store.eq.'y') then
	jst=39
	else if(store.eq.'z') then
	jst=47
	else
	return 1
	end if
c-
c	call openlist(48)
c-
10	ksum=0
	av=0.
	sgm=0.
	ii0=min(ka(1),ka(2))
	ii1=max(ka(1),ka(2))	
	i0=1-ii0
	i1=n_res-ii1

	do 800 i=i0,i1
	 do ii=1,2
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

c-	calculate distance

	curr=f_dist(j_a)
	if(curr.gt.amax.or.curr.lt.amin) goto 800
	ksum=ksum+1
	av=av+curr
	sgm=sgm+curr*curr
	write(48,1001) (text(j_a(ii))(14:27), ii=1,2),
	1 curr
1001	format(1x,2(a14,','),'d= ',f8.3)
	if(ist.gt.0) then
	  igroupa( ksum,1)=j_a(ist)
	  if(jst.gt.0) write(text(j_a(ist))(jst:jst+7),1064) curr
1064	  format(f8.3)
	  end if

800	continue
	if(ist.gt.0) n_groupa(1)=ksum
	if(ksum.eq.0)then
	write(6,1009) ksum
	else
	sumk=1./ksum
	av=av*sumk
	sgm= sqrt( max(sgm*sumk-av*av,0.))
	write(6,1009) ksum, av, sgm
	write(48,1009)ksum, av, sgm
	end if
1009	format('^ab>  ',i6,' records are calculated':
	1', avd=',f8.3,', sgmd=',f8.3)
	call typelist(48)
	return

	entry dfab(*)
c	=============
	n_of_syn=3					!000515
	syntax(1)='syntax:'
	syntax(2)=
	1'dfab [atom_a.s atom_b.s [reg_a.i reg_b.i]'
	syntax(3)=
	1'  [status_a.l status_b.l] [dmin.r dmax.r]]' 

c
c	atoma
c
	do ii=1,2
	if(nword(n_len,txt,ib,ie)) then
	if(ii.gt.1) return 1
	goto 901
	end if
	atoma(ii)=txt(ib:ie)
	ka(ii)=0
	la(ii)=.true.
	end do
	amin=0.
	amax=99.
c
c	ka
c
	do ii=1,2
	  if(nword0(n_len,txt,ib,ie)) return
	  if(ib .le. ie) read(txt(ib:ie),*, err=900) ka(ii)
c-u1065	  format(i<ie-ib+1>)
	  end do
c
c	la
c
	do ii=1,2
	  if(nword0(n_len,txt,ib,ie)) return
	  if(ib .le. ie) read(txt(ib:ib),1066,err=900) la(ii)
1066	  format(l1)
	  end do
c
c	amin, amax
c
        call read_ar(1, amin, *900, *201)
201     call read_ar(1, amax, *900, *202)
202     return
900     return 1

	entry shab
chk	==========
901	write(6,1002) atoma, ka, la, amin, amax
1002	format(' ab>',t10,2a5,t30,2i5,t55,2(l1,1x),t65,2f7.1)
	end

copyright by X. Cai Zhang

	subroutine read_acctxt(def_rdi,*)
chk	======================
	include 'edp_main.inc'
!	use edp_main
	include 'edp_dat.inc'
!	use edp_dat
	
	character*4 rt, r_atom, aa0

	common /accrt/ num_rt,rt(max_rt), irt(max_rt), 
	1 num_r_a, r_atom(max_r_a), rd(max_r_a)
	
	logical nword 
	character*(max_num_chars) txt, txt0
	common /cmm_txt/ n_len,txt,ib,ie
	
	txt0=txt 
	n_len0=n_len
	ib0=ib
	ie0=ie

100	call read_txt(8,*999)

	ie=0
	if(nword(n_len,txt,ib,ie))  goto 100
	if(txt(ib:ie).eq.'dacc') then

	if(nword(n_len,txt,ib,ie)) goto 200
	if(txt(ib:ie).ne.aa0) then	! new residue type
	if(num_rt.ge.max_rt) then
	  write(6,*) 'acc-W> too many res types.'
	  if(verbose .ge. 3 ) write(6,*)
	1 '%EdPDB-I3- increase max_rt in edp_dim.inc (currently'
	1 ,max_rt,')'
	  goto 991
	  end if
	rt(num_rt)= txt(ib:ie)
	irt(num_rt)= num_r_a
	aa0=rt(num_rt)
	num_rt=num_rt+1
	end if

	if(nword(n_len,txt,ib,ie)) goto 200
	if(num_r_a.gt.max_r_a) then
	  write(6,*) 'acc-W> too many atom types.'
	  if(verbose .ge. 3 ) write(6,*)
	1 '%EdPDB-I3- increase max_r_a in edp_dim.inc (currently'
	1 ,max_r_a,')'
 	goto 991
	end if
	r_atom(num_r_a)= txt(ib:ie)
	call read_ar( 1, rd(num_r_a), *200, *200)
	num_r_a=num_r_a+1
	end if
	goto 100

200	write(6,1002) txt	! error message	
1002	format(1x,a78)
	txt(:6)=' error'
	txt(7:ib-4)=' '
	txt(ib-3:)='^'
	write(6,1002) txt	
991	close(8)
	txt=txt0
	n_len=n_len0
	ib=ib0
	ie=ie0
	return 1

999	close(8)
chk	def_rdi is defined by 'dacc * * 1.8' 
	do k=1,num_rt
	  if(rt(k)(1:1) .eq. wildcard) then
	    do kk=irt(k),irt(k+1)-1
	      if(r_atom(kk)(1:1).eq.wildcard) then
	        def_rdi=rd(kk)
		goto 998
	      endif
	    end do
	  endif 
	end do
	if(verbose .ge. 3) write(6,*) 
	1 ' read_acc-I3> default radius=',rd(kk)
998	txt=txt0
	n_len=n_len0
	ib=ib0
	ie=ie0

	irt(num_rt)= num_r_a
	num_rt=num_rt-1
	num_r_a=num_r_a-1
	end


	subroutine acc_index(n,index1,index2,rdi,r_max,*)
chk	====================
	include 'edp_main.inc'
!	use edp_main
	include 'edp_file.inc'
!	use edp_file
	include 'edp_dat.inc'
!	use edp_dat

	logical  open_file
	external open_file
	integer index1(max_atom),index2(max_atom)
	real	rdi(max_atom)
	character*4  rt, r_atom, resj0
!	character*4 resj,atomi, rt, r_atom, resj0

	common /accrt/ num_rt,rt(max_rt), irt(max_rt), 
	1 num_r_a, r_atom(max_r_a), rd(max_r_a)

!	logical L_rt

	character*(108) file_used
	data file_used/ ' ' /

c	character*4 aa0

!	logical nword 
	character*(max_num_chars) txt !, txt0
	common /cmm_txt/ n_len,txt,ib,ie
	data def_rdi/1.8/	! default radius
	
	if(acc_in.ne.file_used) then
	  def_rdi=1.8
	  num_rt=1
	  num_r_a=1
!	  aa0=' '
	  loop_1=0
80	  if(.not.open_file(8,acc_in,'old','.txt')) then
	    loop_1=loop_1+1
	    if(loop_1.gt.2) then
	      write(6,*) 
	1'acc-W> failure in opening the file ['//
	1 acc_in(:ltrim(acc_in))//']'
	      return 1
	      endif
	    acc_in=edp_data(:ltrim(edp_data))//acc_in
	    goto 80
	    endif

c	if(i_pipe.gt.0) write(6,*)	!000427
c	1 '%EdPDB-W- the current command should not ',
c	1 'be followed by a pipe structure.'
	  call read_acctxt(def_rdi,*900)

	  file_used=acc_in
	  endif
	
	resj0=' '
	r_probe= r_max
	r_max=0.
	do ii=1,n
	  i=index1(ii)
	  j=index2(ii)
	  rdi(ii)=read_rdi(i,j,resj0, def_rdi)+r_probe
	  r_max=max(r_max,rdi(ii))
	end do
	return
900	return 1	
	end
	
	function read_rdi(i,j,resj0, def_rdi)
chk	=====================
	include 'edp_main.inc'
!	use edp_main
	include 'edp_file.inc'
!	use edp_file
	include 'edp_dat.inc'
!	use edp_dat

	character*4 resj,atomi, rt, r_atom, resj0

	common /accrt/ num_rt,rt(max_rt), irt(max_rt), 
	1 num_r_a, r_atom(max_r_a), rd(max_r_a)

	logical L_rt
	
	read_rdi=def_rdi
	
	  atomi=atom(i)
	  resj=res(j)
	  L_rt=resj.ne.resj0
	  resj0=resj
	  do k=1,num_rt
	    iw=index(rt(k),wildcard)-1
	    if(iw.lt.0) then  
	      if(resj.eq.rt(k)) goto 500
	    else if(iw .eq. 0) then
	      goto 501   
	    else 
	      if(resj(1:iw).eq.rt(k)(1:iw)) goto 501
	    endif
	  end do
	  if(verbose .ge.1) write(6,1019) 
	1  ' acc-W> undefined residue ['//resj, def_rdi
	  return
	  
501	  if(L_rt) then
	    if(verbose .ge.1) write(6,1019) 
	1 ' acc-W> wildcard ['//rt(k)(1:iw+1)//'] is used for residue ['
	1//resj	
	    L_rt=.false.
	  endif  
500	  k0=irt(k)
	  k1=irt(k+1)-1
	  do k=k0,k1
	    iw=index(r_atom(k),wildcard)-1
	    if(iw.lt.0) then	! there is no wildcard
	      if(atomi.eq.r_atom(k)) goto 600
	    else if(iw.eq.0) then ! absolute wildcard
	      goto 601
	    else ! partial wildcard
	      if(atomi(1:iw).eq.r_atom(k)(1:iw)) goto 601
	    endif
	  end do
	  if(verbose .ge.1) write(6,1019) 
	1 ' acc-W> undefined atom ['//atomi//'] in residue ['//resj, def_rdi
1019	format(a,'].',:,' radius',f4.1,' is used.')
	  return
	  
601	if(verbose .ge.3) write(6,1019) 
	1 ' acc-W3> wildcard ['//r_atom(k)(1:iw+1)//'] is used for atom ['
	1 //atomi//'] in residue ['//resj, rd(k)
600	read_rdi=rd(k)
	end

copyright by X. Cai Zhang

	subroutine analy(key)
chk	================
	include 'edp_main.inc'
!	use edp_main
	include 'edp_file.inc'
!	use edp_file
	include 'edp_dat.inc'
!	use edp_dat

	real max_x,max_y,max_z,max_w,max_b
	real min_x,min_y,min_z,min_w,min_b

	parameter (max_s=256)
!	integer, parameter :: max_s=min(max_atom,256)
	real x0(max_s,5), sum_x1(5), rms_x2(5)
	integer s(max_s)
	character*1 signx,signy,signz

	n_of_syn=2					!000515
	syntax(1)='syntax:'
	syntax(2)='doit'

	k=0
	sum_x=0.
	sum_y=0.
	sum_z=0.
	sum_w=0.
	sum_b=0.
	rms_r=0.
	rms_w=0.
	rms_b=0.
	min_x=9999.
	min_y=9999.
	min_z=9999.
	min_w=9999.
	min_b=9999.
	max_x=-9999.
	max_y=-9999.
	max_z=-9999.
	max_w=-9999.
	max_b=-9999.

	if(key.eq.2) then ! angles
	  do j=1,5
	    rms_x2(j)=0.
	    sum_x1(j)=0.
	    enddo
	  do i=1,n_atom
	    if(lf(i)) then
	      k=k+1
	      if(k.gt.max_s) then
	        write(6,*) 
	1'analyze-W> UNDONE:increase max_s in analy subroutine.'
	        return
	        endif
	      read(text(i)(31:54),1061,err=800) x0(k,1), x0(k,2), x0(k,3)
	      x0(k,4)=w(i)
	      x0(k,5)=b(i)
	      do j=1,5
		tmp=mod(x0(k,j),360.)
	        x0(k,j)= mod(tmp+540.0, 360.0) -180.0
	        t= x0(k,j)
	        rms_x2(j)=rms_x2(j)+ t*t
	        sum_x1(j)=sum_x1(j) + t
	        enddo
	      write(text(i)(31:54),1061,err=800) x0(k,1), x0(k,2), x0(k,3)
	      w(i)=x0(k,4)
	      b(i)=x0(k,5)
	      endif
	    enddo
	  if(k.le.0) then
	    status=-3
		return
	  endif
	  write(6,1003) k
1003	  format(' analyze>',i6,' records included'/
	1    16x,'aver.',6x,'min.',6x,'max.',5x,'range',
	1     7x,'sgm')
	  fk=1./float(k)
	  
	  do j=1,5
	    call sort_rl(k,x0(1,j),s)
	    sgm_x=1.e32
	    do ii=1,k
	      i=s(ii)
	      if(sgm_x.gt. rms_x2(j)-sum_x1(j)*sum_x1(j)*fk) then
	        sgm_x= rms_x2(j)-sum_x1(j)*sum_x1(j)*fk
	        ave_x= sum_x1(j)*fk
	        min_x=x0(i,j)
	        if(ii.le.1) then
	          max_x= x0(s(k),j)
	        else
	          max_x= x0(s(ii-1),j)+360.
	          endif
	        endif
	      rms_x2(j)=rms_x2(j) + 360.*(2.*x0(i,j)+ 360.)
	      sum_x1(j)=sum_x1(j) + 360.
	      enddo
	    write(6,1004) j,ave_x, min_x, max_x, max_x-min_x, 
	1        sqrt(max(sgm_x*fk,0.))
	    enddo
1004	  format(9x,'a',i1,5f10.3)
	  return
800	  write(6,*) 'doit-W> error during reading xyz from the text' 
	  endif

chk	analyze cartician coordinates
	sum_xr2=0.
	sum_yr2=0.
	sum_zr2=0.
	do i=1,n_atom
	  if(lf(i)) then
	    k=k+1
	    read(text(i)(31:54),1061) xi,yi,zi
1061	    format(3f8.3)
	wi=w(i)
	bi=b(i)
	sum_x=sum_x+xi
	sum_y=sum_y+yi
	sum_z=sum_z+zi
	sum_w=sum_w+wi
	sum_b=sum_b+bi
	sum_xr2=sum_xr2+xi*(yi*yi+zi*zi)
	sum_yr2=sum_yr2+yi*(zi*zi+xi*xi)
	sum_zr2=sum_zr2+zi*(xi*xi+yi*yi)
	rms_r=rms_r+ xi*xi+ yi*yi+ zi*zi
	rms_w=rms_w+ wi*wi
	rms_b=rms_b+ bi*bi
	min_x=min(xi,min_x)
	min_y=min(yi,min_y)
	min_z=min(zi,min_z)
	min_w=min(wi,min_w)
	min_b=min(bi,min_b)
	max_x=max(xi,max_x)
	max_y=max(yi,max_y)
	max_z=max(zi,max_z)
	max_w=max(wi,max_w)
	max_b=max(bi,max_b)
	    end if
	  end do
	write( 6,1002) k
	write(48,1002) k
1002	format('!doit>',i6,' atoms included')
	if(k.le.0) then 
	  status=-3
	  return
	else
	  status=k
	endif
	fk=1./float(k)

	sum_xr2=sum_xr2*fk/rms_r
	sum_yr2=sum_yr2*fk/rms_r
	sum_zr2=sum_zr2*fk/rms_r
!	write(*,*) sum_xr2,sum_yr2, sum_zr2
	if(sum_xr2.gt. 1.e-6) then
	 signx='+'
	else if (sum_xr2.lt.-1.e-6) then
	 signx='-'
	else 
	 signx='0'
	endif
	
	if(sum_yr2.gt. 1.e-6) then
	 signy='+'
	else if (sum_yr2.lt.-1.e-6) then
	 signy='-'
	else 
	 signy='0'
	endif
	
	if(sum_zr2.gt. 1.e-6) then
	 signz='+'
	else if (sum_zr2.lt.-1.e-6) then
	 signz='-'
	else 
	 signz='0'
	endif
	
	write( 6,1001) 'x',sum_x*fk,min_x,max_x,max_x-min_x, signx,
	1         'y',sum_y*fk,min_y,max_y,max_y-min_y, signy,
	1         'z',sum_z*fk,min_z,max_z,max_z-min_z, 
	1 sqrt(max((rms_r-(sum_x*sum_x+sum_y*sum_y+sum_z*sum_z)*fk)*fk,0.)),
	1 sqrt(rms_r*fk), signz,
	1         'w',sum_w*fk,min_w,max_w,max_w-min_w,sqrt(rms_w*fk),
	1         'b',sum_b*fk,min_b,max_b,max_b-min_b,
	1 sqrt(max((rms_b-sum_b*sum_b*fk)*fk,0.))
	write(48,1001) 'x',sum_x*fk,min_x,max_x,max_x-min_x, signx,
	1         'y',sum_y*fk,min_y,max_y,max_y-min_y, signy,
	1         'z',sum_z*fk,min_z,max_z,max_z-min_z, 
	1 sqrt(max((rms_r-(sum_x*sum_x+sum_y*sum_y+sum_z*sum_z)*fk)*fk,0.)),
	1 sqrt(rms_r*fk), signz,
	1         'w',sum_w*fk,min_w,max_w,max_w-min_w,sqrt(rms_w*fk),
	1         'b',sum_b*fk,min_b,max_b,max_b-min_b,
	1 sqrt(max((rms_b-sum_b*sum_b*fk)*fk,0.))
1001	format('!',10x,'aver.',6x,'min.',6x,'max.',5x,'range',
	1 2x,'sgm(r,b)',4x,'rms(r,w)',1x,'3rd',
	1 2(/'!',4x,a1,4f10.3,24x,a1),
	1 2(/'!',4x,a1,6f10.3, 4x,a1/'!',4x,a1,4f10.3,10x,f10.3))
	if(verbose.ge.6) then
	  do i=1,n_atom
	    if(lf(i)) then 
	      write(4,1007) text(i)(1:30),
	1       sum_x*fk, sum_y*fk, sum_z*fk, sum_w*fk, sum_b*fk
	      goto 850
	    endif
	  enddo 
	endif
1007	format(A30,3f8.3,2f6.2)     	

850	if(key.eq.0) then
	  sum_x=0.
	  sum_y=0.
	  sum_z=0.
	  do i=1,n_atom
	    if(lf(i)) then
	      sum_x=sum_x +x(i)
	      sum_y=sum_y +y(i)
	      sum_z=sum_z +z(i)
	    endif
	  enddo
	  write(29,1018) 1., 0., 0., 0., 1., 0., 0., 0., 1., 
	1   -sum_x/k,-sum_y/k,-sum_z/k
	  if(verbose.ge.3) write(6,1230) rtn_out(:ltrim(rtn_out))
1230	format(' doit-I3> To apply the transformation, type:'
	1 /'   rtn file ',a)
1018	format(3(3f12.7/),3f12.5)
	endif
	end
chk***	end of analy.for

