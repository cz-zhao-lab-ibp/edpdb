	subroutine date_and_time(pdate,ptime)
chk	========================
	character*10 pdate, ptime
	integer*4 i1,i2,i3, ijk(3)
	
	call idate(i1,i2,i3)
	if(i3.le.999) i3=i3+2000
	write(pdate,1001) i3,i2,i1
1001	format(i4,i2,i2)
	if(i1.le.10) pdate(7:7)='0'
	if(i2.le.10) pdate(5:5)='0'
	
	call itime(ijk)
	write(ptime,1002) ijk
1002	format(3i2)
	if(ijk(1).le.10) ptime(1:1)='0'
	if(ijk(2).le.10) ptime(3:3)='0'
	if(ijk(3).le.10) ptime(5:5)='0'
	end

	function len_trim(string)
chk	=================
	character*(*) string
	len_trim=ltrim(string)
	end
	
	function ltrim(string)
chk	==============
	character*(*) string
	parameter (isp=32, itb=9, icm=44, inl=13, ixx=0)

	i=len(string)
	ltrim=0
	do j=i,1,-1
	  ic=ichar(string(j:j))
	  if( .not. (ic.eq. isp 
	1.or. ic.eq. itb 
	1.or. ic.eq. icm 
	1.or. ic.eq. inl
	1.or. ic.eq. ixx)) then
	    ltrim=j
	    return
	  endif
	enddo
	end

	subroutine search_text(*)
chk	===============
	include 'edp_main.inc'
!	use edp_main
	include 'edp_dat.inc'
!	use edp_dat

	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie
	logical  nword
	logical  ok

	integer k(2)

	n_of_syn=2					!000515
	syntax(1)='syntax:'
	syntax(2)='text  text.s  [p1.i  p2.i] '

c	igr:	id of the basckit group, 0-4
	call dfgroup(igr,*900)
	j1=n_groupa(igr)

	k(1)=0
	k(2)=0

	if(nword( n_len,txt,ib,ie)) return 1

	ic=ichar(delimiter)
	if(ichar(txt(ib:ib)).eq.ic) then
	  ib0=ib+1
	  j=ib+2
	  do while(ichar(txt(j:j)).ne.ic.and.j.le.n_len) 
	    j=j+1
	    enddo
	  ie0=j-1
	  ie=j
	else
	  ib0=ib
	  ie0=ie
	  endif
	length_of_text=ie0-ib0+1

	call read_ai( 2, k, *900, *800)

800	ok=.false.
	if(k(1).le.1) then
	  do j=1,j1
	    i=igroupa( j,igr)
	      if(index(text(i),txt(ib0:ie0)).gt.0) then
	        lf(i)=incl
	        ok=.true.
	        endif
	  enddo
	else 
	  if(k(2).lt.k(1)) k(2)=k(1)+length_of_text-1
	  do j=1,j1
	    i=igroupa( j,igr)
	    if(text(i)(k(1):k(2)).eq.txt(ib0:ie0)) then
	      lf(i)=incl
	      ok=.true.
	      endif
	    enddo
	  endif

	if(.not.ok) then
	  errmsg(:18)=' errmsg: in columns'
	  write(errmsg(19:25),1011) k(1), k(2)
1011	format(i3,'-',i3)	  
	  errmsg(26:)=', no strings matched ['//txt(ib0:ie0)//']'
	  return 1
	  endif
	return
900	return 1
	end
	
	subroutine max_b(*)
chk	===============
	include 'edp_main.inc'
!	use edp_main
	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie
	logical  nword

	n_of_syn=2					!000515
	syntax(1)='syntax:'
	syntax(2)='b max '

c	igr:	id of the basckit group, 0-4
	call dfgroup(igr,*900)
	j1=n_groupa(igr)

	if(nword( n_len,txt,ib,ie)) then
	  x_max=-1.e+32
	  i_max=0
	  do j=1,j1
	    i=igroupa( j,igr)
	    if(lf(i).neqv.incl.and.b(i).gt.x_max) then
	      i_max=i
	      x_max=b(i)
	      endif
	    enddo
	  if(i_max.le.0) then
	    write(6,*) 'maxb-W> UNDONE: no atom selected'
	  else
	    lf(i_max)=incl
	    endif
	else if(txt(ib:ie).eq.'for_each_residue') then
	  ir0=0
	  do j=1,j1
	    ia=igroupa( j,igr)
	    ir=aa_seq(ia)
	    if(ir.ne.ir0) then
	      if(ir0.gt.0) then
	        if(i_max.gt.0) lf(i_max)=incl
	        x_max=-1.e+32
	        i_max=0
	        endif
	      ir0=ir
	      endif
	    if(lf(ia).neqv.incl.and.b(ia).gt.x_max) then
	      i_max=ia
	      x_max=b(ia)
	      endif
	    enddo
	else 
	  return 1
	  endif
	return
900	return 1
	end

	subroutine bwxyz(ii,option0, *)
chk	================
	include 'edp_main.inc'
!	use edp_main
	character*(max_num_chars) txt
	character*1 xyz
	character*2 option
	character*(*) option0
	character*6 s_angles
	data s_angles/'angles'/

	common /cmm_txt/n_len,txt,ib,ie
!	logical  nword

	n_of_syn=2					!000515
	syntax(1)='syntax:'
	syntax(2)='(X,Y,Z,W or B) (<, >, <>, ><, =) cutoff(s).r'

	if(ii.eq.22) then
	  xyz='b'
	else if(ii.eq.23) then
	  xyz='w'
	else if(ii.eq.24) then
	  xyz='x'
	else if(ii.eq.25) then
	  xyz='y'
	else if(ii.eq.26) then
	  xyz='z'
	else
	  return 1
	  endif

	option=option0

	call read_ar( 1, v1, *900, *900)
	if(option.eq.'<>'.or.option.eq.'><') then
	  call read_ar( 1, v2, *900, *900)
	  if(v2.lt.v1) then
	    tmp=v1
	    v1=v2
	    v2=tmp
	    endif
	else if(option.eq.'<'.or.option.eq.'>') then
 	  v2=v1
	else if(option.eq.'=') then
	  tmp=0.0
	  call read_ar( 1, tmp, *900, *100)
100	  tmp=abs(tmp)+1.0e-5	  
	  v2=v1+tmp
	  v1=v1-tmp
	  option='<>'
	else 
	  return 1
	  endif

	ia=match_l1(1,s_angles)

	if(ia.lt.0) then
	  return 1
	else if(ia.eq.1) then
	  if(v2-v1.gt.360.0.or.v1.lt.-360.0.or.v2.gt.360.0) return 1
	  endif 

c	igr:	id of the basckit group, 0-4
	call dfgroup(igr,*900)
	j1=n_groupa(igr)

c	if(.not.nword( n_len,txt,ib,ie)) return 1

	do j=1,j1
	  i=igroupa( j,igr)
	  if( lf(i) .neqv. incl) then
	    if(xyz.eq.'x') then
	      read(text(i)(31:38),1061,err=901) xi
	    else if(xyz.eq.'y') then
	      read(text(i)(39:46),1061,err=901) xi
	    else if(xyz.eq.'z') then
	      read(text(i)(47:54),1061,err=901) xi
	    else if(xyz.eq.'w') then
	      xi=w(i)
	    else if(xyz.eq.'b') then
	      xi=b(i)
	    else 
	      goto 900
	      endif
	    if(option.eq.'<') then
	      if(xi.lt.v1) 	lf(i)=incl
	    else if(option.eq.'>') then
	      if(xi.gt.v1) 	lf(i)=incl
	    else if(ia.eq.0) then
	        if(option.eq.'<>') then
	          if(xi.gt.v1.and.xi.lt.v2) 	lf(i)=incl
	        else if(option.eq.'><') then
	          if(xi.lt.v1.or .xi.gt.v2) 	lf(i)=incl
	          endif
	    else if(ia.eq.1) then
	      xj=mod(xi+360.0,360.)
	      if(option.eq.'<>') then
	        if(xj.gt.v1.and.xj.lt.v2) then
	          lf(i)=incl
	        else 
	          xj=xj-360.0
	          if(xj.gt.v1.and.xj.lt.v2) lf(i)=incl
	          endif
	      else if(option.eq.'><') then
	        if(.not.(xj.gt.v1.and.xj.lt.v2)) then
	          xj=xj-360.0
	          if(.not.(xj.gt.v1.and.xj.lt.v2)) lf(i)=incl
	          endif
	        endif
	      endif
	    endif
	  enddo
	return

1061	format(f8.3)
901	write(6,*) 'EdPDB-W> ABORTED: error during reading x,y,z'
	return
900	return 1
	end

copyright by X. Cai Zhang

	subroutine dfbrg(*)
chk	================
	include 'edp_main.inc'
!	use edp_main
	include 'edp_dat.inc'
!	use edp_dat

	dimension jres(max_res),jatm(max_res)
	character*4 atoma(4), atom_x, atom_y,  
	1 lxy(0:1)*1	!, store*1
	data atoma/4*' '/, lxy/'x','y'/
	
	dimension ka(4), j_a(4), k_a(4), dat(6), iwz(2)
	data iwz/1,4/
	data ka/4*0/
	logical la(4),  la3k, la2k, ist(4),ist1234

	integer wxyz
	
	data la/4*.true./
	data dat/1.0, 4.0, 0.0, 180.0, 0.0, 360.0/
	data wxyz/1/, iskip/0/

	logical nword0, nword
	character*(max_num_chars) txt
	common /cmm_txt/ n_len,txt,ib,ie
	save

	if(verbose.ge.6 ) write(6,1069) 	!021218
1069	format(' bridge-I6> reference: for theta (elevation) angle,'
	1/'  Derewenda ZS, et al.(1995) The occurrence of C-H...O'
	1/'  hudrogen bonds in proteins. J. Mol. Biol., 252:248-62.')

	n_of_syn=5					!000515
	syntax(1)='syntax:'
	syntax(2)=
	1'define bridge atom_w.s atom_x.s atom_y.s atom_z.s '
	syntax(3)=
	1'  [(x,y) reg_w.i (x,y) reg_z.i]' 
	syntax(4)=
	1'  [status_w.l status_x.l status_y.l status_z.l]' 
	syntax(5)=
	1'  [dmin.r dmax.r amin.r amax.r tmin.r tmax.r] '//
     1'[(wxyz, zwxy,theta)] [skip.i]' 

c	atoma_k

	do ii=1,4
	  if(nword(n_len,txt,ib,ie)) then
	    if(ii.eq.1) goto 901
	    return 1
	    endif
	  atoma(ii)=txt(ib:ie)
	  ka(ii)=0
	  la(ii)=.true.
	  enddo
	
c	ka: 	ka(1/2)=0, w (or z) connects with atom_x
c		ka(1/2)=1, w (or z) connects with atom_y
c		ka(3/4)=n, relative registration number.

	do ii=1,2
	  if(nword0(n_len,txt,ib,ie)) return
	  if(ib .le. ie .and. txt(ib:ib).eq.'y' ) 
	1   ka(ii)=1                    ! atom_w connects w/ atom_y
	  if(nword0(n_len,txt,ib,ie)) return
	  if(ib .le. ie ) read(txt(ib:ie),*,err=900) ka(ii+2)
c-u1061	  format(i<ie-ib+1>)
	  enddo

c	la

	do ii=1,4
	  if(nword0(n_len,txt,ib,ie)) return
	  if(ib .le. ie) read(txt(ib:ib),*,err=900) la(ii)
c-u1062	  format(l1)
	  enddo

c	dist.(min,max), angl.(min,max), torsion.(min,max)

	do ii=1,6
	  if(nword0(n_len,txt,ib,ie)) return
	  if(ib .le. ie) read(txt(ib:ie),*,err=900) dat(ii)
c-u1063	  format(f<ie-ib+1>.0)
	  enddo

c	define torsion angle as w-x-y-z or not

	if(nword0(n_len,txt,ib,ie)) return
	if(ib .le. ie) then 
	  if(txt(ib:ie).eq.'wxyz') then
	    wxyz=1
	  else if(txt(ib:ie).eq.'zwxy') then
	    wxyz=2
	  else if(txt(ib:ib+3).eq.'thet') then
	    wxyz=3
	  else
	    return 1
	    endif
	  endif
	
	call read_ai(1,iskip,*900,*899)
899	return
900	return 1
c
	entry shbrg
c	===========
901	if(wxyz .eq. 3) then
	  atom_x='thet'
	else if(wxyz .eq. 2) then
	  atom_x='zwxy'
	else 
	  atom_x='wxyz'
	  endif

	 write(6,1002)  atom_x, (atoma(i),i=1,4), 
	1 (lxy(ka(i)),ka(i+2),i=1,2),
	1 (la(i),i=1,4), (dat(i),i=1,6),iskip
1002	 format(' ',a4,'> ',t10,4a5,t30,2(1x,a1,i3),t40,4(1x,l1),
	1 t50,2f5.1,4f5.0/ t10,'skip',i4,' residues')
	return

	entry brg(*)
c	=========
c	the current version looks for 
c	w(first match in a residue)
c	x(multiple in a residue)
c	y(multiple in a residue)
c	z(first match in a residue)
	
	n_of_syn=2					!000515
	syntax(1)='syntax:'
	syntax(2)='bridge [(w, x, y, z)]' 

	do i=1,4
	  if(atoma(i).eq.' ') then
	    write(6,*) 
	1 'bridge-W> define atom1 -  atom4 first, using define.'
	    return
	    endif
	  enddo

	ist(1)=.false.
	ist(2)=.false.
	ist(3)=.false.
	ist(4)=.false.
	if(.not.nword(n_len,txt,ib,ie)) then
	  do i=ib, ie
	    if(txt(i:i).eq.'w') then
		  ist(1)=.true.
	    else if(txt(i:i).eq.'x') then
		  ist(2)=.true.
	    else if(txt(i:i).eq.'y') then
		  ist(3)=.true.
	    else if(txt(i:i).eq.'z') then
		  ist(4)=.true.
	    else
	      return 1
		  endif
		enddo
	  endif
	ist1234=ist(1).or.ist(2).or.ist(3).or.ist(4)

c-
	ksum=0
	jsum=0
	avd=0.
	ava=0.
	avt=0.
	sgmd=0.
	sgma=0.
	sgmt=0.
	d_max=-1.0E6
	d_min=1.0E6
	a_max=-1.0E6
	a_min=1.0E6
	t_max=-1.0E6
	t_min=1.0E6

	ii0=0
	ii1=0
	do ii=1,2
	  if(ka(ii).eq.1) then
	    ii0=min(ii0,ka(ii+2))
	    ii1=max(ii1,ka(ii+2))
	    endif
	  enddo
	i0=1-ii0
	i1=n_res-ii1

c	list atom_y
	num_y=0
	atom_y=atoma(3)
	itrim=index(atom_y,wildcard)-1
	
	la3k=.not.la(3)
	do 20 i=i0,i1
	  j0=ijk(i)
	  j1=ijk(i+1)-1
	  do j=j0,j1
	   if(la3k.or.lf(j)) then
	     if( (itrim.le.0 .and. atom(j).eq.atom_y) .or.
	1        (itrim.gt.0 .and. 
	1         atom(j)(:itrim).eq.atom_y(:itrim)))  then 
	       num_y=num_y+1
	       jres(num_y)=i
	       jatm(num_y)=j
	       endif
	    end if
	  enddo
20	enddo

c-	select atom_x
	ii0=0
	ii1=0
	do 30 ii=1,2
	  if(ka(ii).eq.0) then
	    ii0=min(ii0,ka(ii+2))
	    ii1=max(ii1,ka(ii+2))
	    endif
30	  enddo
	i0=1-ii0
	i1=n_res-ii1

	atom_x=atoma(2)
	itrim2=index(atom_x,wildcard)-1
	
	la2k=.not.la(2)
	do 800 i=i0,i1
	  j0=ijk(i)
	  j1=ijk(i+1)-1
	  do 700 j700=j0,j1
	    if(la2k.or.lf(j700)) then
	      if( (itrim2.le.0 .and. atom(j700).eq.atom_x) .or.
	1         (itrim2.gt.0 .and. 
	1        atom(j700)(:itrim2).eq.atom_x(:itrim2))) then
	        xx=x(j700)
	        xy=y(j700)
	        xz=z(j700)
	        j_a(2)=j700

c	select atom_w,z if it connects w/ atom_x
	do 300 i300=1,2
	  if(ka(i300).eq.0) then	! ka=0, connect to x; ka=1, connect to y
	    j0=ijk(i+ka(i300+2))
	    j1=ijk(i+ka(i300+2)+1)-1
	    la3k=.not.la(iwz(i300))
	    	
	    atom_y=atoma(iwz(i300))
	    itrim=index(atom_y,wildcard)-1
	    
	    do 200 j200=j0,j1
	      if(la3k .or. lf(j200))  then
	        if( (itrim.le.0 .and. atom(j200).eq.atom_y) .or.
	1           (itrim.gt.0 .and. 
	1            atom(j200)(:itrim).eq.atom_y(:itrim))) then
		  j_a(iwz(i300))=j200 
		  goto 300  	!021218 
		  endif		! match atom names
	        endif		! select on atoms
200	      enddo		! loop all atoms in the residue (see 300)
	    goto 800
	  end if		! select only atom(s) that connected to atom_x
300	enddo			! find atom_w and/or atom_z


c	select atom_y

	do 500 j=1, num_y
	jres_j=jres(j)
	if(abs(jres_j-i).lt.iskip) goto 500
	jatm_j=jatm(j)
	if(j_a(2).eq.jatm_j) goto 500
	dx=xx-x(jatm_j)
	if( abs(dx) .gt. dat(2) ) goto 500
	dy=xy-y(jatm_j)
	if( abs(dy) .gt. dat(2) ) goto 500
	dz=xz-z(jatm_j)
	if( abs(dz) .gt. dat(2) ) goto 500
	r=sqrt( dx*dx+dy*dy+dz*dz)
	if( r.lt.dat(1) .or. dat(2).lt.r ) goto 500
	j_a(3)=jatm_j
c
c	select atom_w,z if it connects w/ atom_y
c
	do 450 ii=1,2
	 if(ka(ii).eq.1) then		! ka=0, connect to x; ka=1, connect to y
	  j0=ijk(jres_j+ka(ii+2))
	  j1=ijk(jres_j+ka(ii+2)+1)-1
	  la3k=.not.la(iwz(ii))
	  	
	  atom_y=atoma(iwz(ii))
	  itrim=index(atom_y,wildcard)-1
	  
	  do 400 jj=j0,j1
	   if(la3k .or. lf(jj))  then
	     if( (itrim.le.0 .and. atom(jj).eq.atom_y) .or.
	1        (itrim.gt.0 .and. 
	1         atom(jj)(:itrim).eq.atom_y(:itrim)))  then
	       j_a(iwz(ii))=jj
	       goto 450				!021218
	       endif   
	     end if
400	  enddo				! match atom names
	  goto 500			! if there is no atom_w or atom_z, go to next atom_y
	 end if				! select atom_w and atom_z only if they are connect to atom_y
450	enddo				! find atom_w and/or atom_z that connect to atom_y
	
	if(j_a(1).eq.j_a(2)) goto 500
	if(j_a(1).eq.j_a(3)) goto 500
	if(j_a(1).eq.j_a(4)) goto 500
	if(j_a(2).eq.j_a(3)) goto 500
	if(j_a(2).eq.j_a(4)) goto 500
	if(j_a(3).eq.j_a(4)) goto 500
c
	angl=f_angl(j_a)
	if(angl.lt.dat(3).or. dat(4).lt.angl) goto 500
c
	if(wxyz .ge.2) then
	  k_a(1)=j_a(4)
	  k_a(2)=j_a(1)
	  k_a(3)=j_a(2)
	  k_a(4)=j_a(3)
	  tor=f_tor(k_a)
	  if(wxyz .eq. 3) tor=asind(sind(angl)*sind(tor))
	else 
	  tor=f_tor(j_a)
	end if

	if(tor.gt.360.) goto 500
	if(tor.lt.dat(5)) tor=tor+360.
	if(tor.lt.dat(5) .or. dat(6).lt. tor) goto 500
c
	ksum=ksum+1
	avd=avd+r
	ava=ava+angl
	avt=avt+tor
	sgmd=sgmd+r*r
	sgma=sgma+angl*angl
	sgmt=sgmt+tor*tor
	d_max=max(d_max, r)
	d_min=min(d_min, r)
	a_max=max(a_max, angl)
	a_min=min(a_min, angl)
	t_max=max(t_max, tor)
	t_min=min(t_min, tor)
	write(48,1001) (text(j_a(ii))(14:27), ii=1,4),
	1 r,angl,tor
	do jst=1,4
	  if(ist(jst)) then	
	    jsum=jsum+1
		igroupa( jsum,1)=j_a(jst)
		endif
	  enddo
1001	format(1x,4(a14,','),'dat=',f5.2,2f5.0)
500	enddo			! loop through potential atom_y

c		  endif		! match atom names
c	        endif		! select on atoms
c200	      enddo		! loop all atoms in the residue (see 300)
cc	    goto 800
c	  end if		! select only atom(s) that connected to atom_x
c300	enddo			! find atom_w and/or atom_z


	        endif		! find atom_x (match names)
	      endif		! find atom_x (select on atoms)
700	    enddo		! loop all atoms in the residues (see 800)
800	enddo			! loop all residues

	if(ist1234) n_groupa(1)=jsum
	if(wxyz .eq. 3) then
	  atom_x='thet'
	else if(wxyz .eq. 2) then
	  atom_x='zwxy'
	else 
	  atom_x='wxyz'
	  endif
	if(ksum.eq.0)then
	  write(6,1009) atom_x,ksum
	else
	  if(ksum .gt. 1) then
	    sumk=1./ksum
	    avd=avd*sumk
	    ava=ava*sumk
	    avt=avt*sumk
	    sgmd=sgmd*sumk-avd*avd
	    sgma=sgma*sumk-ava*ava
	    sgmt=sgmt*sumk-avt*avt
	    if(sgmd .gt. 0.0) sgmd=sqrt(sgmd)
	    if(sgma .gt. 0.0) sgma=sqrt(sgma)
	    if(sgmt .gt. 0.0) sgmt=sqrt(sgmt)
	  else
	    sgmd=0.0
	    sgma=0.0
	    sgmt=0.0
	  endif
	  write( 6,1009) atom_x, ksum, avd, ava, avt,sgmd,sgma,sgmt
	  write(48,1009) atom_x, ksum, avd, ava, avt,sgmd,sgma,sgmt
	  if(verbose .ge.2) then
	    write( 6,1011) atom_x,d_min,a_min,t_min, d_max,a_max,t_max
	    write(48,1011) atom_x,d_min,a_min,t_min, d_max,a_max,t_max
	    endif
	end if
1009	format('^',a4,'> ',i6,' records are calculated':
	1     / '        avdat=(',3f8.3,'), sgmdat=(',3f8.3,')')
1011	format(1x,a4,'-I2> min=(',3f8.3,'),    max=(',3f8.3,')')
	call typelist(48)
	end

copyright by X. Cai Zhang
c	this subroutine selects all
c	the atoms of residue names
c	starting with the chain marker

chk	example input:
chk	1) chain  --list chains
chk	2) chain a b -- select chains a and b
chk	3) chain a-z -- select chains a to z
chk	4) chain a from {ca} -- select Ca atoms from chain a 

	subroutine chain(*)
chk	================
	include 'edp_main.inc'
!	use edp_main
	include 'edp_dat.inc'
!	use edp_dat

	logical ok
	integer n_chain(65:122)

	character*38 cm
	integer ncm

	logical  nword
	character*(max_num_chars) txt, ch*1
	common /cmm_txt/n_len,txt,ib,ie
	save

	n_of_syn=3					!000515
	syntax(1)='syntax:'
	syntax(2)='1) chain [chn_mark_1.x [chn_mark_2.x ...]]' 
	syntax(3)='2) chain chn_mark_1.x-chn_mark_2.x' 
c	syntax(4)='3) chain 1-9'  ! select atoms of a blank chain_mark

chk	igr:    id of the basckit group, 0-4
	call dfgroup(igr,*900)
	j1=n_groupa(igr)

	n=1
	numOfOkayCM=0
100	if(nword(n_len,txt,ib,ie)) then
	  if(n.gt.1) return
chk	  else list the # of ON atom is each chain.
	    do k=65,122
	      n_chain(k)=-1
	      enddo
	    do ir=1,n_res
	      k=ichar(ires(ir)(1:1))
	      if(k.ge.65.and.k.le.122) then
	        n_chain(k)=max(n_chain(k),0)
	        do ia=ijk(ir),ijk(ir+1)-1
	          if(lf(ia)) n_chain(k)=n_chain(k)+1
	          enddo
	        endif
	      enddo
	    kt=0
	    ks=0
	    do k=65,122
	      if(n_chain(k).ge.0) then
	        if(n_chain(k).eq.0) then
	          write(6,1079) char(k),n_chain(k)
1079	format(2x,'chain ',a,' [',i6,']',:,t50,'chain-I3> ',2i6)
	        else
	          kt=kt+1
	          ks=ks+ n_chain(k)
		  if(verbose.lt.3) then
	            write(6,1079) char(k),n_chain(k)
		  else
	            write(6,1079) char(k),n_chain(k),kt,ks
		  endif
	        endif		
	      endif
	    enddo
	  return
	endif

	if(ib.gt.ie) return 1
	cm=txt(ib:ie)
	ncm=ie-ib+1
	call extent_chain_mark(cm, ncm,*900)

	do k=1,ncm
	  ch=cm(k:k)
	  n=n+1

	  ok=.false.
	  do j=1,j1
	    i=igroupa( j,igr)
	    ir=aa_seq(i)
	    if(ires(ir)(1:1).eq.ch) then
	      ok=.true.
	      lf(i)=incl
	      endif
	    enddo

	  if(.not.ok) then
	    do ir=1,n_res
	      if(ires(ir)(1:1).eq.ch) then
	        numOfOkayCM=numOfOkayCM+1
	        goto 100
	        endif
	      enddo
	    if(verbose.ge.1) then
	      write(6,*) 'chain-W> chain ['//ch//'] does not exist'
	      endif
	  else
	    numOfOkayCM=numOfOkayCM+1
	    endif
	  enddo
	if(numOfOkayCM.le.0) goto 900
	goto 100

900	return 1

	entry more_chain(jgr)
chk	================
	do k=65,122
	  n_chain(k)=-1
	  enddo

	do ir=1,n_res
	  k=ichar(ires(ir)(1:1))
	  if(k.ge.65.and.k.le.122) then
	    n_chain(k)=max(n_chain(k),0)
	    do ia=ijk(ir),ijk(ir+1)-1
	      if(Lf(ia)) n_chain(k)= n_chain(k)+1
	      enddo
	    endif
	  enddo
	
	do j=1,n_groupa(jgr)
	  ia=igroupa(j,jgr)
	  ir=aa_seq(ia)
	  k=ichar(ires(ir)(1:1))
	  if(k.ge.65.and.k.le.122) then
	    if(n_chain(k).gt.0) lf(ia)= .true.
	    endif
	  enddo
	end

copyright by X. Cai Zhang

	subroutine punch_sym_o(io,num_sym,sym0)
chk	======================
	include 'edp_main.inc'
!	use edp_main

	real sym0(3,4,num_sym)

	n_of_syn=2					!000515
	syntax(1)='syntax:'
	syntax(2)='symmetry punch file_name.s [o]' 

	if(num_sym.le.0) then
	  write(6,*) '%EdPDB-W- symmetry has not been defined.'
	  return
	  endif

	write(io,1001) num_sym*12
1001	format('.space_group_operators     r  ',i4,'  (3F10.5)')

	write(io,1002) ((( sym0(i,j,isym), i=1,3),j=1,4),isym=1,num_sym)
1002	format(3f10.5)
	end
	
	subroutine chksymm(io,num_sym,sym0,trn)
chk	==================
	include 'edp_main.inc'
!	use edp_main
c                ^^^^^^^ parameter max_symm=48

	real sym0(3,4,max_symm), trn(3,4,max_symm)

	character*32 symtxt, comment(max_symm)       
	integer nsym, mm(max_symm), n6(max_symm)
	integer  k6(max_symm), m6(max_symm)
	character*1 element(-1:1)
	data element/' ',' ','*'/

	real  sym(3,4,max_symm), sym6(3,4,6,max_symm)

	real tmp(3,4), u(3,4)
	data u/1.,0.,0., 0.,1.,0., 0.,0.,1., 0.,0.,0./

	logical mtxsml
	external mtxsml

chk	input symmetry operators
	nsym=num_sym
	do isym=1,nsym
	  call mtxcpy(sym0(1,1,isym), sym(1,1,isym))
	  call get_det(trn(1,1,isym), det, a1, a2, a3, b1, b2, b3)
	  comment(isym)=' '
	  if( abs( det-1.0) .gt. 1.e-3) then
	    write(comment(isym),1007,err=40) det 
1007	    format(' det =',f8.3)
	  else if( 
	1       abs(a1-1.0) .gt. 1.e-3 .or.
	1       abs(a2-1.0) .gt. 1.e-3 .or.
	1       abs(a3-1.0) .gt. 1.e-3 .or.
	1       abs(b1-1.0) .gt. 1.e-3 .or.
	1       abs(b2-1.0) .gt. 1.e-3 .or.
	1       abs(b3-1.0) .gt. 1.e-3 ) then
	    write(comment(isym),1008,err=40) 
1008	    format(' inperfect rotation')
	    endif
40	  enddo

50	if(nsym.le.0) then
	  write(6,*) '%EdPDB-W- symmetry has not been defined.'
	  return
	  endif

chk	check "close"
	do isym=1,nsym
	  call mtxcpy(u, tmp)
	  do j=1,6
	    call mtxmul( sym(1,1,isym), tmp, tmp)
	    call mtxcpy( tmp, sym6(1,1,j,isym) )
	    if( mtxsml(tmp,u) ) then
	      mm(isym)=1
	      n6(isym)=j
	      goto 60
	      endif
	    enddo
	  n6(isym)=0
	  mm(isym)=0
	  comment(isym)=' non-closed element'
60	  enddo

chk	check "duplication"
	do isym=1,nsym-1
	  if(mm(isym).ne.0) then
	    do jsym=isym+1,nsym
	      if(mm(jsym).ne.0) then
	        if(mtxsml(sym(1,1,isym), sym(1,1,jsym))) then
	          mm(jsym)=0
	          write(comment(jsym),1004) isym
1004		  format(' the same as #',i2)
	          endif
	        endif
	      enddo
	    endif
	  enddo
	
chk	check "completeness"
	msym=nsym	      
	do isym=1,msym
	  if(mm(isym).ne.0) then
	do jsym=1,msym
	  if(mm(jsym).ne.0) then
	  call mtxmul( sym(1,1,isym), sym(1,1,jsym), tmp)
	  do ksym=1,msym
	    if(mtxsml( tmp, sym(1,1,ksym) )) goto 80
	    enddo
	  if(nsym.ge.max_symm) then
	    write(6,*) '%EdPDB-W- num_sym => max_symm=',max_symm
	    goto 85
	    endif
	  nsym=nsym+1
	  call mtxcpy(tmp, sym(1,1,nsym))
	  write(comment(nsym),1009) isym, jsym
1009	  format(' s(',i2,')s(',i2,')')
	  goto 50
	    endif
80	  enddo
	    endif
	  enddo

chk	dismiss non-basic elements
	do i=1,nsym
	  k6(i)=n6(i)*max_symm+(max_symm-i)
	  enddo
	call sort_in2(nsym, k6, m6)

	do i=nsym,2,-1
	  isym=m6(i)
	  if(mm(isym).ne.0) then
	    do i1=2,n6(isym)
	      do jsym=1,nsym
	        if(mtxsml(sym(1,1,jsym),sym6(1,1,i1,isym))) then
	          mm(jsym)=0
	          goto 803
	          endif
	        enddo
803	      enddo
	    endif
	  enddo
	    
	do i=nsym,2,-1
	  isym=m6(i)
	  if(mm(isym).ne.0) then
	    do j=nsym,2,-1
	      jsym=m6(j)
	      if(mm(jsym).ne.0.and.
	1       mod(nsym,n6(isym)*n6(jsym)).eq.0) then
	        do i1=1,n6(isym)
	        do j1=1,n6(jsym)
	          call mtxmul(sym6(1,1,i1,isym),sym6(1,1,j1,jsym),tmp) 
	          do ksym=1,nsym
	            if(ksym.ne.isym.and.ksym.ne.jsym.and.
	1             mtxsml(tmp,sym(1,1,ksym))) then
		      mm(ksym)=0
	              goto 801
	              endif
	            enddo
801	          enddo
	          enddo
	        endif
	      enddo
	    endif
	  enddo

	n=nsym
	do i=nsym,2,-1
	  isym=m6(i)
	  if(mm(isym).gt.0) then
	    if(n.eq.1.or. mod(n,n6(isym)).ne.0) then
	      mm(isym)=0
	    else
	      n=n/n6(isym)
	      endif
	    endif
	  enddo

85	if(io.ne.0) write(io,*) 'symmetry initialize'
	do isym=1,nsym
	  call get_symm_string(sym(1,1,isym), symtxt)
	  if(io.ne.0)
	1 write(io,1003) 
	1   symtxt, isym, element(mm(isym)), n6(isym), comment(isym)
	  write(48,1003) 
	1   symtxt, isym, element(mm(isym)), n6(isym), comment(isym)
1003	  format(' symmetry ',a32,'!',i2,a1,' (ord=',i2,')',a32)
c1006	  format(' symmetry ',a32,'#',i2,a1,' (ord=',i2,')',a32)
	  enddo
	call typelist(48)	 
	end

	subroutine mtxmul( a,b,c)
chk	=================
	real a(3,4), b(3,4), c(3,4), t(3,4)
	
	do i=1,3
	do j=1,4
	  t(i,j)=a(i,1)*b(1,j)+ a(i,2)*b(2,j)+ a(i,3)*b(3,j)
	  enddo
	  t(i,4)=t(i,4)+a(i,4)
	  enddo


	do i=1,3
	do j=1,4
	  c(i,j)=t(i,j)
	  enddo
	  enddo
	end

	logical function mtxsml(a,b)
chk	=======================
	real a(12), b(12)
	parameter (eps=1.e-6)
!	real, parameter :: eps=1.e-6

	mtxsml=.false.
	do i=1,9
	  if(abs(a(i)-b(i)).gt.eps) return
	  enddo

	do i=10,12
	  t=mod(abs(a(i)-b(i)),1.)
	  if(t.gt.eps .and. (1.-t).gt.eps) then
	    return
	    endif
	  enddo

	mtxsml=.true.
	end
	
	subroutine mtxcpy(a,b)
chk	=================
	real a(12), b(12)

	do i=1,12
	  b(i)=a(i)
	  enddo	
	end

	subroutine get_symm_string(a, string)
chk	==========================
	character*64 tmp
	character*32 string
	real a(3,4)
	character*2 s(-1:1,3)
	data s/'-x',' ','x','-y',' ','+y','-z',' ','+z'/
	character*6 t(0:11)
	data t/' ','+1/12','+1/6','+1/4','+1/3','+5/12',
	1      '+1/2','+7/12','+2/3','+3/4','+5/6','+11/12'/

	do i=1,3
	do j=1,3
	  if(abs(a(i,j)).gt.1.1) then
	    string=' wrong matrix '
	    return
	    endif
	  enddo
	  enddo

	tmp=
	1 s(nint(a(1,1)),1)//s(nint(a(1,2)),2)//s(nint(a(1,3)),3)//
	1 t(mod(nint(a(1,4)*12.)+120,12))//','//
	1 s(nint(a(2,1)),1)//s(nint(a(2,2)),2)//s(nint(a(2,3)),3)//
	1 t(mod(nint(a(2,4)*12.)+120,12))//','//
	1 s(nint(a(3,1)),1)//s(nint(a(3,2)),2)//s(nint(a(3,3)),3)//
	1 t(mod(nint(a(3,4)*12.)+120,12))

	j=0
	do i=1,64
	  if(tmp(i:i).ne.' ' ) then
	    j=j+1
	    if(j.le.32) string(j:j)=tmp(i:i)
	    if(tmp(i:i).eq.',' ) then
	      j=j+1
	      if(j.le.32) string(j:j)=' '
	      endif
	    endif
	  enddo
	if(j.lt.32) string(j+1:)=' '
	end

	subroutine get_det(a, det, a1, a2, a3, b1, b2, b3)
chk	==================
	real a(3,3)

	det= a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))
	1   +a(2,1)*(a(3,2)*a(1,3)-a(3,3)*a(1,2))
	1   +a(3,1)*(a(1,2)*a(2,3)-a(1,3)*a(2,2))

	a1 = a(1,1)*a(1,1)+ a(1,2)*a(1,2)+ a(1,3)*a(1,3)
	a2 = a(2,1)*a(2,1)+ a(2,2)*a(2,2)+ a(2,3)*a(2,3)
	a3 = a(3,1)*a(3,1)+ a(3,2)*a(3,2)+ a(3,3)*a(3,3)

	b1 = a(1,1)*a(1,1)+ a(2,1)*a(2,1)+ a(3,1)*a(3,1)
	b2 = a(1,2)*a(1,2)+ a(2,2)*a(2,2)+ a(3,2)*a(3,2)
	b3 = a(1,3)*a(1,3)+ a(2,3)*a(2,3)+ a(3,3)*a(3,3)

	end

copyright by X. Cai Zhang

	subroutine mxinv(n,a,ierr)
chk	================
	parameter (nl=100)
!	integer, parameter :: nl=100
	integer ip(nl),jp(nl)
	real a(n,n),b(nl),c(nl), w, z

	ierr=0

	do k=1,n
	  w=0.
	  do j=k,n			! to find the maximum element.
	  do i=k,n			!  in the sub-matrix.
	    if(abs(a(i,j)).gt.abs(w)) then
	      w=a(i,j)
	      i0=i
	      j0=j
	      endif
	    enddo
	    enddo

chk	 set the max element (i0,j0) to be the first element in the sub-matrix.

	  if(abs(w).lt.eps())then		! the determine of the matrix is 0.0
	    ierr=1
	    return 
	    endif

	  w=1./w
	  if(i0.ne.k) then 		! permute the i0, k lines
	    do j=1,n			
	      z=a(i0,j)
	      a(i0,j)=a(k,j)
	      a(k,j)=z
	      enddo
	    endif
	  if(j0.ne.k) then 		! permute the j0, k columes
	    do i=1,n
	      z=a(i,j0)
	      a(i,j0)=a(i,k)
	      a(i,k)=z
	      enddo
	    endif
	  ip(k)=i0
	  jp(k)=j0

	  do j=1,n
	    if(j.eq.k)then
	      b(j)=w
	      c(j)=1.
	    else
	      b(j)=-a(k,j)*w
	      c(j)=a(j,k)
	      endif
	    a(k,j)=0.
	    a(j,k)=0.
	    enddo

	  do j=1,n	
	  do i=1,n
	    if(abs(a(i,j)).lt.eps()) a(i,j)=0.0
	    a(i,j)=a(i,j)+b(j)*c(i)
	    enddo
	    enddo

	  enddo

chk	set the order of the elements of the inverse matrix to the origina one.

	do k=n,1,-1
	  i0=ip(k)
	  j0=jp(k)
	  if(i0.ne.k) then		! permut back the i0,k columes
	    do i=1,n			
	      z=a(i,i0)
	      a(i,i0)=a(i,k)
	      a(i,k)=z
	      enddo
	    endif
	  if(j0.ne.k) then		! permut back the j0,k lines
	    do j=1,n
	      z=a(j0,j)
	      a(j0,j)=a(k,j)
	      a(k,j)=z
	      enddo
	    endif
	  enddo
	end

	function nword(n_len,txt,ib,ie)
chk	==============
! nword return true for a comma (,)
	logical nword, nword0
	character*(*) txt

!	if(ib .gt. n_len) then
!	  nword=.true.
!	else
	  nword=nword0(n_len,txt,ib,ie)
	  if(ib  .gt.ie) nword=.true.
!	endif
	end

	function nword0(n_len,txt,ib,ie)
chk	===============
	logical nword0
	character*(*) txt
	parameter (isp=32, itb=9, icm=44)

	if(ie.le.0) then
	  ib=1
	else
	  ib=ie+1
	  endif

	i=ib
	ic=ichar(txt(i:i))
	if(ic.eq.icm.and.i.gt.1) i=i+1	! belongs to the previous search

	do while (.true.)
	  if(i.gt.n_len) goto 50
	  ic=ichar(txt(i:i))
	  if(   ic.ne.isp .and. ic.ne.itb) goto 100
	  i=i+1
	  enddo
50	nword0=.true.
	return

100	ib=i
	do i=ib,n_len
	ic=ichar(txt(i:i))
	if(  ic.eq.isp .or.ic.eq.icm
	1 .or.ic.eq.itb .or.ic.eq.0) goto 200
	end do
200	ie=i-1
	nword0=.false.
	end

 	subroutine read_ai(n,a,*,*)
chk	==================
	include 'edp_main.inc'
!	use edp_dim

	integer a(n), t
chk	return 1: wrong word
chk	return 2: no word,	does not change the input

	parameter (isp=32, itb=9, icm=44)

	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie

	miss=0
	do j=1,n
	  if(ie.le.0) then
	    ib=1
	  else
c	    ib=ie+2	
	    ib=ie+1
	    endif
	   
	  ic=ichar(txt(ib:ib))
	  if(ic.eq.icm.and.ib.gt.1) ib=ib+1	! belongs to the previous search
	  do i=ib,n_len
	    ic=ichar(txt(i:i))
	    if( ic.ne.isp .and. ic.ne.itb) goto 100
	    end do
	  return 2

100	  ib=i
	  do i=ib,n_len
	    ic=ichar(txt(i:i))
	    if(  ic .eq.isp
	1    .or.ic .eq.icm
	1    .or.ic .eq.itb) goto 200
	  enddo

200	  ie=i-1
	  if(ib  .gt.ie) then
	    miss=miss+1
	  else
	    read(txt(ib:ie),*,err=50) t
	    a(j)=t
	    endif
	  enddo
	if(miss.gt.0) return 2
	return
50	return 1
	end

	subroutine symm(a,iq,sym,ierr)
chk	===============
chk	x'= sym(1,1)*x+sym(1,2)*y+sym(1,3)*z +sym(1,4)
chk	y'= sym(2,1)*x+sym(2,2)*y+sym(2,3)*z +sym(2,4)
chk	z'= sym(3,1)*x+sym(3,2)*y+sym(3,3)*z +sym(3,4)

chk	h'= h*sym(1,1)+k*sym(2,1)+l*sym(3,1)
chk	k'= h*sym(1,2)+k*sym(2,2)+l*sym(3,2)
chk	l'= h*sym(1,3)+k*sym(2,3)+l*sym(3,3)

	include 'edp_dim.inc'
!	use edp_dim

	character*(max_num_chars) errmsg
	character*(max_num_chars) syntax(32)
	common /error_msg/ errmsg, max_err, n_err, syntax, n_of_syn

	character*(*) a
	real sym(3,4)
	logical l1,l2
	integer hkl
	character*1 ai, sp,cm,tb, pl, ng, h,k,l,sx,sy,sz
	data sp,tb,cm,pl,ng/' ',' ',',','+','-'/
	data h,k,l,sx,sy,sz/'h','k','l','x','y','z'/
	tb=char(9)

	n_of_syn=6					!000515
	syntax(1)='syntax:' 
	syntax(2)='1) symmetry -- to show existing operators' 
	syntax(3)='2) symmetry symmetry_operator.s -- e.g. x+1/2,-y,-z' 
	syntax(4)='3) symmetry punch file_name.s [o]' 
	syntax(5)='  to erase existing symmetry operators, type'
	syntax(6)='  cell ,,,,,,,'

c
chk	initializing
c
	ierr=0
	do i=1,3
	do j=1,3
	sym(i,j)=0.
	end do
	end do
	l1=.true.	! true: the element can not be ended.
	l2=.true.	!no deviding symbol
	c=1.
	n=1
	hkl=0

	do 500 i=1,iq
	  ai=a(i:i)
	  jc=ichar(ai)
	  if(65.le.jc.and.jc.le.90) ai=char(jc+32)

!	if((ai.eq.sp.or.ai.eq.tb).and.l1) goto 500
	if((jc.eq.32.or.jc.eq. 9).and.L1) goto 500

!	if(ai.eq.cm.and.l1) goto 600
	if(jc.eq.44.and.l1) goto 600

!	if(ai.eq.cm.or.ai.eq.sp.or.ai.eq.tb) then
	if(jc.eq.44.or.jc.eq.32.or.jc.eq. 9) then
	  if(.not.l2) goto 600
	  if(n.ge.3) goto 900
	  n=n+1
	  l1=.true.
	  goto 500
	  endif

!	if(ai.eq.pl) goto 500
	IF(jc.eq.43) GOTO 500

!	if(ai.eq.ng) then
	IF(jc.eq.45) THEN
	  c=-1.
	  goto 500
	  endif

	if(ai.eq.'/') then
	  l2=.false.
	  goto 500
	  endif

	j=0
	if(ai.eq.h) j=1
	if(ai.eq.k) j=2
	if(ai.eq.l) j=3
	if(j.ne.0) goto 450

	if(hkl.gt.0) goto 600
	if(hkl.eq.0) then
	  hkl=-1
	  sym(1,4)=0.
	  sym(2,4)=0.
	  sym(3,4)=0.
	  endif

	if(ai.eq.sx) j=1
	if(ai.eq.sy) j=2
	if(ai.eq.sz) j=3
	if(j.ne.0) goto 400

	b=0.
	if(ai.eq.'1') b=1.
	if(ai.eq.'2') b=2.
	if(ai.eq.'3') b=3.
	if(ai.eq.'4') b=4.
	if(ai.eq.'5') b=5.
	if(ai.eq.'6') b=6.
	if(b.le.0.) then
	  ierr=1
	  return 
	  endif

	if(l2) then
	 sym(n,4)=c*b
	 c=1.
	else
	 sym(n,4)=sym(n,4)/b
	 l2=.true.
	 endif
	goto 500

400	sym(n,j)=c
	c=1.
	l1=.false.
	goto 500

450	if(hkl.lt.0) goto 600
	hkl=1
	sym(j,n)=c
	c=1.
	l1=.false.
500	continue
600	if(n.lt.3.or.n.eq.3.and.l1) then
	  ierr=1
	  endif
900	end

	subroutine trnsln0(cell,trn)
chk	==================
	real trn(3,3),cell(6)
	logical log
	save

	log=.true.
	goto 100

	entry trnsln1(cell,trn)
chk	=============
	log=.false.

100	a=	cell(1)
	b=	cell(2)
	c=	cell(3)
	alf=	cell(4)
	beta=	cell(5)
	gama=	cell(6)
	if( alf.eq.0.) alf=90.
	if(beta.eq.0.)beta=90.
	if(gama.eq.0.)gama=90.
c	cos_a=cosd( alf +0.00001)	! for the stupid alpha compiler
c	cos_b=cosd(beta +0.00001)
c	cos_g=cosd(gama +0.00001)
c	sin_a=sind( alf +0.00001)
	cos_a=cosd( alf)
	cos_b=cosd(beta)
	cos_g=cosd(gama)
	sin_a=sind( alf)
	c2_a=cos_a*cos_a
	c2_b=cos_b*cos_b
	c2_g=cos_g*cos_g
	c2_f=(1. -c2_a -c2_b -c2_g +2.*cos_a*cos_b*cos_g)/(sin_a*sin_a)
	cos_f= sqrt(c2_f)
	cos_p=-(cos_b - cos_a*cos_g)/sin_a		!000706
c	cos_p= sqrt( max(0., 1. -c2_g -c2_f))

	if(log) then
	trn(1,1)=1./(a*cos_f)
	trn(2,1)=-(cos_g+cos_a*cos_p/sin_a)/(cos_f*b)
	trn(3,1)=cos_p/(sin_a*cos_f*c)
	trn(1,2)=0.
	trn(2,2)=1./b
	trn(3,2)=0.
	trn(1,3)=0.
	trn(2,3)=-cos_a/(sin_a*b)
	trn(3,3)=1./(c*sin_a)
	else
	trn(1,1)=a*cos_f
	trn(2,1)=a*cos_g
	trn(3,1)=-a*cos_p
	trn(1,2)=0.
	trn(2,2)=b
	trn(3,2)=0.
	trn(1,3)=0.
	trn(2,3)=c*cos_a
	trn(3,3)=c*sin_a
	endif
	end

	subroutine trnslnb0(cell,trn,ks)
chk	===================
c---	ks=1	x//a*,       y//b,        z//(a* x b)
c---	ks=2	x//(b x c*), y//b,        z//c*
c---	ks=3	x//(b* x c), y//b*,       z//c
c---	ks=4	x//a*,       y//(c x a*), z//c
c---	ks=5	x//a,        y//(c* x a), z//c*
c---	ks=6	x//a,        y//b*,       z//(a x b*)
c---	ks=7	x//(a-b),    y//(a+b-2c), z//(a+b+c)    for r+ lattice
c---	ks=8	x//(a-c),    y//(2b-a-c), z//(a+b+c)    for r- lattice

	real trn(3,3),cell(6),a(6),tmp(3,3)
	integer j1(6),j2(6),j3(6)
	data j1/1,3,2,1,3,2/
	data j2/2,2,3,3,1,1/
	data j3/3,1,1,2,2,3/
	logical log
	save

	log=.true.
	goto 100

	entry trnslnb1(cell,trn,ks)
chk	==============
	log=.false.

100	continue
	if(ks.eq.7.or.ks.eq.8) then
	  if(  cell(1).ne.cell(2).or.cell(1).ne.cell(3)
	1  .or.cell(4).ne.cell(5).or.cell(4).ne.cell(6)) then
	   write(6,*) 
	1 'trnslnb1-W> unknown xyz->abc alignment'
  	  endif
	  as= sind(0.5*cell(4))*cell(1)
	  ac= cosd(0.5*cell(4))*cell(1)
	  trn(3,1)= sqrt(ac*ac-as*as/3.)
	  trn(3,2)= trn(3,1)
	  trn(3,3)= trn(3,1)
	  if(ks.eq.8) then
	   trn(1,1)= as 
	   trn(2,1)= -as/sqrt(3.)
	   trn(1,2)= 0.
	   trn(2,2)= -2.*trn(2,1)
	   trn(1,3)= -as
	   trn(2,3)= trn(2,1)
	  else ! if(ks.eq.7)
	   trn(1,1)= as 
	   trn(2,1)= as/sqrt(3.)
	   trn(1,2)=-trn(1,1)
	   trn(2,2)= trn(2,1)
	   trn(1,3)= 0.
	   trn(2,3)= -2.*trn(2,1)
	  endif
	  if(log) then
	  call mxinv(3,trn,ierr)
	  if( ierr.ne.0) 
	1  write(6,*)'trnslnb1-W> ERROR in determining the inverse matrix'
	  endif	
	  return
	endif

	k1=j1(ks)
	k2=j2(ks)
	k3=j3(ks)
	a(1)=cell(k1)
	a(2)=cell(k2)
	a(3)=cell(k3)
	a(4)=cell(k1+3)
	a(5)=cell(k2+3)
	a(6)=cell(k3+3)
	if(     log)	call trnsln0(a,tmp)
	if(.not.log)	call trnsln1(a,tmp)

	trn(k1,k1)=tmp(1,1)
	trn(k1,k2)=tmp(1,2)
	trn(k1,k3)=tmp(1,3)
	trn(k2,k1)=tmp(2,1)
	trn(k2,k2)=tmp(2,2)
	trn(k2,k3)=tmp(2,3)
	trn(k3,k1)=tmp(3,1)
	trn(k3,k2)=tmp(3,2)
	trn(k3,k3)=tmp(3,3)
	end

	subroutine xrot(phi,trn)
chk	===============
	real trn(3,3)
	do 10 i=1,3	
	do 10 j=1,3
10	trn(i,j)=0.
	s=sind(phi)
	c=cosd(phi)	

	trn(1,1)=1.
	trn(2,2)=c
	trn(2,3)=-s
	trn(3,2)=s
	trn(3,3)=c
	end

	subroutine yrot(phi,trn)
chk	===============
	real trn(3,3)
	do 20 i=1,3	
	do 20 j=1,3
20	trn(i,j)=0.
	s=sind(phi)
	c=cosd(phi)	

	trn(1,1)=c
	trn(1,3)=s
	trn(3,1)=-s
	trn(3,3)=c
	trn(2,2)=1.
	end

	subroutine zrot(phi,trn)
chk	===============
	real trn(3,3)
	do 30 i=1,3	
	do 30 j=1,3
30	trn(i,j)=0.
	s=sind(phi)
	c=cosd(phi)

	trn(1,1)=c
	trn(1,2)=-s
	trn(2,1)=s
	trn(2,2)=c
	trn(3,3)=1.
	end

	subroutine axbeqc (a,b,c)
chk	=================
	real a(3,3),b(3,3),c(3,3)
	do 10 i=1,3
	do 10 j=1,3
10	c(i,j)= a(i,1)*b(1,j)+a(i,2)*b(2,j)+a(i,3)*b(3,j)
	end

	subroutine uxveqw (u,v,w)
chk	=================
chk	w = u x v, ie. cross product
	real u(3),v(3),w(3)
	w(1)=u(2)*v(3)-u(3)*v(2)
	w(2)=u(3)*v(1)-u(1)*v(3)
	w(3)=u(1)*v(2)-u(2)*v(1)
	end

 	subroutine read_ar(n,a,*,*)
chk	==================
chk	func: read n real numbers from the txt string
chk	return 1: wrong word
chk	return 2: no word
	include 'edp_main.inc'
!	use edp_dim

	real a(n), t

	parameter (isp=32, itb=9, icm=44)

	character*(max_num_chars) txt
	common /cmm_txt/ n_len,txt,ib,ie

	miss=0
	do k=1,n
	  if(ie.le.0) then
	    ib=1
	  else
	    ib=ie+1
	    endif

	  if(ichar(txt(ib:ib)).eq.icm.and.ib.gt.1) ib=ib+1	! belongs to the previous search
	  do i=ib,n_len
	    j=ichar(txt(i:i))
	    if(   j.ne.isp .and.j.ne.itb) goto 100
	    end do
	  return 2

100	  ib=i
	  do i=ib,n_len
	    j=ichar(txt(i:i))
	    if(j.eq.isp.or.j.eq.icm.or.j.eq.itb) goto 200
	  enddo

200	  ie=i-1
	  if(ib  .gt.ie) then
	    miss=miss+1
	  else
	    read(txt(ib:ie),*,err=50) t
	    a(k)=t
	    endif
	  enddo
	if(miss.gt.0) return 2
	return
50	return 1
	end

	real function atand2(x,y)
chk	====================
	real x,y

	if(y.eq.0.0) then
	  if(x.ge.0.0) then
	    atand2=0.0
	  else
	    atand2=180.0
	    endif
	else if(x.eq.0.0) then
	  if(y.gt.0.0) then
	    atand2=90.0
	  else
	    atand2=270.0
	    endif
	else
	  atand2=atand(y/x)
	  endif
	end

	function deps()
chk	=============
	real*8	deps, eps0, eps1

	logical undone
	data undone/.true./

	data eps1/1.0e-6/
	
	if(undone) then
	  eps0=1.
	  do while (eps0+1.0.gt.1.0)
	    eps1=eps0*2.0
	    eps0=eps0/2.0
	    enddo
	  undone=.false.
	  eps1=eps1*1000.0
	  endif
	deps=eps1
	end

	function eps()
chk	=============
	real eps, eps0, eps1

	logical undone
	data undone/.true./

	data eps1/1.0e-6/
	
	if(undone) then
	  eps0=1.
	  do while (eps0+1.0.gt.1.0)
	    eps1=eps0*2.0
	    eps0=eps0/2.0
	    enddo
	  undone=.false.
	  endif
	eps=eps1
	end

chk***	end of chk_lib.for

copyright by X. Cai Zhang
