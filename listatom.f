
	subroutine correlatn(*)
chk	====================
	include 'edp_main.inc'
!	use edp_main
	include 'edp_dat.inc'
!	use edp_dat

	character*1 s(5)
	data s/'x','y','z','w','b'/
	integer col1, col2, col3

	external match_l

	logical put_x

	n_of_syn=2					!000515
	syntax(1)='syntax:' 
	syntax(2)=
	1'correlation group_id.s (x,y,z,w,b) (x,y,z,w,b) [(x,y,z,w,b)]' 

	call find_a_group(igr)

	if( igr.lt.0 .or. igr.eq.999 ) then
	  errmsg=' errmsg: a group name needed'
	  return 1
	  endif

        n_gr=n_groupa(igr)

	n_on=0
	do i=1,n_atom
	  if(lf(i)) n_on=n_on+1
	  enddo

	if(verbose.ge.2 ) write(6,1009)   n_on,cgroup(igr),n_gr
1009	format(' corr-I2>  #of on atoms (n_on)=',i5,
	1 ', #of atoms in group ',a,' (n_gr)=',i5)

	if(n_on.ne.n_gr .or. n_on.eq.0) then
	  errmsg=' errmsg: n_on =/= n_gr or n_on=0'
	  return 1
	  endif

	col1= match_l(5,s)
	if(col1.le.0) return 1

	col2= match_l(5,s)
	if(col2.le.0) return 1

	col3= match_l(5,s)
	if(col3.lt.0) return 1

	sw=0.
	sb=0.
	swb=0.
	sww=0.
	sbb=0.
	jj=0
	kk=0
	do i=1,n_atom
	  if(lf(i)) then
           jj=jj+1
           j=igroupa( jj,igr)
	   vi=get_x(col1,i)
	   vj=get_x(col2,j)
	   if(vi.ne.999.0.and.vj.ne.999.0) then
	    kk=kk+1
	    sw=sw+vi
	    sb=sb+vj
	    swb=swb+vi*vj
	    sww=sww+vi*vi
	    sbb=sbb+vj*vj
	    endif
	   endif
	  enddo

	if(kk.le.2) then
c	  write(6,*) '  corr=',0.0
	  errmsg=' errmsg: too few valid input records'
	  return 1
	  endif

	sn=float(kk)
c	sn=float(n_on)
	sww_swsw=sww-sw*sw/sn
	sbb_sbsb=sbb-sb*sb/sn

	if(abs(sbb_sbsb).le.1.e-6) then
c	  write(6,*) '  corr=',0.0
	  errmsg=' errmsg: colume '//s(col2)//' is featureless.'
	  return 
	  endif

	if(abs(sww_swsw).le.1.e-6)  then
c	  write(6,*) '  corr=',0.0
	  errmsg=' errmsg: colume '//s(col1)//' is featureless.'
	  return 1
	  endif

	write(6,*) '  corr=',(swb-sw*sb/sn)/sqrt(sww_swsw*sbb_sbsb)

	aa= (swb-sw*sb/sn) / sbb_sbsb
	bb= (sw*sbb-swb*sb)/ (sbb_sbsb*sn)

	delta=0.
	jj=0
	do i=1,n_atom
	  if(lf(i)) then
            jj=jj+1
            j=igroupa( jj,igr)
	    d=get_x(col1,i) - (aa*get_x(col2,j)+bb)
	    delta=delta+d*d
	    if(col3.gt.0) then
	      if(put_x(col3,i,d)) then	! put_x return .false. if okey
	        errmsg=' errmsg: ouput uncomplete'
	        return 1
	        endif	
	      endif	
	    endif	
	  enddo
1061	format(f8.3)

	write(6,*) '  s.d.=',sqrt(delta/sn)
	write(6,*) ' ',s(col1),'(on)=',aa,'*',s(col2),'(gr) +',bb
	end

	real function get_x(col,i)
chk	===================
chk	get value in col:col of atom:i 
chk	col range 1-5, correponding to xyzwb

	include 'edp_main.inc'
!	use edp_main
	integer col

	if(col.eq.1) then
	  read(text(i)(31:38),1061,err=900) get_x
	else if(col.eq.2) then
	  read(text(i)(39:46),1061,err=900) get_x
	else if(col.eq.3) then
	  read(text(i)(47:54),1061,err=900) get_x
	else if(col.eq.4) then
	  get_x=w(i)
	else 
	  get_x=b(i)
	  endif
1061	format(f8.3)
	return
900	write(6,*) 'get_x-W> ERROR: xyz are not numbers.'
	get_x=0.0
	end

	logical function put_x(col,i,get_x)
chk	======================
chk	put value in col:col of atom:i 
chk	col range 1-5, correponding to xyzwb
chk	put_x=.fasle. -- successful output
chk	put_x=.true.  -- unsuccessful output

	include 'edp_main.inc'
!	use edp_main
	integer col

	if(col.eq.1) then
	  write(text(i)(31:38),1061,err=900) get_x
	else if(col.eq.2) then
	  write(text(i)(39:46),1061,err=900) get_x
	else if(col.eq.3) then
	  write(text(i)(47:54),1061,err=900) get_x
	else if(col.eq.4) then
	  w(i)=get_x
	else 
	  b(i)=get_x
	  endif
1061	format(f8.3)
	put_x=.false.
	return
900	put_x=.true.
	end

	subroutine jiggle(*)
chk	=================
chk	this idea is stolen from jiggle command of convert of tnt package.
	include 'edp_main.inc'
!	use edp_main

	character*1 s(5)
	data s/'x','y','z','w','b'/

	data iseed/1/

	n_of_syn=2					!000515
	syntax(1)='syntax:' 
	syntax(2)='jiggle (x, y, z, w, b) limit.r [shift.r] '

	lgo=match_l(5,s)
	if(lgo.le.0) return 1

	a1=0.
	call read_ar(1,	a0, *900, *900)
	call read_ar(1,	a1, *900, *800)
800	amin=a1-a0
	amax=a1+a0
	d=amax-amin

	if(lgo.eq.1) then
	  do i=1,n_atom
	    if(lf(i)) then
	      x(i)=x(i)+ran(iseed)*d+amin
	      write( text(i)(31:38),1061,err=910) x(i)
1061	      format(f8.3)
	      end if
	    enddo

	else  if(lgo.eq.2) then
	  do i=1,n_atom
	    if(lf(i)) then
	      y(i)=y(i)+ran(iseed)*d+amin
	      write( text(i)(39:46), 1061,err=910) y(i)
	      end if
	    enddo
	else if(lgo.eq.3) then
	  do i=1,n_atom
	    if(lf(i)) then
	      z(i)=z(i)+ran(iseed)*d+amin
	      write( text(i)(47:54),1061,err=910) z(i)
	      end if
	    enddo
	else if(lgo.eq.4) then
	  do i=1,n_atom
	    if(lf(i)) w(i)=w(i)+ran(iseed)*d+amin
	    enddo
	else
	  do i=1,n_atom
	    if(lf(i)) b(i)=b(i)+ran(iseed)*d+amin
	    enddo
	  endif
	return
900	return 1
910	write(6,*) 'jiggle-W> ABORTED: too large shift, error in output.'
	end

copyright by X. Cai Zhang

	subroutine listzone()
chk	===================
	include 'edp_main.inc'
!	use edp_main
	include 'edp_dat.inc'
!	use edp_dat

	logical first

	character*1 xcode(-1:2)
	data xcode/'R',' ','L','?'/

	parameter (max_xt=100)
!	integer, parameter :: max_xt=100
	real vc(6,max_lv), xt(4,max_xt)

	character*(32) txtv(max_lv)
	logical t(max_lv,max_lv)
	real    a(3,max_lv,max_lv), d3(3)
	integer mkv(max_lv), nkv(max_lv), jkv(max_lv), lkv(max_lv)
	1 ,kkv(3,max_lv), ikv(2,max_lv)
	character*(120) fmtstmt
!	save first
	save	!030422

	ir0=-1
	n=0
	sx=0.0
	sy=0.0
	sz=0.0
	sw=0.0
	sb=0.0
	first=.true.
	if(verbose.lt.1) then
	 io=48
	else
	 io=6
	endif
	
	do ir=1,n_res
	  do i=ijk(ir),ijk(ir+1)-1
            if(lf(i)) then
	      n=n+1
	      sx=sx+x(i)
	      sy=sy+y(i)
	      sz=sz+z(i)
	      sw=sw+w(i)
	      sb=sb+b(i)
	      if(ir.eq.ir0) then
	         !do nothing
	      else if(ir.eq.ir0+1) then
	        ir0=ir
	      else 
	        if(ir0.gt.0) then
	          if(first) write(io,1002)
	          first=.false.
	          n1=n-1
		  write(io,1001) 
	1 ires(ir00),ires(ir0),n1
	1 ,(sx-x(i))/n1, (sy-y(i))/n1, (sz-z(i))/n1
	1 ,(sw-w(i))/n1, (sb-b(i))/n1
	          endif
	        ir00=ir
	        ir0=ir
	        n=1
	        sx= x(i)
	        sy= y(i)
	        sz= z(i)
	        sw= w(i)
	        sb= b(i)
	        endif
	      endif
	    enddo
	  enddo

	if(ir0.gt.0) then
	  if(first) write(io,1002)
	  write(io,1001)
	1 ires(ir00),ires(ir0),n ,sx/n ,sy/n ,sz/n ,sw/n ,sb/n
	endif

	if(io.eq.48) call typelist(io)
1002	format(10x,'zone',t22,' #atom',t39,'center',t56,'<w>',t64,'<b>')
1001	format(' zone ',a5,'- ',a5,'! [',i6,']',5f8.3)

	return 
	
	entry listvector(*)
chk	================
	if(verbose.ge.6 ) write(6,1069) 	!000504
1069	format(' listvector-I6> reference:'
	1/'  Koch, I., F. Kaden, and J. Selbig (1992)'
	1/'  Analysis of protein topolgies by graph theoretical methods.'
	1/'  proteins: Struc. Func. and Genetics 12:314-323.')

	d3(1)=3.0
	d3(2)=6.5
	d3(3)=11.0
	call read_ar(3,d3,*900,*800)
	
800	ir0=-1
	n=0
	kv=1
	do ir=1,n_res
	  do i=ijk(ir),ijk(ir+1)-1
            if(lf(i)) then
	      n=n+1
	      if(n.gt.max_xt) then
	        errmsg=
	1' errmsg: too many atoms in one vector. (increase max_xt)'
		if(verbose.ge.3 ) write(6,'(a,i)') 
	1 ' listvector-I3> current max_xt=',max_xt
	        return 1
	        endif
	      xt(1,n)=x(i)
	      xt(2,n)=y(i)
	      xt(3,n)=z(i)
	      xt(4,n)=w(i)
	      if(ir.eq.ir0) then
	         !do nothing
	      else if(ir.eq.ir0+1) then
	        ir0=ir
	      else 
	        if(ir0.gt.0) then
	          n1=n-1
	          write(txtv(kv),1001) ires(ir00),ires(ir0),n1
	          call find_the_long_axis2(n1,xt,vc(1,kv),*900)
	          ikv(1,kv)=ir00
	          ikv(2,kv)=ir0
		  kv=kv+1
	          if(kv.gt.max_lv) then
	        errmsg=
	1 ' errmsg: too many vectors. (increase max_lv)'
		if(verbose.ge.3 ) write(6,'(a,i)') 
	1 ' listvector-I3> current max_lv=',max_lv
	            return 1
	            endif
	          endif
	        ir00=ir
	        ir0=ir
	        n=1
	        xt(1,n)=x(i)
	        xt(2,n)=y(i)
	        xt(3,n)=z(i)
	        xt(4,n)=w(i)
	        endif
	      endif
	    enddo
	  enddo

	if(ir0.gt.0) then
	  n1=n
	  write(txtv(kv),1001) ires(ir00),ires(ir0),n1
	  call find_the_long_axis2(n1,xt,vc(1,kv),*900)
	  ikv(1,kv)=ir00
	  ikv(2,kv)=ir0
	  endif

	do k=1,kv
	  do m=k+1, kv
	    call v_nayb(k,m,vc,t,a,d3)
	    enddo
	  t(k,k)=.false.
	  enddo
	
	first=.true.
	kcycle=0
40	do k=1,kv
	  mm=0		! mm: # of neighbors of vector k
	  do m=1,kv
	    if(t(k,m)) mm=mm+1
	    enddo
	  if(mm.gt.0) then
	   if(mm.eq.1.or.k.ge.kv-kcycle) then	
c	mm=1, tow strands 
c	mm>1, k.ge.kv too many vectors -> circular.
	    k1=k		! vector k is selected as k1 (starting strand)
	    kcycle=kcycle+1
	    goto 45
	    endif
	   endif
	  enddo
c	call typelist(48)
	return

45	ki=0
	kj=0
50	do k=1,kv
	  if(t(k,k1)) then
	    t(k,k1)=.false.
	    t(k1,k)=.false.
	    if(first) write(6,1104)
	    first=.false.
c	list the first close neighbor of vector k1
	    write(6,1101) txtv(k1)(:28), (a(m,k1,k),m=1,3)
	    ki=ki+1
	    mkv(ki)=k1*kv+(kv-ki)		! order of the vectors
	    jkv(ki)=kj			! direction of the vectors
	    if(a(1,k1,k).gt.90.0) kj=kj+1
1101	format(a28,2x,3f8.3)
1103	format(a28)
1104	format(10x,'zone',t22,' #atom',t34,'angle',t42,'dist.'
	1 ,t50,'center_dist.')
	    k1=k			! next vector
	    goto 50
	    endif
	  enddo
	if(first) write(6,1104)
	write(6,1103) txtv(k1)
	ki=ki+1
	mkv(ki)=k1*kv+(kv-ki)		! order of the vectors
	jkv(ki)=kj			! direction of the vectors

	if(ki.gt.1) then
	  call sort_in2(ki, mkv, nkv)
	  ki1=ki
	  ki=ki-1
	  do k=1,ki
	    k1=nkv(k)
	    k2=nkv(k+1)
	    lkvk=mod(jkv(k2)+jkv(k1)+1,2)
	    k1=mkv(k1)/kv
	    k2=mkv(k2)/kv
	    if(lkvk.ne.0) then ! parrellel, therefore cross
	      lkv(k)= l_or_r_cross(k1,k2, vc, ikv)
	    else ! antiparrallel
	      lkv(k)=0
	      endif
	    if(k1.eq.k2) then ! circular vector barrel
	      kkv(1,k)=0
	      kkv(2,k)=0
	    else
	      kkv(1,k)=nkv(k+1)-nkv(k)
	      kkv(2,k)=ikv(1,k2) -ikv(2,k1) -1
	      endif
	    kkv(3,k+1)=ikv(2,k2) -ikv(1,k2) +1
	    kkv(3,k  )=ikv(2,k1) -ikv(1,k1) +1
	    enddo
!1102	format(/' strand:'   ,<ki1>i4
!	1      /'  key  [   ',<ki>(i3,a1),']'
!	1      /' code  [   ',<ki>(i3,a1),']'
!	1      /' loop  :  ' ,<ki>i4/)
!              123456789012345678901234567890
	fmtstmt='(/" strand:"   ,<n>i4         '//
	1        ' /"  key  [   ",<n>(i3,a1),"]"'//
	1        ' /" code  [   ",<n>(i3,a1),"]"'//
	1        ' /" loop  :  " ,<n>i4/)'
	write(fmtstmt( 17: 19), '(i3)') ki1
	write(fmtstmt( 47: 49), '(i3)') ki 
	write(fmtstmt( 77: 79), '(i3)') ki 
	write(fmtstmt(107:109), '(i3)') ki 

	  write(6,fmtstmt) 
	1  ( kkv(3,k),k=1,ki1),
	1  ( kkv(1,k),xcode(lkv(k)),k=1,ki),
	1  (-kkv(1,k),xcode(lkv(k)),k=1,ki),
	1  ( kkv(2,k),k=1,ki)
	  endif
	goto 40
900	return 1
	end

	function l_or_r_cross(k1,k2,vc,ikv)
chk	=====================
	include 'edp_main.inc'
!	use edp_main
	real vc(6,max_lv)
	integer ikv(2,max_lv)
	real v1(3), v2(3), v3(3)

	v1(1)=vc(4,k1)
	v1(2)=vc(5,k1)
	v1(3)=vc(6,k1)

	v2(1)=vc(1,k2)-vc(1,k1)
	v2(2)=vc(2,k2)-vc(2,k1)
	v2(3)=vc(3,k2)-vc(3,k1)

	ir0=ikv(2,k1)+1
	ir1=ikv(1,k2)-1

	n=0
	xc=0.0
	yc=0.0
	zc=0.0
	wc=0.0
	do ir=ir0,ir1
	  do i=ijk(ir), ijk(ir+1)-1
	    wi=w(i)
	    wc=wc+wi
	    xc=xc+x(i)*wi
	    yc=yc+y(i)*wi
	    zc=zc+z(i)*wi
	    enddo
	  enddo

	if(wc .le. eps() ) then
	  l_or_r_cross=2
	  return
	  endif

	v3(1)=xc/wc - vc(1,k1) 	  
	v3(2)=yc/wc - vc(2,k1) 	  
	v3(3)=zc/wc - vc(3,k1)

	d= v1(1)*(v2(2)*v3(3)-v2(3)*v3(2))
	1 +v1(2)*(v2(3)*v3(1)-v2(1)*v3(3))
	1 +v1(3)*(v2(1)*v3(2)-v2(2)*v3(1))
	
	if(d.gt.0.0) then
	  l_or_r_cross=1
	else
	  l_or_r_cross=-1
	  endif
	  
	end

	subroutine v_nayb(iv1,iv2,vectors,t,a,d3)
chk	=================
	include 'edp_dim.inc'
!	use edp_dim

	real vectors(6,max_lv)
	logical t(max_lv,max_lv)
	real    a(3,max_lv,max_lv), d3(3), f(3)

	theta= acosd(vectors(4,iv1)*vectors(4,iv2)+
	1            vectors(5,iv1)*vectors(5,iv2)+
	1            vectors(6,iv1)*vectors(6,iv2))

	f(1)=vectors(1,iv2)-vectors(1,iv1)
	f(2)=vectors(2,iv2)-vectors(2,iv1)
	f(3)=vectors(3,iv2)-vectors(3,iv1)
	ff=sqrt(f(1)*f(1)+f(2)*f(2)+f(3)*f(3))
	if(ff.gt. eps() ) then
	  f(1)=f(1)/ff
	  f(2)=f(2)/ff
	  f(3)=f(3)/ff
	else
	  f(1)=0.0
	  f(2)=0.0
	  f(3)=1.0
	  endif

	call uxveqw(vectors(4,iv1),vectors(4,iv2),f)
	aa= sqrt( f(1)*f(1)+f(2)*f(2)+f(3)*f(3))
	if( abs( aa) .le. eps() ) then
	  t(iv1,iv2)=.false.
	  t(iv2,iv1)=.false.
	  return
	else
	  f(1)=f(1)/aa
	  f(2)=f(2)/aa
	  f(3)=f(3)/aa
	  endif

	dist=abs( 
	1     f(1)*(vectors(1,iv2)-vectors(1,iv1))+
	1     f(2)*(vectors(2,iv2)-vectors(2,iv1))+
	1     f(3)*(vectors(3,iv2)-vectors(3,iv1)))

	dist1=sqrt( 
	1     (vectors(1,iv2)-vectors(1,iv1))**2+
	1     (vectors(2,iv2)-vectors(2,iv1))**2+
	1     (vectors(3,iv2)-vectors(3,iv1))**2)

	if( d3(1) .le. dist .and. dist .le. d3(2)
	1 .and. dist1 .le. d3(3)) then
	  t(iv1,iv2)=.true.
	  t(iv2,iv1)=.true.
	  a(1,iv1,iv2)= theta
	  a(1,iv2,iv1)= theta
	  a(2,iv1,iv2)= dist
	  a(2,iv2,iv1)= dist
	  a(3,iv1,iv2)= dist1
	  a(3,iv2,iv1)= dist1
	else
	  t(iv1,iv2)=.false.
	  t(iv2,iv1)=.false.
	  endif
	end
 
	subroutine list(*)
chk	===============
	include 'edp_main.inc'
!	use edp_main
	include 'edp_dat.inc'
!	use edp_dat

	common /cmm_txt/ n_len, txt, ib ,ie
	character*(max_num_chars) txt
	logical nword
!	save i0, i1, i2, k, l, i3
	save	!030422

	n_of_syn=4					!000515
	syntax(1)='syntax:' 
	syntax(2)='1) list [ line1.i [line2.i] ]'
	syntax(3)='2) list zone' 
	syntax(4)=
	1'3) list vector [dmin.r, dmax.r, max_center_to_center_dist.r]'

	i1=1
	i2=n_atom
	i3=0

	if(nword(n_len,txt,ib,ie)) then
	  do i=1,n_atom
	    if(lf(i)) i3=i3+1
	    end do
	else if(txt(ib:ie).eq.'zone') then
	  call listzone()
	  return
	else if(txt(ib:ie).eq.'vector') then
	  call listvector(*990)
	  return
	else
	  read( txt(ib:ie),*,err=990) i1
	  if(.not.nword(n_len,txt,ib,ie)) 
	1   read( txt(ib:ie),*, err=990) i2
	  endif

	if(verbose.le.0 ) then
	  write(6,*) 'list-W> UNDONE: verbose <= 0.'
	  return
	  end if


	k=0
	l=1
	goto 20

	entry c_list
!	============
!	if(.not. inter .or.i0.le.0) return
	if(verbose.le.0 .or.i0.le.0) return
	goto 30
20	i0=1
30	if(window_size.le.1) call get_window_size()
	i00=i0
	do ii=i00,n_atom,1
	  i=iorder(ii)
	  if(lf(i)) then
	    k=k+1
	    if(k.gt.i2) return
	    if(k.ge.i1) then
	      if( b(i) .gt. -99.99 .and. b(i) .lt. 999.99 .and.
	1         w(i) .gt. -99.99 .and. w(i) .lt. 999.99 ) then
		write(text(i)(55:66),1011) w(i),b(i)
1011	format(2f6.2)		
	        write(6,1001) text(i)
1001	format(1x,a72)
	      else
	        write(6,'(a)') ' list-W> ERROR in typing the occ/b'
	        end if
	      l=l+1
	      if(mod(l,window_size).eq.0) then
	        if(i3.le.0) then
	          write(6,'(a)') ' list> return for more .....'
	        else if(k.lt.i3) then
	          write(6,1003) (100.*k)/i3
1003	          format(' list>',f5.1,'%, return for more ....')
	        else
	          goto 980 
	          end if
	        i0=ii+1
	        return
	      end if
	    end if
	  end if
	end do
	
	entry e_list
c	============
980	i0=-1
	return

990	return 1
	end 

	subroutine zone(*)
c	=======================
	include 'edp_main.inc'
!	use edp_main
	integer  iseg(max_num_zones),k_sum(max_num_zones)
	character*1 cm

	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie
	logical nword, first, lgo
	dimension ii1(max_num_zones),ii2(max_num_zones)
	logical lf1(max_atom)
!	save first, num_zones, iseg
	save	!030422
	
	n_of_syn=5					!000515
	syntax(1)='syntax:' 
	syntax(2)='1) zone' 
	syntax(3)='2) zone all' 
	syntax(4)='3) zone res_id1.s-res_id2.s' 
	syntax(5)='4) zone res_id1.s [res_id2.s ...]' 

	if(zone_untouch) then
	  first=.true.
	  zone_untouch=.false.
	  endif

	lgo=.false.
	if(.not.nword(n_len,txt,ib,ie)) goto  940

chk	for command	"zone <cr>"
10	if(first) then
	  first=.false.
	  i0=0
	  num_zones=1
	  cm='*'
	  do 920 i=1,n_res
	    j=ijk(i)
	    ie0=22
	    if( nword(72,text(j),ib0,ie0)) goto 999
	    read( text(j)(ib0:ie0),*,err=910) i1
	    if(i1-i0.eq.1.and.text(j)(22:22).eq.cm) goto 920
910	    iseg(num_zones)=i
	    cm=text(j)(22:22)
	    num_zones=num_zones+1
	    if(num_zones.ge.max_num_zones) goto 930
920	    i0=i1
930	  iseg(num_zones)=n_res+1
	  num_zones=num_zones-1
	  endif
	if(lgo) return

	do ii=1,num_zones
	  ksum=0
	  i0=iseg(ii)
	  i1=iseg(ii+1)-1
	  do i=i0,i1
	  do j=ijk(i),ijk(i+1)-1
	    if(lf(j)) ksum=ksum+1
	    end do
	    end do
	  k_sum(ii)=ksum
	  end do
	write(6,1090) 	
	1 (ires(iseg(ii)),ires(iseg(ii+1)-1),k_sum(ii),ii=1,num_zones)
c1090	format(' listzone> ',(t13,3(a5,' - ',a5,'[',i5,'] ')))
1090	format(' zone> ',(t13,3(a5,' - ',a5,'[',i5,'] ')))
	return
999	write(6,*) '%EdPDB-F- unexpected error in zone'
	write(6,*)
	call exit( 4)

940	do i=1,n_atom
	  lf1(i)=.false.
	  enddo

	ib0=ib
	ie0=ie
	call dfgroup(igr,*900)
	j1=n_groupa(igr)
	do j=1,j1
	  i=igroupa( j,igr)
	  lf1(i)=.true.
	  enddo
	ib=ib0
	ie=ie0

	if(txt(ib:ie).eq.'all') then
command	"zone all [from ...]"
	  k=1
	  ii1(1)=1
	  ii2(1)=n_atom
	else
command	"zone .... [from ...]"
	  ie=ib-1
	  call region (k,ii1,ii2,*900)
command	"zone from [...]"
	  if(k.eq.0) then
	    k=1
	    ii1(1)=1
	    ii2(1)=n_atom
	    endif
	  endif

	do j=1,k
	do i=ii1(j),ii2(j)
	 if(lf1(i)) lf(i)=incl
	end do
	end do
	return
c---
900	return 1

	entry code(*)
c	=============
	call find(n_res,ires,i)
	if(i.le.0) return 1
	na = ijk(i+1)-ijk(i)
	n1 = ijk(i)
	call find(na,atom(n1),j)
	if(j.lt.0) return 1
	if(j.eq.0) then
	write(6,1009) ires(i),res(i),i
1009	format(' code> ',2a5,' :',i5: ' , ',a5,' :',i5,' (',i4,')'
	1/t8,' current xyz =',3f8.3)
	else
	j=j+ijk(i)-1
	write(6,1009) ires(i),res(i),i,atom(j),j,iorder(j),x(j),y(j),z(j)
	end if
	
	entry beforecm
c	==============
	if(first.or.zone_untouch) then
	  first=.true.
	  zone_untouch=.false.
	  lgo=.true.
	  goto 10
	  end if
	end

copyright by X. Cai Zhang
	subroutine listatom
chk	===================
	include 'edp_main.inc'
!	use edp_main
	character*3 res_l(max_res),atom_l(max_res)*4
	character*3 resi,atomi*4
	dimension kres(max_res), katom(max_res)
	logical lr,la
	data lr,la/.true.,.true./

	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie
	save	!030422
	
	if(la) then
	la=.false.
	atom_l(1)=atom(1)
	k_atom=1
	n_len=4
	do i=2,n_atom
		txt=atom(i)
		ie=0
		call find(k_atom,atom_l,j)
		if(j.le.0) then
			k_atom=k_atom+1
			atom_l(k_atom)=atom(i)
			if(k_atom.ge.max_res) goto 800
		end if
	end do
	end if

800	do j=1,k_atom
	katom(j)=0
	end do

	do 700 i=1,n_atom
		if(lf(i)) then
		atomi=atom(i)
		do j=1,k_atom
			if(atomi.eq.atom_l(j)) then
				katom(j)=katom(j)+1
				goto 700
			end if
		end do
		write(6,*) '%EdPDB-F- unexpected error in listatom'
		write(6,*)
		call exit( 4)
		end if
700	continue

	write(6,1001) (atom_l(i),min(katom(i),999),i=1,k_atom)
1001	format(' atom> ',(t10,6(a5,'[',i3,']')))
	if(k_atom.eq.max_res) write(6,1002) max_res
1002	format(' atom-W> only',i4,' atom_names have been listed out')
	return

	entry listres
c	=============
	if(lr) then
	lr=.false.
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
			if(k_res.ge.max_res) goto 850
		end if
	end do
	end if

850	do j=1,k_res
	kres(j)=0
	end do

	do 750 ii=1,n_res
	do 740 i=ijk(ii),ijk(ii+1)-1
		if(lf(i)) then
		resi=res(ii)
		do j=1,k_res
			if(resi.eq.res_l(j)) then
				kres(j)=kres(j)+1
				goto 750
			end if
		end do
		write(6,*) '%EdPDB-F- unexpected error in listres'
		call exit( 4) 
		end if
740	continue
750	continue
	write(6,1004) (res_l(i),min(kres(i),999),i=1,k_res)
1004	format(' residue> ',(t10,6(a4,' [',i3,']')))
	end
chk***	the end of listatom.for

copyright by X. Cai Zhang
