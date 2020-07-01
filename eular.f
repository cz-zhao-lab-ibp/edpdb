
	subroutine edp_rms(*)
chk	==================
c
c	overlay the on atoms to the given group.
c	write out a rtn_.txt file.
c
	include 'edp_main.inc'
! 	use edp_main
	include 'edp_file.inc'
!	use edp_file
	include 'edp_dat.inc'
!	use edp_dat

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
	
	logical  open_file
	external open_file
	logical  open_file1
	external open_file1

	is_angle=0
	
	n_of_syn=2					!000515
	syntax(1)='syntax:' 
	syntax(2)='overlay group_id.s [filename.s] [weight] '

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

	do i=1,n_atom
	  if(lf(i)) then
	    jj=jj+1
	    j=igroupa( jj,igr)
	    if(i.ne.j) identical=.false.
	    s(1,jj)=x(i)
	    s(2,jj)=y(i)
	    s(3,jj)=z(i)
	    p(1,jj)=x(j)
	    p(2,jj)=y(j)
	    p(3,jj)=z(j)
	    if(even_weight) then 
	      weight(jj)=1.
	    else
	      weight(jj)=w(i)
	      endif
	    do k=1,3
	      sumx(k) =sumx(k) + (p(k,jj)-s(k,jj))    *weight(jj)
	      sumxx(k)=sumxx(k)+ (p(k,jj)-s(k,jj))**2 *weight(jj)
	      enddo
	    sumw=sumw+ weight(jj)
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
	1  ,sumx(1),sumx(2),sumx(3) 
1018	  format(3(3f12.7/),3f12.5)
	  if(verbose.ge.3) write(6,1230) rtn_out(:ltrim(rtn_out))
1230	format(' overlay-I3> To apply the transformation, type:'
	1 /'   rtn file ',a)

	  rms=sqrt(max(0.0,sgm))
	  shift(1)=sqrt(sumxx(1)+sumxx(2)+sumxx(3))
	  shift(2)=0.
	else
	  rms=-1.0
	  call edp_rms_doit (rms, n_group0,s,p,weight, shift, 29
	1  ,.false., .false.)
	  end if

	close (29)
	if(rms.ge.0.) write(6,1010) rms, shift(1),shift(2)
1010	format(
	1' rms> rms=',e11.4,', shifted  by',e11.4,', rotated by',e11.4)

	if(rtn_out.ne.'rtn_.txt') 
	1 write(6,1011) rtn_out(:ltrim(rtn_out))
1011	format(' rms>  the file ',a,' is created or overwritten.')
	end

	subroutine edp_rms_doit(rms,n,s,p,w, trn, io, rms_only, need_matrix)
c	===========================
c	n -- # of atom pairs
c	s -- array of molecule 1
c	p -- array of molecule 2
c	rms (input) < 0.0 : do quasi mirror
c	rms (output)      : return rms coordinate difference
	implicit real*8  (a-h)		!????
	include 'edp_dat.inc'
!	use edp_dat

!	implicit real*8  (a-h)

	logical rms_only, need_vect, need_matrix

	dimension s(3,n),p(3,n), w(n)
	1	,sc(3),pc(3), f1(3),f2(3),f3(3)
	1	,a(3,3),b(3,3), e(3), v(3,3)

	real rms, s,p, w, real_d, trn(3,4)
	common /cmm_mcs/ is_angle
	real acosd

	if(verbose.ge.6 ) write(6,1069)	!000504
1069	format(' edp_rms-I6> reference:'
	1/'  McLachlan AD. (1979)'
	1/'  Gene duplications in the structural evolution of chymotrypsin.'
	1/'  J Mol Biol., 128(1):49-79.')

	do i=1,3
	  sc(i)=0.
	  pc(i)=0.
	  do j=1,3
	    b(j,i)=0.0
	    enddo
	  enddo

	wt=0.
	do i=1,n
	  wt=wt+w(i)
	  enddo
	ratio=100./wt	! rescale the total weights to 100.

	wt=0.
	do i=1,n
	  w_ratio=w(i)*ratio
	  wt=wt+ w_ratio
	  if(is_angle.eq.0) then
	   do j=1,3
	    sc(j)=sc(j)+s(j,i)*w_ratio
	    pc(j)=pc(j)+p(j,i)*w_ratio
	   enddo
	  endif  
	 enddo
	rn=1./wt

c	ratio1=rn*n
c	ratio2=ratio1*ratio1
c	ratio4=ratio2*ratio2
c	ratio6=ratio4*ratio2

	if(is_angle.eq.0) then
	 do j=1,3
	  scj=sc(j)*rn
	  pcj=pc(j)*rn
	  sc(j)=scj
	  pc(j)=pcj
	  do i=1,n
	    s(j,i)=s(j,i)-scj
	    p(j,i)=p(j,i)-pcj
	  enddo
	 enddo
	endif

	trn(1,1)= sqrt( (sc(1)-pc(1))**2+(sc(2)-pc(2))**2+(sc(3)-pc(3))**2)

	ss=0.
	pp=0.
	do i=1,n
	  w_ratio=w(i)*ratio
	  do j=1,3
	    ss=ss+s(j,i)**2*w_ratio
	    pp=pp+p(j,i)**2*w_ratio
	    do k=1,3
	      b(j,k)=b(j,k)+s(j,i)*p(k,i)*w_ratio
	      enddo
	    enddo
	  enddo

	do i=1,3
	do j=i,3
	  aij=0.0
	  do k=1,3
	    aij=aij+b(j,k)*b(i,k)
	    enddo
	  a(j,i)=aij
	  a(i,j)=aij
	  enddo
	  enddo

	d= b(1,1)*(b(2,2)*b(3,3)-b(2,3)*b(3,2))
	1 +b(1,2)*(b(2,3)*b(3,1)-b(2,1)*b(3,3))
	2 +b(1,3)*(b(2,1)*b(3,2)-b(2,2)*b(3,1))

	if(d.lt.0.) then
	  if(rms.ge.0.and.n.gt.3) then	! test! for the mcs_atm subroutine
	    rms=-1.
	    return
	    endif
	  if(verbose.ge.3 ) then
	    write(6,*) 'rms-W3> quasi mirror symetry may exist.'
	    write(6,*) 
	1 '     only will pure rotation (plus translation) be given.'
	    endif
	  endif

c	rms=-1. :  edp_rms_doit aborted
	rms=-1.

	need_vect=.not.rms_only
	call eigen(a,e,v, need_vect, *900)

	e1=sqrt(abs(e(1)))
	e2=sqrt(abs(e(2)))
	e3=sqrt(abs(e(3)))

	if(d.lt.0.) then
	  rms=(ss+pp-(e1+e2-e3)*2.0)*rn
	else
	  rms=(ss+pp-(e1+e2+e3)*2.0)*rn
	  endif
c	rms=(ss+pp-(e1+e2+e3)  *2.0)*rn
	rms= sqrt(abs(rms))
!	convert the result to the degree unit, if angle rms is used	
	if(is_angle.eq.1) rms=2.0*asind(0.5*rms)
	if(rms_only) return

	if(e2.le.deps()) then
	  write(6,*) '%EdPDB-W-  the matrix can not be determine.'
	  return
	  endif

	do j=1,3
	  fj1=0.0
	  fj2=0.0
	  do i=1,3
	    fj1=fj1+b(i,j)*v(1,i)
	    fj2=fj2+b(i,j)*v(2,i)
	    enddo
	  f1(j)=fj1/e1
	  f2(j)=fj2/e2
	  enddo

	call ixjeqk(f1,f2,f3)

	do 62 i=1,3
	do 62 j=1,3
62	a(j,i)=f1(j)*v(1,i)+f2(j)*v(2,i)+f3(j)*v(3,i)

240	tr1=	pc(1)-a(1,1)*sc(1)-a(1,2)*sc(2)-a(1,3)*sc(3)
	tr2=	pc(2)-a(2,1)*sc(1)-a(2,2)*sc(2)-a(2,3)*sc(3)
	tr3=	pc(3)-a(3,1)*sc(1)-a(3,2)*sc(2)-a(3,3)*sc(3)	

249	if(io.gt.0) then
	  write(io,118) a(1,1),a(1,2),a(1,3)
	2	       ,a(2,1),a(2,2),a(2,3)
	3	       ,a(3,1),a(3,2),a(3,3),tr1,tr2,tr3
	  write(io,*)'  inverse translation:'
	  write(io,118) a(1,1),a(2,1),a(3,1),
	1	        a(1,2),a(2,2),a(3,2),
	2	        a(1,3),a(2,3),a(3,3),
	3   -a(1,1)*tr1-a(2,1)*tr2-a(3,1)*tr3,
	4   -a(1,2)*tr1-a(2,2)*tr2-a(3,2)*tr3,
	5   -a(1,3)*tr1-a(2,3)*tr2-a(3,3)*tr3
118	  format(3(3f12.7/),3f12.5)
	else if(need_matrix) then
	  trn(1,1)=a(1,1)	
	  trn(1,2)=a(1,2)	
	  trn(1,3)=a(1,3)	
	  trn(1,4)=tr1
	  trn(2,1)=a(2,1)	
	  trn(2,2)=a(2,2)	
	  trn(2,3)=a(2,3)	
	  trn(2,4)=tr2
	  trn(3,1)=a(3,1)	
	  trn(3,2)=a(3,2)	
	  trn(3,3)=a(3,3)
	  trn(3,4)=tr3
	  return
	  endif

	d=0.5*(a(1,1)+a(2,2)+a(3,3)-1.)
	if(abs(d).gt.0.99999) then !d=sign(1.0,d)
	  if(d.gt.0) then
	    d=1.0
	  else 
	    d=-1.
	    endif
	  endif
	real_d=d
	trn(2,1)= acosd(real_d)
900	end 

	subroutine vector1(e,a,f)
chk	==================
c	e -- eigen value
c	a -- matrix
c	f -- eigen vector

	implicit real*8  (a-h)
	real*8 a(3,3),f(3)

	do 10 i=1,3
10	a(i,i)=a(i,i)-e
!	write(*,*) a(1,1),a(1,2), a(1,3)
!	write(*,*) a(2,1),a(2,2), a(2,3)
!	write(*,*) a(3,1),a(3,2), a(3,3)

	do 20 i=1,2
	jb=i+1
	do 20 j=jb,3
	d0=a(i,2)*a(j,3)-a(i,3)*a(j,2)
	d1=a(i,3)*a(j,1)-a(i,1)*a(j,3)
	d2=a(i,1)*a(j,2)-a(i,2)*a(j,1)
	fs=sqrt(d0*d0+d1*d1+d2*d2)
	if(fs.gt.deps()) goto 40
20	continue
	write(6,*) 'rms-W> ABORTED: two structures may be identical.'
	return

40	if     (abs(d0).ge.abs(d1).and.abs(d0).ge.abs(d2)) then
	  if(d0.lt.0.0) fs=-fs
	else if(abs(d1).ge.abs(d2).and.abs(d1).ge.abs(d0)) then
	  if(d1.lt.0.0) fs=-fs
	else if(abs(d2).ge.abs(d0).and.abs(d2).ge.abs(d1)) then
	  if(d2.lt.0.0) fs=-fs
	endif
	
	f(1)=d0/fs
	f(2)=d1/fs
	f(3)=d2/fs
	do 60 i=1,3
60	a(i,i)=a(i,i)+e
	end

	subroutine vector2(c,a,f,e,d)
c	=============================
	implicit real*8  (a-h)
	dimension d(3),e(3),f(3), a(3,3)

	do 30 i=1,3
	do 10 j=1,3
10	f(j)=a(j,i)
	f(i)=f(i)-c
	fr=sqrt(f(1)*f(1)+f(2)*f(2)+f(3)*f(3))
	if(fr.gt.deps())goto 40
30	continue
	write(6,*) '%EdPDB-W- aborted in subroutine vector2'
	return

40	do 45 i=1,3
	f(i)=f(i)/fr
45	e(i)=0.0

	do 50 j=1,3
	if(abs(f(j)).lt.deps())then
	e(j)=1.0
	goto 80
	end if
50 	continue

60	fr=1.0/sqrt(f(1)*f(1)+f(2)*f(2))
	e(1)= f(2)*fr
	e(2)=-f(1)*fr

80	call ixjeqk(f,e,d)
	end

	subroutine root(c1,c2,d0,d1,d2,e)
chk	===============
chk	output: e
	implicit real*8  (a-h)
	include 'edp_dat.inc'
!	use edp_dat

!	implicit real*8  (a-h)

	d=c1
	e=c2
	a=((d+d2)*d+d1)*d
	b=a+d0
	d22=d2*2.0

	do i=1,256
	  if(abs((e-d)/e).gt.deps())then
	    d=d-b/(3.0*d*d+d22*d+d1)
	    a=d**3+d2*d*d+d1*d
	    b=a+d0
	    endif
	  ep=e
	  tmp=(e**3+d2*e*e+d1*e-a)
	  if(abs(tmp).lt.deps()) then
	    return 
	    endif
	  e=d-b*(e-d)/tmp
	  if(abs(ep-e).lt.abs(e)) then
	    if(abs((ep-e)/e).lt.deps()) then
	      return
	      endif
	    endif
	  enddo
	if(verbose.ge.3 ) 
	1 write(6,*) '%EdPDB-W3- uncompleted iteration in the subroutine root'
	end

	subroutine ixjeqk(g1,g2,g3)
chk	==============
	implicit real*8  (a-h)

	dimension g1(3), g2(3), g3(3)
	g3(1)=g1(2)*g2(3)-g1(3)*g2(2)
	g3(2)=g1(3)*g2(1)-g1(1)*g2(3)
	g3(3)=g1(1)*g2(2)-g1(2)*g2(1)
	end

copyright by X. Cai Zhang

c	this subroutine converts eulerian or polar angles 
	subroutine std_euler(*)
chk	====================
	include 'edp_main.inc'
!	use edp_main

	common /cmm_symm/ num_sym, sym0(3,4,max_symm), sym(3,4,max_symm)
	
	parameter (num_options=5)
!      integer, parameter :: num_options=5
	character*10 options(num_options)
	data options/ 
 	1 'to_polar',
	1 'to_euler',
	1 'symmetry',
	1 'asymm_red',
	1 'move_to_o'/
c	1 'pair'/

	real tmp1(3,3), ang(3)
		
	n_of_syn=6					!000515
	syntax(1)='syntax:' 
	syntax(2)='euler to_euler'
	syntax(3)='euler to_polar' 
	syntax(4)='euler symmetry symm_num.i' 
	syntax(5)='euler move_to_o res_id.s atom_name.s' 
	syntax(6)='euler asymm [e1.r e2.r e3.r]'
 
	i=match_l1(num_options,options)
	if( i.le.0 ) return 1

	goto (201, 202, 203, 204, 205) i

chk	to_polar
201	write(6,*) 'euler>       '//
	1 'm''(ph'',om'',kp'')=' //
	1 'm( z, y*, z**)'
	do i=1,n_atom
	  if(lf(i)) then
	    call euler_matrix(x(i), y(i), z(i), tmp1)
	    call get_polar(tmp1,ang)
	    x(i)=ang(1)
	    y(i)=ang(2)
	    z(i)=ang(3)
	    write(text(i)(31:54),1052) ang
	    endif
	  enddo
1052	format(3f8.3)
	return

202	write(6,*) 'eular>       '//
	1 'm''(z'',y*'',z**'')=' //
	1 'm( z, y*, z**)'
	do i=1,n_atom
	  if(lf(i)) then
	    call euler_matrix(x(i), y(i), z(i), tmp1)
	    call get_euler(tmp1,ang)
	    x(i)=ang(1)
	    y(i)=ang(2)
	    z(i)=ang(3)
	    write(text(i)(31:54),1052) ang
	    endif
	  enddo
	return

chk	symmery
203	if(num_sym.le.0) then
	  write(6,*)
	1 '%EdPDB-W- UNDONE: symmetry is undefined'
	  return
	  endif

	call read_ai(1,isym,*900,*900)
	if( isym.le.0 .or. isym.gt.num_sym ) return 1
	write(6,1011) isym

1011	format(' euler>       m''(z'',y*'',z**'')=s(',i2,
	1 ')m( z, y*, z**)')

	do i=1,n_atom
	  if(lf(i)) then
	    call euler_matrix(x(i), y(i), z(i), tmp1)
	    call axbeqc(sym(1,1,isym),tmp1,tmp2)
	    call get_euler(tmp2,ang)
	    x(i)=ang(1)
	    y(i)=ang(2)
	    z(i)=ang(3)
	    write(text(i)(31:54),1052) ang
	    endif
	  enddo
	return

chk	asymmetry
204	call euler_asymm(*900)	  
	return

chk	move_to_o
205	call euler_diff(*900)	  
	return

chk	pair
c206	call euler_pair(*900)	  
c	return

900	return 1
	end

	subroutine std_polar(*)	
chk	====================
	include 'edp_main.inc'
!	use edp_main

	common /cmm_symm/ num_sym, sym0(3,4,max_symm), sym(3,4,max_symm)
	
	parameter (num_options=7)
!	integer, parameter :: num_options=7
        character*10 options(num_options)
        data options/ 
 	1 'to_polar',
	1 'to_euler',
	1 'symmetry',
	1 'asymm_red',
	1 'move_to_o',
	1 'srf_red',
	1 'unique'/

	real tmp1(3,3), tmp2(3,3), ang(3)
		
c	logical  nword
c	character*(max_num_chars) txt
c	common /cmm_txt/n_len,txt,ib,ie

	n_of_syn=8					!000515
	syntax(1)='syntax:' 
	syntax(2)='polar to_polar' 
	syntax(3)='polar to_euler' 
	syntax(4)='polar symmetry symm_num.i '
	syntax(5)='polar move_to_o res_id.s atom_name.s '
	syntax(6)='polar asymm [p1.r p2.r p3.r]' 
	syntax(7)='polar srf_red [p1.r p2.r p3.r]' 
	syntax(8)='polar unique delta_angle.r' 

	i=match_l1(num_options,options)
	if( i.le.0 ) return 1

	goto (201, 202, 203, 204, 205, 206,207) i

chk	to_polar
201	write(6,*) 'polar>       '//
	1 'm''(ph'',om'',kp'')=' //
	1 'm( ph, om, kp)'
	do i=1,n_atom
	  if(lf(i)) then
	    call polar_matrix(x(i), y(i), z(i), tmp1)
	    call get_polar(tmp1,ang)
	    x(i)=ang(1)
	    y(i)=ang(2)
	    z(i)=ang(3)
	    write(text(i)(31:54),1052) ang
	    endif
	  enddo
1052	format(3f8.3)
	return

202	write(6,*) 'polar>       '//
	1 'm''(z'',y*'',z**'')=' //
	1 'm( ph, om, kp)'
	do i=1,n_atom
	  if(lf(i)) then
	    call polar_matrix(x(i), y(i), z(i), tmp1)
	    call get_euler(tmp1,ang)
	    x(i)=ang(1)
	    y(i)=ang(2)
	    z(i)=ang(3)
	    write(text(i)(31:54),1052) ang
	    endif
	  enddo
	return

chk	symmery
203	if(num_sym.le.0) then
	  write(6,*)
	1 '%EdPDB-W- UNDONE: symmetry is undefined'
	  return
	  endif

	call read_ai(1,isym,*900,*900)
	if( isym.le.0 .or. isym.gt.num_sym ) return 1
	write(6,1011) isym
1011	format(' polar>       m''(ph'',om'',kp'')=s(',i2,
	1 ')m( ph, om, kp)')

	do i=1,n_atom
	  if(lf(i)) then
	    call polar_matrix(x(i), y(i), z(i), tmp1)
	    call axbeqc(sym(1,1,isym),tmp1,tmp2)
	    call get_polar(tmp2,ang)
	    x(i)=ang(1)
	    y(i)=ang(2)
	    z(i)=ang(3)
	    write(text(i)(31:54),1052) ang
	    endif
	  enddo
	return

chk	asymmetry
204	call polar_asymm(*900)	  
	return

chk	move_to_o
205	call polar_diff(*900)	  
	return

chk	self-rotation function reduce
206	call polar_srf_red(*900)	  
	return

chk	unique
207	call polar_unique(*900)	  
	return

900	return 1
	end

	subroutine polar_asymm(*)
chk	========================
chk	reduce self rotation function solutions.

	include 'edp_main.inc'
!	use edp_main

	common /cmm_symm/ num_sym, sym0(3,4,max_symm), sym(3,4,max_symm)
	
	real tmp(3,3), tmp1(3,3), ang(3), ps(3),pt(3)
		
	if(num_sym.le.0) then
	  write(6,*)
	1 '%EdPDB-W- UNDONE: symmetry is undefined'
	  return
	  endif

	ps(1)=0.0
	ps(2)=0.0
	ps(3)=0.0
	call read_ar(3, ps, *900,*50)

50	do i=1,n_atom
	  if(lf(i)) then
	    call polar_matrix(x(i), y(i), z(i), tmp1)

	    s_min=1.e10
	    do isym=1, num_sym
	      call axbeqc(sym(1,1,isym),tmp1,tmp)
	      call get_polar(tmp,ang)
	      s_curr= polar_diff_func(ps,ang)
	      if(s_curr.lt.s_min) then
	        s_min=s_curr
	        pt(1)=ang(1)
	        pt(2)=ang(2)
	        pt(3)=ang(3)
	        endif
	      enddo

	    x(i)=pt(1)
	    y(i)=pt(2)
	    z(i)=pt(3)
c	    b(i)=s_min
	    write(text(i)(31:54),1052) pt
1052	    format(3f8.3)
	    endif
	  enddo
	return

900	return 1
	end

	subroutine polar_unique(*)
chk	========================
	include 'edp_main.inc'
!	use edp_main

	real tmp(3,3), tmp0(3,3), tmp1(3,3)
		
	call read_ar(1, da, *900,*900)
	da=abs(da)

	do i=1,n_atom
	  if(lf(i)) then
	    call polar_matrix(x(i), y(i),-z(i), tmp0)
	    do j=i+1, n_atom
	      if(lf(j)) then
	        call polar_matrix(x(j), y(j), z(j), tmp1)
	        call axbeqc(tmp1,tmp0,tmp)

	        t=0.5*(tmp(1,1)+tmp(2,2)+tmp(3,3)-1.0)
	        if(t.ge.1.0) then
	          t=0.0
	        else if(t.le.-1.0) then
	          t=180.0
	        else
	          t= abs(acosd(t))
	          endif

	        if(t.le.da) lf(j)=.false.
	        endif
	      enddo
	    endif
	  enddo
	return

900	return 1
	end

	subroutine polar_srf_red(*)
chk	========================
	include 'edp_main.inc'
!	use edp_main

	common /cmm_symm/ num_sym, sym0(3,4,max_symm), sym(3,4,max_symm)
	
	real tmp(3,3), tmp0(3,3), tmp1(3,3), ang(3), ps(3),pt(3)
		
	if(num_sym.le.0) then
	  write(6,*)
	1 '%EdPDB-W- UNDONE: symmetry is undefined'
	  return
	  endif

	ps(1)=0.0
	ps(2)=0.0
	ps(3)=0.0
	call read_ar(3, ps, *900,*50)

50	do i=1,n_atom
	  if(lf(i)) then
	    call polar_matrix(x(i), y(i), z(i), tmp0)

	    s_min=1.e10
	    do jsym=1, num_sym
	      call axbeqc(tmp0,sym(1,1,jsym),tmp1)

	      do isym=1, num_sym
	        call axbeqc(sym(1,1,isym),tmp1,tmp)
	        call get_polar(tmp,ang)
	        s_curr= polar_diff_func(ps,ang)
	        if(s_curr.lt.s_min) then
	          s_min=s_curr
	          pt(1)=ang(1)
	          pt(2)=ang(2)
	          pt(3)=ang(3)
	          endif
	        enddo
	      enddo

	    x(i)=pt(1)
	    y(i)=pt(2)
	    z(i)=pt(3)
c	    b(i)=s_min
	    write(text(i)(31:54),1052) pt
1052	    format(3f8.3)
	    endif
	  enddo
	return

900	return 1
	end

	subroutine euler_asymm(*)
chk	======================
	include 'edp_main.inc'
!	use edp_main

	common /cmm_symm/ num_sym, sym0(3,4,max_symm), sym(3,4,max_symm)
	
	real tmp(3,3), tmp1(3,3), ang(3), ps(3),pt(3)
		
	if(num_sym.le.0) then
	  write(6,*)
	1 '%EdPDB-W- UNDONE: symmetry is undefined'
	  return
	  endif

	ps(1)=0.0
	ps(2)=0.0
	ps(3)=0.0
	call read_ar(3, ps, *900,*50)

50	do i=1,n_atom
	  if(lf(i)) then
	    call euler_matrix(x(i), y(i), z(i), tmp1)

	    s_min=1.e10
	    do isym=1, num_sym
	      call axbeqc(sym(1,1,isym),tmp1,tmp)
	      call get_euler(tmp,ang)
	      s_curr= euler_diff_func(ps,ang)
	      if(s_curr.lt.s_min) then
	        s_min=s_curr
	        pt(1)=ang(1)
	        pt(2)=ang(2)
	        pt(3)=ang(3)
	        endif
	      enddo

	    x(i)=pt(1)
	    y(i)=pt(2)
	    z(i)=pt(3)
c	    b(i)=s_min
	    write(text(i)(31:54),1052) pt
1052	    format(3f8.3)
	    endif
	  enddo
	return

900	return 1
	end

	subroutine polar_matrix(p1,p2,p3, m)
chk	=======================
c       this subroutine calculates the matrix
c       from the polar angles
	real m(3,3)

	al=sind(p2)*cosd(p1)
	am=sind(p2)*sind(p1)
	an=cosd(p2)
	cos_z3=cosd(p3)
	sin_z3=sind(p3)
	m(1,1)=al*al+(am*am+an*an)*cos_z3
	m(1,2)=al*am*(1.-cos_z3)-an*sin_z3
	m(1,3)=al*an*(1.-cos_z3)+am*sin_z3
	m(2,1)=al*am*(1.-cos_z3)+an*sin_z3
	m(2,2)=am*am+(al*al+an*an)*cos_z3
	m(2,3)=am*an*(1.-cos_z3)-al*sin_z3
	m(3,1)=al*an*(1.-cos_z3)-am*sin_z3
	m(3,2)=am*an*(1.-cos_z3)+al*sin_z3
	m(3,3)=an*an+(al*al+am*am)*cos_z3
	end

	subroutine euler_matrix(e1, e2, e3, m)
chk     =======================
	real m(3,3), tmp1(3,3), tmp2(3,3)
	call zrot(e1, m)
	call yrot(e2, tmp1)
	call axbeqc(m, tmp1, tmp2)
	call zrot(e3, tmp1)
	call axbeqc(tmp2, tmp1, m)
	end
 
	real function euler_diff_func(a,b)
chk     =============================
	real a(3), b(3), tmp1(3,3), tmp2(3,3), tmp3(3,3)

	call euler_matrix(a(1), a(2), a(3), tmp1)
	call euler_matrix(-b(3), -b(2),-b(1), tmp2)
	call axbeqc(tmp1,tmp2,tmp3)
	call get_polar(tmp3,c)
	t=0.5*(tmp3(1,1)+tmp3(2,2)+tmp3(3,3)-1.0)
	if(t.ge.1.0) then
	  euler_diff_func=0.0
	else if(t.le.-1.0) then
	  euler_diff_func=180.0
	else
	  euler_diff_func= abs(acosd(t))
	  endif
	end

	real function polar_diff_func(a,b)
chk	=============================
	real a(3), b(3), tmp1(3,3), tmp2(3,3), tmp3(3,3)

	call polar_matrix(a(1), a(2), a(3), tmp1)
	call polar_matrix(b(1), b(2),-b(3), tmp2)
	call axbeqc(tmp1,tmp2,tmp3)
	t=0.5*(tmp3(1,1)+tmp3(2,2)+tmp3(3,3)-1.0)
	if(t.ge.1.0) then
	  polar_diff_func=0.0
	else if(t.le.-1.0) then
	  polar_diff_func=180.0
	else
	  polar_diff_func= abs(acosd(t))
	  endif
	end	

	subroutine polar_diff(*)
chk	====================
	include 'edp_main.inc'
!	use edp_main

	real tmp1(3,3), tmp2(3,3), tmp3(3,3), pt(3)

	call get_atoms(1, i, ierr) 
	if(ierr.ne.0) return 1

	call polar_matrix(x(i), y(i), -z(i), tmp1)

50	do i=1,n_atom
	  if(lf(i)) then
	    call polar_matrix(x(i), y(i), z(i), tmp2)
	    call axbeqc(tmp2, tmp1, tmp3) 
	    call get_polar(tmp3,pt)
	    x(i)=pt(1)
	    y(i)=pt(2)
	    z(i)=pt(3)
	    write(text(i)(31:54),1052) pt
1052	    format(3f8.3)
	    endif
	  enddo
	end

	subroutine euler_diff(*)
chk	====================
	include 'edp_main.inc'
!	use edp_main

	real tmp1(3,3), tmp2(3,3), tmp3(3,3), pt(3)

	call get_atoms(1, i, ierr) 
	if(ierr.ne.0) return 1
	call euler_matrix(-z(i), -y(i), -x(i), tmp1)

50	do i=1,n_atom
	  if(lf(i)) then
	    call euler_matrix(x(i), y(i), z(i), tmp2)
	    call axbeqc(tmp2, tmp1, tmp3) 
	    call get_euler(tmp3,pt)
	    x(i)=pt(1)
	    y(i)=pt(2)
	    z(i)=pt(3)
	    write(text(i)(31:54),1052) pt
1052	    format(3f8.3)
	    endif
	  enddo
	end

c***	end of euler.for

copyright by X. Cai Zhang
