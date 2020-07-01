
	subroutine av_cos(*)
chk	=================
	return 1
	end

	subroutine eigen(a,e,v,need_vect, *)
chk	===============
	implicit real*8 (a-h)

	dimension g1(3),g2(3),g3(3),a(3,3), e(3), v(3,3)
	logical need_vect

	d0=a(1,1)*(a(2,3)**2-a(2,2)*a(3,3))
	1 +a(2,2)* a(3,1)**2-a(1,2)*a(2,3)*a(3,1)*2.
	2 +a(3,3)* a(1,2)**2
	d1=a(1,1)*a(2,2)+a(2,2)*a(3,3)+a(3,3)*a(1,1)
	1 -a(1,2)**2    -a(2,3)**2    -a(3,1)**2
	d2=-a(1,1)-a(2,2)-a(3,3)

	if(abs(d0).lt.deps()) then
	  e3=0.0

	  if(abs(d1).gt.deps()) then       ! 2 order equation
	    d0=sqrt(d2*d2-4.0*d1)
	    e1=.50*(-d2+d0)
	    e2=.50*(-d2-d0)
            goto 170
            endif
	  e2=0.0

	  if(abs(d2).lt.deps()) then  ! 0 order equation
	    write(6,*) 'eigen> det(m)=',d2,'; eps=', deps()
	    return 1
	    endif

	  e1=sqrt(-d2)                ! 1 order equation
	  if(need_vect) call vector1(e1,a,g1)
	  return
	  endif

c	3 order equation
	dc=sqrt(d2*d2-3.*d1)
	if(dc.lt.deps())then     ! trible root
	  e1=(-d2/3.0)		     !
	  e2=e1                 !
	  e3=e1                 !
	  do i=1,3
	    g1(i)=0.
	    g2(i)=0.
	    g3(i)=0.
	    enddo
	  g1(1)=1.
	  g2(2)=1.
	  g3(3)=1.
	  goto 180 		!
	  end if                ! isotropic distribution

	c1=(-d2-dc)/3.
	c2=(-d2+dc)/3.
	fc1=c1**3+d2*c1**2+d1*c1+d0
	fc2=c2**3+d2*c2**2+d1*c2+d0

	if(fc1.lt.deps())then       ! double root
	  e3=c1
	  e2=e3
	  call vector2(c1,a,g1,g2,g3)
	  e1=((a(1,1)*g1(1)+a(1,2)*g1(2)+a(1,3)*g1(3))/g1(1))
	  goto 180
	else if(abs(fc2).lt.deps())then  ! double root
	  e1=c2
	  e2=e1
	  call vector2(c2,a,g3,g1,g2)
	  e3=((a(1,1)*g3(1)+a(1,2)*g3(2)+a(1,3)*g3(3))/g3(1))
	  goto 180
	  end if

	a0=0.
	call root(a0,c1,d0,d1,d2,e3)
	c3=-d2/3.
	fc3=c3**3+d2*c3*c3+d1*c3+d0
	if(abs(fc3).lt.deps())then
	  e2=c3
	else if(fc3.gt.0.0)then
	  call root(c3,c2,d0,d1,d2,e2)
	else
	  call root(c3,c1,d0,d1,d2,e2)
	  end if

	do i=1,256
	  c3=c3+c2*i
	  fc3=c3**3+d2*c3*c3+d1*c3+d0
	  if(fc2.lt.0.0.and.fc3.gt.0.0 .or.
	1   fc2.gt.0.0.and.fc3.lt.0.0 ) goto 160
	  enddo
	write(6,*) 'eigen-W> ABORTED.'
	return 1

160	call root(c3,c2,d0,d1,d2,e1)
170	if(need_vect) then 
	  call vector1(e1,a,g1)
	  call vector1(e2,a,g2)
	  call ixjeqk(g1,g2,g3)
	  endif

180	e(1)=e1
	e(2)=e2
	e(3)=e3

	if(need_vect)  then
	  do i=1,3
	    v(1,i)=g1(i)
	    v(2,i)=g2(i)
	    v(3,i)=g3(i)
	    enddo
	  endif
	end

	subroutine planar(*)
chk	=================
	include 'edp_main.inc'
!	use edp_main

	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie
	dimension ax(3), xx(3,3), f(3), a(3,3), xi(3), af(3,3)
	logical  nword, fix

	n_of_syn=2					!000515
	syntax(1)='syntax:'
	syntax(2)='planar [vector_id.s [FIX]]'
	
	fix=.false.
	call get_vector_id(iv, *100)
	if(.not.nword(n_len,txt,ib,ie)) then
	  if(txt(ib:ie).ne. 'fix') goto 900
	  ax(1)=vectors(1,iv)
	  ax(2)=vectors(2,iv) 
	  ax(3)=vectors(3,iv) 
	  f(1) =vectors(4,iv) 
	  f(2) =vectors(5,iv) 
	  f(3) =vectors(6,iv) 
	  goto 200
	endif

100	do j=1,3
	  ax(j)=0.
	  do k=1,3
	    xx(k,j)=0.
	  enddo
	enddo

	ii= 0
	do  i= 1, n_atom
	 if( lf(i)) then
	   ii= ii+ 1
	   xi(1)= x(i)
	   xi(2)= y(i)
	   xi(3)= z(i)
	   do j=1,3
	     ax(j)= ax(j)+ xi(j)
	     do k=1,3
	       xx(k,j)= xx(k,j)+ xi(k)*xi(j)
	     enddo
	   enddo
	 endif
	enddo

	if( ii .le. 2) then
	  write(6,*) 'planar-W> UNDONE: too few on atoms.'
	  return
	endif
	do j=1,3
	  ax(j)= ax(j)/ii
	  do k=1,3
	    xx(k,j)= xx(k,j)/ii
	  enddo
	enddo

	do j=1,3
	  do k=1,3
	    a(k,j)= xx(k,j)- ax(k)*ax(j)
	  enddo
	enddo

	call uxveqw(a(1,1),a(1,2),af(1,3))
	call uxveqw(a(1,2),a(1,3),af(1,1))
	call uxveqw(a(1,3),a(1,1),af(1,2))

	a1=af(1,1)*af(1,1)+af(2,1)*af(2,1)+af(3,1)*af(3,1)
	a2=af(1,2)*af(1,2)+af(2,2)*af(2,2)+af(3,2)*af(3,2)
	a3=af(1,3)*af(1,3)+af(2,3)*af(2,3)+af(3,3)*af(3,3)
	if(a1.ge. a2 .and. a1.ge. a3) then
	    aa=1./sqrt(a1)
	    f(1)=af(1,1)*aa
	    f(2)=af(2,1)*aa
	    f(3)=af(3,1)*aa
	else if(a2.ge.a3 .and. a2.ge.a1) then
	    aa=1./sqrt(a2)
	    f(1)=af(1,2)*aa
	    f(2)=af(2,2)*aa
	    f(3)=af(3,2)*aa
	else if(a3.ge.a1 .and. a3.ge.a2) then
	    aa=1./sqrt(a3)
	    f(1)=af(1,3)*aa
	    f(2)=af(2,3)*aa
	    f(3)=af(3,3)*aa
	endif

200	aa= 0.
	do i= 1, n_atom
	  if( lf(i)) then
	    w(i)= f(1)*(x(i)-ax(1))+ 
	1         f(2)*(y(i)-ax(2))+ 
	1         f(3)*(z(i)-ax(3))
	    aa= aa+ w(i)*w(i)
	  endif
	enddo

	if( iv .gt. 0 .and. .not.fix ) then
	  vectors(1,iv)= ax(1)
	  vectors(2,iv)= ax(2)
	  vectors(3,iv)= ax(3)
	  vectors(4,iv)= f(1)
	  vectors(5,iv)= f(2)
	  vectors(6,iv)= f(3)
	  vectors(7,iv)= 1.0
	endif
	
	write(6,1001) 
c	1 ( vectors(i,iv),i=1,6), 
	1 sqrt(aa/ii)
1001	format(' planar> the occs of on atoms have been changed to the'
	1     /'         projected distance. rms =',f8.3)
c	1     /'         the plane passes through point (',3f8.3,')' 
c	1     /'                   with a normal vector (',3f8.3,')' ) 
	return
900	return 1
	end	

	subroutine vector_vector(*)
chk	========================
	include 'edp_main.inc'
!	use edp_main

c	character*2 cvector(10)
c	data cvector/'v0','v1','v2','v3','v4','v5','v6','v7','v8','v9'/

	real f(3), g(3), p(3), q(3)

	call find_vector_id(iv1,*900,*900)
	call find_vector_id(iv2,*900,*900)
	call  get_vector_id(iv3,*50)

50	theta= acosd(vectors(4,iv1)*vectors(4,iv2)+
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
	theta1= acosd(vectors(4,iv1)*f(1)+
	1             vectors(5,iv1)*f(2)+
	1             vectors(6,iv1)*f(3))
	theta2= acosd(-vectors(4,iv2)*f(1)
	1             -vectors(5,iv2)*f(2)
	1             -vectors(6,iv2)*f(3))

	call uxveqw(vectors(4,iv1),vectors(4,iv2),f)
	aa= sqrt( f(1)*f(1)+f(2)*f(2)+f(3)*f(3))
	if( abs( aa) .le. eps() ) then
            write(6,*) 'vv-W> UNDONE: '//cvector(iv1)//cvector(iv2)//
	1          ' is ambiguous.'
	    return
	else
	    f(1)=f(1)/aa
	    f(2)=f(2)/aa
	    f(3)=f(3)/aa
	endif

	dist= f(1)*(vectors(1,iv2)-vectors(1,iv1))+
	1     f(2)*(vectors(2,iv2)-vectors(2,iv1))+
	1     f(3)*(vectors(3,iv2)-vectors(3,iv1))

	dist1= (vectors(1,iv2)-vectors(1,iv1))**2
	1     +(vectors(2,iv2)-vectors(2,iv1))**2
	1     +(vectors(3,iv2)-vectors(3,iv1))**2

	if( iv3 .gt. 0) then
c
c	p(t)=p1+v1*t
c	t*(v1xv2)=fx((p1-p2)xf)
c
	  call uxveqw( vectors(4,iv1), vectors(4,iv2), g)
	  do i= 1, 3
	    p(i)= vectors(i,iv2)- vectors(i,iv1)- dist*f(i)
	  enddo
	  call uxveqw( p, vectors(4,iv2), q)
	  do i=1, 3
	    if( abs(g(i)) .gt. eps() ) then
	      aa=q(i)/g(i)
	      goto 100
	    endif
	  enddo
100	  vectors(1,iv3)= vectors(1,iv1)+ aa*vectors(4,iv1)
	  vectors(2,iv3)= vectors(2,iv1)+ aa*vectors(5,iv1)
	  vectors(3,iv3)= vectors(3,iv1)+ aa*vectors(6,iv1)
	  vectors(4,iv3)= f(1)
	  vectors(5,iv3)= f(2)
	  vectors(6,iv3)= f(3)
	  vectors(7,iv3)= 1.0
	endif

	tmp=theta
	if(tmp.gt.90.0)  tmp=180.0-tmp
	tmp1=theta1
	if(tmp1.gt.90.0) tmp1=180.0-tmp1
	tmp2=theta2
	if(tmp2.gt.90.0) tmp2=180.0-tmp2
	write(6,1001) theta, tmp,
	1             theta1,tmp1,
	1             theta2,tmp2,
	1             abs( dist), sqrt(dist1)
1001	format(
	1 ' vv> angle(v^v)      =',f8.2,',     acute_angle=',f8.2
	1/'     angle(v1^p1->p2)=',f8.2,',     acute_angle=',f8.2
	1/'     angle(v2^p2->p1)=',f8.2,',     acute_angle=',f8.2
	1/'     the shortest dist. between the two vectors=',f8.2
	1/'     the dist.  between the two starting points=',f8.2)
	return
900	return 1
	end

	subroutine point_vector(*)
chk	=======================
	include 'edp_main.inc'
!	use edp_main
	dimension ax(3)

	call find_vector_id(iv,*900,*900)

	call get_atoms(1, i, ierr) 
	if(ierr.ne.0) then
	  if(ierr.eq.2) then
	    do i=1,n_atom
	      if(lf(i)) goto 100
	      enddo
	    write(6,*) 
	1 'pv-W> one on atom is needed here to make the connection.'
	    endif
	  return 1
	  endif

c	all the vectors should be normalized.
100	ax(1)= x(i)- vectors(1,iv)
	ax(2)= y(i)- vectors(2,iv)
	ax(3)= z(i)- vectors(3,iv)
	dist= sqrt( ax(1)*ax(1)+ax(2)*ax(2)+ax(3)*ax(3))
	if(dist.le.1.e-12) then
	  write(6,*) 'pv-W> the vector and the point are colinear'
	  return 1
	  endif
	
	theta= (vectors(4,iv)*ax(1)+
	1       vectors(5,iv)*ax(2)+
	1       vectors(6,iv)*ax(3))/dist
	if( theta .gt. 1.) then
	  theta= 0.
	else if( theta .lt. -1.) then
	  theta= 180.
	else
	  theta= acosd( theta)
	endif
	write(6,1001) theta,dist, sind(theta)*dist
1001	format(
	1 ' pv> angle((p-v)^v)=',f8.2
	1,' pv>  distance(p-v)=',f8.2
	1,' pv>      proj_dist=',f8.2)
	return
900	return 1
	end

	subroutine vector_point(*)
chk	=======================
	include 'edp_main.inc'
!	use edp_main
	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie

	call find_vector_id(iv, *900,*900)

	call get_atoms(1, i, ierr) 

	r_length=vectors(7,iv)
	call read_ar(1, r_length, *900,*100)
100	continue

	if(ierr.ne.0) then
	  if(ierr.eq.2) then
	    do i=1,n_atom
	      if(lf(i)) goto 110
	      enddo
	    write(6,*) 
	1'vp-W> one ON atom is needed here to store the result.'
	    endif
	  return 1
	  endif

c	all the vectors should be normalized.
110	x(i)= vectors(1,iv)+vectors(4,iv)*r_length
	y(i)= vectors(2,iv)+vectors(5,iv)*r_length
	z(i)= vectors(3,iv)+vectors(6,iv)*r_length
	write( text(i)(31:54), 1061, err=901) x(i),y(i),z(i)
1061	format(3f8.3)
	return
900	return 1
901	write(6,*) 'vp-W> ERROR in output conversation.'
	end

c****	end of planar.for

copyright by X. Cai Zhang

	subroutine polar(trn,ierr)
chk	==========================
c	this subroutine calculates the polar 
c	(ie. alf,beta and gama/ phi, omega and kappa)	
c	angles from a rotation matrix.

	parameter (eps=1.e-4)
!	real, parameter :: eps=1.e-4
	real trn(3,3), a_polar(3)
	real kappa
	logical local
	save	!030422

	local=.false.
	ierr=1
	d= trn(1,1)*( trn(2,2)*trn(3,3) - trn(3,2)*trn(2,3))
	1 +trn(1,2)*( trn(2,3)*trn(3,1) - trn(3,3)*trn(2,1))
	1 +trn(1,3)*( trn(2,1)*trn(3,2) - trn(3,1)*trn(2,2))

	if( abs(1.-d) .gt. eps) then 
	 write(6,*) 'axis-W> it is not a real rotation. det=', d
	 return
	endif
	ierr=0
	goto 10	

	entry get_polar(trn,a_polar)
c	===============
	local=.true.

c	l=sin(omega)*cos(phi)
c	m=sin(omege)*sin(phi)
c	n=cos(omega)

c	0<=phi<=pi
c	0<=omega<=pi
c	-pi<kappa<=pi

c	ll+(mm+nn)cos_k, 	lm(1-cos_k)-n.sin_k, 	ln(1-cos_k)+m.sin_k
c	lm(1-cos_k)+n.sin_k,	mm+(ll+nn)cos_k,	mn(1-cos_k)-l.sin_k
c	ln(1-cos_k)-m.sin_k,	mn(1-cos_k)+l.sin_k,	nn+(ll+mm)cos_k

10	cos_k= max(-1.,(trn(1,1)+trn(2,2)+trn(3,3)-1.)*0.5)

	if(1.-cos_k.lt.eps) then
	  phi=0.
	  omega=0.
	  kappa=0.

	else if( abs(1.-trn(3,3)) .lt. eps) then
	  phi=0.
	  omega=0.
	  kappa= acosd( cos_k)
	  if( trn(2,1)-trn(1,2) .lt. 0.) kappa= -kappa

	else if( abs(1.-trn(2,2)) .lt. eps) then
	  phi=90.
	  omega=90.
	  kappa= acosd( cos_k)
	  if( trn(1,3)-trn(3,1) .lt. 0.) kappa= -kappa

	else if( abs(1.-trn(1,1)) .lt. eps) then
	  phi=0.
	  omega=90.
	  kappa= acosd( cos_k)
	  if( trn(3,2)-trn(2,3) .lt. 0.) kappa= -kappa

	else if( abs(trn(1,3)+trn(3,1)) .ge. eps) then 
	  phi= atand( (trn(3,2)+trn(2,3))/(trn(1,3)+trn(3,1))) 
	  if( phi .lt. 0.) phi= phi + 180.

	  kappa= acosd( cos_k)

	  if(trn(3,3)- cos_k  .le. eps) then
	    omega=90.0
chk	    sin_k.2m = (T13-T31), where m = sin(omega).sin(phi)
	    if( (trn(1,3)-trn(3,1))*sind(phi) .lt. 0.0) kappa=-kappa
	  else
	    omega= acosd( sqrt( (trn(3,3)-cos_k ) / (1.- cos_k) ))
	    if( (trn(3,1)+trn(1,3))/cosd(phi).lt.0.) omega=180.-omega
	    if( ( trn(2,1)-trn(1,2)) / cosd( omega) .lt. 0.) kappa= -kappa
	    endif
	
	else if( abs(trn(2,3)+trn(3,2)) .ge. eps) then 
	  phi= 90.

	  omega= acosd( sqrt( max(0.,(trn(3,3)-cos_k)) / (1.- cos_k) ))
	  if(trn(2,3).lt.0.) omega=180.- omega

	  kappa= acosd( cos_k)
	  if( abs(cosd(omega)) .ge. eps ) then 
	    if( ( trn(2,1)-trn(1,2)) / cosd( omega) .lt. 0.) kappa= -kappa
	    endif
	else if( abs(trn(1,2)+trn(2,1)) .ge. eps) then 
	  omega=90.

	  phi= acosd( sqrt( (trn(1,1)-cos_k) / (1.- cos_k) ))
	  if( trn(1,2).lt. 0.) phi=180.- phi

	  kappa= acosd( cos_k)
	  if(0. .lt. phi .and. phi .lt. 180.) then
	    if ( trn(1,3)-trn(3,1) .lt. 0.) kappa= -kappa
	  else
	    if( abs(cosd(phi)) .ge. eps) then 
	      if ( (trn(3,2)-trn(2,3)) / cosd(phi) .lt. 0.) kappa= -kappa
	      endif
	    endif

	else
	  write(6,1010) trn
1010	  format(' subroutine polar can not handle the matrix,'/(3e12.4))
	  return
	  endif    

	if(abs(abs(phi)-180.0).lt.0.02) then
	  phi=0.
	  omega=180.-omega
	  kappa=-kappa
	  endif

	if(local) then
	  a_polar(1)=phi
	  a_polar(2)=omega
	  a_polar(3)=kappa
	  return
	  endif
	
	write(6,1106)   phi, omega, kappa
1106	format(  ' the polar angles              are ',3f10.3)
	end
chk***	end of polar.for

copyright by X. Cai zhang

