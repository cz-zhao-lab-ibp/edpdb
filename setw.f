	subroutine doit(*)
chk	===============
	include 'edp_main.inc'
! 	use edp_main
	include 'edp_dat.inc'
!	use edp_dat
	include 'edp_file.inc'
!	use edp_file

	parameter 	(io=48)
!	integer, parameter :: io=48
	character*(max_num_chars) txt
	common /cmm_txt/ n_len,txt,ib,ie
!	logical nword, nword0

	character*6 cnum
	logical L_r2, L_r3, L_r4, L_r5, L_a3, L_a4, L_a5, L_t4, L_t5
	real d12(5), d23(5), d34(5), d45(5)
	1, d13(5), d24(5), d35(5), a13(5), a24(5), a35(5)
	1, t14(5), t25(5)
	character*1 d2c, d3c, d4c, d5c, a3c, a4c, a5c, t4c, t5c

c	n_of_syn=6		
c	syntax(1)='syntax: in one-line' 
c	syntax(2)='doit [gr2.s S/R s2.i d2min.r d2max.r '
c	syntax(3)=' [gr3.s S/R s3.i d3min.r d3max.r A/D a3min.r a3max.r'
c	syntax(4)=' [gr4.s S/R s4.i d4min.r d4max.r A/D a4min.r a4max.r T/H t4min.r t4max.r'
c	syntax(5)=' [gr5.s S/R s5.i d5min.r d5max.r A/D a5min.r a5max.r T/H t5min.r t5max.r'
c	syntax(6)=' ]]]]'
	
	status=-2
	ksum=0
	imax=max_atom
	iskip2=2	!as switches in doit_inp1 
	iskip3=3
	iskip4=4
	iskip5=5
	
	igr2=0
	igr3=0
	igr4=0
	igr5=0
	
	n_group1=0
	do i= 1, n_atom
	  if( lf(i)) n_group1=n_group1+1
	  enddo
	if( n_group1 .le. 0) then
	  write(6,*) 'doit-W>  UNDONE: no on atom exist.'
	  return
	  end if

chk	input group2
	call doit_inp1(igr2,n_group2,L_r2,iskip2,d2c,d2_0,d2_1,*990,*666)
	
chk	input group3
	call doit_inp1(igr3,n_group3,L_r3,iskip3,d3c,d3_0,d3_1,*990,*800)
chk	'a' for angle or 'd' for distance
	call doit_inp2(L_a3,a3c,'a','d',a3_0,a3_1,0.,180.,0.,4.,*990)
	if(L_a3.and.a3_1-a3_0.gt. 360.0) return 1
	
chk	input group4
	call doit_inp1(igr4,n_group4,L_r4,iskip4,d4c,d4_0,d4_1,*990, *800)
	call doit_inp2(L_a4,a4c,'a','d',a4_0,a4_1,0.,180.,0.,4.,*990)
	if(L_a4.and.a4_1-a4_0.gt. 360.0) return 1
chk	'h' for theta or 't' for 'torsion'
	call doit_inp2(L_t4,t4c,'t','h',t4_0,t4_1,-180.,180.,-90.,90.,*990)
	if(t4_1-t4_0.gt. 360.0) return 1
			
chk	input group5
	call doit_inp1(igr5,n_group5,L_r5,iskip5,d5c,d5_0,d5_1,*990, *800)
	call doit_inp2(L_a5,a5c,'a','d',a5_0,a5_1,0.,180.,0.,4.,*990)
	if(L_a5.and.a5_1-a5_0.gt. 360.0) return 1
	call doit_inp2(L_t5,t5c,'t','h',t5_0,t5_1,-180.,180.,-90.,90.,*990)
	if(t5_1-t5_0.gt. 360.0) return 1

800	continue
	if(verbose .ge. 4 ) then
	  write(6,1201) '['//
	1  cgroup(igr2)(:ltrim(cgroup(igr2)))//']'
	1, d2c,iskip2,d2_0,d2_1
	  if(igr3.gt.0) write(6,1201) '['//
	1  cgroup(igr3)(:ltrim(cgroup(igr3)))//']'
	1, d3c,iskip3,d3_0,d3_1, a3c, a3_0,a3_1 
	  if(igr4.gt.0) write(6,1201) '['//
	1  cgroup(igr4)(:ltrim(cgroup(igr4)))//']'
	1, d4c,iskip4,d4_0,d4_1,a4c, a4_0,a4_1,t4c, t4_0,t4_1
	  if(igr5.gt.0) write(6,1201) '['//
	1  cgroup(igr5)(:ltrim(cgroup(igr5)))//']'
	1, d5c,iskip5,d5_0,d5_1,a5c, a5_0,a5_1,t5c, t5_0,t5_1
	  endif

	r2_min=d2_0*d2_0
	r2_max=d2_1*d2_1
	r3_min=d3_0*d3_0
	r3_max=d3_1*d3_1
	r4_min=d4_0*d4_0
	r4_max=d4_1*d4_1
	r5_min=d5_0*d5_0
	r5_max=d5_1*d5_1

	n_scr=0
	do i=2,3
	  d12(i)=0.0
	  d23(i)=0.0
	  d34(i)=0.0
	  d45(i)=0.0
	  d13(i)=0.0
	  d24(i)=0.0
	  d35(i)=0.0
	  a13(i)=0.0
	  a24(i)=0.0
	  a35(i)=0.0
	  t14(i)=0.0
	  t25(i)=0.0
	end do
	  d12(4)=1.e6
	  d23(4)=1.e6
	  d34(4)=1.e6
	  d45(4)=1.e6
	  d13(4)=1.e6
	  d24(4)=1.e6
	  d35(4)=1.e6
	  a13(4)=1.e6
	  a24(4)=1.e6
	  a35(4)=1.e6
	  t14(4)=1.e6
	  t25(4)=1.e6
	  
	  d12(5)=-1.e6
	  d23(5)=-1.e6
	  d34(5)=-1.e6
	  d45(5)=-1.e6
	  d13(5)=-1.e6
	  d24(5)=-1.e6
	  d35(5)=-1.e6
	  a13(5)=-1.e6
	  a24(5)=-1.e6
	  a35(5)=-1.e6
	  t14(5)=-1.e6
	  t25(5)=-1.e6
	  
	do ir1= 1, n_res
	  do 100 i1= ijk(ir1), ijk(ir1+1)-1
	    if( .not.lf(i1)) goto 100
	    do 200 j2=1,n_group2
	      i2= igroupa(j2,igr2)
	      ir2=aa_seq(i2)
	      if( L_r2) then
		if(ir2-ir1 .ne. iskip2) goto 200
	      else ! skip iskip2 residue(s)
	        if( abs(ir2-ir1) .lt. iskip2) goto 200
	      endif

	      call chk_distance(d12(1),i1,i2,d2_1, r2_min, r2_max,*200)
	      
	if(igr3.gt.0) then 
	  do 300 j3=1,n_group3
	      i3= igroupa(j3,igr3)
	      ir3=aa_seq(i3)
	      if( L_r3) then
		if(ir3-ir2 .ne. iskip3) goto 300
	      else ! skip iskip3 residue(s)
	        if( abs(ir3-ir2) .lt. iskip3) goto 300
	      endif

	      call chk_distance(d23(1) ,i2,i3,d3_1, r3_min, r3_max,*300)
	      call chk_angle(a13(1),d13(1),i1,i2,i3,a3_0, a3_1,L_a3,*300)
	      
	if(igr4.gt.0) then 
	  do 400 j4=1,n_group4
	      i4= igroupa(j4,igr4)
	      ir4=aa_seq(i4)
	      if( L_r4) then
		if(ir4-ir3 .ne. iskip4) goto 400
	      else ! skip iskip4 residue(s)
	        if( abs(ir4-ir3) .lt. iskip4) goto 400
	      endif

	      call chk_distance(d34(1),i3,i4,d4_1, r4_min, r4_max,*400)
	      call chk_angle(a24(1),d24(1),i2,i3,i4,a4_0, a4_1,L_a4, *400)
	      call chk_torsion(t14(1),i1,i2,i3,i4,t4_0, t4_1,L_t4,a13(1),*400)

	if(igr5.gt.0) then 
	  do 500 j5=1,n_group5
	      i5= igroupa(j5,igr5)
	      ir5=aa_seq(i5)
	      if( L_r5) then
		if(ir5-ir4 .ne. iskip5) goto 500
	      else ! skip iskip5 residue(s)
	        if( abs(ir5-ir4) .lt. iskip5) goto 500
	      endif

	      call chk_distance(d45(1) ,i4,i5,d5_1, r5_min, r5_max,*500)
	      call chk_angle(a35(1),d35(1),i3,i4,i5,a5_0, a5_1,L_a5,*500)
	      call chk_torsion(t25(1),i2,i3,i4,i5,t5_0,t5_1,L_t5,a24(1),*500)
	      	      
	      call stat1(d12)
	      call stat1(d23)
	      call stat1(d34)
	      call stat1(d45)
	      call stat1(d13)
	      call stat1(d24)
	      call stat1(d35)
	      call stat1(a13)
	      call stat1(a24)
	      call stat1(a35)
	      call stat1(t14)
	      call stat1(t25)
				
	      write(io,1005) 
	1 text(i1)(13:26),text(i2)(13:26),text(i3)(13:26),text(i4)(13:26)
	1,text(i5)(13:26)
	1,d2c//d3c//a3c//d4c//a4c//t4c//d5c//a5c//t5c//'='
	1,d12(1)
	1,d23(1),two21(L_a3,a13(1),d13(1))
	1,d34(1),two21(L_a4,a24(1),d24(1)),t14(1)
	1,d45(1),two21(L_a5,a35(1),d35(1)),t25(1)
	
	      ksum=ksum+1
c	      if( ksum .ge. imax) goto 900
	      if( n_scr+5 .gt. imax) goto 900
	      igroupa(n_scr+1,1)=i1      
	      igroupa(n_scr+2,1)=i2      
	      igroupa(n_scr+3,1)=i3      
	      igroupa(n_scr+4,1)=i4      
	      igroupa(n_scr+5,1)=i5      
	      n_scr=n_scr+5  
500	  end do
	else	      
	      call stat1(d12)
	      call stat1(d23)
	      call stat1(d34)
	      call stat1(d13)
	      call stat1(d24)
	      call stat1(a13)
	      call stat1(a24)
	      call stat1(t14)
	     			      
	      write(io,1004) 
	1 text(i1)(13:26),text(i2)(13:26),text(i3)(13:26),text(i4)(13:26)
	1,d2c//d3c//a3c//d4c//a4c//t4c//'='
	1,d12(1)
	1,d23(1),two21(L_a3,a13(1),d13(1))
	1,d34(1),       two21(L_a4,a24(1),d24(1)),t14(1)

	      ksum=ksum+1
c	      if( ksum .ge. imax) goto 900
	      if( n_scr+4 .gt. imax) goto 900
	      igroupa(n_scr+1,1)=i1      
	      igroupa(n_scr+2,1)=i2      
	      igroupa(n_scr+3,1)=i3      
	      igroupa(n_scr+4,1)=i4      
	      n_scr=n_scr+4   
	endif
400	  end do
	else
	      call stat1(d12)
	      call stat1(d23)
	      call stat1(d13)
	      call stat1(a13)
	     		      
	      write(io,1003) 
	1 text(i1)(13:26),text(i2)(13:26),text(i3)(13:26)
	1,d2c//d3c//a3c//'='
	1,d12(1)
	1,d23(1),two21(L_a3,a13(1),d13(1))

	      ksum=ksum+1
c	      if( ksum .ge. imax) goto 900
	      if( n_scr+3 .gt. imax) goto 900
	      igroupa(n_scr+1,1)=i1      
	      igroupa(n_scr+2,1)=i2      
	      igroupa(n_scr+3,1)=i3      
	      n_scr=n_scr+3   
	endif	      	      
300	  end do
	else		
	      call stat1(d12)
	    	    
	      write(io,1002) text(i1)(13:26),text(i2)(13:26)
	1,d2c//'='
	1,d12(1)
	      ksum=ksum+1
c	      if( ksum .ge. imax) goto 900
	      if( n_scr+2 .gt. imax) goto 900
	      igroupa(n_scr+1,1)=i1      
	      igroupa(n_scr+2,1)=i2 
	      n_scr=n_scr+2   
	endif
200	    end do 
100	  end do
	end do

	n_groupa(1)=n_scr
	if(igr5.gt.0) then
	write( 6,1008)'list', ksum
	1             ,'on', n_group1,
	1      cgroup(igr2), n_group2,
	1      cgroup(igr3), n_group3,
	1      cgroup(igr4), n_group4,
	1      cgroup(igr5), n_group5
	write(io,1008)'list', ksum
	1             ,'on', n_group1,
	1      cgroup(igr2), n_group2,
	1      cgroup(igr3), n_group3,
	1      cgroup(igr4), n_group4,
	1      cgroup(igr5), n_group5
	else if(igr4.gt.0) then
	write( 6,1008)'list', ksum
	1             ,'on', n_group1,
	1      cgroup(igr2), n_group2,
	1      cgroup(igr3), n_group3,
	1      cgroup(igr4), n_group4
	write(io,1008)'list', ksum
	1             ,'on', n_group1,
	1      cgroup(igr2), n_group2,
	1      cgroup(igr3), n_group3,
	1      cgroup(igr4), n_group4
	else if(igr3.gt.0) then
	write( 6,1008)'list', ksum
	1             ,'on', n_group1,
	1      cgroup(igr2), n_group2,
	1      cgroup(igr3), n_group3
	write(io,1008)'list', ksum
	1             ,'on', n_group1,
	1      cgroup(igr2), n_group2,
	1      cgroup(igr3), n_group3
	else 
	write( 6,1008)'list', ksum
	1             ,'on', n_group1,
	1      cgroup(igr2), n_group2
	write(io,1008)'list', ksum
	1             ,'on', n_group1,
	1      cgroup(igr2), n_group2
	endif
	
	if( ksum .gt. 0)then
	    call stat2(d12,ksum,io,d2c,'12')
	  
	  if( igr3 .gt. 0) then
	    call stat2(d23,ksum,io,d3c,'23')
	    call stat2(a13,ksum,io,' ','a123')
	    call stat2(d13,ksum,io,' ','d13')
	  
	  if( igr4 .gt. 0) then
	    call stat2(d34,ksum,io,d4c,'34')
	    call stat2(a24,ksum,io,' ','a234')
	    call stat2(d24,ksum,io,' ','d24')
	    call stat2(t14,ksum,io,t4c,'1234')

	  if( igr5 .gt. 0) then
	    call stat2(d45,ksum,io,d5c,'45')
	    call stat2(a35,ksum,io,' ','a345')
	    call stat2(d35,ksum,io,' ','d35')
	    call stat2(t25,ksum,io,t5c,'2345')

	  endif ! igr5
	  endif ! igr4
	  endif ! igr3
	end if
 
	call typelist(io)
	status=ksum
666	return
1002	format(1x,2(a14,','),a, f7.2)
1003	format(1x,3(a14,','),a,2f7.2,f7.1)
1004	format(1x,4(a14,',')/1x,a,2f7.2,f7.1,f7.2,2f7.1)
1005	format(1x,5(a14,',')/1x,a,2f7.2,f7.1,f7.2,2f7.1,f7.2,2f7.1)
1008	format('!doit> ', 6(a4,':',i5,', '))
1201	format(' doit-I4>',a,t15,a1,':',i3,2f6.1,2(', ',a1,':',2f6.1)) 
1062	format(i6)

900	write(cnum,1062) imax
	errmsg=
	1' errmsg: the number of selection is larger than '//cnum
	status=-3
990	return 1
	end

	real function two21(L,a,b)
!	===================	
	logical L
	if(L) then
	  two21=a
	else
	  two21=b
	endif
	end

	subroutine doit_inp2(L_th,thc, tc, hc, th_0, th_1,t0,t1, h0,h1,*)
chk	=====================
	include 'edp_main.inc'
! 	use edp_main
	
	character*(max_num_chars) txt
	common /cmm_txt/ n_len,txt,ib,ie
	logical nword0
!	logical nword, nword0

	logical L_th
	character*1 thc, tc, hc
	
	if(nword0(n_len,txt,ib,ie) ) then
	  L_th=.true.
	else if(txt(ib:ie).eq.hc) then
	    L_th=.false.
	else if(txt(ib:ie).eq.tc) then
	    L_th=.true.
	else
	    return 1
	endif
		  
	if(L_th) then
	  th_0=t0
	  th_1=t1
	  thc=tc
	else
	  th_0=h0
	  th_1=h1
	  thc=hc
	endif
	call read_ar(1, th_0, *990, *430)
430	call read_ar(1, th_1, *990, *425)
425	if(th_1 .lt. th_0) then
	  tmp= th_1
	  th_1= th_0
	  th_0= tmp
	  endif
	return
990	return 1
	end
	
	subroutine doit_inp1(igr2,n_group2,L_r2,iskip2,d2c,d2_0,d2_1,*,*)
chk	=====================
	include 'edp_main.inc'
! 	use edp_main
	include 'edp_dat.inc'
!	use edp_dat
	include 'edp_file.inc'
!	use edp_file

	character*(max_num_chars) txt
	common /cmm_txt/ n_len,txt,ib,ie
	logical nword0
!	logical nword, nword0

	logical L_r2 

	character*1 d2c 
	
	igr2= match_l( max_gr, cgroup)
	if( igr2 .eq. 0) then
	  if(iskip2.eq.2) call analy(1)
	  return 2
	else if(igr2 .lt. 0) then
	  errmsg=' errmsg: a group name needed'
	  return 1
	  end if
	n_group2= n_groupa(igr2)
	if( n_group2.le.0) then
	  write(6,*) 
	1 'doit-W> UNDONE: define group ['//cgroup(igr2)//'] first.'
	  return
	  end if
	  
chk	's' for skip or 'r' for 'regest'
	if(nword0(n_len,txt,ib,ie) ) then
	  L_r2=.true.
	else if(txt(ib:ie).eq.'s') then
	    L_r2=.false.
	else if(txt(ib:ie).eq.'r') then
	    L_r2=.true.
	else 
	    return 1
	  endif
	  
chk	input iskip
	if(L_r2) then
	  iskip2=0
	  d2c='r'
	  d2_0=2.0  
	else
	  iskip2=1
	  d2c='s'
	  d2_0=3.0
	endif
	call read_ai(1, iskip2, *990, *250)
	if(.not.L_r2) iskip2=abs(iskip2)

chk	input d2_0, d2_1
250 	d2_1=0.0
	call read_ar(1, d2_0, *990, *246)
246	call read_ar(1, d2_1, *990, *245)
245	if(d2_1.lt.d2_0) then
	  r2_max=d2_1
	  d2_1=d2_0
	  d2_0=r2_max
	  endif
	if(d2_0.lt.0.) return 1
	return
990	return 1
	end
	
	subroutine stat1(a)
chk	================
	real a(5)
	a(2)=a(2)+a(1)
	a(3)=a(3)+a(1)*a(1)
	a(4)=min(a(4),a(1))
	a(5)=max(a(5),a(1))
	end
	
	subroutine stat2(a,k,io,c1,c2)
	character*(*) c1,c2
	character*(5) c12
	real a(5)
c	if(k.gt.0) then
	  a(2)=a(2)/k
	  a(3)=sqrt(max(0., a(3)/k-a(2)*a(2)))
c	endif

	i=len(c1)
	j=len(c2)
	c12(:i)=c1
	c12(i+1:)=c2
	write( 6,1018) c12,(a(i),i=2,5)
	write(io,1018) c12,(a(i),i=2,5)
1018	format('!doit> av+sgm+min+max(',a5,')',4f9.3) 
	end
	
	subroutine chk_torsion(t,i1,i2,i3,i4, t_min, t_max,L_t4, a123,*)
 	include 'edp_main.inc'
!	use edp_main
	logical L_t4
	
	a1=x(i2)-x(i1)
	b1=x(i3)-x(i2)
	c1=x(i4)-x(i3)
	a2=y(i2)-y(i1)
	b2=y(i3)-y(i2)
	c2=y(i4)-y(i3)
	a3=z(i2)-z(i1)
	b3=z(i3)-z(i2)
	c3=z(i4)-z(i3)
	axb1=a2*b3-a3*b2
	axb2=a3*b1-a1*b3
	axb3=a1*b2-a2*b1
	bxc1=b2*c3-b3*c2
10	bxc2=b3*c1-b1*c3
	bxc3=b1*c2-b2*c1
	ab=axb1*axb1+axb2*axb2+axb3*axb3
	if(ab.lt.1.e-6) return 1
	bc=bxc1*bxc1+bxc2*bxc2+bxc3*bxc3
	if(bc.lt.1.e-6) return 1
	vol=b1*(axb2*bxc3-axb3*bxc2)+	
	1   b2*(axb3*bxc1-axb1*bxc3)+	
	1   b3*(axb1*bxc2-axb2*bxc1)
	f_tor= (axb1*bxc1+axb2*bxc2+axb3*bxc3)/sqrt(ab*bc)
	if(abs(f_tor).ge.1.) then
	  c1=c1+0.001			! force the result away from single point
	  goto 10
	end if
	
	t= sign(acosd(f_tor),vol)
	if(.not.L_t4) then		! L_t4: true-torsion, false-theta
	 t=asind(sind(t)*sind(a123))
	 endif

	if(t .lt. t_min) then
	  t= t+ 360.0
	else if(t .gt. t_max) then
	  t= t- 360.0
	endif
	if(t .lt. t_min .or. t .gt. t_max) return 1
	end
 	
	subroutine chk_angle(a,d,i1,i2,i3,a_min, a_max,L,*)
 	include 'edp_main.inc'
!	use edp_main
	logical L
	
	if(L) then ! check angle
	  x1=x(i1)-x(i2)
	  y1=y(i1)-y(i2)
	  z1=z(i1)-z(i2)
	  x2=x(i3)-x(i2)
	  y2=y(i3)-y(i2)
	  z2=z(i3)-z(i2)
	  r=sqrt((x1*x1+y1*y1+z1*z1)*(x2*x2+y2*y2+z2*z2))
	  if(r.lt. eps()*10.0) return 1
	  a= acosd((x1*x2+y1*y2+z1*z2)/r)
	  if(a .lt. a_min) then
	    a= a+ 360.0
	  else if(a .gt. a_max) then
	    a= a- 360.0
	  endif
	  if(a .lt. a_min .or. a .gt. a_max) return 1
	  x2=x(i3)-x(i1)
	  y2=y(i3)-y(i1)
	  z2=z(i3)-z(i1)
	  d=sqrt(x2*x2+y2*y2+z2*z2)
	else ! check distance
	  x2=x(i3)-x(i1)
	  y2=y(i3)-y(i1)
	  z2=z(i3)-z(i1)
	  d=sqrt(x2*x2+y2*y2+z2*z2)
	  if(d .lt. a_min .or. d .gt. a_max) return 1
	  x1=x(i1)-x(i2)
	  y1=y(i1)-y(i2)
	  z1=z(i1)-z(i2)
	  x2=x(i3)-x(i2)
	  y2=y(i3)-y(i2)
	  z2=z(i3)-z(i2)
	  r=sqrt((x1*x1+y1*y1+z1*z1)*(x2*x2+y2*y2+z2*z2))
	  if(r.lt. eps()*10.0) return 1
	  a= acosd((x1*x2+y1*y2+z1*z2)/r)	
	endif	
	end
	
	subroutine chk_distance(dr,i,j,dmax, rsmin, rsmax, *)
chk	=======================
 	include 'edp_main.inc'
! 	use edp_main

	dx=x(i)-x(j)
	if( abs(dx) .gt. dmax) return 1
	dy=y(i)-y(j)
	if( abs(dy) .gt. dmax) return 1
	dz=z(i)-z(j)
	if( abs(dz) .gt. dmax) return 1
	dr=dx*dx+dy*dy+dz*dz
	if(dr .lt. rsmin .or. rsmax .lt. dr) return 1
	dr=sqrt(dr)
	end
	
	subroutine setw(*)
chk	===============
 	include 'edp_main.inc'
!	use edp_main

	character*(108) file_name

	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie
	logical nword 

	dimension rad(max_atom)
	character*2 ss1,ss2,ss3,ss4
	character*4 symb_ca
	data ss1,ss2,ss3/'+w','-w','*w'/
	ss4='/w'
	file_name='*'

	if(nword(n_len,txt,ib,ie)) then
	  write(6,*) 'setw> w=j_res*0.1'
	  do j=1,n_res
	    do i=ijk(j),ijk(j+1)-1
	      if(lf(i)) w(i)=j*0.1
	      end do
	    end do
	  return
	else if(txt(ib:ie).eq.'x') then
	 do i=1,n_atom
	  if(lf(i)) 
	1   read( text(i), 1001, err=900) rad(i)
	 enddo
	else if(txt(ib:ie).eq.'y') then
	 do i=1,n_atom
	  if(lf(i)) 
	1  read( text(i), 1002, err=900) rad(i)
	 enddo
	else if(txt(ib:ie).eq.'z') then
	 do i=1,n_atom
	  if(lf(i)) 
	1  read( text(i), 1003, err=900) rad(i)
	 enddo
	else if(txt(ib:ie).eq.'b') then
	 do i=1,n_atom
	  if(lf(i)) rad(i)=b(i)
	 enddo
	else if(txt(ib:ie).eq.'avw') then	!010330
	 do i=1,n_res
	  kk=0
	  av=0.
	  do j=ijk(i), ijk(i+1)-1
	   if(lf(j)) then
	    av = av+w(j)
	    kk=kk+1
	   end if
	  end do
	  if(kk.gt.0) then
	   av=av/kk
	   do j=ijk(i), ijk(i+1)-1
	    if(lf(j)) rad(j)=av
	   end do
	  end if
	 end do
	else if(txt(ib:ie).eq.'ca') then	!010330
	 call get_ca(symb_ca)
	 do i=1,n_res
	  do j0=ijk(i), ijk(i+1)-1
	   if(atom(j0).eq.symb_ca) goto 151
	  end do
	  j0=ijk(i)
151	  continue
	  do j=ijk(i), ijk(i+1)-1
	    if(lf(j)) rad(j)=w(j0)
	  end do
	 end do
	else
	 read( txt(ib:ie), *, err=800) wv
	 do i=1,n_atom
	  if(lf(i)) rad(i)=wv
	 enddo
	endif

100	if(nword(n_len,txt,ib,ie)) then
	  if(file_name.ne.'*') call read_vdw( file_name,rad,*901)	 
	  do i=1,n_atom
	   if(lf(i)) w(i)=rad(i)
	  enddo
	else if(txt(ib:ie).eq.ss1(1:min(len(ss1),ie-ib+1))) then
	  if(file_name.ne.'*') call read_vdw( file_name,rad,*901)	 
	  do i=1,n_atom
	   if(lf(i)) w(i)=rad(i)+w(i)
	  enddo
	else if(txt(ib:ie).eq.ss2(1:min(len(ss2),ie-ib+1))) then
	  if(file_name.ne.'*') call read_vdw( file_name,rad,*901)	 
	  do i=1,n_atom
	   if(lf(i)) w(i)=rad(i)-w(i)
	  enddo
	else if(txt(ib:ie).eq.ss3(1:min(len(ss3),ie-ib+1))) then
	  if(file_name.ne.'*') call read_vdw( file_name,rad,*901)	 
	  do i=1,n_atom
	   if(lf(i)) w(i)=rad(i)*w(i)
	  enddo
	else if(txt(ib:ie).eq.ss4(1:min(len(ss4),ie-ib+1))) then
	  if(file_name.ne.'*') call read_vdw( file_name,rad,*901)	 
	  do i=1,n_atom
	   if(lf(i)) then 
	    if(w(i).eq.0.) goto 900
	    w(i)=rad(i)/w(i)
	   endif
	  enddo
	else 
	  goto 901
	endif
	return
1001	format(30x,f8.3)
1002	format(38x,f8.3)
1003	format(46x,f8.3)

800	file_name=txt(ib:ie)
	goto 100

900	errmsg=' errmgs: error during converting '
901	return 1
	end
