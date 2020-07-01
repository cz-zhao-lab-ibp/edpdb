c	this subroutine calculates the harker peaks
c	  of a given fractional coordinates (the first (two) on atom(s)).
c
	subroutine aniso(*)
	include 'edp_main.inc'
	real 	a_array(6), v(3,3)
	real*8 	a(3,3), e(3)
	character*(6) my_vec(3)
	
	status=-2
	call read_ar(6,a_array, *900,*900)
	if(a_array(1).le.0.0.or.a_array(2).le.0.0.or.a_array(3).le.0.0) 
	1 goto 900
	scale=1.0
	call read_ar(1,scale, *900,*40)
40	if(scale.le.0.0) goto 900
	call get_vector_id(iv1 ,*10)
10	call get_vector_id(iv2 ,*20)
20	call get_vector_id(iv3 ,*30)	
30	call get_atoms(1,id_atom,ierr)
	if(ierr.eq.1) then
	  goto 900
	else if(ierr.eq. 2) then
	  id_atom=0
	endif	
	status=-3
	a(1,1)=a_array(1)*scale
	a(1,2)=a_array(6)*scale
	a(1,3)=a_array(5)*scale
	a(2,1)=a_array(6)*scale
	a(2,2)=a_array(2)*scale
	a(2,3)=a_array(4)*scale
	a(3,1)=a_array(5)*scale
	a(3,2)=a_array(4)*scale
	a(3,3)=a_array(3)*scale
	call eigen(a,e,v,.true.,*900)
	if( iv1.gt.0) then
	  my_vec(1)='['//cvector(iv1)//']'
	  if(id_atom.gt.0) then
	    vectors(1,iv1)= x(id_atom)
	    vectors(2,iv1)= y(id_atom)
	    vectors(3,iv1)= z(id_atom)
	  else
	    vectors(1,iv1)= 0.0
	    vectors(2,iv1)= 0.0
	    vectors(3,iv1)= 0.0
	    endif
	  vectors(4,iv1)= v(1,1)
	  vectors(5,iv1)= v(1,2)
	  vectors(6,iv1)= v(1,3)
	  vectors(7,iv1)= e(1)
	else
	  my_vec(1)='    v:'
	  endif
	if( iv2.gt.0) then
	  my_vec(2)='['//cvector(iv2)//']'
	  if(id_atom.gt.0) then
	    vectors(1,iv2)= x(id_atom)
	    vectors(2,iv2)= y(id_atom)
	    vectors(3,iv2)= z(id_atom)
	  else
	    vectors(1,iv2)= 0.0
	    vectors(2,iv2)= 0.0
	    vectors(3,iv2)= 0.0
	    endif
	  vectors(4,iv2)= v(2,1)
	  vectors(5,iv2)= v(2,2)
	  vectors(6,iv2)= v(2,3)
	  vectors(7,iv2)= e(2)
	else
	  my_vec(2)='    v:'
	  endif
	if( iv3.gt.0) then
	  my_vec(3)='['//cvector(iv3)//']'
	  if(id_atom.gt.0) then
	    vectors(1,iv3)= x(id_atom)
	    vectors(2,iv3)= y(id_atom)
	    vectors(3,iv3)= z(id_atom)
	  else
	    vectors(1,iv3)= 0.0
	    vectors(2,iv3)= 0.0
	    vectors(3,iv3)= 0.0
	    endif
	  vectors(4,iv3)= v(3,1)
	  vectors(5,iv3)= v(3,2)
	  vectors(6,iv3)= v(3,3)
	  vectors(7,iv3)= e(3)
	else
	  my_vec(3)='    v:'
	  endif
	aniso_factor=min(e(1),e(2),e(3))/max(e(1),e(2),e(3))
	write(6,1002) 
	1 my_vec(1), v(1,1),v(1,2),v(1,3), e(1),
	1 my_vec(2), v(2,1),v(2,2),v(2,3), e(2),(e(1)+e(2)+e(3))/3.0,
	1 my_vec(3), v(3,1),v(3,2),v(3,3), e(3),aniso_factor
1002	format(' aniso>',
	1(t9,a6,3f8.3,', e1: ',e10.3/),
	1(t9,a6,3f8.3,', e2: ',e10.3,',   <e>: ',f8.3/),
	1(t9,a6,3f8.3,', e3: ',e10.3,', aniso: ',f8.3))
	if(id_atom.gt.0) w(id_atom)=aniso_factor
	status=0
	return
900	return 1
	end

	subroutine harker(*)
chk	=================
	include 'edp_main.inc'
!	use edp_main

	common /cmm_symm/ num_sym, sym(3,4,max_symm)
	1 ,jnk1(3,4,max_symm)
	character*5 scross
	data  scross/'cross'/

	character*50 out
	character*4 s(0:24)
	data s/'    ','***1','***2','***3','+1/6','***5','+1/4',
	1      '***7','+1/3','***9','**10','**11','+1/2','**13',
	1      '**14','**15','+2/3','**17','+3/4','**19','+5/6',
	1      '**21','**22','**23','**24'/
	character*3 xyzl(-3:3,3)
	data xyzl/'-3x','-2x',' -x','   ','  x',' 2x',' 3x',
	1         '-3y','-2y',' -y','   ',' +y','+2y','+3y',
	1         '-3z','-2z',' -z','   ',' +z','+2z','+3z'/
	real d(3,4), p(4,2), q(4), sc(3)
		
	logical  nword
	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie

	n_of_syn=2					!000515
	syntax(1)='syntax:' 
	syntax(2)=
	1'harker [grid_a.i grid_b.i grid_c.i] [symm_1.i [symm_2.i]] [cross]' 
c---
	isym0=1
	jsym0=1
	isym1=num_sym
	jsym1=num_sym
	sc(1)=1.
	sc(2)=1.
	sc(3)=1.
	lcross=0

	call read_ar(3,sc, *900,*40)
40	if(sc(1).le.0.0 .or.sc(2).le.0.0 .or.sc(3).le.0.0) goto 900
	
	call read_ai(1, isym1, *900,*55)
	if(isym1.le.0.or.isym1.gt.num_sym) goto 900
	isym0=isym1
55	call read_ai(1, jsym1, *900,*60)
	if(jsym0.le.0.or.jsym1.gt.num_sym) goto 900
	jsym0=jsym1

60	if(.not.nword(n_len,txt,ib,ie)) then
	  if(txt(ib:ie).ne.scross(1:min(len(scross),ie-ib+1))) goto 900
	  lcross=1
	  endif

100	j=0
	do i=1,n_atom
	  if(lf(i)) then
	    j=j+1
	    p(1,j)=x(i)/sc(1)
	    p(2,j)=y(i)/sc(2)
	    p(3,j)=z(i)/sc(3)
	    p(4,j)=w(i)
	    if(j.gt.lcross) goto 200
	    endif
	  enddo
	  write(6,*) 'harker-W> UNDONE: select one or two records first.'
	  return
900	  return 1

200	do isym= isym0, isym1
	do jsym= jsym0, jsym1
	  if(isym.eq.jsym.and.lcross.eq.0) goto 500
	  if(lcross.eq.0) then
	    do i=1,3
	     do j=1,4
	      d(i,j)=sym(i,j,jsym)-sym(i,j,isym)
	      enddo
	     d(i,4)=mod(d(i,4),1.)
	     if(d(i,4).lt.0.) d(i,4)=1.+ d(i,4)
	     enddo
	   write(out,1001) 
	1 ((xyzl(nint(d(i,j)),j),j=1,3),s(nint(d(i,4)*24.)),i=1,3)
1001	format(3(3a3,a4,', '))
	   j=1
	   do i=1,50
	     if(out(i:i).eq.'+'.and.(out(j:j).eq.','.or.j.eq.1)) out(i:i)=' ' 
	     if(out(i:i).eq.','.and.out(j:j).eq.',') then
	       out(i-1:i-1)='0'
	       j=i-1
	       endif
	     if(out(i:i).eq.',') then
	       j=j+1
	       out(i:i)=' '
	       out(j:j)=','
	       endif	      
	     if(out(i:i).ne.' ') j=i
	     enddo

	   do i=1,3
	     q(i)= mod(d(i,1)*p(1,1)+d(i,2)*p(2,1)+d(i,3)*p(3,1)+d(i,4),1.)
	     if(q(i).lt.0.) q(i)=q(i)+1.
	     q(i)=q(i)*sc(i)
	     enddo
	   if(q(1).eq.0..and.q(2).eq.0..and.q(3).eq.0.) then
	     write(48,1051) out
1051	     format(2x,a50)
	   else
	     write(48,1052) out,q
1052	     format(1x,a50,4f7.2)
	     endif
	else
	  do i=1,3
	    q(i)= mod(sym(i,1,jsym)*p(1,2) +sym(i,2,jsym)*p(2,2)
	1            +sym(i,3,jsym)*p(3,2) +sym(i,4,jsym)
	1           -(sym(i,1,isym)*p(1,1) +sym(i,2,isym)*p(2,1)
	1            +sym(i,3,isym)*p(3,1) +sym(i,4,isym)), 1.)
	    if(q(i).lt.0.) q(i)=q(i)+1.
	    q(i)=q(i)*sc(i)
	    enddo
	  write(48,1002) jsym,isym,q
1002	  format(' symm(',i2,')*x2-symm(',i2,')x1 = ',t52,4f7.2)
	  endif 
500	 enddo
	 enddo
	call typelist(48)
	end
c----	end of harker.for

copyright by X. Cai Zhang
