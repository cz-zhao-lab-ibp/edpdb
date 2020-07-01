	subroutine update(*)
chk	=================
	include 'edp_main.inc'
!	use edp_main
	include 'edp_dat.inc'
!	use edp_dat

        common /main_para/ tolower
	1 ,jnk(m_main_para-1)
        logical tolower
  
	character*6 options(4)
	data options/'xyz','b','weight','text'/
	character*(108) file_name
	character*(32) form		!,skip
	character*(80) tmp

        character*(max_num_chars)  txt
        common /cmm_txt/n_len,txt,ib,ie
        logical nword

	logical  open_file1
	external open_file1

	integer jt(2)
	equivalence (jt1,jt(1)),(jt2,jt(2))

	n_of_syn=5					!000515
	syntax(1)='syntax:' 
	syntax(2)='update ' 
	syntax(3)='update (xyz, w, b) input_filename.s fortran_format.s' 
	syntax(4)=
	1'update t column_1.i column_2.i input_filename.s fortran_format.s' 
	syntax(5)='  [jump_after_string.s]'

	call find1(4,options,i_option)
	if(i_option.le.0) then

	j=0
	do i=1,n_atom
	  if(lf(i)) then
	    read(text(i)(31:66),1009,err=900) x(i), y(i), z(i), w(i), b(i)
1009	format(3f8.3,2f6.2)
	    j=j+1
	    endif
	  enddo
	  goto 200	
	else if(i_option.gt.4) then
	  return 1
	endif
	
	if(i_option.eq.4) then
	  call read_ai(2,jt,*900,*900)
	  if(jt1.le.0 .or. jt2.gt.72 .or. jt1.gt.jt2) return 1
	  endif

	if(.not. open_file1(29,file_name,'old','.txt')) return 1
	
	if(nword(n_len,txt,ib,ie)) then 
	  form='*'
	  if(verbose.ge.2 ) write(6,1005)
1005	  format(' update-I2> free format will be used.')
	else
	  ie=n_len
	  ic=ichar(delimiter)
	  if( ic .eq. ichar(txt(ib:ib))) then !use '(...)' to input the user prefored output-format
	    ib= ib+1
	    do j=ib, n_len
	      if( ichar(txt(j:j)).eq.ic) then
	        ie=j-1
	        txt(j:j)=' '  !000303, it may screw up the history file.
	        goto 80
	        endif
	      enddo
	    endif
80	  form= txt(ib:ie)
	  if(verbose.ge.2 ) write(6,1004) form
1004	  format(' update-I2> the input format will be: ',a)
	  endif

	if(.not.nword(n_len,txt,ib,ie)) then 
	  ie=n_len
	  ic=ichar(delimiter)
	  if( ic .eq. ichar(txt(ib:ib)) ) then !use '(...)' to input the skip string
	    ib= ib+1
	    do j=ib, n_len
	      if( ichar(txt(j:j)).eq.ic) then
	        ie=j-1
	        txt(j:j)=' '  !000303, it may screw up the history file.
	        goto 84
	        endif
	      enddo
	    endif
84	 if(ie.ge.ib) then
	  num_skip_lines=0
85	  read(29,1002,end=920) tmp
	  num_skip_lines=num_skip_lines+1
	  if(verbose.ge.4 ) write(6,*) tmp
	  if(index(tmp,txt(ib:ie)).le.0) goto 85
	  if(verbose.ge.2 ) write(6,1006) num_skip_lines, txt(ib:ie)
1002	  format(a)
1006	  format(' update-I2>',i4,
	1' lines are skipped until the first [',a,'] is met.')
	  endif
!	  write(*,*) 'ok1'
	else
!	  write(*,*) 'ok2'
	 endif

	j=0
	do i=1,n_atom
	  if(lf(i)) then
	    if(i_option.eq.1) 	then		! update xyz
	      if(form.ne.'*') then
	         read(29,form,err=900,end=910) x(i), y(i), z(i)
chk			^ 
chk		the format must have correct syntax, if on an alpha machine.
	      else
	         read(29,*,err=900,end=910) x(i), y(i), z(i)
	         endif
	      write(text(i)(31:54),1001,err=900) x(i),y(i),z(i)
1001	      format(3f8.3)
	    else if(i_option.eq.2) 	then	! update b-factor
	      if(form.ne.'*') then
	        read(29,form,err=900,end=910) b(i)
	      else
	        read(29,*,err=900,end=910) b(i)
	        endif
	    else if(i_option.eq.3) 	then	! update occ
	      if(form.ne.'*') then
	        read(29,form,err=900,end=910) w(i)
	      else
	        read(29,*,err=900,end=910) w(i)
	        endif
	    else if(i_option.eq.4) 	then	! update text
	      if(form.ne.'*') then
	        read(29,form,err=900,end=910)  text(i)(jt1:jt2)
	      else
	        read(29,'(a)',err=900,end=910) text(i)(jt1:jt2) 
	        endif
	      if(tolower) then
	        do jc=jt1,jt2
	          ic=ichar(text(i)(jc:jc))
	          if(65.le.ic.and.ic.le.90) text(i)(jc:jc)=char(ic+32)
	          enddo
	        endif
	      endif
	    j=j+1
	    endif
	  enddo
200	if(verbose.ge.1 ) write(6,1007) j
1007	format(' update>     ',i6,' records have been updated.')
	return
900	errmsg=' errmsg: format error'
	return 1
920	errmsg=' errmsg: the starting text string not found'
	return 1
910	write(6,1003) j
1003	format(' update-W> only',i6,' records have been updated.',
	1 ' end of file detected.')
	end

	subroutine extract(*)
c	=====================
	include 'edp_main.inc'
! 	use edp_main

	character*(max_num_chars) txt
	common /cmm_txt/ n_len,txt,ib,ie
	logical nword
	logical lf1(max_atom)

	n_of_syn=2					!000515
	syntax(1)='syntax:' 
	syntax(2)='load group2.s | extract group1 sub_group1 '

	call dfgroup(igr,*901)

	ie=0
	if( nword(n_len,txt,ib,ie))  goto 901	! get ie back

	igr1= match_l( max_gr, cgroup)
	if( igr1 .le. 0) then
	  errmsg=' errmsg: a group name is needed'
	  return 1
	end if
	n_group1= n_groupa(igr1)
	if( n_group1.le.0) then
	  write(6,*) 
	1 'extract-W> UNDONE: define group ['//cgroup(igr1)//'] first.'
	  return
	else if(n_group1.ne.n_groupa(igr)) then
	  write(6,*) 
	1 'extract-W> UNDONE: num(',cgroup(igr1),')=/=num(',cgroup(igr),
	1        ').'
	  return
	  end if

	igr2= match_l( max_gr, cgroup)

	if( igr2 .le. 0 .or. igr2 .eq. igr1) then
	  errmsg=' errmsg: another group name needed'
	  return 1
	else if( n_groupa(igr2) .le.0) then
	  write(6,*) 
	1'extract-W> UNDONE: define group ['//cgroup(igr2)//'] first.'
	  return
	  end if
	n_group2= n_groupa(igr2)

	do i=1,n_atom
	  lf1(i)=.false.
	  enddo
	j1=n_groupa(igr2)
	do j=1,j1
	  i=igroupa( j,igr2)
	  lf1(i)=.true.
	  enddo
	  
	n=0
	j1=n_groupa(igr1)
	do j=1,j1
	  i=igroupa( j,igr1)
	  if(lf1(i)) then
	    i=igroupa( j,igr)
	    lf(i)=incl
	    n=n+1
	    endif
	  enddo
	if(n.le.0) write(6,*) 
	1 'extract-W> UNDONE: the overlap between ',
	1 cgroup(igr1),' & ',cgroup(igr2),
	1 ' is zero.'
	return

901	errmsg= ' errmsg: wrong group/zone information'
	return 1
	end

	subroutine dock(*)
chk	===============
 	include 'edp_main.inc'
 	include 'edp_dat.inc'
!	use edp_main

	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie
	logical nword 

	dimension rad(max_atom), isort(max_atom)

	call read_ar(4,rad, *901, *901)
	if(rad(1).le.rad(2)) then 
	  zc=rad(1)
	  ze=rad(2)
	else 
	  zc=rad(2)
	  ze=rad(1)
	endif
	
	if(rad(3).gt.0.0) then 
	  dz=rad(3)
	else 
	  goto 901
	endif
	
	if(rad(4).gt. 0.0) then
	  shell=rad(4)
	else
	  goto 901
	endif
		
200	if(zc.le.ze) then
	 j=0
	 do i=1,n_atom
	  if(lf(i)) then
	    j=j+1
	    rad(j)=sqrt(x(i)**2+y(i)**2+(z(i)-zc)**2)
	    if(verbose .ge. 6) w(i)= rad(j)
	  endif
	 enddo
	 if(j.le.2) goto 901
	 call sort_rl(j,rad,isort)
	 
	 r_min=rad(isort(1))
	 n=0
	 do i=1,j
	  if((rad(isort(i))-r_min).le.shell) then
	    n=n+1
	  else 
	    write(6,1002) r_min, n
	    zc=zc+dz
	    goto 200
	  endif
	 enddo
	endif
	return
	
1002	format(' dock> r=',f8.1,', # of atoms in touch=',i6)

901	return 1
	end

c***	end of pair.for

copyright by X. Cai Zhang
