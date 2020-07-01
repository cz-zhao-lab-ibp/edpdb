	subroutine find1(n_lab,label,i)
chk	================
	include 'edp_dim.inc'
!	use edp_dim

	character*(*) label(n_lab)
	character*16  tmp

	character*(max_num_chars)  txt
	common /cmm_txt/n_len,txt,ib,ie
	logical nword0
	parameter (max_string_length=16)
!	integer, parameter :: max_string_length=16
	integer nk(max_string_length)
c----
	i=0
	if(nword0(n_len,txt,ib,ie)) then
	  ib=ie
	  return
	  endif

	if(len(label(n_lab)).gt.max_num_chars) then
	  write(6,*) '%EdPDB-F- error in the subroutine find1'
	  call exit(4)
	  endif
	
	l=min(ie-ib+1, len(label(1)))
	if(l.le.0) return
	k=0
	tmp=txt(ib:ie)
	do i=1,n_lab
	  if(tmp.eq.label(i)(:l)) then
	    if(index(label(i),' ').eq.l+1) then
	      return
	      endif
	    if(k.lt.max_string_length) k=k+1
	    nk(k)=i
	    end if
	  end do

	if(k.eq.1) then
	  i=nk(1)
	  return
	else if (k.gt.1) then
	  write(6,1001) (label(nk(j)),j=1,k)
	  i=0
	  return
	  end if
1001	format(' EdPDB-W> ambiguous command verb -- supply more characters'
	1/(8x,5a12))
100	i=-1
	end

	subroutine find(n_lab,label,i)
chk	===============
	include 'edp_dim.inc'
!	use edp_dim

	character*(*) label(n_lab)
	character*16  tmp

	character*(max_num_chars)  txt
	common /cmm_txt/n_len,txt,ib,ie
	logical nword

	i=0
	if(nword(n_len,txt,ib,ie)) return

	if(len(label(n_lab)).gt.max_num_chars) then
	  write(6,*) '%EdPDB-F- error in subroutine find'
	  call exit(4)
	  endif

	tmp=txt(ib:ie)
	do i=1,n_lab
	if(tmp.eq.label(i)) return
	end do
100	i=-1
	end

	function match_l1(n_lab,label)
chk	=================
	include 'edp_dim.inc'
!	use edp_dim

	character*(*) label(n_lab)
	character*16  tmp

	character*(max_num_chars)  txt
	common /cmm_txt/n_len,txt,ib,ie
	logical nword

	parameter (max_string_length=16)
!	integer, parameter :: max_string_length=16
	dimension nk(max_string_length)
c----
	match_l1=0
	if(nword(n_len,txt,ib,ie)) return

	if(len(label(n_lab)).gt.max_num_chars) then
	  write(6,*) '%EdPDB-F- error in  function match_l1'
	  call exit(4)
	  endif

	l=min(ie-ib+1,len(label(1)))
	if(l.le.0) return
	k=0
	tmp=txt(ib:ie)
	do i=1,n_lab
	if(tmp.eq.label(i)(:l)) then
	  if(index(label(i),' ').eq.l+1) then
	    match_l1=i
	    return
	  endif
	if(k.lt.max_string_length) k=k+1
	nk(k)=i
	end if
	end do
	if(k.eq.1) then
	 match_l1=nk(1)
	 return
	else if (k.gt.1) then
	write(6,1001) (label(nk(i)),i=1,k)
	match_l1=0
	return
	end if
1001	format(' EdPDB-W> ambiguous command verb -- supply more characters'/
	1    ,(8x,5a10))
100	match_l1=-1
	end

	function match_l(n_lab,label)
chk	================
	include 'edp_dim.inc'
!	use edp_dim

	character*(*) label(n_lab)
	character*16 tmp

	character*(max_num_chars)  txt
	common /cmm_txt/n_len,txt,ib,ie
	logical nword

	match_l=0
	if(nword(n_len,txt,ib,ie)) return

	if(len(label(n_lab)).gt.max_num_chars) then
	  write(6,*) '%EdPDB-F- error in the function match_l'
	  call exit(4)
	  endif

	tmp=txt(ib:ie)
	do i=1,n_lab
	  match_l=match_l+1
	  if(tmp.eq.label(match_l)) return
	  end do
100	match_l=-1
	end

	function match_id(n_lab,label)
chk	=================
	include 'edp_dim.inc'
!	use edp_dim

	character*(*) label(n_lab)
	character*16 tmp
	data id_save/0/

	character*(max_num_chars)  txt
	common /cmm_txt/n_len,txt,ib,ie
	logical nword

	if( nword(n_len,txt,ib,ie)) then
	  match_id=0
	else
	  if(len(label(n_lab)).gt.max_num_chars) then
	    write(6,*) '%EdPDB-F- error in the function match_id'
	    call exit(4)
	    endif
	  tmp=txt(ib:ie)
	  do match_id=1,n_lab
	    if( tmp .eq. label(match_id)) then
	      id_save= match_id
	      return
	    endif
	  enddo
	  match_id= -1
	  i= index( tmp, '+')
	  if( i .gt. 0) then
            ibi= ib+ i
	    read(txt(ibi:ie),*, err=900) j
	    if( i .eq. 1) then
	      match_id= id_save+ j
              if( match_id .gt. n_lab) match_id=-1
	    else
	      tmp=txt(ib:ibi-2)
	      do jj=1,n_lab
	        if( tmp .eq. label(jj)) then
	          id_save= jj
	          match_id= id_save+ j
                  if( match_id .gt. n_lab) match_id=-1
	          return
	        endif
	      enddo
	    endif
	  endif
	endif
900	end

	subroutine read_txt(io,*)
c	===================
c	get rid of the biginning '*' 
c	convert the string to lower case if 'tolower' is .true.
	include 'edp_dim.inc'
!	use edp_dim

	logical starting_input

	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib0,ie0

	common /main_para/ tolower
	1 ,jnk1(m_main_para-1)
	logical tolower	

	if(io.gt.0) then
	  n_len=0
	  read(io,'(a)',end=900) txt
	  n_len=ltrim(txt)
	  if(n_len.gt.max_num_chars) write(6,*) 
	1  '%EdPDB-W- the input line is truncated to'
	1 ,max_num_chars
	1 ,' chars.'
	  ie0=0
	  n_len=min(n_len,max_num_chars)
	  endif

	starting_input=.false.

	do ic=1,n_len
	  if(starting_input.or.
	1 txt(ic:ic).ne.'*'.and.txt(ic:ic).ne.' ') then 
	    starting_input=.true.
	    jc=ichar(txt(ic:ic))
	    if(tolower.and.(65.le.jc.and. jc.le.90)) txt(ic:ic)=char(jc+32)
	  else 	! ie. (.not.starting_input.and.txt(ic:ic).eq.'*') 
	    txt(ic:ic)=' '
	    endif
	  enddo
	return
900	return 1
	end

copyright by X. Cai Zhang


	subroutine init_group0()
chk	======================
	include 'edp_main.inc'
!	use edp_main

	if (n_groupa(0) .ne. n_atom) then
	  do j=1,n_atom
	    igroupa(j,0)=j
	    enddo
	  n_groupa(0)=n_atom
	  endif
	end

	subroutine chk_pipe()
chk	===================
	include 'edp_main.inc'
!	use edp_main
	include 'edp_dat.inc'
!	use edp_dat

        character*(max_num_chars) txt
        common /cmm_txt/n_len,txt,ib,ie

	common /main_para/ tolower, jou, echo_L, junk_i0, incl0, incl1,
	1 jjou, nest, n_length,			!9
	1 jnk_main_para(m_main_para-9)

	integer q1, q2

!	save j0,j1,j2,j3, n_length_prev
	save

	j1=index(txt,'|')
	if(j1.le.0) j1=9999
	j2=index(txt,';')
	if(j2.le.0) j2=9999
	j3=index(txt,'}')
	if(j3.le.0) j3=9999
	j0=min(j1,j2,j3)

	n_length_prev=n_length
	n_len0=index(txt(:n_length),'!')
	if(n_len0.le.0) n_len0=n_length

!000601
	q1=index(txt(:n_length),delimiter)
	if(q1.ge.1) then
	  q2=index(txt(q1+1:n_length),delimiter)
	else
	  q2=0
	  endif

	if(q1.lt.j0.and.j0.lt.q2+q1) then
	  i_pipe=-1
	else if(j0.gt.n_len0) then	!no  pipe
	  i_pipe=-1
	else if(index(txt(:n_length),'from').gt.0) then
	  i_pipe=-1		! from-prase and piping are incompatable.
	else			! pipe
	  i_pipe=1
	  call init_group0()
	  j3=index(txt,'{')
	  if(j3.le.0) j3=9999
	  j4=index(txt,'}')
	  if(j4.le.0) j4=9999
	  j0=min(j0,j3,j4)
	  n_len=j0-1
c	  if(n_len.le.0) goto 100
	  endif
	return

	entry piping(j_dollar_1)
chk	============
	n_delta=n_length-n_length_prev
	n_length_prev=n_length
	if(n_delta.ne.0.and. j0.gt. j_dollar_1) then
	 j0=j0+n_delta
	 j1=j1+n_delta
	 j2=j2+n_delta
	 j3=j3+n_delta
	endif

	n_len0=index(txt(:n_length),'!')
	if(n_len0.le.0) n_len0=n_length

100	continue
	if(j0.gt.n_len0) then
	  i_pipe=-1
	  return		! pipe completed
	  endif
	i_pipe=1

	do j=1,j0
	  txt(j:j)=' '
	  enddo
	if(j0.eq.j1) then	! '|': logical 'and'
	  j=0
	  do jj=1,n_atom
	    if(lf(jj)) then
	      j=j+1
	      igroupa(j,0)=jj
	      lf(jj)=.false.
	      endif
	    enddo
	  n_groupa(0)=j
	else if(j0.eq.j2) then	! ';': logical 'or'
	  call init_group0()
	else if(j0.eq.j3) then	! '{': copy ON to 0000 and initialize
	  j=0
	  do jj=1,n_atom
	    if(lf(jj)) then
	      j=j+1
	      igroupa(j,2)=jj
	      lf(jj)=.false.
	      endif
	    enddo
	  n_groupa(2)=j
          call init_group0
	else if(j0.eq.j4) then	! '}': copy back from 0000 to ON 
	  j0=n_groupa(2)
	  do j=1,j0
	    lf(igroupa(j,2))=.true.
	    enddo
	  endif

	j1=index(txt,'|')
	if(j1.le.0) j1=9999
	j2=index(txt,';')
	if(j2.le.0) j2=9999
	j3=index(txt,'{')
	if(j3.le.0) j3=9999
	j4=index(txt,'}')
	if(j4.le.0) j4=9999
	j0=min(j1,j2,j3,j4)

	if(j0.le.n_len0) then	
	  n_len=j0-1
	else 
	  n_len=n_len0
	  i_pipe=0
	  endif
	return
	end

chk**** end of find1.for file
