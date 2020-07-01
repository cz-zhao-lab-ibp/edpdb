	subroutine true_or_false(edp_if,*)
chk	========================
	include 'edp_main.inc'
!	use edp_main
	include 'edp_dat.inc'
!	use edp_dat

	logical edp_if

	integer np, npp
	character*8  pp
	character*72 pn
	common /ud_param/ npp(max_param), pp(max_param), 
	1 np(max_param), pn(max_param)

!	character*1  pn1
!	character*72 input

	logical tolower
	common /main_para/ tolower
	1 ,jnk(m_main_para-1)

	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie
	logical nword

	parameter (n_items=5)
	character*8 items(n_items)
	data items/'==', '^=', '>=', '<=','/='/

	parameter (n_opts=3)
	character*8 p_opts(n_opts)
	data  p_opts/'-s', '-r', '-i'/

	n_of_syn=2					!000515
	syntax(1)='syntax:' 
	syntax(2)=
	1'if ( [-(s,r,i)] parameter_name.s operator.s value.s ) command.s' 

	ie0=ie
	n_len0=n_len

	i0=index( txt, '(')
	if(i0.le.0)  goto 900
	ie=i0
	i1=index( txt, ')')
	if(i1.le.i0)  goto 900
	n_len=i1-1

	call find1(n_opts, p_opts, i_opt)
	if(i_opt.lt.0) then
	  ie=i0
	  i_opt=1
	  endif

	i_par= match_l( max_param, pp)
	if( i_par .lt. 0) goto 900

	call find1(n_items, items, i_log)
	if(i_log.le.0) goto 900

	edp_if=.false.
	if( nword(n_len,txt,ib,ie) ) goto 900
	if( i_opt.eq.1) then	! '-s: string'
c	  copy from setcm.for		
	  ic=ichar(delimiter)
	  if(ic.eq. ichar(txt(ib:ib)) ) then	! using ' ' to input parameter
	    ib=ib+1
	    do j=ib, n_len
	      if(ichar(txt(j:j)).eq.ic) then
	        txt(j:j)=' ' 	!000303, it may screw up the history file
		goto 2890
	        endif
	      enddo
	    j=n_len+1
2890	    ie=j-1
	    endif

	  if( i_log.eq.1) then	! '=='
	    if(pn(i_par)(:np(i_par)) .eq. txt(ib:ie)) edp_if=.true.
	  else if( i_log.eq.2.or.i_log.eq.5) then	! '^=' or '/=' 
	    if(pn(i_par)(:np(i_par)) .ne. txt(ib:ie)) edp_if=.true.
	  else if( i_log.eq.3) then	! '>='
	    if(pn(i_par)(:np(i_par)) .ge. txt(ib:ie)) edp_if=.true.
	  else if( i_log.eq.4) then	! '<='
	    if(pn(i_par)(:np(i_par)) .le. txt(ib:ie)) edp_if=.true.
            endif
	else if( i_opt.eq.2) then	! '-r'
	    read(pn(i_par)(:np(i_par)),*,err=900) rv
	    read(txt(ib:ie),*,err=900) rv1
	    if( i_log.eq.1) then
	      if( rv.eq. rv1) edp_if=.true.
	    else if( i_log.eq.2) then
	      if( rv.ne. rv1) edp_if=.true.
	    else if( i_log.eq.3) then
	      if( rv.ge. rv1) edp_if=.true.
	    else if( i_log.eq.4) then
	      if( rv.le. rv1) edp_if=.true.
	      endif
	  else if( i_opt.eq.3) then	! '-i'
	    read(pn(i_par)(:np(i_par)),*,err=900) iv
	    read(txt(ib:ie),*,err=900) iv1
	    if( i_log.eq.1) then
	      if( iv.eq. iv1) edp_if=.true.
	    else if( i_log.eq.2) then
	      if( iv.ne. iv1) edp_if=.true.
	    else if( i_log.eq.3) then
	      if( iv.ge. iv1) edp_if=.true.
	    else if( i_log.eq.4) then
	      if( iv.le. iv1) edp_if=.true.
	      endif
	  else 
	    goto 900
	  endif 

	ie=n_len+1
	n_len=n_len0
	return 

900	n_len=n_len0
	ie=ie0
	return 1
	end

	subroutine param(*,*)
chk	=================
chk	return	: normal
chk	return 1: error
chk	return 2: some wrong, but may be okey

	include 'edp_main.inc'
!	use edp_main
	include 'edp_dat.inc'
!	use edp_dat

	integer np, npp
	character*8  pp
	character*72 pn
	common /ud_param/ npp(max_param), pp(max_param), 
	1 np(max_param), pn(max_param)

	character*1  pn1
	character*72 input

	character*6 sexit, sreport
	data sexit,sreport/'exit','report'/

	parameter (n_items=10)
!	integer, parameter :: n_items=10
	character*8 items(n_items)
	data items/'entry', 'atom', 'residue', 'chain', 'id', 
	1 'x', 'y', 'z', 'w', 'b'/

	character*(64) prompt_string

	logical tolower
	common /main_para/ tolower
	1 ,jnk(m_main_para-1)

	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie
	logical nword, junk

	n_of_syn=9					!000515
	syntax(1)='syntax:' 
	syntax(2)=
	1'1) parameter [pn.s]' 
	syntax(3)=
	1'2) parameter pn.s = [value.s]' 
	syntax(4)=
	1'3) parameter pn.s ? [prompt_string.s] [default_value.s] '
	syntax(5)='    [(exit, report) ]' 
	syntax(6)=
	1'4) parameter pn.s (+, -) step_size.i limit.i [(exit, report)]' 
	syntax(7)=
	1'5) parameter pn.s group_id.s '
	syntax(8)=
	1'    (entry, atom, residue, chain, id, x, y, z, w, b)'
	syntax(9)='    [(exit, report)]' 

	iv= match_l( max_param, pp)
	if( iv .lt. 0) then
	  do iv=max_reserved_param+1, max_param
	    if(np(iv).le.0) then
	      pp(iv)=txt(ib:ie)
	      npp(iv)=ie-ib+1
	      goto 1280
	      endif
	    enddo
	  errmsg=' errmsg: Too many paramters have been defined.'
	  return 1 
	else if(iv.eq.0) then
	  do iv=1,max_param
	    if(np(iv).gt.0) then
	      write(6,1089) pp(iv),pn(iv)(:np(iv))
	      endif
	    enddo
	  return  
	  endif

1280	if( nword(n_len,txt,ib,ie) ) then		! parameter P1
	  if(np(iv).gt.0) then
	    write(6,1089)  pp(iv),pn(iv)(:np(iv))
	  else
	    write(6,1089)  pp(iv),'null'
	    endif
1089	    format(' param> ',a8,' is defined as ',a)
	else if(txt(ib:ie).eq.'=') then			! parameter p1 =
	  if( nword(n_len,txt,ib,ie))  then		! parameter p1 = <cr>
	    np(iv)=0
	    return
	    endif

c	  copy from setcm.for		
	  ic=ichar(delimiter)				! parameter p1 = something
	  if(ic.eq. ichar(txt(ib:ib)) ) then	! using ' ' to input parameter
	    ib=ib+1
	    do j=ib, n_len
	      if(ichar(txt(j:j)).eq.ic) then
		txt(j:j)=' '	!000303, it may screw up the history file 
		goto 2890
	        endif
	      enddo
	    j=n_len+1
2890	    ie=j-1
	    endif
	  np(iv)=ie-ib+1
	  pn(iv)=txt(ib:ie)

	else if(txt(ib:ie).eq.'?') then			! parameter p1 ? ...
	  if( nword(n_len,txt,ib,ie)) then		! parameter p1 ? ,
	    prompt_string='input '
	    ns=6
	  else if( txt(ib:ib) .eq. delimiter ) then	! parameter P1 ? 'something'
	    ib=ib+1
	    do j=ib, n_len
	      if( txt(j:j) .eq. delimiter ) goto 2896
	      end do
	    j=n_len+1
2896	    ie=j-1
	    prompt_string=txt(ib:ie)
	    ns=ie-ib+1
	  else						! parameter P1 ? something
	    prompt_string=txt(ib:ie)
	    ns=ie-ib+1
	    endif

	  if( .not.nword(n_len,txt,ib,ie) ) then	! parameter p1 ? something default
	    if( txt(ib:ib) .eq. delimiter ) then	! parameter p1 ? something 'default'
	      ib=ib+1
	      do j=ib, n_len
	        if( txt(j:j) .eq. delimiter ) goto 2897
	        end do
	      j=n_len+1
2897	      ie=j-1
	      endif
	    prompt_string(ns+1:)='['//txt(ib:ie)//'] '
	    ns=ns+ie-ib+4
	    pn(iv)=txt(ib:ie)
	    np(iv)=ie-ib+1
	  else if(np(iv).gt.0) then			! parameter p1 ? something ,
	    prompt_string(ns+1:)='['//pn(iv)(:np(iv))//'] '
	    ns=ns+np(iv)+3
	    endif
							! parameter P1 ? string default
	  call smg_read(input,prompt_string,min(64,ns),n_input,*2892)

	  if(n_input.gt.0) then				! got something
	    np(iv)=n_input
	    do ic=1,n_input
              jc=ichar(input(ic:ic))
              if(tolower.and.(65.le.jc.and.jc.le.90)) then
	        pn(iv)(ic:ic)=char(jc+32)
	      else
	        pn(iv)(ic:ic)=char(jc)
	        endif
	      enddo
	    endif

	else if(txt(ib:ie).eq.'+') then			! parameter p1 + ...
	  iv1=np(iv)
	  pn1=pn(iv)(1:1)
	  if( (pn1.lt.'0'.or.pn1.gt.'9') 
	1 .and.pn1.ne.'-'
	1 .and.pn1.ne.'+') goto 2101
	  read( pn(iv)(:iv1),*,err=2101) iv2
	  if( nword(n_len,txt,ib,ie)) goto 2151
	  read( txt(ib:ie),*,err=2151) iv3
	  if( nword(n_len,txt,ib,ie))  goto 2151
	  read( txt(ib:ie),*,err=2151) iv4
	  iv2=iv2+iv3
	  if(iv2.gt.iv4) goto 2892
	  write( pn(iv),*) iv2
	  ie_v=0
	  if( nword(72,pn(iv),ib_v,ie_v)) goto 2151
	  np(iv)=ie_v-ib_v+1
	  pn(iv)(:np(iv))=pn(iv)(ib_v:ie_v)
	  return
2101	  iv2=ichar(pn(iv)(1:1))
	  if( nword(n_len,txt,ib,ie)) goto 2151
	  read( txt(ib:ie),*,err=2151) iv3
	  if( nword(n_len,txt,ib,ie))  goto 2151
	  iv4=ichar(txt(ib:ib))
	  iv2=iv2+iv3
	  if(iv2.gt.iv4) goto 2892
	  pn(iv)(1:1)=char(iv2)

	else if(txt(ib:ie).eq.'-') then			! parameter p1 - ...
	  iv1=np(iv)
	  pn1=pn(iv)(1:1)
	  if( (pn1.lt.'0'.or.pn1.gt.'9') 
	1 .and.pn1.ne.'-'
	1 .and.pn1.ne.'+') goto 2102
	  read( pn(iv)(:iv1),*,err=2102) iv2
	  if( nword(n_len,txt,ib,ie)) goto 2151
	  read( txt(ib:ie), *, err=2151) iv3
	  if( nword(n_len,txt,ib,ie)) goto 2151
	  read( txt(ib:ie), *, err=2151) iv4
	  iv2=iv2-iv3
	  if(iv2.lt.iv4) goto 2892
	  write( pn(iv),*) iv2
	  ie_v=0
	  if( nword(72,pn(iv),ib_v,ie_v)) goto 2151
	  np(iv)=ie_v-ib_v+1
	  pn(iv)(:np(iv))=pn(iv)(ib_v:ie_v)
	  return
2102	  iv2=ichar(pn(iv)(1:1))
	  if( nword(n_len,txt,ib,ie)) goto 2151
	  read( txt(ib:ie),*,err=2151) iv3
	  if( nword(n_len,txt,ib,ie))  goto 2151
	  iv4=ichar(txt(ib:ib))
	  iv2=iv2+iv3
	  if(iv2.gt.iv4) goto 2892
	  pn(iv)(1:1)=char(iv2)

	else				! parameter p1 group_id
	  ie=ib-1
	  igr=match_l(max_gr,cgroup)
	  if(igr.gt.0) then
	    ng=n_groupa(igr)
	  else
	    ng=0
	    endif
	  if(ng.le.0) then
	    junk=nword(n_len,txt,ib,ie)
	    goto 2892
	  else
	    ia=igroupa(1,igr)
	    ir=aa_seq(ia)
	    call find1(n_items, items, i)
	    if(i.le.0) goto 2151
	    if( i.eq. 1) then		! entry
	      pn(iv)=text(ia)(8:11)
	      np(iv)=4
	    else if( i .eq. 2) then	! atom
	      pn(iv)=atom(ia)
	      np(iv)=4
	    else if( i .eq. 3) then	! residue
	      pn(iv)=res(ir)
	      np(iv)=3
	    else if( i .eq. 4) then	! chain_mark
	      pn(iv)=text(ia)(22:22)
	      np(iv)=1
	    else if( i .eq. 5) then	! id
	      pn(iv)=ires(ir)
	      np(iv)=5
	    else if( i .eq. 6) then	! x
	      pn(iv)=text(ia)(31:38)
	      np(iv)=8
	    else if( i .eq. 7) then	! y
	      pn(iv)=text(ia)(39:46)
	      np(iv)=8
	    else if( i .eq. 8) then	! z
	      pn(iv)=text(ia)(47:54)
	      np(iv)=8
	    else if( i .eq. 9) then	! w
	      write(pn(iv),'(f6.2)',err=2891) w(ia)
	      np(iv)=6
	    else if( i .eq. 10) then	! b
	      write(pn(iv),'(f6.2)',err=2891) b(ia)
	      np(iv)=6
	    else 
	      return 1
	      endif
	    n_groupa(igr)=n_groupa(igr)-1
	    if(ng.gt.1) then
	      do i=2,ng
	        igroupa(i-1,igr)=igroupa(i,igr)
	        enddo
	      endif
	    endif 
	  endif 
	return

2891	errmsg=' errmsg: error during read the numberic value'
	return 1
2892	status= -3 
	if( nword(n_len,txt,ib,ie)) then
	 errmsg= ' errmsg: parameter offlimit'
	else if(index(sreport,txt(ib:ie)).eq.1 ) then
	 errmsg= ' errmsg: parameter offlimit'
	else if(index(sexit,txt(ib:ie)).eq.1 ) then
	 return 2
	endif
2151	return 1
	end 

c1004	 format(q,a)
c***	end of param.for

copyright by X. Cai Zhang
