chkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkc
c                                                                       c
	program EdPDB                                                    !c
c       Written in otc. 1987 when the author worked                     c
c         at the Institute of Biophysics, Academia Sinica, China        c
c       Modified extensively in 1990-1993,                              c 
c         at Institute of Molecular Biology,                            c
c         University of Oregon, USA                                     c
c                                                                       c
c       Copyright by                                                    c
c                                               X. Cai Zhang, Ph.D.     c
c                                                                       c
chkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkc

chk	unit  1 pdb_file.pdb (1st)                  (i)
chk	unit  3 pdb_file.pdb (2nd)                  (i)
chk	unit  2 new_file.edp 	                    (o)
chk	unit  4 new_file.pdb 	                    (o)
chk	unit  7 sequ_file.seq 	and	
chk	        sequ_alignment block file         	(i/o)
chk	unit  8 'edp_data/acc','edp_data/pdbstd'	(o)
chk	unit  9 pdb_file.pdb (after the 2nd) 		(i)
chk	unit 11--18 preprepared.edp                 (i)
chk	unit 19 symmetry .edp file                  (i)
chk	unit 21 new .edp file                       (o)
chk	unit 29 old rtnmtx                          (i)
chk	unit 22 new rtnmtx                          (o)
chk	unit 34 new kinemage vector list file (no longer used)
chk	unit 48 scr_file.out                        (o)
chk	unit in_unit	sys$input                   (i)

!	verbose 0 shut up
!			1 normal, including to list .out file
!			2 describe output in more details
!			3 warning (-W-), information (-I-) and command syntax
!			4 describe input line interpretation, talky output
!			5 more talky output
!			6 list references, more talky output
  
	include 'edp_main.inc'
!	use edp_main
	include 'edp_file.inc'
!	use edp_file
	include 'edp_dat.inc'
!	use edp_dat

	logical  cli_s_present
	character*(max_num_chars) input_line
	common /cmm_cli/ input_line
 
	character*8 s_quit
	character*(64) lcs(2)

	data s_quit/ 'quit'/

	character*12 vp

	integer np, npp
	character*8  pp
	character*72 pn
	common /ud_param/ npp(max_param), pp(max_param),
	1 np(max_param), pn(max_param)
  
	logical	nword
        character*(max_num_chars) txt,fln
        common /cmm_txt/n_len,txt,ib,ie

        common /cmm_symm/ num_symm
	1 ,sym0(3,4,max_symm) ,sym(3,4,max_symm)

	character*(108) file
	logical inter0, tolower, incl0, incl1
	integer echo_l
!	character*8 blink,normal	!010322
	character*1 blink,normal
	character*2 junk_bl
	character*12 prompt
	character*4 err_control

	common /main_para/ tolower, jou, echo_L, inter0, incl0, incl1,
	1 jjou, nest, n_length,
	1 blink,normal, junk_bl, prompt,n_prompt,
	1 err_control, file
c	1 ,jnk_main_para(m_main_para-17)     

	logical making_edp
	character*(108) tmp_edp
	common /create_edp/ making_edp, tmp_edp

	character*(120) fmtstmt
	character pdate*12,ptime*10
	character*132 txt11

	logical  open_file
	external open_file

chk	initialize
	do i=1,max_reserved_param
	  if(i.lt.10) then
	    npp(i)=2
	    write(pp(i),1034) i
	  else
	    npp(i)=3
	    write(pp(i),1035) i
	    endif
1034	  format('p',i1)
1035	  format('p',i2)
	  enddo

	do i= 1,max_param
	  np(i)=0
	  enddo

	do i=1,max_vector
	  cvector(i)=' '
	  enddo

	call date_and_time(pdate,ptime)

chk	initialization
	call okay(lcs)
	lcs(1)(11:14)=version

	pdb_in='?'  ;  pdb_in2='?' ; pdb_in3='?'; pdb_out='?'
	rtn_in='?'  ;  rtn_out='?' ; edp_in='?';  edp_out='?'
	scr_out='?' ;  acc_in='?'  ; std_in='?';  seq_out='?'
	update_in='?'; log_out='?'
		
	call get_window_size()
	if(window_size.le.1) window_size=window_size0
	delimiter=char(39)		! character '
	wildcard ='*'
	wildcard1='%'
	num_symm=-1
	n_prompt=1
	prompt='*'
	n_err=0
	max_err=1024
	err_control='quit'
	n_of_syn=0
	cgroup(1)='scr'
	cgroup(2)='0000'
	verbose= 1
	status=0
	making_edp=.false.
	zone_untouch=.true.; atom_untouch=.true.; residue_untouch=.true.
	call set_parameter(2,'sg',2,'p1',2)
	call set_parameter(3,'edp_data',8,edp_data,ltrim(edp_data))


chk	open input pdb file.
	tolower=.false.
	txt=input_line
	n_len=ltrim(txt)
	call read_txt(0,*999)
	ie=0
	if(nword(n_len,txt,ib,ie)) goto 100
	ib0= ib
	fln=input_line(:ltrim(input_line))
	call add_to_title(fln)
	n_fn=ltrim(input_line)

	n0=n0_from_fln(fln)

	n1=index(fln(n0:),'.')
	if(n1.le.0) then
	  n1=n_fn
	else
	  n1=n0+n1-2
	  endif
	pdb_in=fln(:n_fn)
	if(.not.open_file(1,pdb_in,'old','.pdb')) then
	  write(6,*) 
	1'%EdPDB-E- failure in opening the pdb file ['//
	1 pdb_in(:ltrim(pdb_in))//']'
	  goto 101
	  endif
	file=fln(n0:n1)

	cm1='_'
	cm2=' '
	ncm1=1
	ncm2=0
	if(.not. cli_s_present('more input?')) goto 110

	tolower=.true.
	txt=input_line
	n_len=max_num_chars
	call read_txt(0,*999)

	if(txt(1:2).eq.'-i') goto 110
	if(txt(1:2).eq.'-c') then
	  if(.not.cli_s_present('chain names?')) goto 101 ! wrong input '-c <eol>'

	tolower=.true.
	txt=input_line
	n_len=max_num_chars
	call read_txt(0,*999)

	  n_len=72
	  ie=0

	  if(nword(n_len,txt,ib,ie)) goto 101

	  n_len=ie
	  call read_txt(0,*999)
	  cm1=txt(ib:ie)
	  ncm1=ie-ib+1
	  call set_parameter(4,'cm1',3,cm1,ncm1)
	  if(.not.cli_s_present('pdb file2?')) goto 110
	  endif

chk	open file 2 if existing
	tolower=.false.
	txt=input_line
	n_len=max_num_chars
	call read_txt(0,*999)
 
	ie=0
	if(nword(n_len,txt,ib,ie)) goto 100
	n_len=ie
	call read_txt(0,*999)
	if(txt(1:2).eq.'-i')  goto 110

	pdb_in2=txt(:n_len)
	if(.not. open_file(3,pdb_in2, 'old','.pdb')) then
	  write(6,*) 
	1'%EdPDB-E- failure in opening the pdb file ['//
	1 pdb_in2(:ltrim(pdb_in2))//']'
	   goto 101
	   endif

	cm2='_'
	ncm2=1
	if(.not. cli_s_present('more input?')) goto 110

	tolower=.true.
	txt=input_line
	n_len=max_num_chars
	call read_txt(0,*999)

	if(txt(1:2).eq.'-i') goto 110
	if(txt(1:2).eq.'-c') then
	  if(.not.cli_s_present('chain names?')) goto 101 ! wrong input '-c <eol>'

	tolower=.true.
	txt=input_line
	n_len=max_num_chars
	call read_txt(0,*999)

	  n_len=72
	  ie=0
	  if(nword(n_len,txt,ib,ie)) goto 101
	  n_len=ie
	  call read_txt(0,*999)
	  cm2=txt(ib:ie)
	  ncm2=ie-ib+1
	  call set_parameter(5,'cm2',3,cm2,ncm2)
	  if(.not.cli_s_present('edp file?')) goto 110
	  endif

chk	open .edp history file
110	edp_out=file
	if(.not. open_file( 2,edp_out,
	1 'new_or_unknown',
	1 '.edp')) then
	  write(6,*) 
	1 '%EdPDB-F- failure in opening a new history file ['//
     1 edp_out(:ltrim(edp_out))//']'
	  goto 101
	  endif

chk	open scratch file
	scr_out=file
	if(.not. open_file(48,scr_out,
	1 'new_or_unknown','.out')) then
	  write(6,*) 
	1 '%EdPDB-F- failure in opening a new scratch file ['//
	1 scr_out(:ltrim(scr_out))//']'
	  write(6,*) 
	  call my_stop() !exit( 4)
	  endif

	if(cm2.eq.' ') then
	  if(cm1.eq.'_') then
	    write( 2,1062)  lcs, pdb_in(:ltrim(pdb_in))
	1 ,' ',' ',' ',' ',' ' 
 	    write(48,1062)  lcs, pdb_in(:ltrim(pdb_in))
	1 ,' ',' ',' ',' ',' ' 
          else 
	    write( 2,1062)  lcs, pdb_in(:ltrim(pdb_in))
	1 ,'-c',cm1(:ltrim(cm1)),' ',' ',' '
	    write(48,1062)  lcs, pdb_in(:ltrim(pdb_in))
	1 ,'-c',cm1(:ltrim(cm1)),' ',' ',' '
	  endif
	else
	  if(cm1.eq.'_') then
	    if(cm2.eq.'_') then
	      write( 2,1062)  lcs, pdb_in(:ltrim(pdb_in))
	1 ,' ',' ',pdb_in2(:ltrim(pdb_in2)),' ',' ' 
	      write(48,1062)  lcs, pdb_in(:ltrim(pdb_in))
	1 ,' ',' ',pdb_in2(:ltrim(pdb_in2)),' ',' ' 
	    else
	      write( 2,1062)  lcs, pdb_in(:ltrim(pdb_in))
	1 ,' ',' ',pdb_in2(:ltrim(pdb_in2)),'-c',cm2(:ltrim(cm2)) 
	      write(48,1062)  lcs, pdb_in(:ltrim(pdb_in))
	1 ,' ',' ',pdb_in2(:ltrim(pdb_in2)),'-c',cm2(:ltrim(cm2)) 
	    endif
          else 
	    if(cm2.eq.'_') then
	      write( 2,1062)  
	1 lcs, pdb_in(:ltrim(pdb_in)),'-c',cm1(:ltrim(cm1))
	1 ,pdb_in2(:ltrim(pdb_in2)),' ',' ' 
	      write(48,1062)  
	1 lcs, pdb_in(:ltrim(pdb_in)),'-c',cm1(:ltrim(cm1))
	1 ,pdb_in2(:ltrim(pdb_in2)),' ',' ' 
	    else
	      write( 2,1062)  
	1 lcs, pdb_in (:ltrim(pdb_in )),'-c',cm1(:ltrim(cm1))
	1     ,pdb_in2(:ltrim(pdb_in2)),'-c',cm2(:ltrim(cm2)) 
	      write(48,1062)  
	1 lcs, pdb_in (:ltrim(pdb_in )),'-c',cm1(:ltrim(cm1))
	1     ,pdb_in2(:ltrim(pdb_in2)),'-c',cm2(:ltrim(cm2)) 
	    endif
	  endif
	endif 
	
 	write( 2,1064) 
	1 ptime(1:2),ptime(3:4),ptime(5:6), pdate(1:4),pdate(5:6),pdate(7:8)
 	write(48,1064) 
	1 ptime(1:2),ptime(3:4),ptime(5:6), pdate(1:4),pdate(5:6),pdate(7:8)

1062	format(2('!',a/),
	1 '!% edpdb ',6(a,' '))
1064	format('! at ',a,':',a,':',a,' on ',a,'-',a,'-',a)

	n_fn=n_len-ib0+1
	fln= txt(ib0:n_len)
	iadd=1
	do igr=0, max_gr
	  n_groupa(igr)=0
	  enddo
	call set_parameter(1,'file1',5,file,n1-n0+1)
	errmsg=' errmsg: invalid parameter'

	call lib_s_getjpi(txt,jlen)
	inter0=txt(:jlen).eq.'INTERACTIVE'

!chk	open smg$xxxx system service	
	call smg_open
	inter=inter0
	if(inter) then							!010319 
!	  blink =char(7)//char(27)//'[5m'
!	  normal=char(27)//'[0m'
	  blink=char(7)
	  normal=' '
	else
	  blink =' '
	  normal=' '
	  write(6,*) 
	  endif

	jou=10
	tolower=.true.
	txt=input_line
	n_len=max_num_chars
	call read_txt(0,*999)

	if(txt(1:2).eq.'-i') then !'-i'

	  if(.not.cli_s_present('edp file?')) then
	    edp_in='edpini.edp'
	  else

	tolower=.false.
	txt=input_line
	n_len=max_num_chars
	call read_txt(0,*999)

	    edp_in=txt(:ltrim(txt))
	    endif

	  jou=11
	  if(.not. open_file(11, edp_in, 'old','.edp')) then
	    if(.not. open_file(11, 
	1 edp_data(:ltrim(edp_data))//edp_in,
	1 'old','.edp')) then
	      write(6,*) 
	1'%EdPDB-W- the file ['//
	1 edp_in(:ltrim(edp_in))//'] does not exist.'
	      close(11)
	      jou=10
	      endif
	    endif
	  endif !'-i'

	tolower=.true.
	txt=' '
	n_len=0	  

2	incl0=.true.
	incl1=.true.
	call readf(ierr)
	if( ierr.ne.0) then
	  write(6,*) '%EdPDB-F- error in reading the pdb file(s)'		
	  write(6,*)
	  call my_stop() !exit(4)
	  endif
1001	  format(a)

149	n_length= max(1,min(max(n_length,n_len), max_num_chars))
	write(2,1001)	txt(:n_length)
	write(48,1001)	'!'//txt(:n_length)
	if( 
	1  ((jou.gt.10 .and. echo_L .gt. 1) 
	1  .and. (txt(1:1) .eq. '!' .and. txt(1:2) .ne. '!!') 
	1  .and. (verbose .ge. 1) 
	1  ) .or. (verbose .ge. 12)
	1  ) write(6,'(a)') ' '//prompt(:n_prompt)//txt(:n_length)

150	continue
	if(jou.gt.10) then
	  read(jou,'(a)',end=200) txt
	  n_length=ltrim(txt)
	else	! input unit:sys$input (in_unit)
	  call smg_read(txt,prompt,n_prompt,n_length,*999)
	  endif

	if(n_length.le.0) then 			! "return for more ...."
	  call c_list
	  call c_typelist
	  goto 150
	  endif
	call e_list
	call e_typelist

	n_length= max(1,min(n_length, max_num_chars))
	n_len=n_length
	
	call read_txt(0,*999)	! convert to lower case

chk	!parameter substitution
	if(.not.making_edp .and. index(txt(1:n_len),'$(').gt.0) then

	  j_substitution=0
501	  if(j_substitution .gt. max_param ) then	
	    errmsg=' errmsg: recursive parameter substitution.'
	    goto 305
	    endif
	  do jp=1,max_param
	    npjp=np(jp)
	    if(npjp.gt.0) then
	      vp='$('//pp(jp)(1:npp(jp))//')'
	      ips=npp(jp)+3
	      ib   =index(txt,vp(1:ips))
	      if(ib.gt.0) then
	        ie   =min(max_num_chars,max(n_length,ib+ips))
		n_length=min(max_num_chars,n_length-ips+npjp)
		txt(ib:n_length)=pn(jp)(1:npjp)//txt(ib+ips:ie)
	        j_substitution=j_substitution + 1
		goto 501
	  	endif
	      endif
	    enddo

	  endif

	  n_length= max(1,min(n_length, max_num_chars))
	  n_len=index(txt(:n_length),'!') -1 
	  if(n_len.eq.0) goto 149	! starting from '!'
	  if(n_len.lt.0) n_len=n_length


	if(jou.le.10) then	! from the terminal
	  if(echo_l.ge.2 ) write(6,1001) ' !'//txt(:n_length)
	else			! from a macro
	  if(echo_l.ge.1) write(6,1001) 
	1   ' '//prompt(:n_prompt)//txt(:n_length)
	  endif

	nest=0	
	call get_command(*300,*305)
	goto 150

100	write(6,*) '%EdPDB-F- file open failure'
101	write(6,*)
	write(6,*) '%EdPDB-I- version: '//version
	write(6,*) '%EdPDB-I- usage:   edpdb pdb_file '
	write(6,*) '          [-c chainnames] '//
	1'[2nd_pdb_file [-c chainnames]]  [-i edp_file]'
	call my_stop() !exit(4) 

200	close (jou) 
	jou=jou-1
	goto 150

300	if(status.eq.0) status=-1
	fmtstmt='(<n>x,a:/5x,a)'
	write(fmtstmt(2:4),'(i3)') min(max(1,(ib+ie)/2+n_prompt),79)
	write(6,fmtstmt) blink//'^'//normal//' invalid command'
	goto 307

305	if(status.eq.0) status=-2
	if (echo_L.eq.0.and.jou.gt.10) 
	1 write(6,1001) ' '//prompt(:n_prompt)//txt(:n_length)
	do i=60,1,-1
	  if( errmsg(i:i) .ne. ' ') goto 306
	  enddo
!1019	format(<min(max(1,(ib+ie)/2+n_prompt-2),79)>x,a:/5x,a)
306	fmtstmt='(<n>x,a:/5x,a)'
	write(fmtstmt(2:4),'(i3)') min(max(1,(ib+ie)/2+n_prompt),79)
	write(6,fmtstmt) blink//'^'//normal,errmsg(:i)
	if(verbose.ge.3 ) then
        if(syntax(1)(1:7) .ne. 'syntax:') write(6,*) 
	1 '%EdPDB-I3> syntax:'
	  do i=1,n_of_syn
	     write(6,*) syntax(i)(:ltrim(syntax(i)))
	     enddo
	  endif
	n_of_syn=0
	errmsg=' errmsg: invalid parameter'
	backspace (2)
307	n_err=n_err+1
	if(n_err.gt.max_err) then
	  if(jou.le.10) then

308	    write(6,*) '%EdPDB-E- abort for too many errors'
	    write(6,*)
	    call my_stop() !exit(4)
	  else if(err_control.eq.s_quit) then
	    write(6,*) '%EdPDB-E- abort for too many errors'
	    write(6,*)
	    call my_stop() !exit(4)
	    endif
	  goto 200
	else if(max_err-n_err.le.max(10,max_err/20) ) then
	  write(6,*) '%EdPDB-W- ',blink//'reset maxerr'//normal
	1 ,'(e.g. setenv maxerr)'
	  endif 
	goto 150

999	end


	subroutine error_handling_2151(*)
c	A dirty fix to make pipe work for multiple commands 
	include 'edp_main.inc'
!	use edp_main
	include 'edp_file.inc'
!	use edp_file
	include 'edp_dat.inc'
!	use edp_dat

	logical	nword
        character*(max_num_chars) txt,fln
        common /cmm_txt/n_len,txt,ib,ie
	character*(120) fmtstmt

	character*8 s_quit
	data s_quit/ 'quit'/

	character*(108) file
	logical inter0, tolower, incl0, incl1
	integer echo_l
!	character*8 blink,normal	!010322
	character*1 blink,normal
	character*2 junk_bl
	character*12 prompt
	character*4 err_control

	common /main_para/ tolower, jou, echo_L, inter0, incl0, incl1,
	1 jjou, nest, n_length,
	1 blink,normal, junk_bl, prompt,n_prompt,
	1 err_control, file
c	1 ,jnk_main_para(m_main_para-17)     

305	if(status.eq.0) status=-2
	if (echo_L.eq.0.and.jou.gt.10) 
	1 write(6,1001) ' '//prompt(:n_prompt)//txt(:n_length)
	do i=60,1,-1
	  if( errmsg(i:i) .ne. ' ') goto 306
	  enddo
!1019	format(<min(max(1,(ib+ie)/2+n_prompt-2),79)>x,a:/5x,a)
306	continue
	if(verbose.ge.3 ) then
	fmtstmt='(<n>x,a:/5x,a)'
	write(fmtstmt(2:4),'(i3)') min(max(1,(ib+ie)/2+n_prompt),79)
	write(6,fmtstmt) blink//'^'//normal,errmsg(:i)
        if(syntax(1)(1:7) .ne. 'syntax:') write(6,*) 
	1 '%EdPDB-I3> syntax:'
	  do i=1,n_of_syn
	     write(6,*) syntax(i)(:ltrim(syntax(i)))
	     enddo
	  endif
	n_of_syn=0
	errmsg=' errmsg: invalid parameter'
c	backspace (2)
307	n_err=n_err+1
	if(n_err.gt.max_err) return 1
1001	format(a) 
	end

chk**** end of edpdb.for file
