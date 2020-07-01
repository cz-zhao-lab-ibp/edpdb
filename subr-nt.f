chkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkc
c								                                      c
c	This is a version of EDPDB for NT computer,                       c
c 	  modified from a unix version of EDPDB.			              c
c									                                  c
c 	To compile edpdb,  see the instruction in readme.txt              c
c									                                  c
c	To run EDPDB, one needs to include the following into one's 	  c
c	  .cshrc file.							                          c
c									                                  c
c	# variables for EDPDB						                      c
c	setenv  edp_data   ~/edpdb/data					                  c
c	setenv  hlplib     $edp_data/hlplib				                  c
c	setenv  edphlp     $edp_data/edpdb_v95a.hlp			              c
c	alias	edpdb      ~/edpdb/edpdb_v96a				              c
c									                                  c
chkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkc

chk	the following subroutines are NT specific.
chk	^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	subroutine vu_open_file(io,file,def_type,ierr)
chk     =======================
	include 'edp_file.inc'
!	use edp_file

	character*(108) file
	character*(*)  def_type

c-vms	call s_open_file(io,file,'new',type,ierr)
	call s_open_file(io,file,'unknown',def_type,ierr)
	do i=1,108
	  if( 
	1  (pdb_in(i:i).eq.' '.or.pdb_in(i:i).eq.'.') 
	1   .and.
	1  (file(i:i).eq.' '.or.file(i:i).eq.'.') 
	1  ) then
	    ierr=4
	    close (io)
	    return
	  else if(pdb_in(1:i).ne.file(1:i)) then
	    return
	    endif
	  enddo
	end

 
	function n0_from_fln(fln)
chk     ====================
	include 'edp_dim.inc'
!	use edp_dim

	character*(max_num_chars) fln
!    for Window NT, find the last '\' in the file name
!    for unix, find the last '/' in the file name
!    for vax, find the last ']' in the file name

	n0=1 
!    for Window NT                                                           
	do while (index(fln(n0:),'\') > 0)
	  n0=n0+index(fln(n0:),'\')
	  enddo

!    for UWIN (Window version of unix)
	if(n0 == 1 ) then
	 do while (index(fln(n0:),'/') > 0)
	  n0=n0+index(fln(n0:),'/')
	 enddo
	endif 

	n0_from_fln=n0
	end
 
	subroutine open_file0(io,file_name,stt,dft,fn,*)
chk	=====================
	character*(*) stt
	character*(*) dft
	character*(108) file_name
	character*(*)  fn
	character*(32)  stt0

	if(stt.eq.'new_or_unknown') then
	  stt0='unknown'
	else
	  stt0=stt
	  endif

	if(stt0.eq.'old') then
	 open(io,file=trim(file_name),form='formatted',status=stt0,
	1   readonly, err=900)
	else if(stt0.eq.'unknown') then
	   open(io,file=trim(file_name),form='formatted',status=stt0,
	1  err=900)
	  fn=trim(file_name)
	else
	   open(io,file=trim(file_name),form='formatted',status=stt0,
	1   carriagecontrol='list', err=900)
c	 if(index(file_name,' ')-1.le.len(fn)) then
	  fn=trim(file_name)
	  endif
	return
900	return 1
	end

	subroutine my_stop
chk	=================
	read(*,*)
	call exit(0)
	end

	logical function cli_s_present(string)	! to mimic VMS function
chk	==============================
c	system input: iargc()
c	input:  string is not used anymore.
c	output: input
c	bookkeeping: done, num_arg

	include 'edp_dat.inc'
!	use edp_dat
	include 'edp_dim.inc'
!	use edp_dim

	character*(*) string
	character*(max_num_chars) input
	common /cmm_cli/ input
	integer num_arg
	data num_arg/0/
	logical done
	data done/.false./
	integer*2 istatus

 	cli_s_present=.false.
	if(done) return
	
	if(iargc_lc().le.0) then
	  done=.true.
	  return
	  endif

	num_arg=num_arg+1
	cli_s_present= (num_arg <= iargc_lc())
	done=.not.cli_s_present
	if(cli_s_present) call getarg_lc(num_arg,input,istatus)
	end

	integer function iargc_lc()
chk	================
	use DFLIB
	use DFPORT
	include 'edp_dat.inc'
!	use edp_dat
	include 'edp_dim.inc'
!	use edp_dim

	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie

	character*(max_num_chars) input
	integer*2 istatus,jj
	character*(72) argc(8)
	data argc/8*' '/

	integer jargc
	data jargc/-1/

	iargc_lc=jargc
	if(jargc >= 1) return

	jargc=min(8,iargc())
	if(jargc >= 1) then
	  do j=1,jargc
	    jj=j
	    call getarg(jj,argc(j), istatus)
	  enddo
	  iargc_lc=jargc
	  return
	endif

	call change_dir('cwd.txt')
	write(6,1002)
1002	format('$%edpdb ')
	read(5,'(a)') txt
	n_len=ltrim(txt)
	
	jargc=0
	ib=0
	ie=0
	do while(.true.)
	if(.not.nword(n_len,txt,ib,ie)) then
	  if(jargc >= 8) then
	    write(6,*)'argc-W> too many input items'
	    call my_stop()
	    endif
	  jargc=jargc+1
	  argc(jargc)=txt(ib:ie)
	else
	  ib=0
	  ie=0
	  iargc_lc=jargc
	  return
	endif
	end do
	return

	entry igetarg(i,input)
	igetarg=0
	input=argc(i)
	end

	subroutine getarg_lc(i,input,istatus)
	include 'edp_dim.inc'
!	use edp_dim
	character*(max_num_chars) input

	integer*2 istatus

	istatus=igetarg(i,input)
	end

	subroutine smg_open
chk	===================
	use DFLIB
	implicit integer (e,s)

	include 'edp_dat.inc'
!	use edp_dat
	include 'edp_dim.inc'
!	use edp_dim

	character*(max_num_chars)  txt
	character*(*) prompt
	save kbid

	logical tolower
	integer echo_L
	common /main_para/ tolower, jou, echo_L, inter0
	1 ,jnk_main_para(m_main_para-4)     

	if(inter0) then 
	  smg_ok= 1	!smg_create_virtual_keyboard(kbid,,,,20) 
	  if( smg_ok.ne.1) then
	    write(6,*) '%EDPDB-F- smg_open: error'
            call exit(4)
	    endif
	  echo_L=1
	else
	  echo_L=2
	  endif
	return

	entry smg_read(txt,prompt,n_prompt,n_length,*)
chk	==============
	if(inter0) then 
c	  smg_ok= smg_read_composed_line
c	1     (kbid,,txt,prompt,n_prompt,n_length,)
c	  if(smg_ok.ne.1) then
c	    write(*,*) '%EDPDB-F- read_smg: error.'
c	    call exit(4)
c	    endif
	write(6,1001) prompt(:n_prompt)
	n_length=GETSTRQQ(txt)
1001	format(a\)
c1001	format(1x,a\)
	else
	  read(in_unit,'(a)',end=999) txt
	  n_length=ltrim(txt)
	  endif
	return
999	return 1
	end

	subroutine lib_s_getjpi(txt,jlen) 
chk	=======================
chk	!minic a vms subroutine lib$get_jpi
	character*(*) txt

chk	use "% man fstat" to check fstat 
	integer fstat
	integer statb(12)

c	i=fstat(5,statb)
c	if(i.ne.0) then
c	  write(6,1001) i
c1001	format(' %EDPDB-F- stdin (channel 5) open failure (code:',i2,')')
c	  call exit(4)
c	  endif

c	if(statb(8).eq.0) then 
	  txt='INTERACTIVE'
	  jlen=11
c	else
c	  txt='BATCH'
c	  jlen=5
c	  endif
	end

	subroutine lib_s_erase_page(v1,v2)	!minic vms lib$erase_page()
chk	=====================
c	call system('clear')
	write(6,*) 'this command is deactivated in this NT version'
	end

	subroutine lib_s_spawn0
chk	=======================
! this subroutine does not work on PC
	use DFLIB
   	character*(*) v1
	INTEGER(4) I, errnum

!	I = SYSTEM('')
!	If (I == -1) then
!		errnum = ierrno( )
!		print *, 'Error ', errnum
!	end if
!	type *,'%EDPDB-W- currently, the system subroutine does not work.'
	return

	entry lib_s_spawn1(v1)
chk	=====================
	goto 100

	entry lib_s_spawn2(v1)
chk	==================
100	continue
	if(ltrim(v1).gt.0) then
	I = SYSTEM(v1)
	If (I .eq. -1) then
		errnum = ierrno( )
		print *, 'Error ', errnum
	end if
	else
!	write(*,*) i
	If (I .eq. -1) then
		errnum = ierrno( )
		print *, 'Error ', errnum
	end if
	  endif
	I = SYSTEM('')
!	write(*,*) 
!	1'%EDPDB-W- currently, the system subroutine does not work.'
	end

	subroutine okay(lcs)
chk	===============
	use DFLIB
	include 'edp_file.inc'
!	use edp_file
	include 'edp_dat.inc'
!	use edp_dat
	include 'edp_dim.inc'
!	use edp_dim

	character*(64) lcs(2)
	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie
	logical  cli_s_present

	integer(2) oldcur, istat

	do i=1,2
	  do j=1,64
	    k=ichar(guess(i)(j:j))
	    if(65.le.k.and.k.le.90.or.97.le.k.and.k.le.122)
	1     guess(i)(j:j)=char(187-k)
	    enddo
	    lcs(i)=guess(i)
	  enddo

	call get_dir(edp_data)
	if(index(edp_data,'/') > 0 ) then		!for UWIN (Window version of unix)
	 edp_data= trim(edp_data)//'data/'
	else !	if(index(edp_data,'\') > 0 ) then ! for WindowNT
	 edp_data= trim(edp_data)//'data\'
	endif

	edp_data0=edp_data
				
	call set_window()

	in_unit=5
	if(cli_s_present('FNAME')) n_len=0 			!unix, nt
	end

	subroutine set_window()
!	=====================
! USE statements to include routine interfaces
	use dflib
	use dfwin
	include 'edp_main.inc'
!	use edp_main

! Data declarations
!	Type(T_CONSOLE_SCREEN_BUFFER_INFO) cinfo
	Type(T_COORD) wpos
	Type(T_SMALL_RECT) sr
	Type(T_CONSOLE_CURSOR_INFO) cursor

	integer fhandle
	logical lstat
	character*(64) title
	character*(40) fln

	fhandle = GetStdHandle(STD_OUTPUT_HANDLE)

! Executable code to set console window size
	sr.top    =   0
	sr.left   =   0
	sr.bottom =   39 ! <= console buffer height -1
	sr.right  =   89 ! <= console buffer width  -1

! Executable code to set console buffer size
	wpos.x = 90		! columns >= console window width
	wpos.y = 1024	! lines   >= console window height
	lstat = SetConsoleScreenBufferSize(fhandle, wpos)
 
! Executable code to set console window size
	lstat = SetConsoleWindowInfo(fhandle, .true., sr)

	cursor.DwSize = 96
	cursor.bVisible = .true.
	lstat = SetConsoleCursorInfo(fhandle, cursor)

	title='edpdb '//version
	lstat = SetConsoleTitle(title)
	return

	entry add_to_title(fln)
!	==================
	title='EdPDB ('//version//')--'//fln
	lstat = SetConsoleTitle(title)
	end

	subroutine get_window_size()
!	==========================
! USE statements to include routine interfaces
	use dflib
	use dfwin
	include 'edp_dat.inc'
!	use edp_dat

! Data declarations
	integer fhandle
	logical lstat

	Type(T_CONSOLE_SCREEN_BUFFER_INFO) conbuf
	type (T_COORD)        dwSize
      type (T_SMALL_RECT)   srWindow

	fhandle = GetStdHandle(STD_OUTPUT_HANDLE)

! Executable code to get console buffer size
	lstat = GetConsoleScreenBufferInfo(fhandle, conbuf)
	window_size= max(1,(conbuf.srWindow.bottom-conbuf.srWindow.top)-2)	
c	 write (*,*) "window size= ",window_size
c	 write (*,*) "Window coordinates= ", conbuf.srWindow
c	 write (*,*) "Buffer size= ", conbuf.dwSize

	end

	subroutine get_dir(dir)
chk	=====================	
	use DFLIB
	character*(108) program_name
	character*(*) dir

	call getarg(0,program_name)
	n0=n0_from_fln(program_name)-1
	program_name=program_name(1:n0)
	if(index(program_name,'Debug\') >= 1 ) n0=n0-6
	if(index(program_name,'Debug/') >= 1 ) n0=n0-6
	dir=program_name(1:n0)
	end

	subroutine change_dir(file)
chk	=====================	
	use DFPORT

	integer(4) istatus
	character*(108) newdir, cwd, cwd0, program_dir
	character*(*) file

	logical  open_file
	external open_file

	ISTATUS = GETCWD (cwd0)
	call get_dir(program_dir)

	if(open_file(7,trim(program_dir)//trim(file),'unknown','') )then
!	  read(7,*,err=100, end=100)
	 do
	  read(7,'(a)',err=100, end=100) cwd
	  if(cwd(1:11) == 'setenv cwd ') cwd=trim(cwd(12:))
	  if( chdir(cwd) == 0) goto 10
	 enddo
	endif

100	cwd=cwd0
	if( chdir(cwd) /= 0) cwd='????'
      
10	write(6,1001) trim(cwd) 
1001	format(' %default directory name [',a,'] : ',$)
	read(5,'(a)') newdir
	j_nd=ltrim(newdir)
	if(j_nd <= 0 ) newdir=cwd
	if(newdir(1:4) == '????') newdir=cwd0

	ISTATUS = CHDIR(newdir)
	select case (istatus)
       case (enoent)
		write(6,*) 'The directory'//trim(newdir)//' does not exist'
		goto 10
       case (enotdir)
 		write(6,*) trim(newdir)//' is not a directory'
		goto 10
       case default
 		if(j_nd >= 1 ) 
	1	write(6,*)'Default directory successfully changed.'
		write(6,*)
	    rewind(7)
		write(7,'(a)') 
	1	'!This file stores/updates the name of working-dirctory.'
	    write(7,'(a)') 'setenv cwd '//trim(newdir)
		close(7)
	end select
	end

	subroutine sgi_backspace(io)
! sgi f90 compiler does not allow writing after an end-of-file. 
	end

chk***  end of nt_subrtn.for
