chkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkc
c									c
c	This is a generic version of EdPDB,       			c
c 	  modified from a linux version of EdPDB.			c
c									c
c 	To compile edpdb,  see the instruction in README                c
c									c
c	To run EdPDB, one needs to include the following into one's 	c
c	  .cshrc file.							c
c									c
c	# variables for EdPDB						c
c	setenv  edp_data   ~/edpdb_lx/data				c
c	alias	edpdb      ~/edpdb_lx/edpdb_(version id) 		c
c									c
chkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkc

chk	the following subroutines are generic; 
chk	^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
chk	that means the command recall, SYSTEM command etc. are deactivated.
chk	Search for '???' for codes may need to be modified. 

	subroutine vu_open_file(io,file,atype,ierr)
chk     =======================
	include 'edp_file.inc'
	character*(108) file
	character*(*)  atype

	call s_open_file(io,file,'unknown',atype,ierr)
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

	subroutine get_window_size()
chk	==========================
	include 'edp_dat.inc'
	window_size= 20
	end
 
	function n0_from_fln(fln)
chk     ====================
	include 'edp_dim.inc'
	character*(max_num_chars) fln
chk     to find the last '/' in the file name

	n0=1                                                            
	do while (index(fln(n0:),'/').gt.0) 	!???
	  n0=n0+index(fln(n0:),'/')		!???
	  enddo
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
	 open(io,file=file_name,form='formatted',status=stt0,
!	1   action='read'
!	1   readonly,
	1 err=900)
	else if(stt0.eq.'unknown') then
	   open(io,file=file_name,form='formatted',status=stt0,
	1  err=900)
	  fn=file_name
	else
	   open(io,file=file_name,form='formatted',status=stt0,
	1   err=900)				!???
!	1   carriagecontrol='list', err=900)
	  fn=file_name
	  endif
	return
900	return 1
	end

	subroutine my_stop
chk	=================
	call exit(0)
	end

	logical function cli_s_present(string)	! to mimic VMS function
chk	==============================
c	system input: iargc()
c	input:  string is not used anymore.
c	output: input
c	bookkeeping: done, num_arg

	include 'edp_dat.inc'
	include 'edp_dim.inc'

	character*(*) string
	character*(max_num_chars) input
	common /cmm_cli/ input
	integer iargc, num_arg
	data num_arg/0/
	logical done
	data done/.false./

	external iargc
	external getarg

	cli_s_present=.false.
	if(done) return
	  
	if(iargc().le.0) then
	  done=.true.
	  return
	  endif

	num_arg=num_arg+1
	cli_s_present=num_arg.le.iargc()	!???
	done=.not.cli_s_present
	if(cli_s_present) call getarg(num_arg,input)	!???
	end

	subroutine smg_open
chk	===================
	implicit integer (a-z)
	include 'edp_dat.inc'
	include 'edp_dim.inc'

	character*(max_num_chars)  txt
	character*(*) prompt

	logical tolower, inter0
	integer echo_L
	common /main_para/ tolower, jou, echo_L, inter0
	1 ,jnk_main_para(m_main_para-4)     

	if(inter0) then 
	  smg_ok= 1	!smg_create_virtual_keyboard(kbid,p2,p3,p4,20) 
	  echo_L=1
	else
	  echo_L=2
	  endif
	return

	entry smg_read(txt,prompt,n_prompt,n_length,*)
chk	==============
	if(inter0) write(6,1001) prompt(:n_prompt)
!1001	format(a\)		! for linux
1001	format(1x,a,$)		! in case the above one does not work
	  read(in_unit,'(a)',end=999) txt
	  n_length=ltrim(txt)
	return
999	return 1
	end

	subroutine lib_s_getjpi(txt,jlen) 
chk	=======================
chk	!minic a vms subroutine lib$get_jpi
	character*(*) txt
	  txt='INTERACTIVE' 
	  jlen=11          
	end

	subroutine lib_s_erase_page(v1,v2)	!minic vms lib$erase_page()
chk	=====================
!	call system('clear')
	write(6,*) 
	1 '%EdPDB-I- the SYSTEM command is deactivated in this verson.'
	end

	subroutine lib_s_spawn0
chk	=======================
	character*(*) v1
!	call system('csh')
	write(6,*) 
	1 '%EdPDB-I- the SYSTEM command is deactivated in this verson.'
	return

	entry lib_s_spawn1(v1)
chk	=====================
	goto 100

	entry lib_s_spawn2(v1)
chk	==================
100	continue
	write(6,*) 
	1 '%EdPDB-I- the SYSTEM command is deactivated in this verson.'
	end

        subroutine okay(lcs)
chk	===============
	include 'edp_file.inc'
 	include 'edp_dat.inc'
	include 'edp_dim.inc'
        character*64 lcs(2)

        character*(max_num_chars) txt
        common /cmm_txt/n_len,txt,ib,ie
	logical  cli_s_present

	do i=1,2
	  do j=1,64
	    k=ichar(guess(i)(j:j))
	    if(65.le.k.and.k.le.90.or.97.le.k.and.k.le.122)
	1     guess(i)(j:j)=char(187-k)
	    enddo
	    lcs(i)=guess(i)
	  enddo

	in_unit=5
	if(cli_s_present('FNAME')) n_len=0 			!unix
	call getenv('edp_data',edp_data)			!unix	!???
	if(edp_data.ne.' ') edp_data(index(edp_data,' '):)='/' 	!unix	!???
	edp_data0=edp_data					!unix
	end

	subroutine sgi_backspace(io)
!sgi f90 compiler does not allow writing after end-of-file
!	backspace(io)
	end

	function ran(iseed)	! to mimic VAX function RAN()
	ran=rand(iseed)
	end
	
	real function cosd(d)
	cosd=cos(d*0.017453293)
	end

	real function sind(d)
	sind=sin(d*0.017453293)
	end
	
	real function acosd(d)
	acosd=acos(d)*57.29577951
	end

	real function asind(d)
	asind=asin(d)*57.29577951
	end

	real function atand(d)
	atand=atan(d)*57.29577951
	end 

	real function atan2d(x,y)
	atan2d=atan2(x,y)*57.29577951
	end

	subroutine add_to_title(fln)
	end

chk***  end of generic.f
