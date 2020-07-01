chkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkc
c									c
c	This is a version of EdPDB for SGI unix computer, 		c
c 	  modified from an openVMS version of EdPDB.			c
c									c
c 	To compile edpdb,  see the instruction in readme.txt            c
c									c
c	To run EdPDB, one needs to include the following into one's 	c
c	  .cshrc file.							c
c									c
c	# variables for EdPDB						c
c	setenv  edp_data   ~/edpdb/data					c
c	setenv  hlplib     $edp_data/hlplib				c
c	setenv  edphlp     $edp_data/edpdb_v95a.hlp			c
c	alias	edpdb      ~/edpdb/edpdb_v96a				c
c									c
chkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkchkc

chk	the following subroutines are unix specific.
chk	^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	
	subroutine vu_open_file(io,file,type0,ierr)
chk     =======================
	include 'edp_file.inc'
!	use edp_file
	character*(108) file
	character*(*)  type0

c-vms	call s_open_file(io,file,'new',type,ierr)
	call s_open_file(io,file,'unknown',type0,ierr)
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
!	use edp_dat
	integer lines, cols
	external lines, cols
c	logical first_run
c	data first_run/.true./

!	write(*,*) lines(), cols()
c	if(first_run) then
		window_size= max(window_size0,lines()-2)
c		first_run=.false.
c	else
c		window_size= max(1,lines()-2)
c	endif
	end
 
	function n0_from_fln(fln)
chk     ====================
	include 'edp_dim.inc'
!	use edp_dim
	character*(max_num_chars) fln
chk     to find the last '/' in the file name

	n0=1                                                            
	do while (index(fln(n0:),'/').gt.0)
	  n0=n0+index(fln(n0:),'/')
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
!	1   action='read', err=900)
	1   readonly, err=900)
	else if(stt0.eq.'unknown') then
	   open(io,file=file_name,form='formatted',status=stt0,
	1  err=900)
	  fn=file_name
	else
	   open(io,file=file_name,form='formatted',status=stt0,
	1   err=900)
!	1   carriagecontrol='list', err=900)
c	 if(index(file_name,' ')-1.le.len(fn)) then
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
!	use edp_dat
	include 'edp_dim.inc'
!	use edp_dim
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
	cli_s_present=num_arg.le.iargc()
	done=.not.cli_s_present
	if(cli_s_present) call getarg(num_arg,input)
	end

	subroutine smg_open
chk	===================
	implicit integer (e,p,s)
	include 'edp_dat.inc'
!	use edp_dat
	include 'edp_dim.inc'
!	use edp_dim
	character*(max_num_chars)  txt
	character*(*) prompt
	save kbid

	logical tolower, inter0
	integer echo_L
	common /main_para/ tolower, jou, echo_L, inter0
	1 ,jnk_main_para(m_main_para-4)     

	if(inter0) then 
	  smg_ok= smg_create_virtual_keyboard(kbid,p2,p3,p4,20) 
	  if( smg_ok.ne.1) then
	    write(6,*) '%EdPDB-F- smg_open: error'
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
	  write(6,*)	! for sgi-unix
	  smg_ok= smg_read_composed_line
	1     (kbid,p2,txt,prompt,n_prompt,n_length)
	  if(smg_ok.ne.1) then
	    write(6,*) '%EdPDB-F- read_smg: error.'
            call exit(4)
	    endif
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

	i=fstat(5,statb)
	if(i.ne.0) then
	  write(6,1001) i
1001	format(' %EdPDB-F- stdin (channel 5) open failure (code:',i2,')')
	  call exit(4)
	  endif

	if(statb(8).ne.0) then 
	  txt='BATCH'
	  jlen=5
	else
	  txt='INTERACTIVE' !if "fstat" does not work, keep thest two statements
	  jlen=11           ! and deactivate everything else in this subroutine.
	  endif
	end

	subroutine lib_s_erase_page(v1,v2)	!minic vms lib$erase_page()
chk	=====================
	call system('clear')
	end

	subroutine lib_s_spawn0
chk	=======================
	character*(*) v1
	call system('csh')
	return

	entry lib_s_spawn1(v1)
chk	=====================
	goto 100

	entry lib_s_spawn2(v1)
chk	==================
100	if(len(v1).gt.0) then
	  call system(v1)
	else
	  call system('csh')
	  endif
	end

        subroutine okay(lcs)
chk	===============
	include 'edp_file.inc'
!        use edp_file
 	include 'edp_dat.inc'
!       use edp_dat
	include 'edp_dim.inc'
!        use edp_dim
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
	call getenv('edp_data',edp_data)			!unix
	if(edp_data.ne.' ') edp_data(index(edp_data,' '):)='/' 	!unix
	edp_data0=edp_data					!unix
	end

	subroutine sgi_backspace(io)
!sgi f90 compiler does not allow writing after end-of-file
	backspace(io)
	end

	subroutine add_to_title(fln)
	end
chk***  end of unixsubrtn.for
