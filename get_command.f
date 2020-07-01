	subroutine get_command(*,*)
!return 1: invalid command
!return 2: invalid parameter

	include 'edp_main.inc'
!	use edp_main
	include 'edp_file.inc'
!	use edp_file
	include 'edp_dat.inc'
!	use edp_dat

	character*8 s_save,s_exit,s_quit,s_rms,s_scale,s_wait,s_morph
	character*6 show
	character*28 cm_1, cm_2, cm_3

	character*(108) file1
	character*(64) goto_label

	integer np, npp
	character*8  pp
	character*72 pn
	common /ud_param/ npp(max_param), pp(max_param), np(max_param),
	1 pn(max_param)

	data s_exit,s_quit,s_save,s_rms,s_scale, s_wait, s_morph !&
	1 /'exit','quit','save','rms','scale','wait','morph'/

	logical	nword, nword0
	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie

	common /cmm_symm/ num_symm, jnk1(3,4,max_symm) ,jnk2(3,4,max_symm)

	character*(108) file
	logical inter0, tolower, incl0, incl1
	integer echo_L 
	character*1 blink,normal, junk_bl*2
	character*12 prompt, prompt_save
	save prompt_save

	character*4 err_control, istatus

	common /main_para/ tolower, jou, echo_L, inter0, incl0, incl1, !&
	1 jjou, nest, n_length,	!&
	1 blink,normal, junk_bl, prompt,n_prompt,     	!&
	1 err_control, file		
!   ,jnk_main_para(m_main_para-17)

	logical making_edp
	character*(108) tmp_edp
	character*(64)  eof, pm_inline, pm_inline0
	common /create_edp/ making_edp, tmp_edp, eof, neof
	data j_sub/0/

	integer scan_times
	logical jerr, edp_if
	real mw

	character*30 edp_data_save
	integer echo_L_save, verbose_save, window_size_save, max_err_save, n_err_save
	logical tolower_save, inter_save, save_all, notSaveYet
	character*(1)  delimiter_save, wildcard_save, wildcard1_save 
	character*(4) err_control_save

	data notSaveYet/.true./	
	save edp_data_save
	1,echo_L_save, verbose_save, window_size_save, max_err_save, n_err_save
	1,tolower_save, inter_save, save_all
	1,delimiter_save, wildcard_save, wildcard1_save 
	1,err_control_save

	character*12 stterm(10)
	data stterm/'verbose','edp_data','window_size','echo','tolower', !&
	1 'interactive','delimiter','wildcard','prompt','maxerr'/

!	'define'
	character*8 dfterm(8)
	data dfterm/'ca','main','ab','abc','abcd','bridge', 'newxyz','residue'/

	logical  open_file
	external open_file
	logical  open_file1
	external open_file1
	integer  threading, w2t
                                     
	parameter	(n_lab=134)
	character*16 label(n_lab) 
	data label/'atom','residue','initialize','exit','quit',				!5
	1 'cell','symmetry','append','close','zone','main',				!11
	1 'side','ca','include','exclude','list','help','write','analyze',		!19
	1 'setb','setw','b','w','x','y','z','distance','nayb',	     	    		!28
	1 'reset','abc','abcd','naybr','rtn','dfmain','dfca',				!35
	1 'dfab','dfabc','dfabcd','ab','shdf','dfbrg','bridge','clear','pause',		!44
	1 'swap','more','clique','dfres','switchwb','group','code',			!51	
	1 'load','setc','rewind','interactive','seta','avb','blank',			!58
	1 'newxyz','dfnewxyz','setr','seti','sett','sequence','sumw','file',		!66
	1 'set','read','mmi','mmir','system','axis','jiggle','sete',			!74
	1 'correlation','access','mmig','diff','overlay','rmsw','prompt',		!81
	1 'match','sort','ratio','planar','match3d','align3d',
	1 'permute','parameter','maxerr',
	1 'pattern','shape','volume','alias','movecenter','chain','harker',		!97
	1 'euler','polar','vm','extract','momentinertia','update',			!103
	1 'closer','setenv','operator','lattice','text','return','goto',		!110
	1 'if','vector','mkfile','mw','match1d','seq2pdb','copy',			!117
	1 'sdistance','snayb','snaybr','define','threading','w2t','doit',		!124
	1 'pdist','poverlay','pclique','pdoit','psdist','sclique','eigen',		!134
	1 'remark','dock','touch'/

!copy right by Xj Cai Zhang
	
	if(making_edp) then
	  if(txt(1:neof).ne.eof(1:neof)) then
	    write(21,'(a)') txt(:n_length)
	  else
	    close(21)
	    making_edp=.false.
	    endif
	    goto 150
	  endif

	edp_if=.false.
	n_length= max(1,min(max(n_length,n_len), max_num_chars))

	if(nest.le.0) then
	  incl=incl0
	  incl1=incl0
	else
	  incl=.true.
	  endif

20	i_plus=index(txt,'--')
	if(i_plus.gt.0) then
	  txt(i_plus:i_plus+1)=' +'
	  goto 20
	  endif
	call chk_pipe()
	if(i_pipe.ge.0) then
	  write(2,1001)  txt(1:n_length)
	  write(48,1001) '!!'//txt(1:n_length)
	  endif
! get the inline parameter, i.e. the text string after the last '$=',
! if exists in the original input line
	j_dollar=0
	k_dollar=index(txt(:n_length),'$=')
	do while (k_dollar.ge.1)
	 j_dollar=j_dollar+k_dollar
	 k_dollar=index(txt(j_dollar+1:n_length),'$=')
	enddo 

	if(j_dollar.gt.0) then
	  pm_inline0=txt(j_dollar+2:n_length)
	  if( index(pm_inline0,'|').le.0
 	1 .and.index(pm_inline0,';').le.0
 	1 .and.index(pm_inline0,'{').le.0
 	1 .and.index(pm_inline0,'}').le.0    ) then
	   j_sub0=ltrim(pm_inline0)
	   if(j_sub0.gt.0) then
	    txt(j_dollar:n_length)=' '
	    if(pm_inline0(1:1).eq.delimiter) then
	     j20=index(pm_inline0(2:j_sub0),delimiter)
		 if(j20.gt.0) then
	      pm_inline0=pm_inline0(2:j20)
		  j_sub0=j20-1
	     endif
	    endif
	  endif
	 else
	  j_sub0=0
	 endif
	else
	  j_sub0=0	
	endif

!	user-defining-keyword
!	ierr=999 -- empty
!	ierr=0   -- substitution
!	ierr<0   -- define udk
!	ierr>0   -- non udk

	loop_1=0
502	continue
	loop_1=loop_1+1
	if(loop_1.gt.20) return 2	! safety loop

	call udk(n_length,ierr, blink//'WARNING:'//normal)

	if(ierr.eq.999) then
	  ib=max(1,ib)
	  if(txt(ib:ib).eq.',') return 1 ! error: a line starting w/ ','
	  goto 150

	else if(ierr.eq.0) then   !substitution
	  ie=0
	  goto 502

	else if(ierr.lt.0) then	! define keyword
	  i=-1
	  goto 2300
	  
	else if(ierr.gt.1) then 
	  return 2
	endif

	ie=0
	ib=0
503	call find1(n_lab,label,i)
	
	if(i.eq.0) then			 ! i=0 : no word
	  if(ib.gt.ie)	return 1 ! error: a line starting w/ ','
	  if(ib.le.0)	return
	  if(txt(ib:ib).eq.'!') then
	    write(2,1001)  txt(:n_length)
	    write(48,1001) '!!'//txt(:n_length)
	    endif
	  return			! ambiguous command

	else if(i.lt.0) then			! i<0 : non-edpdb-keyword
	  if(txt(ib:ib).eq.'!') return
	  if(txt(ie:ie).eq.':') goto 150
	  if(nest.gt.0) return 1
	  if(txt(ib:ib).eq.'.'.or.txt(ib:ib).eq.'@') then
	    if(i_pipe.gt.0) then
	      errmsg=' errmsg: illegal piping order'
	      goto 2151
	      endif
	    if(ib.eq.ie) then		! there is no file name
	      errmsg=' errmsg: file open failure'
	      goto 2151
	      endif
	    if(jou.ge.18) then
	      write(6,*) '%EdPDB-F- ',blink//'^'//normal,
	1	  ' too many .edp type files'
	      if(verbose.ge.3 ) 
	1	  write(6,*)'%EdPDB-I3- currently max_num_edp_files= 18'
	      write(6,*)
	      call exit( 4)		! exit to prevent infinite loop
	      endif
	    file1 =txt(ib+1:ie)
	
!	    defind parameter $(p1),$(p2), .... $(p8)
	    num_new_p=0
505	    if(.not.nword0(n_len,txt,ib,ie)) then
	      num_new_p=num_new_p+1
	      if(ib.gt.ie) goto 505
	      if(num_new_p.gt.max_reserved_param) then
	        write(6,*)
	1 '%EdPDB-W- too many parameters. only p1 - p8 will be used.'
	        goto 511
	        endif
	      if(ib .le. ie) then
	        ic=ichar(delimiter)
	        if(ichar(txt(ib:ib)) .eq. ic) then	! starting  with '
	          ib=ib+1
	          do j=ib, n_len
	            if(ichar(txt(j:j)).eq.ic) then
	              txt(j:j)=' '	!000303, it may screw the history file
	              goto 510
	              endif
	            enddo
	          j=n_len+1
510	          ie=j-1
	          endif
	        pn(num_new_p)=txt(ib:ie)
	        np(num_new_p)=ie-ib+1
	        if(txt(ib:ie).eq.
	1		'$('//pp(num_new_p)(1:npp(num_new_p))//')') goto 2151
	        endif
	      goto 505
	      endif

511	    jou=jou+1
	    edp_in=file1
	    if(.not.open_file(jou,edp_in,'old','.edp')) then
	      edp_in=edp_data(:ltrim(edp_data))//edp_in
	      if(.not.open_file(jou,edp_in,'old','.edp')) then
	edp_in=edp_data(:ltrim(edp_data))//'symmetry/'//
	1 edp_in(ltrim(edp_data)+1:)
	       if(.not.open_file(jou,edp_in,'old','.edp')) then
	        write(6,*) '    ',blink//'^'//normal,'file open failure'
	        jou=jou-1
	        endif
	       endif
	      endif
	    goto 150
	    endif	! @file_name
	  if (index(txt(:n_len), '???').ge.1) goto 3017
	  return 1	! invalid command
	  endif		! i < 0

2300	continue
! get the inline parameter if exists in a (user-defined) command.
	j_dollar=index(txt(:n_len),'$=')
	if(j_dollar.gt.0) then
	  i_dollar=n_len
	  k_dollar=index(txt(j_dollar+1:n_len),';')
	  if(k_dollar.gt.0) i_dollar=min(i_dollar,j_dollar+k_dollar-1)
	  k_dollar=index(txt(j_dollar+1:n_len),'|')
	  if(k_dollar.gt.0) i_dollar=min(i_dollar,j_dollar+k_dollar-1)
	  k_dollar=index(txt(j_dollar+1:n_len),'$=')
	  if(k_dollar.gt.0) i_dollar=min(i_dollar,j_dollar+k_dollar-1)
	  pm_inline=txt(j_dollar+2:i_dollar)
	  j_sub=ltrim(pm_inline)
	  if(j_sub.gt.0) then
	    txt(j_dollar:i_dollar)=' '
		if(pm_inline(1:1).eq.delimiter) then
	     j20=index(pm_inline(2:j_sub),delimiter)
		 if(j20.gt.0) then
	      pm_inline=pm_inline(2:j20)
		  j_sub=j20-1
	     endif
	    endif
	   endif
	else
	  pm_inline=pm_inline0
	  j_sub=j_sub0
	endif
	j_1st_dollar=n_length
	if(j_sub.gt.0) then
	 do j_loop=1,32
	  j_dollar=index(txt(:n_len),'$')
	  if(j_dollar.le.0) goto 2301
	  j_1st_dollar=min(j_1st_dollar,j_dollar)
	  txt(j_dollar+j_sub:)=txt(j_dollar+1:)
	  txt(j_dollar:j_dollar+j_sub-1)=pm_inline(:j_sub)
	  n_len=n_len+j_sub-1
	  n_length=n_length+j_sub-1
	 enddo
	endif

	
2301	continue

!	i>0 deal with edpdb-keyword
!	in case of save,quit,rewind,return, goto 	!,recall 
	if(i.eq.4 .or. i.eq.5.or.i.eq.54.or. i.eq.109.or.i.eq.110) then 
	  write(2,1001)	'!!'//txt(:n_length)
	else
	  if(i_pipe.lt.0 .and. nest.le.0 .and. .not.edp_if)  then 
	    write(2,1001)	txt(:n_length)
	  else 
	    write(2,1001)	'!!'//txt(:n_length)
	    endif
	  endif
	write(48,1001)	'!!'//txt(:n_length)
!	if(verbose.ge.4 .and. echo_L < 2 ) write(6,'(a)') ' ;'//txt(1:n_length)
	if(i.le.0) goto 150
	
	n_of_syn=0
	write(istatus,'(i4)') status
	call set_parameter(6,'status',3,istatus,4)
	status=0
	
	goto (  2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,
	1	2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,
	1	2021,2022,2023,2024,2025,2026,2027,2028,2029,2030,
	1	2031,2032,2033,2034,2035,2036,2037,2038,2039,2040,
	1	2041,2042,2043,2044,2045,2046,2047,2048,2049,2050,
	1	2051,2052,2053,2054,2055,2056,2057,2058,2059,2060,
	1	2061,2062,2063,2064,2065,2066,2067,2068,2069,2070,
	1	2071,2072,2073,2074,2075,2076,2077,2078,2079,2080,
	1	2081,2082,2083,2084,2085,2086,2087,2088,2089,2090,
	1	2091,2092,2093,2094,2095,2096,2097,2098,2099,2100,
	1	2101,2102,2103,2104,2105,2106,2107,2108,2109,2110,
	1	2111,2112,2113,2114,2115,2116,2117,2118,2119,2120,
	1	2121,2122,2123,2124,2125,2126,2127,2128,2129,2130,          
	1	2131,2132,2133,2134)  i          

149	write(6,1149) 
1149	format(60x,'<<< dull keyword')
150	continue
	if (i_pipe.lt.0) return	
	call piping(j_1st_dollar)
	if (i_pipe.lt.0) return
	incl=incl0
	goto 502

2151	continue
	if (i_pipe.lt.0) goto 2153 
	call error_handling_2151(*2153)
	call piping(j_1st_dollar)
	if (i_pipe.lt.0) return 
	incl=incl0
	goto 502
2153	return 2 


2152	if(syntax(1)(1:7) .ne. 'syntax:') write(6,*) 'syntax:'
      do i=1,n_of_syn
	  write(6,*) syntax(i)(:ltrim(syntax(i)))
	enddo
	goto 150	

2001	syntax(1)='1) ATOM' 
	syntax(2)='2) ATOM atom_1.s [atom_2.s ... ]' 
	syntax(3)='3) ATOM EXCEPT atom_1.s [atom_2.s ... ]'
	n_of_syn=3; if (index(txt(:n_len),'???').ge.1) goto 2152
	 
	call pick(*2151)
	goto 150

2002	syntax(1)='1) RESIDUE'
	syntax(2)='2) RESIDUE res_type1.s [res_type2.s ... ]' 
	syntax(3)='3) RESIDUE EXCEPT res_type1.s [res_type2.s ... ]'
	n_of_syn=3; if (index(txt(:n_len),'???').ge.1) goto 2152
	 
	call pickr(*2151)
	goto 150

2003	syntax(1)='INITIALIZE '
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	if(.not.nword(n_len,txt,ib,ie)) goto 2151
	call init
	incl0=.true.
	goto 150

2004	syntax(1)='EXIT [file_name.s] [(COS, HEADER, title.s)]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call writef(file,ierr)
	if(ierr.ne.0) goto 2151
	write(4,'(a3)') 'END'
	goto 3051

2005	syntax(1)='QUIT [SAVE]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	if( nword0(n_len,txt,ib,ie) ) goto 3051
	if(ib  .gt.ie .or. txt(ib:ie).ne.
	1 s_save(1:min(len(s_save),ie-ib+1))) goto 2151
	goto 999

3051	close(2,status='delete',err=30511)
30511	close(48,status='delete',err=30512)
30512	if(verbose.ge.6 ) write(6,1069) 	!000504
1069	format( !&
	1' EdPDB-I6> references:', !&
	1'  Zhang, X. and Matthews, B.(1995).'/ !&
	1'  EdPDB: a multifunctional tool for protein strucure analysis.'/ !&
	1'  J. Appl. Cryst. 28:624-630.'/ !&
	1'  Also see the EdPDB web site at'/ !&
	1'   <http://omega.omrf.ouhsc.edu/zhangc/edpdb/edpdb.html>')
	goto 999

2006	syntax(1)='1) CELL' 
	syntax(2)='2) CELL [a b c alpha beta gamma] [convention#] '
	n_of_syn=2; if (index(txt(:n_len),'???').ge.1) goto 2152

	call dfcell(*2151)
	goto 150

2007	syntax(1)='1) SYMMETRY' 
	syntax(2)='2) SYMMETRY symmetry_operator' 
	syntax(3)='3) SYMMETRY PUNCH file_name [O] '
	syntax(4)='4) SYMMETRY INITIALIZE '
	n_of_syn=4; if (index(txt(:n_len),'???').ge.1) goto 2152

	call dfsymm(*2151)
	goto 150

2008	syntax(1)='APPEND [comment.s] '
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call append(ierr)
	if(ierr.ne.0) goto 2151
	goto 150

2009	syntax(1)='CLOSE [io.i]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	i1=4
	call read_ai(1, i1, *2151,*3091)
3091	if(i1.eq.4) then
	  pdb_out='?'
	  write(4,'(a3)') 'END'
	  endif
	close (i1)
	goto 150

2010	syntax(1)='1) ZONE' 
	syntax(2)='2) ZONE res_id1.s [res_id2.s ...] '
	n_of_syn=2; if (index(txt(:n_len),'???').ge.1) goto 2152

	call zone(*2151)
	goto 150

2011	syntax(1)='MAIN'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call main_chain(i-10,*2151,*20)
	goto 150

2012	syntax(1)='SIDE'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call main_chain(i-10,*2151,*20)
	goto 150

2013	syntax(1)='CA'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call main_chain(i-10,*2151,*20)
	goto 150

2014	syntax(1)='1) INCLUDE'
	syntax(2)='2) INCLUDE another_selection_command.s'
	n_of_syn=2; if (index(txt(:n_len),'???').ge.1) goto 2152

	if(.not.nword(n_len,txt,ib,ie)) then
	  ie=ib-1
	  incl=.true.
	  incl1=.true.
	  backspace (2)
	  goto 503
	  endif
	incl0=.true.
	goto 150

2015	syntax(1)='1) EXCLUDE'
	syntax(2)='2) EXCLUDE another_selection_command.s'
	n_of_syn=2; if (index(txt(:n_len),'???').ge.1) goto 2152

	if(.not.nword(n_len,txt,ib,ie)) then
	  ie=ib-1
	  incl=.false.
	  incl1=.false.
	  backspace (2)
	  goto 503
	  endif
	incl0=.false.
	goto 150

2016	syntax(1)='1) LIST [ n1 [n2] ]'
	syntax(2)='2) LIST ZONE' 
	syntax(3)='3) LIST VECTOR [dmin, dmax, max_center_to_center_dist]'
	n_of_syn=3; if (index(txt(:n_len),'???').ge.1) goto 2152

	call list(*2151)
	goto 150

2017	syntax(1)='1) HELP ; or ???'
	syntax(2)='2) command.s ???'
	n_of_syn=2; if (index(txt(:n_len),'???').ge.1) goto 2152

3017	write(48,1023) version, label
1023	format(' EdPDB> available commands (in version ',a,'):'/(4x,4a16))
	call typelist(48)
	goto 150

2018	syntax(1)=
	1'WRITE file_name.s [(COS, HEADER, BLANK, title.s )] [format.s]' 
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call writef('tmp_.pdb',ierr)
	if(ierr.ne.0) goto 2151
	goto 150

2019	syntax(1)='syntax:'
	syntax(2)='analyze [ANGLE]'
	n_of_syn=2; if (index(txt(:n_len),'???').ge.1) goto 2152

	if(nword(n_len,txt,ib,ie)) then
	  call analy(1)
	else if(txt(ib:ie).eq.'angle') then
          call analy(2)
	  endif
	goto 150

2020	syntax(1)='1) SET B' 
	syntax(2)='2) SET B [b.r]'
	n_of_syn=2; if (index(txt(:n_len),'???').ge.1) goto 2152

	call setb(*2151)
	goto 150

2021	syntax(1)='1) SET WEIGHT' 
	syntax(2)=
	1'2) SET WEIGHT (wv.r, X, Y, Z, B, AVW, CA, file_name.s) '//
     1'[(+W, -W, *W, /W)]'
	n_of_syn=2; if (index(txt(:n_len),'???').ge.1) goto 2152

	call setw(*2151)
	goto 150

2022	syntax(1)='1) B (<, >, <>) cutoff(s) [ANGLE]' 
	syntax(2)='2) B MAX '
	n_of_syn=2; if (index(txt(:n_len),'???').ge.1) goto 2152

	if(nword(n_len,txt,ib,ie)) goto 2151
	if(txt(ib:ie).eq.'max'.and.i.eq.22) then
	  call max_b(*2151)
	  goto 150
	  endif

	call bwxyz(i,txt(ib:ie),*2151)
	goto 150

2023	syntax(1)='W (<, >, <>) cutoff(s) [ANGLE]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	if(nword(n_len,txt,ib,ie)) goto 2151
	call bwxyz(i,txt(ib:ie),*2151)
	goto 150

2024	syntax(1)='X (<, >, <> and ><) cutoff(s) [ANGLE] '
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	if(nword(n_len,txt,ib,ie)) goto 2151
	call bwxyz(i,txt(ib:ie),*2151)
	goto 150

2025	syntax(1)='Y (<, >, <> and ><) cutoff(s) [ANGLE] '
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	if(nword(n_len,txt,ib,ie)) goto 2151
	call bwxyz(i,txt(ib:ie),*2151)
	goto 150

2026	syntax(1)='Z (<, >, <> and ><) cutoff(s) [ANGLE] '
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	if(nword(n_len,txt,ib,ie)) goto 2151
	call bwxyz(i,txt(ib:ie),*2151)
	goto 150

2027	syntax(1)= !&
	1'DISTANCE group_id.s dmin.r dmax.r [skip#.i max_output#.i] '//
     1'[(LOAD, COPY, MARK)]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call dist(*2151)
	goto 150

2028	syntax(1)='1) NAYB radius [res_id [atom_name]]'
	syntax(2)='2) NAYB radius CENTER x, y, z '
	n_of_syn=2; if (index(txt(:n_len),'???').ge.1) goto 2152

	call nayb (*2151)
	goto 150

2029	syntax(1)='RESET'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call readf(ierr)
	if( ierr.ne.0) then
	  write(6,1001) '%EdPDB-F- error during read the pdb file(s)'
	  write(6,*)
	  call exit( 4)
	  endif
	ncm3=0
	goto 150

2030	syntax(1)='ABC [(A, B, C) [(X, Y, Z)]] '
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call abc(*2151)
	goto 150

2031	syntax(1)='ABCD [(A, B, C, D}, [(X, Y, Z)]] '
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call abcd(*2151)
	goto 150

2032	syntax(1)='NAYBR radius.r res_id.s'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call naybr(*2151)
	goto 150


2033	syntax(1)= !&
	1'RTN main_option (parameters) [(SAVE, MULT, INVE) [file_name(s)]] '
	syntax(2)='RTN ABCD res_a.s [atom_a.s] res_b.s [atom_b.s] '
	syntax(3)='   res_c.s [atom_c.s] res_d.s [atom_d.s] tor_ang.r'
	syntax(4)='RTN AXIS vector_id.s rotation_angle.r [translation.r]'
	syntax(5)='RTN CENTER '
	syntax(6)='RTN DEORTH grid_a.r grid_b.r grid_c.r'
	syntax(7)='RTN EZXZ e1.r e2.r e3.r [t1.r t2.r t3.r]'
	syntax(8)='RTN EZYZ e1.r e2.r e3.r [t1.r t2.r t3.r]'
	syntax(9)='RTN FILE file_name.s'
	syntax(10)= !&
	1'RTN MATCH res_id1.s [atom_1.s] res_id2.s [atom_2.s] '//
     1'phi.r, omega.r, kappa.r'
	syntax(11)= !&
	1'RTN MATRIX r11.r, r12.r, r13.r, ... r33.r, t1.r, t2.r, t3.r' 
	syntax(12)='RTN ORTHOG grid_a.r grid_b.r grid_c.r'
	syntax(13)= !&
	1'RTN OVERLAY res_id1.s [atom_11.s atom_12.s atom_14.s] '//
     1'[reg_11.i reg_12.i reg_13.i]'
	syntax(15)= !&
	1'   [res_id2.s [atom_21.s atom_22.s atom_23.s] '
     1//'[reg_21.i reg_22.i reg_23.i]]'
	syntax(16)='RTN POLAR phi.r omega.r kappa.r [t1.r t2.r t3.r]'
	syntax(17)='RTN SYMMETRY [symm_#.i] [fx.r fy.r fz.r]'
	syntax(18)='RTN V_ALIGN vector_id1.s  vector_id2.s'
	syntax(19)='RTN REMARK remark_id.s  remark_#.i'
	n_of_syn=19; if (index(txt(:n_len),'???').ge.1) goto 2152

	jjou=jou
	call edp_rtn(*2151)
	goto 150

2034	syntax(1)='1) DEFINE MAIN'
	syntax(2)='2) DEFINE MAIN atom_1.s [atom_2.s ...]'
	n_of_syn=2; if (index(txt(:n_len),'???').ge.1) goto 2152

	call dfmain(*2151)
	goto 150

2035	syntax(1)='1) DEFINE CA'
	syntax(2)='2) DEFINE CA atom_name.s '
	n_of_syn=2; if (index(txt(:n_len),'???').ge.1) goto 2152

	call dfca(*2151)
	goto 150

2036	syntax(1)='1) DEFINE AB' 
	syntax(2)= !&
	1'2) DEFINE AB atom_a.s atom_b.s [reg_a.i reg_b.i] '//
     1'[status_a.L status_b.L] [Dmin.r Dmax.r]'
	n_of_syn=2; if (index(txt(:n_len),'???').ge.1) goto 2152

	call dfab(*2151)
	goto 150

2037	syntax(1)='1) DEFINE ABC' 
	syntax(2)='2) DEFINE ABC atom_a.s atom_b.s atom_c.s '
	syntax(3)= !&
	1'   [reg_a.i reg_b.i reg_c.i] [status_a.L status_b.L status_c.L]'
     1//' [Dmin.r Dmax.r]'
	n_of_syn=3; if (index(txt(:n_len),'???').ge.1) goto 2152

	call dfabc(*2151)
	goto 150

2038	syntax(1)='1) DEFINE ABCD' 
	syntax(2)='2) DEFINE ABCD atom_a.s atom_b.s atom_c.s atom_d.s '
	syntax(3)= !&
	1'   [reg_a.i ... reg_d.i] [status_a.L ... status_d.L] '//
     1'[Dmin.r Dmax.r]'
	n_of_syn=3; if (index(txt(:n_len),'???').ge.1) goto 2152

	call dfabcd(*2151)
	goto 150

2039	syntax(1)='AB [(A, B) [(X, Y, Z)]]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call ab(*2151)
	goto 150

2040	syntax(1)='SHDF (obsolete) '
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	show=' '
	if(.not.nword(n_len,txt,ib,ie)) show=txt(ib:ie) 
	jerr=show.eq.'*'
	if(show.eq.' ' .or. show.eq. wildcard ) then
		if(incl) then
			write(6,*) 'incl> the selection switch is on'
			else
			write(6,*) 'excl> the deselection switch is on'
		endif
	endif
	if(jerr .or. show.eq.'main') 	call shmain
	if(jerr .or. show.eq.'ca')	call shca
	if(jerr .or. show.eq.'ab') 	call shab
	if(jerr .or. show.eq.'abc') 	call shabc
	if(jerr .or. show.eq.'abcd') 	call shabcd
	if(jerr .or. show.eq.'brg') 	call shbrg
	if(jerr .or. show.eq.'res') 	call shdfres
	if(jerr .or. show.eq.'newx') 	call shnewx
	if(jerr .or. show.eq.'cell') 	call shcell
	if(jerr .or. show.eq.'symm') 	call shsymm
	goto 150

2041	syntax(1)='DEFINE BRIDGE atom_w.s atom_x.s atom_y.s atom_z.s '
	syntax(2)= !&
	1'  [Rw.s reg_w.i Rz.s reg_z.i] [status_w.L status_x.L '//
     1'status_y.L status_z.L]'
	syntax(3)= !&
	1'  [Dmin.r Dmax.r Amin.r Amax.r Tmin.r Tmax.r] '//
     1'[(WXYZ, ZWXY)] [skip.i]'
	n_of_syn=3; if (index(txt(:n_len),'???').ge.1) goto 2152

	call dfbrg(*2151)
	goto 150

2042	syntax(1)='BRIDGE [(W, X, Y, Z)]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

 	call brg(*2151)
	goto 150

2043	syntax(1)='clear '
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152
	call lib_s_erase_page(1,1)
	goto 150

2044	syntax(1)='PAUSE'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	if(.not.inter.or.jou.lt.10) goto 150
	call smg_read(txt
	1,'pause> type <cr> to continue or q to quit ... '
	1,46, n_len, *206)
	call read_txt(0,*206)
	ie=0
	if(.not.nword(n_len,txt,ib,ie)) then
	  if(label(5)(1:min(len(label(5)),ie-ib+1)).eq.txt(ib:ie)) 
	1  goto 206
	  endif
	goto 150

206	do i=jou,11,-1
	  close (i)
	  enddo
	jou=10
	goto 150

2045	syntax(1)='1) SWAP'
	syntax(2)='2) SWAP group_id '
	n_of_syn=2; if (index(txt(:n_len),'???').ge.1) goto 2152

	call swap(*2151)
	goto 150

2046	syntax(1)='1) MORE [i0.i [i1.i]]' 
	syntax(2)='2) MORE CHAIN '
	n_of_syn=2; if (index(txt(:n_len),'???').ge.1) goto 2152

	call more(*2151)
	goto 150

2047	syntax(1)= !&
	1'CLIQUE group_id.s min_clique.i rms_cutoff.r eps.r '//
	1'max_#_cliques.i [mch_char.i]' 
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call mcs_atm_n(*2151)
	goto 150

2048	syntax(1)='1) DEFINE RESIDUE' 
	syntax(2)='2) DEFINE RESIDUE res_type.s [:id.s] [(atom_names.s)]'
	n_of_syn=2; if (index(txt(:n_len),'???').ge.1) goto 2152

	call dfres(*2151)
	goto 150

2049	syntax(1)='SWITCHWB '
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call switchwb
	goto 150

2050	syntax(1)='1) GROUP'
	syntax(2)='2) GROUP group_id '
	n_of_syn=2; if (index(txt(:n_len),'???').ge.1) goto 2152

	call group(*2151)
	goto 150

2051	syntax(1)='CODE (obsolete)'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call code(*2151)
	goto 150

2052	syntax(1)='LOAD group_id.s'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call load(*2151)
	goto 150

2053	syntax(1)='SET CHAIN [chain_name.s]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call setcm(*2151)
	goto 150

2054	syntax(1)='1) REWIND' 
	syntax(2)='2) REWIND (EDP, SCR, PDB)'
	n_of_syn=2; if (index(txt(:n_len),'???').ge.1) goto 2152

	if(nword(n_len, txt, ib, ie)) then
	 if(jou.le.10) then
	  errmsg= ' errmsg: there is no .edp file being read.'
	  goto 2151 
	 else
	  rewind (jou)
	  n_err=n_err+1
	  if(n_err.gt.max_err) then
	    close (jou)
	    jou=jou-1
	    errmsg=' errmsg: too many error/rewind in one file'
	    return 2
	    endif
	 endif
	else if(txt(ib:ie).eq.'out' .or. txt(ib:ie).eq.'scr') then
	  rewind (48)
	  call typelist(-48)
	else if(txt(ib:ie).eq.'edp') then
	  rewind (2)
	else if(txt(ib:ie).eq.'pdb') then
	  rewind (4)
	else 
	goto 2151
	endif
	goto 150

2055	syntax(1)='1) SETENV INTERACTIVE' 
	syntax(2)='2) SETENV INTERACTIVE ON' 
	syntax(3)='3) SETENV INTERACTIVE OFF'
	n_of_syn=3; if (index(txt(:n_len),'???').ge.1) goto 2152

	if(nword(n_len,txt,ib,ie)) then
	  inter=inter0
	else if(txt(ib:ie).eq.'no' .or. txt(ib:ie).eq.'off') then
	  inter=.false.
	else if(txt(ib:ie).eq.'yes' .or. txt(ib:ie).eq.'on') then
	  inter=inter0
	  iadd=1
	  call typelist(48)
	else
	  goto 2151
	endif
	if(inter) then
!		blink =char(7)//char(27)//'[5m'			!vms/vax 
!		normal=char(27)//'[0m'
		blink =char(7)
		normal=' '
	else
		blink =' '
		normal=' '
	endif
	goto 150

2056	syntax(1)='SET ATOM atom_name.s '
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call seta(*2151)
	goto 150

2057	syntax(1)='AVB (X, Y, Z, W)'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call avb(*2151,0)
	goto 150

2058	syntax(1)='BLANK'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call blank()
	goto 150

2059	syntax(1)='NEWXYZ [(A, B, C)]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call newx(*2151)
	goto 150

2060	syntax(1)='1) DEFINE NEWXYZ'
	syntax(2)= !&
	1'2) DEFINE NEWXYZ atom_a.s atom_b.s atom_c.s [reg_a.i '//
     1'reg_b.i reg_c.i]'
	syntax(3)= !&
	1'   [status_a.L status_b.L status_c.L] [distance.r] [angle.r] '//
     1'[torsion_angle.r]'
	n_of_syn=3; if (index(txt(:n_len),'???').ge.1) goto 2152

	call dfnewx(*2151)
	goto 150

2061	syntax(1)='SET RESIDUE residue_type.s'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call setr(*2151)
	goto 150

2062	syntax(1)='1) SET ID' 
	syntax(2)='2) SET ID [new_res_#.i [incr_#.i]] '
	n_of_syn=2; if (index(txt(:n_len),'???').ge.1) goto 2152

	call seti(*2151)
	goto 150

2063	syntax(1)='SET TEXT text_string.s column_1.i column_2.i '
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call sett(*2151)
	goto 150

2064	syntax(1)='SEQUENCE [file_name.s] [format.s] [(C,R)]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call sequ(file,*2151)
	goto 150

2065	syntax(1)='SUMW (X, Y, Z, B)'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call avb(*2151,1)
	goto 150

2066	syntax(1)='1) FILE '
	syntax(2)='2) FILE OUT file_name.s'
	syntax(3)='3) FILE LOG file_name.s'
	n_of_syn=3; if (index(txt(:n_len),'???').ge.1) goto 2152


	if(.not.nword(n_len,txt,ib,ie)) then
	  if(txt(ib:ie).eq.'out' .or. txt(ib:ie).eq.'scr') then
	    scr_out='?'
	    if(.not.open_file1(48, scr_out,'new','.out')) goto 998
	    call typelist(-48)
	  else if(txt(ib:ie).eq.'log') then
	    log_out='?'
	    if(.not.open_file1(6, log_out,'unknown','.log')) goto 998
	  else 
	    goto 2151
	    endif
	  goto 150

998	  errmsg=' errmsg: error opening a new file'
	  goto 2151
	  endif

	if(ncm1.le.0 .or. cm1.eq.'_') then
	  cm_1='?'
	else
	  cm_1=cm1
	endif
	if(ncm2.le.0 .or. cm2.eq.'_') then
	  cm_2='?'
	else
	  cm_2=cm2
	endif
	if(ncm3.le.0 .or. cm3.eq.'_') then
	  cm_3='?'
	else
	  cm_3=cm3
	endif
	write(6,10661) version	!, 	!&
	write(6,10662) pdb_in (:ltrim(pdb_in )) 
	1 //' ['//cm_1(:ltrim(cm_1))//']'
	if(pdb_in2(1:1).ne.'?'.or. verbose.ge.3) 
	1 write(6,10662) pdb_in2(:ltrim(pdb_in2))
	1 //' ['//cm_2(:ltrim(cm_2))//']'
	if(pdb_in3(1:1).ne.'?'.or. verbose.ge.3) 
	1 write(6,10662) pdb_in3(:ltrim(pdb_in3))
	1 //' ['//cm_3(:ltrim(cm_3))//']'
	if(pdb_out(1:1).ne.'?'.or. verbose.ge.3)  
	1 write(6,10663) pdb_out(:ltrim(pdb_out))
	if(rtn_in(1:1).ne.'?'.or. verbose.ge.3)  
	1 write(6,10664) rtn_in (:ltrim(rtn_in))
	if(rtn_out(1:1).ne.'?'.or. verbose.ge.3) 
	1 write(6,10665) rtn_out(:ltrim(rtn_out))
	if(edp_in(1:1).ne.'?'.or. verbose.ge.3) 
	1 write(6,10666) edp_in (:ltrim(edp_in))
	write(6,10667) edp_out(:ltrim(edp_out))
	write(6,10668) scr_out(:ltrim(scr_out))
	if(log_out(1:1).ne.'?'.or. verbose.ge.3) 
	1 write(6,10669) log_out(:ltrim(log_out))
	if(acc_in(1:1).ne.'?'.or. verbose.ge.3) 
	1 write(6,10670) acc_in (:ltrim(acc_in))
	if(std_in(1:1).ne.'?'.or. verbose.ge.3) 
	1 write(6,10671) std_in (:ltrim(std_in))
	if(seq_out(1:1).ne.'?'.or. verbose.ge.3) 
	1 write(6,10672) seq_out(:ltrim(seq_out))
	if(update_in(1:1).ne.'?'.or. verbose.ge.3) 
	1 write(6,10673) update_in(:ltrim(update_in)) 
10661	format( !&
	1' program:'         	,t16,'EdPDB (',a,')')
10662	format(' pdb_in (.pdb)'	,t16,a	)
10663	format(' pdb_out(.pdb)'	,t16,a	)
10664	format(' rtn_in'	,t16,a	)
10665	format(' rtn_out'	,t16,a	)
10666	format(' edp_in (.edp)'	,t16,a	)
10667	format(' edp_out(.edp)'	,t16,a	)
10668	format(' scr_out(.out)'	,t16,a	)
10669	format(' log_out(.log)'	,t16,a	)
10670	format(' acc_in'	,t16,a	)
10671	format(' std_in'	,t16,a	)
10672	format(' seqence(.seq)'	,t16,a	)
10673	format(' update_in'	,t16,a	) 
	goto 150

2067	syntax(1)='SET field_name parameters'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call set_x(*2151)
	goto 150

2068	syntax(1)='READ file_name.s [mark.s] [INITIALIZE]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call readf1(*2151)		! read one more PDB file.
	goto 150

2069	syntax(1)='1) SNAYB radius.r [res_id.s [atom_name.s]]' 
	syntax(2)='2) SNAYB radius.r CENTER x.r, y.r, z.r '
	n_of_syn=2; if (index(txt(:n_len),'???').ge.1) goto 2152

	call snayb(*2151)
	goto 150

2070	syntax(1)='SNAYBR radius.r res_id.s'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call snaybr(*2151)
	goto 150

2071	syntax(1)='SYSTEM [system_command.s]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	if( nword(n_len,txt,ib,ie)  )  then
	  call lib_s_spawn0
	else if(txt(ib:ie).eq.s_wait(1:min(len(s_wait),ie-ib+1))) then
	  if( nword(n_len,txt,ib,ie)  )  then
	    call lib_s_spawn0
	  else
	    call lib_s_spawn1 (txt(ib:n_len))
	    endif
	else  
	  call lib_s_spawn2 (txt(ib:n_len))
	endif
	goto 150

2072	syntax(1)='AXIS file_name.s [vector_id.s, [axis_id.s]] '
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call axisdist(*2151)
	goto 150

2073	syntax(1)='JIGGLE (X, Y, Z, W, B) limit.r [shift.r] '
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call jiggle(*2151)
	goto 150

2074	syntax(1)='SET ENTRY'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call sete(*2151)
	goto 150	

2075	syntax(1)='CORRELATION grp_id.s (X,Y,Z,W,B) (X,Y,Z,W,B) '//
	1'[(X,Y,Z,W,B)] '
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call correlatn(*2151)
	goto 150

2076	syntax(1)= !&
	1'ACCESS [grp_id.s] [ISOLATED,BURIED] [r_probe.r] [zstep.r]'//
	1' [file_name.s]' 
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call acc(*2151)
	goto 150

2077	syntax(1)='MMIG (obsolete; SDISTANCE is used)'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call sdistance(*2151)
	goto 150

2078	syntax(1)='DIFF group_id.s [(RMS,MORPH [lamda.r],'//
	1' SCALE [scale1.r, scale2.r])]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	igr=match_l( max_gr, cgroup)
	if( igr .le. 0) then
	  errmsg= ' errmsg: a group name is needed here'
	  goto 2151
	endif
	if(nword0(n_len,txt,ib,ie)) then
	  call diff_g_o(igr,0,1.0,1.0)
	else if(ib .le. ie .and. txt(ib:ie).eq.s_rms(1:min(len(s_rms) 
	1,ie-ib+1))) then
	  call diff_g_o(igr,1,1.0,1.0)
	else if(ib .le. ie .and. txt(ib:ie).eq.s_scale(1:min(len(s_scale)
	1,ie-ib+1))) then
	  scale1=1.
	  scale2=1.
	  call read_ar(1, scale1, *2151, *2781)
	  call read_ar(1, scale2, *2151, *2781)
2781	  call diff_g_o(igr,0,scale1,scale2)
	else if(ib .le. ie .and. txt(ib:ie).eq.s_morph(1:min(len(s_morph)
	1,ie-ib+1))) then
	  scale1=1.
	  scale2=0.
	  call read_ar(1, scale1, *2151, *2782)
2782	  scale2=scale1-1.0
	  call diff_g_o(igr,0,scale1,scale2)
	else
	  goto 2151
	endif
	goto 150

2079	syntax(1)='OVERLAY group_id.s [file_name.s] [WEIGHT]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call edp_rms(*2151)
	goto 150

2080	syntax(1)='RMSW (X, Y, Z, B)'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call avb(*2151,2)
	goto 150

2081	syntax(1)='SETENV PROMPT prompt_character_string.s' 
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	if( nword(n_len,txt,ib,ie)  )  goto 2151
	ic= ichar(delimiter)
	if( ic .eq. ichar( txt(ib:ib)) ) then	! using '...' to input prompt
	  ib=ib+1
	  do j=ib, n_len
	    if( ichar(txt(j:j)) .eq. ic) then
	      txt(j:j)=' '	!000303, it may screw the history file.
	      ie= j-1
	      goto 2811
	      endif
	    end do
	  ie= n_len
	  endif
2811	n_prompt= min( 9, ie-ib+1)
	if(n_prompt.le.0) then
	  n_prompt=1
	  prompt='*'
	else
	  prompt(:n_prompt)= txt(ib:ie)
	  endif
	goto 150

2082	syntax(1)=
	1'load grp_id2.s | MATCH grp_id1.s sub_grp_id1.s [fmtstmt.s]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call edp_match(*2151)
	goto 150

2083	syntax(1)='SORT '
	syntax(2)='SORT [(-B, B, -W, W)]'
	syntax(3)='SORT [(DFRES, SWAP, LOAD)]'
 	n_of_syn=3; if (index(txt(:n_len),'???').ge.1) goto 2152

	call order(*2151)
	goto 150

2084	syntax(1)='RATIO group_id.s [scale.r] [def_value.r]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	igr=match_l( max_gr, cgroup)
	if( igr .le. 0) then
	  errmsg= ' errmsg: a group name is needed here'
	  goto 2151
	endif
	
	scale=1.
	dfr=999.
	call read_ar( 1, scale, *2151, *2831)
2831	call read_ar( 1, dfr, *2151, *2832)
2832	call ratio_g_o(igr,scale,dfr)
	goto 150

2085	syntax(1)='PLANAR vector_id.s'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call planar(*2151)
	goto 150

2086	syntax(1)='MATCH3D grp_id.s min_clique_num.i max_rms.r '//
	1'[file_name.s] [NONSEQU]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call match3d(*2151)
	goto 150

2087	syntax(1)='load grp_id2 | ALIGN3D grp_id1 sub_grp_id1 '
	syntax(2)='  [distance_cutoff [penalty_to_break]] '
	n_of_syn=2; if (index(txt(:n_len),'???').ge.1) goto 2152

	call align3d(*2151)
	goto 150

2088	syntax(1)='1) PERMUTE' 
	syntax(2)='2) PERMUTE i0.i i1.i shift.i'
	n_of_syn=2; if (index(txt(:n_len),'???').ge.1) goto 2152

	call permut(*2151,*2016)
	goto 150

2089	syntax(1)='1) PARAMETER [Pn]' 
	syntax(2)='2) PARAMETER Pn = [value] '
	syntax(3)='3) PARAMETER Pn ? [prompt_string] [default_value] '//
	1'[(EXIT, REPORT)]' 
	syntax(4)='4) PARAMETER Pn (+, -) step_size limit [(EXIT, REPORT)] '
	syntax(5)='5) PARAMETER Pn.s group_id.s (ENTRY, ATOM, RESIDUE, '//
	1'CHAIN, ID, X, Y, Z, W, B)'
	syntax(6)='	    [(EXIT, REPORT)]'
	n_of_syn=6; if (index(txt(:n_len),'???').ge.1) goto 2152

	call param(*2151,*2891)
	goto 150
2891	if(jou.le.10) then
	  write(6,*) '%EdPDB-F- wrong parameters or top level.'
	  call exit( 4) 
          endif

3202	close (jou)
	jou=jou-1
	goto 150

2090	syntax(1)='SETENV MAXERR max_err.i [(EXIT, QUIT)]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	if( nword(n_len,txt,ib,ie)  )  then
	  write(6,*) 'EdPDB>  max_number_of_error=',max_err
	  write(6,*) 'EdPDB> current_num_of_error=',n_err
	  goto 150
	  endif
	read( txt(ib:ie),*,err=2151) max_err
!1063	format(i<ie-ib+1>)
	n_err=0
	if( nword(n_len,txt,ib,ie)  )  then
	 err_control=s_exit
	else if(txt(ib:ie).eq.s_exit(1:min(len(s_exit),ie-ib+1))) then
	 err_control=s_exit
	else if(txt(ib:ie).eq.s_quit(1:min(len(s_quit),ie-ib+1))) then
	 err_control=s_quit
	else 
	 goto 2151
	endif
	goto 150

2091	syntax(1)='load grp_id.s | PATTERN pattern_string.s [position.i]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call pattern(*2151)
	goto 150

2092	syntax(1)='SHAPE search_radius.r res_id.s [atom_name.s]'
	syntax(2)= !&
	1'  max_RT.r [probe_radius.r] [file_name.s] [random_seed.i] [box(6).r]'
	n_of_syn=2; if (index(txt(:n_len),'???').ge.1) goto 2152

	call edp_shape(*2151)
	goto 150

2093	syntax(1)='VOLUME [probe_radius.r] [zstep.r] [file_name.s]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call volume(*2151)
	goto 150

2094	syntax(1)='1) ALIAS' 
	syntax(2)='2) ALIAS key.s'
	syntax(3)='3) ALIAS key.s := value.s ;or ALIAS key.s value.s ' 
	syntax(4)='4) ALIAS key.s :=' 
	n_of_syn=4; if (index(txt(:n_len),'???').ge.1) goto 2152

	call udk1(*2151)
	goto 150

2095	syntax(1)='MOVECENTER [file_name.s] [fx1.r fy1.r fz1.r '//
	1'[fx2.r fy2.r fz2.r]]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call move2o(*2151)
	goto 150

2096	syntax(1)='1) CHAIN' 
	syntax(2)='2) CHAIN [chn_mark_1.s chn_mark_2.s ...]'
	n_of_syn=2; if (index(txt(:n_len),'???').ge.1) goto 2152

	call chain(*2151)
	goto 150

2097	syntax(1)='HARKER [grid_a.r grid_b.r grid_c.r] [symm_#1.i '//
	1'[symm_#2.i]] [CROSS]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call harker(*2151)
	goto 150

2098	syntax(1)='1) EULER TO_EULER' 
	syntax(2)='2) EULER TO_POLAR'
	syntax(3)='3) EULER SYMMETRY symm_#.i'
	syntax(4)='4) EULER MOVE_TO_O res_id.s, atom_name.s'
	syntax(5)='5) EULER ASYMM [e1.r, e2.r, e3.r] '
	n_of_syn=5; if (index(txt(:n_len),'???').ge.1) goto 2152

	call std_euler(*2151)
	goto 150

2099	syntax(1)='1) POLAR TO_POLAR'
	syntax(2)='2) POLAR TO_EULER'
	syntax(3)='3) POLAR MOVE_TO_O res_id.s atom_name.s' 
	syntax(4)='4) POLAR ASYMM [p1.r, p2.r, p3.r]' 
	syntax(5)='5) POLAR SRF_RED [p1.r, p2.r, p3.r]' 
	syntax(6)='6) POLAR UNIQUE delta_angle.r '
	n_of_syn=6; if (index(txt(:n_len),'???').ge.1) goto 2152

	call std_polar(*2151)
	goto 150

2100	syntax(1)='VM [molecular_weight_in_kd.r] [num_symm_op.i]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call read_ar( 1, mw, *2151, *3001)
	if(mw.le.0.0) goto 2151
! assume that if the molecular weight number is smaller than 1000,
! we are talking about kilo-dalton.
	if(mw.le.5000.0) mw=mw*1000.0	 
	goto 3002
3001	call get_mw(mw)
	if(mw.le.0.0) goto 2151
3002	call vm(mw,*2151)
	goto 150

2101	syntax(1)='load grp_id2.s | EXTRACT grp_id1.s sub_grp_id1.s'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call extract(*2151); goto 150

2102	syntax(1)='MOMENTINERTIA [file_name.s] [vector_id1.s] '//
	1'[vector_id2.s] [vector_id3.s]' 
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call rotmomentum(*2151)
	goto 150

2103	syntax(1)='1) UPDATE (XYZ, W, B) file_name.s fortran_format.s' 
	syntax(2)='2) UPDATE T column_1.i column_2.i file_name.s '
	syntax(3)='    	fortran_format.s [jump_after_string]'
	n_of_syn=3; if (index(txt(:n_len),'???').ge.1) goto 2152

	call update(*2151)
	goto 150

2104	syntax(1)='CLOSER grp_id1.s grp_id2.s dmax.r'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call closer(*2151)
	goto 150

2105	syntax(1)= '1)  SETENV -- display the current values'
	syntax(2)= '2)  SETENV -(S,R) variable.s -- save, restore'
	syntax(3)= '3)  SETENV DELIMITER [delimiter.x] -- e.g. # ' 
	syntax(4)= '4)  SETENV ECHO [(0,1,2)[ '
	syntax(5)= '5)  SETENV EDP_DATA [directory_name.s] '
	syntax(6)= '6)  SETENV MAXERR [max_err.i [(exit, quit)]]'
	syntax(7)= '7)  SETENV TOLOWER [(on, off)]'
	syntax(8)= '8)  SETENV VERBOSE [level.i]' 
	syntax(9)= '9)  SETENV WILDCARD [wildcard.s] -- e.g. + '
	syntax(10)='10) SETENV WINDOW_SIZE [#OfLinesPerWindow.i]' 
	syntax(11)='11) SETENV PROMPT [text.s]' 
	n_of_syn=11; if (index(txt(:n_len),'???').ge.1) goto 2152

	if(nword(n_len,txt,ib,ie)) then
	  write(6,*) 	'delimiter  := ',delimiter
	  write(6,*) 	'echo       := ',echo_L
	  write(6,*) 	'edp_data   := ',edp_data(:ltrim(edp_data))
	  if(inter) then
	    write(6,*) 	'interactive:= on'
	  else
	    write(6,*) 	'interactive:= off'
	    endif
	  write(6,*) 	'maxerr     := ',max_err,' ',err_control
	1 ,' (curr#err=',n_err,')'
	  if(tolower) then
	    write(6,*) 	'tolower    := on'
	  else
	    write(6,*) 	'tolower    := off'
	    endif
	  write(6,*) 	'verbose    := ',verbose
	  write(6,*) 	'wildcard   := ',wildcard//' '//wildcard1
	  write(6,*) 	'window_size:= ',window_size
	  write(6,*) 	'prompt     := ',prompt(:n_prompt)
	  goto 150
	  endif

	if(txt(ib:ie).eq.'-s' .or. notSaveYet) then
	  ie_setenv=ie
	  ib_setenv=ib
	  save_all=nword(n_len,txt,ib,ie)
	  if(notSaveYet) save_all=.true.
 
	  if(txt(ib:ie).eq. stterm(1)(:1+ie-ib).or.save_all) 
	1 verbose_save=verbose
	  if(txt(ib:ie).eq. stterm(2)(:1+ie-ib).or.save_all) 
	1 edp_data_save=edp_data
	  if(txt(ib:ie).eq. stterm(3)(:1+ie-ib).or.save_all) 
	1 window_size_save=window_size
	  if(txt(ib:ie).eq. stterm(4)(:1+ie-ib).or.save_all) 
	1 echo_L_save=echo_L
	  if(txt(ib:ie).eq. stterm(5)(:1+ie-ib).or.save_all) 
	1 tolower_save=tolower
	  if(txt(ib:ie).eq. stterm(6)(:1+ie-ib).or.save_all) 
	1 inter_save=inter
	  if(txt(ib:ie).eq. stterm(7)(:1+ie-ib).or.save_all) 
	1 delimiter_save=delimiter
	  if(txt(ib:ie).eq. stterm(8)(:1+ie-ib).or.save_all) then
	    wildcard_save  =wildcard
	    wildcard1_save =wildcard1
	    endif
	  if(txt(ib:ie).eq. stterm(9)(:1+ie-ib).or.save_all) then
	    prompt_save   =prompt
	    n_prompt_save =n_prompt
	    endif
	  if(txt(ib:ie).eq. stterm(10)(:1+ie-ib).or.save_all ) then 
	    max_err_save=max_err
	    n_err_save=n_err
	    err_control_save=err_control
	    endif
	  notSaveYet=.false.
	  if(txt(ib_setenv:ie_setenv).ne.'-s') then 
	    ib=ib_setenv
	    ie=ie_setenv
	  else 
	    if( ib.gt.ie) goto 150
	    endif
	  endif

	if(txt(ib:ie).eq.stterm(1)(:1+ie-ib)) then
	  call read_ai(1,verbose,*2151,*3050)
	  goto 150
3050	  verbose= 1 
	  goto 150
	else if(txt(ib:ie).eq.stterm(2)(:1+ie-ib)) then
	  if(nword(n_len,txt,ib,ie)) then
	    edp_data=edp_data0
	  else
	    edp_data=txt(ib:n_len)
	    endif
	  call set_parameter(3,'edp_data',8,edp_data,
	1 index(edp_data,'/ ')-1)
	  goto 150
	else if(txt(ib:ie).eq.stterm(3)(:1+ie-ib)) then
	  call read_ai(1,window_size,*2151,*4051)
	  goto 150
4051	  window_size=window_size0
	  goto 150
	else if(txt(ib:ie).eq.stterm(4)(:1+ie-ib)) then
	  call read_ai(1,echo_L,*2151,*2151)
	  goto 150
	else if(txt(ib:ie).eq.stterm(5)(:1+ie-ib)) then
	  if(nword(n_len,txt,ib,ie)) then
	    tolower=.true.
	  else if(txt(ib:ie).eq.'on') then
	    tolower=.true.
	  else if(txt(ib:ie).eq.'off') then
	    tolower=.false.
	    endif
	  goto 150
	else if(txt(ib:ie).eq.stterm(6)(:1+ie-ib) )  then! 'interactive'
	  goto 2055
	else if(txt(ib:ie).eq.stterm(7)(:1+ie-ib)) then
	  if(nword(n_len,txt,ib,ie)) then
	    delimiter=char(39)
	  else
	    delimiter=txt(ib:ib)
	    endif
	  goto 150
	else if(txt(ib:ie).eq.stterm(8)(:1+ie-ib)) then
	  if(nword(n_len,txt,ib,ie)) then
	    wildcard  ='*'
	    wildcard1 ='%'
	  else
	    wildcard =txt(ib:ib)
	    if(nword(n_len,txt,ib,ie)) then
	      wildcard1='%'
	    else
	      wildcard1=txt(ib:ib)
	      endif
	    endif
	  goto 150

	else if(txt(ib:ie).eq.stterm(9)(:1+ie-ib) ) then ! 'prompt'
	  goto 2081		!010221

	else if(txt(ib:ie).eq.stterm(10)(:1+ie-ib)) then ! 'maxerr'
	  call read_ai(1,max_err,*2151,*3052)
	  goto 3053
3052	  max_err=1024
3053	  n_err=0
	  if( nword(n_len,txt,ib,ie)  )  then
	    err_control=s_quit
	  else if(txt(ib:ie).eq.s_exit(1:min(len(s_exit),ie-ib+1))) then
	    err_control=s_exit
	  else if(txt(ib:ie).eq.s_quit(1:min(len(s_quit),ie-ib+1))) then
	    err_control=s_quit
	  else 
	    goto 2151
	    endif
	  goto 150

	else if(txt(ib:ie).eq.'-r') then
	  if(notSaveYet) then
	    errmsg=' errmsg: UNDONE: you should save them (using -s) first '
	    goto 2151
	    endif
	  if(nword(n_len,txt,ib,ie)) then
	    save_all=.true.
	  else
	    save_all=.false.
	    endif

	  if(txt(ib:ie).eq.stterm(1)(:1+ie-ib) .or.save_all) 
	1 verbose=verbose_save
	  if(txt(ib:ie).eq.stterm(2)(:1+ie-ib) .or.save_all) 
	1 edp_data=edp_data_save
	  if(txt(ib:ie).eq.stterm(3)(:1+ie-ib) .or.save_all) 
	1 window_size=window_size_save
	  if(txt(ib:ie).eq.stterm(4)(:1+ie-ib) .or.save_all) 
	1 echo_L=echo_L_save
	  if(txt(ib:ie).eq.stterm(5)(:1+ie-ib) .or.save_all) 
	1 tolower=tolower_save
	  if(txt(ib:ie).eq.stterm(6)(:1+ie-ib) .or.save_all) 
	1 inter=inter_save
	  if(txt(ib:ie).eq.stterm(7)(:1+ie-ib) .or.save_all) 
	1 delimiter=delimiter_save
	  if(txt(ib:ie).eq.stterm(8)(:1+ie-ib) .or.save_all) then
	    wildcard  =wildcard_save
	    wildcard1 =wildcard1_save
	    endif
	  if(txt(ib:ie).eq.stterm(9)(:1+ie-ib)  .or.save_all ) then ! 'prompt'
	    prompt=prompt_save
	    n_prompt=n_prompt_save
	    endif
	  if(txt(ib:ie).eq.stterm(10)(:1+ie-ib) .or.save_all ) then ! 'maxerr'
	    max_err=max_err_save
	    n_err=n_err_save
	    err_control=err_control_save
	    endif
	  goto 150
	  endif

	goto 2151

2106	goto 2007	!operator.eq.symmetry

2107	goto 150	!lattice --  not used

2108	syntax(1)='TEXT text_string.s [t1.i [t2.i]] '
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call search_text(*2151)
	goto 150

2109	syntax(1)='RETURN'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152
	
	goto 3202

2110	syntax(1)='GOTO label.s'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	if(jou.le.10) then
	  errmsg=
	1 ' errmsg: a goto statement can not be used at the top level.'
	  goto 2151
	  endif

	if(nword(n_len,txt,ib,ie)) then
	  errmsg=' errmsg: a label is required.'
	  goto 2151
	  endif

	goto_label=txt(ib:ie)//':'
	ngoto_label=ie-ib+2		! 990713
	ib0=ib
	ie0=ie
	scan_times=0
	do while( scan_times .le. 1 )
	  scan_times = scan_times + 1
	  do while(.true.)
	    read(jou,'(a)',end=3101) txt
	    n_len=ltrim(txt)
	    call read_txt(0,*150)		! convert to lower cases
	    if(index(txt,goto_label(:ngoto_label)) .gt. 0) then	
	      ie=0
	      if(.not.nword(n_len,txt,ib,ie)) then
	        if(txt(ib:ie).eq.goto_label(:ngoto_label)) then
!	          write(2,1001)		'!!'//trim(txt)
!	          write(48,1001)	'!!'//trim(txt)
!	          if(echo_L.ge.1) write(6,1001)' '//prompt(:n_prompt)//trim(txt)
	          n_err=n_err+1
              backspace(jou)
	          goto 150
	          endif
	        endif
	      endif
	    enddo
3101	  if(scan_times .le. 1) rewind(jou)
	  enddo
	errmsg=' errmsg: the label statement can not be found.'
	close (jou)
	jou=jou-1
	ib=ib0
	ie=ie0
	goto 2151

2111	syntax(1)=
	1 'IF ( [-(S,R,I)] parameter_name operator value ) command.s'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call true_or_false(edp_if,*2151)
	if(edp_if)  then
	  do j=1,ie
	    txt(j:j)=' '
	    enddo
	  goto 502
	  endif
	goto 150

2112	syntax(1)='1) vector by_atom res_1.s atom_1.s res_2.s atom_2.s' 
	syntax(2)=
	1'2) vector by_xyz vector_id.s p1.r p2.r p3.r r1.r r2.r r3.r '//
	1'[length.r]' 
	syntax(3)='3) vector delete vector_id.s '
	syntax(4)='4) vector list [vector_id.s] '
	syntax(5)='5) vector pv vector_id.s [res_id.s [atom_name.s]]'
	syntax(6)='6) vector vp vector_id.s [res_id.s [atom_name.s]] '//
	1'[length.r]' 
	n_of_syn=6; if (index(txt(:n_len),'???').ge.1) goto 2152

	call list_vector(*2151)
	goto 150

2113	syntax(1)='MKFILE filename.s [eof.s]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	if(.not.open_file1(21, tmp_edp,'unknown','.edp')) then
	  errmsg=' errmsg: error in opening the file'
	  goto 2151
	  endif
	making_edp=.true.
	if(.not.nword(n_len,txt,ib,ie)) then
	  neof=ie-ib+1
	  if(neof.ge.1) then
	    eof=txt(ib:ie)
	    goto 150
	    endif
	  endif
	eof='eof'
	neof=3
	goto 150

2114	syntax(1)='MW [molecular_weight_in_Da, delta]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call read_ar( 1, mw, *2151, *21141)
	if(mw.le.0.0) goto 2151
	call read_ar( 1, delta, *2151, *2151)
	call match_mw(mw,delta)
	goto 150
21141	call get_mw(mw)
	goto 150

2115	syntax(1)='load grp_id2.s | MATCH1D grp_id1.s sub_grp_id1.s '
	syntax(2)=
	1'  [score_threshold.r [list_threshold(1).r [list_threshold(2).r, '
	syntax(3)=
	1'  [extension_penalty.r, [num_of_random_trial.i '//
     1'[random_seed.i]]]]]] '
	n_of_syn=3; if (index(txt(:n_len),'???').ge.1) goto 2152
	
	if(threading('match1d').ge.1 ) goto 2151 
	goto 150

2116	syntax(1)='SEQ2PDB [file_name.s] [format.s]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call seq2pdb(*2151)
	goto 150

2117	syntax(1)='COPY group_name.s [(X,Y,Z,W,B), T [col1.i [ col2.i]]]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call copy(*2151)
	goto 150

2118	syntax(1)='SDISTANCE group_id.s dmin.r dmax.r '//
	1'[(LOAD, MOVE, PUNCH_ALL)]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call sdistance(*2151)
	goto 150

2119	goto 2069	!snayb.eq.mmi

2120	goto 2070	!snaybr.eq.mmir

2121	syntax(1)='1) DEFINE -- display the current values'
	syntax(2)='2) DEFINE CA -- pseudo ca atom'
	syntax(3)='3) DEFINE MAIN -- pseudo main chain atoms'
	syntax(4)='4) DEFINE AB -- bond type'
	syntax(5)='5) DEFINE ABC -- bond angle type'
	syntax(6)='6) DEFINE ABCD -- torsional angle type'
	syntax(7)='7) DEFINE BRIDGE -- see bridge'
	syntax(8)='8) DEFINE NEWXYZ -- see newxyz'
	syntax(9)=
	1'9) DEFINE RESIDUE -- for ordering atoms using sort/dfres'
	n_of_syn=9; if (index(txt(:n_len),'???').ge.1) goto 2152

	if(nword(n_len,txt,ib,ie)) then
	  call shca
	  call shmain
	  call shab
	  call shabc
	  call shabcd
	  call shbrg
	  call shdfres
	  call shnewx
 	  goto 150
 	  endif

	if(txt(ib:ie).eq.dfterm(1)(:1+ie-ib)) then
	  call dfca(*2151)
	else 	if(txt(ib:ie).eq.dfterm(2)(:1+ie-ib)) then
	  call dfmain(*2151)
	else 	if(txt(ib:ie).eq.dfterm(3)(:1+ie-ib)) then
	  call dfab(*2151)
	else 	if(txt(ib:ie).eq.dfterm(4)(:1+ie-ib)) then
	  call dfabc(*2151)
	else 	if(txt(ib:ie).eq.dfterm(5)(:1+ie-ib)) then
	  call dfabcd(*2151)
	else 	if(txt(ib:ie).eq.dfterm(6)(:1+ie-ib)) then
	  call dfbrg(*2151)
	else 	if(txt(ib:ie).eq.dfterm(7)(:1+ie-ib)) then
	  call dfnewx(*2151)
	else 	if(txt(ib:ie).eq.dfterm(8)(:1+ie-ib)) then
	  call dfres(*2151)
	else 	
	  goto 2151
	  endif
	goto 150

2122	syntax(1)='load group2.s | THREAD  group1.s sub_group1.s '
	syntax(2)='  [score_threshold.r [scale1.r, [scale2.r,'
	syntax(3)='  [extension_penalty.r, [num_of_random_trial.i, '//
	1'[seed.i]]]]]]' 
	n_of_syn=3; if (index(txt(:n_len),'???').ge.1) goto 2152

	if(threading('thread').ge.1) goto 2151
	goto 150

2123	syntax(1)='W2T [scale.r] [score_file.s]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	if(w2t().ge.1) goto 2151
	goto 150
	
2124	syntax(1)='syntax:' 
	syntax(2)= 'doit [gr2.s S/R s2.i d2min.r d2max.r '
	syntax(3)= ' [gr3.s S/R s3.i d3min.r d3max.r A/D a3min.r a3max.r'
	syntax(4)= ' [gr4.s S/R s4.i d4min.r d4max.r A/D a4min.r a4max.r' !&
	1         //' T/H t4min.r t4max.r'
	syntax(5)= ' [gr5.s S/R s5.i d5min.r d5max.r A/D a5min.r a5max.r' !&
	1         //' T/H t5min.r t5max.r]]]]'
	n_of_syn=5; if (index(txt(:n_len),'???').ge.1) goto 2152

	call doit(*2151)
	goto 150

2125	syntax(1)= 'pdist gr2.s d2min.r d2max.r '
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call dist_p(*2151)
	goto 150

2126	syntax(1)='poverlay group_id.s [file_name.s] [WEIGHT]'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call edp_rms_p(*2151)
	goto 150

2127	syntax(1)= !&
	1'pclique group_id.s min_clique.i rms_cutoff.r eps.r '//
	1'max_#_cliques.i [mch_char.i]' 
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call mcs_atm_p(*2151)
	goto 150
	
2128	syntax(1)='syntax:'
	syntax(2)='pdoit'
	n_of_syn=2; if (index(txt(:n_len),'???').ge.1) goto 2152

	call analy(2)
	goto 150

2129	syntax(1)= 'psdist gr2.s d2min.r d2max.r '
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call sdist_p(*2151)
	goto 150

2130	syntax(1)= !&
	1'sclique group_id.s min_clique.i eps.r '//
	1'max_#_cliques.i [mch_char.i]' 
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call mcs_atm_sp(*2151)
	goto 150
	
2131	syntax(1)='aniso u11.r u22.r u33.r u23.r u31.r u12.r '
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call aniso(*2151)
	goto 150

2132	syntax(1)='remark text.s'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152


	if(pdb_out.eq.'?') then
	  write(6,*)
	1'remark-W> UNDONE: open an output pdb file first'
	  return
	  endif
	write(4,1132) txt(ie+1:n_len)
	goto 150
1132	format('REMARK',a)	

2133	syntax(1)='dock z1.r z2.r dz.r shell.r'
	n_of_syn=1; if (index(txt(:n_len),'???').ge.1) goto 2152

	call dock(*2151)
	goto 150

2134	syntax(1)='touch box_margin.r (def=9.0; min=1.4) probe_margin.r (def=1.0; min=0.2)' 
	syntax(2)='      [file_name.s] '
	n_of_syn=2; if (index(txt(:n_len),'???').ge.1) goto 2152

	call edp_touch(*2151)
	goto 150

1001	format(a)

999	call my_stop
	end !subroutine get_command

!copyright by X. Cai Zhang
