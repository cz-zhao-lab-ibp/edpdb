	subroutine vm(mw,*)
chk	=============
c	this subroutine calculates the Matthews' coefficient.
	include 'edp_dim.inc'
!	use edp_dim
	include 'edp_dat.inc'
!	use edp_dat


	character*(max_num_chars) errmsg
	character*(max_num_chars) syntax(32)
	common /error_msg/ errmsg, max_err, n_err, syntax, n_of_syn
	real mw, cell(12)
	common /cmm_symm/ num_symm, junk(24)
	equivalence 
	1  ( a,cell(1)), ( b,cell(2)),( c,cell(3))
	1 ,(al,cell(4)), (be,cell(5)),(ga,cell(6))

!	character*(max_num_chars) txt
!	common /cmm_txt/n_len,txt,ib,ie
!	logical nword, nword0

	real junk1(3,4)	!,junk2(3,4)
	character*32 junk_txt
	external no_symm
	logical  no_symm

	if(verbose.ge.6 ) write(6,1069) 	!000504
1069	format(' vm-I6> reference: Matthews, B.W. (1968)'
	1/'  Solvent content of protein crystals.J. Mol. Biol., 33(2):491-7.'
	1/'  And assume water density = 1 gm/cm^3.')

	isymm=-3   ! to get cell dimension isymm<=-3   5-sep-1991
	call get_trn( isymm,cell,junk1,junk_txt)
	if(isymm.le.-3) return 1

	num_op=0
	call read_ai(1,num_op, *900, *101)
101	if(num_op.le.0) then 
	 if(num_symm.le.0) then
	  if(no_symm()) then	
	    write(6,*) 'vm-W> UNDONE: symmetry is undefined'
	    return
	  endif
	 endif
	 num_op=num_symm
	endif

	cos_a=cosd(al)
	cos_b=cosd(be)
	cos_g=cosd(ga)
	sin_a=sind(al)
	sin_b=sind(be)
	sin_g=sind(ga)
	sin2_a=sin_a*sin_a
	sin2_b=sin_b*sin_b
	sin2_g=sin_g*sin_g

	volume=a*b*c*sqrt(1.-cos_a**2-cos_b**2-cos_g**2
	1                   +2.*cos_a*cos_b*cos_g) 
chk	end of cell volume calculation
	
	write(6,1001) mw/1000.0, volume,num_op
1001	format(' molecular weight (kDa) =',f10.1/
	1      ' cell volume      (A^3) =',g10.3/
	1      ' num of symmetry opt    =',i10/)

	vmth=volume/(mw*num_op)
	i=1
	do while(vmth/i .ge. 1.2 .and. i.le.20)
	  percent_sol=(1.-1.23*i/vmth)*100.0
	  num_water=int(volume/num_op*percent_sol*0.00033)	! assume water density 1 g/cm^3 -> 0.033 water mol/A^3
	  write(6,1088) vmth/i,i, nint(percent_sol), num_water
1088	    format(' VM=',f10.1,' A^3/Da for',i3,
	1' protein(s)/asu, (Vsolv=',i4,'%,',i5,' waters/asu)')
	  i=i+1
	  enddo
	return
900	return 1
	end

	subroutine init
chk	===============
	include 'edp_main.inc'
!	use edp_main
	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie
	logical ierr
	logical lf1( max_res)
	integer im(2)
	save	!030422

	n_of_syn=2					!000515
	syntax(1)='syntax:' 
	syntax(2)='initialize'

	do i=1,n_atom
	lf(i)=.false.
	end do
	return

	entry more(*)
c	=============
	n_of_syn=3					!000515
	syntax(1)='syntax:'
	syntax(2)='1) more [i0.i [i1.i]]' 
	syntax(3)='2) more chain '

	ierr=index(txt(:n_len),'from').le.0
	call dfgroup(igr,*101)

	im(1)=0
	im(2)=-999
	call read_ai(2,im, *103, *102)
102	if(im(2).eq.-999) then
	  im(2)=im(1)
!	else if(im(2).lt.im(1)) then
!	  goto 101
	  endif

	do j= 1, n_res
	  lf1(j)= .false.
	  enddo
	
	im1=min(im(1),im(2))
	im2=max(im(1),im(2))
	j0=n_res
	j1=1

	if(.not.ierr) then
	  k1=n_groupa(igr)
	  j0=aa_seq(igroupa( 1,igr))+im(1)
	  j1=aa_seq(igroupa(k1,igr))+im(2)
	  do k=1,k1
	    j=aa_seq(igroupa(k,igr))
	    i0= max(  1, j+ im1,    j0)
	    i1= min( j1, j+ im2, n_res)
	    do i= i0, i1
	      lf1(i)= .true.
	      end do
	    enddo
	else
	  do 100 j=1,n_res
	    do i=ijk(j),ijk(j+1)-1
	      if( lf(i)) goto 90
	      end do
            goto 100
90	    j0=min(j,j0)
	    j1=max(j,j1) 
	    i0= max(        1, j+ im1, j0+im(1))
	    i1= min( j1+im(2), j+ im2,    n_res)
	    do i= i0, i1
	      lf1(i)= .true.
	      end do
100	    enddo
	  endif

	do i=1,n_atom
	  lf(i)=.false.
	  enddo
	do j=1,n_res
	  if(lf1(j)) then
	    i0=ijk(j)
	    i1=ijk(j+1)-1
	    do i=i0,i1
	      lf(i)=.true.
	      enddo
	    endif
	  enddo
	return

103	if(txt(ib:ie).ne.'chain') goto 101
	call more_chain(igr)
	return
101	return 1

	entry setb(*)
c	=============
	call read_ar(1,bv,*201,*202)
 	do i=1,n_atom
	  if(lf(i)) b(i)=bv
	end do
	return
201	return 1
202	write(6,*) 'setb> b=<b>'
	do i=1,n_res
	  kk=0
	  av=0.
	  do j=ijk(i), ijk(i+1)-1
	   if(lf(j)) then
	    av = av+b(j)
	    kk=kk+1
	   end if
	  end do
	  if(kk.gt.0) then
	   av=av/kk
	   do j=ijk(i), ijk(i+1)-1
	    if(lf(j)) b(j)=av
	   end do
	  end if
	end do
	return

	entry blank
c	===========
	n_of_syn=2					!000515
	syntax(1)='syntax:'
	syntax(2)='blank -- the same as set text 31 54 " "'

	do i=1, n_atom
	if(lf(i)) text(i)(31:54)=' '
	end do
	end

	subroutine writef(file0,ierr)
chk	=================
	include 'edp_main.inc'
!	use edp_main
	include 'edp_file.inc'
!	use edp_file
	include 'edp_dat.inc'
!	use edp_dat

	parameter (n_lab=3)	!	integer, parameter :: n_lab=3
	character*(16) label(n_lab)
	data label/'header','cos','blank'/

	character*(*) file0

	logical ap, same_file
c--
	character*3 scos, junk
	character*12 junk1
	character*6 sheader
	data scos/'cos'/, sheader/'header'/
	character*(108) file
	character*(32) form, output_format

	character*(80) header 
	real trnb(3,3), cell(12)

	character*(72) texti
	data output_format/'(A72)'/

	character*32 option
	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie
	logical nword	!, nword0
	
	character pdate*10,ptime*10

	integer np, npp
	character*8  pp
	character*72 pn
	common /ud_param/ npp(max_param), pp(max_param),
	1 np(max_param), pn(max_param)
	save	!030422

	n_of_syn=2					!000515
	syntax(1)='syntax:' 
	syntax(2)=
	1'write filename.s [(cos, header, blank, title)] [fortran_format.s]' 

	ierr=1
	file=file0
	pdb_out='?'
	call vu_open_file(4, file,'.pdb',ierr) 
	if(ierr.ge.3) goto 100		! open failure
	if(ierr.eq.2) goto 103		! blank
chk	if(ierr.eq.1) 			! comma ,

	n_l=ltrim(file)

50	same_file=.false.
	n_l0=index(file,'.pdb')-1
	if(n_l0.gt.0) then
	  if(file0 .eq. file(:n_l0)) same_file=.true.
	  endif

	ierr=1
	pdb_out=file
	ap=.false.

	goto 2000
c--
	entry append(jerr)
chk	============
	n_of_syn=2					!000515
	syntax(1)='syntax:' 
	syntax(2)='append [comment.s]' 

	if(pdb_out.eq.'?') goto 100
	same_file=.false.
	jerr=1

2000	form= output_format
	
	call find1(n_lab,label,i)
	
	option='unknown'
	if(i.eq.0) then
	  if(.not.same_file) then
	    if(ap) then
	      option='blank'
	    else
	      option='cos_or_blank'
	      endif
	  else
	    option='copy_header'
	    endif
	else if(i.eq.1) then
	  option='copy_header'
	else if(i.eq.2) then
	  option='cos'
	else if(i.eq.3) then
	  option='blank'
	else if(i.lt.0) then
	  option='user_input'
	else
	  errmsg=' errmsg: unknown error in writef'
	  return
	  endif

	if(option.eq.'copy_header') then
	  write(6,*) '{list of headers copied from the original file:'
	  rewind (1)
60	  read(1,'(a80)',err=102) header
	  if(header(1:5).ne.'ATOM '.and.header(1:4).ne.'HETA') then 
	    write(4,1051,err=101) header
1051	    format(a80)
	    write(6,1019) header
1019	    format(1x,a79)
	    goto 60
	    endif
	  write(6,*) 'end_of_list}'
	else if(option.eq.'blank') then
	  goto 75
	else if(option.eq.'cos_or_blank'.or.option.eq.'cos') then
	  call date_and_time(pdate,ptime)
1003	format('REMARK FILENAME="',a,'"'
	1     /'REMARK TIME:',a,'/',a,'/',a,' DATE: ',a,'/',a,'/',a	!,'       created by user: '
	1     /'REMARK EdPDB VERSION:',a4)

	  isymm=-3   ! to get cell dimension isymm<=-3   5-sep-1991
	  call get_trn( isymm,cell,trnb,junk)
	  if(isymm.le.-3) then
	    if(option.eq.'cos_or_blank') then
	write(4,1003) pdb_out(:ltrim(pdb_out)),
	1ptime(1:2),ptime(3:4),ptime(5:6), pdate(1:4),pdate(5:6),pdate(7:8)
	1, version
	if(verbose.ge.6 ) write(6,1004)
1004	format(
	1' writef-I6> no cell parameter is written to the output pdb file')
	      goto 75
	    else 
	      return
	      endif
	    endif
	  junk1=char(ichar(pn(max_reserved_param+2)(1:1))-32)//' '//
	1 pn(max_reserved_param+2)(2:np(max_reserved_param+2))
	  write(4,1003) pdb_out(:ltrim(pdb_out)),
	1ptime(1:2),ptime(3:4),ptime(5:6), pdate(1:4),pdate(5:6),pdate(7:8)
	1, version
	  write(4,1009,err=101)  (cell(i),i=1,6), junk1,
	1   1.,0.,0., 0.,1.,0., 0.,0.,1., ((trnb(i,j),j=1,3),i=1,3)
1009	  format('CRYST1', 3F9.3, 3F7.2, 1x, A12
	1/ 'ORIGX1',4X,3F10.6,8X,'0.00000',
	1/ 'ORIGX2',4X,3F10.6,8X,'0.00000',
	1/ 'ORIGX3',4X,3F10.6,8X,'0.00000',
	1/ 'SCALE1',4X,3F10.6,8X,'0.00000',
	1/ 'SCALE2',4X,3F10.6,8X,'0.00000',
	1/ 'SCALE3',4X,3F10.6,8X,'0.00000')
	else if(option.eq.'user_input') then
	  ic=ichar(delimiter)
	  if( ic .eq. ichar(txt(ib:ib))) then	! using '....' to input header
	    ib= ib+1
	    do j=ib, n_len
	      if( ichar(txt(j:j)).eq.ic) then
	        ie=j-1
		txt(j:j)=' '	!000303, it may screw up the .out file
                goto 70
	        endif
	      enddo
            ie=n_len
            endif	  
70	  if( ib .le. ie) write(4,1053) txt(ib:ie)
	  ie=ie+1
1053	format(a)
	else
	  goto 101
	  endif

chk	get format
75	if(nword(n_len,txt,ib,ie)) goto 90	! blank
	ie=n_len
	ic=ichar(delimiter)
	if( ichar(txt(ib:ib)).eq. ic) then !use '(...)' to input the user prefored output-format
	  ib= ib+1
	  do j=ib, n_len
	    if( ichar(txt(j:j)).eq.ic) then
	      txt(j:j)=' '	!000303, it may screw up the history file
	      ie=j-1
	      goto 80
	      endif
	    enddo
          endif	  
80	if(ib.ge.ie) then
	  errmsg=' errmsg: wrong format'
	  return
	  endif
	form= txt(ib:ie)
	nf=ie-ib+1
	write(6,*) form(:nf)//' format will be used.'

90	n=0
	do ii=1,n_atom
	  i=iorder(ii)
	  if(lf(i)) then
	    texti=text(i)
	    do ic=1,72
	      jc=ichar(texti(ic:ic)) 
	      if(97.le.jc.and.jc.le.122) texti(ic:ic)=char(jc-32)
	      enddo
c	    text(i)=texti
	    write(texti(55:66),1011,err=101) w(i),b(i)
1011	format(2f6.2)	    
	    write(4,form,err=101) texti
	    n=n+1
	    endif
	  enddo
	if(ap) then
	  jerr=0
	  write(6,1001) file(1:n_l),n
1001	  format(' writef-I> append on file : ' ,a,i6,' lines')
	else
	  ierr=0
	  write(6,1002) file(1:n_l),n
1002	  format('         new file : ' ,a,i6,' lines')
	  ap=.true.	
	  endif
	return
100	errmsg=' errmsg: error during openning the file'
	close(4)
	return
101	errmsg=' errmsg: error during writing the file'
	return
102	errmsg=' errmsg: error during reading the original file'
	return
103	errmsg=' errmsg: file specification required. [, for old-name]'
	close(4)
	return
	end
c--
	subroutine readf (ierr)
chk	================
	include 'edp_main.inc'
!	use edp_main
	include 'edp_file.inc'
!	use edp_file

	parameter (max_tmp=max_atom+256)
!	integer, parameter :: max_tmp=max_atom+256

 	character*1 jcm,kcm, nres*3
 	character*1 curr_cm, subs_cm
	character*4 key, curr_id
	character*5 atom_tmp, kres, jres
	character*10 keywords(1)
	data keywords/'initialize'/

	character*72 txt0,tmp(max_tmp)
	equivalence (tmp(1)(:1),text(1)(:1))
	equivalence (txt0(:1),key(:1))
	logical nreset
	data kcm/' '/,nreset/.true./

	logical tolower
	common /main_para/ tolower
	1 ,jnk(m_main_para-1)

	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie
	logical nword, nword0

	logical  open_file1
	external open_file1
	save	!030422

	ierr=0
!	write(6,*)

!	format('atom',3x,i4,2x,a4,a3,1x,a1,1x,i3,4x,3f8.3,2f6.2)
1001	format(12x,            a5,a3,1x,a1,      8x,3f8.3,2f6.2)

	je=0
	call read_records(1, max_tmp,tmp,je)

	call extent_chain_mark(cm1, ncm1,*901)
	call extent_chain_mark(cm2, ncm2,*901)
	icm=0
	curr_cm='_'

	do jj=1,je
!	  if(tmp(jj)(1:4).eq.'HETA') tmp(jj)(1:4)='ATOM'
 	  if(tmp(jj)(1:4).eq.'ATOM'.or.tmp(jj)(1:4).eq.'HETA') then
	    if(tmp(jj)(22:22).ne.curr_cm .or.
	1      tmp(jj)(23:26).lt.curr_id) then
	      curr_cm=tmp(jj)(22:22)
	      icm=mod(icm,ncm1)+1
	      subs_cm=cm1(icm:icm)
	      endif
	    curr_id=tmp(jj)(23:26)
	    if(subs_cm.ne.'_') tmp(jj)(22:22)=subs_cm	     
	    endif
	  enddo

	if(ncm2.gt.0) then
	  j1=je+1
	  call read_records(3,max_tmp,tmp,je)
c<cm2
	icm=0
	curr_cm='_'
	do jj=j1,je
!	  if(tmp(jj)(1:4).eq.'HETA') tmp(jj)(1:4)='ATOM'
	  if(tmp(jj)(1:4).eq.'ATOM'.or.tmp(jj)(1:4).eq.'HETA') then
	    if(tmp(jj)(22:22).ne.curr_cm .or.
	1      tmp(jj)(23:26).lt.curr_id) then
	      curr_cm=tmp(jj)(22:22)
	      icm=mod(icm,ncm2)+1
	      subs_cm=cm2(icm:icm)
	      endif
	    curr_id=tmp(jj)(23:26)
	    if(subs_cm.ne.'_') tmp(jj)(22:22)=subs_cm	     
	    endif
	  enddo
cm2>
	  endif
	n_res=0
	n_atom=0
	j1=1
	cm3=' '
	ncm3=0
	goto 50

	entry readf1(*)
chk	============
	n_of_syn=2					!000515
	syntax(1)='syntax:' 
	syntax(2)='read filename.s [mark.s] [initialize]' 

	pdb_in3='?'
	if(.not. open_file1( 9, pdb_in3, 'old','.pdb')) then
	  errmsg=' errmsg: error during opening the PDB file'
	  return 1
	  endif
	pdb_in3=txt(ib:ie)

	cm3='_'
	ncm3=1
	if(.not.nword0(n_len,txt,ib,ie)) then
	  if(ib .le. ie) then
	    cm3=txt(ib:ie)
	    ncm3=ie-ib+1
	    endif
	  call find1(1,keywords,i)
	  if(i.lt.0) then
	    close(9)
	    return 1
	  else if(i.eq.1) then	! initialize
	    n_atom=0
	    n_res=0
	    n_groupa(0)=0
	    endif
	  endif
	call extent_chain_mark(cm3, ncm3,*901)

	je=n_atom
	j1=je+1
	call read_records(9,max_tmp,tmp,je)

	icm=0
	curr_cm='_'
	do jj=j1,je
!	  if(tmp(jj)(1:4).eq.'HETA') tmp(jj)(1:4)='ATOM'
!	  write(6,*) tmp(jj)
!	  write(6,*) tmp(jj)(22:22),curr_cm
	  if(tmp(jj)(1:4).eq.'atom'.or.tmp(jj)(1:4).eq.'heta'
	1.or.tmp(jj)(1:4).eq.'ATOM'.or.tmp(jj)(1:4).eq.'HETA') then
	    if(tmp(jj)(22:22).ne.curr_cm) then
	      curr_cm=tmp(jj)(22:22)
	      icm=mod(icm,ncm3)+1
	      subs_cm=cm3(icm:icm)
!	      write(6,*) subs_cm,icm
	      endif
	    if(subs_cm.ne.'_') tmp(jj)(22:22)=subs_cm	     
	    endif
	  enddo
	write(6,*) 
	goto 50

901	errmsg=' errmsg: invalid chain mark string'
	close(9)
	return 1

50	num_res=	n_res
	num_atom=	n_atom

	kres= '****'
	im=2
	do 100 jj=j1,je
	  txt0=tmp(jj)
	  do ic=1,72
	    jc=ichar(txt0(ic:ic))
	    if((ic.le.6.or.tolower).and.(65.le.jc.and.jc.le.90)) 
	1     txt0(ic:ic)=char(jc+32)
	    enddo

	  if(key.ne.'atom'.and.key.ne.'heta') goto 100
	  if( num_atom.ge.max_atom-1) goto 400 
	  num_atom=num_atom+1
	  text(num_atom)=txt0

c1001	format(12x,            a5,a3,1x,a1,      8x,3f8.3,2f6.2)
	  read( txt0,1001, err=300) 
	1   atom_tmp,nres,jcm,x(num_atom),y(num_atom),z(num_atom)
	1  ,w(num_atom),b(num_atom)
	  if(atom_tmp(1:1).eq.' ') atom_tmp(1:4)=atom_tmp(2:5)
	  if(atom_tmp(4:4).ne.' ') then
	    if(atom_tmp(3:3).eq.' ') atom_tmp(3:3)='_'
	    if(atom_tmp(2:2).eq.' ') atom_tmp(2:2)='_'
	    endif
	  atom(num_atom)=atom_tmp(1:4)
	  ie=22
	  if(nword(30,txt0,ib,ie)) goto 350
          jres=txt0(ib:ie)
	  if(jres.ne.kres.or.jcm.ne.kcm)   then
	    kres=jres
	    kcm=jcm
	    if(num_res.ge.max_res-1) then 
	      num_atom=num_atom-1
	      goto 400
	    endif
	    num_res=num_res+1
	    res(num_res)=nres
	    if(jcm.eq.' ') then
	      ires(num_res)=jres
	    else
	      ires(num_res)=jcm(:1)//jres
	      endif
	    ijk(num_res)=num_atom
	    endif
	  lf(num_atom)=.false.
	  iorder(num_atom)=num_atom
	  aa_seq(num_atom)=num_res

	  if(nreset) then 
	    if(num_atom.gt.1000) then
	      jp=nint(100.0*float(num_atom)/float(je))
	      if(jp.ge.ir) then
	        do while (jp.ge.ir) 
	           ir=ir+10
	           enddo
c	        call monitor(jp,'read> ','%')
	        endif
	    else if(mod(num_atom,im).eq.0) then
	      im=im*2
c	      call monitor(num_atom,'read> ',' ')
	      ir=10
	      endif
	    endif
100	  enddo


101	continue
	if(num_atom.gt.0) then
	  write(6,1006)  num_atom-n_atom, num_res-n_res
1006	  format(
	1  ' read>',i6,' atoms,',i6,' res.') !  type help for help')	
c	1  '+read>',i6,' atoms,',i6,' res.') !  type help for help')	
	  nreset=.false.
	  ijk(num_res+1)=num_atom+1
	  n_atom=	num_atom
	  n_res=	num_res
	  zone_untouch		=.true.
	  atom_untouch		=.true.
	  residue_untouch	=.true.
	  close(9)
	  return
	  endif
	if(cm3.ne.' ') then
	  errmsg=' errmsg: zero record is read in'
	  n_err=max_err+1
	  close(9)
	  return 1
	  endif

300	write(6,1008) txt0
1008	format(' read> ',a)
c1008	format('+read> ',a)
350	if(cm3.eq.' ') ierr=1		! initial reading
	return

400	write(6,1009) 
	1 ' read-W> max_num_atom=',max_atom-1,', max_num_res=',max_res-1
	write(6,*) 
1009	format(a,i7,a,i7)
	goto 101
!	if(cm3.eq.' ') ierr=1		! initial reading
	end

	subroutine read_records(io,max_tmp,tmp,j_end)
chk	=======================
	character*72 tmp(max_tmp)
!	character*80 tmp(max_tmp)	! expanding text
	rewind (io)
	do while(.true.)
	  j_end= j_end+1
10	  read(io,1007,end=20,err=10) tmp(j_end)
	  if(j_end.ge.max_tmp) return
	  enddo
1007	  format(a72)
!1007	  format(a80)		! expanding text
20	j_end=j_end-1
	end

	subroutine extent_chain_mark(cm, ncm,*)
chk	============================
chk	extent the chain mark string.
chk	eg. from a-i to abcdefghi

	character*28 cm

100	i=index(cm,'-')
	if(i.le.0) return
	if(i.eq.1) return 1
	ic=ichar(cm(i-1:i-1))
	if(.not.(
	1  (97.le.ic.and.ic.le.122) .or. 
	1  (65.le.ic.and.ic.le.90)  .or.
	1  (48.le.ic.and.ic.le.57) ))return 1
	jc=ichar(cm(i+1:i+1))
	if(.not.(
	1  (97.le.jc.and.jc.le.122) .or. 
	1  (65.le.jc.and.jc.le.90)  .or.
	1  (48.le.jc.and.jc.le.57) ))return 1

	cm(i:)=cm(i+1:)
	ncm=ncm-1
	do j=ic+1,jc-1
	  if(ncm.ge.28) return 1
	  cm(i:)=char(j)//cm(i:)
	  ncm=ncm+1
	  i=i+1
	  enddo
	goto 100
	end

chk*** end of readf.for

copyright by X. Cai Zhang

	function open_file(io,fn,stt,dft)
chk	===================
chk	stt= 'old'		-- readonly
chk	     'unknown','new' 	-- normal

	include 'edp_dat.inc'
!	use edp_dat
	include 'edp_file.inc'
!	use edp_file
	include 'edp_dim.inc'
!	use edp_dim

	logical open_file
	character*(*) stt
	character*(*) dft
	character*(108) file_name
	character*(*)  fn

	character*(max_num_chars) txt
	common /cmm_txt/ n_len,txt,ib,ie
!	logical nword

	if(len(fn).le.0) then
	  open_file=.false.
	  return
	  endif

	file_name= fn(:ltrim(fn))
	close (io)
	if(index(file_name,'.').le.0)  then
	  if(ltrim(file_name).le.0) then
	    open_file=.false.
	    return
	    endif
	  file_name= fn(:ltrim(fn))//dft
	  endif
	
		call open_file0(io,file_name,stt,dft,fn,*900)

	open_file=.true.
	if(verbose.ge.2 ) then 
	  if(stt.eq.'old') then 
	    write(6,1501) file_name(:ltrim(file_name)) 
	  else
	    write(6,1503) file_name(:ltrim(file_name)) 
	    endif
	  endif 
1501	format(' open_file-I2> [',a,'] is opened for readonly.')
1503	format(' open_file-I2> [',a,'] is opened for writing.')
	return

900	open_file=.false.
	if(verbose.ge.2 ) write(6,1502) file_name(:ltrim(file_name)) 
1502	format(' open_file-W2> [',a,'] can not be opened')
	close(io)
	end

	function open_file1(io,file_name,stt,dft)
chk	===================
chk	stt= 'old'		-- readonly
chk	    'unknown','new' 	-- normal

	include 'edp_dim.inc'
!	use edp_dim

	logical open_file
	logical open_file1
	character*(*) stt
	character*(*) dft
	character*(108) file_name

	character*(max_num_chars) txt
	common /cmm_txt/ n_len,txt,ib,ie
	logical nword

	if(.not.nword(n_len,txt,ib,ie)) then
	  open_file1=open_file(io,txt(ib:ie),stt,dft)
	  file_name=txt(ib:ie)
	else ! no file_name or comma
	  open_file1=open_file(io,file_name,stt,dft)
	  endif
	end

	subroutine s_open_file(io,file_name,stt,dft,ierr)
chk	======================
chk	stt= 'old'		-- readonly
chk	    'unknown','new' 	-- normal
	include 'edp_dim.inc'
!	use edp_dim

	logical open_file

	character*(*) stt
	character*(*) dft
	character*(108) file_name

	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie
	logical nword0
	
	if(nword0(n_len,txt,ib,ie)) then
	  ierr=2
	  if(stt.eq.'new') return
	else if(ib.gt.ie) then	! comma
	  ierr=1
	else
	  ierr=0
	  file_name=txt(ib:ie)
	  endif

	if(.not. open_file(io, file_name, stt, dft)) ierr=3
	end

copyright by X. Cai Zhang

	subroutine regionr (k,i1,i2,*)
c	==============================
	include 'edp_main.inc'
!	use edp_main
	dimension i1(max_num_zones), i2(max_num_zones)
	logical  followed	
	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie
c
	k=0
	followed=.false.
500	i= match_id1(n_res,ires,followed)
	if(i.lt.0) then
	  return 1
	else if(i.gt.0) then
	  k=k+1
	  if(k.gt.max_num_zones) then
	    errmsg=' errmsg: too many zones'
	    return 1
	    endif
	  i1(k)=i
	  if(followed) then 
	    i= match_id1(n_res,ires,followed)
	    if(followed.or.i.le.0) return 1
	    endif
	  i2(k)=i
	  if(i2(k).lt.i1(k)) return 1
	  goto 500
	  endif
	end

	subroutine region (k,i1,i2,*)
c	==============================
	include 'edp_main.inc'
!	use edp_main
	dimension i1(max_num_zones),i2(max_num_zones)
	logical  followed	
	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie
c
	k=0
	followed=.false.
500	i= match_id1(n_res,ires,followed)

	if(i.lt.0) then
	  errmsg=' errmsg: illegal residue id'
	  return 1
	else if(i.gt.0) then
	  k=k+1
	  if(k.gt.max_num_zones) then
	    errmsg=' errmsg: too many zones'
	    return 1
	    endif
	  i1(k)=ijk(i)
	  if(followed) then 
	    i= match_id1(n_res,ires,followed)
	    if(followed.or.i.le.0) return 1
	    endif
	  i2(k)=ijk(i+1)-1
	  if(i2(k).lt.i1(k)) then
	    errmsg=' errmsg: illegal id order'
	    return 1
	    endif
	  goto 500
	  endif
	end

	function match_id1(n_lab,label,followed)
chk	==================
	include 'edp_dim.inc'
!	use edp_dim

	character*(*) label(n_lab)
	character*10 b
	data id_save/0/
	logical followed

	character*(max_num_chars)  txt
	common /cmm_txt/n_len,txt,ib,ie
	logical nword1
	
	if(id_save.le.0 .or. .not.followed ) then
	  m1=1
	else
	  m1=id_save
	  endif

	if( nword1(n_len,txt,ib,ie,followed)) then
	  match_id1=0			! no word
	  return
	  endif

	b=txt(ib:ie)
	do match_id1=m1, n_lab
	  if( b .eq. label(match_id1)) then
	    id_save= match_id1	! e.g. return value for 'a10'
	    return
	    endif
	  enddo

	if(b.eq.'last') then 	
	  match_id1=n_lab
	  id_save= match_id1		! e.g. return value for 'last'
	  return
	  endif

	if(b.eq.'first') then
	  match_id1=1
	  id_save= match_id1		! e.g. return value for 'first'
	  return
	  endif

	match_id1= -1
	i= index( b, '+')
	if( i .le. 0)  return		! error, illegal ID
	  
	ibi= ib+ i
	if(ibi.gt.ie) return		! error, illegal ID, e.g. 'a+'

	read( txt(ibi:ie), *, err=900) j
	if( i .eq. 1) then				
	  match_id1= id_save+ j		! e.g. '+10'
          if( match_id1 .gt. n_lab) match_id1=-1
	  return
	  endif

	b=txt(ib:ibi-2)
	if(b.eq.'_') then		! e.g. change '_+10' to  '10'
	  b=txt(ibi:ie)
	else			! e.g. change 'a+10' to 'a10'
	  b=txt(ib:ibi-2)//txt(ibi:ie)
	  endif
	do jj=1,n_lab
	  if( b .eq. label(jj)) then
	    id_save= jj-j		! set registor	of 'a10' as 10
	    match_id1= jj		! and return the value for 'a10'
	    return
	    endif
	  enddo
900	end

c--
	function nword1(n_len,txt,ib,ie,followed)
	logical nword1, followed
	character*(*) txt
	parameter (isp=32, itb=9, icm=44, ihf=45)

	if(ie.le.0) then
	  ib=1
	else
	  ib=ie+1
	end if

	i=ib
	do while(i.le.n_len)
	  j=ichar(txt(i:i))
	  if(j.eq.isp.or.j.eq.itb.or.j.eq.icm.or.j.eq.ihf) then
	    i=i+1	
	    goto 90
	  else 
 	    goto 100
	  endif
90	enddo
	nword1=.true.
	return

100	ib=i
	do i=ib,n_len
	  j=ichar(txt(i:i))
	  if(j.eq.isp.or.j.eq.itb.or.j.eq.icm.or.j.eq.ihf) goto 200
	  enddo

200	ie=i-1
	nword1=.false.
	followed=.false.
	i0=i
	do i=i0,n_len
	  j=ichar(txt(i:i))
	  if(j.eq.isp.or.j.eq.itb.or.j.eq.icm) then 
	    goto 210
 	  else if(j.eq.ihf) then 
	    followed=.true.
	    return
	  endif
 	  return
210	enddo
	end

copyright by X. Cai Zhang

	subroutine rotmomentum(*)
c	======================
	implicit real* 8 (a-h)
	include 'edp_main.inc'
!	use edp_main
	include 'edp_file.inc'
!	use edp_file

!	implicit real* 8 (a-h)
	real xc(3), pc(3), i1,i2,i3, soln(3)
	dimension a(3,3)

	logical  open_file1
	external open_file1

c	logical  punch_vec
c	integer pseudo_id

	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie

	n_of_syn=3					!000515
	syntax(1)='syntax:' 
	syntax(2)='momentinertia [filename.s]' 
	syntax(3)='   [vector_id1.s] [vector_id2.s] [vector_id3.s]' 

	rtn_out='rtn_.txt'
	if(.not. open_file1(29,rtn_out,'unknown','.txt')) return 1
	  
	iv1=0
	iv2=0
	iv3=0
	call get_vector_id(iv1, *200)
	call get_vector_id(iv2, *200)
	call get_vector_id(iv3, *200)

200	do i=1,3
	do j=1,3
	  a(j,i)=0.
	  enddo
	  pc(i)=0.
	  enddo
	i1=0.
	i2=0.
	i3=0.

	j1=0
	sum_w=0.
	rmax=0.
	do i=1,n_atom
	  if(lf(i)) then
	    sum_w=sum_w + w(i)
	    rmax=max(rmax,abs(x(i)),abs(y(i)),abs(z(i)))
	    if(j1.le.0) j1=i
	    j2=i
	    endif
	  enddo
	if(sum_w.le.0.) then 
	  write(6,*) 'momentinertia-W> UNDONE:  mwt.le.0.'
	  return
	  endif

	scale=1. /(sum_w*rmax*rmax)
	sum_w=sum_w*scale

	do i=1,n_atom
	  if(lf(i)) then
	    wi=w(i)*scale
	    pc(1)=pc(1) + x(i)*wi
	    pc(2)=pc(2) + y(i)*wi
	    pc(3)=pc(3) + z(i)*wi
	    endif
	  enddo
	pc(1)= pc(1)/sum_w
	pc(2)= pc(2)/sum_w
	pc(3)= pc(3)/sum_w

chk	this loop takes time to run
	do i=1, n_atom
	  if(lf(i)) then
	    xc(1)=x(i)-pc(1)
	    xc(2)=y(i)-pc(2)
	    xc(3)=z(i)-pc(3)
	    wi=w(i)*scale
	    do j=1,3
	    do k=1,3
	      a(k,j)=a(k,j) - xc(k)*xc(j)*wi
	      enddo
	      enddo
	    i1=i1+ (y(i)*y(i)+z(i)*z(i))*wi
	    i2=i2+ (z(i)*z(i)+x(i)*x(i))*wi
	    i3=i3+ (x(i)*x(i)+y(i)*y(i))*wi
	    endif
	  enddo

	a3=a(1,1)+a(2,2)+a(3,3)
	a(1,1)=a(1,1)-a3
	a(2,2)=a(2,2)-a3
	a(3,3)=a(3,3)-a3

	write(6,1001) 	i1/scale, sqrt(i1/sum_w),
	1         	i2/scale, sqrt(i2/sum_w),
	1         	i3/scale, sqrt(i3/sum_w)
1001	FORMAT(' the current   momentinertia and gyration radii:'
	1     /'    ix=',e12.4,', rx=',f8.3
	1     /'    iy=',e12.4,', ry=',f8.3
	1     /'    iz=',e12.4,', rz=',f8.3)

	xc(1)=x(j2)-x(j1)
	xc(2)=y(j2)-y(j1)
	xc(3)=z(j2)-z(j1)
	call rotm(a,pc,sum_w,scale,xc,soln)

	if(iv1.gt.0) then
	  vectors(1,iv1)=pc(1)
	  vectors(2,iv1)=pc(2)
	  vectors(3,iv1)=pc(3)
	  vectors(4,iv1)=a(3,1)
	  vectors(5,iv1)=a(3,2)
	  vectors(6,iv1)=a(3,3)
	  vectors(7,iv1)= soln(3)
	  endif
	if(iv2.gt.0) then
	  vectors(1,iv2)=pc(1)
	  vectors(2,iv2)=pc(2)
	  vectors(3,iv2)=pc(3)
	  vectors(4,iv2)=a(2,1)
	  vectors(5,iv2)=a(2,2)
	  vectors(6,iv2)=a(2,3)
	  vectors(7,iv2)= soln(2)
	  endif
	if(iv3.gt.0) then
	  vectors(1,iv3)=pc(1)
	  vectors(2,iv3)=pc(2)
	  vectors(3,iv3)=pc(3)
	  vectors(4,iv3)=a(1,1)
	  vectors(5,iv3)=a(1,2)
	  vectors(6,iv3)=a(1,3)
	  vectors(7,iv3)= soln(1)
	  endif	  
	return
c900	return 1
	end

	subroutine rotm(a,pc,sum_w,scale,xv,soln)
chk	===============
	implicit real*8 (a-h)
	include 'edp_file.inc'
!	use edp_file
	include 'edp_dat.inc'
!	use edp_dat
	
	dimension a(3,3),e(3),v(3,3)

	real pc(3), xv(3), soln(3)

	call eigen(a,e,v,.true.,*900)
	
	soln(1)=sqrt(e(1)/sum_w)
	soln(2)=sqrt(e(2)/sum_w)
	soln(3)=sqrt(e(3)/sum_w)
	write(6,1001) 
	1    e(1)/scale, soln(1),
	1    e(2)/scale, soln(2),
	1    e(3)/scale, soln(3)
1001	format(' the principle momentinertia and gyration radii:'
	1     /'    i1=',e12.4,', r1=',f8.3
	1     /'    i2=',e12.4,', r2=',f8.3
	1     /'    i3=',e12.4,', r3=',f8.3)

	if(xv(1)*v(3,1)+ xv(2)*v(3,2)+ xv(3)*v(3,3) .lt.0.) then
	  v(2,1)=-v(2,1)
	  v(2,2)=-v(2,2)
	  v(2,3)=-v(2,3)
	  v(3,1)=-v(3,1)
	  v(3,2)=-v(3,2)
	  v(3,3)=-v(3,3)
	  endif
c	xv(1)=v(3,1)
c	xv(2)=v(3,2)
c	xv(3)=v(3,3)

	a(1,1)=v(1,1)
	a(1,2)=v(1,2)
	a(1,3)=v(1,3)
	a(2,1)=v(2,1)
	a(2,2)=v(2,2)
	a(2,3)=v(2,3)
	a(3,1)=v(3,1)
	a(3,2)=v(3,2)
	a(3,3)=v(3,3)

249	write(29,118)a(1,1),a(1,2),a(1,3),
	1            a(2,1),a(2,2),a(2,3),
	1            a(3,1),a(3,2),a(3,3),
	1   -(a(1,1)*pc(1)+a(1,2)*pc(2)+a(1,3)*pc(3)),
	1   -(a(2,1)*pc(1)+a(2,2)*pc(2)+a(2,3)*pc(3)),
	1   -(a(3,1)*pc(1)+a(3,2)*pc(2)+a(3,3)*pc(3))

	write(29,*)'  matrix of diagonalization of momentinertia.'
118	format(3(3f12.7/),3f12.5)
	if(verbose.ge.3) write(6,1230) rtn_out(:ltrim(rtn_out))
1230	format(' rtn-I3> To apply the transformation, type:'
	1 /'   rtn file ',a)
	close(29)
	return
900	end 

copyright by X. Cai Zhang
	subroutine edp_rtn(*)
chk	==================
c
c	this program is for coordinate translation 
c
	include 'edp_main.inc'
!	use edp_main
	include 'edp_file.inc'
!	use edp_file
	include 'edp_dat.inc'
!	use edp_dat

	real e(3),p1(3),p2(3), eu(12)
	1 ,arr(3,3)  ,arr2(3,3)
	1 ,trn0(3,3) ,trn1(3,4) ,trn2(3,4)

	integer atoms_id(4)

	character*32 tmp  ,yes*1
	character*10 curr_date, curr_time
	character*32 symm_txt
	character*(108) rtnfil
	character*6 key1
	character*4 key2
	character*80 header

	character*4 ssave,sinv,smult
	data ssave,sinv,smult/'save','inve','mult'/

	logical  open_file
	external open_file
	logical  open_file1
	external open_file1

 	character*(max_num_chars) txt
	logical nword0
	common /cmm_txt/n_len,txt,ib,ie

        logical tolower
        common /main_para/ tolower, jou0, echo_L, inter0, incl0, incl1
	1 ,jou
	1 ,jnk1(m_main_para-7)

	integer pass
	logical fcell
	data fcell,iii/.true.,1/

	parameter (num_options=17)
!	integer, parameter :: num_options=16
	character*12 options(num_options)
	data options/
	1 'file','symmetry','ezxz','ezyz','deorth','orthog','axis',	!7
	1 'overlay','polar','match','abcd','center','matrix',		!13
	1 'c2sp','sp2c','v_align','remark'/				!17

	n_of_syn=12					!000515
	syntax(1)='syntax:' 
	syntax(2)='rtn abcd res_a.s atom_a.s res_b.s atom_b.s' 
	syntax(3)='  res_c.s atom_c.s res_d.s atom_d.s torsion_angle.r' 
	syntax(4)='rtn axis vector_id.s rotation_angle.r [translation.r]' 
	syntax(5)='rtn (deorth,orthog) grid_a.r grid_b.r grid_c.r'
	syntax(6)='rtn (ezxz,ezyz,polar) e1.r e2.r e3.r [t1.r t2.r t3.r]' 
	syntax(7)='rtn file filename.s' 
	syntax(8)='rtn overlay res_id1.s'
 	syntax(9)=
	1'  [atom11.s atom12.s atom13.s] [reg11.s reg12.s reg13.s]'
	syntax(10)=
	1'  [res_id2.s [atom21.s atom22.s atom23.s] [reg21.s reg22.s reg23.s]]' 
	syntax(11)='rtn symmetry [symm_#.i] [fx.r fy.r fz.r]' 
	syntax(12)='rtn remark (BIOMT,SMTRY) Op#' 

	errmsg=' errmsg: insufficient or wrong input'
	rtnfil='rtn_.txt'

	pass=match_l1(num_options,options)
	if( pass.le.0) return 1

	trn1(1,4)=0.
	trn1(2,4)=0.
	trn1(3,4)=0.
	goto (100 ,200 ,300 ,300 ,500 ,600 ,700, 800,900,2000
	1 ,2100,2200,2300,2400,2500,2600,2700)  pass

100	rtnfil='rtn_.txt'
	if(.not.open_file1(29, rtnfil, 'old', '.txt')) return 1
	rtn_in=rtnfil
	read(29,*,err=105, end=105) 
	1  ((trn1(i,j),j=1,3),i=1,3), (trn1(i,4),i=1,3)
	close (29)
	rtnfil=' '
	goto 40

105	errmsg=' errmsg: error during reading '//rtnfil
	return 1

200	eu(1)=0.0
	eu(2)=0.0	!a-shift
	eu(3)=0.0	!b-shift
	eu(4)=0.0	!c-shift

	call read_ar(4,eu,*901,*201)
201	if(nint(eu(1)).eq.0) then	! non standard translation.
	  call nomorlize(eu(2))
	  return
	  endif

	jsymm=nint(abs(eu(1)))
	isymm=jsymm-1
	call get_trn( isymm,trn1,trn0,symm_txt)
	if(isymm.ne.jsymm) return 1	
	if(eu(1).lt.0.) then			! inverse matrix
	  call mxinv(3,trn1,ierr)
	  if(ierr.ne.0) then
	    errmsg=' errmsg: error during inversing the matrix.'
	    return 1 
	    endif
	  e1=-trn1(1,1)*trn1(1,4)-trn1(1,2)*trn1(2,4)-trn1(1,3)*trn1(3,4)
	  e2=-trn1(2,1)*trn1(1,4)-trn1(2,2)*trn1(2,4)-trn1(2,3)*trn1(3,4)
	  e3=-trn1(3,1)*trn1(1,4)-trn1(3,2)*trn1(2,4)-trn1(3,3)*trn1(3,4)
	  trn1(1,4)=e1
	  trn1(2,4)=e2
	  trn1(3,4)=e3
	  endif

	trn1(1,4)= trn1(1,4)
	1 +trn0(1,1)*eu(2) +trn0(1,2)*eu(3) +trn0(1,3)*eu(4)
	trn1(2,4)= trn1(2,4)
	1 +trn0(2,1)*eu(2) +trn0(2,2)*eu(3) +trn0(2,3)*eu(4)
	trn1(3,4)= trn1(3,4)
	1 +trn0(3,1)*eu(2) +trn0(3,2)*eu(3) +trn0(3,3)*eu(4)
	goto 40
c---
300	do i=1,6
	  eu(i)=0.
	  enddo
	call read_ar(6,eu,*901,*301)
301	z1=eu(1)
	z2=eu(2)
	z3=eu(3)
	trn1(1,4)=eu(4)
	trn1(2,4)=eu(5)
	trn1(3,4)=eu(6)

	call zrot(z1,trn1)
	if( pass.eq.3)call xrot(z2,arr2)
	if( pass.eq.4)call yrot(z2,arr2)
	call zrot(z3,trn0)
	call axbeqc(trn1,arr2,arr)
	call axbeqc(arr,trn0,trn1)
c***
	goto 40

c	default: anstrongs along the cell edges.
500	call read_ar(3,eu,*901,*901)
	isymm=-2
	call get_trn( isymm,trn0,trn1,symm_txt)
	if(isymm.le.-2) return

	s1=eu(1)
	s2=eu(2)
	s3=eu(3)
	write(6,1099) 
	1 ' rtn> translate from orthergonal'
	2,' to frn.(or grd.)_coordinates'
 
	do i=1,3
	trn1(1,i)=trn1(1,i)*s1
	trn1(2,i)=trn1(2,i)*s2
	trn1(3,i)=trn1(3,i)*s3
	end do
	goto 40

c	default: anstrongs along the cell edges.
600	call read_ar(3,eu,*901,*901)
	isymm=-1
	call get_trn( isymm,trn0,trn1,symm_txt)
	if(isymm.le.-1) return

	s1=1./eu(1)
	s2=1./eu(2)
	s3=1./eu(3)
	write(6,1099) 
	1 ' rtn> translate from frn.(or grd.)_coordinates'
	2,' to orthergonal'
 
	do i=1,3
	trn1(i,1)=trn1(i,1)*s1
	trn1(i,2)=trn1(i,2)*s2
	trn1(i,3)=trn1(i,3)*s3
	end do
	goto 40

700	call rtn_axis(trn1,*901)
	goto 40

800	call ncac(trn1, *901)
	goto 40

901	return 1

900	do i=1,6
	  eu(i)=0.
	  enddo
	call read_ar(6,eu,*901,*902)
902	z1=eu(1)
	z2=eu(2)
	z3=eu(3)
	trn1(1,4)=eu(4)
	trn1(2,4)=eu(5)
	trn1(3,4)=eu(6)
	al=sind(z2)*cosd(z1)
	am=sind(z2)*sind(z1)
	an=cosd(z2)
	trn1(1,1)=al*al+(am*am+an*an)*cosd(z3)
	trn1(1,2)=al*am*(1.-cosd(z3))-an*sind(z3)
	trn1(1,3)=al*an*(1.-cosd(z3))+am*sind(z3)
	trn1(2,1)=al*am*(1.-cosd(z3))+an*sind(z3)
	trn1(2,2)=am*am+(al*al+an*an)*cosd(z3)
	trn1(2,3)=am*an*(1.-cosd(z3))-al*sind(z3)
	trn1(3,1)=al*an*(1.-cosd(z3))-am*sind(z3)
	trn1(3,2)=am*an*(1.-cosd(z3))+al*sind(z3)
	trn1(3,3)=an*an+(al*al+am*am)*cosd(z3)
	goto 40

2000	call get_atoms(2, atoms_id, ierr) 
	if(ierr.ne.0) return 1
	ia=atoms_id(1)
	ja=atoms_id(2)
	do i=1,3
	  eu(i)=0.
	  end do
	call read_ar(3,eu,*901,*2001)
2001	z1=eu(1)
	z2=eu(2)
	z3=eu(3)
	al=sind(z2)*cosd(z1)
	am=sind(z2)*sind(z1)
	an=cosd(z2)
	trn1(1,1)=al*al+(am*am+an*an)*cosd(z3)
	trn1(1,2)=al*am*(1.-cosd(z3))-an*sind(z3)
	trn1(1,3)=al*an*(1.-cosd(z3))+am*sind(z3)
	trn1(2,1)=al*am*(1.-cosd(z3))+an*sind(z3)
	trn1(2,2)=am*am+(al*al+an*an)*cosd(z3)
	trn1(2,3)=am*an*(1.-cosd(z3))-al*sind(z3)
	trn1(3,1)=al*an*(1.-cosd(z3))-am*sind(z3)
	trn1(3,2)=am*an*(1.-cosd(z3))+al*sind(z3)
	trn1(3,3)=an*an+(al*al+am*am)*cosd(z3)
	if(ierr.ne.0) return 1
	trn1(1,4)=x(ja) - (trn1(1,1)* x(ia)+trn1(1,2)* y(ia)+trn1(1,3)* z(ia))
	trn1(2,4)=y(ja) - (trn1(2,1)* x(ia)+trn1(2,2)* y(ia)+trn1(2,3)* z(ia))
	trn1(3,4)=z(ja) - (trn1(3,1)* x(ia)+trn1(3,2)* y(ia)+trn1(3,3)* z(ia))
	goto 40

2100	call skrew0(trn1,*901)
	goto 40

2200	e(1)=0.
	e(2)=0.
	e(3)=0.
	k=0
	do i=1,n_atom
	  if(lf(i)) then
	    e(1)=e(1)+x(i)
	    e(2)=e(2)+y(i)
	    e(3)=e(3)+z(i)
	    k=k+1
	    endif
	  enddo
	if(k.le.0) then
	  errmsg=' errmsg: no ON atom exists.'
	  return 1
	  endif
	do i=1,3
	do j=1,3
	  trn1(i,j)=0.
	  enddo
	  trn1(i,i)=1.
	  trn1(i,4)=-e(i)/k
	  enddo
	goto 40

2300	ie0=ie
	call read_ar(12,eu,*23001,*23001)
	trn1(1,1)=eu(1)
	trn1(1,2)=eu(2)
	trn1(1,3)=eu(3)
	trn1(2,1)=eu(4)
	trn1(2,2)=eu(5)
	trn1(2,3)=eu(6)
	trn1(3,1)=eu(7)
	trn1(3,2)=eu(8)
	trn1(3,3)=eu(9)
	trn1(1,4)=eu(10)
	trn1(2,4)=eu(11)
	trn1(3,4)=eu(12)
	goto 40

23001	ie=ie0
	goto 100

2400	do i=1,n_atom
	  if(lf(i)) then
	    cx=x(i)
	    cy=y(i)
	    cz=z(i)
	    x(i)=sqrt(cx*cx+cy*cy+cz*cz)
	    if(x(i).le.1.e-6) then
	      y(i)=0.0
	      z(i)=0.0
	    else
	      y(i)=acosd(cz/x(i))
	      r0=sqrt(cx*cx+cy*cy)
              if(r0.le.1.e-6) then
	        z(i)=0.0
	      else
	        z(i)=atand2(cx,cy)		!?????
	        endif
	      endif
	    endif
	  enddo
	return

2500	do i=1,n_atom
	  if(lf(i)) then
	    cx=x(i)
	    cy=y(i)
	    cz=z(i)
	    r0=sind(cy)*cx
	    x(i)=cosd(cz)*r0
	    y(i)=sind(cz)*r0
	    z(i)=cosd(cy)*cx
	    endif
	  enddo
	return

2600	call find_vector_id(iv1,*900,*900)
	call find_vector_id(iv2,*900,*900)
	call ixjeqk4(vectors(4,iv2),vectors(4,iv1),p1,*900)
	call ixjeqk4(vectors(4,iv1),p1,p2,*900)
	do i=1,3
	  trn1(1,i)=p1(i)
	  trn1(2,i)=p2(i)
	  trn1(3,i)=vectors(3+i,iv1)
	  enddo
	do i=1,3
	  trn1(i,4)= -trn1(i,1)*vectors(1,iv1)
	1            -trn1(i,2)*vectors(2,iv1)
	1            -trn1(i,3)*vectors(3,iv1)
	  enddo
	goto 40

2700	if(nword0(n_len,txt,ib,ie)) return 1
	do ic=ib,ie
	      jc=ichar(txt(ic:ic)) 
	      if(97.le.jc.and.jc.le.122) txt(ic:ic)=char(jc-32)
	      enddo
	key1=txt(ib:ie)
	call read_ai(1,key3,*901,*901)
	write(key2,'(i4)') key3
c	write(6,*) key1,key2
	  rewind (1)
2760	  read(1,'(a80)',err=2702,end=2703) header
	  if(header(1:6).eq.'REMARK') then 
c	write(6,*) header(1:23)
	  if(header(14:18).eq.key1) then
	  if(header(20:23).eq.key2) then
	  if(header(19:19).eq.'1') then
		read(header,'(23x,3f10.1,5x,f10.1)',err=2702) (trn1(1,i),i=1,4)
	  else if(header(19:19).eq.'2') then 
		read(header,'(23x,3f10.1,5x,f10.1)',err=2702) (trn1(2,i),i=1,4)
	  else if(header(19:19).eq.'3') then
		read(header,'(23x,3f10.1,5x,f10.1)',err=2702) (trn1(3,i),i=1,4)
		goto 40
	  else 
		goto 2702
	  endif
		endif
		endif
		endif
	goto 2760
2702	write(6,*) 'error during reading the matrix from the PDB file.'
	return 1 
2703	write(6,*) 'I can not find the keyword ',key1,' with opertor # ',key2
	return 1 

40	write(6,1031) 
	1 'x', (trn1(1,i),i=1,4),
	1 'y', (trn1(2,i),i=1,4),
	1 'z', (trn1(3,i),i=1,4)
1031	format(' rtn> ',
	1    3(t8,a1,'new=',f8.3,'*xold+',f8.3,'*yold+',f8.3,'*zold+',
	2      f8.3,:/))

	call trnxyz(trn1,*901)

	if(nword0(n_len,txt,ib,ie)) return

	if(ib  .gt.ie .or. index(ssave,txt(ib:ie)).eq.1
	1    .or.txt(ib:ib).eq.'y') then
	    yes='y'
	else if( index(sinv,txt(ib:ie)).eq.1) then
	    yes='i'
	else if( index(smult,txt(ib:ie)).eq.1) then
	    yes='m'
	else
	    return 1
	    endif

	if(yes.eq.'m') then
	  rtn_in='rtn_.txt'
	  if(.not.open_file1(29, rtn_in, 'old', '.txt')) goto 901
	  read(29,*,err=105, end=105) 
	1   ((arr(i,j),j=1,3),i=1,3), (p1(i),i=1,3)
	  close (29)

	  do i=1,3
	    p2(i)=trn1(i,1)*p1(1)+trn1(i,2)*p1(2)+trn1(i,3)*p1(3)+trn1(i,4)
	    enddo
	  call axbeqc(trn1,arr,trn0)
	  do i=1,3
	    do j=1,3
	      trn1(i,j)=trn0(i,j)
	      enddo
	    trn1(i,4)=p2(i)
	    enddo
	  yes='y'
	  endif

	rtn_out='rtn_.txt'
	if(.not.open_file1(29, rtn_out, 'unknown', '.txt')) goto 901

	do i=1,4
	do j=1,3
	  trn2(j,i)=trn1(j,i)
	  enddo
	  enddo
	call mxinv(3,trn1,ierr)

	if(ierr.eq.0) then
	  tmp(1:10)=' inv:o.k. '
	  e(1)=-trn1(1,1)*trn1(1,4)-trn1(1,2)*trn1(2,4)-trn1(1,3)*trn1(3,4)
	  e(2)=-trn1(2,1)*trn1(1,4)-trn1(2,2)*trn1(2,4)-trn1(2,3)*trn1(3,4)
	  e(3)=-trn1(3,1)*trn1(1,4)-trn1(3,2)*trn1(2,4)-trn1(3,3)*trn1(3,4)
	else
	  tmp(1:10)=' inv:err  '
	  endif

	iii=0
	fcell=yes.eq.'y'
854	if(fcell) then
	  do i=1,3
	  do j=1,4
	    if(abs(trn2(i,j)).lt.eps()) trn2(i,j)=0.0
	    enddo
	    enddo
	  write(29,1023) ((trn2(i,j),j=1,3),i=1,3),(trn2(i,4),i=1,3)
	  if(iii.ge.1) !990618
	1   write(29,1025) ((trn1(i,j),i=1,3),j=1,3),e
	else
	  do i=1,3
	  do j=1,3
	    if(abs(trn1(i,j)).lt.eps()) trn1(i,j)=0.0
	    enddo
	    if(abs(e(i)).lt.eps()) e(i)=0.0
	    enddo
	  write(29,1023) ((trn1(i,j),j=1,3),i=1,3),e
	  if(iii.ge.1) 	!990618
	1   write(29,1025) ((trn2(i,j),i=1,3),j=1,3),(trn2(i,4),i=1,3)
	  endif

	if(iii.eq.0) then
	  iii=iii+1
	  fcell=.not.fcell
	  call date_and_time(curr_date, curr_time)
	  tmp(21:30)=curr_time
	  tmp(11:20)=curr_date
	  write(29,1099) 
	1' program EdPDB--rtn, date: '//
	1tmp(11:14)//'-'//tmp(15:16)//'-'//tmp(17:18)//', time: '//
	1tmp(21:22)//':'//tmp(23:24)//':'//tmp(25:26)
	  goto 854
	  end if
	if(verbose.ge.3) write(6,1230) rtn_out(:ltrim(rtn_out))
1230	format(' rtn-I3> To apply the transformation, type:'
	1 /'   rtn file ',a)
	close(29)

1023	format(3e16.8)
1025	format(' OMAT'/(3f12.6))
1099	format(a,a)
	end

        subroutine ixjeqk4(g1,g2,g3,*)
chk     ==============
        real g1(3), g2(3), g3(3)
        g3(1)=g1(2)*g2(3)-g1(3)*g2(2)
        g3(2)=g1(3)*g2(1)-g1(1)*g2(3)
        g3(3)=g1(1)*g2(2)-g1(2)*g2(1)
	r=sqrt(g3(1)**2+ g3(2)**2+ g3(3)**2)
	if(r.le.eps()) return 1
	do i=1,3
	  g3(i)=g3(i)/r
	  enddo 
        end

	subroutine rtn_axis(trn,*)
chk	===================
	include 'edp_main.inc'
!	use edp_main
	real trn(3,4), trn0(3,3), trn1(3,3), arr(3,3)
	real cx(3), cy(3), cz(3)
		
	call find_vector_id(iv,*900,*900)

	call read_ar(1,phi,*900,*900)
	t=0.0
	call read_ar(1,t  ,*900,*800)

800	call zrot(phi,trn)

	cz(1)=vectors(4,iv)
	cz(2)=vectors(5,iv)
	cz(3)=vectors(6,iv)

	if(cz(1).eq.0.) then
	  cy(1)=1.
	  cy(2)=0.
	  cy(3)=0.
	else
	  r=1./sqrt(cz(1)*cz(1)+cz(2)*cz(2))
	  cy(1)=cz(2)*r
	  cy(2)=-cz(1)*r
	  cy(3)=0.
	  end if
	call uxveqw(cy,cz,cx)
	do i=1,3
	  trn0(1,i)=cx(i)
	  trn0(2,i)=cy(i)
	  trn0(3,i)=cz(i)
	  trn1(i,1)=cx(i)
	  trn1(i,2)=cy(i)
	  trn1(i,3)=cz(i)
	  enddo
	call axbeqc(trn1,trn,arr)
	call axbeqc(arr,trn0,trn)
	
	cz(1)=vectors(1,iv)
	cz(2)=vectors(2,iv)
	cz(3)=vectors(3,iv)

	do i=1,3
	  trn(i,4)=cz(i)
	1         -trn(i,1)*cz(1)
	1         -trn(i,2)*cz(2)
	1         -trn(i,3)*cz(3)
	1        +trn1(i,3)*t
	  enddo
	return
900	return 1
	end

	subroutine nomorlize(center)
chk	====================		(930605)
chk	convert individual atom to one box cell, centered at the "center" 
chk	(in fractional)
chk	both the input and output coordinates are cartetional. 
	include 'edp_main.inc'
!	use edp_main
	real center(3), trn0(3,3), trn1(3,3),junk(3,4)
	character*32 junk_txt

	isymm=-2	! get matrix from cartetianal to fractional
	call get_trn( isymm,junk,trn0, junk_txt)
	if(isymm.le.-2) return

	isymm=-1	! get matrix from fractional to cartetianal
	call get_trn( isymm,junk,trn1, junk_txt)
	if(isymm.le.-1) return

	sx=center(1)-0.5
	sy=center(2)-0.5
	sz=center(3)-0.5

	do i=1, n_atom
	  if(lf(i)) then
	    cx=x(i)
	    cy=y(i)
	    cz=z(i)
	    ca= mod(trn0(1,1)*cx+ trn0(1,2)*cy+ trn0(1,3)*cz -sx,1.)
	    cb= mod(trn0(2,1)*cx+ trn0(2,2)*cy+ trn0(2,3)*cz -sy,1.)
	    cc= mod(trn0(3,1)*cx+ trn0(3,2)*cy+ trn0(3,3)*cz -sz,1.)
	    ca= mod(ca+1.,1.) +sx
	    cb= mod(cb+1.,1.) +sy
	    cc= mod(cc+1.,1.) +sz
	    x(i)= trn1(1,1)*ca+ trn1(1,2)*cb+ trn1(1,3)*cc 
	    y(i)= trn1(2,1)*ca+ trn1(2,2)*cb+ trn1(2,3)*cc 
	    z(i)= trn1(3,1)*ca+ trn1(3,2)*cb+ trn1(3,3)*cc
	    write(text(i)(31:54),1051,err=900) x(i),y(i),z(i)
1051	    format(3f8.3)
	    endif
	  enddo
	return
900	write(6,*) 'nomorlize-W> ABORTED: Output error.'
	end 

	subroutine get_atoms(howmany,atoms_id,ierr)
chk	====================
chk	ierr = 1 wrong residue_id 
chk	       2 blank
chk	       0 okay.

	include 'edp_main.inc'
	include 'edp_dat.inc'
	character atom_name*4
	integer howmany, atoms_id(howmany)

	character*(max_num_chars) txt
	logical nword0
	common /cmm_txt/n_len,txt,ib,ie

	ierr=1
	iatom=1

40	call find(n_res,ires,i)
	if(i.eq.0) then
	  ierr=2
	  return
	else if(i.lt.0) then
	  errmsg=' errmsg: wrong residue_id'
	  return 
	  end if

	if(nword0(n_len,txt,ib,ie)) goto 901

	if(ib  .gt.ie) then
	  j=ijk(i)
	  if(verbose.ge.4)    write(6,*)
	1 'get_atom-I4> the first atom of residue '//ires(i)//'is used.'
	  goto 50
	  endif

	atom_name=txt(ib:ie)
	do j=ijk(i),ijk(i+1)-1
	  if(atom(j).eq.atom_name) goto 50
	  end do
	errmsg=' errmsg: atom_name not found'
	return 

901	if(howmany .eq. 1) then
	  j=ijk(i)
	  if(verbose.ge.6) write(6,*) 
	1 'get_atom-I6> the first atom of residue '//ires(i)//'is used.'
	  goto 50
	else
	  errmsg=' errmsg: an atom_name is required'
	endif
902	return

50	atoms_id(iatom)=j
	iatom=iatom+1
	if(iatom.le.howmany) goto 40
	ierr=0
	end

	subroutine skrew0(trn,*)
chk	=================
	include 'edp_main.inc'
!	use edp_main
	include 'edp_dat.inc'
!	use edp_dat

	character*(max_num_chars) txt
	real trn(3,4), cx(3), cy(3), cz(3), p1(3), p2(3),
	1 trn1(3,3), trn0(3,3), arr(3,3) 
	integer atoms_id(4)
	common /cmm_txt/n_len,txt,ib,ie

	call get_atoms(4, atoms_id, ierr)
	if(ierr.ne.0) return 1
	curr=f_tor(atoms_id)
	z1=curr
	call read_ar(1,z1,*900,*50)

50	if(verbose.ge.2 ) then 
	  write(6,1001) curr,z1
1001	  format(' rtn-I2> the torsional angle will be set from ',
	1 f8.3,' to ',f8.3,' (degrees).')
	  endif

	phi= z1-curr
c	s=0.

	ia=atoms_id(2)
	ja=atoms_id(3)

	p1(1)=x(ia)
	p1(2)=y(ia)
	p1(3)=z(ia)
	p2(1)=x(ja)
	p2(2)=y(ja)
	p2(3)=z(ja)

	call zrot(phi,trn)
	cz(1)=p2(1)-p1(1)
	cz(2)=p2(2)-p1(2)
	cz(3)=p2(3)-p1(3)
	r=cz(1)*cz(1)+cz(2)*cz(2)+cz(3)*cz(3)
	if(r.le.1e-4) return
	  r=1./sqrt(r)
	  cz(1)=cz(1)*r
	  cz(2)=cz(2)*r
	  cz(3)=cz(3)*r
	  if(cz(1).eq.0.) then
	    cy(1)=1.
	    cy(2)=0.
	    cy(3)=0.
	  else
	    r=1./sqrt(cz(1)*cz(1)+cz(2)*cz(2))
	    cy(1)=cz(2)*r
	    cy(2)=-cz(1)*r
	    cy(3)=0.
	    end if
	  call uxveqw(cy,cz,cx)
	  do 100 i=1,3
	    trn0(1,i)=cx(i)
	    trn0(2,i)=cy(i)
	    trn0(3,i)=cz(i)
	    trn1(i,1)=cx(i)
	    trn1(i,2)=cy(i)
100	    trn1(i,3)=cz(i)
	call axbeqc(trn1,trn,arr)
	call axbeqc(arr,trn0,trn)
	do 200 i=1,3
200	trn(i,4)=p1(i)-trn(i,1)*p1(1)-trn(i,2)*p1(2)-trn(i,3)*p1(3)
c	1 +trn1(i,3)*s
	return
900	return 1
	end

	subroutine trnxyz(trn1,*)
chk	=================
	include 'edp_main.inc'
!	use edp_main
	real trn1(3,4)
	integer num_atom

	num_atom=0

	do i=1,n_atom
	  if(lf(i)) then
	    num_atom= num_atom +1 
	    xi=x(i)
	    yi=y(i)
	    zi=z(i)
	    x(i)= trn1(1,4) +trn1(1,1)*xi +trn1(1,2)*yi +trn1(1,3)*zi
	    y(i)= trn1(2,4) +trn1(2,1)*xi +trn1(2,2)*yi +trn1(2,3)*zi
	    z(i)= trn1(3,4) +trn1(3,1)*xi +trn1(3,2)*yi +trn1(3,3)*zi
	    write(text(i)(31:54),1051,err=900) x(i),y(i),z(i)
1051	    format(3f8.3)
	    endif
	  enddo

	if(num_atom.le.0 ) then
	  write(6,1099) 
	1 ' rtn-W> no on atoms is found.'
	1,' only is the matrix calculated.'
	  status= 0
	else 
	  status= num_atom
	endif
	return

900	write(6,1099) 
	1 ' rtn-W> ERROR in writing in the pdb format.'
	1,' for safety, enter reset.'
1099	format(a,a)
	return 1
	end

copyright by X. Cai Zhang

	subroutine vector(e,b,f,*)
chk	=================
	dimension a(3,3),f(3), b(3,3)
	data eps/1.e-4/

	do 10 i=1,3
	a(1,i)=b(1,i)
	a(2,i)=b(2,i)
	a(3,i)=b(3,i)
10	a(i,i)=b(i,i)-e

	do 20 i=1,2
	do 20 j=i+1,3
	  d0=a(i,2)*a(j,3)-a(i,3)*a(j,2)
	  d1=a(i,3)*a(j,1)-a(i,1)*a(j,3)
	  d2=a(i,1)*a(j,2)-a(i,2)*a(j,1)
	  fs=sqrt(d0*d0+d1*d1+d2*d2)
	  if(fs.gt.eps) goto 40
20	  enddo
	return 1

40	f(1)=d0/fs
	f(2)=d1/fs
	f(3)=d2/fs
	end

	subroutine axis03 (r, f,tt1,tt2,tt3,screw_length,i_axis, *)
chk	=================
c	transform: 	x'=rx+t
c	axis:		x+(t.f)f=rx+t => (i-r)x=t-(t.f)f
	real r(3,4),f(3)

	screw_length= r(1,4)*f(1) +r(2,4)*f(2) +r(3,4)*f(3)

	fmax=0.
	if(i_axis.eq.0) then
	  do i=1,3
	    if(fmax.lt.abs(f(i))) then
	      fmax=abs(f(i))
	      j=i
	      endif
	    enddo
	else
	  j=i_axis
	  endif

	if(j.eq.1) then
	  a11=1.-r(2,2)
	  a12=  -r(2,3)
	  a21=  -r(3,2)
	  a22=1.-r(3,3)
	  b1=r(2,4)-screw_length*f(2)
	  b2=r(3,4)-screw_length*f(3)
	  d=a11*a22-a21*a12
	  if(d.eq.0.) return 1
	  tt1=0.
	  tt2=(b1*a22-b2*a12)/d
	  tt3=(b2*a11-b1*a21)/d
	else if(j.eq.2) then
	  a11=1.-r(3,3)
	  a12=  -r(3,1)
	  a21=  -r(1,3)
	  a22=1.-r(1,1)
	  b1=r(3,4)-screw_length*f(3)
	  b2=r(1,4)-screw_length*f(1)
	  d=a11*a22-a21*a12
	  if(d.eq.0.) return 1
	  tt2=0.
	  tt3=(b1*a22-b2*a12)/d
	  tt1=(b2*a11-b1*a21)/d
	else if(j.eq.3) then
	  a11=1.-r(1,1)
	  a12=  -r(1,2)
	  a21=  -r(2,1)
	  a22=1.-r(2,2)
	  b1=r(1,4)-screw_length*f(1)
	  b2=r(2,4)-screw_length*f(2)
	  d=a11*a22-a21*a12
	  if(d.eq.0.) return 1
	  tt3=0.
	  tt1=(b1*a22-b2*a12)/d
	  tt2=(b2*a11-b1*a21)/d
	  endif
	end
chk***	end of rtnout.for

copyright by X. Cai Zhang
