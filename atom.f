copyright by X. Cai Zhang

	subroutine pick(*)
chk	===============
	include 'edp_main.inc'
!	use edp_main
	include 'edp_dat.inc'
!	use edp_dat

	logical  except, lgo
	character*4  atomi, atomj, key(max_num_zones)

	logical  nword
	character*(max_num_chars) txt, txt1
	common /cmm_txt/n_len,txt,ib,ie

	parameter (max_atom_list=max_res)
!	integer, parameter :: max_atom_list=max_res
	character*4 atom_l(max_atom_list)
	character*6 sexcept
	integer  katom(max_atom_list), nk(max_num_zones), mk(max_num_zones)
	data sexcept/'except'/
	save
	
	n_of_syn=4					!000515
	syntax(1)='syntax:'
	syntax(2)='1) atom' 
	syntax(3)='2) atom atom_1.s [atom_2.s ... ]' 
	syntax(4)='3) atom except atom_1.s [atom_2.s ... ]'

	if(atom_untouch) then
	  k_atom=0
	  atom_untouch=.false.
	  endif

	ie0=ie
	if(k_atom.le.0) then
	  lgo=.false.
          n_len1=n_len
	  txt1=txt
	  goto 800
	  endif

400	ie=ie0
	if(nword(n_len,txt,ib,ie)) goto 850	! list atom information
        ie=ie0

c	igr:	id of the basckit group, 0-max_groups
	call dfgroup(igr,*901)

	jj=0
	except=.false.
601	call find(k_atom, atom_l, kk)
	if(kk) 6011, 6014, 6013
6011	if(jj.eq.0 .and. sexcept(1:min(len(sexcept),ie-ib+1)).eq. 
	1 txt(ib:ie) ) then
	  except=.true.
	  goto 601
	  end if

	if( index( txt(ib:ie), wildcard) .gt. 0 ) then
	  if( txt(ib:ie) .eq. '*') then
	    j1=n_groupa(igr)
	    do j=1,j1
	      lf(igroupa(j,igr))=incl
	      enddo
	      return
	    endif
	else if( index( txt(ib:ie), wildcard1) .le. 0 ) then
	  if(verbose.ge.2 ) then
	    write(6,'(a)') 
	1' atom-W2> atom ['//txt(ib:ie)//'] does not exist'
	    endif
	  goto 601
	  endif

6013	jj=jj+1
	if(jj.gt.max_num_zones) goto 700
	key(jj)=txt(ib:ie)
	nk(jj) =max(0,index(key(jj), wildcard)-1)
	mk(jj) =index( txt(ib:ie), wildcard1)
	goto 601

6014	if(jj.lt.1) goto 900

	if(except) then
	  assign 602 to igo
	else
	  assign 6012 to igo
	  end if
 
	j1=n_groupa(igr)
	do 602 j=1,j1
	  i=igroupa(j,igr)
	  if(lf(i).eqv.incl) goto 602
	  atomi=atom(i)
	  do ii=1,jj
	    atomj=key(ii)
	    if(mk(ii).gt.0) then
	      do kk=mk(ii),4
	        if(atomj(kk:kk).eq.'%'.and.atomi(kk:kk).ne.' ') 
	1         atomj(kk:kk)=atomi(kk:kk)
	        enddo
	      endif
	    njj=nk(ii)
	    if(njj.eq.0) then
	      if( atomi.eq.atomj ) goto igo
	    else 
	      if( atomi(:njj).eq. atomj(:njj) ) goto igo
	      endif
	    enddo
	  if(.not.except) goto 602
6012	  lf(i)=incl
602	  continue
	return
c---
700	errmsg=' errmsg: too many atom names in the input list'
	return 1

900	errmsg=' errmsg: unrecoganizible atom name(s)'
	return 1

901	errmsg=' errmsg: wrong group/zone information'
	return 1

	entry beforea
c	=============
	lgo=.true.
	txt1=txt
	ie1=ie
	n_len1=n_len
800	atom_l(1)=atom(1)
	k_atom=1
	n_len=4
	do i=2,n_atom
	  atomi=atom(i)
	  j=ltrim(atomi)
	  do jj=1,j
	    if(atomi(jj:jj).eq.' ' ) atomi(jj:jj)='_'
	    enddo
	  atom(i)=atomi
	  txt=atomi
	  ie=0
	  call find(k_atom,atom_l,j)
	  if(j.le.0) then
	    k_atom=k_atom+1
	    atom_l(k_atom)=atom(i)
	    if(k_atom.ge.max_atom_list) then
	      write(6,*) ' warning: too many atom types.'
	      if(verbose.ge.3 ) write(6,*)
	1'%EdPDB-I3- increase max_atom_list in edp_dim.inc (currently'
	1,max_atom_list,')'
 	      goto 810
	      endif
	    end if
	  end do
810	txt=txt1
	n_len=n_len1
	if(.not.lgo) goto 400
	ie=ie1
	return

850	do j=1,k_atom
	katom(j)=0
	end do

	do 750 i=1,n_atom
	  if(lf(i)) then
	    atomi=atom(i)
	    do j=1,k_atom
	      if(atomi.eq.atom_l(j)) then
	        katom(j)=katom(j)+1
		goto 750
		end if
	      end do
	    write(6,*) '%EdPDB-F- unexpected error in listatom'
	    write(6,*)
	    call exit( 4)
	  endif
750	continue

	write(6,1004) (atom_l(i),min(katom(i),999),i=1,k_atom)
1004	format(' atom> ',(t10,6(a5,'[',i3,']')))
	if(k_atom.eq.max_atom_list) write(6,1002) max_atom_list
1002	format(' atom-W> only',i4,' atom_names have been listed out')
	end

	subroutine dfgroup(igr,*)
chk	==================
	include 'edp_main.inc'
!	use edp_main
	include 'edp_dat.inc'
!	use edp_dat

c	character*1 string1

	logical  nword
	character*(max_num_chars) txt
	common /cmm_txt/n_len,txt,ib,ie

	integer i1(max_num_zones), i2(max_num_zones)

	if (i_pipe.ge.0) then
	  igr=0
	  return
	  endif

	ie0=ie
	i=index(txt(:n_len),'from') 
	if( i.gt.1 ) then
	  ie=i+3
	else 
	  goto 200
	  endif
	ie1=ie
	in=i-2

	igr=match_l(max_gr,cgroup)
	if(igr.gt.0) then
	    if(.not.nword(n_len,txt,ib,ie)) return 1
	    if( n_groupa(igr).le.0) return 1
	    n_len=in
	    ie=ie0
	    ib=ie
	    return
	else if(igr.eq.0) then
	    goto 200
	    endif

	if(txt(ib:ib).eq.'{') then
	  call get_zone(ie0,in,*900)
	  n_len=in
	  ie=ie0
	  ib=ie
	  igr=0
	  return
	  endif

	ie=ie1
	call region(k,i1,i2,*900)
	n_len=in
	ie=ie0
	ib=ie
	if(k.eq.0) then
	    goto 200
	else if(k.gt.max_num_zones) then
	    write(6,*) 'dfgroupr-W> too many zones.'
	    if(verbose.ge.3 ) write(6,*)
	1'%EdPDB-I3- increase max_num_zones in edp_dim.inc (currently'
	1 ,max_num_zones,')'
	    return 1
	    endif

	n_g=0
	do j=1,k
	do i=i1(j),i2(j)
	    n_g=n_g+1
	    igroupa(n_g,0)=i
	    enddo 
	    enddo
	n_groupa(0)=n_g
	igr=0
	return	  

200	if(index(txt(:n_len),'from').gt.1) then
	  errmsg=' errmsg: from where ?'
	  return 1
	  endif
c	if(n_groupa(0).ne. n_atom) then
c	  n_groupa(0)=n_atom
c	  do i=1,n_atom
c	    igroupa(i,0)=i
c	    enddo
c	  endif
	call init_group0()
	igr=0
	return
900	continue
	return 1
	end

	subroutine get_zone(ie_s,in_s,*)
chk	===================
	include 'edp_main.inc'
!	use edp_main

	logical lf_s(max_atom)
	integer ie0(max_num_zones), in0(max_num_zones), ie1(max_num_zones)
	character*(max_num_chars) txt_s(max_num_zones)

        character*(max_num_chars) txt
        common /cmm_txt/n_len,txt,ib,ie

	common /main_para/ tolower, jou, echo_L, junk_i0, incl0, incl1
	1 ,jjou ,nest
	1 ,jnk1(m_main_para-8)

	logical junk_i0, echo_L, tolower, incl0, incl1
 
	n_of_syn=2					!000515
	syntax(1)='syntax:' 
	syntax(2)='zone [res_id1.s [res_id2.s ... ]] '

	if(nest.lt.0) return 1
	nest=nest+1
	if(nest.gt.max_num_zones) then
	  errmsg= 'errmsg: too many nested commands'
	  goto 900
	  endif

	nc=0
	do i=ib,n_len
	  if(txt(i:i).eq.'{') then
	    nc=nc+1
	  else if(txt(i:i).eq.'}') then
	    nc=nc-1
	    endif
	  if(nc.eq.0) then
	    ie=i
	    goto 200
	    endif
	  enddo
	errmsg=' errmsg: wrong nested commands'
	goto 900

200	txt_s(nest)=txt(1:in_s)
	in0(nest)=in_s
	ie0(nest)=ie_s
	ie1(nest)=ib
	if(nest.eq.1) then
	  do i=1,n_atom
	    lf_s(i)=lf(i)
	    lf(i)=.false.
	    enddo
	  endif

	txt=txt(ib+1:ie-1)
	n_len=ie-ib-1
	if(i_pipe.ge.0) then
	  errmsg=' errmsg: PIPE system and FROM_subselection'//
	1 ' are mutually exclusive.'
	  return 1
	  endif

	call get_command(*901,*900)

	n_g=0
	do ii=1,n_atom
	    if(lf(ii)) then 
	      n_g=n_g+1
	      igroupa(n_g,0)=ii
	      lf(ii)=.false.
	      endif
	  enddo
	n_groupa(0)=n_g

	txt=txt_s(nest)
	in_s=in0(nest)
	ie_s=ie0(nest)
	if(nest.eq.1) then
	  do i=1,n_atom
	    lf(i)=lf_s(i)
	    enddo
	  endif
	nest=nest-1
	if(nest.le.0) incl=incl1
	return

901	errmsg=' errmsg: invalid command'
900	do i=1,nest
	  ie=ie+ie1(i)
	  enddo
	ib=ie
	do i=1,n_atom
	  lf(i)=lf_s(i)
	  enddo
	nest=-1

	return 1
	end
c***	the end of atom.for

copyright by X. Cai Zhang
