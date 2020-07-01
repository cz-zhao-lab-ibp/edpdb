
	SUBROUTINE NEEDLEMAN_1d(RATIO,JTRY) 
chk	====================
	include 'edp_dim.inc'
!	use edp_dim
	include 'edp_dat.inc'
!	use edp_dat

	PARAMETER (MAX_T=  max_rt)
	parameter (max_n=  max_l)
	parameter (max_2n=max_n*2)
!      integer, PARAMETER :: MAX_T=  max_rt 
!      integer, parameter :: max_n=  max_L
!      integer, parameter :: max_2n=max_n*2

	integer ia0 ,ib0 ,ic0 
	common /cmm_1d/ n1, n2, ad (max_n) ,bd(max_n)
	1 ,ia0(max_n) ,ib0(max_n) ,ic0(max_n)

	character*1 ss, pp
	common /cmm_tbl_1d/ kp, ss(max_t), sc(max_t,max_t)
	1 ,mp ,pp(max_t), dd(max_t)

	common /cmm_ctl_1d/ nr, iseed, pe, n_end0, sl0, sl1, sl2
	1 ,igr0, igr1, igr2, blk_file

	DIMENSION SM(max_n,max_n), TK_A(max_n)
	INTEGER*2 ic1(MAX_2N), ic2(MAX_2N)
	INTEGER*2 IGO(max_n,max_n),  KGO_A(max_n)

	if(verbose.ge.6 ) write(6,1069) 	!000504
1069	format(' needleman-I6> reference:'
	1/'  Needleman SB, Wunsch CD.'
	1/'  A general method applicable to the search for similarities '
	1/'    in the amino acid sequence of two proteins.'
	1/'  J Mol Biol. 1970 Mar;48(3):443-53.')

	n_end= min(n_end0,n1-1,n2-1)
	do i=1,max_2n
	  ic1(i)=max_t
	  ic2(i)=max_t
	  enddo
C
	IBJ=ib0(1)
	DO I= 1, N1
	SM(I,1)= SC(ia0(I),IBJ)
	IGO(I,1)=0
	END DO

	IAI=ia0(1)
	DO J= 2, N2
	SM(1,J)= SC(IAI, ib0(J))
	IGO(1,J)=0
	END DO
C
	IAI=ia0(2)
	IBJ=ib0(2)
	SM(2,2)=SM(1,1)+SC(IAI,IBJ)
	IGO(2,2)=0
C
	BDJ=BD(1)
	PEBDJ=PE*BDJ
	TL=-1024.
	DO I=3,N1
		JGO=0
		SMAX=SM(I-1,1)

		STMP=SM(I-2,1)+BDJ
		IF(SMAX.LT.STMP) THEN
		  SMAX=STMP
		  JGO=2-I
		END IF
		
		IF(SMAX.LT.TL) THEN
		  JGO=LGO
		  SMAX=TL
		  TL=TL+PEBDJ
		ELSE IF(TL.LT.STMP) THEN
		  TL=STMP+PEBDJ
		  LGO=2-I
		ELSE
		  TL=TL+PEBDJ
		END IF
	
		SM(I,2)=SMAX+SC(ia0(I),IBJ)
		IGO(I,2)=JGO
	END DO

	BDJ=AD(1)
	PEBDJ=PE*BDJ
	TL=-1024.
	DO I=3,N2
		JGO=0
		SMAX=SM(1,I-1)

		STMP=SM(1,I-2)+BDJ
		IF(SMAX.LT.STMP) THEN
		  SMAX=STMP
		  JGO= I-2
		END IF
		
		IF(SMAX.LT.TL) THEN
		  JGO=LGO
		  SMAX=TL
		  TL=TL+PEBDJ
		ELSE IF(TL.LT.STMP) THEN
		  TL=STMP+PEBDJ
		  LGO=I-2
		ELSE
		  TL=TL+PEBDJ
		END IF
	
		SM(2,I)=SMAX+SC(IAI,ib0(I))
		IGO(2,I)=JGO
	END DO
C
	DO I=3,N1
		TK_A(I)=-1024.
	END DO
	DO J=3,N2
		J1=J-1
		J2=J-2
		J3=J-3
		BDJ=BD(J)
		PEBDJ=PE*BDJ
		IBJ=ib0(J)
		
		TL=-1024.
		DO I=3,N1
			I1=I-1
			I2=I-2
			I3=I-3
			ADI=AD(I)
			PEADI=PE*ADI

			JGO=0
			SMAX=SM(I1, J1)

			STMP=SM(I2,J1)+BDJ
			IF(SMAX.LT.STMP) THEN
			  SMAX=STMP
			  JGO=-I2
			END IF
		
			IF(SMAX.LT.TL) THEN
			  JGO=LGO
			  SMAX=TL
			  TL=TL+PEBDJ
			ELSE IF(TL.LT.STMP) THEN
			  TL=STMP+PEBDJ
			  LGO=-I2
			ELSE
			  TL=TL+PEBDJ
			END IF

			STMP=SM(I1,J2)+ADI
			IF(SMAX.LT.STMP) THEN
			  SMAX=STMP
			  JGO=J2		! why not -j2
			END IF
		
			TK=TK_A(I)
			KGO=KGO_A(I)
			IF(SMAX.LT.TK) THEN
			  JGO=KGO
			  SMAX=TK
			  TK=TK+PEADI
			ELSE IF(TK.LT.STMP) THEN
			  TK=STMP+PEADI
			  KGO=J2		! why not -j2
			ELSE
			  TK=TK+PEADI
			END IF
			TK_A(I)=TK
			KGO_A(I)=KGO

			SM(I,J)=SMAX+SC(ia0(I),IBJ)
			IGO(I,J)=JGO
		END DO
	END DO

C
	SMAX=SM(N1,N2)
	IMAX=N1
	JMAX=N2
	I1=N1- n_end
	DO I=I1,N1
		IF(SMAX.LT.SM(I,N2)) THEN
			IMAX=I
			SMAX=SM(I,N2)
		END IF
	END DO

	J1=N2- n_end
 	DO J= J1, N2
		IF(SMAX.LT.SM(N1,J)) THEN
			JMAX=J
			SMAX=SM(N1,J)
		END IF
	END DO

	ratio=SMAX/MIN(N1,N2)

	IF(JTRY.GT.0) RETURN
	if(verbose.ge.2 ) WRITE(6,1006) SMAX, RATIO
1006	format(' match1d-I2>   Score= ', F8.3,', Ratio= ',F8.3)

c	IF(JTRY.GT.0) RETURN
	K=0
	IF(JMAX.LT.N2) THEN
		IMAX=N1
		DO J=N2, JMAX+1, -1
			K=K+1
			IC2(K)=ib0(J)
		END DO
	ELSE 
		DO I=N1, IMAX+1, -1
			K=K+1
			IC1(K)=ia0(I)
		END DO
	END IF

C	Search maximun score

	kgap1=0		!970211
	kgap2=0

	DO WHILE(IMAX.GE.1.AND.JMAX.GE.1)
		K=K+1
		IC1(K)=ia0(IMAX)
		IC2(K)=ib0(JMAX)
		JGO=IGO(IMAX,JMAX)
		IF(JGO) 100, 200, 300

100	KGAP2=KGAP2+1
	IS=IMAX-1
	IMAX=1-JGO
	DO I=IS,IMAX,-1
	K=K+1
	IC1(K)=ia0(I)
	END DO
	GOTO 200

300	KGAP1=KGAP1+1
	JS=JMAX-1
	JMAX=JGO+1
	DO J= JS,JMAX,-1
	K=K+1
	IC2(K)=ib0(J)
	END DO

200	IMAX=IMAX-1
	JMAX=JMAX-1
	END DO

	IF(IMAX.GT.1) THEN
		DO I=IMAX,1,-1
			K=K+1
			IC1(K)=ia0(I)
		END DO
	ELSE IF(JMAX.GT.1) THEN
		DO J=JMAX,1,-1
			K=K+1
			IC2(K)=ib0(J)
		END DO
	END IF

500	if(verbose.ge.2 ) WRITE(6,1005) KGAP1,KGAP2
1005	FORMAT(' match1d-I2>   # of gaps = ',i3,', ',i3)
	CALL output_ndlmn_1d(IC1,IC2,K,'match1d')
	END


        subroutine score_table_1d(file,*)
chk     ======================
	include 'edp_main.inc'
!	use edp_main
	include 'edp_dat.inc'
!	use edp_dat

	logical fo,fo2
	data	fo,fo2/.false.,.false./
	PARAMETER 	(max_t=  max_rt )
	parameter 	(max_n=  max_l)
	parameter 	(max_2n=max_n*2)
!	integer, PARAMETER :: max_t=  max_rt 
!	integer, parameter :: max_n=  max_L
!	integer, parameter :: max_2n=max_n*2

	character*(*) file
 	character*1 ss, pp
	common /cmm_tbl_1d/ kp, ss(max_t), sc(max_t,max_t)
	1 ,mp ,pp(max_t), dd(max_t)

	character*(120) fmtstmt, curr_line

        character*(108) curr_file
        logical  open_file
        external open_file
	save

	if(fo) return

        curr_file=file
80      if(.not.open_file(79,curr_file,'old','.txt')) then
          if(curr_file .ne. file) goto 81
          curr_file=edp_data(:ltrim(edp_data))//curr_file
          goto 80
          endif

	do i=1,max_t
	  fmtstmt='(3x,a1,<i>f4.1)'
	  write(fmtstmt(8:10),'(i3)') i
	  read(79,1001,end=100) curr_line
	  if(curr_line(1:1).ne.'!') then 
	    read(curr_line,fmtstmt,err=110) ss(i), (sc(i,j),j=1,i)
	    if(ss(i).eq.' ') goto 100
	  endif
          enddo
1001    format(a120)

100     kp=i-1
        do i=1,kp
        do j=1,i
          sc(j,i)=sc(i,j)       !if(sc(j,i).eq.0.) 
          enddo
          enddo
	close(79)
	fo=.true. 
        return

110	if(verbose.ge.9 ) then
	  write(6,'(a,i)') 
	1 ' score_table-I9> # of lines read in =',i
	  write(6,'(a)') 
	1 ' score_table-I9> format= ['//fmtstmt(:ltrim(fmtstmt))//']' 
	  do ii=1,i
	    write(6,fmtstmt) ss(ii), (sc(ii,j),j=1,i)
	  enddo
	endif 
	errmsg= 
	1' errmsg: error during read score_table'
	return 1

81      errmsg=
	1' errmsg: failure in opening the score file ['//file//']'
        return 1
        
        entry penalty_table_1d(file)
chk     ======================
	if(fo2) return
	fo2=.true.
        curr_file=file
90      if(.not.open_file(79,curr_file,'old','.txt')) then
          if(curr_file .ne. file) goto 91
          curr_file=edp_data(:ltrim(edp_data))//curr_file
          goto 90
          endif

	do i=1,max_t
          read(79,1002,end=105,err=104) pp(i), dd(i)
	  if(pp(i).eq.' ') goto 105
          enddo
1002    format(3x,a1,f8.3)
105     mp=i-1
	close(79)  
        return

104	if(verbose.ge.9 ) then
	  write(6,'(a,i)') 
	1 'penalty_table-I9> # of lines read in =',i
	  do ii=1,i
	    write(6,1002) pp(ii), dd(ii)
	  enddo
	endif 
	write(6,*) 'match1d-W> error during read penalty_table'

91      do i=1,max_t
          pp(i)=ss(i)
          dd(i)=-3.0
          enddo
        end

	subroutine match_ab
chk	===================
	include 'edp_dim.inc'
!	use edp_dim

	PARAMETER (MAX_T=  max_rt) 
	parameter (max_n=  max_l)
	parameter (max_2n=max_n*2)

!      integer, PARAMETER :: MAX_T=  max_rt 
!      integer, parameter :: max_n=  max_L
!      integer, parameter :: max_2n=max_n*2

	integer ia0, ib0, ic0
	common /cmm_1d/ num_aa(2), ad (max_n) ,bd(max_n)
	1 ,ia0(max_n) ,ib0(max_n), ic0(max_n)

	character*1 ss, pp
	common /cmm_tbl_1d/ kp, ss(max_t), sc(max_t,max_t)
	1 ,mp ,pp(max_t), dd(max_t)

	INTEGER*2 iaa(max_n)

	LOGICAL LF(max_n)

	common /cmm_ctl_1d/ nr, iseed, pe, n_end, sl0, sl1, sl2
	1 ,igr0, igr1, igr2, blk_file

	n1=num_aa(1)
	n2=num_aa(2)
	AV=0.
	SGM=0.
	DO J=0,NR
 	CALL NEEDLEMAN_1d(RATIO,J)
	IF(J.GT.0) THEN
	AV=AV+RATIO
	SGM=SGM+RATIO*RATIO
	ELSE
	RATIO_1ST=RATIO
	END IF
	IF(J.EQ.NR) GOTO 800

	DO I=1,N1
	LF(I)=.FALSE.
	iaa(I)=ia0(I)
	END DO

	DO I=1, N1
	N=NINT(RAN(ISEED)*N1)+1
210	IF(N.GT.N1) N=N-N1
	IF(LF(N)) THEN
		N=N+1
		GOTO 210
	END IF
	LF(N)=.TRUE.
	ia0(I)=iaa(N)
	END DO
	END DO

800	IF(NR.GT.2) THEN
	RN=1./FLOAT(NR)
	AV=AV*RN
	SGM=SQRT(MAX(0.,SGM*RN-AV*AV))
	WRITE(6,1007) NR,AV, SGM	
	IF(SGM.GT.0.) WRITE(6,1008) (RATIO_1ST-AV)/SGM
1007	FORMAT(
	1/' match1d-I> # of RANDOM TRY= ',i8
	1/' match1d-I> AV_RATIO       = ',F8.3
	1/' match1d-I> SIGMA          = ',F8.3)
1008	FORMAT( 
	1 ' match1d-I> SIGNIFICANCE   = ',f8.3,' sgm')
	END IF
	END 
	
	subroutine match1d(*)
chk	==================
c       unit 79: penalty_table file,    format(3x,a1,f8.3)
c       	 score_table file,      format(3x,a1,<max_t>f4.1)

	include 'edp_main.inc'
!        use edp_main
	include 'edp_dat.inc'
!        use edp_dat

	PARAMETER  (MAX_T=  max_rt )
	parameter  (max_n=  max_l)
	parameter  (max_2n=max_n*2)

!      integer, PARAMETER :: MAX_T=  max_rt 
!      integer, parameter :: max_n=  max_L
!      integer, parameter :: max_2n=max_n*2

        integer ia0, ib0, ic0
	common /cmm_1d/ num_aa(2), ad (max_n,2)
	1 ,ia0(max_n,2), ic0(max_n)

        character*1 ss, pp
        common /cmm_tbl_1d/ kp, ss(max_t), sc(max_t,max_t)
	1 ,mp ,pp(max_t), dd(max_t)

        common /cmm_ctl_1d/ nr, iseed, pe, n_end, sl0, sl1, sl2
	1 ,igr0, igr1, igr2, blk_file

	character*1 s1(max_n), si

	character*5 j60(2,20)
	integer*2 ic1(max_2n), ic2(max_2n)
	character*1 c1(max_2n), c2(max_2n), c3(max_2n)
c	character*1 cj, ck

	logical nword 
	character*(max_num_chars) txt
	common /cmm_txt/ n_len,txt,ib,ie

	logical open_file, blk_file
	external open_file


	n_of_syn=4					!000515
	syntax(1)='syntax:' 
	syntax(2)='load group2.s | match1d group1.s sub_group1.s '
	syntax(3)=
	1'  [score_threshold.r [list_threshold_1.r, [list_threshold_2.r,'
	syntax(4)=
	1'  [extension_penalty.r, [num_of_random_trial.i, [seed.i]]]]]]' 
	
	call dfgroup(igr0,*901)
	n_group0=n_groupa(igr0)
	if(n_group0.gt.max_n) then
	  write(errmsg,'(a,i5,a)') 
	1 ' errmsg: too many atoms to select.(max_n=',max_n,')'
	  return 1
	  endif

	igr1= match_l( max_gr, cgroup)
	if( igr1 .le. 0) then
	  errmsg=' errmsg: a group name is needed'
	  return 1
	  end if
	n_group1= n_groupa(igr1)
	if( n_group1.le.0) then
	  write(6,*) 
	1 'match1d-W> UNDONE: define group ['//cgroup(igr1)//'] first.'
	  return
	else if(n_group1.gt.max_n) then
	  write(errmsg,'(a,i5,a)') 
	1 ' errmsg: too many atoms in group ['//cgroup(igr1)//
     1 '] (max_n=',max_n,')'
	  return 1
	  end if

	call find_a_group(ng2)
	if( ng2 .le. 0) return 1

	igr2=ng2

        NR=0
        ISEED=1
        PE=0.1
        SL0=1.5
        SL1=1.0
        SL2=0.5
        N_END=1000
        io=0

5005    call read_ar(1,sl0      ,*901,*5000)   !  default      1.5
        call read_ar(1,sl1      ,*901,*5000)   !               1.0
	call read_ar(1,sl2      ,*901,*5000)   !               0.5
        call read_ar(1,pe       ,*901,*5000)   !               0.1
	call read_ai(1,nr       ,*901,*5000)   !	       0
        call read_ai(1,iseed    ,*901,*5000)   !               1
        call read_ai(1,n_end    ,*901,*5000)   !               1000
        
5000	nr=max(0,min(99,nr))
        n_end=max(0,n_end)

	call   score_table_1d('score_seq.txt',*901)
	call penalty_table_1d('penalty_table.txt')

	ks=1
	call get_sequence(igr0,s1)

	num_aa(1)=n_group0
200	do i=1,num_aa(ks)
	  si=s1(i)
	  do j=1,kp
	    if(si.eq.ss(j))goto 220
	    enddo
	  errmsg=' errmsg: a symbol unmatched with score_seq: '//si
	  return 1
220	  ia0(i,ks)=j
	  do j=1,mp
	    if(si.eq.pp(j)) goto 240
	    enddo
	  errmsg=' errmsg: a symbol unmatched with penalty_table: '//si
	  return 1
240	  ad(i,ks)=dd(j)
	  enddo
	if(ks.eq.1) then	
	  ks=2
	  call get_sequence(igr1,s1)
	  num_aa(2)=n_group1
	  goto 200
	  endif

	if(verbose.ge.4 ) then
	  write(6,1016) sl0, sl1, sl2, pe, nr,iseed
	  WRITE(6,1006) num_aa
	  endif
1016	format(
	1 ' match1d-I4> threshold for pair-selection=',f6.1
	1/' match1d-I4> threshold for listing as (:)=',f6.1
	1/' match1d-I4> threshold for listing as (.)=',f6.1
	1/' match1d-I4>            extension penalty=',f6.1
	1/' match1d-I4>           # of random trials=',I6
	1/' match1d-I4>     seed for random trial(s)=',I6)
1006	format(
	1 ' match1d-I4>     Length= ',i3, ',',i4)

	if(num_aa(1).gt.3 .and. num_aa(2).gt.3) then
	  blk_file=.false.			!000413
	  if(verbose.ge.6 ) then 
	    blk_file=open_file(37, 'seq_.txt', 'new', '')
	    if(blk_file) write(37,1031)
	  endif
	  call match_ab
	endif
1031	format('>Seq-1'/'>Seq_2'/'*')
	return
901	return 1
	end

	subroutine output_ndlmn_1d(IC1,IC2,K,switch)
chk	=====================
	include 'edp_main.inc'
!	use edp_main
	include 'edp_dat.inc'
!	use edp_dat

	PARAMETER (MAX_T=  max_rt )
	PARAMETER (max_n=  max_l)
	PARAMETER (max_2n=max_n *2)

!      integer, PARAMETER :: MAX_T=  max_rt 
!      integer, parameter :: max_n=  max_L
!      integer, parameter :: max_2n=max_n*2

        integer ia0, ib0, ic0
	common /cmm_1d/ num_aa(2), ad (max_n,2)
	1 ,ia0(max_n,2), ic0(max_n)

        character*1 ss, pp
        common /cmm_tbl_1d/ kp, ss(max_t), sc(max_t,max_t)
	1 ,mp ,pp(max_t), dd(max_t)

        common /cmm_ctl_1d/ nr, iseed, pe, n_end, sl0, sl1, sl2
	1 ,igr0, igr1, igr2, blk_file

	character*5 j60(2,20)
	integer*2 ic1(max_2n), ic2(max_2n)
	character*1 c1(max_2n), c2(max_2n), c3(max_2n)

	logical blk_file, open_file
	character*(*) switch

	n_g=0
	i1=0
	i2=0
	do k0=k,1,-1
	  j1 = ic1(k0)
	  j2 = ic2(k0)
	  if(j1.ne.max_t) i1=i1+1
	  if(j2.ne.max_t) i2=i2+1
	  if( j1.ne.max_t .and. j2.ne. max_t) then
	    if( sc(j1,j2) .ge. sl0) then
	      n_g=n_g+1
	      lf(igroupa(i1,igr0)) = incl
	      w(igroupa(i1,igr0)) = sc(j1,j2)
	      igroupa(n_g,igr2) = igroupa(i2,igr1)
	      endif
	    endif
	  enddo

	n_groupa(igr2)= n_g

	if(verbose.ge.4 ) write(6,1018) switch(:ltrim(switch)),n_g
1018	format(' ',a,'-I4>  selected matches= ',I3)

	ss(max_t)='-'

	ii=0		! number of identical residues
	k60=0		! number of output characters
	j=0
	i1=0
	i2=0
	do k0=1,20
	  j60(1,k0)='     '
	  j60(2,k0)='     '
	  enddo

	do k0=k,1,-1
	  j=j+1
	  j1 = ic1(k0)
	  j2 = ic2(k0)
	  c1(j) = ss(j1)
	  c2(j) = ss(j2)
	  if(j1.ne.max_t) i1=i1+1
	  if(j2.ne.max_t) i2=i2+1
	  sc12=sc(j1,j2)

	  if(j1.eq.j2) then
	    c3(j)='|'
	    ii=ii+1
	  else if(sc12.gt.sl1) then
	    c3(j)=':'
	  else if(sc12.gt.sl2) then
	    c3(j)='.'
	  else
	    c3(j)=' '
	    endif

	  if(mod(j,60).eq.1) then
	    k60=k60+1
	    j60(1,k60)=ires(aa_seq(igroupa(1,igr0)))
	    j60(2,k60)=ires(aa_seq(igroupa(1,igr1)))
	    if(i1.gt.0) j60(1,k60)=ires(aa_seq(igroupa(i1,igr0)))
	    if(i2.gt.0) j60(2,k60)=ires(aa_seq(igroupa(i2,igr1)))
	    endif

	if(verbose.gt.6.and..not.blk_file) then
	    blk_file= open_file(37, 'seq_.txt', 'unknown', '')
	    if(blk_file) write(37,1031)
	  endif

	  if(blk_file) write(37,1032) c1(j),c2(j) 
	  enddo
	  if(blk_file) write(37,1033)
1031	format('>Seq-1'/'>Seq_2'/'*') 
1032	format(2a1)
1033	format('*')

	if(verbose.ge.4 ) write(6,1008) switch(:ltrim(switch)),ii
!1008	format(' match1d> identical matches= ',I3)
1008	format(' ',a,'-I4> identical matches= ',I3)

	k60=k/60
	j1=0
	do i=1,k60
	j1=i*60
	j0=j1-59
	write (48,1010)
	write (48,1011) j60(1,i),(c1(j),j=j0,j1)
	write (48,1009)          (c3(j),j=j0,j1)
	write (48,1011) j60(2,i),(c2(j),j=j0,j1)
!	if(i_pipe.ge.0 .and. verbose.ge.4 ) then 	!090528
!	  write (6,1014)
!	  write (6,1011) j60(1,i),(c1(j),j=j0,j1)
!	  write (6,1009)          (c3(j),j=j0,j1)
!	  write (6,1011) j60(2,i),(c2(j),j=j0,j1)
!	  endif
	end do

	if(k-k60*60.gt.0) then
	j0=j1+1
 	j1=k
	write (48,1010)
	write (48,1011) j60(1,i),(c1(j),j=j0,j1)
	write (48,1009)          (c3(j),j=j0,j1)
	write (48,1011) j60(2,i),(c2(j),j=j0,j1)
	write (48,*)
!	if(i_pipe.ge.0 .and. verbose.ge.4 ) then	!090528
!	  write (6,1014) 
!	  write (6,1011) j60(1,i),(c1(j),j=j0,j1)
!	  write (6,1009)          (c3(j),j=j0,j1)
!	  write (6,1011) j60(2,i),(c2(j),j=j0,j1)
!	  write (6,*)
!	  endif
	end if

1009	format(10x, 60a1)
1010	format(10x, 6(9x,','))
1014	format(              19x,',',5(9x,','))
!1014	format(' match1d-I4>',7x,',',5(9x,','))
1011	format(4x,a5,1x,60a1)
	call typelist(48)
	end
  
	subroutine seq2pdb(*)	!991228
chk	==================
	include 'edp_main.inc'
!        use edp_main
	include 'edp_dat.inc'
!        use edp_dat
	include 'edp_file.inc'
!        use edp_file

	character*(108) file_seq
	character*(32) form, astar 
	character*(1) seq2(max_res)

	logical nword 
	character*(max_num_chars) txt
	common /cmm_txt/ n_len,txt,ib,ie
	logical  open_file1
	external open_file1

	!external trim; character trim
	
	n_of_syn=2					!000515
	syntax(1)='syntax:' 
	syntax(2)='seq2pdb [filename.s] [fortran_format.s] '

	if(pdb_out.eq.'?') then
	  write(6,*)
	1'seq2pdb-W> UNDONE: open an output pdb file first'
	  return
	  endif

	file_seq='?'
	if(.not.open_file1(7, file_seq(:ltrim(file_seq))
	1 , 'old', '.seq')) goto 901

!	form='(5(1x,10a1))'
	form='(70a1)'	! for FASTA format, 121004
	if(nword(n_len,txt,ib,ie)) goto 305
	ic=ichar(delimiter)
	if(ichar(txt(ib:ib)).eq. ic ) then
	  do j=ib+1,n_len
	    if(ichar(txt(j:j)).eq.ic) then
	      ib=ib+1
	      ie=j-1
	      txt(j:j)=' '	!000303, it may screw up the history file
	      goto 303
	      endif
	    enddo
	  endif
303	form=txt(ib:ie)
305	if(verbose.ge.2 ) write(6,'(a)') 
	1 ' seq2pdb-I2> input format= '//form//' (default: FASTA)'

	seq_out=file_seq(:ltrim(file_seq))
	ns1=1
	read(7,'(a)',err=903,end=306) astar
	if(astar(1:1).eq.'*') then
	  read(astar(2:),*) ns1		!????
	else if(astar(1:1).eq.'>') then	! fasta 
	  if(verbose.ge.2 ) write(6,*) astar
	else
	  rewind (7)
	  endif
	read(7,form,err=903,end=306) (seq2(ns0),ns0=1,max_res)
306	ns0=ns0-1
	call seq2pdb_output(seq2,ns0, ns1, *904)
!	write(6,*) ns0, ns1		! what is ns1 for? 
	close (7)
	return

901	close (7)
	errmsg=' errmsg: FAILURE in opening the sequence file ['
	1//file_seq(:ltrim(file_seq))//']'
	return 1
903	close (7)
	errmsg=' errmsg: ERROR in reading the sequence file ['
	1//file_seq(:ltrim(file_seq))//']'
	return 1
904	close (7)
	return 1
	end

	subroutine copy(*)
chk	==================
c
c	copy xzy, w and b from a given group to the ON group.
c
	include 'edp_main.inc'
!	use edp_main
	include 'edp_dat.inc'
!	use edp_dat
	
	logical nword 
	character*(max_num_chars) txt
	common /cmm_txt/ n_len,txt,ib,ie

	logical L_x, L_y, L_z, L_w, L_b, L_t

	integer match_l
	external match_l

	n_of_syn=2					!000515
	syntax(1)='syntax:' 
	syntax(2)='copy group_id.s [(x,y,z,w,z), t t1.i t2.i]' 

	status=-2

	igr=match_l( max_gr, cgroup )
	if( igr .le. 0) return 1

	n_group0=0
	do i=1,n_atom
	 if(lf(i)) n_group0=n_group0+1
	 enddo
	n_group1=n_groupa(igr)

	if(verbose.ge.2 ) write(6,1009) n_group0,cgroup(igr),n_group1
1009	format(' copy-I2>  #of on atoms (na)=',i5,
	1 ', #of atoms in group ',a,' (ng)=',i5)

	if(n_group0.ne.n_group1.or.n_group0.eq.0) then
	  errmsg= ' errmsg: na =/= ng or na=0'
	  return 1
	  endif

	L_t=.false.
	L_x=.true.
	L_y=.true.
	L_z=.true.
	L_w=.true.
	L_b=.true.
	if(.not.nword(n_len,txt,ib,ie).and.ib .le.ie) then
	  if(index(txt(ib:ie),'x').le.0) L_x=.false.
	  if(index(txt(ib:ie),'y').le.0) L_y=.false.
	  if(index(txt(ib:ie),'z').le.0) L_z=.false.
	  if(index(txt(ib:ie),'w').le.0) L_w=.false.
	  if(index(txt(ib:ie),'b').le.0) L_b=.false.
	  if(index(txt(ib:ie),'t').gt.0) then 
	    L_t=.true. 
	    i0=1
	    i1=30
            call read_ai(1,i0      ,*900,*100)   ! 
	    i1=i0
            call read_ai(1,i1      ,*900,*100)   ! 
100	    if(i1.lt.i0) then
	      i2=i1
	      i1=i0
	      i0=i2
	    endif
	    if(i0.le.0.or.i1.gt.72) return 1
	  endif
	endif

	jj=0
	if(L_t) then
	do i=1,n_atom
	  if(lf(i)) then
	    jj=jj+1
	    j=igroupa( jj,igr)
	    text(i)(i0:i1)= text(j)(i0:i1)
	    endif
	  enddo
	else 
	 do i=1,n_atom
	  if(lf(i)) then
	    jj=jj+1
	    j=igroupa( jj,igr)
	    if(L_x) x(i)=x(j)
	    if(L_y) y(i)=y(j)
	    if(L_z) z(i)=z(j)
	    if(L_w) w(i)=w(j)
	    if(L_b) b(i)=b(j)
	    write(text(i)(31:54),'(3f8.3)',err=800) x(i),y(i),z(i)
	  endif
	 enddo
	endif
	status=jj
	return
800	errmsg=' errmsg: ABORTED: error during writing the coordinates'
	return 1
900	return 1
	end

	subroutine set_parameter
chk	=========================
	1 (ith, name, name_length, content, cont_length)

	include 'edp_main.inc'
!	use edp_main
	integer np, npp
	character*8  pp
	character*72 pn
	common /ud_param/ npp(max_param), pp(max_param),
	1 np(max_param), pn(max_param)

	character*(*) name, content
	integer   name_length, cont_length

	i=max_reserved_param+ith
	pp(i) =name
	npp(i)=name_length
	pn(i) =content
	np(i) =cont_length
	end

chk	end of copy subroutine

copyright by X. Cai Zhang, 2000
