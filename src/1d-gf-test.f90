program GF_1d_driver

  use modGate
  use modConvParGF, only: convParGFDriver, initModConvParGF &
                         , modConvParGF_initialized &
                         , readGFConvParNML, OUTPUT_SOUND
  use modConstants, only : c_cp, c_alvl
   implicit none

   integer,parameter :: mxp=1   ,myp=1, mtp =1 , nmp = 1 , mynum=1
   integer,parameter :: mgmxp=mxp,mgmyp=myp
   integer,parameter :: maxiens    = 3 !cloud spectral size
   integer,parameter :: &
     ims=1   ,ime=1 ,jms=1  ,jme=1  ,kms=1 &  
    ,its=1  ,ite=1  ,jts=1  ,jte=1  ,kts=1   

   integer :: itime1
   logical :: land
   !- here are the place for data related with the gate soundings
   !- soundings arrays
   integer ::jk, nruns, version, klon_local,klev_local,kte,kme, mzp,mgmzp

   !- this for the namelist gf.inp
   namelist /run/ runname, runlabel, rundata,version, land , klev_sound 

   !- for grads output
   integer :: nrec,nvx,nvar,nvartotal,klevgrads(0:300),int_byte_size,n1,n2,n3
   real    :: real_byte_size
   logical :: init_stat
   
   !C: Carlos Renato testes: 
   integer :: icr, icrf
   real :: dtlt, time, tkmin

   integer, allocatable, dimension(:)       :: flip
   integer, allocatable, dimension(:,:)     :: kpbl,do_this_column
   integer, allocatable, dimension(:,:,:)   :: kstabi4d_tmp, kstabm4d_tmp
   integer, allocatable, dimension(:,:,:)   :: ierr4d_tmp, jmin4d_tmp ,klcl4d_tmp 
   integer, allocatable, dimension(:,:,:)   :: k224d_tmp, kbcon4d_tmp, ktop4d_tmp
   !- 1d
   real   , allocatable, dimension(:)       :: fscav_int
   !- 2d
   real   , allocatable, dimension(:,:)     :: dx2d, lons, lats, aot500, tke_pbl, col_sat, glon, glat
   real   , allocatable, dimension(:,:)     :: stochastic_sig, temp2m, xland, sfc_press,rtgt
   real   , allocatable, dimension(:,:)     :: conprr,lightn_dens, rh_dicy_fct
   real   , allocatable, dimension(:,:)     :: sflux_r, sflux_t, topt, var2d
   real   , allocatable, dimension(:,:)     :: cnv_frc,aa0,aa1,aa2,aa3
   real   , allocatable, dimension(:,:)     :: aa1_bl,aa1_cin,tau_bl,tau_ec
   real   , allocatable, dimension(:,:)     :: wlpool, aa1_adv, aa1_radpbl
   !-- 3d
   real   , allocatable, dimension(:,:,:)   :: qexcp, hexcp 
   real   , allocatable, dimension(:,:,:)   :: advf_t, sgsf_t, sgsf_q, thsrc, rtsrc
   real   , allocatable, dimension(:,:,:)   :: clsrc, nlsrc, nisrc, usrc, vsrc, src_buoy
   real   , allocatable, dimension(:,:,:)   :: revsu_gf, prfil_gf,var3d_agf
   real   , allocatable, dimension(:,:,:)   :: var3d_bgf,var3d_cgf,var3d_dgf
   real   , allocatable, dimension(:,:,:)   :: zm3d, zt3d, dm3d, up, vp, wp
   real   , allocatable, dimension(:,:,:)   :: temp, press, rvap, buoy_exc, gsf_t, gsf_q
   real   , allocatable, dimension(:,:,:)   :: cprr4d_tmp, xmb4d_tmp, edt4d_tmp
   real   , allocatable, dimension(:,:,:)   :: pwav4d_tmp, sigma4d_tmp
   !- 4d
   real   , allocatable, dimension(:,:,:,:) :: clwup5d_tmp, tup5d_tmp, conv_cld_fr5d_tmp
   real   , allocatable, dimension(:,:,:,:) :: dd_massdetr5d_tmp ,zup5d_tmp, zdn5d_tmp
   real   , allocatable, dimension(:,:,:,:) :: prup5d_tmp, prdn5d_tmp, dd_massentr5d_tmp
   real   , allocatable, dimension(:,:,:,:) :: pcup5d_tmp, up_massentr5d_tmp, up_massdetr5d_tmp 
   real   , allocatable, dimension(:,:,:,:) :: src_chem
   real   , allocatable, dimension(:,:,:,:) :: tracer
   real   , allocatable, dimension(:,:,:,:) :: mp_ice, mp_liq, mp_cf 
   real   , allocatable, dimension(:,:,:,:) :: sub_mpqi, sub_mpql, sub_mpcf

   p_use_gate = .true.

   !------------------- simulation begins  ------------------
   !
   !- reads namelists
   open(15,file='gf.inp',status='old',form='formatted')    
    read(15,nml=run)
   close(15)

   modConvParGF_initialized = .false.
   init_stat = initModConvParGF()
   !-- read the GF namelist
   call readGFConvParNML(mynum)

   if(trim(rundata) == "gate.dat" .or. trim(rundata) == "GATE.dat" ) then
     klon_local=p_klon
     klev_local=p_klev  
     output_sound = 0
   else
     output_sound = 1
     klon_local=1
!    klon_local=100 !40
     klev_local=klev_sound
   endif

   mzp = klev_local
   kte = mzp
   kme = mzp
   mgmzp=mzp

    print *, "mxp, myp, mzp               : ", mxp, myp, mzp
    print *, "mtp, nmp, itime1, maxiens   : ", mtp, nmp, itime1, maxiens
    print *, "ims, ime, jms, jme, kms, kme: ", ims,ime, jms,jme, kms,kme
    print *, "its, ite, jts, jte, kts, kte: ", its,ite, jts,jte, kts,kte
    
   allocate(flip(mzp),fscav_int(mtp))
   !
   allocate(kpbl(mxp,myp),do_this_column(mxp,myp))
   allocate(wlpool(mxp,myp))
   allocate(dx2d(mxp,myp), lons(mxp,myp), lats(mxp,myp), aot500(mxp,myp), tke_pbl(mxp,myp), col_sat(mxp,myp))
   allocate(glon(mxp,myp), glat(mxp,myp))
   allocate(stochastic_sig(mxp,myp), temp2m(mxp,myp), xland(mxp,myp), sfc_press(mxp,myp))
   allocate(conprr(mxp,myp),lightn_dens(mxp,myp), rh_dicy_fct(mxp,myp))
   allocate(sflux_r(mxp,myp), sflux_t(mxp,myp), topt(mxp,myp), var2d(mxp,myp))
   allocate(cnv_frc(mxp,myp),aa0(mxp,myp),aa1(mxp,myp),aa2(mxp,myp),aa3(mxp,myp))
   allocate(aa1_bl(mxp,myp),aa1_cin(mxp,myp),tau_bl(mxp,myp),tau_ec(mxp,myp))
   allocate(aa1_adv(mxp,myp), aa1_radpbl(mxp,myp), rtgt(mxp,myp))
  
   !
   allocate(advf_t(mzp,mxp,myp), sgsf_t(mzp,mxp,myp), sgsf_q(mzp,mxp,myp), thsrc(mzp,mxp,myp), rtsrc(mzp,mxp,myp))
   allocate(clsrc(mzp,mxp,myp), nlsrc(mzp,mxp,myp), nisrc(mzp,mxp,myp), usrc(mzp,mxp,myp), vsrc(mzp,mxp,myp), src_buoy(mzp,mxp,myp))
   allocate(revsu_gf(mzp,mxp,myp), prfil_gf(mzp,mxp,myp),var3d_agf(mzp,mxp,myp))
   allocate(var3d_bgf(mzp,mxp,myp),var3d_cgf(mzp,mxp,myp),var3d_dgf(mzp,mxp,myp))
   allocate(zm3d(mzp,mxp,myp), zt3d(mzp,mxp,myp), dm3d(mzp,mxp,myp), up(mzp,mxp,myp), vp(mzp,mxp,myp), wp(mzp,mxp,myp))
   allocate(temp(mzp,mxp,myp), press(mzp,mxp,myp), rvap(mzp,mxp,myp), buoy_exc(mzp,mxp,myp), gsf_t(mzp,mxp,myp), gsf_q(mzp,mxp,myp))
   allocate(qexcp(mzp,mxp,myp),hexcp(mzp,mxp,myp))

   allocate(cprr4d_tmp(mxp,myp,maxiens), xmb4d_tmp(mxp,myp,maxiens), edt4d_tmp(mxp,myp,maxiens))
   allocate(pwav4d_tmp(mxp,myp,maxiens), sigma4d_tmp(mxp,myp,maxiens))
   allocate(kstabi4d_tmp(mxp,myp,maxiens), kstabm4d_tmp(mxp,myp,maxiens))
   allocate(ierr4d_tmp(mxp,myp,maxiens), jmin4d_tmp(mxp,myp,maxiens) ,klcl4d_tmp(mxp,myp,maxiens))
   allocate(k224d_tmp(mxp,myp,maxiens), kbcon4d_tmp(mxp,myp,maxiens), ktop4d_tmp(mxp,myp,maxiens))
    !
   allocate(clwup5d_tmp(mxp,myp,mzp,maxiens), tup5d_tmp(mxp,myp,mzp,maxiens), conv_cld_fr5d_tmp(mxp,myp,mzp,maxiens))
   allocate(dd_massdetr5d_tmp(mxp,myp,mzp,maxiens) ,zup5d_tmp(mxp,myp,mzp,maxiens), zdn5d_tmp(mxp,myp,mzp,maxiens))
   allocate(prup5d_tmp(mxp,myp,mzp,maxiens), prdn5d_tmp(mxp,myp,mzp,maxiens), dd_massentr5d_tmp(mxp,myp,mzp,maxiens))
   allocate(pcup5d_tmp(mxp,myp,mzp,maxiens), up_massentr5d_tmp(mxp,myp,mzp,maxiens), up_massdetr5d_tmp(mxp,myp,mzp,maxiens))
   allocate(src_chem(mtp,mzp,mxp,myp), tracer(mxp,myp,mzp,mtp))
   allocate(mp_ice(nmp,mzp,mxp,myp), mp_liq(nmp,mzp,mxp,myp), mp_cf(nmp,mzp,mxp,myp))
   allocate(sub_mpqi(nmp,mzp,mxp,myp), sub_mpql(nmp,mzp,mxp,myp), sub_mpcf(nmp,mzp,mxp,myp))

   !--- allocation      
   allocate(cupout(0:p_nvar_grads))
   do nvar=0,p_nvar_grads
        allocate(cupout(nvar)%varp(klon_LOCAL,KLEV_LOCAL))
        allocate(cupout(nvar)%varn(3))
        cupout(nvar)%varp(:,:)=0.0
        cupout(nvar)%varn(:)  ="xxxx"
   enddo
   if(.not. p_use_gate) then
       !print*,"====================================================================="
       !print*, "use_gate logical flag must be true to run in 1-d, model will stop"
       !print*,"====================================================================="
       stop "use_gate flag"
   endif


   !  
   !- reads gate soundings                
   IF(trim(RUNDATA) == "GATE.dat") THEN
   open(7,file="GATE.dat",form="formatted",STATUS="OLD")
     read(7,*)
     do jl=1,p_klon
     	read(7,*)
     	!z(m)  p(hpa) t(c) q(g/kg) u  v (m/s) w(pa/s) q1 q2 !!!qr (k/d) advt(k/d) advq(1/s)
     	do jk=p_klev,1,-1
     	read(7,*)pgeo(jl,jk),ppres(jl,jk),ptemp(jl,jk),pq(jl,jk),pu(jl,jk),pv(jl,jk),pvervel(jl,jk), &
     		       zq1(jl,jk),zq2(jl,jk),zqr(jl,jk),zadvt(jl,jk),zadvq(jl,jk)			    
     	!print*,"GATE=",jl,jk,pgeo(jl,jk),zadvq(jl,jk)
     	end do
     enddo
   close(7)
   ENDIF


   !-
   !-
   !- general  initialization ---------------------------------------
   flip      = 1      !integer
   dx2d      = 22000. !meters
   dtlt      = 450.   !seconds
   time      = 0.
   zt3d      = 0.
   zm3d      = 0.
   dm3d      = 1.
   lons      = 0.
   lats      = 0.
   !
   RTGT     (:,:) = 1.     !don't change this
   aot500   (:,:) = 0.1    ! #
   temp2m   (:,:) = 303.   ! Kelvin
   sflux_r  (:,:) = 700./(1.15*c_alvl) !(kg/kg/s)
   sflux_t  (:,:) = 100./(1.15*c_cp) !(K/s)
   CONPRR   (:,:) = 0.

   topt     (:,:) = 0.
   sfc_press(:,:) = 1000.
   kpbl     (:,:) = 5
   wlpool   (:,:) = 5.
   if(land) then
     xland  (:,:) = 0. !land
   else
     xland  (:,:) = 1. !ocean
   endif

   ! up     (:,:,:)= 1.
   ! vp     (:,:,:)= 1.
   ! theta  (:,:,:)= 300. 
   ! thp    (:,:,:)= 300. 
   ! pp     (:,:,:)= 1000.
   ! pi0    (:,:,:)= 0.1
   ! wp     (:,:,:)= 0.
   ! rv     (:,:,:)= 0.001
   ! rv2    (:,:,:)= 0.001
   ! rtp    (:,:,:)= 0.001
   ! tend_pt(:,:,:)= 0.
   ! sub_qils      = 0.
   ! sub_qlls      = 0.
   ! sub_cfls      = 0.
   ! ls_ice        = 0.
   ! ls_liq        = 0.
   ! ls_clfrac     = 0.
   ! advf_t        = 0.
   ! rh_dicy_fct   = 0.
   ! thsrc  (:,:,:)= 0.
   ! rtsrc  (:,:,:)= 0.
   ! clsrc  (:,:,:)= 0.
   ! nlsrc  (:,:,:)= 0.
   ! nisrc  (:,:,:)= 0.
   ! usrc   (:,:,:)= 0.
   ! vsrc   (:,:,:)= 0.
   ! rcp      (:,:,:)= 0.
   ! tkep     (:,:,:)= tkmin
   ! sc_aer     =0.
   ! sc_chem    =0.
   ! src_aer    =0.
   ! src_chem   =0.
   ! buoy_exc   =0.
   ! qexp       =0.
   ! hexcp      =0.
   ! col_sat    =1.
   ! lightn_dens=0.
   ! revsu_gf   =0.
   ! pfil_gf    =0.
   ! var3d_agf  =0.
   ! var3d_bgf  =0.
   ! var3d_cgf  =0.
   ! var3d_dgf  =0.
   ! var2d      =0.
   ! zkbcon     =0.
   stochastic_sig =1.
   do_this_column =1 !integer
   tke_pbl        =tkmin
!- end of  initialization ---------------------------------------


!- big loop on the gate soundings

   icrf=1
   do icr=1, icrf     !CR: loop para dar volume de execucao no codigo      
   
   do jl=1,klon_LOCAL !klon=number of soundings
   !do jl=1,1 !klon=number of soundings

     TIME=TIME+DTLT
     !IF(TIME/86400. > 2.) CYCLE

     !write(0,*) "############ Sounding:",jl!,TIME/86400.
     !grid_length= float(jl)*1000.
     
  !-initialization 
    ierr4d_tmp         =0               
    jmin4d_tmp         =0 
    klcl4d_tmp         =0
    k224d_tmp          =0 
    kbcon4d_tmp        =0 
    ktop4d_tmp         =0 
    kstabi4d_tmp       =0 
    kstabm4d_tmp       =0

    xmb4d_tmp          =0. 
    cprr4d_tmp         =0. 
    edt4d_tmp          =0. 
    pwav4d_tmp         =0. 
    sigma4d_tmp        =0.
    pcup5d_tmp         =0. 
    up_massentr5d_tmp  =0.        
    up_massdetr5d_tmp  =0.
    dd_massentr5d_tmp  =0.
    dd_massdetr5d_tmp  =0.
    zup5d_tmp          =0.
    zdn5d_tmp          =0. 
    prup5d_tmp         =0. 
    prdn5d_tmp         =0. 
    clwup5d_tmp        =0. 
    tup5d_tmp          =0.  
    conv_cld_fr5d_tmp  =0.                          
    REVSU_GF           =0.

   !if(JL .ne. 40) cycle
   !print*," ====================================================================="
   !print*,"Sounding =",jl

    print*, "Processando GF"
    CALL convParGFDriver(mxp,myp,KLEV_LOCAL,mtp ,nmp, time, itime1 &
                      ,ims,ime, jms,jme, kms,kme   &
                      ,its,ite, jts,jte, kts,kte   &
                      ,flip        &
                      ,fscav_int   &
                      ,mynum       &
                      ,dtlt        &
                      ,dx2d        &
                      ,stochastic_sig &
                      ,zm3d        &
                      ,zt3d        &
                      ,dm3d        &
                      !--- sfc inputs 
                      ,lons        &
                      ,lats        &
                      ,aot500      &
                      ,temp2m      &
                      ,sflux_r     &
                      ,sflux_t     &
                      ,qexcp       &
                      ,hexcp       & 
                      ,wlpool      &
                      ,topt        &
                      ,xland                 &
                      ,sfc_press             &
                      ,kpbl                  &
                      ,tke_pbl               &
                      !--- atmos state
                      ,col_sat&
                      ,up     &
                      ,vp     &
                      ,wp     &
                      ,temp   &
                      ,press  &
                      ,rvap   &
                      ,mp_ice &
                      ,mp_liq &
                      ,mp_cf  &
                      ,rvap   &
                      !--- atmos composition state
                      ,TRACER   & !- note: uses GEOS-5 data structure
                      !---- forcings---
                      ,buoy_exc &
                      , gsf_t   & ! forcing for theta adv+rad
                      , gsf_q   & ! forcing for rv    adv
                      ,advf_t   &
                      ,sgsf_t   & ! forcing for theta pbl
                      ,sgsf_q   & ! forcing for rv    pbl
                      !---- output ----
                      ,conprr     &
                      ,lightn_dens& 
                      ,rh_dicy_fct&
                      ,thsrc      & ! temp tendency
                      ,rtsrc      & ! rv tendency
                      ,clsrc      & ! cloud/ice  mass   mix ratio tendency
                      ,nlsrc      & ! cloud drop number mix ratio tendency
                      ,nisrc      & ! ice        number mix ratio tendency
                      ,usrc       & ! u tendency
                      ,vsrc       & ! v tendency
                      ,sub_mpqi    & 
                      ,sub_mpql    & 
                      ,sub_mpcf    & 
                      ,src_buoy    &
                      ,src_chem    & ! tracer tendency
                      ,revsu_gf    &
                      ,prfil_gf    & 
                      !
                      ,do_this_column    &
                      ,ierr4d_tmp        & 
                      ,jmin4d_tmp        &
                      ,klcl4d_tmp        &
                      ,k224d_tmp         &
                      ,kbcon4d_tmp       &
                      ,ktop4d_tmp        &
                      ,kstabi4d_tmp      &
                      ,kstabm4d_tmp      &
                      ,cprr4d_tmp        &
                      ,xmb4d_tmp         &
                      ,edt4d_tmp         &
                      ,pwav4d_tmp        &
                      ,sigma4d_tmp       &
                      ,pcup5d_tmp        &
                      ,up_massentr5d_tmp &
                      ,up_massdetr5d_tmp &
                      ,dd_massentr5d_tmp &
                      ,dd_massdetr5d_tmp &
                      ,zup5d_tmp         &
                      ,zdn5d_tmp         &
                      ,prup5d_tmp        &
                      ,prdn5d_tmp        &
                      ,clwup5d_tmp       &
                      ,tup5d_tmp         &
                      ,conv_cld_fr5d_tmp &
                      !-- for debug/diagnostic
                      ,aa0,aa1,aa1_adv,aa1_radpbl,aa1_bl,aa2,aa3,AA1_CIN,tau_bl,tau_ec  &
                      ,var2d,var3d_agf,var3d_bgf,var3d_cgf,var3d_dgf)
                   
   enddo ! loop over gate soundings				     
   enddo
   !!  
   !
   !-- output
   print*,"writing grads control file:',trim(runname)//'.ctl"
   !
   !number of variables to be written
   nvartotal=0
   do nvar=0,p_nvar_grads
     if(cupout(nvar)%varn(1) .ne. "xxxx") nvartotal=nvartotal+1
     if(cupout(nvar)%varn(3)  ==  "3d"  ) klevgrads(nvar)=KLEV_LOCAL-1
     if(cupout(nvar)%varn(3)  ==  "2d"  ) klevgrads(nvar)=1
   enddo
  !- binary file 
   inquire (iolength=int_byte_size) real_byte_size  ! inquire by output list
   print*, 'opening grads file:',trim(runname)//'.gra'
   open(19,file= '../dataout_1_coluna/'//trim(runname)//'.gra',form='unformatted',&
           access='direct',status='replace', recl=int_byte_size*(klon_LOCAL))
   nrec=0
   do nvar=0,p_nvar_grads
       if(cupout(nvar)%varn(1) .ne. "xxxx") then
        do jk=1,klevgrads(nvar)
          nrec=nrec+1
          write(19,REC=nrec) real((cupout(nvar)%varp(:,jk)),4)
        enddo
       endif
   enddo

   close (19)

   !-setting vertical dimension '0' for 2d var
   where(klevgrads==1)klevgrads=0
   !- ctl file
   open(20,file='../dataout_1_coluna/'//trim(runname)//'.ctl',status='unknown')
   write(20,2001) '^'//trim(runname)//'.gra'
   write(20,2002) 'undef -9.99e33'
   write(20,2002) 'options'!byteswapped' ! zrev'
   write(20,2002) 'title '//trim(runlabel)
   write(20,2003) 1,0.,1. ! units m/km
   write(20,2004) klon_LOCAL,1.,1.

   IF(trim(RUNDATA) == "GATE.dat") THEN
     write(20,2005) KLEV_LOCAL-1,(ppres(1,jk),jk=1,KLEV_LOCAL-1)
   ELSE
    n1 = KLEV_LOCAL/3
    write(20,2005) KLEV_LOCAL-1,(cupout(0)%varp(1,jk),jk=1,n1)
    n2 = n1 + KLEV_LOCAL/3
    write(20,2009)            (cupout(0)%varp(1,jk),jk=n1+1,n2)
    write(20,2009)            (cupout(0)%varp(1,jk),jk=n2+1,KLEV_LOCAL-1)
   ENDIF
   
   write(20,2006) 1,'00:00Z01JAN2000','1mn'
   write(20,2007) nvartotal
   do nvar=0,p_nvar_grads
    if(cupout(nvar)%varn(1) .ne. "xxxx") then
     write(20,2008) cupout(nvar)%varn(1)(1:len_trim(cupout(nvar)%varn(1)))&
                   ,klevgrads(nvar),cupout(nvar)%varn(2)(1:len_trim(cupout(nvar)%varn(2)))
    endif
   enddo
  
   write(20,2002) 'endvars'
   close(20)
 
  2001 format('dset ',a)
  2002 format(a)
  2003 format('xdef ',i4,' linear ',2f15.3)
  2004 format('ydef ',i4,' linear ',2f15.3)

  2005 format('zdef ',i4,' levels ',200f10.2)
  2009 format(200f10.2)

  2006 format('tdef ',i4,' linear ',2a15)
  2007 format('vars ',i4)
  2008 format(a10,i4,' 99 ',a40)!'[',a8,']')
  2055 format(60f7.0)
   133 format (1x,F7.0)

END PROGRAM GF_1d_driver

