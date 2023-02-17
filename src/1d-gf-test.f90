program GF_1d_driver

  USE modGate
  use modConvParGF  , only: modConvParGFDriver,icumulus_gf, closure_choice, p_deep, p_shal, p_mid &
      ,use_scale_dep,dicycle,tau_deep,tau_mid,hcts                       &
      ,use_tracer_transp, use_tracer_scaven,use_memory,convection_tracer &
      ,use_flux_form,use_tracer_evap,downdraft,use_fct                   &
      ,use_rebcb, vert_discr, satur_calc, clev_grid, apply_sub_mp, alp1  &
      ,sgs_w_timescale, lightning_diag, tau_ocea_cp,  tau_land_cp        &
      ,autoconv, bc_meth,overshoot,use_wetbulb                           &
      ,c1,c0_deep, qrc_crit,lambau_deep,lambau_shdn,c0_mid               &
      ,cum_max_edt_land  ,cum_max_edt_ocean, cum_hei_down_land           &
      ,cum_hei_down_ocean,cum_hei_updf_land, cum_hei_updf_ocean          &
      ,use_momentum_transp,cum_entr_rate                                 &
      ,p_nmp, p_lsmp, p_cnmp,moist_trigger,frac_modis,max_tq_tend  &
      ,cum_fadj_massflx, cum_use_excess, cum_ave_layer, adv_trigger      &
      ,use_smooth_prof, output_sound,use_cloud_dissipation      &
      ,use_smooth_tend,GFConvparInit,beta_sh,c0_shal                   &
      ,use_linear_subcl_mf,cap_maxs,liq_ice_number_conc,alpha_adv_tuning &
      ,sig_factor,lcl_trigger, rh_dicycle, add_coldpool_prop,cum_t_star  &
      ,add_coldpool_clos,mx_buoy1, mx_buoy2, cum_t_star,cum_zuform       &
      ,add_coldpool_diff, modConvParGF_initialized, initModConvParGF


  implicit none

  integer,parameter :: mzp=200 ,mxp=1   ,myp=1

  character(len=*), parameter :: procedureName = '1d-gf-test' 

  integer,parameter :: &
       ngrid=1 ,ngrids_cp=1     & 
     ,iens  =1   , mynum=1 ,npatch=1        &
     ,i0    =1   , j0=1                     &    
     ,ia    =1   , ja=1 ,iz=1 ,jz=1         ! 1-d column
     
  integer,parameter :: mgmxp=mxp &
                      ,mgmyp=myp &
                      ,mgmzp=mzp

  integer,parameter :: &
     ims=1   ,ime=mxp ,jms=1   ,jme=myp ,kms=1, kme=mzp  &   
    ,its=ia  ,ite=iz  ,jts=ja  ,jte=jz  ,kts=1   

  
  integer, parameter ::   maxiens    = 3 !cloud spectral size
     
   integer, parameter :: n_aer=1
   real     :: fscav(n_aer)

   real, parameter ::      &
       rgas    = 287.,    &
       cp      = 1004.,   &
       rm      = 461.,    &
       p00     = 1.e5,    &
       tcrit   = 273.15,  &
       g       = 9.80,    &
       cpor    = cp/rgas, &
       xl      = 2.5e6,   &
       akmin   = 1.0,     &
       tkmin   = 1.e-5
     
   real,    dimension(kms:kme)  ::  zmn,ztn
   integer, dimension(kms:kme)  :: flip

   real,  dimension(kms:kme ,  ims:ime , jms:jme ) :: &
                                                          up,   &
                                                          vp,   &
                                                          wp,   &
                                                          rv,   &
                                                          rtp,  &
                                                          theta   ,& 
                                                          thp     ,& 
                                                          pp      ,& 
                                                          pi0     ,& 
                                                          dn0     ,&
                                                          tend_pt ,&
                                                          rcp,tkep,zt3d,zm3d,dm3d,rv2,&
                                                          buoy_exc, qexp, hexcp 
           
   real,  dimension(p_nmp,kms:kme ,  ims:ime , jms:jme ) ::  ls_ice,  &
                                                           ls_liq,  &
                                                           ls_clfrac

   real, dimension(ims:ime , jms:jme,npatch)  :: patch_area


   real, dimension( ims:ime , jms:jme ) ::  aot500 , temp2m , rtgt &
                                           ,sflux_r, sflux_t, topt &
                       ,rshort , xland,sfc_press&
                       ,lons,lats,col_sat,cnv_frc, tke_pbl,wlpool

   real :: dtlt,  time

   real, dimension( ims:ime , jms:jme ) ::      &
                          conprr,xmb_deep,                     &
                          apr_gr,apr_w,apr_mc,apr_st,apr_as,    &
                          xmb_shallow,err_deep ,mass_flux_dp ,mass_flux_sh  ,mass_flux_md     
   real, dimension( ims:ime , jms:jme ) ::        &
                  weight_gr,weight_w,weight_mc,weight_st,weight_as &
                 ,aa0,aa1,aa2,aa3,aa1_bl,aa1_cin,tau_bl,tau_ec &
                 ,xmbdn,lightn_dens,var2d,zkbcon,stochastic_sig, rh_dicy_fct, aa1_adv, aa1_radpbl
              
   real,  dimension(ims:ime, kms:kme ,  jms:jme, n_aer) :: sc_aer,sc_chem
   real,  dimension(n_aer, kms:kme , ims:ime ,  jms:jme) :: src_aer,src_chem

!
   real, dimension(kms:kme , ims:ime ,  jms:jme )   ::   &  
          thsrc  & ! temp tendency         
         ,rtsrc  & ! rv tendency           
         ,clsrc  & ! cloud/ice tendency        
         ,nisrc  & ! ice # tendency        
         ,nlsrc  & ! cloud # tendency        
         ,usrc   & ! u tendency            
         ,vsrc   & ! v tendency            
         ,mup    & ! updraft mass flux         
         ,lsfth  & ! forcing for theta deep    
         ,advf_t     &
         ,lsfrt      & ! forcing for rv deep       
         ,lsfth_sh   & ! forcing for theta shallow 
         ,lsfrt_sh   &   ! forcing for rv shallow    
         ,src_buoy   & ! lifting tendency from downdrafts
         ,revsu_gf   &   
         ,pfil_gf    & ! ice_or_liq convective_precipitation flux: kg m2 s-1 (deep only)
         ,var3d_agf  & ! dummy 3-d var for output
         ,var3d_bgf  & ! dummy 3-d var for output
         ,var3d_cgf  & ! dummy 3-d var for output
         ,var3d_dgf    ! dummy 3-d var for output

   integer :: itime1

   real, dimension(p_nmp,kms:kme , ims:ime ,  jms:jme )   ::   &  
              sub_qils & ! subsidence transport applied to grid-scale ice mix ratio
             ,sub_qlls & ! subsidence transport applied to grid-scale cloud mix ratio
             ,sub_cfls   ! subsidence transport applied to grid-scale cloud fraction
!
!- for convective transport-start
  integer, dimension(mgmxp,mgmyp,maxiens,ngrids_cp) ::         &
           ierr4d  		     &  	     
	      ,jmin4d  		     & 
	      ,kdet4d  		     & 
	      ,k224d	         & 
	      ,kbcon4d 		     & 
	      ,ktop4d  		     & 
	      ,kpbl4d  		     & 
          ,klcl4d            &
	      ,kstabi4d		     & 
	      ,kstabm4d          &
          ,do_this_column		   

   real,dimension(mgmxp,mgmyp,maxiens,ngrids_cp) :: &
           cprr4d            & 
	      ,xmb4d		     & 
	      ,edt4d		     & 
	      ,sigma4d		     & 
	      ,pwav4d		     	      
   real,dimension(mgmzp,mgmxp,mgmyp,maxiens,ngrids_cp) ::&
 	       pcup5d 		     & 
          ,up_massentr5d	 &        
	      ,up_massdetr5d	 &
	      ,dd_massentr5d	 &
	      ,dd_massdetr5d	 &
 	      ,zup5d		     &
	      ,zdn5d   		     & 
	      ,prup5d  		     & 
	      ,prdn5d  		     & 
	      ,clwup5d 		     & 
	      ,tup5d             &
          ,conv_cld_fr5d                  		      
  
  real,    dimension(ims:ime , jms:jme ) :: grid_length
  integer, dimension(ims:ime , jms:jme ) ::  kpbl
  logical :: land
 !- here are the place for data related with the gate soundings
 !- soundings arrays
   integer ::jk, nruns, version, klon_local,klev_local,kte

!- this for the namelist gf.inp
  namelist /run/ runname, runlabel, rundata,version, land , klev_sound 

!- for grads output
   integer :: nrec,nvx,nvar,nvartotal,klevgrads(0:300),int_byte_size,n1,n2,n3
   real    :: real_byte_size
   logical :: init_stat
   
   
!C: Carlos Renato testes: 
   integer :: icr, icrf

!Informa que não houve inicialização do GF5
modConvParGF_initialized = .false.
! Chama a rotina que inicializa as variáveis do módulo
init_stat = initModConvParGF()
!------------------- simulation begins  ------------------
!
!- reads namelists
!- namelist 1 
   open(15,file='gf.inp',status='old',form='formatted')	  
    read(15,nml=run)
   close(15)

!- namelist 2 
  call GFConvparInit(mynum)


!- print the namelist
  print*,"           "
  print*,"------------ namelist for:  ",trim(runname)
  print*, 'rundata           ' , trim(rundata)	  
  print*, 'icumulus_gf       ' , icumulus_gf	  
  print*, 'closure_choice    ' , closure_choice   
  print*, 'clev_grid         ' , clev_grid	  
  print*, 'use_rebcb         ' , use_rebcb	  
  print*, 'vert_discr        ' , vert_discr	  
  print*, 'satur_calc        ' , satur_calc	  
  print*, 'autoconv          ' , autoconv	  
  print*, 'bc_meth           ' , bc_meth	
  print*, 'overshoot         ' , overshoot	  
  print*, 'use_wetbulb       ' , use_wetbulb	  
  print*, 'hei_down_land     ' , real(cum_hei_down_land  ,4)
  print*, 'hei_down_ocean    ' , real(cum_hei_down_ocean ,4)  
  print*, 'hei_updf_land     ' , real(cum_hei_updf_land  ,4)	    
  print*, 'hei_updf_ocean    ' , real(cum_hei_updf_ocean ,4)  
  print*, 'max_edt_land      ' , real(cum_max_edt_land 	 ,4) 
  print*, 'max_edt_ocean     ' , real(cum_max_edt_ocean	 ,4) 
  print*, 'cum_use_excess    ' , cum_use_excess	  
  print*, 'c0_deep           ' , real(c0_deep	         ,4)	  
  print*, 'c0_mid            ' , real(c0_mid 	         ,4)	  
  print*, 'c0_shal           ' , real(c0_shal 	         ,4)	  
  print*, 'c1                ' , real(c1		     ,4)	  
  print*, 'qrc_crit          ' , real(qrc_crit	 	 ,4) 
  print*, 'lambau_deep       ' , real(lambau_deep	 ,4)	  
  print*, 'lambau_shdn       ' , real(lambau_shdn	 ,4)	  
  print*, 'use_tracer_transp ' , use_tracer_transp
  print*, 'use_tracer_scaven ' , use_tracer_scaven
  print*, 'use_flux_form     ' , use_flux_form    
  print*, 'use_fct           ' , use_fct	  
  print*, 'use_tracer_evap   ' , use_tracer_evap  
  print*, 'use_momentum_transp', use_momentum_transp   
  print*, 'downdraft         ' , downdraft	  
  print*, 'sgs_w_timescale   ' , sgs_w_timescale  
  print*, 'apply_sub_mp      ' , apply_sub_mp	  
  print*, 'alp1              ' , real(alp1	 ,4)		  
  print*, 'lightning_diag    ' , lightning_diag   
  print*, 'use_scale_dep     ' , USE_SCALE_DEP    
  print*, 'dicycle	         ' , DICYCLE	  
  print*, 'land 	         ' , LAND		  
  print*, 'CUM_ENTR 	     ' , real(cum_entr_rate	 ,4)
  print*, 'FRAC_MODIS	     ' , FRAC_MODIS
  print*, 'MOIST_TRIGGER     ' , MOIST_TRIGGER	  
  print*, 'ADV_TRIGGER       ' , ADV_TRIGGER	
  print*, 'tau_deep,tau_mid  ' , real(tau_deep,4),real(tau_mid,4)
  print*, 'use_smooth_prof   ' , use_smooth_prof
  print*, 'use_cloud_dissipation', real(use_cloud_dissipation,4)
  print*, 'beta_sh 	         ' , real(beta_sh,4)
  print*, 'use_linear_subcl_mf',use_linear_subcl_mf
  print*, 'cap_maxs          ',real(cap_maxs,4)
  print*,"========================================"
  print*,"           "

   
   IF(trim(RUNDATA) == "GATE.dat") THEN
     KLON_LOCAL=p_klon
     KLEV_LOCAL=p_klev	
     OUTPUT_SOUND = 0
   ELSE
     OUTPUT_SOUND = 1
     KLON_LOCAL=1
!    KLON_LOCAL=100 !40
     KLEV_LOCAL=klev_sound
   ENDIF
!--- allocation      
   allocate(cupout(0:p_nvar_grads))
   do nvar=0,p_nvar_grads
        allocate(cupout(nvar)%varp(klon_LOCAL,KLEV_LOCAL))
        allocate(cupout(nvar)%varn(3))
        cupout(nvar)%varp(:,:)=0.0
        cupout(nvar)%varn(:)  ="xxxx"
   enddo
   print*,"USE_GATE=",p_use_gate
   if(.not. p_use_gate) then
       print*,"====================================================================="
       print*, "use_gate logical flag must be true to run in 1-d, model will stop"
       print*,"====================================================================="
       stop "use_gate flag"
   endif
!  
!- reads gate soundings                
   IF(trim(RUNDATA) == "GATE.dat") THEN
   print*,"reading GATE soundings"
   open(7,file="GATE.dat",form="formatted",STATUS="OLD")
     read(7,*)
     do jl=1,p_klon
     	read(7,*)
     	!z(m)  p(hpa) t(c) q(g/kg) u  v (m/s) w(pa/s) q1 q2 !!!qr (k/d) advt(k/d) advq(1/s)
     	do jk=p_klev,1,-1
     	read(7,*)pgeo(jl,jk),ppres(jl,jk),ptemp(jl,jk),pq(jl,jk),        &
     		 pu(jl,jk),pv(jl,jk),pvervel(jl,jk), &
     		 zq1(jl,jk),zq2(jl,jk),zqr(jl,jk),zadvt(jl,jk),&
     		 zadvq(jl,jk)			    
     	!print*,"GATE=",jl,jk,pgeo(jl,jk),zadvq(jl,jk)
     	end do
     enddo
   close(7)
   ENDIF
!-
!-
!- general  initialization ---------------------------------------

   grid_length = 22000. ! meters
   dtlt=450. !seconds
   time=0.
   kte=KLEV_LOCAL
   	
   != nz
   flip        = 1 !integer
   zmn	       = 0.
   ztn	       = 0.
   zt3d        = 0.
   zm3d        = 0.
   dm3d        = 0.
   != nx,ny
   RTGT    (:,:) = 1.     !don�t change this
   aot500  (:,:) = 0.1    ! #
   temp2m  (:,:) = 303.   ! Kelvin
   sflux_r (:,:) = 700./(1.15*xl) !(kg/kg/s)
   sflux_t (:,:) = 100./(1.15*cp) !(K/s)
   cnv_frc (:,:) = 0.
   lons=0.
   lats=0.
   CONPRR  (:,:) = 0.
   apr_GR  (:,:) = 0.
   apr_W   (:,:) = 0.
   apr_MC  (:,:) = 0.
   apr_ST  (:,:) = 0.
   apr_AS  (:,:) = 0.
   xmb_deep(:,:) = 0.
   err_deep(:,:) = 0.
   xmb_shallow(:,:) = 0.
   weight_GR  (:,:) = 0.
   weight_W   (:,:) = 0.
   weight_MC  (:,:) = 0.
   weight_ST  (:,:) = 0.
   weight_AS  (:,:) = 0.
   topt       (:,:) = 0.
   rshort     (:,:) = 0.
   sfc_press  (:,:) = 1000.
   kpbl        (:,:) = 5
   MASS_FLUX_DP(:,:) = 0.  
   MASS_FLUX_SH(:,:) = 0.  
   wlpool      (:,:) = 5.
   != nx,ny,npatch
   patch_area (:,:,:)= 1.!don�t change this
   if(land) then
     xland (:,:) = 0. !land
   else
     xland (:,:) = 1. !ocean
   endif

   dn0    (:,:,:)= 1.
   up     (:,:,:)= 1.
   vp     (:,:,:)= 1.
   theta  (:,:,:)= 300. 
   thp    (:,:,:)= 300. 
   pp     (:,:,:)= 1000.
   pi0    (:,:,:)= 0.1
   wp     (:,:,:)= 0.
   rv     (:,:,:)= 0.001
   rv2    (:,:,:)= 0.001
   rtp    (:,:,:)= 0.001
   tend_pt(:,:,:)= 0.
   sub_qils      = 0.
   sub_qlls      = 0.
   sub_cfls      = 0.
   ls_ice        = 0.
   ls_liq        = 0.
   ls_clfrac     = 0.
   advf_t        = 0.
   rh_dicy_fct   = 0.
   thsrc  (:,:,:)= 0.
   rtsrc  (:,:,:)= 0.
   clsrc  (:,:,:)= 0.
   nlsrc  (:,:,:)= 0.
   nisrc  (:,:,:)= 0.
   usrc   (:,:,:)= 0.
   vsrc   (:,:,:)= 0.
   mup    (:,:,:)= 0.
   lsfth  (:,:,:)= 0.
   lsfrt  (:,:,:)= 0.
   lsfth_sh (:,:,:)= 0.
   lsfrt_sh (:,:,:)= 0.
   rcp      (:,:,:)= 0.
   tkep     (:,:,:)= tkmin
   sc_aer     =0.
   sc_chem    =0.
   src_aer    =0.
   src_chem   =0.
   buoy_exc   =0.
   qexp       =0.
   hexcp      =0.
   col_sat    =1.
   lightn_dens=0.
   revsu_gf   =0.
   pfil_gf    =0.
   var3d_agf  =0.
   var3d_bgf  =0.
   var3d_cgf  =0.
   var3d_dgf  =0.
   var2d      =0.
   zkbcon     =0.
   stochastic_sig =1.
   do_this_column =1 !integer
   tke_pbl        =tkmin
!- end of  initialization ---------------------------------------


!- big loop on the gate soundings

   icrf=1
   do icr=1, icrf     !CR: loop para dar volume de execucao no codigo
      print*, 'CR: ', procedureName, ' loop-bombando:', icr
      
      !do jl=1,klon_LOCAL !klon=number of soundings
         do jl=1,1 !klon=number of soundings


         TIME=TIME+DTLT
         !IF(TIME/86400. > 2.) CYCLE

         !CR:     write(0,*) "############ Sounding:",jl!,TIME/86400.
         !grid_length= float(jl)*1000.
     
     
 
         !-initialization 
         ierr4d         =0               
         jmin4d         =0 
         klcl4d         =0

         k224d          =0 
         kbcon4d        =0 
         ktop4d         =0 
         kstabi4d       =0 
         kstabm4d       =0   
         xmb4d          =0. 
         cprr4d         =0. 
         edt4d          =0. 
         pwav4d         =0. 
         sigma4d        =0.
         pcup5d         =0. 
         up_massentr5d  =0.        
         up_massdetr5d  =0.
         dd_massentr5d  =0.
         dd_massdetr5d  =0.
         zup5d          =0.
         zdn5d          =0. 
         prup5d         =0. 
         prdn5d         =0. 
         clwup5d        =0. 
         tup5d          =0.  
         conv_cld_fr5d  =0.                		      
         REVSU_GF       =0.

         !if(JL .ne. 40) cycle
         !CR:       print*," ====================================================================="
         !CR:       print*,"Sounding =",jl

   
         !CR:   print*, 'CR: 1d-gf-test.f90: chamando modConvParGFDriver', icr
         print*, 'CR: ', procedureName, ' chamando modConvParGFDriver', jl
         CALL modConvParGFDriver(mxp,myp,KLEV_LOCAL,n_aer,p_nmp, time, itime1 &
            ,ims,ime, jms,jme, kms,kme                        & 
            ,its,ite, jts,jte, kts,kte                        & 
	         ,flip        &
            , FSCAV      &
            ,mynum       &
            ,dtlt        &
            ,grid_length & 
	         ,stochastic_sig &
	         ,zm3d	       & !zmn(:,ngrid) 
	         ,zt3d        & !ztn(:,ngrid) 
	         ,dm3d        &
            ,lons        &
	         ,lats        &
            ,aot500      & ! aot at 500nm
	         ,temp2m      & ! 2m-temp
  	         ,sflux_r     & !turb_g(ngrid)%sflux_r  
            ,sflux_t     & !turb_g(ngrid)%sflux_t
            ,qexp        &
            ,hexcp       & 
            ,wlpool      &
	         ,topt	       & !grid_g(ngrid)%topt 
	         ,xland       & ! land fc
	         ,sfc_press   &
            ,KPBL        &  
            ,tke_pbl     &
            ,col_sat     &
            ,up	   & !basic_g(ngrid)%up      
            ,vp	   & !basic_g(ngrid)%vp      
            ,wp	   & !basic_g(ngrid)%wp      
            ,theta   & !basic_g(ngrid)%theta   
            ,pp	   & !basic_g(ngrid)%pp      
	         ,rv	   & !basic_g(ngrid)%rv      
	         ,ls_ice	   &
	         ,ls_liq	   &
            ,ls_clfrac   &
	         ,rv2	       & !basic_g(ngrid)%rv      
	         ,SC_CHEM     &
	         ,buoy_exc    &
            ,lsfth       & ! forcing for theta deep    cuforc_g(ngrid)% lsfth  
	         ,lsfrt 	   & ! forcing for rv deep       cuforc_g(ngrid)% lsfrt  
            ,advf_t      &
	         ,lsfth_sh    & ! forcing for theta shallow cuforc_sh_g(ngrid)%lsfth
	         ,lsfrt_sh    & ! forcing for rv shallow    cuforc_sh_g(ngrid)%lsfrt
	         ,CONPRR      & !cupout_g(ngrid)%CONPRR 
            ,LIGHTN_DENS &
            ,rh_dicy_fct &
            ,THSRC	   & ! temp tendency		   g3d_g(ngrid)%THSRC	   
            ,RTSRC	   & ! rv tendency		   g3d_g(ngrid)%RTSRC	   
            ,CLSRC	   & ! cloud/ice tendency	   g3d_g(ngrid)%CLSRC	   
            ,NLSRC       &
            ,NISRC       &
            ,USRC	       & ! U tendency		   g3d_g(ngrid)%USRC	   
            ,VSRC	       & ! V tendency		   g3d_g(ngrid)%VSRC	   
            ,SUB_QILS    & 
            ,SUB_QLLS    & 
            ,SUB_CFLS	   &   
            ,SRC_BUOY    &
	         ,SRC_CHEM    &
            ,REVSU_GF    & 
	         ,PFIL_GF     & 
	         ,do_this_column&
            ,ierr4d	    &
            ,jmin4d	    &
            ,klcl4d	    &
            ,k224d	    &
            ,kbcon4d      &
            ,ktop4d	    &
            ,kstabi4d     &
            ,kstabm4d     &
            ,cprr4d	    &
            ,xmb4d	    &
            ,edt4d	    &
            ,pwav4d	    &
	         ,sigma4d      &
            ,pcup5d	    &
            ,up_massentr5d&
            ,up_massdetr5d&
            ,dd_massentr5d&
            ,dd_massdetr5d&
            ,zup5d	    &
            ,zdn5d	    & 
            ,prup5d	    & 
            ,prdn5d	    & 
            ,clwup5d      & 
            ,tup5d	    & 
            ,conv_cld_fr5d&				       
	         !-- for debug/diagnostic
            ,AA0,AA1,AA1_ADV,AA1_RADPBL,AA1_BL,AA2,AA3,AA1_CIN,TAU_BL,TAU_EC  &
            ,VAR2d,VAR3d_aGF,VAR3d_bGF,VAR3d_cGF,VAR3d_dGF&
         )



      enddo ! loop over gate soundings		 !klon=number of soundings		     
   enddo    !CR: loop bomba

END PROGRAM GF_1d_driver

