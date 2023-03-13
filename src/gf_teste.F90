program  gf_test
   use modConvParGF, only: convParGFDriver

   implicit none 

   integer :: l_unit, icnt
   logical :: first_read, ok, exists
   real :: confrq, local_time
   character(len=4) :: ctime

   integer :: its,ite, jts,jte, kts,kte, mynum
   integer :: mxp,myp,mzp,mtp ,nmp, itime1, maxiens
   integer :: ims,ime, jms,jme, kms,kme
   real    :: time, dtlt
   !
   integer, allocatable, dimension(:)       :: flip
   integer, allocatable, dimension(:,:)     :: kpbl,do_this_column
   integer, allocatable, dimension(:,:,:)   :: kstabi4d_tmp, kstabm4d_tmp
   integer, allocatable, dimension(:,:,:)   :: ierr4d_tmp, jmin4d_tmp ,klcl4d_tmp 
   integer, allocatable, dimension(:,:,:)   :: k224d_tmp, kbcon4d_tmp, ktop4d_tmp
   !  
   real   , allocatable, dimension(:)       :: fscav_int
   real   , allocatable, dimension(:,:)     :: dx2d, lons, lats, aot500, tke_pbl, col_sat
   real   , allocatable, dimension(:,:)     :: stochastic_sig, temp2m, xland, sfc_press
   real   , allocatable, dimension(:,:)     :: conprr,lightn_dens, rh_dicy_fct
   real   , allocatable, dimension(:,:)     :: sflux_r, sflux_t, topt, var2d
   real   , allocatable, dimension(:,:)     :: cnv_frc,aa0,aa1,aa2,aa3
   real   , allocatable, dimension(:,:)     :: aa1_bl,aa1_cin,tau_bl,tau_ec
   real   , allocatable, dimension(:,:,:)   :: advf_t, sgsf_t, sgsf_q, thsrc, rtsrc
   real   , allocatable, dimension(:,:,:)   :: clsrc, nlsrc, nisrc, usrc, vsrc, src_buoy
   real   , allocatable, dimension(:,:,:)   :: revsu_gf, prfil_gf,var3d_agf
   real   , allocatable, dimension(:,:,:)   :: var3d_bgf,var3d_cgf,var3d_dgf
   real   , allocatable, dimension(:,:,:)   :: zm3d, zt3d, dm3d, up, vp, wp
   real   , allocatable, dimension(:,:,:)   :: temp, press, rvap, buoy_exc, gsf_t, gsf_q
   real   , allocatable, dimension(:,:,:)   :: cprr4d_tmp, xmb4d_tmp, edt4d_tmp
   real   , allocatable, dimension(:,:,:)   :: pwav4d_tmp, sigma4d_tmp
   real   , allocatable, dimension(:,:,:,:) :: clwup5d_tmp, tup5d_tmp, conv_cld_fr5d_tmp
   real   , allocatable, dimension(:,:,:,:) :: dd_massdetr5d_tmp ,zup5d_tmp, zdn5d_tmp
   real   , allocatable, dimension(:,:,:,:) :: prup5d_tmp, prdn5d_tmp, dd_massentr5d_tmp
   real   , allocatable, dimension(:,:,:,:) :: pcup5d_tmp, up_massentr5d_tmp, up_massdetr5d_tmp 
   real   , allocatable, dimension(:,:,:,:) :: src_chem
   real   , allocatable, dimension(:,:,:,:) :: tracer
   real   , allocatable, dimension(:,:,:,:) :: mp_ice, mp_liq, mp_cf 
   real   , allocatable, dimension(:,:,:,:) :: sub_mpqi, sub_mpql, sub_mpcf

   first_read = .true.
   ok = .true.
   icnt = 0
   local_time = 0

   do while(ok)
      icnt = icnt+1
      write(ctime,fmt="(I4.4)") int(local_time)
      inquire (file="gf_dataIn-"//ctime//".bin", exist=exists)
      if (.not. exists) then
         print *,"Arquivo gf_dataIn-"//ctime//".bin nao encontrado. Encerrando!"
         stop "!!!"
      end if
      print *, icnt,"Abrindo e lendo arquivo gf_dataIn-"//ctime//".bin"
      open(newunit = l_unit,file = "gf_dataIn-"//ctime//".bin",ACCESS = "stream", action="read", status="old")

      read(l_unit) confrq
      read(l_unit) mxp,myp,mzp,mtp ,nmp, time, itime1, maxiens
      read(l_unit) ims,ime, jms,jme, kms,kme
      read(l_unit) its,ite, jts,jte, kts,kte

      local_time = local_time + confrq

      if(first_read) then
         first_read = .false.
         print *, "mxp, myp, mzp               : ", mxp, myp, mzp
         print *, "mtp, nmp, itime1, maxiens   : ", mtp, nmp, itime1, maxiens
         print *, "ims, ime, jms, jme, kms, kme: ", ims,ime, jms,jme, kms,kme
         print *, "its, ite, jts, jte, kts, kte: ", its,ite, jts,jte, kts,kte
         
         allocate(flip(mzp),fscav_int(mtp))
         !
         allocate(kpbl(mxp,myp),do_this_column(mxp,myp))
         allocate(dx2d(mxp,myp), lons(mxp,myp), lats(mxp,myp), aot500(mxp,myp), tke_pbl(mxp,myp), col_sat(mxp,myp))
         allocate(stochastic_sig(mxp,myp), temp2m(mxp,myp), xland(mxp,myp), sfc_press(mxp,myp))
         allocate(conprr(mxp,myp),lightn_dens(mxp,myp), rh_dicy_fct(mxp,myp))
         allocate(sflux_r(mxp,myp), sflux_t(mxp,myp), topt(mxp,myp), var2d(mxp,myp))
         allocate(cnv_frc(mxp,myp),aa0(mxp,myp),aa1(mxp,myp),aa2(mxp,myp),aa3(mxp,myp))
         allocate(aa1_bl(mxp,myp),aa1_cin(mxp,myp),tau_bl(mxp,myp),tau_ec(mxp,myp))
         !
         allocate(kstabi4d_tmp(mxp,myp,maxiens), kstabm4d_tmp(mxp,myp,maxiens))
         allocate(ierr4d_tmp(mxp,myp,maxiens), jmin4d_tmp(mxp,myp,maxiens) ,klcl4d_tmp(mxp,myp,maxiens))
         allocate(k224d_tmp(mxp,myp,maxiens), kbcon4d_tmp(mxp,myp,maxiens), ktop4d_tmp(mxp,myp,maxiens))
         allocate(advf_t(mzp,mxp,myp), sgsf_t(mzp,mxp,myp), sgsf_q(mzp,mxp,myp), thsrc(mzp,mxp,myp), rtsrc(mzp,mxp,myp))
         allocate(clsrc(mzp,mxp,myp), nlsrc(mzp,mxp,myp), nisrc(mzp,mxp,myp), usrc(mzp,mxp,myp), vsrc(mzp,mxp,myp), src_buoy(mzp,mxp,myp))
         allocate(revsu_gf(mzp,mxp,myp), prfil_gf(mzp,mxp,myp),var3d_agf(mzp,mxp,myp))
         allocate(var3d_bgf(mzp,mxp,myp),var3d_cgf(mzp,mxp,myp),var3d_dgf(mzp,mxp,myp))
         allocate(zm3d(mzp,mxp,myp), zt3d(mzp,mxp,myp), dm3d(mzp,mxp,myp), up(mzp,mxp,myp), vp(mzp,mxp,myp), wp(mzp,mxp,myp))
         allocate(temp(mzp,mxp,myp), press(mzp,mxp,myp), rvap(mzp,mxp,myp), buoy_exc(mzp,mxp,myp), gsf_t(mzp,mxp,myp), gsf_q(mzp,mxp,myp))
         allocate(cprr4d_tmp(mxp,myp,maxiens), xmb4d_tmp(mxp,myp,maxiens), edt4d_tmp(mxp,myp,maxiens))
         allocate(pwav4d_tmp(mxp,myp,maxiens), sigma4d_tmp(mxp,myp,maxiens))
         !
         allocate(clwup5d_tmp(mxp,myp,mzp,maxiens), tup5d_tmp(mxp,myp,mzp,maxiens), conv_cld_fr5d_tmp(mxp,myp,mzp,maxiens))
         allocate(dd_massdetr5d_tmp(mxp,myp,mzp,maxiens) ,zup5d_tmp(mxp,myp,mzp,maxiens), zdn5d_tmp(mxp,myp,mzp,maxiens))
         allocate(prup5d_tmp(mxp,myp,mzp,maxiens), prdn5d_tmp(mxp,myp,mzp,maxiens), dd_massentr5d_tmp(mxp,myp,mzp,maxiens))
         allocate(pcup5d_tmp(mxp,myp,mzp,maxiens), up_massentr5d_tmp(mxp,myp,mzp,maxiens), up_massdetr5d_tmp(mxp,myp,mzp,maxiens))
         allocate(src_chem(mtp,mzp,mxp,myp), tracer(mxp,myp,mzp,mtp))
         allocate(mp_ice(nmp,mzp,mxp,myp), mp_liq(nmp,mzp,mxp,myp), mp_cf(nmp,mzp,mxp,myp))
         allocate(sub_mpqi(nmp,mzp,mxp,myp), sub_mpql(nmp,mzp,mxp,myp), sub_mpcf(nmp,mzp,mxp,myp))
      endif

      read(l_unit) flip
      read(l_unit) fscav_int 
      read(l_unit) mynum     
      read(l_unit) dtlt     
      read(l_unit) dx2d     
      read(l_unit) stochastic_sig
      read(l_unit) zm3d   
      read(l_unit) zt3d   
      read(l_unit) dm3d   
      read(l_unit) lons   
      read(l_unit) lats   
      read(l_unit) aot500 
      read(l_unit) temp2m 
      read(l_unit) sflux_r 
      read(l_unit) sflux_t 
      read(l_unit) topt    
      read(l_unit) xland                 
      read(l_unit) sfc_press             
      read(l_unit) kpbl                  
      read(l_unit) tke_pbl               
      read(l_unit) col_sat
      read(l_unit) up     
      read(l_unit) vp     
      read(l_unit) wp     
      read(l_unit) temp   
      read(l_unit) press  
      read(l_unit) rvap   
      read(l_unit) mp_ice 
      read(l_unit) mp_liq 
      read(l_unit) mp_cf  
      read(l_unit) rvap   
      read(l_unit) TRACER   
      read(l_unit) buoy_exc 
      read(l_unit)  gsf_t   
      read(l_unit)  gsf_q   
      read(l_unit) advf_t   
      read(l_unit) sgsf_t   
      read(l_unit) sgsf_q   
      read(l_unit) conprr  
      read(l_unit) lightn_dens
      read(l_unit) rh_dicy_fct
      read(l_unit) thsrc      
      read(l_unit) rtsrc      
      read(l_unit) clsrc      
      read(l_unit) nlsrc      
      read(l_unit) nisrc      
      read(l_unit) usrc       
      read(l_unit) vsrc       
      read(l_unit) sub_mpqi     
      read(l_unit) sub_mpql     
      read(l_unit) sub_mpcf     
      read(l_unit) src_buoy    
      read(l_unit) src_chem   
      read(l_unit) revsu_gf    
      read(l_unit) prfil_gf     
      read(l_unit) do_this_column
      read(l_unit) ierr4d_tmp    
      read(l_unit) jmin4d_tmp    
      read(l_unit) klcl4d_tmp    
      read(l_unit) k224d_tmp     
      read(l_unit) kbcon4d_tmp   
      read(l_unit) ktop4d_tmp    
      read(l_unit) kstabi4d_tmp  
      read(l_unit) kstabm4d_tmp  
      read(l_unit) cprr4d_tmp    
      read(l_unit) xmb4d_tmp     
      read(l_unit) edt4d_tmp     
      read(l_unit) pwav4d_tmp    
      read(l_unit) sigma4d_tmp   
      read(l_unit) pcup5d_tmp        
      read(l_unit) up_massentr5d_tmp 
      read(l_unit) up_massdetr5d_tmp 
      read(l_unit) dd_massentr5d_tmp 
      read(l_unit) dd_massdetr5d_tmp 
      read(l_unit) zup5d_tmp         
      read(l_unit) zdn5d_tmp         
      read(l_unit) prup5d_tmp        
      read(l_unit) prdn5d_tmp        
      read(l_unit) clwup5d_tmp       
      read(l_unit) tup5d_tmp         
      read(l_unit) conv_cld_fr5d_tmp 
      read(l_unit) AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_EC
      read(l_unit) VAR2d,VAR3d_aGF,VAR3d_bGF,VAR3d_cGF,VAR3d_dGF
      close(l_unit)

      !Chamada da GF
      CALL convParGFDriver(mxp,myp,mzp,mtp ,nmp, time, itime1 &
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
                     ,sflux_r &
                     ,sflux_t &
                     ,topt    &
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
                     ,conprr  &
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
                     ,AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_EC  &
                     ,VAR2d,VAR3d_aGF,VAR3d_bGF,VAR3d_cGF,VAR3d_dGF &
                     )


      open(newunit = l_unit,file = "gf_dataOut-"//ctime//".txt",form = "formatted", action="write", status="replace")
      write(l_unit,*) mxp,myp,mzp,mtp ,nmp, time, itime1
      write(l_unit,*) ims,ime, jms,jme, kms,kme
      write(l_unit,*) its,ite, jts,jte, kts,kte
      write(l_unit,*) cuparm_g(ngrid)%conprr
      write(l_unit,*) g3d_g(ngrid)%lightn_dens
      write(l_unit,*) g3d_g(ngrid)%thsrc      
      write(l_unit,*) g3d_g(ngrid)%rtsrc      
      write(l_unit,*) g3d_g(ngrid)%clsrc      
      write(l_unit,*) g3d_g(ngrid)%nlsrc      
      write(l_unit,*) g3d_g(ngrid)%nisrc      
      write(l_unit,*) g3d_g(ngrid)%usrc       
      write(l_unit,*) g3d_g(ngrid)%vsrc       
      write(l_unit,*) src_buoy
      write(l_unit,*) src_chem
      write(l_unit,*) revsu_gf
      write(l_unit,*) prfil_gf
      write(l_unit,*) VAR3d_aGF
      write(l_unit,*) VAR3d_bGF
      write(l_unit,*) VAR3d_cGF
      write(l_unit,*) VAR3d_dGF
      close(l_unit)

   end do




end program gf_test