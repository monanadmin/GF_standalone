module ConvPar_GF_GEOS5
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !- This convective parameterization is build to attempt                      !
   !  a smooth transition to cloud resolving scales as proposed                 !
   !  by Arakawa et al (2011, ACP). The scheme is  described                    !
   !  in the paper Grell and Freitas (ACP, 2014).                               !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !- Implemented in GEOS5 GCM by Saulo Freitas (July 2016)                     !
   !- Use the following references for this implementation:                     !
   !- Freitas et al (2018, JAMES/AGU, https://doi.org/10.1029/2017MS001251)     !
   !- Freitas et al (2021, GMD/EGU,   https://doi.org/10.5194/gmd-14-5393-2021) !
   !- Please, contact Saulo Freitas (saulo.r.de.freitas@gmail.com) for comments !
   !- questions, bugs, etc.                                                     !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !- Adapted for BRAMS 6.0 by Saulo Freitas (November 2021)                    !
   !- Refactoring by Luiz Flavio Rodrigues at 20 December 2021 (Monday)         !
   !  Keywords using ; are separeted, some loops receives exit instead goto,    !
   !  The identation was fixed and all keywords are lowercase                   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   use module_gate
   !USE MAPL
   use MAPL_ConstantsMod ! - only for GATE soundings
   !
   use Henrys_law_ConstantsMod, only: get_HenrysLawCts
   !.. USE GTMP_2_GFCONVPAR, only : GTMP_2_GFCONVPAR_interface

   !use node_mod, only: mynum

   implicit none
   !- physical constants
   real, parameter ::    &
      rgas    = 287.,    & ! J K-1 kg-1
      cp      = 1004.,   & ! J K-1 kg-1
      rv      = 461.,    & ! J K-1 kg-1
      p00     = 1.e5,    & ! hPa
      tcrit   = 258.,    & ! K
      g       = MAPL_GRAV,&! m s-2
      cpor    = cp/rgas, &
      xlv     = 2.5e6,   & ! J kg-1
      akmin   = 1.0,     & ! #
      tkmin   = 1.e-5,   & ! m+2 s-2
      ccnclean= 250.,    & ! # cm-3
      T_0     = 273.16,  & ! K
      T_ice   = 235.16,  & ! K
      xlf     = 0.333e6, & ! latent heat of freezing (J kg-1)
      max_qsat= 0.5        ! kg/kg

   private
   public  gf_geos5_interface, maxiens, icumulus_gf, closure_choice, deep, shal, mid &
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
      ,zero_diff , nmp, lsmp, cnmp,moist_trigger,frac_modis,max_tq_tend  &
      ,cum_fadj_massflx, cum_use_excess, cum_ave_layer, adv_trigger      &
      ,use_smooth_prof, evap_fix,output_sound,use_cloud_dissipation      &
      ,use_smooth_tend,GF_convpar_init,beta_sh,c0_shal                   &
      ,use_linear_subcl_mf,cap_maxs,liq_ice_number_conc,alpha_adv_tuning &
      ,sig_factor,lcl_trigger, rh_dicycle, add_coldpool_prop             &
      ,add_coldpool_clos,mx_buoy1, mx_buoy2, cum_t_star,cum_zuform       &
      ,add_coldpool_diff

   public GF_GEOS5_DRV,make_DropletNumber ,make_IceNumber,fract_liq_f    &
      ,use_gustiness, use_random_num, dcape_threshold,coldpool_start 
   !
   logical :: wrtgrads = .false.
   integer :: nrec = 0, ntimes = 0
   real    :: int_time = 0.
   !-
   !- plume spectral size
   integer, parameter  :: maxiens = 3, deep=1 ,shal=2 , mid = 3
   character(len=10),parameter,dimension(maxiens)  :: cumulus_type = (/ &
       'deep      ' &
      ,'shallow   ' &
      ,'mid       ' &
      /)
   !- number of microphysics schemes in the host model
   !integer, parameter  :: nmp = 2, lsmp = 1, cnmp = 2 ! --- for GEOS-5
    integer, parameter  :: nmp = 2, lsmp = 1, cnmp = 2 ! --- for BRAMS

   !------------------- namelist variables
   !-- plume to be activated (1 true, 0 false): deep, shallow, congestus
   integer, dimension(maxiens) :: icumulus_gf = (/1,1,1/)

   !-- choice for the closures:
   !--  deep   : 0 ensemble (all)          , 1 GR, 4 ll omega, 7 moist conv, 10 PB
   !--  shallow: 0 ensemble (all)          , 1 Wstar, 4 heat-engine, 7 BLQE, 10 TKE-based
   !--  mid    : 0 ensemble (Wstar/BLQE/PB), 1 Wstar, 2 BLQE, 3 PB, 4 PB_BL
   integer, dimension(maxiens) :: closure_choice = (/0,  7,  3/) ! deep, shallow, congestus
   integer, parameter :: shall_closures = 12 

   !-- gross entraiment rate: deep, shallow, congestus
   real,    dimension(maxiens) :: cum_entr_rate = (/&
       1.00e-4  & !deep
      ,2.00e-3  & !shallow
      ,9.00e-4  & !mid
      /)
   !-- zu updraft format : deep, shallow, congestus
   integer, dimension(maxiens) :: cum_zuform = (/20,20,20/) ! for deep: 10 or 20, congestus: 20

   integer :: USE_TRACER_TRANSP = 1 != 0/1     - default 1

   integer :: USE_TRACER_SCAVEN = 2 != 0/1/2/3 - default 2

   integer :: USE_FLUX_FORM     = 1 != 1/2/3   - default 1

   integer :: USE_FCT           = 1 != 0/1     - default 1 (only for USE_FLUX_FORM     = 2)

   integer :: USE_TRACER_EVAP   = 1 != 0/1     - default 1 (only for USE_TRACER_SCAVEN > 0)

   integer :: USE_MEMORY        = -1 != -1/0/1/2 .../10    !-

   integer :: CONVECTION_TRACER = 0 != 0/1:  turn ON/OFF the "convection" tracer

   integer :: ADD_COLDPOOL_PROP = -1 != -1,0,1,2 add coldpool propagation

   integer :: ADD_COLDPOOL_CLOS = 0 ! add the kinetic energy at leading of the gust front
   
   integer :: ADD_COLDPOOL_DIFF = 0 ! add vert/horizontal diffusion to the cold pool propaga

   integer :: USE_SCALE_DEP     = 1 != 0/1:  scale dependence flag, default = 1

   integer :: DICYCLE           = 1 != 0/1/2:  diurnal cycle closure, default = 1
                                    != 2 uses Qadv closure (Becker et al 2021)

   integer :: RH_DICYCLE        = 0 ! controls of RH on the diurnal cycle (see Tian et al 2022 GRL)

   integer :: CLEV_GRID         = 1 != 0/1/2: interpolation method to define environ state at the
                                    != cloud levels (at face layer), default = 0
                                    != clev_grid = 0 default method
                                       != clev_grid = 1 interpolation method based on Tiedtke (1989)
                                       != clev_grid = 2 for GATE soundings only

   integer :: USE_REBCB         = 1 != 0/1: turn ON/OFF rainfall evap below cloud base, default = 0

   integer :: VERT_DISCR        = 1 != 0/1: 1=new vert discretization, default = 0

   integer :: SATUR_CALC        = 1 != 0/1: 1=new saturation specific humidity calculation, default = 0

   integer :: SGS_W_TIMESCALE   = 0 != 0/1: vertical velocity for tau_ecmwf, default = 0

   integer :: LIGHTNING_DIAG    = 0 != 0/1: do LIGHTNING_DIAGgnostics based on Lopez (2016, MWR)

   integer :: APPLY_SUB_MP      = 0 != 0/1: subsidence transport applied the to grid-scale/anvil ice/liq mix
                                    !=      ratio and cloud fraction

   real    :: ALP1              = 1 != 0/0.5/1: apply subsidence transport of LS/anvil cloud fraction using
                                    !=          time implicit discretization

   integer :: USE_WETBULB       = 0 != 0/1

                                    != boundary condition determination for the plumes
   integer :: BC_METH           = 1 != 0: simple arithmetic mean around the source level
                                    != 1: mass weighted mean around the source level

   real    :: OVERSHOOT         = 0.!= 0, 1

   integer :: AUTOCONV          = 1     != 1, 3 or 4 autoconversion formulation: (1) Kessler,
                                        !  (3) Kessler with temp dependence, (4) Sundvisqt
   real    ::  C0_DEEP          = 1.5e-3!= default= 3.e-3   conversion rate (cloud to rain, m-1) - for deep      plume
   real    ::  C0_MID           = 1.5e-3!= default= 2.e-3   conversion rate (cloud to rain, m-1) - for congestus plume
   real    ::  C0_SHAL          = 0.    != default= 0.e-3   conversion rate (cloud to rain, m-1) - for shallow   plume
   real    ::  QRC_CRIT         = 2.e-4 != default= 2.e-4   kg/kg
   real    ::  C1               = 0.0   != default= 1.e-3   conversion rate (cloud to rain, m-1) - for the 'C1d' detrainment approach

   integer :: USE_MOMENTUM_TRANSP = 1   != 0/1:  turn ON/OFF conv transp of momentum
   real    ::  LAMBAU_DEEP        = 0.0 != default= 2.0 lambda parameter for deep/congestus convection momentum transp
   real    ::  LAMBAU_SHDN        = 2.0 != default= 2.0 lambda parameter for shallow/downdraft convection momentum transp

   integer :: DOWNDRAFT           = 1   != 0/1:  turn ON/OFF downdrafts, default = 1


   real    ::  TAU_DEEP           = 1800.  != deep      convective timescale
   real    ::  TAU_MID            = 900.   != congestus convective timescale

   real    :: MAX_TQ_TEND         = 100.   != max T,Q tendency allowed (100 K/day)

   integer :: ZERO_DIFF           = 0      != to get the closest solution of the stable version Dec 2019 for single-moment

   integer :: USE_SMOOTH_PROF     = 1      != 1 makes the normalized mass flux, entr and detraiment profiles smoother

   integer :: USE_SMOOTH_TEND     = 0      != 0 => OFF, > 0 produces smoother tendencies (e.g.: for 1=> makes average between k-1,k,k+1)
   !---                                              deep, shallow, congestus
   real,   dimension(maxiens) :: CUM_HEI_DOWN_LAND =(/0.30,  0.20,  0.20/)!= [0.2,0.8] height of the max Z Downdraft , default = 0.50
   real,   dimension(maxiens) :: CUM_HEI_DOWN_OCEAN=(/0.30,  0.20,  0.20/)!= [0.2,0.8] height of the max Z Downdraft , default = 0.50
   real,   dimension(maxiens) :: CUM_HEI_UPDF_LAND =(/0.35,  0.10,  0.10/)!= [0.2,0.8] height of the max Z Updraft   , default = 0.35
   real,   dimension(maxiens) :: CUM_HEI_UPDF_OCEAN=(/0.35,  0.10,  0.10/)!= [0.2,0.8] height of the max Z Updraft   , default = 0.35
   real,   dimension(maxiens) :: CUM_MAX_EDT_LAND  =(/0.40,  0.00,  0.40/)!= maximum evap fraction allowed over the land  ,default= 0.9
   real,   dimension(maxiens) :: CUM_MAX_EDT_OCEAN =(/0.30,  0.00,  0.30/)!= maximum evap fraction allowed over the ocean ,default= 0.9

   real,   dimension(maxiens) :: CUM_FADJ_MASSFLX  =(/1.00,  1.00,  1.00/)!= multiplicative factor for tunning the mass flux at cloud base
                                                                          != default = 1.0
   real,   dimension(maxiens) :: CUM_AVE_LAYER     =(/50.,   30.,   50. /)!= layer depth for average the properties
                                                                          != of source air parcels (mbar)
   real,   dimension(maxiens) :: CUM_T_STAR        =(/1.,    10.,   10. /) != scale temperature for the diurnal cycle closure

   integer,dimension(maxiens) :: CUM_USE_EXCESS    =(/1,     1,     1   /)!= use T,Q excess sub-grid scale variability

   integer :: MOIST_TRIGGER  = 0   != relative humidity effects on the cap_max trigger function
   integer :: FRAC_MODIS     = 1   != use fraction liq/ice content derived from MODIS/CALIPO sensors
   integer :: ADV_TRIGGER    = 3   != 1 => Kain (2004),  2 => moisture adv trigger (Ma & Tan, 2009, Atmos Res), 3 => dcape trigger
   real    :: dcape_threshold= 70. != CAPE time rate threshold for ADV_TRIGGER = 3 (J kg^-1 hr^-1)
                                   != typical range is [-200,200] J/kg/hr, Wu et all (2007) recomends ~ 70 J/kg/hr
                                   != 55 J/kg/hr is indicated for the Amazon basin (Song&Zhang 2017)
   integer :: LCL_TRIGGER    = 1   != greater than zero, activates the LCL trigger which requires the lcl height
                                   != be lower than the pbl height, only for shallow convection
   
   integer :: EVAP_FIX       = 1   != fix total evap > column rainfall
   integer :: OUTPUT_SOUND   = 0   != outputs a "GEOS" vertical profile for the GF stand alone model

   real    :: tau_ocea_cp    = 2.*3600. != 
   real    :: tau_land_cp    = 2.*3600. != 
   real    :: mx_buoy1       = (cp*5.0 + xlv*2.e-3)*0.025  !=   250.5 J/kg
   real    :: mx_buoy2       = (cp*10. + xlv*4.e-3)        != 20004.0 J/kg: temp exc=10 K, q deficit=4 g/kg (=> mx_buoy ~ 20 kJ/kg)

   real    :: use_cloud_dissipation = 0.   != to acccount for the cloud dissipation at the decayment phase
   integer :: use_gustiness         = 0    != not in use
   real    :: use_random_num        = 0.   != stochastic pertubation for the height of maximum Zu

   real    :: beta_sh               = 2.2  != only for shallow plume
   integer :: use_linear_subcl_mf   = 1    != only for shallow plume
   real    :: cap_maxs              = 50.  != max distance (hPa) the air parcel is allowed to go up looking for the LFC
   integer :: liq_ice_number_conc   = 1    != include drop/ice number mixing ratio convective tendencies
   real    :: alpha_adv_tuning      = 0.8  != tuning parameter for the Becker et al (2021) closure
   real    :: col_sat_adv_threshold = 0.94 != suppress Qadv closure for col_sat_adv > col_sat_adv_threshold
   real    :: sig_factor            = 0.22 != exponential factor for the sigma determination (orig = 0.1)
   !------------------- internal variables  -------------------
   integer, parameter :: ON = 1, OFF = 0 !=  ON/OFF integer paremeters

   real    ::  HEI_DOWN_LAND     != [0.2,0.8] height of the max Z Downdraft , default = 0.50
   real    ::  HEI_DOWN_OCEAN    != [0.2,0.8] height of the max Z Downdraft , default = 0.50
   real    ::  HEI_UPDF_LAND     != [0.2,0.8] height of the max Z Updraft   , default = 0.35
   real    ::  HEI_UPDF_OCEAN    != [0.2,0.8] height of the max Z Updraft   , default = 0.35
   real    ::  MAX_EDT_LAND      != default= 0.9 - maximum evap fraction allowed over the land
   real    ::  MAX_EDT_OCEAN     != default= 0.9 - maximum evap fraction allowed over the ocean
   real    ::  FADJ_MASSFLX      != default= 1.0 - multiplicative factor for the mass flux at cloud base
   real    ::  T_star            != Scale Temperature for the DC closure
   real    ::  AVE_LAYER         != layer depth for average the properties of source air parcels (mbar)
   real    ::  c0                != autoconversion constant
   integer ::  USE_EXCESS        != default= 1   - use T,Q excess sub-grid scale variability

   !-- General internal controls for the diverse options in GF
   logical, parameter :: ENTRNEW        = .true.  != new entr formulation

   logical, parameter :: COUPL_MPHYSICS = .true.  != coupling with cloud microphysics (do not change to false)

   logical, parameter :: MELT_GLAC      = .true.  != turn ON/OFF ice phase/melting

   logical, parameter :: FEED_3DMODEL   = .true.  != set "false" to not feedback the AGCM with the
                                                  != heating/drying/transport conv tendencies
   logical            :: USE_C1D        = .false. != turn ON/OFF the 'c1d' detrainment approach, don't change this.

   logical            :: FIRST_GUESS_W  = .false. != use it to calculate a 1st guess of the updraft vert velocity

   integer, parameter :: aeroevap = 1             !=rainfall evaporation (1) orig  - (2) mix orig+new - (3) new

   integer, parameter ::               &
      maxens  = 1,  & ! 1  ensemble one on cap_max
      maxens2 = 1,  & ! 1  ensemble two on precip efficiency
      maxens3 = 16, & !16 ensemble three done in cup_forcing_ens16 for G3d
      ensdim  = maxens*maxens2*maxens3,&
      ens4    = 1
   !
  
   !- proportionality constant to estimate pressure
   !- gradient of updraft (Zhang and Wu, 2003, JAS) => REAL, PARAMETER ::    pgcon=-0.55
   real, parameter ::    pgcon= 0.0
   !- numerical constraints
   real, parameter ::       &
      xmbmaxshal  =  0.05,    &  ! kg/m2/s
      mintracer   =  tiny(1.),&  ! kg/kg - tiny(x)
      smallerQV   =  1.e-16  ,&  ! kg/kg
      PI = 3.1415926536

   integer, parameter :: MAX_NSPEC=200
   integer           ,dimension(MAX_NSPEC)    :: ind_chem
   character(len=100),dimension(MAX_NSPEC)    ::  CHEM_NAME
   integer           ,dimension(MAX_NSPEC)    ::  CHEM_NAME_MASK,CHEM_NAME_MASK_EVAP
   real              ,dimension(MAX_NSPEC)    ::  CHEM_ADJ_AUTOC
   integer :: ispc_CO
   type Hcts_vars
      real :: hstar,dhr,ak0,dak
   end type Hcts_vars
   type (Hcts_vars), allocatable :: Hcts(:)

   integer :: whoami_all, JCOL,itime1_in
   real    :: time_in

contains


   !---------------------------------------------------------------------------------------------------
   subroutine GF_GEOS5_INTERFACE(mxp,myp,mzp,mtp,ITRCR,LONS,LATS,DT_MOIST          &
      ,T, PLE, PLO, ZLE, ZLO, PK,  U, V, OMEGA , KH      &
      ,TH1, Q1, U1, V1 ,QLCN ,QICN,QLLS,QILS, CNPCPRATE  &
      ,CNV_MF0, CNV_PRC3, CNV_MFD, CNV_DQLDT ,ENTLAM     &
      ,CNV_MFC, CNV_UPDF, CNV_CVW, CNV_QC, CLCN,CLLS     &
      ,QV_DYN_IN,PLE_DYN_IN,U_DYN_IN,V_DYN_IN,T_DYN_IN   &
      ,RADSW   ,RADLW ,DQDT_BL  ,DTDT_BL                 &
      ,FRLAND  ,AREA  ,USTAR ,TSTAR ,QSTAR ,T2M ,Q2M     &
      ,TA      ,QA    ,SH    ,EVAP  ,PHIS                &
      ,KPBLIN                                            &
      ,STOCHASTIC_SIG, SIGMA_DEEP, SIGMA_MID             &
      ,DQDT_GF,DTDT_GF,MUPDP,MUPSH,MUPMD,MDNDP           &
      ,MFDP,MFSH,MFMD,ERRDP,ERRSH,ERRMD                  &
      ,AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_EC      &
      ,DTDTDYN,DQVDTDYN                                  &
      ,NCPL, NCPI, CNV_NICE, CNV_NDROP,CNV_FICE,CLDMICRO &
      ,TRACER,FSCAV,CNAMES,QNAMES,DTRDT_GF               &
      ,RSU_CN,REV_CN, PFI_CN, PFL_CN                     &
      ,TPWI,TPWI_star,LIGHTN_DENS                        &
      ,VAR3d_a,VAR3d_b,VAR3d_c,VAR3d_d                   &
       !,CNV_TR &!,VAR2d,ZKBCON
      )

      implicit none
      !INCLUDE "mpif.h"
      character(len=*),intent(in) :: CLDMICRO !set two-moment microphysics

      integer ,intent(in) :: mxp,myp,mzp,mtp,ITRCR

      real                             ,intent(in)   :: DT_moist

      real   ,dimension(mxp,myp)       ,intent(in)   :: FRLAND ,AREA ,USTAR ,TSTAR ,QSTAR &
         ,T2M ,Q2M ,TA ,QA ,SH ,EVAP ,PHIS  &
         ,KPBLIN,LONS,LATS &
         ,STOCHASTIC_SIG
      real   ,dimension(mxp,myp)       ,intent(out)  :: SIGMA_DEEP, SIGMA_MID

      real   ,dimension(mxp,myp,0:mzp) ,intent(in)   :: PLE,ZLE,PLE_DYN_IN
      real   ,dimension(mxp,myp,mzp)   ,intent(in)   :: T ,U ,V ,ZLO ,PLO ,PK ,OMEGA, KH       &
         ,RADSW  ,RADLW  ,DQDT_BL  ,DTDT_BL      &
         ,QV_DYN_IN,U_DYN_IN,V_DYN_IN,T_DYN_IN   &
         ,DTDTDYN,DQVDTDYN

      real   ,dimension(mxp,myp,mzp)   ,intent(inout):: TH1,Q1,U1,V1,QLCN ,QICN, NCPL, NCPI    &
         ,QLLS,QILS,CLLS
      real   ,dimension(mxp,myp,mzp)   ,intent(inout):: RSU_CN,REV_CN,VAR3d_a,VAR3d_b,VAR3d_c,VAR3d_d

      real   ,dimension(mxp,myp,0:mzp) ,intent(inout):: PFI_CN, PFL_CN

      !   REAL   ,DIMENSION(mxp,myp,mzp)   ,INTENT(INOUT):: CNV_TR
      real   ,dimension(mxp,myp,mzp)                 :: CNV_TR

      real   ,dimension(mxp,myp,0:mzp) ,intent(out)  :: CNV_MFC
      real   ,dimension(mxp,myp,mzp)   ,intent(out)  :: CNV_MF0 ,CNV_PRC3,CNV_MFD,CNV_DQLDT  &
         ,CNV_UPDF, CNV_CVW, CNV_QC,CLCN,ENTLAM&
         ,CNV_NICE, CNV_NDROP, CNV_FICE
      !-for debug purposes
      real   ,dimension(mxp,myp,mzp)   ,intent(inout):: DQDT_GF,DTDT_GF,MUPDP,MDNDP,MUPSH,MUPMD ,DTRDT_GF
      real   ,dimension(mxp,myp)       ,intent(inout):: MFDP,MFSH,MFMD,ERRDP,ERRSH,ERRMD
      real   ,dimension(mxp,myp)       ,intent(inout):: AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_EC &
         ,TPWI, TPWI_star

      real   ,dimension(MXP,MYP)       ,intent(out)  :: CNPCPRATE ,LIGHTN_DENS
      real   ,dimension(MXP,MYP)                     :: VAR2d,ZKBCON

      real   ,dimension(mxp,myp,mzp,itrcr) ,intent(inout)   :: TRACER  !=XHO in grid_moist_comp.f90
      real   ,dimension(itrcr)             ,intent(in   )   :: FSCAV
      character(len=*)  ,dimension(mtp)    ,intent(in   )   :: CNAMES,QNAMES

      real   ,dimension(mxp,myp,mzp)                        :: frct_liq,AA1_ADV,AA1_RADPBL

      integer      :: ims,ime, jms,jme, kms,kme,    &
         its,ite, jts,jte, kts,kte,    &
         mynum

      real,  dimension(mzp , mxp, myp ) :: up       &
         ,vp       &
         ,wp       &
         ,rvap     &
         ,temp     &
         ,press    &
         ,zm3d     &
         ,zt3d     &
         ,dm3d     &
         ,curr_rvap&
         ,buoy_exc &
         ,khloc



      real,  dimension(mzp , mxp, myp ) ::          &
           gsf_t   & ! grid-scale forcing for temp
         , gsf_q   & ! advection forcing for rv
         ,advf_t   & ! advection forcing for temp
         ,sgsf_t   & ! sub-grid scale forcing for temp
         ,sgsf_q   & ! sub-grid scale forcing for rv
         ,SRC_T    & ! temp tendency      from convection
         ,SRC_Q    & ! rv tendency        from convection
         ,SRC_CI   & ! cloud/ice tendency from convection
         ,SRC_U    & ! U tendency         from convection
         ,SRC_V    & ! V tendency         from convection
         ,SRC_NI   & ! Ice     number tendency from convection
         ,SRC_NL   & ! Droplet number tendency from convection
         ,SRC_BUOY & ! buoyancy tendency from downdrafts
         ,REVSU_GF & ! evaporation_or_sublimation of_convective_precipitation kg kg-1 s-1
         ,PRFIL_GF & ! ice_or_liq convective_precipitation flux: kg m2 s-1 (deep only)
         ,VAR3d_aGF& ! dummy 3-d var for output
         ,VAR3d_bGF& ! dummy 3-d var for output
         ,VAR3d_cGF& ! dummy 3-d var for output
         ,VAR3d_dGF& ! dummy 3-d var for output
         ,qexcp    & ! placeholder for Q   ex from cold pool param
         ,hexcp      ! placeholder for MSE ex from cold pool param


      real,  dimension(nmp, mzp , mxp, myp ) ::     &
          mp_ice   &
         ,mp_liq   &
         ,mp_cf

      real,  dimension(nmp, mzp , mxp, myp ) ::     &
          SUB_MPQI & ! subsidence transport applied to ice mix ratio
         ,SUB_MPQL & ! subsidence transport applied to cloud mix ratio
         ,SUB_MPCF   ! subsidence transport applied to cloud fraction


      real,  allocatable, dimension(:,:,:,:) :: SRC_CHEM ! tracer mixing ratio tendencies from the parameterized convection

      real,  dimension(mtp)       :: FSCAV_INT
      character(len=100)          :: AER_CHEM_MECH

      real,    dimension(mxp,myp) :: CONPRR
      real                        :: time=0.
      real,    dimension(mxp,myp) ::  aot500  ,temp2m  ,sfc_press &
         ,sflux_r ,sflux_t ,topt      &
         ,xland   ,dx2d    ,water_bud &
         ,col_sat, tke_pbl,rh_dicy_fct&
         ,wlpool
      integer, dimension(mxp,myp) :: kpbl,do_this_column
      integer, dimension(mzp)     :: flip
      integer :: k,i,j,iens,ispc, itime1=0
      !- for convective transport
      integer, dimension(mxp,myp,maxiens) ::  &
          ierr4d                       &
         ,jmin4d                       &
         ,klcl4d                       &
         ,k224d                        &
         ,kbcon4d                      &
         ,ktop4d                       &
         ,kstabi4d                     &
         ,kstabm4d

      real,dimension(mxp,myp,maxiens)     ::  &
          cprr4d                       &
         ,xmb4d                        &
         ,edt4d                        &
         ,pwav4d                       &
         ,sigma4d
      real,dimension(mxp,myp,mzp,maxiens) ::  &
          pcup5d                       &
         ,up_massentr5d                &
         ,up_massdetr5d                &
         ,dd_massentr5d                &
         ,dd_massdetr5d                &
         ,zup5d                        &
         ,zdn5d                        &
         ,prup5d                       &
         ,prdn5d                       &
         ,clwup5d                      &
         ,tup5d                        &
         ,conv_cld_fr5d

      !-----------local var in GEOS-5 data structure
      real,  dimension(mxp, myp, mzp) :: T1, ZLO_N,PLO_N,PK_N,TH_N,Q_N,T_N,U_N,V_N,DM
      real,  dimension(mxp,myp,0:mzp) :: ZLE_N,PLE_N
      integer :: status,alloc_stat ,wantrank=-99999
      real    :: tem1,dz,air_dens, src_cnvtr,snk_cnvtr,dz_int,tau_cp
      character(len=10) :: ENV_SETTING='DEFAULT'! ! 'CURRENT'/'DEFAULT'
      integer, parameter :: itest=1!3
      real :: RL, RI, disp_factor,x1,x2

      !--- to reproduce model behavior when using single-moment and version X0039_p5/f525_p5_fp of Dec 2019
      if(ZERO_DIFF == 1) then
         CUM_USE_EXCESS(:)  = 0
         FIRST_GUESS_W     = .false.
      endif

      !-- some initialization
      do_this_column =0
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
      CNV_NDROP      =0.
      CNV_NICE       =0.
      CNV_FICE       =0.
      SRC_NI         =0.
      SRC_NL         =0.
      SRC_T          =0.
      SRC_Q          =0.
      SRC_CI         =0.
      SRC_U          =0.
      SRC_V          =0.
      CNPCPRATE      =0.
      SUB_MPQI       =0.
      SUB_MPQL       =0.
      SUB_MPCF       =0.
      LIGHTN_DENS    =0.
      rh_dicy_fct    =0.
      SRC_BUOY       =0.
      REVSU_GF       =0.
      PRFIL_GF       =0.
      VAR3d_aGF      =0.
      VAR3d_bGF      =0.
      VAR3d_cGF      =0.
      VAR3d_dGF      =0.
      VAR2d          =0.
      !-
      !---temporary settings for debugging purposes
      !- special setting for SCM runs
      if(mxp==1 .and. myp==1 .and. maxval(T2m) < 1.e-6) return

      !- special setting for SCM runs
      if(mxp>1 .and. myp>1) wrtgrads = .false.
      !call mpi_comm_rank(MPI_COMM_WORLD,WHOAMI_ALL,status)
      !----
      if(use_gate) stop "use_gate must be false for GEOS5 runs"
      if( .not. use_gate .and. wrtgrads) call alloc_grads_arr(1,mzp,1,jl)
      !--------------------------------------------------------

      !- time counter
      ntimes = ntimes + 1
      !if(ntimes==1 .and. WHOAMI_ALL == 0) print *,'==> Applying GF convection scheme '
      mynum = -999
      call set_index_loops( ims,ime, jms,jme, kms,kme,    &
         its,ite, jts,jte, kts,kte,    &
         mxp,myp,mzp                   )

      !- define the vector "flip" to invert the z-axis orientation
      call flipz(flip,mzp)
      !
      if(.not.allocated(SRC_CHEM)) then
         allocate(SRC_CHEM(mtp, mzp, mxp, myp),stat=alloc_stat) !- tendency from convection
         if(alloc_stat==0) SRC_CHEM=0.0
      endif
      if(USE_TRACER_TRANSP==1) then
         AER_CHEM_MECH='GOCART' !in the future set as intent(in) from MoistGridComp
         call interface_aerchem(mtp,itrcr,aer_chem_mech, cnames,qnames, fscav, fscav_int)
      endif

      !
      !- 2-d input data
      aot500  (:,:) = 0.1  ! #
      !- as moist is called before surface, at the 1st time step all arrays
      !- from surface are zero
      if(maxval(T2m) < 1.e-6) then
         temp2m   (:,:) = T  (:,:,mzp) ! Kelvin
      else
         temp2m   (:,:) = T2M(:,:) ! or TA(:,:) ! Kelvin
      endif
      !- moisture flux from sfc
      sflux_r  (:,:) = EVAP(:,:) ! kg m-2 s-1
      !- sensible heat flux (sh) comes in W m-2, below it is converted to K m s-1
      !-(air_dens_sfc = ple(:,:,mzp)/( 287.04*TA(:,:)*(1.+0.608*QA(:,:)))))
      sflux_t  (:,:) = SH  (:,:) /(1004. * ple(:,:,mzp)/(287.04*T(:,:,mzp)*(1.+0.608*Q1(:,:,mzp)))) ! K m s-1
      !- topography height  (m)
      topt     (:,:) = PHIS(:,:)/MAPL_GRAV
      !- land/ocean fraction: land if < 1 ,ocean if = 1
      xland    (:,:) = 1.0-FRLAND(:,:)
      !
      !- grid length for the scale awareness
      dx2d(:,:) = sqrt(AREA(:,:)) ! meters
      !- special setting for SCM runs
      if(mxp==1 .and. myp==1) dx2d = 100000.

      !-pbl heigth index
      do j=1,myp
         do i=1,mxp
            if (nint(KPBLIN(i,j)) /= 0) then
               kpbl(i,j) = max(1, flip(min( nint(KPBLIN(i,j)), mzp)))
            else
               kpbl(i,j) = 1
            endif
            tke_pbl(i,j) = tkmin ! dummy 
         enddo
      enddo
      !
      !- 3-d input data
      !- current temperature T1 (after dyn+every physics process called before moist)
      T1 = PK * TH1
      !- any var with index "1" (and omega and pk) are already updated with dynamics
      !  tendencies and everything else (from physics) that was called before moist
      !
      if(trim(env_setting)=='CURRENT') then
         PLO_N = PLO
         T_N   = T1
         DM = ( PLE(:,:,1:mzp)-PLE(:,:,0:mzp-1) )*(1./MAPL_GRAV)
         !- 1st setting: enviromental state is the one already modified by dyn + physics
         do j=1,myp
            do i=1,mxp
               do k=1,mzp
                  temp  (k,i,j) = T1    (i,j,flip(k))
                  press (k,i,j) = PLO   (i,j,flip(k))*100. !Pa
                  rvap  (k,i,j) = Q1    (i,j,flip(k))!
                  up    (k,i,j) = U1    (i,j,flip(k))! already @ A-grid (m/s)
                  vp    (k,i,j) = V1    (i,j,flip(k))! already @ A-grid (m/s)
                  wp    (k,i,j) = OMEGA (i,j,flip(k))! Pa/s
                  zt3d  (k,i,j) = ZLO   (i,j,flip(k))! mid -layer level
                  zm3d  (k,i,j) = ZLE   (i,j,flip(k))! edge-layer level
                  dm3d  (k,i,j) = DM    (i,j,flip(k))
                  khloc (k,i,j) = KH    (i,j,flip(k))
                  curr_rvap(k,i,j) = Q1 (i,j,flip(k))!current rvap

                  mp_ice(lsmp,k,i,j) = QILS  (i,j,flip(k))
                  mp_liq(lsmp,k,i,j) = QLLS  (i,j,flip(k))
                  mp_cf (lsmp,k,i,j) = CLLS  (i,j,flip(k))
                  mp_ice(cnmp,k,i,j) = QICN  (i,j,flip(k))
                  mp_liq(cnmp,k,i,j) = QLCN  (i,j,flip(k))
                  mp_cf (cnmp,k,i,j) = CLCN  (i,j,flip(k))

               enddo
            enddo
         enddo
         !- sfc pressure (Pa)
         sfc_press(:,:) = PLE(:,:,mzp)
         !- Grid and sub-grid scale forcings for convection
         do j=1,myp
            do i=1,mxp
               do k=1,mzp
                  gsf_t (k,i,j) = 0.
                  gsf_q (k,i,j) = 0.
                  sgsf_t (k,i,j) = RADSW  (i,j,flip(k))+ RADLW(i,j,flip(k)) + DTDT_BL(i,j,flip(k))
                  sgsf_q (k,i,j) = DQDT_BL(i,j,flip(k))
                  advf_t (k,i,j) = 0.
               enddo
            enddo
         enddo
         !.. !- Tracer transport/scavenging section
         if(USE_TRACER_TRANSP==1) then
            where(TRACER<=0.0) TRACER = mintracer
         endif
      elseif(trim(env_setting)=='DEFAULT') then
         !-2nd setting: environmental state is that one before any tendency
         !- is applied (i.e, at begin of each time step).
         !- Get back the model state, heights and others variables at time N
         !- (or at the beggining of current time step)
         !- In physics, the state vars (T,U,V,PLE) are untouched and represent the
         !- model state after dynamics phase 1. But, "Q" is modified by physics, so
         !- depending on what was called before this subroutine, "Q" may be already
         !- changed from what it was just after dynamics phase 1. To solve this issue,
         !- "Q" just after dynamics is saved in the var named "QV_DYN_IN" in "GEOS_AgcmGridComp.F90".
         Q_N   =  QV_DYN_IN
         T_N   =   T_DYN_IN
         U_N   =   U_DYN_IN
         V_N   =   V_DYN_IN
         PLE_N = PLE_DYN_IN

         DM = ( PLE_N(:,:,1:mzp)-PLE_N(:,:,0:mzp-1) )*(1./MAPL_GRAV)
         !DM = ( PLE  (:,:,1:mzp)-PLE  (:,:,0:mzp-1) )*(1./MAPL_GRAV)
         call get_vars(mzp,mxp,myp,Q_N,T_N,PLE_N,ZLE_N,ZLO_N,PLO_N,PK_N,TH_N)
         !
         do j=1,myp
            do i=1,mxp
               do k=1,mzp
                  !
                  temp       (k,i,j) = T_N   (i,j,flip(k))! (K)
                  press      (k,i,j) = PLO_N (i,j,flip(k))*100.! (Pa) @ mid-layer level
                  rvap       (k,i,j) = Q_N   (i,j,flip(k))! water vapor mix ratio
                  up         (k,i,j) = U_N   (i,j,flip(k))! already @ A-grid (m/s)
                  vp         (k,i,j) = V_N   (i,j,flip(k))! already @ A-grid (m/s)
                  wp         (k,i,j) = OMEGA (i,j,flip(k))! (Pa/s)
                  zt3d       (k,i,j) = ZLO_N (i,j,flip(k))! mid -layer level (m)
                  zm3d       (k,i,j) = ZLE_N (i,j,flip(k))! edge-layer level (m)
                  dm3d       (k,i,j) = DM    (i,j,flip(k))
                  khloc      (k,i,j) = KH    (i,j,flip(k))
                  curr_rvap  (k,i,j) = Q1    (i,j,flip(k)) ! current rvap (dyn+phys)
                  mp_ice(lsmp,k,i,j) = QILS  (i,j,flip(k))
                  mp_liq(lsmp,k,i,j) = QLLS  (i,j,flip(k))
                  mp_cf (lsmp,k,i,j) = CLLS  (i,j,flip(k))
                  mp_ice(cnmp,k,i,j) = QICN  (i,j,flip(k))
                  mp_liq(cnmp,k,i,j) = QLCN  (i,j,flip(k))
                  mp_cf (cnmp,k,i,j) = CLCN  (i,j,flip(k))
               enddo
            enddo
         enddo
         !- sfc pressure (Pa)
         sfc_press(:,:) = PLE_N(:,:,mzp)
         !- Grid and sub-grid scale forcings for convection
         do j=1,myp
            do i=1,mxp
               do k=1,mzp
                  advf_t (k,i,j) = DTDTDYN (i,j,flip(k))
                  gsf_t  (k,i,j) = DTDTDYN (i,j,flip(k)) + RADLW(i,j,flip(k)) + RADSW  (i,j,flip(k))
                  gsf_q  (k,i,j) = DQVDTDYN(i,j,flip(k))

                  sgsf_t (k,i,j) = DTDT_BL (i,j,flip(k))
                  sgsf_q (k,i,j) = DQDT_BL (i,j,flip(k))
               enddo
            enddo
         enddo
         !
         !- Tracer transport/scavenging section
         if(USE_TRACER_TRANSP==1) then
            where(TRACER<=0.0) TRACER = mintracer
         endif
      else
         stop 'unknown env_setting at convpar_gf_geos5.F90'
      endif

      if(CONVECTION_TRACER==1) then
         do j=1,myp
            do i=1,mxp
               !--- saturation CWV
               col_sat(i,j) = TPWI(i,j)/(1.e-6+TPWI_star(i,j))
               col_sat(i,j) = min(1.,max(0.,col_sat(i,j)))
               !--- temporary to hold only CWV in mm
               ! col_sat(i,j) = TPWI(i,j)
               !
               do k=1,mzp
                  buoy_exc(k,i,j)=CNV_TR (i,j,flip(k))
               enddo
            enddo
         enddo
         qexcp = 0. 
         hexcp = 0. 
         wlpool= 0.
        !print*,"buoy_exc1=",maxval(buoy_exc),minval(buoy_exc)
      else
         buoy_exc = 0.0
      endif

      !- call the driver routine to apply the parameterization
      call GF_GEOS5_DRV(mxp,myp,mzp,mtp,nmp, time, itime1   &
         ,ims,ime, jms,jme, kms,kme   &
         ,its,ite, jts,jte, kts,kte   &
         ,flip        &
         ,fscav_int   &
         ,mynum       &
         ,dt_moist    &
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
         ,xland       &
         ,sfc_press   &
         ,kpbl        &
         ,tke_pbl     &
         !--- atmos state
         ,col_sat     &
         ,up          &
         ,vp          &
         ,wp          &
         ,temp        &
         ,press       &
         ,rvap        &
         ,mp_ice      &
         ,mp_liq      &
         ,mp_cf       &
         ,curr_rvap   &
         !--- atmos composition state
         ,TRACER      & !- note: uses GEOS-5 data structure
         !---- forcings---
         ,buoy_exc    &
         ,gsf_t       &
         ,gsf_q       &
         ,advf_t      &
         ,sgsf_t      &
         ,sgsf_q      &
         !---- output ----
         ,conprr      &
         ,lightn_dens &
         ,rh_dicy_fct &
         ,src_t       &
         ,src_q       &
         ,src_ci      &
         ,src_nl      &
         ,src_ni      &
         ,src_u       &
         ,src_v       &
         ,sub_mpqi    &
         ,sub_mpql    &
         ,sub_mpcf    &
         ,src_buoy    &
         ,src_chem    &
         ,revsu_gf    &
         ,prfil_gf    &
         !
         !
         ,do_this_column&
         ,ierr4d       &
         ,jmin4d       &
         ,klcl4d       &
         ,k224d        &
         ,kbcon4d      &
         ,ktop4d       &
         ,kstabi4d     &
         ,kstabm4d     &
         ,cprr4d       &
         ,xmb4d        &
         ,edt4d        &
         ,pwav4d       &
         ,sigma4d      &
         ,pcup5d       &
         ,up_massentr5d&
         ,up_massdetr5d&
         ,dd_massentr5d&
         ,dd_massdetr5d&
         ,zup5d        &
         ,zdn5d        &
         ,prup5d       &
         ,prdn5d       &
         ,clwup5d      &
         ,tup5d        &
         ,conv_cld_fr5d&
         !-- for debug/diagnostic
         ,AA0,AA1,AA1_ADV,AA1_RADPBL,AA1_BL,AA2,AA3,AA1_CIN,TAU_BL,TAU_EC &
         ,VAR2d,VAR3d_aGF,VAR3d_bGF,VAR3d_cGF,VAR3d_dGF)


      if(FEED_3DMODEL)then
         !--- vertical fraction of liq/ice water phases
         if(FRAC_MODIS==1 .and. icumulus_gf(DEEP)==ON) then
            do j=1,myp
               do i=1,mxp
                  do k=1,mzp
                     frct_liq(i,j,k) = fract_liq_f(tup5d(i,j,flip(k),DEEP))
                  enddo
               enddo
            enddo
         else
            do j=1,myp
               do i=1,mxp
                  do k=1,mzp
                     frct_liq(i,j,k) = fract_liq_f(T(i,j,k))
                  enddo
               enddo
            enddo
         endif
         !-- update GEOS-5 model state with the feedback from cumulus convection
         !- to include the tendencies from the convection,  update the vars th1,q1,v1 and u1
         do j=1,myp
            do i=1,mxp
               if(do_this_column(i,j) == 0) cycle
               !- conv precip rate: mm/s = kg m-2 s-1
               CNPCPRATE (i,j) =  CONPRR (i,j)

               if(ITEST==0) CNPCPRATE(i,j) =  0.

               !
               do k=1,mzp ! in the future, limit the vertical loop to ktop (DO k=mzp,flip(ktop),-1)?
                  !- convert from d temp/dt to d theta/dt using PK => d theta/dt = (1/pk)*d temp/dt
                  !- (think if PK must be the current one _or_ at the begin of the time step
                  TH1(i,j,k) = TH1(i,j,k) + DT_moist * SRC_T(flip(k),i,j) / PK(i,j,k)

                  Q1 (i,j,k) = Q1 (i,j,k) + DT_moist * SRC_Q(flip(k),i,j)
                  Q1 (i,j,k) = max(smallerQV, Q1 (i,j,k))
                  !
                  !- simple splitting of condensate tendency into liq/ice phases
                  !- these are 'anvil' mixing ratio and not 'grid-scale' mix ratio
                  !- (the convective source will be applied in progno_cloud, routine "consrc")
                  ! QLCN (i,j,k) = QLCN (i,j,k) + DT_moist * SRC_CI(flip(k),i,j) * frct_liq(i,j,k)
                  ! QICN (i,j,k) = QICN (i,j,k) + DT_moist * SRC_CI(flip(k),i,j) * (1.0-frct_liq(i,j,k))

                  if(ITEST==3) then
                     !- simple splitting of condensate tendency into liq/ice phases
                     !- these are 'grid-scale' mix ratio
                     !- (the convective source will be set to zero, see below)
                     QLLS (i,j,k) = QLLS (i,j,k) + DT_moist * SRC_CI(flip(k),i,j) *  frct_liq(i,j,k)
                     QILS (i,j,k) = QILS (i,j,k) + DT_moist * SRC_CI(flip(k),i,j) *  (1.-frct_liq(i,j,k))
                  endif
               enddo
               !--- sublimation/evaporation tendencies (kg/kg/s)
               do k=1,mzp
                  !--- sublimation/evaporation tendencies (kg/kg/s)
                  RSU_CN (i,j,k) = REVSU_GF(flip(k),i,j)* (1.0-frct_liq(i,j,k))
                  REV_CN (i,j,k) = REVSU_GF(flip(k),i,j)*  frct_liq(i,j,k)
                  !--- preciptation fluxes (kg/kg/s)
                  PFI_CN (i,j,k) = PRFIL_GF(flip(k),i,j)* (1.0-frct_liq(i,j,k))
                  PFL_CN (i,j,k) = PRFIL_GF(flip(k),i,j)*  frct_liq(i,j,k)
               enddo

            enddo
         enddo
         !-----
         if(USE_MOMENTUM_TRANSP > 0) then
            do j=1,myp
               do i=1,mxp
                  if(do_this_column(i,j) == 0) cycle
                  do k=1,mzp
                     U1 (i,j,k) = U1 (i,j,k) + DT_moist * SRC_U(flip(k),i,j)
                     V1 (i,j,k) = V1 (i,j,k) + DT_moist * SRC_V(flip(k),i,j)
                  enddo
               enddo
            enddo
         endif


         if(APPLY_SUB_MP == 1) then
            do j=1,myp
               do i=1,mxp
                  if(do_this_column(i,j) == 0) cycle
                  do k=1,mzp
                     QLLS (i,j,k) = QLLS (i,j,k) + DT_moist * SUB_MPQL(LSMP,flip(k),i,j)
                     QILS (i,j,k) = QILS (i,j,k) + DT_moist * SUB_MPQI(LSMP,flip(k),i,j)
                     CLLS (i,j,k) = CLLS (i,j,k) + DT_moist * SUB_MPCF(LSMP,flip(k),i,j)
                     QLCN (i,j,k) = QLCN (i,j,k) + DT_moist * SUB_MPQL(CNMP,flip(k),i,j)
                     QICN (i,j,k) = QICN (i,j,k) + DT_moist * SUB_MPQI(CNMP,flip(k),i,j)
                     CLCN (i,j,k) = CLCN (i,j,k) + DT_moist * SUB_MPCF(CNMP,flip(k),i,j)
                  enddo
               enddo
            enddo
         endif

         if(USE_TRACER_TRANSP==1) then

            do j=1,myp
               do i=1,mxp
                  if(do_this_column(i,j) == 0) cycle
                  do k=1,mzp
                     !-special array for output of CO tendency
                     DTRDT_GF(i,j,k)=SRC_CHEM(ispc_CO,flip(k),i,j)

                     !- update tracer mass mixing ratios
                     do ispc=1,mtp

                        TRACER(i,j,k,ispc)=TRACER(i,j,k,ispc)+ DT_moist * SRC_CHEM(ispc,flip(k),i,j) * CHEM_NAME_MASK(ispc)!

                        !-- final check for negative tracer mass mixing ratio
                        TRACER(i,j,k,ispc)=max(mintracer, TRACER(i,j,k,ispc))
                     enddo
                  enddo
               enddo
            enddo
           !-- final check for negative tracer mass mixing ratio
           !where (TRACER < mintracer) TRACER = mintracer

         endif

         !--- for lightning flash rate
         !--- ice/liq precip fluxes
         if(icumulus_gf(deep) == ON) then
            do j=1,myp
               do i=1,mxp
                  if(ierr4d(i,j,deep) .ne. 0) cycle
                  ZKBCON(i,j) = ZLE(i,j,flip(kbcon4d(i,j,DEEP)))
                  do k=mzp,1,-1
                     PFI_CN (i,j,k) = prfil_gf(flip(k),i,j)* (1.0-frct_liq(i,j,k))
                     PFL_CN (i,j,k) = prfil_gf(flip(k),i,j)*  frct_liq(i,j,k)
                  enddo
               enddo
            enddo
         endif

         !--- for dummy output
         do j=1,myp
            do i=1,mxp
               do k=1,mzp
                  VAR3d_a  (i,j,k) = VAR3d_aGF(flip(k),i,j)
                  VAR3d_b  (i,j,k) = VAR3d_bGF(flip(k),i,j)
                  VAR3d_c  (i,j,k) = VAR3d_cGF(flip(k),i,j)
                  VAR3d_d  (i,j,k) = VAR3d_dGF(flip(k),i,j)
               enddo
            enddo
         enddo


         do IENS=1, maxiens
            if(icumulus_gf(IENS) == ON) then
               do j=1,myp
                  do i=1,mxp
                     if(ierr4d(i,j,IENS) .ne. 0) cycle
                     do k=mzp,flip(ktop4d(i,j,IENS))-1,-1

                        DZ       = -( ZLE(i,j,k) - ZLE(i,j,k-1) )
                        air_dens = 100.*plo_n(i,j,k)/(287.04*T_n(i,j,k)*(1.+0.608*Q_n(i,j,k)))

                        !- special treatment for CNV_DQLDT: 'convective_condensate_source',  UNITS     = 'kg m-2 s-1',
                        !- SRC_CI contains contributions from deep, shallow,... . So, do not need to be accumulated over  CNV_DQLDT
                        !- note the SRC_CI has different array structure (k,i,j) _not_ (i,j,k)
                        if(ITEST /= 3) &
                           CNV_DQLDT(i,j,k)=  SRC_CI(flip(k),i,j)*DZ * air_dens !units: kg[w]/(kg[air] s) * m * kg[air]/m^3 = kg[w]/(m^2 s)

                        !CNV_DQLDT(i,j,k)=  SRC_CI(flip(k),i,j)*DM(i,j,k)     !units: kg[w]/(kg[air] s) * m * kg[air]/m^3 = kg[w]/(m^2 s)
                        !
                        !-'detraining_mass_flux',UNITS     = 'kg m-2 s-1',
                        CNV_MFD (i,j,k) = CNV_MFD (i,j,k) + ( up_massdetr5d(i,j,flip(k),IENS) )


                        !-'cloud_base_mass_flux',    units = 'kg m-2 s-1',                                  &
                        CNV_MF0 (i,j,k) = CNV_MF0 (i,j,k) + zup5d(i,j,flip(k),IENS)

                        !-convective mass flux [kg m^{-2} s^{-1}] - with downdrafts
                        ! CNV_MFC (i,j,k) = CNV_MFC (i,j,k) +  ( zup5d(i,j,flip(k),IENS) + edt4d(i,j,IENS)*  zdn5d(i,j,flip(k),IENS))

                        !---only updraft
                        CNV_MFC (i,j,k) = CNV_MFC (i,j,k) + zup5d(i,j,flip(k),IENS)

                        if(zup5d(i,j,flip(k),IENS) > 1.0e-6) then
                           !-'entrainment parameter',  UNITS     ='m-1',
                           ENTLAM  (i,j,k) =  ENTLAM   (i,j,k) + (up_massentr5d(i,j,flip(k),IENS)/(DZ*zup5d(i,j,flip(k),IENS)))

                           !-'updraft_vertical_velocity',            UNITS     = 'hPa s-1',
                           CNV_CVW (i,j,k) = -0.2 ! hPa/s =>  4 m/s
                        endif

                        !-'grid_mean_convective_condensate', UNITS     ='kg kg-1'
                        CNV_QC  (i,j,k) = CNV_QC  (i,j,k) + clwup5d(i,j,flip(k),IENS)
                        !
                        !
                        !~ !----------------------------------------------------------------------------------------------------
                        !- not using progno-cloud to calculate the precip from the convective column
                        !- if CNV_PRC3 will be send to progno-cloud, set CNPCPRATE = zero
                        !-'convective_precipitation_from_GF',UNITS     = 'kg m-2 s-1',
                        !- JAN/17/2017 : the units above are wrong. The correct are kg[precip water]/kg[air]
                        CNV_PRC3(i,j,k) = CNV_PRC3(i,j,k) +                 (prup5d(i,j,flip(k),IENS) + &
                           edt4d(i,j,IENS)* prdn5d(i,j,flip(k),IENS) ) &
                           * DT_moist/(DZ*air_dens)

                        !-'updraft_areal_fraction',
                        if(zup5d(i,j,flip(k),IENS) > 1.0e-6) CNV_UPDF(i,j,k) = 0.033
                        !----------------------------------------------------------------------------------------------------
                        if(ITEST==2) then
                           !-'updraft_areal_fraction',
                           if(zup5d(i,j,flip(k),IENS) > 1.0e-6) then
                              CNV_UPDF(i,j,k) = 0.033
                           else
                              CNV_UPDF(i,j,k) = 0.0
                           endif

                           !-'convective_cloud_area_fraction', adimensional
                           !- Tiedtke formulation
                           CLCN(i,j,k) = CLCN(i,j,k) + (1.0-CLCN(i,j,k))*(up_massdetr5d(i,j,flip(k),IENS) &
                              * DT_moist/(DZ*air_dens)+CNV_UPDF(i,j,k))

                           !- Chab&Bechtold 2002/2005 formulation
                           ! CLCN(i,j,k) = CLCN(i,j,k) + (1.0-CLCN(i,j,k))*(conv_cld_fr5d(i,j,flip(k),IENS)+CNV_UPDF(i,j,k))
                           !
                           CLCN(i,j,k) = max(0.,min(CLCN(i,j,k),0.99))
                        endif
                     !----------------------------------------------------------------------------------------------------
                     enddo
                   !print*,"iens=",iens,maxval(conv_cld_fr5d(i,j,:,IENS)),minval(conv_cld_fr5d(i,j,:,IENS));call flush(6)
                  enddo
               enddo
            endif
         enddo
      endif

      if(adjustl(CLDMICRO) =="2MOMENT") then
         !- we adjust convective cloud condensate and number here
         do j=1,myp
            do i=1,mxp
               do k=1,mzp

                  !---obsolete
                  !tem1 = T(i,j,k)
                  !RL =   10.0  + (12.0*(283.0- tem1)/40.0)
                  !RL =   min(max(RL, 10.0), 30.0)*1.e-6
                  !RI =   100.0 + (80.0*(tem1- 253.0)/40.0)
                  !RI =   min(max(RI, 20.0), 250.0)*1.e-6
                  !tem1 = 1.- (tem1 - 235.0) /38.0
                  !tem1 =  min(max(0.0, tem1), 1.0)
                  ! make up some "number" sources. In the future this should depend explicitly on the convective mphysics
                  !disp_factor =  10.0 ! used to account somehow for the size dist
                  !SRC_NL(flip(k),i,j) = SRC_CI(flip(k),i,j)* (1.0-tem1) /(1.333 * MAPL_PI*RL*RL*RL*997.0*disp_factor)
                  !SRC_NI(flip(k),i,j)= SRC_CI(flip(k),i,j) * tem1 /(1.333 * MAPL_PI *RI*RI*RI*500.0*disp_factor)
                  !CNV_FICE (i, j, k)   =   tem1
                  !---obsolete

                  tem1 = frct_liq(i,j,k)

                  QLCN (i,j,k) = QLCN (i,j,k) + DT_moist * SRC_CI(flip(k),i,j) * (1.0-tem1)
                  QICN (i,j,k) = QICN (i,j,k) + DT_moist * SRC_CI(flip(k),i,j) * tem1

                  NCPL (i,j,k) = NCPL (i,j,k) + DT_moist * SRC_NL(flip(k),i,j)
                  NCPI (i,j,k) = NCPI (i,j,k) + DT_moist * SRC_NI(flip(k),i,j)

                  DZ       = -( ZLE(i,j,k) - ZLE(i,j,k-1) )
                  air_dens = 100.*PLO_n(i,j,k)/(287.04*T_n(i,j,k)*(1.+0.608*Q_n(i,j,k)))

                  if (CNV_MFD (i,j,k)  .gt. 0.) then
                     CNV_NICE (i,j,k)=   SRC_NI(flip(k),i,j)*DZ * air_dens/CNV_MFD (i,j,k)
                     CNV_NDROP(i,j,k)=   SRC_NL(flip(k),i,j)*DZ * air_dens/CNV_MFD (i,j,k)
                  endif
               enddo
            enddo
         enddo
         !--- special section for convective cloud fraction
         do iens=1,maxiens
            if(icumulus_gf(iens) .ne. ON) cycle
            do j=1,myp
               do i=1,mxp
                  if(ierr4d(i,j,IENS) .ne. 0) cycle
                  do k=mzp,flip(ktop4d(i,j,IENS))-1,-1

                     DZ       = -( ZLE(i,j,k) - ZLE(i,j,k-1) )
                     air_dens = 100.*plo_n(i,j,k)/(287.04*T_n(i,j,k)*(1.+0.608*Q_n(i,j,k)))
                     CLCN(i,j,k) = CLCN(i,j,k) + (1.0-CLCN(i,j,k))*(up_massdetr5d(i,j,flip(k),iens) &
                        * dt_moist/(dz*air_dens)+CNV_UPDF(i,j,k))
                     CLCN(i,j,k) = max(0.,min(CLCN(i,j,k),0.99))
                  enddo
               enddo
            enddo
         enddo
      endif !2 moment
      !
      !--- cold pool/"convection tracer"
      if(CONVECTION_TRACER==1) then
         DTRDT_GF=0. !temporary   for output only
         do j=1,myp
            do i=1,mxp

               tau_cp=FRLAND(i,j)*tau_land_cp + (1.0-FRLAND(i,j))*tau_ocea_cp

               do k=1,mzp

                  !- sink term (exp decay 1h)
                  snk_cnvtr =  dt_moist * abs(CNV_TR (i,j,k))/tau_cp

                  !DZ          = -( ZLE(i,j,k) - ZLE(i,j,k-1) )
                  !air_dens = 100.*PLO_N(i,j,k)/(287.04*T_n(i,j,k)*(1.+0.608*Q_n(i,j,k)))
                  !
                  !- downdraft convective mass flux [kg m^{-2} s^{-1}]
                  ! iens =?
                  !src_cnvtr =  edt4d(i,j,iens)*zdn5d(i,j,flip(k),iens)

                  !- downdraft detrainment mass flux [kg m^{-2} s^{-1}]
                  ! iens =?
                  !src_cnvtr =  edt4d(i,j,iens)*dd_massdetr5d(i,j,flip(k),iens)


                  !- source term
                  !- downdraft detrainment of buoyancy [ J/kg s^{-1}]
                  !- negative sign => source for updraft lifting
                  src_cnvtr = - dt_moist * min(0.,SRC_BUOY(flip(k),i,j))

                  !- 'continuity' equation = ADV + SRC - SINK
                  CNV_TR (i,j,k) = CNV_TR (i,j,k)  + src_cnvtr  - snk_cnvtr

                  !temporary for output only
                  DTRDT_GF(i,j,k)=SRC_BUOY(flip(k),i,j)

               enddo
            enddo
         enddo
        !print*,"buoy_exc2=",maxval(SRC_BUOY),minval(SRC_BUOY)
      endif

      !
      if(maxval(icumulus_gf)>0) then
         do IENS=1, maxiens
            if(icumulus_gf(IENS) == ON .and. IENS== DEEP) then
               do j=1,myp
                  do i=1,mxp
                     if(ierr4d(i,j,DEEP) /= 0) cycle
                     MFDP      (i,j)      =   xmb4d(i,j,DEEP)
                     SIGMA_DEEP(i,j)      = sigma4d(i,j,deep)
                     MUPDP     (i,j,1:mzp)=zup5d(i,j,flip(1):flip(mzp):-1,DEEP)
                     MDNDP     (i,j,1:mzp)=zdn5d(i,j,flip(1):flip(mzp):-1,DEEP) * edt4d(i,j,IENS)
                  enddo
               enddo
            elseif(icumulus_gf(IENS) == ON .and. IENS== SHAL) then
               do j=1,myp
                  do i=1,mxp
                     if(ierr4d(i,j,SHAL) /= 0) cycle
                     MFSH (i,j)      =xmb4d(i,j,SHAL)
                     MUPSH(i,j,1:mzp)=zup5d(i,j,flip(1):flip(mzp):-1,SHAL)
                  enddo
               enddo
            elseif(icumulus_gf(IENS) == ON .and. IENS== MID) then
               do j=1,myp
                  do i=1,mxp
                     if(ierr4d(i,j,MID) /= 0) cycle
                     MFMD      (i,j)      = cprr4d (i,j,MID) ! xmb4d(i,j,MID)temporary saving for mid precip
                     SIGMA_MID (i,j)      = sigma4d(i,j,MID )
                     MUPMD     (i,j,1:mzp)=zup5d(i,j,flip(1):flip(mzp):-1,MID)
                  enddo
               enddo
            endif
         enddo
         !- for output purposes, ierr=0 (convection is ON) will be changed to 1
         where(ierr4d==0)ierr4d=1
         where(ierr4d>1)ierr4d=0
         do j=1,myp
            do i=1,mxp
               DQDT_GF(i,j,1:mzp)=SRC_Q(flip(1):flip(mzp):-1,i,j)!note SQR_Q is (k,i,j)
               DTDT_GF(i,j,1:mzp)=SRC_T(flip(1):flip(mzp):-1,i,j)!note SQR_T is (k,i,j)

               ERRDP(i,j)=float(ierr4d(i,j,DEEP))
               ERRSH(i,j)=float(ierr4d(i,j,SHAL))
               ERRMD(i,j)=float(ierr4d(i,j,MID ))
            enddo
         enddo
      endif
      !
      if( allocated(src_chem))  deallocate(src_chem,stat=alloc_stat) !tendency   from convection

      !- for debugging purposes only
      if(.not. use_gate .and. wrtgrads) call alloc_grads_arr(1,mzp,2,jl)

   end subroutine GF_GEOS5_INTERFACE
   !---------------------------------------------------------------------------------------------------

   subroutine GF_GEOS5_DRV(mxp,myp,mzp,mtp,nmp, time, itime1 &
      ,ims,ime, jms,jme, kms,kme               &
      ,its,ite, jts,jte, kts,kte               &
      ,flip                                    &
      ,FSCAV                                   &
      ,mynum                 &
      ,dt                    &
      ,dx2d                  &
      ,stochastic_sig        &
      ,zm                    &
      ,zt                    &
      ,dm                    &

      ,lons                  &
      ,lats                  &
      ,aot500                &
      ,temp2m                &
      ,sflux_r               &
      ,sflux_t               &
      ,qexcp                 &
      ,hexcp                 &
      ,wlpool                &      
      ,topt                  &
      ,xland                 &
      ,sfc_press             &
      ,kpbl                  &
      ,tke_pbl               &

      ,col_sat               &
      ,u                     &
      ,v                     &
      ,w                     &
      ,temp                  &
      ,press                 &
      ,rvap                  &
      ,mp_ice                &
      ,mp_liq                &
      ,mp_cf                 &
      ,curr_rvap             &
      ,TRACER                &!-note: uses GEOS-5 data structure

      !---- forcings---
      ,buoy_exc              &
      ,rthften               &! gsf_t
      ,rqvften               &! gsf_q
      ,rth_advten            &!advf_t
      ,rthblten              &!sgsf_t
      ,rqvblten              &!sgsf_q
      !---- output ----
      ,conprr                &
      ,lightn_dens           &
      ,rh_dicycle_fct        &
      ,rthcuten              &
      ,rqvcuten              &
      ,rqccuten              &
      ,rnlcuten              &
      ,rnicuten              &
      ,rucuten               &
      ,rvcuten               &
      ,sub_mpqi              &
      ,sub_mpql              &
      ,sub_mpcf              &
      ,rbuoycuten            &
      ,rchemcuten            &
      ,revsu_gf              &
      ,prfil_gf              &
      !
      ,do_this_column        &
      ,ierr4d                &
      ,jmin4d                &
      ,klcl4d                &
      ,k224d                 &
      ,kbcon4d               &
      ,ktop4d                &
      ,kstabi4d              &
      ,kstabm4d              &
      ,cprr4d                &
      ,xmb4d                 &
      ,edt4d                 &
      ,pwav4d                &
      ,sigma4d               &
      ,pcup5d                &
      ,up_massentr5d         &
      ,up_massdetr5d         &
      ,dd_massentr5d         &
      ,dd_massdetr5d         &
      ,zup5d                 &
      ,zdn5d                 &
      ,prup5d                &
      ,prdn5d                &
      ,clwup5d               &
      ,tup5d                 &
      ,conv_cld_fr5d         &
      !-- for debug/diagnostic
      ,AA0,AA1,AA1_ADV,AA1_RADPBL,AA1_BL,AA2,AA3,AA1_CIN,TAU_BL,TAU_EC &
      ,VAR2d,VAR3d_aGF,VAR3d_bGF,VAR3d_cGF,VAR3d_dGF)

      implicit none
      !include "mpif.h"
      !------------------------------------------------------------------------
      integer, intent(in) :: ims,ime, jms,jme, kms,kme,    &
         its,ite, jts,jte, kts,kte,    &
         mynum,mzp,mxp,myp,mtp,nmp, itime1

      real,    intent(in) :: DT, time

      integer, intent(in), dimension(mzp) :: flip

      real:: FSCAV(mtp)

      real,    dimension(kts:kte,its:ite,jts:jte), intent(in)  ::       &
         zm,        &
         zt,        &
         u,         &
         v,         &
         w,         &
         rvap,      &
         temp,      &
         press,     &
         dm,        &
         curr_rvap, &
         buoy_exc,  &
         qexcp,     &
         hexcp 

      integer, dimension(its:ite,jts:jte), intent(in) :: kpbl

      !-- intent (in)
      real,    dimension(its:ite,jts:jte) :: topt ,aot500 ,temp2m ,sfc_press &
         ,sflux_r ,sflux_t,xland,lons,lats,dx2d,col_sat,stochastic_sig,tke_pbl,wlpool

      real,    dimension(kts:kte,its:ite,jts:jte), intent(in) :: rthften    &
         ,rqvften    &
         ,rth_advten &
         ,rthblten   &
         ,rqvblten

      real,    dimension(its:ite,jts:jte),         intent(out  ) ::   CONPRR,LIGHTN_DENS
      real,    dimension(its:ite,jts:jte),         intent(inout) ::   rh_dicycle_fct

      real,    dimension(kts:kte,its:ite,jts:jte), intent(out) :: &
          rthcuten   &
         ,rqvcuten   &
         ,rqccuten   &
         ,rnlcuten   &
         ,rnicuten   &
         ,rucuten    &
         ,rvcuten    &
         ,rbuoycuten &
         ,revsu_gf   &
         ,prfil_gf   &
         ,var3d_agf  &
         ,var3d_bgf  &
         ,var3d_cgf  &
         ,var3d_dgf


      real,    dimension(nmp,kts:kte,its:ite,jts:jte), intent(in)  :: &
          mp_ice  &
         ,mp_liq  &
         ,mp_cf

      real,    dimension(nmp,kts:kte,its:ite,jts:jte), intent(out) :: &
          sub_mpqi   &
         ,sub_mpql   &
         ,sub_mpcf

      !-***** TRACER has different data structure   (i,j,k,ispc) *********
      real,    dimension(its:ite,jts:jte,kts:kte,mtp), intent(in )  :: TRACER
      !-***** rchemcuten uses the GF data structure (ispc,k,i,j) *********
      real,    dimension(mtp,kts:kte,its:ite,jts:jte), intent(out)  :: rchemcuten

      integer, dimension(its:ite,jts:jte), intent(inout) :: do_this_column

      !- for convective transport and cloud/radiation (OUT)
      integer,dimension(mxp,myp,maxiens)::  &
         !  integer, dimension(its:ite,jts:jte,maxiens) , INTENT(OUT) ::    &
          ierr4d                    &
         ,jmin4d                    &
         ,klcl4d                    &
         ,k224d                     &
         ,kbcon4d                   &
         ,ktop4d                    &
         ,kstabi4d                  &
         ,kstabm4d

      real,dimension(mxp,myp,maxiens)::     &
         !   real,dimension(its:ite,jts:jte,maxiens)    , INTENT(OUT) ::    &
          cprr4d                    &
         ,xmb4d                     &
         ,edt4d                     &
         ,pwav4d                    &
         ,sigma4d
      real,dimension(mxp,myp,mzp,maxiens):: &
         !   real,dimension(its:ite,jts:jte,kts:kte,maxiens), INTENT(OUT) ::    &
          pcup5d                    &
         ,up_massentr5d             &
         ,up_massdetr5d             &
         ,dd_massentr5d             &
         ,dd_massdetr5d             &
         ,zup5d                     &
         ,zdn5d                     &
         ,prup5d                    &
         ,prdn5d                    &
         ,clwup5d                   &
         ,tup5d                     &
         ,conv_cld_fr5d
      !--for debug
      real   ,dimension(mxp,myp)  ,intent(inout)  :: AA0,AA1,AA1_ADV,AA1_RADPBL &
                                                    ,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_EC,var2d

      !----------------------------------------------------------------------
      ! LOCAL VARS

      ! basic environmental input includes
      ! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
      ! convection for this call only and at that particular gridpoint
      !

      real,   dimension (kts:kte,its:ite,jts:jte)  :: Tpert_h,Tpert_v

      real,   dimension (its:ite,jts:jte) ::  rtgt

      real,   dimension (its:ite,kts:kte) ::               &
          zo,temp_old,qv_old,PO,US,VS,rhoi,phil            &
         ,temp_new_dp,qv_new_dp,temp_new_sh,qv_new_sh,z2d  &
         ,tkeg,rcpg,dhdt,temp_new_md,qv_new_md             &
         ,temp_new_BL,qv_new_BL,dm2d,temp_tendqv,qv_curr   &
         ,buoy_exc2d,revsu_gf_2d,prfil_gf_2d,var3d_Agf_2d  &
         ,var3d_Bgf_2d,temp_new,qv_new,Tpert_2d            &
         ,temp_new_adv,qv_new_adv


      real,   dimension (its:ite,kts:kte,maxiens) ::       &
         outt,outq,outqc,outu,outv,outbuoy,outnliq,outnice

      real,   dimension (mtp,its:ite,kts:kte)         :: se_chem
      real,   dimension (mtp,its:ite,kts:kte,maxiens) :: out_chem

      real,   dimension (nmp,its:ite,kts:kte)         :: mpqi,mpql,mpcf
      real,   dimension (nmp,its:ite,kts:kte,maxiens) :: outmpqi,outmpql,outmpcf

      real,   dimension (its:ite)   :: ter11, xlandi,pbl,zws,ccn,psur &
         ,ztexec,zqexec,h_sfc_flux,le_sfc_flux,tsur&
         ,xlons,xlats,fixout_qv,cum_ztexec,cum_zqexec


      real,   dimension (its:ite,kts:kte,1:ens4)      ::  omeg

      real,   dimension (kts:kte) :: min_tend,distance
      integer,dimension (its:ite) :: kpbli,last_ierr

      integer :: i,j,k,kr,n,itf,jtf,ktf,ispc,zmax,status

      real :: dp,dq,exner, dtdt,PTEN,PQEN,PAPH,ZRHO,PAHFS,PQHFL,ZKHVFL,PGEOH
      real :: fixouts,dt_inv

      real,   dimension(mxp,myp,-1:5) :: dummy_precip
      integer :: imemory,irun,jlx,kk,kss,plume,ii_plume


      !----------------------------------------------------------------------
      !-do not change this
      itf=ite
      ktf=kte-1
      jtf=jte
      int_time = int_time + dt
      WHOAMI_ALL=mynum
      time_in   = time
      itime1_in = itime1
      !----------------------------------------------------------------------
      if(abs(C1) > 0.) USE_C1D = .true.

      !-- for the moisture adv trigger
      if(ADV_TRIGGER == 2) then
         call prepare_temp_pertubations(kts,kte,ktf,its,ite,itf,jts,jte,jtf,dt,xland,topt,zm &
            ,temp,rqvften,rthblten,rthften,Tpert_h,Tpert_v &
            ,AA1_CIN,AA1_BL)
         !-- for output only
         AA0(:,:) = Tpert_h(2, :,:)
         AA1(:,:) = Tpert_h(5,:,:)
         AA2(:,:) = Tpert_h(10,:,:)
         AA3(:,:) = Tpert_v(2,:,:)
      endif

      !-- big loop over j dimension
      do j = jts,jtf
         JCOL = J

         !-- initialization
         do I= its,itf
            rtgt(i,j)=1.0
         enddo
         do i= its,itf
            ztexec   (i) = 0.0
            zqexec   (i) = 0.0
            last_ierr(i) = -999
            fixout_qv(i) = 1.0
            !
            conprr     (i,j) = 0.0
            lightn_dens(i,j) = 0.0
            var2d      (i,j) = 0.0
            !--- (i,k)
            revsu_gf_2d (i,:) = 0.0
            prfil_gf_2d (i,:) = 0.0
            var3d_agf_2d(i,:) = 0.0
            var3d_bgf_2d(i,:) = 0.0
            Tpert_2d    (i,:) = 0.0
            !
            temp_tendqv(i,:) = 0.0
            !- tendencies (w/ maxiens)
            outt   (i,:,:)=0.0
            outu   (i,:,:)=0.0
            outv   (i,:,:)=0.0
            outq   (i,:,:)=0.0
            outqc  (i,:,:)=0.0
            outnice(i,:,:)=0.0
            outnliq(i,:,:)=0.0
            outbuoy(i,:,:)=0.0
         enddo

         if(APPLY_SUB_MP == 1) then
            do i= its,itf
               !- tendencies (w/ nmp and maxiens)
               outmpqi(:,i,:,:)=0.0
               outmpql(:,i,:,:)=0.0
               outmpcf(:,i,:,:)=0.0
            enddo
         endif
         do i= its,itf
            omeg (i,:,:)=0.0
         enddo
         !-- for the moisture adv trigger (Ma and Tan, AR 2009)
         if(ADV_TRIGGER == 2) then
            do k=kts,kte
               do i=its,itf
                  Tpert_2d(i,k)= Tpert_h(k,i,j) + Tpert_v(k,i,j) !- don't use "kr" here
               enddo
            enddo
         endif

         if(USE_TRACER_TRANSP==1) then
            out_chem = 0.0
         endif
         !
         if(autoconv == 2) then
            do I= its,itf
               ccn(i) = max( 100., ( 370.37*(0.01+MAX(0.,aot500(i,j))))**1.555 )
            enddo
         else
            do I= its,itf
               ccn(i) = 100.
            enddo
         endif

         do i=its,itf

            xlandi(i) = xland(i,j)!flag < 1 para land
                                  !flag  =1 para water
            psur  (i) = sfc_press(i,j)*1.e-2 ! mbar
            tsur  (i) = temp2m(i,j)
            ter11 (i) = max(0.,topt(i,j))
            kpbli (i) = kpbl(i,j)
            xlons (i) = lons(i,j)*180./3.14159
            xlats (i) = lats(i,j)*180./3.14159
         enddo

         do k=kts,ktf
            do i=its,itf
               kr=k   !+1   !<<<< only kr=k (the input was already converted to the BRAMS vertical grid,
                            !                see cup_grell3.f90 routine)

               !- heigths, current pressure, temp and water vapor mix ratio
               zo      (i,k)  = zt(kr,i,j)*rtgt(i,j)+topt(i,j)
               po      (i,k)  = press(kr,i,j)*1.e-2 !mbar
               temp_old(i,k)  = temp(kr,i,j)

               qv_old  (i,k)  = rvap     (kr,i,j) ! @ begin of the timestep
               qv_curr (i,k)  = curr_rvap(kr,i,j) ! current (after dynamics + physical processes called before GF)

               !- air density, TKE and cloud liq water mixing ratio
               rhoi    (i,k)  = 1.e2*po (i,k)/( 287.04*temp_old(i,k)*(1.+0.608*qv_old(i,k)))
               tkeg    (i,k)  = tkmin
               rcpg    (i,k)  = 0.

               !- wind velocities
               us      (i,k)  =  u (kr,i,j)
               vs      (i,k)  =  v (kr,i,j)
               omeg    (i,k,:)=  w (kr,i,j)
               dm2d    (i,k)  =  dm(kr,i,j)
               !omeg   (i,k,:)= -g*rho(kr,i,j)*w(kr,i,j)
               !-buoyancy excess
               buoy_exc2d(i,k)= buoy_exc(kr,i,j)
               !- temp/water vapor modified only by advection
               temp_new_ADV(i,k)= temp_old(i,k)  +  (rth_advten(kr,i,j) )*dt
               qv_new_ADV(i,k)=   qv_old(i,k)  +  (rqvften   (kr,i,j) )*dt
            enddo
         enddo

         if(APPLY_SUB_MP == 1) then
            do k=kts,ktf
               do i=its,itf
                  kr=k   !+1   !<<<< only kr=k
                  !- microphysics ice and liq mixing ratio, and cloud fraction of the host model
                  !- (only subsidence is applied)
                  mpqi   (:,i,k) = mp_ice  (:,kr,i,j) ! kg/kg
                  mpql   (:,i,k) = mp_liq  (:,kr,i,j) ! kg/kg
                  mpcf   (:,i,k) = mp_cf   (:,kr,i,j) ! 1
               enddo
            enddo
         endif
         if(USE_TRACER_TRANSP==1) then
            do k=kts,kte
               do i=its,itf
                  kr=k !+1
                  !- atmos composition
                  do ispc=1,mtp
                     se_chem(ispc,i,k) = max(mintracer, TRACER(i,j,flip(kr),ispc))
                  enddo
               enddo
            enddo
         endif
         !- pbl  (i) = depth of pbl layer (m)
         !- kpbli(i) = index of zo(i,k)
         !call get_zi_gf(j,its,ite,kts,kte,its,itf,ktf,ierrs,kpbli,pbl,&
         !             tkeg,rcpg,zo,ter11,tkmin)
         do i=its,itf
            pbl  (i)  = zo(i,kpbli(i)) - topt(i,j)
             !print*,"PBL=",kpbli(i),zo(i,kpbli(i)),topt(i,j),pbl  (i)
         enddo

         !- begin: for GATE soundings-------------------------------------------
         !- this section is intended for model developments only and must
         !- not be used for normal runs.
         if(USE_GATE) then
            if(CLEV_GRID == 0) stop "use_gate requires CLEV_GRID 1 or 2"
            if(USE_TRACER_TRANSP==1) then
               ispc_CO=1
               if( .not. allocated(Hcts)) allocate(Hcts(mtp))
               CHEM_NAME_MASK (:) = 1
               !--- dummy initization FSCAV
               do i=1,mtp
                  !FSCAV(i) = 0.1  !km^-1

                  FSCAV(i) = 1.e-5  !km^-1
                  Hcts(i)%hstar  = 0.0 !8.300e+4! 2.4E+3 !59.
                  Hcts(i)%dhr    = 0.0 !7400.   !5000.  !4200.
                  Hcts(i)%ak0    = 0.0
                  Hcts(i)%dak    = 0.0
                   ! H2O2      0.00000      8.300e+4    7400.00000       0.00000       0.00000
                   ! HNO3      0.00000      2.100e+5    8700.00000       0.00000       0.00000
                   ! NH3       0.00000      59.00000    4200.00000       0.00000       0.00000
                   ! SO2       0.00000      2.400e+3    5000.00000       0.00000       0.00000
               enddo
               do i=its,itf
                  se_chem(1:mtp,i,kts:kpbli(i)-1) = 1.+1.e-6
                  do k=kpbli(i),kte
                     se_chem(1:mtp,i,k) = 1.*exp(-(max(0.,0.9*float(k-kpbli(i)))/float(kpbli(i))))+1.e-6
                  enddo
                  do k=kts+1,kte-1
                     se_chem(1:mtp,i,k) = 1./3. *( se_chem(1:mtp,i,k) + se_chem(1:mtp,i,k-1) + se_chem(1:mtp,i,k+1))
                  enddo
               enddo
            endif

            !--- only for GATE soundingg
            if(trim(RUNDATA) == "GATE.dat") then
               jlx= jl
               !jlx= 1 ! to run with only one soundings
               !jlx= 42 ! to run with only one soundings

               do i=its,itf
                  do k=kts,kte
                     po       (i,k) = 0.5*(ppres(jlx,k)+ppres(jlx,min(kte,k+1)))
                     temp_old (i,k) = ptemp(jlx,k)+273.15
                     qv_old   (i,k) = pq(jlx,k)/1000.
                     us       (i,k) = pu(jlx,k)
                     vs       (i,k) = pv(jlx,k)
                     omeg     (i,k,:)=pvervel(jlx,k)
                     phil     (i,k) = pgeo(jlx,k)*g   !geo
                     rhoi     (i,k) = 1.e2*po(i,k)/(rgas*temp_old(i,k))
                  enddo

                  do k=kts,kte
                     mpql     (:,i,k) = 0.
                     mpql     (:,i,k) = 0.
                     mpcf     (:,i,k) = 0.
                     if(po(i,k) > 900. .or. po(i,k)<300.) cycle
                     pqen  =  exp((-3.e-5*(po(i,k)-550.)**2))
                     pten  =  min(1., (max(0.,(temp_old(i,k)-T_ice))/(T_0-T_ice))**2)
                     mpql  (:,i,k) =3.*pqen* pten
                     mpqi  (:,i,k) =3.*pqen*(1.- pten)
                     mpcf  (:,i,k) = (mpqi  (:,i,k)+mpql  (:,i,k))*100.
                  enddo

                  do k=kts,kte
                     zo       (i,k) = 0.5*(phil(i,k)+phil(i,min(kte,k+1)))/g    !meters
                  enddo
                  ter11(i)  = phil(i,1)/g  ! phil is given in g*h.
                  psur (i)  = ppres(jlx,1)
                  tsur (i)  = temp2m(i,j) !temp_old(i,1)
                  kpbli(i)  = 5
                  pbl  (i)  = zo(i,kpbli(i))
                  zws  (i)  = 1.0 ! wstar
                  do k=kts,ktf
                     temp_new(i,k) = temp_old(i,k) + dt *(zadvt(jlx,k)+zqr(jlx,k))/86400.
                     qv_new  (i,k) = qv_old  (i,k) + dt * zadvq(jlx,k)

                     temp_new_dp (i,k) = temp_old(i,k) + dt *(zadvt(jlx,k)+zqr(jlx,k))/86400.
                     qv_new_dp   (i,k) = qv_old  (i,k) + dt * zadvq(jlx,k)

                     temp_new_md (i,k) = temp_new_dp(i,k)
                     qv_new_md   (i,k) = qv_new_dp  (i,k)
                     temp_new_bl (i,k) = temp_new_dp(i,k)
                     qv_new_bl   (i,k) = qv_new_dp  (i,k)
                     temp_new_adv(i,k) = temp_old(i,k) + dt * zadvt(jlx,k)/86400.
                     qv_new_adv  (i,k) = qv_old  (i,k) + dt * zadvq(jlx,k)
                  enddo
               enddo
            endif
         endif !- end:   for GATE soundings-------------------------------------------
         !
         !- get execess T and Q for source air parcels
         do i=its,itf
            pten = temp_old(i,1)
            pqen = qv_old  (i,1)
            paph = 100.*psur(i)
            zrho = paph/(287.04*(temp_old(i,1)*(1.+0.608*qv_old(i,1))))
            !- sensible and latent sfc fluxes for the heat-engine closure
            h_sfc_flux (i)=zrho*cp *sflux_t(i,j)!W/m^2
            le_sfc_flux(i)=zrho*xlv*sflux_r(i,j)!W/m^2
            !
            !- local le and h fluxes for W*
            pahfs=-sflux_t(i,j) *zrho*1004.64  !W/m^2
            pqhfl=-sflux_r(i,j)                !kg/m^2/s
            !- buoyancy flux (h+le)
            zkhvfl= (pahfs/1004.64+0.608*pten*pqhfl)/zrho ! K m s-1
            !- depth of 1st model layer
            !- (zo(1)-top is ~ 1/2 of the depth of 1st model layer, => mult by 2)
            pgeoh =  2.*( zo(i,1)-topt(i,j) )*g ! m+2 s-2
            !-convective-scale velocity w*
            !- in the future, change 0.001 by ustar^3
            zws(i) = max(0.,0.001-1.5*0.41*zkhvfl*pgeoh/pten) ! m+3 s-3

            if(zws(i) > tiny(pgeoh)) then
               !-convective-scale velocity w*
               zws(i) = 1.2*zws(i)**.3333
               !- temperature excess
               ztexec(i)     = max(0.,-1.5*pahfs/(zrho*zws(i)*1004.64)) ! K
               !print*,"exce1=",pahfs,zrho,ztexec(i),zws(i),pgeoh,zo(i,1),topt(i,j)
               !call flush(6)
               !- moisture  excess
               zqexec(i)     = max(0.,-1.5*pqhfl/(zrho*zws(i)))        !kg kg-1
            endif   ! zws > 0
            !if(ztexec(i) > 1.) print*,"T",ztexec(i),h_sfc_flux(i)
            !if(zqexec(i)*1000 > 0.5) then 
            ! print*,"Q",1000*zqexec(i),le_sfc_flux(i),(ztexec(i)*cp+xlv*zqexec(i))/1000.
            !endif

            !
            !- zws for shallow convection closure (Grant 2001)
            !- depth of the pbl
            pgeoh = pbl(i)*g
            !-convective-scale velocity W* (m/s)
            zws(i) = max(0.,0.001-1.5*0.41*zkhvfl*pgeoh/pten)
            zws(i) = 1.2*zws(i)**.3333
         enddo
         !
         !------ CALL CUMULUS PARAMETERIZATION
         !

         do ii_plume = 1, maxiens
            if(ii_plume == 1) then
               plume = shal
               c0 = c0_shal
            endif
            if(ii_plume == 2) then
               plume = deep
               c0 = c0_deep
            endif
            if(ii_plume == 3) then
               plume = mid
               c0 = c0_mid
            endif

            hei_down_land  =  cum_hei_down_land  (plume)
            hei_down_ocean =  cum_hei_down_ocean (plume)
            hei_updf_land  =  cum_hei_updf_land  (plume)
            hei_updf_ocean =  cum_hei_updf_ocean (plume)
            max_edt_land   =  cum_max_edt_land   (plume)
            max_edt_ocean  =  cum_max_edt_ocean  (plume)
            fadj_massflx   =  cum_fadj_massflx   (plume)
            use_excess     =  cum_use_excess     (plume)
            ave_layer      =  cum_ave_layer      (plume)
            T_star         =  cum_t_star         (plume)
            !print*,"plume=",plume,shal,mid,deep

            if(icumulus_gf(plume) /= ON ) cycle

            !-- set minimum/max for excess of T and Q
            if(use_excess == 0) then
               cum_ztexec(:)= 0.
               cum_zqexec(:)= 0.
            elseif (use_excess == 1) then
               cum_ztexec(:)= ztexec(:)
               cum_zqexec(:)= zqexec(:)
            elseif (use_excess == 2) then
               do i=its,itf
                  cum_zqexec(i)=min(5.e-4, max(1.e-4,zqexec(i)))! kg kg^-1
                  cum_ztexec(i)=min(0.5,   max(0.2  ,ztexec(i)))! Kelvin
               enddo
            else
               do i=its,itf
                  if(xlandi(i) > 0.98) then ! ocean
                     cum_zqexec(i)=min(8.e-4, max(5.e-4,zqexec(i)))! kg kg^-1
                     cum_ztexec(i)=min(1.,    max(0.5  ,ztexec(i)))! Kelvin
                  else                      ! land
                     cum_ztexec(i)= ztexec(i)
                     cum_zqexec(i)= zqexec(i)
                  endif
               enddo
            endif
            !
            !--- replace 'q' and 't' excess in case of use of the cold pool scheme
            !
            if(convection_tracer == 1 .and. use_gustiness == 1 ) then           
            !---- only for deep plume
               if(plume == deep) then
                  k = 2 ! surface in brams
                  do i=its,itf
                   cum_ztexec(i)= (hexcp(k,i,j) - qexcp(k,i,j)*xlv)/cp
                   cum_zqexec(i)= qexcp(k,i,j)
                  enddo
               endif
            endif
            !
            !
            !-- shallow convection
            !
            if(plume == shal) then
               do i=its,itf
                  do k=kts,ktf
                     kr=k!+1 <<<<
                     if(use_gate) then
                        dhdt    (i,k)= cp*(temp_new_dp(i,k)-temp_old(i,k))+xlv*(qv_new_dp(i,k)-qv_old(i,k))
                        temp_new(i,k)= temp_new_dp(i,k)
                        qv_new  (i,k)= qv_new_dp  (i,k)
                     else

                        temp_new(i,k)=temp_old(i,k) + (rthblten(kr,i,j)+rthften(kr,i,j))*dt
                        qv_new  (i,k)=  qv_old(i,k) + (rqvblten(kr,i,j)+rqvften(kr,i,j))*dt
                        qv_new  (i,k)= max(smallerqv,qv_new  (i,k))

                        !- only pbl forcing changes moist static energy
                        dhdt(i,k)= cp   *(rthblten(kr,i,j)) +  xlv  *(rqvblten(kr,i,j))

                        !- all forcings change moist static energy
                        dhdt(i,k)=dhdt(i,k) + cp*rthften(kr,i,j) + xlv*rqvften(kr,i,j)

                     endif
                  enddo
               enddo
            endif
            !
            !--- deep convection
            if(plume == deep) then

               if(use_gate) then
                  do k=kts,ktf
                     do i=its,itf
                        temp_new(i,k) = temp_new_dp(i,k)
                        qv_new  (i,k) = qv_new_dp  (i,k)
                     enddo
                  enddo
               endif
               if(.not. use_gate) then
                  do k=kts,ktf
                     do i=its,itf
                        kr=k!+1 <<<<

                        temp_new    (i,k)= temp_old(i,k)  +  (rthblten(kr,i,j) + rthften (kr,i,j))*dt
                        qv_new    (i,k)=   qv_old(i,k)  +  (rqvblten(kr,i,j) + rqvften (kr,i,j))*dt

                        temp_new_BL (i,k)= temp_old(i,k)  +  (rthblten(kr,i,j) )*dt
                        qv_new_BL (i,k)=   qv_old(i,k)  +  (rqvblten(kr,i,j) )*dt

                        if(DICYCLE==100) then
                           temp_new(i,k) = temp_new(i,k) + outt(i,k,mid)*dt + outt(i,k,shal)*dt
                           qv_new  (i,k) =   qv_new(i,k) + outq(i,k,mid)*dt + outq(i,k,shal)*dt
                        endif

                     enddo
                  enddo
               endif
            endif
            !
            !--- mid/congestus type convection
            if(plume == mid) then

               if(use_gate) then
                  do k=kts,ktf
                     do i=its,itf
                        temp_new(i,k) = temp_new_dp(i,k)
                        qv_new  (i,k) = qv_new_dp  (i,k)
                     enddo
                  enddo
               endif
               !last_ierr(:)= ierr4d(:,j,deep)

               if(.not. use_gate) then
                  do i=its,itf
                     do k=kts,ktf
                        kr=k!+1 <<<<

                        temp_new(i,k)=temp_old(i,k) + (rthblten(kr,i,j)+rthften(kr,i,j))*dt
                        qv_new  (i,k)=  qv_old(i,k) + (rqvblten(kr,i,j)+rqvften(kr,i,j))*dt
                        qv_new  (i,k)= max(smallerqv,qv_new  (i,k))

                        !- only pbl forcing changes moist static energy
                        dhdt(i,k)= cp   *(rthblten(kr,i,j)) +  xlv  *(rqvblten(kr,i,j))

                        !- all forcings change moist static energy
                        dhdt(i,k)=dhdt(i,k) + cp*rthften(kr,i,j) + xlv*rqvften(kr,i,j)

                        !- temp/water vapor modified only by bl processes
                        temp_new_BL(i,k)= temp_old(i,k)  +  (rthblten(kr,i,j) )*dt
                        qv_new_BL  (i,k)= qv_old  (i,k)  +  (rqvblten(kr,i,j) )*dt

                     enddo
                  enddo
               endif
            endif
            !

            call CUP_gf(its,ite,kts,kte, itf,ktf, mtp, nmp, FSCAV  &
               ,cumulus_type  (plume)            &
               ,closure_choice(plume)            &
               ,cum_entr_rate (plume)            &
               ,cum_use_excess(plume)            &
               !- input data
               ,dx2d          (:,j)              &
               ,stochastic_sig(:,j)              &
               ,col_sat       (:,j)              &
               ,tke_pbl       (:,j)              &
               ,rh_dicycle_fct(:,j)              &
               ,wlpool        (:,j)              &
               ,dt                               &
               ,kpbli                            &
               ,cum_ztexec                       &
               ,cum_zqexec                       &
               ,ccn                              &
               ,rhoi                             &
               ,omeg                             &
               ,temp_old                         &
               ,qv_old                           &
               ,ter11                            &
               , h_sfc_flux                      &
               ,le_sfc_flux                      &
               ,xlons                            &
               ,xlats                            &
               ,xlandi                           &
               ,temp_new                         &
               ,qv_new                           &
               ,temp_new_BL                      &
               ,qv_new_BL                        &
               ,temp_new_ADV                     &
               ,qv_new_ADV                       &
               ,zo                               &
               ,po                               &
               ,tsur                             &
               ,psur                             &
               ,us                               &
               ,vs                               &
               ,dm2d                             &
               ,se_chem                          &
               ,zws                              &
               ,dhdt                             &
               ,buoy_exc2d                       &
               ,mpqi                             &
               ,mpql                             &
               ,mpcf                             &
               ,last_ierr            (:)         &
               !output data
               ,outt                 (:,:,plume) &
               ,outq                 (:,:,plume) &
               ,outqc                (:,:,plume) &
               ,outu                 (:,:,plume) &
               ,outv                 (:,:,plume) &
               ,outnliq              (:,:,plume) &
               ,outnice              (:,:,plume) &
               ,outbuoy              (:,:,plume) &
               ,outmpqi            (:,:,:,plume) &
               ,outmpql            (:,:,:,plume) &
               ,outmpcf            (:,:,:,plume) &
               ,out_chem           (:,:,:,plume) &
               !- for convective transport
               ,ierr4d               (:,j,plume) &
               ,jmin4d               (:,j,plume) &
               ,klcl4d               (:,j,plume) &
               ,k224d                (:,j,plume) &
               ,kbcon4d              (:,j,plume) &
               ,ktop4d               (:,j,plume) &
               ,kstabi4d             (:,j,plume) &
               ,kstabm4d             (:,j,plume) &
               ,cprr4d               (:,j,plume) &
               ,xmb4d                (:,j,plume) &
               ,edt4d                (:,j,plume) &
               ,pwav4d               (:,j,plume) &
               ,sigma4d              (:,j,plume) &
               ,pcup5d             (:,j,:,plume) &
               ,up_massentr5d      (:,j,:,plume) &
               ,up_massdetr5d      (:,j,:,plume) &
               ,dd_massentr5d      (:,j,:,plume) &
               ,dd_massdetr5d      (:,j,:,plume) &
               ,zup5d              (:,j,:,plume) &
               ,zdn5d              (:,j,:,plume) &
               ,prup5d             (:,j,:,plume) &
               ,prdn5d             (:,j,:,plume) &
               ,clwup5d            (:,j,:,plume) &
               ,tup5d              (:,j,:,plume) &
               ,conv_cld_fr5d      (:,j,:,plume) &
               !-- for debug/diag
               ,AA0(:,j),AA1(:,j),AA1_ADV(:,j),AA1_RADPBL(:,j),AA2(:,j),AA3(:,j) &
               ,AA1_BL(:,j),AA1_CIN(:,j),TAU_BL(:,j),TAU_EC(:,j) &
               !-- for diag
               ,lightn_dens  (:,j)               &
               ,var2d        (:,j)               &
               ,revsu_gf_2d                      &
               ,prfil_gf_2d                      &
               ,var3d_agf_2d                     &
               ,var3d_bgf_2d                     &
               ,Tpert_2d                         &
               )

         enddo !- plume

         !--- reset ierr4d to value different of zero in case the correspondent
         !--- plume (shalllow, congestus, deep) was not actually used
         do n=1,maxiens
            if(icumulus_gf(n) == OFF ) ierr4d (:,j,n) = -99
         enddo

         do i=its,itf
            do_this_column(i,j) = 0
            loop1:  do n=1,maxiens
               if(ierr4d (i,j,n) == 0 ) then
                  do_this_column(i,j) = 1
                  exit loop1
               endif
            enddo loop1
         enddo
         !----------- check for negative water vapor mix ratio
         do i=its,itf
            if(do_this_column(i,j) == 0) cycle
            do k = kts,ktf
               temp_tendqv(i,k)= outq (i,k,shal) + outq (i,k,deep) + outq (i,k,mid )
            enddo

            do k = kts,ktf
               distance(k)= qv_curr(i,k) + temp_tendqv(i,k) * dt
            enddo

            if(minval(distance(kts:ktf)) < 0.) then
               zmax   =  MINLOC(distance(kts:ktf),1)

               if( abs(temp_tendqv(i,zmax) * dt) <  mintracer) then
                  fixout_qv(i)= 0.999999
                !fixout_qv(i)= 0.
               else
                  fixout_qv(i)= ( (smallerQV - qv_curr(i,zmax))) / (temp_tendqv(i,zmax) *dt)
               endif
               fixout_qv(i)=max(0.,min(fixout_qv(i),1.))
            endif
         enddo
         !------------ feedback
         !-- deep convection
         do i=its,itf
            if(do_this_column(i,j) == 0) cycle
            cprr4d(i,j,deep) =  cprr4d(i,j,deep)* fixout_qv(i)
            cprr4d(i,j,mid)  =  cprr4d(i,j,mid) * fixout_qv(i)
            cprr4d(i,j,shal) =  cprr4d(i,j,shal)* fixout_qv(i)
            CONPRR(i,j)= (cprr4d(i,j,deep) + cprr4d(i,j,mid) + cprr4d(i,j,shal))
            CONPRR(i,j)= max(0.,CONPRR(i,j))
         enddo

         !-- deep + shallow + mid convection
         do i = its,itf
            if(do_this_column(i,j) == 0) cycle
            do k = kts,kte
               kr=k!+1
               !- feedback the tendencies from convection
               RTHCUTEN (kr,i,j)= (outt (i,k,shal) + outt (i,k,deep) + outt (i,k,mid )) *fixout_qv(i)

               RQVCUTEN (kr,i,j)= (outq (i,k,shal) + outq (i,k,deep) + outq (i,k,mid )) *fixout_qv(i)

               RQCCUTEN (kr,i,j)= (outqc(i,k,shal) + outqc(i,k,deep) + outqc(i,k,mid )) *fixout_qv(i)

               REVSU_GF (kr,i,j)= revsu_gf_2d(i,k)*fixout_qv(i) !-- already contains deep and mid amounts.

               !---these arrays are only for the deep plume mode
               PRFIL_GF (kr,i,j)= prfil_gf_2d (i,k)*fixout_qv(i) !-- ice/liq prec flux of the deep plume
               !VAR3d_aGF(kr,i,j)= var3d_gf_2d(i,k)               !-- vertical velocity of the deep plume
               VAR3d_aGF(kr,i,j)= outt (i,k,mid )*fixout_qv(i)   !--
               VAR3d_bGF(kr,i,j)= outq (i,k,mid )*fixout_qv(i)   !--

               if(icumulus_gf(shal) == OFF) then
                  VAR3d_cGF(kr,i,j)= outqc (i,k,deep)*fixout_qv(i)  !--
                  VAR3d_dGF(kr,i,j)= outqc (i,k,mid )*fixout_qv(i)  !--
               else
                  VAR3d_cGF(kr,i,j)= outt (i,k,shal)*fixout_qv(i)   !--
                  VAR3d_dGF(kr,i,j)= outq (i,k,shal)*fixout_qv(i)   !--
               endif

            enddo
         enddo
         if(USE_MOMENTUM_TRANSP > 0) then
            do i = its,itf
               if(do_this_column(i,j) == 0) cycle
               do k = kts,kte
                  kr=k!+1
                  RUCUTEN (kr,i,j)= (outU(i,k,deep)+outU(i,k,mid)+outU(i,k,shal)) *fixout_qv(i)
                  RVCUTEN (kr,i,j)= (outV(i,k,deep)+outV(i,k,mid)+outV(i,k,shal)) *fixout_qv(i)
               enddo
            enddo
         endif

         if(APPLY_SUB_MP == 1) then
            do i = its,itf
               if(do_this_column(i,j) == 0) cycle
               do k = kts,kte
                  kr=k!+1
                  SUB_MPQL (:,kr,i,j)= (outmpql(:,i,k,deep)+outmpql(:,i,k,mid)+outmpql(:,i,k,shal)) *fixout_qv(i)
                  SUB_MPQI (:,kr,i,j)= (outmpqi(:,i,k,deep)+outmpqi(:,i,k,mid)+outmpqi(:,i,k,shal)) *fixout_qv(i)
                  SUB_MPCF (:,kr,i,j)= (outmpcf(:,i,k,deep)+outmpcf(:,i,k,mid)+outmpcf(:,i,k,shal)) *fixout_qv(i)
               enddo
            enddo
         endif

         if(LIQ_ICE_NUMBER_CONC == 1) then
            do i = its,itf
               if(do_this_column(i,j) == 0) cycle
               do k = kts,kte
                  kr=k!+1
                  RNICUTEN (kr,i,j)= (outnice(i,k,shal) + outnice(i,k,deep) + outnice(i,k,mid )) *fixout_qv(i)
                  RNLCUTEN (kr,i,j)= (outnliq(i,k,shal) + outnliq(i,k,deep) + outnliq(i,k,mid )) *fixout_qv(i)
               enddo
            enddo
         endif

         if(USE_TRACER_TRANSP==1) then
            do i = its,itf
               if(do_this_column(i,j) == 0) cycle
               do k = kts,kte
                  kr=k!+1
                  RCHEMCUTEN (:,kr,i,j)= (out_CHEM(:,i,k,deep) +out_CHEM(:,i,k,mid)+out_CHEM(:,i,k,shal)) *fixout_qv(i)
               enddo
            enddo

            !- constrain positivity for tracers
            do i = its,itf
               if(do_this_column(i,j) == 0) cycle

               do ispc=1,mtp
                  if(CHEM_NAME_MASK (ispc) == 0 ) cycle

                  do k=kts,ktf
                     distance(k)= se_chem(ispc,i,k) + RCHEMCUTEN(ispc,k,i,j)* dt
                  enddo

                  !-- fixer for mass of tracer
                  if(minval(distance(kts:ktf)) < 0.) then
                     zmax   =  MINLOC(distance(kts:ktf),1)

                     if( abs(RCHEMCUTEN(ispc,zmax,i,j)*dt) <  mintracer) then
                        fixouts= 0.999999
                      !fixouts= 0.
                     else
                        fixouts=  ( (mintracer - se_chem(ispc,i,zmax))) / (RCHEMCUTEN(ispc,zmax,i,j)*dt)
                     endif
                     if(fixouts > 1. .or. fixouts <0.)fixouts=0.

                     RCHEMCUTEN(ispc,kts:ktf,i,j)=fixouts*RCHEMCUTEN(ispc,kts:ktf,i,j)
                  endif
               enddo
            enddo
         endif

         if(CONVECTION_TRACER==1) then
            do i = its,itf
               if(do_this_column(i,j) == 0) cycle
               do k = kts,kte
                  kr=k!+1
                  RBUOYCUTEN (kr,i,j)= (outbuoy(i,k,deep)+outbuoy(i,k,mid)+outbuoy(i,k,shal)) *fixout_qv(i)
                  !print*,"RBUOYCUTEN", RBUOYCUTEN (kr,i,j),outbuoy(i,k,deep),&
                  !             outbuoy(i,k,mid),outbuoy(i,k,shal)
               enddo
            enddo
         
           !----- for output only
           !if(use_gustiness==1 .or. use_gustiness ==2 ) then 
               !print*,"H-T",cp*1.1*maxval(sflux_t(:,j))&
               !            ,xlv*1.1*maxval(sflux_r(:,j)),maxval(ztexec),maxval(zqexec)
               !sflux_t(:,j) = ztexec(:)
               !sflux_r(:,j) = zqexec(:)
           !endif 
         endif

         !-----memory
         !AA3(:,j)=cprr4d(:,j,deep) *fixout_qv(:)
         !AA2(:,j)=cprr4d(:,j,mid)  *fixout_qv(:)
      enddo

   end subroutine GF_GEOS5_DRV
   !---------------------------------------------------------------------------------------------------

   subroutine CUP_gf(its,ite,kts,kte ,itf,ktf, mtp, nmp , FSCAV   &
      ,cumulus           &
      ,ichoice           &
      ,entr_rate_input   &
      ,use_excess        &
      !input data
      ,dx                &
      ,stochastic_sig    &
      ,col_sat           &
      ,tke_pbl           &
      ,rh_dicycle_fct    &
      ,wlpool            &
      ,dtime             &
      ,kpbl              &
      ,ztexec            &
      ,zqexec            &
      ,ccn               &
      ,rho               &
      ,omeg              &
      ,t                 &
      ,q                 &
      ,z1                &
      , h_sfc_flux       &
      ,le_sfc_flux       &
      ,xlons             &
      ,xlats             &
      ,xland             &
      ,tn                &
      ,qo                &
      ,tn_bl             &
      ,qo_bl             &
      ,tn_adv            &
      ,qo_adv            &
      ,zo                &
      ,po                &
      ,tsur              &
      ,psur              &
      ,us                &
      ,vs                &
      ,dm2d              &
      ,se_chem           &
      ,zws               &
      ,dhdt              &
      ,buoy_exc          &
      ,mpqi              &
      ,mpql              &
      ,mpcf              &
      ,last_ierr         &
      !output data
      ,outt              &
      ,outq              &
      ,outqc             &
      ,outu              &
      ,outv              &
      ,outnliq           &
      ,outnice           &
      ,outbuoy           &
      ,outmpqi           &
      ,outmpql           &
      ,outmpcf           &
      ,out_chem          &
      !- for convective transport
      ,ierr              &
      ,jmin              &
      ,klcl              &
      ,k22               &
      ,kbcon             &
      ,ktop              &
      ,kstabi            &
      ,kstabm            &
      ,pre               &
      ,xmb               &
      ,edto              &
      ,pwavo             &
      ,sig               &
      ,po_cup            &
      ,up_massentro      &
      ,up_massdetro      &
      ,dd_massentro      &
      ,dd_massdetro      &
      ,zuo               &
      ,zdo               &
      ,pwo               &
      ,pwdo              &
      ,qrco              &
      ,tup               &
      ,clfrac            &
      !- for convective transport-end
      !- for debug/diag
      ,AA0_,AA1_,AA1_ADV_,AA1_RADPBL_,AA2_,AA3_,AA1_BL_,AA1_CIN_,TAU_BL_,TAU_EC_   &
      ,lightn_dens       &
      ,var2d             &
      ,revsu_gf          &
      ,prfil_gf          &
      ,var3d_agf         &
      ,var3d_bgf         &
      ,Tpert             &
      )
      implicit none

      !-local settings
      logical, parameter:: USE_INV_LAYERS=.true.

      character*(*),intent(in) :: cumulus
      integer      ,intent(in) :: itf,ktf,its,ite,kts,kte,ichoice,use_excess,mtp, nmp
      integer      ,intent(inout),  dimension (its:ite) ::   kpbl,last_ierr
      !
      ! outtem = output temp tendency (per s)
      ! outq   = output q tendency (per s)
      ! outqc  = output qc tendency (per s)
      ! pre    = output precip
      real,    dimension (its:ite,kts:kte) ,intent (inout)   ::       &
          outu,outv,outt,outq,outqc,outbuoy,revsu_gf,prfil_gf,var3d_agf ,var3d_bgf &
         ,outnliq ,outnice

      real,    dimension (its:ite)         ,intent (out  )   ::       &
         pre,sig,lightn_dens,var2d

      !
      ! basic environmental input includes
      ! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
      ! convection for this call only and at that particular gridpoint
      !
      real,    dimension (its:ite,kts:kte)       ,intent (inout   )    ::   &
         dhdt,rho,t,po,us,vs,tn,dm2d,buoy_exc,tn_bl,tn_adv
      real,    dimension (its:ite,kts:kte,1:ens4),intent (inout)    ::   &
         omeg
      real,    dimension (its:ite,kts:kte)       ,intent (inout)    ::   &
         q,qo,Tpert,qo_bl,qo_adv
      real,    dimension (its:ite)               ,intent (inout   )    ::   &
         ccn,Z1,PSUR,xland,xlons,xlats, h_sfc_flux,le_sfc_flux,tsur,dx

      real,    dimension (its:ite)               ,intent (in   )    ::   &
         col_sat,stochastic_sig,tke_pbl
      real,    dimension (its:ite)               ,intent (inout)    ::   &
         zws,ztexec,zqexec, rh_dicycle_fct,wlpool
      real                                       ,intent (in   )    ::   &
         dtime,entr_rate_input
      real,    dimension (nmp,its:ite,kts:kte)   ,intent (inout)   ::   &
         mpqi,mpql,mpcf
      real,    dimension (nmp,its:ite,kts:kte)   ,intent (inout)   ::   &
         outmpqi,outmpql,outmpcf

      !
      ! local ensemble dependent variables in this routine
      real,    dimension (its:ite,1:maxens2) ::                         &
         edtc
      real,    dimension (its:ite,1:ensdim)        ::                   &
         xf_ens,pr_ens
      !
      !*******the following are your basic environmental
      !          variables. They carry a "_cup" if they are
      !          on model cloud levels (staggered). They carry
      !          an "o"-ending (z becomes zo), if they are the forced
      !          variables. They are preceded by x (z becomes xz)
      !          to indicate modification by some typ of cloud
      !
      ! z           = heights of model levels
      ! q           = environmental mixing ratio
      ! qes         = environmental saturation mixing ratio
      ! t           = environmental temp
      ! p           = environmental pressure
      ! he          = environmental moist static energy
      ! hes         = environmental saturation moist static energy
      ! z_cup       = heights of model cloud levels
      ! q_cup       = environmental q on model cloud levels
      ! qes_cup     = saturation q on model cloud levels
      ! t_cup       = temperature (Kelvin) on model cloud levels
      ! p_cup       = environmental pressure
      ! he_cup = moist static energy on model cloud levels
      ! hes_cup = saturation moist static energy on model cloud levels
      ! gamma_cup = gamma on model cloud levels
      !
      !
      ! hcd = moist static energy in downdraft
      ! zd normalized downdraft mass flux
      ! dby = buoancy term
      ! entr = entrainment rate
      ! zd   = downdraft normalized mass flux
      ! entr= entrainment rate
      ! hcd = h in model cloud
      ! bu = buoancy term
      ! zd = normalized downdraft mass flux
      ! gamma_cup = gamma on model cloud levels
      ! qcd = cloud q (including liquid water) after entrainment
      ! qrch = saturation q in cloud
      ! pwd = evaporate at that level
      ! pwev = total normalized integrated evaoprate (I2)
      ! entr= entrainment rate
      ! z1 = terrain elevation
      ! entr = downdraft entrainment rate
      ! jmin = downdraft originating level
      ! kdet = level above ground where downdraft start detraining
      ! psur    = surface pressure
      ! z1      = terrain elevation
      ! pr_ens  = precipitation ensemble
      ! xf_ens  = mass flux ensembles
      ! massfln = downdraft mass flux ensembles used in next timestep
      ! omeg    = omega from large scale model
      ! mconv   = moisture convergence from large scale model
      ! zd      = downdraft normalized mass flux
      ! zu      = updraft normalized mass flux
      ! dir     = "storm motion"
      ! mbdt    = arbitrary numerical parameter
      ! dtime   = dt over which forcing is applied
      ! iact_gr_old = flag to tell where convection was active
      ! kbcon       = LFC of parcel from k22
      ! k22         = updraft originating level
      ! ichoice     = flag if only want one closure (usually set to zero!)
      ! dby  = buoancy term
      ! ktop = cloud top (output)
      ! xmb  = total base mass flux
      ! hc   = cloud moist static energy
      ! hkb  = moist static energy at originating level
      logical :: keep_going
      character*128 :: ierrc(its:ite)

      real,    dimension (its:ite,kts:kte) ::                           &
         entr_rate_2d,mentrd_rate_2d                                    &
         , he, hes, qes, z,heo,heso,qeso,zo, zu,zd                       &
         ,xhe,xhes,xqes,xz,xt,xq                                         &
         , qes_cup, q_cup, he_cup, hes_cup, z_cup, p_cup, gamma_cup, t_cup &
         ,qeso_cup,qo_cup,heo_cup,heso_cup,zo_cup,       gammao_cup,tn_cup &
         ,xqes_cup,xq_cup,xhe_cup,xhes_cup,xz_cup                        &
         ,xt_cup,hcot,evap_bcb                                           &
         ,dby,hc,clw_all                                                 &
         ,dbyo,qco,qrcdo,hcdo,qcdo,dbydo,hco                             &
         ,xdby,xzu,xzd,xhc,cupclw,pwo_eff,                               &

         ! cd  = detrainment function for updraft
         ! cdd = detrainment function for downdraft
         ! dellat = change of temperature per unit mass flux of cloud ensemble
         ! dellaq = change of q per unit mass flux of cloud ensemble
         ! dellaqc = change of qc per unit mass flux of cloud ensemble

         cd,cdd,dellah,dellaq,dellat,dellaqc,dsubq,dsubh,dellabuoy,      &
         u_cup,v_cup,uc,vc,ucd,vcd,dellu,dellv,                          &
         up_massentr,up_massdetr,dd_massentr,dd_massdetr,                &
         subten_H,subten_Q,subten_T

      ! aa0 cloud work function for downdraft
      ! edt = epsilon
      ! aa0     = cloud work function without forcing effects
      ! aa1     = cloud work function with forcing effects
      ! xaa0    = cloud work function with cloud effects
      ! edt     = epsilon

      real,    dimension (its:ite) ::                             &
         edt,edtx,aa1,aa0,xaa0,hkb,                               &
         hkbo,xhkb,qkb, pwevo,bu,bud,cap_max,xland1,vshear,       &
         cap_max_increment,psum,psumh,sigd,mconv,rescale_entrain,entr_rate,mentrd_rate

      integer,    dimension (its:ite) ::                         &
         kzdown,kdet,kb, ierr2,ierr3,kbmax
      integer :: iloop,nall,iedt,nens,nens3,ki,I,K,KK,iresult,nvar,nvarbegin
      integer :: jprnt,k1,k2,kbegzu,kdefi,kfinalzu,kstart,jmini,imid,k_free_trop

      real :: day,dz,dzo,radius,entrd_rate,                         &
         zcutdown,depth_min,zkbmax,z_detr,zktop,               &
         massfld,dh,trash,frh,xlamdd,radiusd,frhd,effec_entrain
      real :: detdo1,detdo2,entdo,dp,subin,detdo,entup,             &
         detup,subdown,entdoj,entupk,detupk,totmas,min_entr_rate
      real :: tot_time_hr,beta,env_mf,env_mf_p,env_mf_m
      real :: dts,denom,denomU

      real,    dimension (its:ite,1:maxens3) ::  xff_mid
      real,    dimension (kts:kte)   :: dummy1,dummy2
      integer :: iversion,bl=1,fa=2,step
      real :: umean

      real, dimension (its:ite)         :: aa0_bl,aa1_bl,tau_bl,tau_ecmwf,wmean,aa1_fa,aa1_tmp,hkbo_x &
         ,aa2,aa3,cin0,cin1,edtmax,edtmin,aa1_lift,aa_tmp,aa_ini,aa_adv&
         ,daa_adv_dt,wlpool_bcon
      real, dimension (its:ite,kts:kte) :: tn_x, qo_x, qeso_x, heo_x, heso_x,zo_cup_x &
         ,qeso_cup_x,qo_cup_x, heo_cup_x,heso_cup_x,po_cup_x&
         ,gammao_cup_x,tn_cup_x,hco_x,DBYo_x,u_cup_x,v_cup_x

      real, dimension (its:ite,kts:kte) :: xhe_x,xhes_x,xt_x,xq_x,xqes_x, &
         xqes_cup_x,xq_cup_x,xhe_cup_x,xhes_cup_x,gamma_cup_x,xt_cup_x
      real, dimension (its:ite)         :: xaa0_x,xk_x

      real, dimension(its:ite) :: xf_dicycle,mbdt,xf_coldpool
      real :: C_up, E_dn,G_rain,trash2, pgc,bl2dp,trash3, ke_mx
      character(len=2) :: cty

      real   :: dsubh_aver,dellah_aver,x_add,cap_max_inc,tlll,plll,rlll,tlcl,plcl,dzlcl,zlll
      integer:: start_k22,start_level(its:ite)
      real,    dimension (its:ite,kts:kte) ::  dtempdz
      integer, dimension (its:ite,kts:kte) ::  k_inv_layers
      integer :: ipr=0,jpr=0,fase

      real,    dimension (its:ite,kts:kte) ::  vvel2d,tempco,tempcdo
      real,    dimension (its:ite        ) ::  vvel1d, x_add_buoy

      real,    dimension (its:ite,kts:kte) ::  p_liq_ice,melting_layer,melting
      real,    dimension (its:ite,kts:kte) ::  c1d

      real,    dimension (its:ite)         ::lambau_dn,lambau_dp
      real,    dimension (its:ite,kts:kte) ::  up_massentru,up_massdetru,&
         dd_massentru,dd_massdetru

      real,    dimension (its:ite) :: q_wetbulb,t_wetbulb
      integer :: i_wb=0

      real,    dimension (its:ite) :: col_sat_adv,Q_adv,alpha_adv,aa1_radpbl,aa1_adv
      integer, dimension (its:ite) :: ierr_dummy

      real,    dimension (its:ite) :: p_cwv_ave,cape, depth_neg_buoy,frh_bcon &
         ,check_sig,random,rh_entr_factor

      real,    dimension (its:ite,kts:kte) :: prec_flx,evap_flx,qrr

      real,    dimension (its:ite,shall_closures) :: xff_shal

      !- atmos composition arrays
      real, dimension (mtp ),               intent (in)    ::   fscav
      !real, dimension (mtp,its:ite,kts:kte),intent (in)    ::   se_chem
      real, dimension (mtp,its:ite,kts:kte),intent (inout) ::   se_chem
      real, dimension (mtp,its:ite,kts:kte),intent (inout) ::   out_chem

      !-locals
      integer :: ispc,kmp,istep,lstep
      real, dimension (mtp,its:ite,kts:kte) ::  se_cup_chem,sc_up_chem,sc_dn_chem,pw_up_chem,pw_dn_chem
      real, dimension (mtp,its:ite)         ::  tot_pw_up_chem,tot_pw_dn_chem
      real, dimension (mtp,kts:kte)         ::  trcflx_in,sub_tend,ddtr,ddtr_upd,zenv_diff,fp_mtp,fm_mtp
      real, dimension (mtp)                 ::  min_tend_chem,dummy_chem,delS_up,delS_dn,env_sub,outchem1
      real, dimension (nmp,its:ite,kts:kte) ::  dellampqi,dellampql,dellampcf
      real, dimension (kts:kte)             ::  aa,bb,cc,ddu,ddv,ddh,ddq,fp,fm
      real, dimension (nmp,kts:kte)         ::  dd
      real, dimension (its:ite,kts:kte)     ::  massflx,zenv,rho_hydr,alpha_H,alpha_Q
      real                                  ::  evap_(mtp),wetdep_(mtp),trash_(mtp),trash2_(mtp) &
         ,massi,massf,dtime_max,evap,wetdep,residu_(mtp)

      !----------------------------------------------------------------------
      integer, dimension (its:ite)      ,intent (inout)  :: &
          ierr              &
         ,jmin              &
         ,klcl              &
         ,k22               &
         ,kbcon             &
         ,ktop              &
         ,kstabi            &
         ,kstabm

      real,  dimension (its:ite)        ,intent (inout)  :: &
          xmb               &
         ,edto              &
         ,pwavo
      real,  dimension (its:ite,kts:kte),intent (inout)  :: &
          po_cup            &
         ,up_massentro      &
         ,up_massdetro      &
         ,dd_massentro      &
         ,dd_massdetro      &
         ,zuo               &
         ,zdo               &
         ,pwo               &
         ,pwdo              &
         ,qrco              &
         ,tup               &
         ,clfrac
      !----------------------------------------------------------------------
      !-- debug/diag
      real,  dimension (its:ite)        ,intent (inout)  :: &
         aa0_,aa1_,aa2_,aa3_,aa1_bl_,aa1_cin_,tau_bl_,tau_ec_, AA1_RADPBL_,AA1_ADV_
      real,  dimension (its:ite,kts:kte) :: dtdt,dqdt
      real :: s1,s2,q1,q2,rzenv,factor,CWV,entr_threshold,resten_H,resten_Q,resten_T
      integer :: status
      real :: alp0,beta1,beta2,dp_p,dp_m,delt1,delt2,delt_Tvv,wkf,ckf,wkflcl(its:ite),rcount
      real :: min_deep_top,min_shall_top
      character(200) :: lixo
      integer :: X_kte,X_k,X_i,X_jcol
      real,    dimension (its:ite)           :: X_dx,X_stochastic_sig
      real,    dimension (kts:kte,8)       ::  tend2d
      real,    dimension (8)               ::  tend1d
      real,    dimension (its:ite,8)       ::  check_cons_I,check_cons_F

      !-- only for debug (atmos composition)
      real, allocatable, dimension (:,:,:),save    ::   se_chem_update
          !----------------------------------------------------------------------


      !--only for debug
      if(use_gate) then
         if( .not. allocated(se_chem_update)) allocate(se_chem_update(3,its:ite,kts:kte))
         if(jl==1) then
         !    se_chem_update(1,:,:) = mpql (lsmp,:,:)
         !    se_chem_update(2,:,:) = mpqi (lsmp,:,:)
         !    se_chem_update(3,:,:) = mpcf (lsmp,:,:)
         else
         !    mpql (lsmp,:,:)= se_chem_update(1,:,:)
         !    mpqi (lsmp,:,:)= se_chem_update(2,:,:)
         !    mpcf (lsmp,:,:)= se_chem_update(3,:,:)
         endif
      endif
      !--only for debug

      !
      !--- maximum depth (mb) of capping inversion (larger cap = no convection)
      !
      if(ZERO_DIFF==1 .or. MOIST_TRIGGER==0) then
         if(cumulus == 'deep'   ) then
            cap_max_inc=20.
         endif ! cap_maxs=50.
         if(cumulus == 'mid'    ) then
            cap_max_inc=10.
         endif ! cap_maxs=50.
         if(cumulus == 'shallow') then
            cap_max_inc=25.
         endif ! cap_maxs=50.
      else
         if(cumulus == 'deep'   ) then
            cap_max_inc=90.
         endif !--- test cap_maxs=10.  ; cap_max_inc=50
         if(cumulus == 'mid'    ) then
            cap_max_inc=90.
         endif !--- test cap_maxs=10.  ; cap_max_inc=50
         if(cumulus == 'shallow') then
            cap_max_inc=10.
         endif !--- test cap_maxs=25.  ; cap_max_inc=50
      endif
      !
      !--- lambda_U parameter for momentum transport
      !
      if(cumulus == 'deep'   ) then
         lambau_dp (:) = lambau_deep
         lambau_dn (:) = lambau_shdn
      endif
      if(cumulus == 'mid'    ) then
         lambau_dp (:) = lambau_shdn
         lambau_dn (:) = lambau_shdn
      endif
      if(cumulus == 'shallow') then
         lambau_dp (:) = lambau_shdn
         lambau_dn (:) = lambau_shdn
      endif

      if(pgcon .ne. 0.) then
         lambau_dp (:) = 0.
         lambau_dn (:) = 0.
      endif

      do i=its,itf
         kbmax  (i) = 1
         kstabm (i) = ktf-1
         ierr2  (i) = 0
         ierr3  (i) = 0
         xland1 (i) = xland(i) ! 1.
         cap_max(i) = cap_maxs
         ierrc  (i) = "ierrtxt"
         aa0    (i) = 0.0
         aa1    (i) = 0.0
         aa2    (i) = 0.0
         aa3    (i) = 0.0
         aa1_bl (i) = 0.0
         aa1_fa (i) = 0.0
         aa0_bl (i) = 0.0
         q_adv  (i) = 0.0
         aa1_radpbl(i) = 0.0
         aa1_adv   (i) = 0.0
         alpha_adv (i) = 0.0
         cin1   (i) = 0.0
         xk_x   (i) = 0.0
         edt    (i) = 0.0
         edto   (i) = 0.0
         tau_bl (i) = 0.0
         q_wetbulb (i) = 0.0
         t_wetbulb (i) = 0.0
         tau_ecmwf (i) = 0.0
         xf_dicycle(i) = 0.0
         x_add_buoy(i) = 0.0
         xf_coldpool(i) = 0.0
         wlpool_bcon(i) = 0.0
         z     (i,:) = zo(i,:)
         xz    (i,:) = zo(i,:)
         hcdo  (i,:) = 0.0
         cupclw(i,:) = 0.0
         qrcdo (i,:) = 0.0
         hcot  (i,:) = 0.0
         c1d   (i,:) = 0.0
         xf_ens(i,:) = 0.0
         pr_ens(i,:) = 0.0
         evap_bcb(i,:) = 0.0
         cap_max_increment(i)=cap_max_inc
      enddo
      !
      !---  create a real random number in the interval [-use_random_num, +use_random_num]
      !
      if( cumulus == 'deep' .and. use_random_num > 1.e-6) then
         call gen_random(its,ite,use_random_num,random)
      else
         random = 0.0
      endif
      !
      !
      !--- max/min allowed value for epsilon (ratio downdraft base mass flux/updraft
      !    base mass flux
      !--  note : to make the evaporation stronger => increase "edtmin"
      if( cumulus == 'shallow') then
         edtmin(:)=0.0
         edtmax(:)=0.0
      endif
      if(cumulus == 'mid'     ) then
         do i=its,itf
            if(xland(i) > 0.99 ) then !- over water
               edtmin(i)=0.1
               edtmax(i)=MAX_EDT_OCEAN  !
            else!- over land
               edtmin(i)=0.1
               edtmax(i)=MAX_EDT_LAND  !
            endif
         enddo
         if(c0_mid < 1.e-8) edtmin(:)=0.0
      endif
      if(cumulus == 'deep'    ) then
         do i=its,itf
            if(xland(i) > 0.99 ) then !- over water
               edtmin(i)=0.1
               edtmax(i)=MAX_EDT_OCEAN  !
            else!- over land
               edtmin(i)=0.1
               edtmax(i)=MAX_EDT_LAND  !
            endif
         enddo
      endif
      !
      !--- minimum depth (m), clouds must have
      !
      if(cumulus == 'deep'                         ) depth_min=1000.
      if(cumulus == 'mid' .or. cumulus == 'shallow') depth_min=500.
      !
      !--- max height(m) above ground where updraft air can originate
      !
      if(cumulus == 'deep'                         ) zkbmax=4000.
      if(cumulus == 'mid' .or. cumulus == 'shallow') zkbmax=3000.
      !
      !--- height(m) above which no downdrafts are allowed to originate
      !
      zcutdown=3000.
      !
      !--- depth(m) over which downdraft detrains all its mass
      !
      z_detr=1000.
      if(cumulus == 'deep'                         ) z_detr= 1000.
      if(cumulus == 'mid' .or. cumulus == 'shallow') z_detr= 300.

      !
      !--- mbdt ~ xmb * timescale
      !
      do i=its,itf
         mbdt(i)= 0.1!*dtime*xmb_nm1(i)
         !mbdt(i)= 100.*(p_cup(i,kbcon(i))-p_cup(i,kbcon(i)+1))/(g*dtime)
         !mbdt(i)= 0.1*mbdt(i)
      enddo
      !
      !--- environmental conditions, FIRST HEIGHTS
      !--- calculate moist static energy, heights, qes
      !
      call cup_env(z ,qes ,he ,hes ,t ,q ,po,z1 ,psur,ierr,-1,itf,ktf,its,ite, kts,kte)
      call cup_env(zo,qeso,heo,heso,tn,qo,po,z1, psur,ierr,-1,itf,ktf,its,ite, kts,kte)

      !
      !--- outputs a model sounding for the stand-alone code (part 1)
      !
      if(OUTPUT_SOUND == 1) then
         call SOUND(1,cumulus,int_time,dtime,ens4,itf,ktf,its,ite, kts,kte,xlats,xlons,jcol,whoami_all  &
            ,z ,qes ,he ,hes ,t ,q ,po,z1 ,psur,zo,qeso,heo,heso,tn,qo,us,vs ,omeg,xz            &
            ,h_sfc_flux,le_sfc_flux,tsur, dx,stochastic_sig,zws,ztexec,zqexec, xland             &
            ,kpbl,k22,klcl,kbcon,ktop,aa0,aa1,sig,xaa0,hkb,xmb,pre,edto                          &
            ,zo_cup,dhdt,rho,zuo,zdo,up_massentro,up_massdetro,outt, outq,outqc,outu,outv)
      endif
      !
      !--- environmental values on cloud levels
      !
      call cup_env_clev(t,qes,q,he,hes,z,po,qes_cup,q_cup,he_cup, &
         us,vs,u_cup,v_cup,                                     &
         hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur,tsur,         &
         ierr,z1,itf,ktf,its,ite, kts,kte)

      call cup_env_clev(tn,qeso,qo,heo,heso,zo,po,qeso_cup,qo_cup, heo_cup,   &
         us,vs,u_cup,v_cup,                                           &
         heso_cup,zo_cup,po_cup,gammao_cup,tn_cup,psur,tsur,          &
         ierr,z1,itf,ktf,its,ite, kts,kte)
      !
      !--- get air density at full layer (model levels) by hydrostatic balance (kg/m3)
      !
      do i=its,itf
         rho_hydr(i,:) = 0.0
         if(ierr(i) /= 0)cycle
         do k=kts,ktf
            rho_hydr(i,k)=100.*(po_cup(i,k)-po_cup(i,k+1))/(zo_cup(i,k+1)-zo_cup(i,k))/g
             !print*,"rhohidr=",k,rho_hydr(i,k),po_cup(i,k+1),zo_cup(i,k+1)
         enddo
      enddo
      !
      !--- partition between liq/ice cloud contents
      !
      call get_partition_liq_ice(ierr,tn,z1,zo_cup,po_cup,p_liq_ice,melting_layer,&
         itf,ktf,its,ite,kts,kte,cumulus)
      !
      do i=its,itf
         if(ierr(i) /= 0) cycle
         do k=kts,ktf
               if(zo_cup(i,k).gt.zkbmax+z1(i))then
                  kbmax(i)=k
                  exit
               endif
          enddo
          !--- level where detrainment for downdraft starts
          do k=kts,ktf
               if(zo_cup(i,k).gt.z_detr+z1(i))then
                  kdet(i)=k
                  exit
               endif
          enddo
      enddo
      !
      !--- determine level with highest moist static energy content - k22
      !
      if(cumulus == 'shallow') then
         start_k22 = 1
      else
         start_k22 = 2
      endif
      k22(:)=kts
      do i=its,itf
         if(ierr(i) /= 0) cycle
         k22(i)=maxloc(HEO_CUP(i,start_k22:kbmax(i)+1),1)+start_k22-1
         k22(i)=max(k22(i),start_k22)
         if(cumulus == 'shallow') then
            k22(i)=min(2,k22(i))

            if(K22(I).gt.KBMAX(i))then
               ierr(i)=2
               ierrc(i)="could not find k22"
            endif

         else

            if(k22(i) > kbmax(i))then
               !- let's try k22=start_k22 for the cases k22>kbmax
               k22(i)= start_k22
               cycle
            endif

         endif
      enddo
      !
      !
      !-- get the pickup of ensemble ave prec, following Neelin et al 2009.
      !
      call precip_cwv_factor(itf,ktf,its,ite,kts,kte,ierr,tn,po,qo,po_cup,cumulus,p_cwv_ave)
      !
      !------- determine LCL for the air parcels around K22
      !
      do i=its,itf
         klcl(i) = k22(i) ! default value
         if(ierr(i) == 0)then
            !tlll, rlll,plll - temp, water vapor and pressure of the source air parcel
            x_add = max(0.,zqexec(i))
            call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),q_cup (i,kts:kte),rlll,k22(i),x_add)
            x_add = max(0.,ztexec(i))
            call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),t_cup (i,kts:kte),tlll,k22(i),x_add)
            call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),p_cup (i,kts:kte),plll,k22(i))

            call get_lcl(tlll,100.*plll,rlll,tlcl,plcl,dzlcl)

            !-get LCL
            if(dzlcl >= 0.) then ! LCL found (if dzlcl<0 => not found)
               call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),z_cup (i,kts:kte),zlll,k22(i))
               loop0:  do k=kts,ktf
                         if(z_cup(i,k).gt.zlll+dzlcl)then
                             klcl(i)=max(k,k22(i))
                             exit loop0
                         endif
               enddo loop0
               klcl(i)=min(klcl(i),ktf-4)
            endif
         endif
         !write(12,111)'MDlcl',tlcl,plcl,dzlcl,klcl(i),ierr(i)
         !111      format(1x,A5,3F10.2,2i4)
      enddo
      !
      !-- check if LCL height is below PBL height to allow shallow convection
      !
      if(lcl_trigger > 0 .and. cumulus == 'shallow')then
         do i=its,itf
            if(ierr(i) /= 0) cycle
            if(klcl(i) > max(1,kpbl(i)-lcl_trigger)) then
                  ierr(i)=21
                  ierrc(i)='for shallow convection:  LCL height < PBL height'
            endif
         enddo
         !print*,"LCL",maxval(klcl),minval(klcl),maxval(kpbl),minval(kpbl)
      endif
      !
      !------- trigger function based on Kain (JAS 2004) 
      !
      if(ADV_TRIGGER==1 .and. cumulus /= 'shallow') then
         wkf = 0.02 ! m/s
         do i=its,itf
            if(ierr(i) /= 0) cycle

            k     = klcl(i)
            dzlcl = z_cup(i,k)-z1(i)
            ckf = wkf
            if(dzlcl .le. 2.e+3) ckf = wkf * dzlcl/2000.
            wkflcl(i) =-(omeg(i,max(kts,k-1),1)/rho(i,max(k-1,kts)) + &
                         omeg(i,k           ,1)/rho(i,k           ) + &
                         omeg(i,k+1         ,1)/rho(i,k+1         ) )/(3.*g)

            !...check to see if cloud is buoyant using fritsch-chappell trigger
            !...function described in kain and fritsch (1992)...w0avg is an
            !...aproximate value for the running-mean grid-scale vertical
            !...velocity, which gives smoother fields of convective initiation
            !...than the instantaneous value...formula relating temperature
            !...perturbation to vertical velocity has been used with the most
            !...success at grid lengths near 25 km.  for different grid-lengths,
            !...adjust vertical velocity to equivalent value for 25 km grid
            !...length, assuming linear dependence of w on grid length...
            if(dx(i) >= 25.e+3) then
               wkflcl(i) = wkflcl(i)*dx(i)/25.e3 - ckf
            else
               wkflcl(i) = wkflcl(i) - ckf
            endif
            !--- think about letting wkflcl <0 => Tpert<0 =>prevent convection in subsidence areas
            ! wkflcl = max(wkflcl,0.) ! -- only positive.
            !
            !-- Kain (2004) Eq. 1           
            Tpert(i,kts   ) = max(0., 4.64*abs(wkflcl(i))**(1./3.)*sign(1.,wkflcl(i)))
            Tpert(i,kts+1:) = Tpert(i,kts)
           
            !DTLCL=4.64*WKL**0.33,  WSIGNE = signal of GDT
            !GDT=G*DTLCL*(ZLCL-Z0(LC))/(TV0(LC)+TVEN)
            !WLCL=1.+.5*WSIGNE*SQRT(ABS(GDT)+1.E-10)    !<< velocity at LCL             
         enddo
         !print*," TPERT=",maxval(Tpert),maxval(wkflcl),maxval(klcl);call flush(6)
      endif
      !
      !--- define entrainment/detrainment profiles for updrafts
      !
      !- initial entrainment/detrainment
      entr_rate   (:  ) = entr_rate_input
      min_entr_rate     = entr_rate_input * 0.1
      !
      !
      !-- cold pool parameterization and convective memory
      !
      if(convection_tracer == 1 .and. cumulus == 'deep') then
            if(USE_MEMORY >= 0) then 
                do i=its,itf
                   if(ierr(i) /= 0) cycle
                   !x_add_buoy(i) = min(mx_buoy2, maxval(buoy_exc(i,kts:klcl(i))))
                
                   call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte), &
                                     buoy_exc (i,kts:kte),x_add_buoy (i),kts)
                                   ! buoy_exc (i,kts:kte),x_add_buoy (i),klcl(i))

                enddo
                !print*,"BU=",maxval(x_add_buoy),minval(x_add_buoy)
            endif 

            !-- avoid extra-buoyancy where rained before
            if(USE_MEMORY == 1 .or. USE_MEMORY == 12) then
               where(AA2_ > 10./3600.) 
                  x_add_buoy = 0.0
                  wlpool     = 0.0
               end where
            endif
             !-- avoid extra-buoyancy where rained before
            if(USE_MEMORY == 4 .or. USE_MEMORY == 14) then
                do i=its,itf
                   if(ierr(i) /= 0) cycle
                   if(AA2_(i) > 1.e-6 .and. x_add_buoy(i) < 1000. .and.  x_add_buoy(i) > 250.) then 
                     x_add_buoy(i) = 0.0
                     wlpool    (i) = 0.0
                     ierr      (i) = 100
                   endif 
                enddo 
               !where(AA2_ > 1.e-6. .and. x_add_buoy < 1000.) 
               !   x_add_buoy = 0.0
               !   wlpool     = 0.0
               !   ierr       = 100
               !end where
            endif
             
            if(USE_MEMORY == 2 .or. USE_MEMORY == 12 .or. USE_MEMORY == 14) then  
                !- initial entrainment/detrainment
                entr_rate(:)  = entr_rate_input       !* 2.0
                min_entr_rate = entr_rate_input * 0.1 !* 2.0
                 
                do i=its,itf !-- reduce entr rate, where cold pools exist
                   if(ierr(i) /= 0) cycle
                  !entr_rate(i) = max(0.1, 1.-coldpool_start(x_add_buoy(i))) * entr_rate(i) 
                  !entr_rate(i) = max(0.5, 1.-coldpool_start(x_add_buoy(i))) * entr_rate(i) 
                   entr_rate(i) = max(0.7, 1.-coldpool_start(x_add_buoy(i))) * entr_rate(i) 
                  !entr_rate(i) = max(0.8, 1.-coldpool_start(x_add_buoy(i))) * entr_rate(i) 
               enddo
                !print*,"ENT",1000.*maxval(entr_rate),1000.*minval(entr_rate)&
                !            ,coldpool_start(maxval((x_add_buoy(:)))),coldpool_start(minval((x_add_buoy(:))))
            endif

            if( USE_MEMORY == 3 .or. ADD_COLDPOOL_CLOS >= 1)  then ! increase capmax
                do i=its,itf
                   if(ierr(i) /= 0) cycle
                   cap_max(i) = cap_max(i) + coldpool_start(x_add_buoy(i)) * 35.
                enddo
            endif
            if( ADD_COLDPOOL_CLOS == 3)  then ! increase x_add_buoy
                do i=its,itf
                   if(ierr(i) /= 0) cycle
                    x_add_buoy(i) = x_add_buoy(i) + 0.5*wlpool(i)**2
                enddo
            endif

            !--- temporary for output
            if(use_gustiness == 0) AA3_(:)=x_add_buoy(:) 
            if(use_gustiness == 4 ) then 
                 AA3_(:)=x_add_buoy(:) ; x_add_buoy(:) = 0.0
            endif
            
            !-- using ztexc and zqexc as perturbation:
            if (use_gustiness == 1 .or. use_gustiness == 2) then 
               AA3_(:)=cp*ztexec(:)+xlv*zqexec(:)
               x_add_buoy(:) = 0. 
            endif 
            
      endif 
      !
      !--- determine the entrainment dependent on environmental moist (here relative humidity)
      !--- also the controls of RH on the diurnal cycle (see Tian et al 2022 GRL)
      if(cumulus == 'deep') &
          call rh_controls(itf,ktf,its,ite,kts,kte,ierr,tn,po,qo,qeso,po_cup,cumulus,rh_entr_factor, &
                           rh_dicycle_fct,entr_rate_input, entr_rate ,xlons,dtime)         
      !
      !
      !--- determine the vertical entrainment/detrainment rates, the level of convective cloud base -kbcon-
      !--- and the scale dependence factor (sig).
      !
      do i=its,itf
         entr_rate_2d(i,:) = entr_rate  (i)
         cd          (i,:) = entr_rate  (i)
         if(ierr(i) /= 0) cycle

         if(cumulus /= 'shallow') then

            do k=kts,ktf

               frh = min(qo_cup(i,k)/qeso_cup(i,k),1.)
               !-------------------------------------------
               if(ENTRNEW) then
                  !- v 2
                  if(k >= klcl(i)) then
                    !entr_rate_2d(i,k)=entr_rate(i)*(1.3-frh)*(qeso_cup(i,k)/qeso_cup(i,klcl(i)))**3
                     entr_rate_2d(i,k)=entr_rate(i)*(1.3-frh)*(qeso_cup(i,k)/qeso_cup(i,klcl(i)))**1.25
                  else
                     entr_rate_2d(i,k)=entr_rate(i)*(1.3-frh)
                  endif
                  cd(i,k)=0.75e-4*(1.6-frh)
                  entr_rate_2d(i,k) = max(entr_rate_2d(i,k),min_entr_rate)
               else
                  !- v 1
                  entr_rate_2d(i,k)=max(entr_rate(i)*(1.3-frh)*max(min(1.,(qeso_cup(i,k)&
                                   /qeso_cup(i,klcl(i)))**1.25),0.1),1.e-5)
                  if(cumulus == 'deep') cd(i,k)=1.e-2*entr_rate(i)
                  if(cumulus == 'mid' ) cd(i,k)=0.75*entr_rate_2d(i,k)
               endif
            enddo

         else

            do k=kts,ktf
               frh = min(qo_cup(i,k)/qeso_cup(i,k),1.)
               !entr_rate_2d(i,k)=entr_rate(i)*(1.3-frh)*max(min(1.,(qeso_cup(i,max(k,klcl(i)))&
               !                                                    /qeso_cup(i,klcl(i)))**3) ,0.1)
               entr_rate_2d(i,k)=entr_rate(i)*(1.3-frh)*max(min(1.,(qeso_cup(i,max(k,klcl(i)))&
                                                                   /qeso_cup(i,klcl(i)))**1) ,0.1)

               ! entr_rate_2d(i,k)=entr_rate(i)*(1.3-frh)*(min(z(i,klcl(i))/z(i,k),1.))
               ! entr_rate_2d(i,k) = max(entr_rate_2d(i,k),min_entr_rate)
               !print*,"ent=",k,real(z(i,k),4),real(min(z(i,klcl(i))/z(i,k),1.),4),real(entr_rate_2d(i,k)*1000.,4)

               cd(i,k)=0.75*entr_rate_2d(i,k)!+0.5e-3
            enddo
         endif
      enddo
      !
      !--- start_level
      !
      start_level(:)=  KLCL(:)
      !start_level(:)=  KTS
      !
      !--- determine the moist static energy of air parcels at source level
      !
      do i=its,itf
         if(ierr(i) /= 0)cycle
         x_add = (xlv*zqexec(i)+cp*ztexec(i)) +  x_add_buoy(i)
         call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),he_cup (i,kts:kte),hkb (i),k22(i),x_add,Tpert(i,kts:kte))
         call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),heo_cup(i,kts:kte),hkbo(i),k22(i),x_add,Tpert(i,kts:kte))
         !print*,"xi=",i, xlv*zqexec(i)+cp*ztexec(i) , x_add_buoy(i), hkbo(i)
      enddo
      !
      !--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
      !
      call cup_cloud_limits(cumulus,ierrc,ierr,cap_max_increment,cap_max,heo_cup,heso_cup,qo_cup   &
         ,qeso_cup,po,po_cup,zo_cup,heo,hkbo,qo,qeso,entr_rate_2d,hcot,k22,kbmax &
         ,klcl,kbcon,ktop,depth_neg_buoy,frh_bcon,Tpert,start_level              &
         ,use_excess,zqexec,ztexec,x_add_buoy,xland,itf,ktf,its,ite, kts,kte)
      !
      !--- SCALE DEPENDENCE FACTOR (SIG), version new
      !
      if( USE_SCALE_DEP == 0  .or.  cumulus == 'shallow') then
         sig(:)=1.
      else
         do i=its,itf
            sig(i) = 0.
            if(ierr(i) /= 0) cycle
            !--original
            !sig(i) = 1.0-0.9839*exp(-0.09835.  *(dx(i)/1000.))

            !-- for similar curve as in IFS/EC, use sig_factor = 0.22
            sig(i) = 1.0 - exp(-sig_factor*(dx(i)/1000.))

            !print*,"sig=",sig(i),dx(i), sig_factor
            if (stochastic_sig(i) /= 1.0) then
               sig(i) = sig(i)**(stochastic_sig(i)*MAX(0.9,0.9*sig(i)))
            endif
            sig(i)= max(0.001,min(sig(i),1.))
         enddo
         !print*,'sig',maxval(sig),minval(sig),maxval(dx),minval(dx)
      endif
      !
      !--- define entrainment/detrainment profiles for downdrafts
      !
      if(ENTRNEW) then
         mentrd_rate(:) = entr_rate(:)*0.3
      else
         mentrd_rate(:) = entr_rate(:)
      endif
      do i=its,itf
         cdd(i,kts:kte) = mentrd_rate(i)
      enddo
      !- scale dependence factor
      sigd(:)  = 1.
      if( DOWNDRAFT == 0) sigd(:) = 0.0
      !
      !--- update hkb/hkbo in case of k22 is redefined in 'cup_kbon'
      !
      do i=its,itf
         if(ierr(i) /= 0) cycle
         x_add = (xlv*zqexec(i)+cp*ztexec(i)) +  x_add_buoy(i)
         call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),he_cup (i,kts:kte),hkb (i),k22(i),x_add,Tpert(i,kts:kte))
         call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),heo_cup(i,kts:kte),hkbo(i),k22(i),x_add,Tpert(i,kts:kte))
      enddo
      !
      !--- increase detrainment in stable layers
      !
      call cup_minimi(HEso_cup,Kbcon,kstabm,kstabi,ierr,itf,ktf,its,ite, kts,kte)
      !
      !--- option for using the inversion layers as a barrier for the convection development
      !
      if(cumulus == 'mid') then
         if(USE_INV_LAYERS) then
            !
            !--- get inversion layers
            call get_inversion_layers(cumulus,ierr,psur,po_cup,tn_cup,zo_cup,k_inv_layers,&
               dtempdz,itf,ktf,its,ite, kts,kte)
            do i=its,itf
               if(ierr(i) /= 0)cycle
               ktop(i) = min(ktop(i),k_inv_layers(i,mid))
              !print*,"ktop=",ktop(i),k_inv_layers(i,mid)
             enddo
         endif
         !
         !-- check if ktop is above 450hPa layer for mid convection
         !
         do i=its,itf
            if(ierr(i) /= 0)cycle
            !print*,"sta=",Kbcon(i),kstabm(i),kstabi(i),p_cup(i,ktop(i)),z_cup(i,kstabi(i))
            if(po_cup(i,ktop(i)) < 450.) then
               ierr(i)=25
               ierrc(i)='mid convection with cloud top above 450 hPa (~ 7km asl)'
            endif
         enddo
         !-- check if ktop is below 750hPa layer for mid convection
         !
         do i=its,itf
            if(ierr(i) /= 0)cycle
            if(po_cup(i,ktop(i)) > 750.) then
               ierr(i)=55
               ierrc(i)='ktop too low for mid'
            endif
         enddo
      endif

      if(cumulus == 'shallow') then
         if(USE_INV_LAYERS) then
            call get_inversion_layers(cumulus,ierr,psur,po_cup,tn_cup,zo_cup,k_inv_layers,&
               dtempdz,itf,ktf,its,ite, kts,kte)
            do i=its,itf
               if(ierr(i) /= 0)cycle
               ktop(i) = min(ktop(i),k_inv_layers(i,shal))

            enddo
         endif

         !--- Check if ktop is above 700hPa layer for shallow convection
         do i=its,itf
            if(ierr(i) /= 0)cycle
            min_shall_top=700.
            !if(icumulus_gf(mid) == 0) min_shall_top=500.
            if(po_cup(i,ktop(i)) < min_shall_top) then
               ierr(i)=26
               ierrc(i)='shallow convection wit h cloud top above min_shall_top hPa' 
            endif
         enddo
      endif

      !
      !-- last checks for ktop (deep and mid)
      !
      if(zero_diff == 1) then
         do i=its,itf
            if(ierr(i) /= 0)cycle
            if(cumulus /= 'shallow') then
               k=5
            else
               k=2
            endif
            if(ktop(i) < kbcon(i)+k)then
               ierr(i)=5
               ierrc(i)='ktop too low'
            endif
         enddo
      !      ELSE
      !       DO i=its,itf
      !        if(ierr(i) /= 0)cycle
      !        if( (z_cup(i,ktop(i))-z_cup(i,kbcon(i))) < depth_min ) then
      !            ierr(i)=5
      !            ierrc(i)='cloud depth too small'
      !        endif
      !       ENDDO
      endif

      do i=its,itf
         if(ktop(i) <= kbcon(i))then
            ierr(i)=5
            ierrc(i)='ktop too small'
         endif
      enddo

      !---not-zero-diff-APR-06-2020
      if(zero_diff == 0) then
         if(cumulus == 'deep') then
            min_deep_top=500.
            if(icumulus_gf(mid) == 0) min_deep_top=750.
            do i=its,itf
               if(ierr(i) /= 0)cycle
               if(po_cup(i,ktop(i)) > min_deep_top) then
                  ierr(i)=55
                  ierrc(i)='ktop too low for deep'
               endif
            enddo
         endif
      endif
      !---not-zero-diff-APR-06-2020


      !
      !-- avoid double-counting with shallow scheme (deep and mid)
      !
      do i=its,itf
         if(ierr(i) /= 0)cycle
         if(last_ierr(i) == 0) then
             !--- if 'mid' => last was 'shallow'
            ! if(cumulus == 'mid' .and. po_cup(i,ktop(i)) > 700.) then
            !   ierr(i)=27
            !   ierrc(i)='avoiding double-counting shallow and mid'
            ! endif
            !--- if 'mid' => last was 'shallow'
            if(cumulus == 'mid' ) then
               ierr(i)=27
               ierrc(i)='avoiding double-counting deep and mid'
            endif
         endif
      enddo
      !
      !--- determine the normalized mass flux profile for updraft
      !
      do i=its,itf
         zuo(i,:)=0.
         if(ierr(i) /= 0) cycle
         call get_zu_zd_pdf(trim(cumulus),trim(cumulus)//"_up",ierr(i),k22(i),ktop(i),zuo(i,kts:kte),kts,kte,ktf  &
            ,kpbl(i),k22(i),kbcon(i),klcl(i),po_cup(i,kts:kte),psur(i),xland(i),random(i))
      enddo

      do i=its,itf
         if(ierr(i) /= 0) cycle
         xzu(i,:)= zuo(i,:)
         zu (i,:)= zuo(i,:)
      enddo
      !
      ! calculate mass entrainment and detrainment
      !
      call get_lateral_massflux(itf,ktf, its,ite, kts,kte,min_entr_rate              &
         ,ierr,ktop,zo_cup,zuo,cd,entr_rate_2d,po_cup          &
         ,up_massentro, up_massdetro ,up_massentr, up_massdetr &
         ,cumulus,kbcon,k22,kpbl,up_massentru,up_massdetru,lambau_dp)
      uc  =0.
      vc  =0.
      hc  =0.
      hco =0.
      do i=its,itf
         if(ierr(i).eq.0)then
            do k=kts,start_level(i)
               hc (i,k) =hkb (i)
               hco(i,k) =hkbo(i)
               !-get uc and vc as average between layers below k22
               call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),u_cup(i,kts:kte),uc(i,k),k22(i))
               call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),v_cup(i,kts:kte),vc(i,k),k22(i))
            enddo
         endif
      enddo
      !
      !--- 1st guess for moist static energy and dbyo (not including ice phase)
      !
      do i=its,itf
         if(ierr(i) /= 0) cycle
         do k=start_level(i)  +1,ktop(i) + 1  ! mass cons option
            denom=(zu(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))
            if(denom > 0.0)then
               hco(i,k)=(hco(i,k-1)*zuo(i,k-1)-.5*up_massdetro(i,k-1)*hco(i,k-1)+ &
                  up_massentro(i,k-1)*heo(i,k-1))  / denom
               if(k==start_level(i)+1) then
                  x_add = (xlv*zqexec(i)+cp*ztexec(i)) +  x_add_buoy(i)
                  hco(i,k)= hco(i,k) + x_add*up_massentro(i,k-1)/denom
               endif
            else
               hco(i,k)= hco(i,k-1)
            endif
         enddo
         do k=ktop(i)+2,ktf
            hco (i,k)=heso_cup(i,k)!=heo_cup(i,k)
         enddo
      enddo
      !
      !--- Get buoyancy of updrafts
      !
      call get_buoyancy(itf,ktf, its,ite, kts,kte,ierr,klcl,kbcon,ktop,hco,heo_cup,heso_cup,dbyo,zo_cup)

      !--- get "c1d" profile ----------------------------------------
      if(cumulus == 'deep' .and. USE_C1D) then
         do i=its,itf
            if(ierr(i) /= 0) cycle
            c1d(i,kbcon(i)+1:ktop(i)-1)=abs(c1)
         enddo
      endif

      if(FIRST_GUESS_W .or. AUTOCONV == 4) then
         call cup_up_moisture_light(cumulus,start_level,klcl,ierr,ierrc,zo_cup,qco,qrco,pwo,pwavo,hco,tempco,xland   &
            ,po,p_cup,kbcon,ktop,cd,dbyo,clw_all,t_cup,qo,GAMMAo_cup,zuo          &
            ,qeso_cup,k22,qo_cup,ZQEXEC,use_excess,rho,up_massentr,up_massdetr    &
            ,psum,psumh,c1d,x_add_buoy,1,itf,ktf,ipr,jpr,its,ite, kts,kte         )

         call cup_up_vvel(vvel2d,vvel1d,zws,entr_rate_2d,cd,zo,zo_cup,zuo,dbyo,GAMMAo_CUP,tn_cup &
            ,tempco,qco,qrco,qo,start_level,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte,wlpool,wlpool_bcon,1)
      endif

      !
      !--- calculate moisture properties of updraft
      call cup_up_moisture(cumulus,start_level,klcl,ierr,ierrc,zo_cup,qco,qrco,pwo,pwavo,hco,tempco,xland   &
         ,po,p_cup,kbcon,ktop,cd,dbyo,clw_all,t_cup,qo,GAMMAo_cup,zuo,qeso_cup &
         ,k22,qo_cup,ZQEXEC,use_excess,ccn,rho,up_massentr,up_massdetr,psum    &
         ,psumh,c1d,x_add_buoy,vvel2d,vvel1d,zws,entr_rate_2d                  &
         ,1,itf,ktf,ipr,jpr,its,ite, kts,kte                                   )

      do i=its,itf
         if(ierr(i) /= 0) cycle
         cupclw(i,kts:ktop(i)+1)=qrco(i,kts:ktop(i)+1)
      enddo


      !
      !--- get melting profile
      !
      call get_melting_profile(ierr,tn_cup,po_cup, p_liq_ice,melting_layer,qrco    &
         ,pwo,edto,pwdo,melting                               &
         ,itf,ktf,its,ite, kts,kte, cumulus                   )

      !
      !--- updraft moist static energy + momentum budget
      !
      !--- option to produce linear fluxes in the sub-cloud layer.
      if(cumulus == 'shallow' .and. use_linear_subcl_mf == 1) then
         do i=its,itf
            if(ierr(i) /= 0) cycle
            call get_delmix(cumulus,kts,kte,ktf,xland(i),start_level(i),po(i,kts:kte) &
               ,he_cup (i,kts:kte), hc (i,kts:kte))
            call get_delmix(cumulus,kts,kte,ktf,xland(i),start_level(i),po(i,kts:kte) &
               ,heo_cup(i,kts:kte), hco(i,kts:kte))
         enddo
      endif

      do i=its,itf
         if(ierr(i) /= 0) cycle

         do k=start_level(i)+1 , ktop(i)+1  ! mass cons option
            denom =(zu(i,k-1)-.5*up_massdetr (i,k-1)+up_massentr (i,k-1))
            denomU=(zu(i,k-1)-.5*up_massdetru(i,k-1)+up_massentru(i,k-1))
            if(denom > 0.0 .and. denomU >0.0)then

               hc (i,k)=(hc (i,k-1)*zu (i,k-1)-.5*up_massdetr (i,k-1)*hc (i,k-1) + &
                  up_massentr (i,k-1)*he (i,k-1))/ denom

               hco(i,k)=(hco(i,k-1)*zuo(i,k-1)-.5*up_massdetro(i,k-1)*hco(i,k-1)+ &
                  up_massentro(i,k-1)*heo(i,k-1))/ denom
               if(k==start_level(i)+1) then
                  x_add = (xlv*zqexec(i)+cp*ztexec(i)) +  x_add_buoy(i)
                  hco(i,k)= hco(i,k) + x_add*up_massentro(i,k-1)/denom
                  hc (i,k)= hc (i,k) + x_add*up_massentr (i,k-1)/denom
               endif
               !assuming zuo=zu,up_massdetro=up_massdetr, ...
               !(zuo(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))

               uc(i,k)=(uc(i,k-1)*zu(i,k-1)-.5*up_massdetru(i,k-1)*uc(i,k-1)+     &
                  up_massentru(i,k-1)*us(i,k-1)      &
                  -pgcon*.5*(zu(i,k)+zu(i,k-1))*(u_cup(i,k)-u_cup(i,k-1))) /  denomU

               vc(i,k)=(vc(i,k-1)*zu(i,k-1)-.5*up_massdetru(i,k-1)*vc(i,k-1)+     &
                  up_massentru(i,k-1)*vs(i,k-1)      &
                  -pgcon*.5*(zu(i,k)+zu(i,k-1))*(v_cup(i,k)-v_cup(i,k-1))) /  denomU

            else
               hc (i,k)= hc (i,k-1)
               hco(i,k)= hco(i,k-1)
               uc (i,k)= uc (i,k-1)
               vc (i,k)= vc (i,k-1)
            endif
            !---meltglac-------------------------------------------------
            !- includes glaciation effects on HC,HCO
            !                    ------ ice content --------
            !print*,"H=",hc (i,k),(1.-p_liq_ice(i,k))*qrco(i,k)*xlf,hc (i,k)+(1.-p_liq_ice(i,k))*qrco(i,k)*xlf
            hc (i,k)= hc (i,k)+(1.-p_liq_ice(i,k))*qrco(i,k)*xlf
            hco(i,k)= hco(i,k)+(1.-p_liq_ice(i,k))*qrco(i,k)*xlf

         enddo
         !
         do k=ktop(i)+2,ktf
            hc  (i,k)= hes_cup(i,k)!= he_cup(i,k)
            uc  (i,k)=   u_cup(i,k)
            vc  (i,k)=   v_cup(i,k)
            hco (i,k)=heso_cup(i,k)!=heo_cup(i,k)
            zu  (i,k)=0.
            zuo (i,k)=0.
         enddo
      enddo
      !
      !--- Get buoyancy of updrafts
      !
      call get_buoyancy(itf,ktf, its,ite, kts,kte,ierr,klcl,kbcon,ktop,hc, he_cup, hes_cup,  dby, z_cup)
      call get_buoyancy(itf,ktf, its,ite, kts,kte,ierr,klcl,kbcon,ktop,hco,heo_cup,heso_cup,dbyo,zo_cup)
      !---
      !
      if(convection_tracer == 1 .and. SGS_W_TIMESCALE == 1 .and. cumulus == 'deep') then 
            !--- compute vertical velocity
            !
            !call cup_up_vvel(vvel2d,vvel1d,zws,entr_rate_2d,cd,zo,zo_cup,zuo,dbyo,GAMMAo_CUP,tn_cup &
            !                ,tempco,qco,qrco,qo,start_level,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte&
            !                ,wlpool,wlpool_bcon,2)
             wlpool_bcon(:) = wlpool(:)
            !--- trigger function based on KE > CIN
            if( add_coldpool_clos == 2 )then
               call cup_up_aa0(cin1,zo_cup,zuo,dbyo ,GAMMAo_CUP   ,tn_cup  ,k22,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte,'CIN')
               do i=its,itf
                   if (ierr(i) /= 0) cycle
                   ke_mx = 0.5 * max( wlpool_bcon(i)**2, zws(i)**2) + 1.e-6
                   if (ke_mx < abs( min(cin1(i), 0.) )) ierr(i)=500
               enddo
            endif 
      endif
      !
      if(.not. FIRST_GUESS_W) then
         !--- calculate in-cloud/updraft air temperature for vertical velocity
         !
         do i=its,itf
            if(ierr(i) == 0)then
               do k=kts,ktf
                  tempco (i,k) = (1./cp)*(hco (i,k)-g*zo_cup(i,k)-xlv*qco (i,k))
               enddo
               tempco (i,kte)=tn_cup(i,kte)
            else
               tempco (i,:)  =tn_cup(i,:)
            endif
         enddo
         !
         !--- vertical velocity
         !
         call cup_up_vvel(vvel2d,vvel1d,zws,entr_rate_2d,cd,zo,zo_cup,zuo,dbyo,GAMMAo_CUP,tn_cup &
            ,tempco,qco,qrco,qo,start_level,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte,wlpool,wlpool_bcon,1)
      endif

      !---- new rain
      !
      !--- calculate rain mixing ratio in updrafts
      !
      !       call cup_up_rain(cumulus,klcl,kbcon,ktop,k22,ierr,xland         &
      !                       ,zo_cup,qco,qrco,pwo,pwavo,po,p_cup,t_cup,tempco&
      !                       ,zuo,up_massentr,up_massdetr,vvel2d,rho         &
      !                       ,qrr                                            &
      !                       ,itf,ktf,its,ite, kts,kte)

      !---- new rain
      !
      !
      !--- DOWNDRAFT section
      !
      do i=its,itf
         kzdown(i)=0
         if(ierr(i).eq.0)then
            zktop=(zo_cup(i,ktop(i))-z1(i))*.6
            zktop=min(zktop+z1(i),zcutdown+z1(i))
            do k=kts,ktf
               if(zo_cup(i,k).gt.zktop)then
                  kzdown(i)=k
                  go to 37
               endif
            enddo
         endif
37    continue
      enddo
      !
      !--- DOWNDRAFT ORIGINATING LEVEL - JMIN
      !
      call cup_minimi(heso_cup,k22,kzdown,jmin,ierr,itf,ktf,its,ite, kts,kte)

      call get_jmin(cumulus,itf,ktf,its,ite, kts,kte,ierr,kdet,ktop,kbcon,jmin,ierrc  &
         ,beta,depth_min,heso_cup,zo_cup,melting_layer)
      !
      !--- this calls routine to get downdrafts normalized mass flux
      !
      do i=its,itf
         zd(i,:)=0.
         if(ierr(i) /= 0) cycle
         call get_zu_zd_pdf(trim(cumulus),"DOWN",ierr(i),kdet(i),jmin(i),zdo(i,:),kts,kte,ktf&
            ,kpbl(i),k22(i),kbcon(i),klcl(i),po_cup(i,kts:kte),psur(i),xland(i),random(i))
      enddo
      !
      !---  calls routine to get lateral mass fluxes associated with downdrafts
      !
      call get_lateral_massflux_down(trim(cumulus),itf,ktf, its,ite, kts,kte &
         ,ierr,jmin,zo_cup,zdo,xzd,zd,cdd,mentrd_rate_2d      &
         ,dd_massentro,dd_massdetro ,dd_massentr, dd_massdetr &
         ,cumulus,mentrd_rate,dd_massentru,dd_massdetru,lambau_dn)
      !
      !---  calls routine to get wet bulb temperature and moisture at jmin
      !
      if(USE_WETBULB == 1 .and. cumulus /= 'shallow' ) then
         do i=its,itf
            if(ierr(i) /= 0) cycle
            k = jmin(i)
            call get_wetbulb(jmin(i),qo_cup(i,k),t_cup(i,k),po_cup(i,k),q_wetbulb(i),t_wetbulb(i))
          !print*,"wb       =",jmin,qo_cup(i,k),t_cup(i,k),q_wetbulb(i),t_wetbulb(i)
          !print*,"evap/cool=",q_wetbulb(i)-qo_cup(i,k),t_wetbulb(i)-t_cup(i,k)
         enddo
      endif
      !
      !
      !--- downdraft moist static energy + moisture budget
      !
      do i=its,itf
         hcdo (i,:)= heso_cup(i,:)
         ucd  (i,:)=    u_cup(i,:)
         vcd  (i,:)=    v_cup(i,:)
         dbydo(i,:)= 0.
      enddo

      do i=its,itf
         bud(i)=0.
         if(ierr(i)/= 0 .or. cumulus == 'shallow')cycle
         i_wb=0
         !--for future test
         if(use_wetbulb==1) then
            !--option 1
            !hcdo(i,jmin(i))=cp*t_wetbulb(i)+xlv*q_wetbulb(i)+zo_cup(i,jmin(i))*g
            !--option 2
            hcdo(i,jmin(i))=0.5*(cp*t_wetbulb(i)+xlv*q_wetbulb(i)+zo_cup(i,jmin(i))*g + hc(i,jmin(i)))
            i_wb=1
         endif

         dbydo(i,jmin(i))=hcdo(i,jmin(i))-heso_cup(i,jmin(i))
         bud(i)=dbydo(i,jmin(i))*(zo_cup(i,jmin(i)+1)-zo_cup(i,jmin(i)))

         do ki=jmin(i) - i_wb ,kts,-1!do ki=jmin(i)-1,1,-1
            denom = zdo(i,ki+1)-0.5*dd_massdetro(i,ki)+dd_massentro(i,ki)
            denomU= zdo(i,ki+1)-0.5*dd_massdetru(i,ki)+dd_massentru(i,ki)
             !-tmp fix for denominator being zero
            if(denom > 0.0 .and. denomU >0.0)then
               dzo=zo_cup(i,ki+1)-zo_cup(i,ki)

               ucd(i,ki)=(ucd(i,ki+1)*zdo(i,ki+1)-.5*dd_massdetru(i,ki)*ucd(i,ki+1)+ &
                  dd_massentru(i,ki)*us (i,ki)    &
                  -pgcon*zdo(i,ki+1)*(us(i,ki+1)-us(i,ki)))   /  denomU
               vcd(i,ki)=(vcd(i,ki+1)*zdo(i,ki+1)-.5*dd_massdetru(i,ki)*vcd(i,ki+1)+        &
                  dd_massentru(i,ki)*vs (i,ki)           &
                  -pgcon*zdo(i,ki+1)*(vs(i,ki+1)-vs(i,ki)))   /  denomU

               hcdo(i,ki)=(hcdo(i,ki+1)*zdo(i,ki+1)-.5*dd_massdetro(i,ki)*hcdo(i,ki+1)+     &
                  dd_massentro(i,ki)*heo(i,ki))  /denom

               dbydo(i,ki)=hcdo(i,ki)-heso_cup(i,ki)
               !if(i.eq.ipr)write(0,*)'ki,bud = ',ki,bud(i),hcdo(i,ki)
               bud(i)=bud(i)+dbydo(i,ki)*dzo
            else
               ucd (i,ki)= ucd(i,ki+1)
               vcd (i,ki)= vcd(i,ki+1)
               hcdo(i,ki)=hcdo(i,ki+1)
            endif
         enddo
         if(bud(i).gt.0)then
            ierr(i)=7
            ierrc(i)='downdraft is not negatively buoyant '
         endif
      enddo
      !
      !--- calculate moisture properties of downdraft
      !
      call cup_dd_moisture(cumulus,ierrc,zdo,hcdo,heso_cup,qcdo,qeso_cup,     &
         pwdo,qo_cup,zo_cup,dd_massentro,dd_massdetro,jmin,ierr,gammao_cup, &
         !--test    pwevo,bu,qrcdo,qo,heo,t_cup,1,t_wetbulb,q_wetbulb,qco,pwavo,       &
         pwevo,bu,qrcdo,qo,heo,tn_cup,1,t_wetbulb,q_wetbulb,qco,pwavo,      &
         itf,ktf,its,ite, kts,kte)
      !
      !
      !--- calculate workfunctions for updrafts
      !
      call cup_up_aa0(aa0,z_cup ,zu ,dby  ,GAMMA_CUP    ,t_cup   ,k22,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte)
      call cup_up_aa0(aa1,zo_cup,zuo,dbyo ,GAMMAo_CUP   ,tn_cup  ,k22,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte)

      do i=its,itf
         if(ierr(i) /= 0) cycle
         if(aa1(i).eq.0.)then
            ierr(i)=17
            ierrc(i)="cloud work function zero"
         endif
      enddo

      !
      !--- Implements Becker et al (2021) closure, part 1
      !
      if ( (DICYCLE==2 .or. DICYCLE==3) .and. cumulus == 'deep') then

        do ki = 1,2
         if(DICYCLE==2 .and. ki==2) cycle 

         if(ki==1) then
          !-- get the cloud work function for updrafts associated only with RAD + PBL
          tn_x = t + tn - tn_adv
          qo_x = q + qo - qo_adv

          !-- to check => aa1_radpbl=aa1
          !! tn_x = tn
          !! qo_x = qo
         endif

         if(ki==2) then
            !-- get the cloud work function for updrafts associated only with Qv-advection
           !tn_x = tn_adv  ! orig 
           tn_x = t        ! v2
           qo_x = qo_adv
         endif

         ierr_dummy=ierr

         call cup_env(zo,qeso_x,heo_x,heso_x,tn_x,qo_x,po,z1,psur,ierr_dummy,-1,itf,ktf, its,ite, kts,kte)
         call cup_env_clev(tn_x,qeso_x,qo_x,heo_x,heso_x,zo,po,qeso_cup_x,qo_cup_x,heo_cup_x,us,vs   &
            ,u_cup_x,v_cup_x,heso_cup_x,zo_cup_x,po_cup_x,gammao_cup_x,tn_cup_x,psur,tsur  &
            ,ierr_dummy,z1,itf,ktf,its,ite, kts,kte)

         !--- get MSE
         do i=its,itf
            if(ierr_dummy(i) /= 0) cycle
            x_add = (xlv*zqexec(i)+cp*ztexec(i)) +  x_add_buoy(i)
            call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),heo_cup_x(i,kts:kte),hkbo_x(i),k22(i),x_add,Tpert(i,kts:kte))
            hco_x (i,kts:start_level(i)) = hkbo_x(i)

            do k=start_level(i)+1 , ktop(i)+1  ! mass cons option
               denom =(zu(i,k-1)-.5*up_massdetr (i,k-1)+up_massentr (i,k-1))
               if(denom > 0.0 )then

                  hco_x(i,k)=(hco_x(i,k-1)*zuo(i,k-1)-.5*up_massdetro(i,k-1)*hco_x(i,k-1)+ &
                     up_massentro(i,k-1)*heo_x(i,k-1))/ denom
                  if(k==start_level(i)+1) then
                     x_add = (xlv*zqexec(i)+cp*ztexec(i)) +  x_add_buoy(i)
                     hco_x(i,k)= hco_x(i,k) + x_add*up_massentro(i,k-1)/denom
                  endif
               else
                  hco_x(i,k)= hco_x(i,k-1)
               endif
               !- includes glaciation effects on HCO_X
               hco_x(i,k)= hco_x(i,k)+(1.-p_liq_ice(i,k))*qrco(i,k)*xlf
            enddo
            hco_x (i,ktop(i)+2:ktf)=heso_cup_x(i,ktop(i)+2:ktf)
         enddo

         call get_buoyancy(itf,ktf, its,ite, kts,kte,ierr_dummy,klcl,kbcon,ktop &
            ,hco_x,heo_cup_x,heso_cup_x,dbyo_x,zo_cup_x)

         if(ki==1) then  ! RAD+PBL only
            call cup_up_aa0(aa1_radpbl,zo_cup_x,zuo,dbyo_x,GAMMAo_CUP_x,tn_cup_x   &
                           ,k22,klcl,kbcon,ktop,ierr_dummy,itf,ktf,its,ite, kts,kte)
            !-- get AA1_ADV
            !aa1_adv = aa1 + aa0 - aa1_radpbl

         endif

         if(ki==2) & ! ADV of Qv only
         call cup_up_aa0(aa1_adv,zo_cup_x,zuo,dbyo_x,GAMMAo_CUP_x,tn_cup_x   &
            ,k22,klcl,kbcon,ktop,ierr_dummy,itf,ktf,its,ite, kts,kte)
        
         !Observe that : 
         !aa1 ~ aa0 + (aa1_radpbl-aa0) + (aa1_adv-aa0)

        enddo ! ki
      endif
      !
      !--- calculate CIN for updrafts
      !
      ! call cup_up_aa0(cin0,z_cup ,zu ,dby  ,GAMMA_CUP    ,t_cup   ,k22,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte,'CIN')
      ! call cup_up_aa0(cin1,zo_cup,zuo,dbyo ,GAMMAo_CUP   ,tn_cup  ,k22,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte,'CIN')
      !
      !----trigger function for KE+CIN < 0 => no convection
      !
      ! IF( DICYCLE>1) THEN
      !   do i=its,itf
      !     if(ierr(i) /= 0) cycle
      !     print*,"cin=",cin0(i),0.5*zws(i)**2, omeg(i,kpbl(i),1)/(-g*rho(i,kpbl(i)))
      !call flush(6)
      !     if(cin0(i) + 0.5*zws(i)**2 < 0.)then !think about including the grid scale vertical velocity at KE calculation
      !          ierr(i)=19
      !          ierrc(i)="CIN negat"
      !     endif
      !    enddo
      !  ENDIF
      !
      !
      !--- calculate in-cloud/updraft and downdraft air temperature for vertical velocity
      !
      do i=its,itf
         if(ierr(i) == 0)then
            do k=kts,ktf
               tempcdo(i,k) = (1./cp)*(hcdo(i,k)-g*zo_cup(i,k)-xlv*qcdo(i,k))
            enddo
         else
            tempcdo(i,:)=tn_cup(i,:)
         endif
      enddo
      !
      !
      !--- diurnal cycle section
      !
      !--- Bechtold et al 2008 time-scale of cape removal
      !
      if(cumulus=='deep') then
               tau_ecmwf(:)=tau_deep;   wmean(:) = 3. !  mean vertical velocity m/s
      else
               tau_ecmwf(:)=tau_mid ;   wmean(:) = 3. 
      endif
      !--- we shall let all scale dependence on the sig parameter
      !tau_ecmwf(:)= tau_ecmwf(:) * (1. + 1.66 * (dx(:)/(125*1000.)))! dx must be in meters

      if(SGS_W_TIMESCALE == 1 .and. cumulus=='deep') then

         do i=its,itf
            if(ierr(i) /= 0) cycle
            !- mean vertical velocity based on integration of vertical veloc equation
            wmean(i) = min(max(vvel1d(i),3.),20.)

            !- time-scale cape removal from Bechtold et al. 2008
            tau_ecmwf(i)=( zo_cup(i,ktop(i))- zo_cup(i,kbcon(i)) ) / wmean(i)
            !tau_ecmwf(i)= min(10800., max(720.,tau_ecmwf(i)))
            tau_ecmwf(i)= min(10800., max(1000.,tau_ecmwf(i)))
         enddo
!====
!  do i=its,itf
!      if(ierr(i) /= 0) cycle
!      print*,'tauec',wlpool_bcon(i),vvel1d(i), wmean(i) ,tau_ecmwf(i)
!      call flush(6)
!  enddo
!====
      endif

      !
      !--- Implements the Bechtold et al (2014) and Becker et al (2021) closures
      !
      do i=its,itf
         if(ierr(i) /= 0) cycle
         !- over water
         !   umean= 2.0+sqrt(0.5*(US(i,1)**2+VS(i,1)**2+US(i,kbcon(i))**2+VS(i,kbcon(i))**2))
         !   tau_bl(i) = (zo_cup(i,kbcon(i))- z1(i)) /umean
         !- over land
         !   tau_bl(i) = tau_ecmwf(i)
         !-----------
         umean= 2.0+sqrt(0.5*(US(i,1)**2+VS(i,1)**2+US(i,kbcon(i))**2+VS(i,kbcon(i))**2))
         !--                    - over land -            -          over ocean       s   -
         tau_bl(i)=(1.-xland(i))*tau_ecmwf(i) + xland(i)*(zo_cup(i,kbcon(i))- z1(i)) /umean
      enddo

      if( dicycle<=3 .and. cumulus == 'deep') then

         !-- calculate "pcape" or equivalent cloud work function from the BL forcing only
         iversion=0
         call cup_up_aa1bl(iversion,aa1_bl,aa1_fa,aa1,t,tn,q,qo,dtime,po_cup,zo_cup,zuo,dbyo,GAMMAo_CUP,tn_cup, &
            rho,klcl,kpbl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte,                          &
            xland,ztexec,xlons,xlats, h_sfc_flux,le_sfc_flux,tau_bl,tau_ecmwf,t_star,cumulus,tn_bl,qo_bl  )

         do i=its,itf
            if(ierr(i) /= 0) cycle
            aa1_bl(i) = (aa1_bl(i)/T_star) * tau_bl(i) ! units J/kg
           !aa1_bl(i) = (aa1_bl(i)/T_star) * tau_bl(i) - cin1(i)
            aa1_bl(i) = min(2000., abs(aa1_bl(i)))*sign(1.,aa1_bl(i))
         enddo
         !
         !--- Adds Becker et al (2021) closure, part 2
         !
         if(dicycle==2) then
            call get_Qadv(cumulus,itf,ktf,its,ite,kts,kte,ierr,dtime,q,qo,qo_adv,po,po_cup &
               ,qeso, q_adv,col_sat_adv,alpha_adv,tau_bl,zo_cup,kbcon,ktop)
         endif
      
         if(dicycle==3) then
            do i=its,itf
                if(ierr(i) /= 0) cycle
                aa1_adv(i) =  (aa1_adv(i) - aa0(i)) * tau_bl(i)/dtime
            enddo
         endif

      !
      !--- Implements the Zhang(2002) closure
      !
      elseif( dicycle==4 .and. cumulus == 'deep') then

         !- T and Q profiles modified only by RAD+ADV tendencies
         do i=its,itf
            if ( ierr(i) /= 0 ) cycle
            tn_x(i,kts:ktf) = tn(i,kts:ktf)-tn_bl(i,kts:ktf)+t(i,kts:ktf)
            qo_x(i,kts:ktf) = qo(i,kts:ktf)-qo_bl(i,kts:ktf)+q(i,kts:ktf)
         enddo

         !--- calculate moist static energy, heights, qes, ... only by free troposphere tendencies
         call cup_env(zo,qeso_x,heo_x,heso_x,tn_x,qo_x,po,z1, &
            psur,ierr,-1,itf,ktf, its,ite, kts,kte)
         !--- environmental values on cloud levels only by FT tendencies
         call cup_env_clev(tn_x,qeso_x,qo_x,heo_x,heso_x,zo,po,qeso_cup_x,qo_cup_x,heo_cup_x,           &
            us,vs,u_cup,v_cup,                                    &
            heso_cup_x,zo_cup,po_cup,gammao_cup_x,tn_cup_x,psur,tsur,  &
            ierr,z1,itf,ktf,its,ite, kts,kte)
         !--- this is (DT_ve/Dt)_adv+rad
         do i=its,itf
            if(ierr(i) /= 0) cycle
            aa3(i)=0.
            do k=max(kbcon(i),kts+1),ktop(i)
               dp=- (log(100.*po(i,k))-log(100.*po(i,k-1))) !no units
               aa3(i)=aa3(i)- (tn_cup_x(i,k)*(1.+0.608*qo_cup_x(i,k)) - &
                  t_cup   (i,k)*(1.+0.608*q_cup   (i,k)) ) *dp/dtime !units = K/s
               !print*,"tve=",k,aa3(i),tn_cup_x(i,k)*(1.+0.608*qo_cup_x(i,k)),&
               !               t_cup (i,k)*(1.+0.608*q_cup   (i,k)),dp
            enddo
         enddo
         do i=its,itf
            if(ierr(i)/= 0)cycle
            !- this is (DCAPE_env/Dt)_adv+rad
            !aa1_bl(i) = -aa3(i)
            !- Zhang threshold:  65 J/kg/hour => 65/(Rd *3600)= 63 10^-6 K/s
            aa1_bl(i) = aa3(i)-(63.e-6)!*1.5
            !print*,"dcape_env=",aa3(i),aa1_bl(i)
            if(xland(i) > 0.90 ) aa1_bl(i)=1.4*aa1_bl(i) !- over water
         enddo
         !--- this is (DT_ve/Dt)_cu
         do i=its,itf
            dtdt(i,:)=0.
            dqdt(i,:)=0.
            if(ierr(i) /= 0) cycle
            do k=max(kbcon(i),kts+1),ktop(i)
               dp=100.*(po_cup(i,k+1)-po_cup(i,k))
               RZenv=0.5*(zuo(i,k+1)+zuo(i,k) - (zdo(i,k+1)+zdo(i,k))*edto(i))
               S2= cp*tn_cup_x(i,k+1) + g*zo_cup(i,k+1)
               S1= cp*tn_cup_x(i,k  ) + g*zo_cup(i,k  )
               Q2= qo_cup_x(i,k+1)
               Q1= qo_cup_x(i,k  )

               dqdt(i,k)=         -RZenv*(Q2-Q1)*g/dp
               dtdt(i,k)= -(1./cp)*RZenv*(S2-S1)*g/dp

               dqdt(i,k)= dqdt(i,k)+(up_massdetro(i,k)*0.5*(qco (i,k+1)+qco (i,k)-(Q2+Q1)) &
                  + edto(i)*dd_massdetro(i,k)*0.5*(qcdo(i,k+1)+qcdo(i,k)-(Q2+Q1)))*g/dp

               dtdt(i,k)= dtdt(i,k)+(up_massdetro(i,k)*0.5*(tempco  (i,k+1) + tempco  (i,k)-  &
                  (tn_cup_x(i,k+1) + tn_cup_x(i,k))) &
                  + edto(i)*dd_massdetro(i,k)*0.5*(tempcdo (i,k+1) + tempcdo (i,k)-  &
                  (tn_cup_x(i,k+1) + tn_cup_x(i,k))) &
                  )*g/dp
               !print*,"dtdt=",k, dtdt(i,k),zuo(i,k+1),zdo(i,k+1),dqdt(i,k)
            enddo
            xk_x(i)=0.
            do k=max(kbcon(i),kts+1),ktop(i)
               dp=-(log(100.*po_cup(i,k+1))-log(100.*po_cup(i,k)))      ! no units here
               xk_x(i)=xk_x(i)+   ( (1.+0.608*qo_x(i,k))*dtdt(i,k) + &
                  0.608*tn_x(i,k) *dqdt(i,k) )*dp !  units=K m/Pa s2
               !=> aa3/xk_x will have units of kg/m2/s for the mass flux at cloud base.
               !print*,"xk_x=",k, xk_x(i),dtdt(i,k),dqdt(i,k)
            enddo
         enddo
      endif
      !
      !--- Trigger function based on Xie et al (2019)
      !
       if(ADV_TRIGGER == 3 .and. cumulus == 'deep') then
     ! if(ADV_TRIGGER == 3) then
         daa_adv_dt=0.
         do step=1,2
            !--- calculate moist static energy, heights, qes, ... only by ADV tendencies
            if(step==1) then
               tn_x = t
               qo_x = q
            else
               tn_x = tn_adv
               qo_x = qo_adv
            endif
            call cup_env(zo,qeso_x,heo_x,heso_x,tn_x,qo_x,po,z1,psur,ierr,-1,itf,ktf, its,ite, kts,kte)
            call cup_env_clev(tn_x,qeso_x,qo_x,heo_x,heso_x,zo,po,qeso_cup_x,qo_cup_x,heo_cup_x,us,vs   &
               ,u_cup_x,v_cup_x,heso_cup_x,zo_cup_x,po_cup_x,gammao_cup_x,tn_cup_x,psur,tsur  &
               ,ierr,z1,itf,ktf,its,ite, kts,kte)

            !--- get MSE
            do i=its,itf
               if(ierr(i) /= 0) cycle
               call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),heo_cup_x(i,kts:kte),hkbo_x(i),k22(i))
               hco_x (i,kts:start_level(i)) = hkbo_x(i)

               do k = start_level(i) +1,ktop(i)+1

                  denom=(zuo(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))
                  if(denom > 0.0)then
                     hco_x (i,k)=(hco_x(i,k-1)*zuo(i,k-1)-.5*up_massdetro(i,k-1)*hco_x(i,k-1) + &
                        up_massentro(i,k-1)*heo_x(i,k-1))/ denom
                  else
                     hco_x (i,k)=hco_x(i,k-1)
                  endif
               enddo
               hco_x (i,ktop(i)+2:ktf)=heso_cup_x(i,ktop(i)+2:ktf)
            enddo
            call get_buoyancy(itf,ktf, its,ite, kts,kte,ierr,klcl,kbcon,ktop &
               ,hco_x,heo_cup_x,heso_cup_x,dbyo_x,zo_cup_x)
            !--- get cloud work function
            aa_tmp=0.
            call cup_up_aa0(aa_tmp,zo_cup_x,zuo,dbyo_x,GAMMAo_CUP_x,tn_cup_x   &
               ,k22,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte)

            if(step==1) aa_ini=aa_tmp ! cloud work function initial
            if(step==2) aa_adv=aa_tmp ! cloud work function modified by advection tendencies
         enddo
         !
         do i=its,itf
            if(ierr(i) /= 0) cycle

            daa_adv_dt(i)=(aa_adv(i)-aa_ini(i))/dtime

            !print*,"daa_adv_dt J. kg-1 hr-1=",daa_adv_dt(i)*3600.
            ! call flush(6)

            if( daa_adv_dt(i) > dcape_threshold/3600. .and. aa_ini(i) > 0.) cycle !
            ierr(i)=90
            ierrc(i) = "dcape trigger not satisfied"

         enddo
         !--- only for output
         AA0_(:) = daa_adv_dt(:)*3600. ! J/kg/hour

      endif
      !
      !---
      !
      !--- DETERMINE DOWNDRAFT STRENGTH IN TERMS OF WINDSHEAR
      !
      call cup_dd_edt(cumulus,ierr,us,vs,zo,ktop,kbcon,edt,po,pwavo, &
         pwo,ccn,pwevo,edtmax,edtmin,maxens2,edtc,psum,psumh, &
         ccnclean,rho,aeroevap,itf,ktf,ipr,jpr,its,ite, kts,kte,vshear)

      do iedt=1,maxens2
         do i=its,itf
            if(ierr(i).eq.0)then
               edto(i)=sigd(i)*edtc(i,iedt)
               edt (i)=edto(i)
            endif
         enddo
         !
         !--- get the environmental mass flux
         !
         do i=its,itf
            zenv(i,:) = 0.0
            if(ierr(i) /= 0) cycle
            zenv(i,:) = zuo(i,:)-edto(i)*zdo(i,:)
         enddo
         !
         !--- check mass conservation
         !
         do i=its,itf
            if(ierr(i) /= 0) cycle
            do k=kts,ktop(i)
               ! these three are only used at or near mass detrainment and/or entrainment levels
               entupk=0.
               detupk=0.
               entdoj=0.
               ! detrainment and entrainment for downdrafts
               detdo=edto(i)*dd_massdetro(i,k)
               entdo=edto(i)*dd_massentro(i,k)
               ! entrainment/detrainment for updraft
               entup=up_massentro(i,k)
               detup=up_massdetro(i,k)
               ! subsidence by downdrafts only
               subin=-zdo(i,k+1)*edto(i)
               subdown=-zdo(i,k)*edto(i)
               if(k.eq.ktop(i))then
                  detupk=zuo(i,ktop(i))
                  subin=0.
                  subdown=0.
                  detdo=0.
                  entdo=0.
                  entup=0.
                  detup=0.
               endif
               totmas=subin-subdown+detup-entup-entdo+ &
                  detdo-entupk-entdoj+detupk+zuo(i,k+1)-zuo(i,k)
               if(abs(totmas).gt.1.e-6)then
                  write(6,*)'**mass cons: k,ktop,zo(ktop),totmas,subin,subdown,detup,entup,detdo,entdo,entupk,detupk'
                  write(6,123)'mass*1.e+6',k,ktop(i),zo(i,ktop(i)),totmas*1.e+6,subin*1.e+6,subdown*1.e+6,detup*1.e+6,entup*1.e+6&
                     ,detdo*1.e+6,entdo*1.e+6,entupk*1.e+6,detupk*1.e+6
123               format(1X,A11,2i5,10e12.5)
                 ! call error_fatal ( 'totmas .gt.1.e-6' )
               endif
            enddo   ! k
         enddo


         !
         !--- change per unit mass that a model cloud would modify the environment
         !
         !--- 1. in bottom layer
         !
         dellu     =0.
         dellv     =0.
         dellah    =0.
         dellat    =0.
         dellaq    =0.
         dellaqc   =0.
         dellabuoy =0.
         subten_H  =0.
         subten_Q  =0.
         subten_T  =0.

         !
         !----------------------------------------------  cloud level ktop
         !
         !- - - - - - - - - - - - - - - - - - - - - - - - model level ktop-1
         !      .               .                 .
         !      .               .                 .
         !      .               .                 .
         !      .               .                 .
         !      .               .                 .
         !      .               .                 .
         !
         !----------------------------------------------  cloud level k+2
         !
         !- - - - - - - - - - - - - - - - - - - - - - - - model level k+1
         !
         !----------------------------------------------  cloud level k+1
         !
         !- - - - - - - - - - - - - - - - - - - - - - - - model level k
         !
         !----------------------------------------------  cloud level k
         !
         !      .               .                 .
         !      .               .                 .
         !      .               .                 .
         !      .               .                 .
         !      .               .                 .
         !      .               .                 .
         !      .               .                 .
         !      .               .                 .
         !      .               .                 .
         !      .               .                 .
         !
         !----------------------------------------------  cloud level 3  _cup
         !
         !- - - - - - - - - - - - - - - - - - - - - - - - model level 2
         !
         !----------------------------------------------  cloud level 2  _cup
         !
         !- - - - - - - - - - - - - - - - - - - - - - - - model level 1

         !
         !------------------------------------------------------------------------------------------
         if(VERT_DISCR == 0) then

            do i=its,itf
               if(ierr(i) /= 0) cycle
               do k=kts,ktop(i)

                  dp=100.*(po_cup(i,k)-po_cup(i,k+1))

                  dellu(i,k) =-(zuo(i,k+1)*(uc (i,k+1)-u_cup(i,k+1) ) -             &
                                zuo(i,k  )*(uc (i,k  )-u_cup(i,k  ) ) )*g/dp        &
                              +(zdo(i,k+1)*(ucd(i,k+1)-u_cup(i,k+1)) -              &
                                zdo(i,k  )*(ucd(i,k  )-u_cup(i,k  )) )*g/dp*edto(i)

                   dellv(i,k) =-(zuo(i,k+1)*(vc (i,k+1)-v_cup(i,k+1) ) -             &
                                 zuo(i,k  )*(vc (i,k  )-v_cup(i,k  ) ) )*g/dp        &
                               +(zdo(i,k+1)*(vcd(i,k+1)-v_cup(i,k+1) ) -             &
                                 zdo(i,k  )*(vcd(i,k  )-v_cup(i,k  ) ) )*g/dp*edto(i)

               enddo   ! k
            enddo

            do i=its,itf
               trash  = 0.0
               trash2 = 0.0
               if(ierr(i).eq.0)then
                  do k=kts,ktop(i)

                     dp=100.*(po_cup(i,k)-po_cup(i,k+1))

                     dellah(i,k) =-(zuo(i,k+1)*(hco (i,k+1)-heo_cup(i,k+1) ) -                 &
                                    zuo(i,k  )*(hco (i,k  )-heo_cup(i,k  ) ) )*g/dp            &
                                  +(zdo(i,k+1)*(hcdo(i,k+1)-heo_cup(i,k+1) ) -                 &
                                    zdo(i,k  )*(hcdo(i,k  )-heo_cup(i,k  ) ) )*g/dp*edto(i)

                     !---meltglac-------------------------------------------------
                     dellah(i,k) = dellah(i,k) + xlf*((1.-p_liq_ice(i,k))*0.5*(qrco(i,k+1)+qrco(i,k)) &
                                 - melting(i,k))*g/dp

                     !-- for output only
                     subten_H(i,k) = -(zuo(i,k+1)*(-heo_cup(i,k+1)) - zuo(i,k)*(-heo_cup(i,k)))*g/dp       &
                                     +(zdo(i,k+1)*(-heo_cup(i,k+1)) - zdo(i,k)*(-heo_cup(i,k)))*g/dp*edto(i)

                      !- check H conservation
                     trash2 = trash2+ (dellah(i,k))*dp/g

                      !-- take out cloud liquid/ice water for detrainment
                     detup=up_massdetro(i,k)
                     if( cumulus == 'mid'  .or. cumulus == 'shallow') then

                        dellaqc(i,k) = detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp

                     elseif(cumulus == 'deep') then

                        if(.not. USE_C1D) then

                           dellaqc(i,k) = detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp

                        elseif(c1>0.0) then
                           if(k == ktop(i)) then
                              dellaqc(i,k) = detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp
                           else
                              dz=zo_cup(i,k+1)-zo_cup(i,k)
                              dellaqc(i,k) = zuo(i,k)*c1d(i,k)*qrco(i,k)*dz/dp*g
                           endif
                        else
                           if(k == ktop(i)) then
                              dellaqc(i,k) = detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp
                           else
                              dz=zo_cup(i,k+1)-zo_cup(i,k)
                              dellaqc(i,k) = ( zuo(i,k)*c1d(i,k)*qrco(i,k)*dz/dp*g  + &
                                             detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp )*0.5
                           endif
                        endif
                     endif
                     !
                     !---
                     G_rain=  0.5*(pwo (i,k)+pwo (i,k+1))*g/dp
                     E_dn  = -0.5*(pwdo(i,k)+pwdo(i,k+1))*g/dp*edto(i) ! pwdo < 0 and E_dn must > 0
                     !
                     !print*,"eva=",k,pwdo(i,k),E_dn,zdo(i,k  ),G_rain
                     !
                     !-- condensation source term = detrained + flux divergence of
                     !-- cloud liquid/ice water (qrco) + converted to rain

                     C_up = dellaqc(i,k)+(zuo(i,k+1)* qrco(i,k+1) -       &
                        zuo(i,k  )* qrco(i,k  )  )*g/dp + G_rain

                     !-- water vapor budget
                     !-- = flux divergence z*(Q_c - Q_env)_up_and_down &
                     !--   - condensation term + evaporation
                     dellaq(i,k) =-(zuo(i,k+1)*(qco (i,k+1)-qo_cup(i,k+1) ) -                 &
                                    zuo(i,k  )*(qco (i,k  )-qo_cup(i,k  ) ) )*g/dp            &
                                  +(zdo(i,k+1)*(qcdo(i,k+1)-qo_cup(i,k+1) ) -                 &
                                    zdo(i,k  )*(qcdo(i,k  )-qo_cup(i,k  ) ) )*g/dp*edto(i)    &
                                  - C_up + E_dn

                     !-- for output only
                     subten_Q(i,k) =-(zuo(i,k+1)*(-qo_cup(i,k+1)) - zuo(i,k)*(-qo_cup(i,k)))*g/dp       &
                                    +(zdo(i,k+1)*(-qo_cup(i,k+1)) - zdo(i,k)*(-qo_cup(i,k)))*g/dp*edto(i)

                     !- check water conservation liq+condensed (including rainfall)
                     trash= trash+ (dellaq(i,k)+dellaqc(i,k)+ G_rain-E_dn)*dp/g

                     !---
                     dellabuoy(i,k) = edto(i)*dd_massdetro(i,k)*0.5*(dbydo(i,k+1)+dbydo(i,k))*g/dp
                     !---

                     !write(3,*)'=>H= ',k,real(trash2,4),real(dellah(i,k),4)
                     !write(4,*)'=>W= ',k,real(trash,4),real(dellaq(i,k),4)
                  enddo   ! k
                !--- test only with double precision:
                !write(0,*)'=>H/W-FINAL= ',real(trash2,4),real(trash,4),k22(i),kbcon(i),ktop(i)
                !if(abs(trash)>1.e-6 .or. abs(trash2) > 1.e-6) then
                !    write(0,*)'=> not water mass or H cons for deep= ',i,trash,trash2
                !    !stop 33
                !endif

               endif

            enddo

         elseif(VERT_DISCR == 1) then

            !---- convective transport of momentum
            if(alp1 == 0.) then !-- fully time explicit
               do i=its,itf
                  if(ierr(i) /= 0) cycle
                  do k=kts,ktop(i)
                     dp=100.*(po_cup(i,k)-po_cup(i,k+1))

                     dellu(i,k) =-(zuo(i,k+1)*(uc (i,k+1)-u_cup(i,k+1) ) -             &
                                   zuo(i,k  )*(uc (i,k  )-u_cup(i,k  ) ) )*g/dp        &
                                 +(zdo(i,k+1)*(ucd(i,k+1)-u_cup(i,k+1)) -              &
                                   zdo(i,k  )*(ucd(i,k  )-u_cup(i,k  )) )*g/dp*edto(i)

                     dellv(i,k) =-(zuo(i,k+1)*(vc (i,k+1)-v_cup(i,k+1) ) -             &
                                   zuo(i,k  )*(vc (i,k  )-v_cup(i,k  ) ) )*g/dp        &
                                 +(zdo(i,k+1)*(vcd(i,k+1)-v_cup(i,k+1) ) -             &
                                   zdo(i,k  )*(vcd(i,k  )-v_cup(i,k  ) ) )*g/dp*edto(i)
                  enddo   ! k
               enddo

            elseif (alp1 > 0.) then              !-- time alp0*explict + alp1*implicit + upstream

               alp0=1.-alp1
               do i=its,itf
                  if(ierr(i) /= 0) cycle
                  do k=kts,ktop(i)+1
                     fp(k) = 0.5*(zenv(i,k)+abs(zenv(i,k)))
                     fm(k) = 0.5*(zenv(i,k)-abs(zenv(i,k)))
                  enddo

                  do k=kts,ktop(i)
                     dp=100.*(po_cup(i,k)-po_cup(i,k+1))

                     beta1 = dtime*g/dp
                     aa(k) =    alp1*beta1*fm(k)
                     bb(k) = 1.+alp1*beta1*(fp(k)-fm(k+1))
                     cc(k) =   -alp1*beta1*fp(k+1)

                     ddu(k) = us(i,k)-( zuo(i,k+1)*uc (i,k+1)-zuo(i,k)*uc (i,k) )*beta1 + &
                                      ( zdo(i,k+1)*ucd(i,k+1)-zdo(i,k)*ucd(i,k) )*beta1*edto(i)

                     ddu(k) = ddu(k) + alp0*beta1*(-fm(k)*us(i,max(kts,k-1)) +(fm(k+1)-fp(k))*us(i,k) +fp(k+1)*us(i,k+1))


                     ddv(k) = vs(i,k)-( zuo(i,k+1)*vc (i,k+1)-zuo(i,k)*vc (i,k) )*beta1 + &
                                      ( zdo(i,k+1)*vcd(i,k+1)-zdo(i,k)*vcd(i,k) )*beta1*edto(i)

                     ddv(k) = ddv(k) + alp0*beta1*(-fm(k)*vs(i,max(kts,k-1)) +(fm(k+1)-fp(k))*vs(i,k) +fp(k+1)*vs(i,k+1))

                    !print*,"trX4=",k,aa(k),bb(k),cc(k)!, 1.+alp1*beta1*zenv(i,k  ), -alp1*beta1*zenv(i,k+1)
                  enddo
                  call tridiag (ktop(i),aa (kts:ktop(i)),bb (kts:ktop(i)),cc (kts:ktop(i)),ddu (kts:ktop(i)))
                  dellu(i,kts:ktop(i))=(ddu(kts:ktop(i))-us(i,kts:ktop(i)))/dtime

                  call tridiag (ktop(i),aa (kts:ktop(i)),bb (kts:ktop(i)),cc (kts:ktop(i)),ddv (kts:ktop(i)))
                  dellv(i,kts:ktop(i))=(ddv(kts:ktop(i))-vs(i,kts:ktop(i)))/dtime
               enddo
            endif


            !--- convective transport of MSE and Q/Qc
            !if(USE_FLUX_FORM == 1) then
            do i=its,itf
               if(ierr(i) /= 0) cycle

               !--- moist static energy : flux form + source/sink terms + time explicit
               !
               !   if(use_fct == 0 .or. adjustl(cumulus) == 'shallow') then
               if(use_fct == 0 ) then

                  do k=kts,ktop(i)
                     dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                     dellah(i,k) =-(zuo(i,k+1)*(hco (i,k+1)-heo_cup(i,k+1) ) -          &
                                    zuo(i,k  )*(hco (i,k  )-heo_cup(i,k  ) ) )*g/dp     &
                                   +(zdo(i,k+1)*(hcdo(i,k+1)-heo_cup(i,k+1) ) -          &
                                     zdo(i,k  )*(hcdo(i,k  )-heo_cup(i,k  ) ) )*g/dp*edto(i)

                     dellah(i,k) = dellah(i,k) + xlf*((1.-p_liq_ice(i,k))* &
                                 0.5*(qrco(i,k+1)+qrco(i,k)) - melting(i,k))*g/dp

                     !--- for output only
                     subten_H(i,k) = -(zuo(i,k+1)*(-heo_cup(i,k+1)) - zuo(i,k)*(-heo_cup(i,k)))*g/dp       &
                                     +(zdo(i,k+1)*(-heo_cup(i,k+1)) - zdo(i,k)*(-heo_cup(i,k)))*g/dp*edto(i)
                  enddo   ! k

               else

                  !-- FCT scheme for the subsidence transport: d(M_env*S_env)/dz
                  sub_tend (1,:) = 0. ! dummy array
                  trcflx_in(1,:) = 0. ! dummy array
                  massflx  (i,:) = 0.
                  dtime_max      = dtime

                  do k=kts,ktop(i)
                     dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                     trcflx_in (1,k) =-(zuo(i,k)  -edto(i)*zdo(i,k))*heo_cup(i,k) !* xmb(i)
                     massflx   (i,k) =-(zuo(i,k)  -edto(i)*zdo(i,k))        !* xmb(i)
                     dtime_max=min(dtime_max,.5*dp)
                  enddo
                  call fct1d3 (ktop(i),kte,dtime_max,po_cup(i,:),heo(i,:),massflx(i,:),trcflx_in(1,:),sub_tend(1,:))

                  do k=kts,ktop(i)
                     dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                     dellah(i,k) =-( zuo(i,k+1)*hco (i,k+1) - zuo(i,k)*hco (i,k) )*g/dp      &
                                  +( zdo(i,k+1)*hcdo(i,k+1) - zdo(i,k)*hcdo(i,k) )*g/dp*edto(i)

                     dellah(i,k) = dellah(i,k) + xlf*((1.-p_liq_ice(i,k))* &
                                 0.5*(qrco(i,k+1)+qrco(i,k)) - melting(i,k))*g/dp
                     !- update with subsidence term from the FCT scheme
                     dellah(i,k) = dellah(i,k) + sub_tend(1,k)
                     !--- for output only
                     subten_H(i,k) = sub_tend(1,k)
                  enddo   ! k
               endif
            enddo

            !     elseif(USE_FLUX_FORM == 2) THEN
            !
            !        !- flux form + source/sink terms + time explicit + upstream with anti-diffusion step (Smolarkiewicz 1983)
            !        alp0=1.
            !        do i=its,itf
            !          if(ierr(i) /= 0) cycle
            !    do istep=1,-1, -2
            !
            !      if(istep == 1) then
            !         ddu(:) = heo(i,:)
            !     do k=kts,ktop(i)+1
            !       fp(k) = 0.5*(zenv(i,k)+abs(zenv(i,k)))
            !       fm(k) = 0.5*(zenv(i,k)-abs(zenv(i,k)))
            !     enddo
            !       else
            !         ddu(kts:ktop(i)+1) = heo(i,kts:ktop(i)+1) + dellah(i,kts:ktop(i)+1)*dtime
            !         zenv_diff(1,kts) = 0.
            !         do k=kts,ktop(i)+1
            !        dp = 100.*(po_cup(i,k)-po_cup(i,k+1))
            !        zenv_diff (1,k+1) = 1.06* ( dp*abs(zenv(i,k+1))/g - dtime*zenv(i,k+1)**2 )/dp/g &
            !                * (ddu(k+1) - ddu(k)) /(ddu(k+1) + ddu(k) + 1.e-16)
            !         enddo
            !         do k=kts,ktop(i)+1
            !        fp(k) = 0.5*(zenv_diff(1,k)+abs(zenv_diff(1,k)))
            !        fm(k) = 0.5*(zenv_diff(1,k)-abs(zenv_diff(1,k)))
            !         enddo
            !       endif
            !       do k=kts,ktop(i)
            !         dp=100.*(po_cup(i,k)-po_cup(i,k+1))
            !         beta1 = dtime*g/dp
            !         ddh(k) = ddu(k) + alp0*beta1*( -fm(k)*ddu(max(kts,k-1)) + (fm(k+1)-fp(k))*ddu(k) + fp(k+1)*ddu(k+1) )
            !       enddo
            !
            !       dellah(i,kts:ktop(i)+1)=(ddh(kts:ktop(i)+1)-heo(i,kts:ktop(i)+1))/dtime
            !
            !     enddo
            !
            !     do k=kts,ktop(i)
            !       dp=100.*(po_cup(i,k)-po_cup(i,k+1))
            !       beta1 = g/dp
            !
            !       ddh(k) =  -( zuo(i,k+1)*hco (i,k+1) - zuo(i,k)*hco (i,k) )*beta1       &
            !            +( zdo(i,k+1)*hcdo(i,k+1) - zdo(i,k)*hcdo(i,k) )*beta1*edto(i)
            !
            !       ddh(k) = ddh(k) + xlf*((1.-p_liq_ice(i,k))* &
            !              0.5*(qrco(i,k+1)+qrco(i,k)) - melting(i,k))*beta1
            !
            !       dellah(i,k) =  dellah(i,k) + ddh(k)
            !
            !     enddo
            !       enddo
            !     endif
            !-------------------------
            !--- water vapor + condensates : flux form + source/sink terms + time explicit
            do i=its,itf

               if(ierr(i) /= 0) cycle
               !
               do k=kts,ktop(i)
                  dp=100.*(po_cup(i,k)-po_cup(i,k+1))

                   !-- take out cloud liquid/ice water for detrainment
                  detup=up_massdetro(i,k)
                  if( cumulus == 'mid'  .or. cumulus == 'shallow') then

                     dellaqc(i,k) = detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp

                  elseif(cumulus == 'deep') then

                     if(.not. USE_C1D) then

                        dellaqc(i,k) = detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp

                     elseif(c1>0.0) then
                        if(k == ktop(i)) then
                           dellaqc(i,k) = detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp
                        else
                           dz=zo_cup(i,k+1)-zo_cup(i,k)
                           dellaqc(i,k) = zuo(i,k)*c1d(i,k)*qrco(i,k)*dz/dp*g
                        endif
                     else
                        if(k == ktop(i)) then
                           dellaqc(i,k) = detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp
                        else
                           dz=zo_cup(i,k+1)-zo_cup(i,k)
                           dellaqc(i,k) = ( zuo(i,k)*c1d(i,k)*qrco(i,k)*dz/dp*g  + &
                              detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp )*0.5
                        endif
                     endif
                  endif
                   !
                  !---
                  G_rain=  0.5*(pwo (i,k)+pwo (i,k+1))*g/dp
                  E_dn  = -0.5*(pwdo(i,k)+pwdo(i,k+1))*g/dp*edto(i) ! pwdo < 0 and E_dn must > 0
                  !
                  !-- condensation source term = detrained + flux divergence of
                  !-- cloud liquid/ice water (qrco) + converted to rain

                  C_up = dellaqc(i,k)+(zuo(i,k+1)*qrco(i,k+1) - zuo(i,k)* qrco(i,k))*g/dp + G_rain

                  !-- water vapor budget
                  !-- = flux divergence z*(Q_c - Q_env)_up_and_down  - condensation term + evaporation
                  dellaq(i,k) =-(zuo(i,k+1)*qco (i,k+1) - zuo(i,k)*qco (i,k))*g/dp     &
                               +(zdo(i,k+1)*qcdo(i,k+1) - zdo(i,k)*qcdo(i,k))*g/dp*edto(i)    &
                               - C_up + E_dn

                   !--- source of cold pools
                  dellabuoy(i,k)=edto(i)*dd_massdetro(i,k)*0.5*(dbydo(i,k+1)+dbydo(i,k))*g/dp

               enddo
               !        if(use_fct == 0 .or. adjustl(cumulus) == 'shallow') then
               if(use_fct == 0 ) then
                  do k=kts,ktop(i)
                     dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                     sub_tend(1,k) =-(zuo(i,k+1)*(-qo_cup(i,k+1)) - zuo(i,k)*(-qo_cup(i,k)))*g/dp       &
                                    +(zdo(i,k+1)*(-qo_cup(i,k+1)) - zdo(i,k)*(-qo_cup(i,k)))*g/dp*edto(i)
                  enddo
               else
                  !-- FCT scheme for the subsidence transport: d(M_env*S_env)/dz
                  sub_tend (1,:) = 0. ! dummy array
                  trcflx_in(1,:) = 0. ! dummy array
                  massflx  (i,:) = 0.
                  dtime_max      = dtime

                  do k=kts,ktop(i)
                     dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                     trcflx_in (1,k) =-(zuo(i,k)  -edto(i)*zdo(i,k))*qo_cup(i,k) !* xmb(i)
                     massflx   (i,k) =-(zuo(i,k)  -edto(i)*zdo(i,k))             !* xmb(i)
                     dtime_max=min(dtime_max,.5*dp)
                  enddo
                  call fct1d3 (ktop(i),kte,dtime_max,po_cup(i,:),qo(i,:),massflx(i,:),trcflx_in(1,:),sub_tend(1,:))
               endif

               !--- add the contribuition from the environ subsidence
               dellaq(i,kts:ktop(i)) = dellaq(i,kts:ktop(i)) + sub_tend(1,kts:ktop(i))

                !--- for output only
               subten_Q (i,kts:ktop(i)) = sub_tend(1,kts:ktop(i))

               !     do k=kts,ktop(i)
               !     print*,"delq=",use_fct,k,dellaq(i,k) , sub_tend(1,k)
               !     enddo

               !- check H and water conservation liq+condensed (including rainfall)
               trash  = 0.
               trash2 = 0.0
               do k=kts,ktop(i)
                  dp     = 100.*(po_cup(i,k)-po_cup(i,k+1))
                  G_rain =  0.5*(pwo (i,k)+pwo (i,k+1))*g/dp
                  E_dn   = -0.5*(pwdo(i,k)+pwdo(i,k+1))*g/dp*edto(i)
                  trash  = trash + (dellaq(i,k) + dellaqc(i,k)+ G_rain-E_dn)*dp/g
                  trash2 = trash2+  dellah(i,k)*g/dp + xlf*((1.-p_liq_ice(i,k))*0.5*(qrco(i,k+1)+qrco(i,k)) &
                     - melting(i,k))*g/dp
               enddo   ! k
               !--- test only with double precision:
               !write(0,*)'=>H/W-FINAL= ',real(trash2,4),real(trash,4),k22(i),kbcon(i),ktop(i)
               !if(abs(trash)>1.e-6 .or. abs(trash2) > 1.e-6) then
               !    write(0,*)'=> not water mass or H cons for deep= ',i,trash,trash2
               !    !stop 33
               !endif

            enddo


         endif ! vertical discretization formulation

         !
         !--- apply environmental subsidence on grid-scale ice and liq water contents, and cloud fraction (Upwind scheme)
         !
         if(APPLY_SUB_MP == 1) then
            dellampqi =0.
            dellampql =0.
            dellampcf =0.

            do i=its,itf
               if(ierr(i) /= 0) cycle
               do k=kts,ktop(i)
                  dp=100.*(po_cup(i,k  )-po_cup(i,k+1))

                  !--- apply environmental subsidence on grid-scale/anvil ice and liq water contents (Upwind scheme)
                  !
                  env_mf   = - 0.5* (zenv(i,k+1) + zenv(i,k))
                  env_mf_m = min(env_mf,0.)*g/dp
                  env_mf_p = max(env_mf,0.)*g/dp

                  dellampqi(:,i,k) = - (  env_mf_m*(mpqi(:,i,k+1)-mpqi(:,i,k))  +        &
                                          env_mf_p*(mpqi(:,i,k  )-mpqi(:,i,max(k-1,kts))))
                  dellampql(:,i,k) = - (  env_mf_m*(mpql(:,i,k+1)-mpql(:,i,k))  +         &
                                          env_mf_p*(mpql(:,i,k  )-mpql(:,i,max(k-1,kts))))

                  !--- apply environmental subsidence on grid-scale/anvil cloud fraction
                  dellampcf(:,i,k) = - (  env_mf_m*(mpcf(:,i,k+1)-mpcf(:,i,k))  +         &
                                          env_mf_p*(mpcf(:,i,k  )-mpcf(:,i,max(k-1,kts))))
               enddo

               !--- apply environmental subsidence on grid-scale and anvil cloud fraction using time implicit/explict method
               if(alp1 > 0.) then
                  alp0=1.0-alp1
                  do k=kts,ktop(i)
                     dp=100.*(po_cup(i,k  )-po_cup(i,k+1))
                     env_mf   = - 0.5* (zenv(i,k+1) + zenv(i,k))
                     env_mf_m = min(env_mf,0.)*g/dp
                     env_mf_p = max(env_mf,0.)*g/dp

                     beta1 = -env_mf_m
                     beta2 = -env_mf_p

                     aa(k) =    alp1*beta2             ! coef of f(k-1,t+1),
                     bb(k) = 1.+alp1*beta1-alp1*beta2  ! coef of f(k  ,t+1),
                     cc(k) =   -alp1*beta1             ! coef of f(k+1,t+1),

                     !-- this is the rhs of the discretization
                     dd(:,k) = (1.-alp0*beta1+alp0*beta2)*mpcf(:,i,k  ) +&       ! coef of  f(k  ,t),
                        alp0*beta1 *mpcf(:,i,k+1) -&       ! coef of  f(k+1,t),
                        alp0*beta2 *mpcf(:,i,max(kts,k-1)) ! coef of  f(k-1,t),
                  enddo
                  do kmp =1,nmp
                      !-- this routine solves the problem: aa*f(k-1,t+1) + bb*f(k,t+1) + cc*f(k+1,t+1) = dd
                     call tridiag (ktop(i),aa(kts:ktop(i)), bb(kts:ktop(i)), cc(kts:ktop(i)), dd(kmp,kts:ktop(i)))

                     dellampcf(kmp,i,kts:ktop(i)) = dd(kmp,kts:ktop(i))-mpcf(kmp,i,kts:ktop(i))
                  enddo
               endif
            enddo
         endif
         !
         !
         !--- make the smoothness procedure
         !
         if(USE_SMOOTH_TEND > 0) then
            do i=its,itf
               if(ierr(i) /= 0) cycle
               tend2d=0.

               do k=kts,ktop(i)
                  rcount = 1.e-8
                  tend1d=0.
                  do kk= max(kts,k-USE_SMOOTH_TEND),min(ktop(i),k+USE_SMOOTH_TEND)
                     dp         = (po_cup(i,kk)-po_cup(i,kk+1))
                     rcount     = rcount     +  dp
                     tend1d(1)  = tend1d(1)  +  dp* DELLAH  (i,kk)
                     tend1d(2)  = tend1d(2)  +  dp* DELLAQ  (i,kk)
                     tend1d(3)  = tend1d(3)  +  dp* DELLAQC (i,kk)
                     tend1d(4)  = tend1d(4)  +  dp* DELLU   (i,kk)
                     tend1d(5)  = tend1d(5)  +  dp* DELLV   (i,kk)
                  enddo
                  tend2d(k,1:5)  = tend1d(1:5) /rcount
               enddo
               !--- get the final/smoother tendencies
               do k=kts,ktop(i)
                  DELLAH  (i,k) = tend2d(k,1)
                  DELLAQ  (i,k) = tend2d(k,2)
                  DELLAQC (i,k) = tend2d(k,3)
                  DELLU   (i,k) = tend2d(k,4)
                  DELLV   (i,k) = tend2d(k,5)
               enddo
            enddo
         endif ! USE_SMOOTH_TEND == 1
         !
         !--- using dellas, calculate changed environmental profiles
         !
         do k=kts,ktf
            do i=its,itf
               dellat(i,k)=0.
               if(ierr(i) /= 0) cycle
               !
               XHE(I,K)=(DELLAH(I,K)             )*MBDT(i)+HEO(I,K)
               XQ (I,K)=(DELLAQ(I,K)+DELLAQC(i,k))*MBDT(i)+QO(I,K)
               if(XQ(I,K).le.0.)XQ(I,K)=1.e-08

               !- do not feed dellat with dellaqc if the detrainment of liquid water
                !- will be used as a source for cloud microphysics
               if(COUPL_MPHYSICS) then
                  DELLAT(I,K)=(1./cp)*(DELLAH(I,K)-xlv*DELLAQ(I,K))
               else
                  !---meltglac-------------------------------------------------
                  DELLAT (I,K)=(1./cp)*( DELLAH(I,K) - xlv*(DELLAQ(I,K) + DELLAQC(i,k))*(1.+(xlf/xlv)*(1.-p_liq_ice(i,k))))
                  !DELLAT (I,K)=(1./cp)*( DELLAH(I,K)  -xlv*(DELLAQ(I,K) + DELLAQC(i,k)))

                  !-adding dellaqc to dellaq:
                  DELLAQ (I,K)= DELLAQ(I,K)+DELLAQC(I,K)
                  DELLAQC(I,K)= 0.0
               endif
               !---meltglac-------------------------------------------------
               XT(I,K)=((1./cp)*DELLAH(i,k)-(xlv/cp)*(DELLAQ(i,k)+DELLAQC(i,k)*(1.+(xlf/xlv)*(1.-p_liq_ice(i,k)))))*MBDT(i) &
                  + TN(I,K)
               !XT(I,K)=((1./cp)*DELLAH(i,k)-(xlv/cp)*(DELLAQ(i,k)+DELLAQC(i,k)))*MBDT(i)+TN(I,K)

               !--- temp tendency due to the environmental subsidence
               subten_T(i,k)=(1./cp)*(subten_H(i,k)-xlv*subten_Q(i,k))

            enddo
         enddo
         do i=its,itf
            if(ierr(i) /= 0) cycle
            !XHKB(I)=(dsubh(i,k22(i))+DELLAH(I,K22(i)))*MBDT+HKBO(I)
            XHE(I,ktf)=HEO(I,ktf)
            XQ (I,ktf)=QO(I,ktf)
            XT (I,ktf)=TN(I,ktf)
            if(XQ(I,ktf).le.0.)XQ(I,ktf)=1.e-08
         enddo
         !- new way for defining XHKB
         do i=its,itf
            if(ierr(i) /= 0)cycle
            !XHKB(I)= DELLAH(I,K22(i))*MBDT+HKBO(I)
            !-note that HKBO already contains the contribuition from
            !-ztexec and zqexec
            call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),DELLAH (i,kts:kte),DELLAH_aver,k22(i))
            XHKB(I)= DELLAH_aver*MBDT(i) + HKBO(I)
         enddo
         !
         !--- calculate moist static energy, heights, qes
         !
         call cup_env(xz,xqes,xhe,xhes,xt,xq,po,z1,psur,ierr,-1,itf,ktf,its,ite,kts,kte)
         !
         !--- environmental values on cloud levels
         !
         call cup_env_clev(xt,xqes,xq,xhe,xhes,xz,po,xqes_cup,xq_cup, xhe_cup   &
            ,us,vs,u_cup,v_cup,xhes_cup,xz_cup,po_cup,gamma_cup    &
            ,xt_cup,psur,tsur,ierr,z1,itf,ktf,its,ite, kts,kte)
         !
         !--- static control
         !
         !--- moist static energy inside cloud
         !
         do i=its,itf
            xhc(i,:)=0.
            if(ierr(i)/= 0)cycle
            do k=kts,start_level(i) !k22(i)
               xhc(i,k)=xhkb(i)
            enddo
         enddo
         !
         !--- option to produce linear fluxes in the sub-cloud layer.
         if(cumulus == 'shallow' .and. use_linear_subcl_mf == 1) then
           do i=its,itf
             if(ierr(i) /= 0) cycle
             call get_delmix(cumulus,kts,kte,ktf,xland(i),start_level(i),po(i,kts:kte) &
                            ,xhe_cup (i,kts:kte), xhc (i,kts:kte))
           enddo
         endif
         do i=its,itf
            if(ierr(i)/=0)cycle
            do k=start_level(i)  +1,ktop(i)+1  ! mass cons option
               denom= (xzu(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))
               if(denom==0.0) then
                  xhc(i,k)= xhc(i,k-1)
               else
                  xhc(i,k)=(xhc(i,k-1)*xzu(i,k-1)-.5*up_massdetro(i,k-1)*xhc(i,k-1)+ &
                     up_massentro(i,k-1)*xhe(i,k-1)) / denom
                  if(k==start_level(i)+1)  then
                     x_add = (xlv*zqexec(i)+cp*ztexec(i)) +  x_add_buoy(i)
                     xhc(i,k)= xhc(i,k) + x_add*up_massentro(i,k-1)/denom
                  endif
               endif
               !
               !- include glaciation effects on XHC
               !                                   ------ ice content --------
               xhc (i,k)= xhc (i,k)+ xlf*(1.-p_liq_ice(i,k))*qrco(i,k)
            enddo
            do k=ktop(i)+2,ktf
               xHC (i,k)=xhes_cup(i,k)
               xzu (i,k)=0.
            enddo
         enddo
         call get_buoyancy(itf,ktf, its,ite, kts,kte,ierr,klcl,kbcon,ktop,xhc,xhe_cup,xhes_cup,xdby,xz_cup)
         !
         !--- workfunctions for updraft
         !
         call cup_up_aa0(xaa0,xz_cup,xzu,xdby,GAMMA_CUP,xt_cup, k22,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte)

         do nens=1,maxens
            do i=its,itf
               if(ierr(i) /= 0) cycle
               !~ xaa0_ens(i,nens)=xaa0(i)
               do k=kts,ktop(i)
                  do nens3=1,maxens3
                     if(nens3.eq.7)then
                        !--- b=0
                        pr_ens(i,nens3)=pr_ens(i,nens3) +pwo(i,k)+edto(i)*pwdo(i,k)
                     !--- b=beta
                     else if(nens3.eq.8)then
                        pr_ens(i,nens3)=pr_ens(i,nens3) +pwo(i,k)+edto(i)*pwdo(i,k)
                     !--- b=beta/2
                     else if(nens3.eq.9)then
                        pr_ens(i,nens3)=pr_ens(i,nens3) +pwo(i,k)+edto(i)*pwdo(i,k)
                     else
                        pr_ens(i,nens3)=pr_ens(i,nens3) +pwo(i,k)+edto(i)*pwdo(i,k)
                     endif
                  enddo
               enddo
               if(pr_ens(i,7).lt.1.e-6 .and. c0_mid > 0. .and.  cumulus /= 'shallow' )then
                  ierr(i)=18
                  ierrc(i)="total normalized condensate too small"
                  do nens3=1,maxens3
                     pr_ens(i,nens3)=0.
                  enddo
               endif
               do nens3=1,maxens3
                  if(pr_ens(i,nens3) < 1.e-5) pr_ens(i,nens3)=0.
               enddo
            enddo
         enddo
         !
         !--- LARGE SCALE FORCING
         !
         do i=its,itf
            ierr2(i)=ierr(i)
            ierr3(i)=ierr(i)
         enddo
         !
         !--- calculate cloud base mass flux
         !
         if(cumulus == 'deep') &
            call cup_forcing_ens_3d(itf,ktf,its,ite, kts,kte,ens4,ensdim,ichoice,maxens,maxens2,maxens3 &
            ,ierr,ierr2,ierr3,k22,kbcon,ktop         &
            ,xland1,aa0,aa1,xaa0,mbdt,dtime          &
            ,xf_ens,mconv,qo                         &
            ,po_cup,omeg,zdo,zuo,pr_ens,edto         &
            ,tau_ecmwf,aa1_bl,xf_dicycle, xk_x       &
            ,alpha_adv,Q_adv,aa1_radpbl,aa1_adv,wlpool_bcon,xf_coldpool)

         if(cumulus == 'mid') &
            call cup_forcing_ens_3d_mid(aa0,aa1,xaa0,mbdt,dtime,ierr                         &
            ,po_cup,ktop,k22,kbcon,kpbl,ichoice,maxens,maxens3    &
            ,itf,ktf,its,ite, kts,kte,tau_ecmwf,aa1_bl,xf_dicycle &
            ,dhdt,xff_mid,zws,hc,hco,he_cup,heo_cup,wlpool_bcon,xf_coldpool)

         if(cumulus == 'shallow') then
            call cup_up_cape(cape,z,zu,dby,gamma_cup,t_cup,k22,kbcon,ktop,ierr    &
               ,tempco,qco,qrco,qo_cup,itf,ktf,its,ite, kts,kte)

            call cup_forcing_ens_3d_shal(itf,ktf,its,ite,kts,kte,dtime,ichoice    &
               ,ierrc,ierr,klcl,kpbl,kbcon,k22,ktop      &
               ,xmb,tsur,cape,h_sfc_flux,le_sfc_flux,zws &
               ,po, hco, heo_cup,po_cup,t_cup,dhdt,rho   &
               ,xff_shal,xf_dicycle,tke_pbl,wlpool_bcon,xf_coldpool)
         endif


         !do k=kts,ktf
         !  do i=its,itf
         !    if(ierr(i) /= 0) cycle
         !    pwo_eff(i,k)=pwo(i,k)+edto(i)*pwdo(i,k)
         !  enddo
         !enddo
         !
         !--- get the net precipitation at surface
         !
         do i=its,itf
            if(ierr(i) == 0) then
               pwo_eff(i,:)=pwo(i,:)+edto(i)*pwdo(i,:)
            else
               pwo_eff(i,:)=0.
            endif
         enddo


      enddo
      !
      !--- Include kinetic energy dissipation converted to heating
      !
      call ke_to_heating(itf,ktf,its,ite, kts,kte,ktop,ierr,po_cup,us,vs,dellu,dellv,dellat)
      !
      !--- FEEDBACK
      !
      call cup_output_ens_3d(cumulus,xff_shal,xff_mid,xf_ens,ierr,dellat,dellaq,            &
         dellaqc,outt, outq,outqc,zuo,pre,pwo_eff,xmb,ktop,           &
         maxens2,maxens,ierr2,ierr3,                                  &
         pr_ens,maxens3,ensdim,sig,xland1,                            &
         ichoice,ipr,jpr,itf,ktf,its,ite, kts,kte,                    &
         xf_dicycle,outu,outv,dellu,dellv,dtime,po_cup,kbcon,         &
         dellabuoy,outbuoy,                                           &
         dellampqi,outmpqi,dellampql,outmpql,dellampcf,outmpcf,nmp,   &
         rh_dicycle_fct,xf_coldpool,wlpool_bcon)

      !
      !
      !--- get the net precipitation flux (after downdraft evaporation)
      call get_precip_fluxes(cumulus,klcl,kbcon,ktop,k22,ierr,xland,pre,xmb  &
         ,pwo,pwavo,edto,pwevo,pwdo,t_cup,tempco          &
         ,prec_flx,evap_flx,itf,ktf,its,ite, kts,kte)

      !
      !--- rainfall evap below cloud base
      !
      if(use_rebcb == 1)                                                               &
         call rain_evap_below_cloudbase(cumulus,itf,ktf, its,ite, kts,kte,ierr,kbcon,ktop&
         ,xmb,psur,xland,qo_cup,t_cup                      &
         ,po_cup,qes_cup,pwavo,edto,pwevo,pwo,pwdo         &
         ,pre,prec_flx,evap_flx,outt,outq,outbuoy,evap_bcb)


      !
      !--- includes effects of the remained cloud dissipation into the enviroment
      !
      if(use_cloud_dissipation >= 0.)                                              &
         call cloud_dissipation(cumulus,itf,ktf, its,ite, kts,kte,ierr,kbcon,ktop    &
         ,dtime,xmb,xland,qo_cup,qeso_cup,po_cup,outt,outq     &
         ,outqc,zuo,vvel2d,rho_hydr,qrco,sig,tempco,qco,tn_cup &
         ,heso_cup,zo)

      !
      !--- get the total (deep+congestus) evaporation flux for output (units kg/kg/s)
      !
      do i=its,itf
         if(ierr(i) /= 0) cycle
         do k=kts,ktop(i)
            dp=100.*(po_cup(i,k)-po_cup(i,k+1))
             !--- add congestus and deep plumes, and convert to kg/kg/s
            revsu_gf(i,k) = revsu_gf(i,k) + evap_flx(i,k)*g/dp
         enddo
      enddo
      !
      !
      !--- get lightning flashes density (parameterization from Lopez 2016, MWR)
      !
      if( lightning_diag == 1 .and. trim(cumulus) == 'deep') then
         call cup_up_cape(cape,z,zu,dby,gamma_cup,t_cup,k22,kbcon,ktop,ierr &
            ,tempco,qco,qrco,qo_cup,itf,ktf,its,ite, kts,kte)

         call cup_up_lightning(itf,ktf,its,ite, kts,kte, ierr, kbcon,ktop,xland,cape &
            ,zo,zo_cup,t_cup,t,tempco,qrco,po_cup,rho,prec_flx     &
            ,lightn_dens)
      endif
      !
      !--- for outputs (only deep plume)
      !
      if( trim(cumulus) == 'deep') then
         do i=its,itf
            if(ierr(i) /= 0) cycle
            var2d(i) = p_cwv_ave(i)
            do k=kts,ktop(i)+1
               prfil_gf  (i,k) = prec_flx(i,k)
               var3d_agf (i,k) = vvel2d  (i,k)
            enddo
         enddo
      endif
      !
      !--- for tracer convective transport / outputs
      !
      do i=its,itf
         if(ierr(i) /= 0) cycle
         do k=kts,ktf
            !clwup5d     (i,k) = qrco (i,k) !ice/liquid water
            !tup         (i,k) = (1./cp)*(hco(i,k)-g*zo_cup(i,k)-xlv*qco(i,k))!in-updraft temp
            tup          (i,k) = tempco(i,k) !in-updraft temp
         enddo
         tup (i,kte) = t_cup(i,kte)
      enddo

      !--- convert mass fluxes, etc...
      do i=its,itf
         if(ierr(i) /= 0) cycle
         pwavo       (i)   = xmb(i)*pwavo       (i)
         pwevo       (i)   = xmb(i)*pwevo       (i)
         zuo         (i,:) = xmb(i)*zuo         (i,:)
         zdo         (i,:) = xmb(i)*zdo         (i,:)
         pwo         (i,:) = xmb(i)*pwo         (i,:)
         pwdo        (i,:) = xmb(i)*pwdo        (i,:)
         up_massentro(i,:) = xmb(i)*up_massentro(i,:)
         up_massdetro(i,:) = xmb(i)*up_massdetro(i,:)
         dd_massentro(i,:) = xmb(i)*dd_massentro(i,:)
         dd_massdetro(i,:) = xmb(i)*dd_massdetro(i,:)
         zenv        (i,:) = xmb(i)*zenv        (i,:)
      enddo

      !--for output only.
      do i=its,itf
         subten_Q    (i,:) = xmb(i)*subten_Q    (i,:)
         subten_H    (i,:) = xmb(i)*subten_H    (i,:)
         subten_T    (i,:) = xmb(i)*subten_T    (i,:)
      enddo

      !
      !--- outputs a model sounding for the stand-alone code (part 2)
      !
      if(OUTPUT_SOUND == 1) then
         call SOUND(2,cumulus,int_time,dtime,ens4,itf,ktf,its,ite, kts,kte,xlats,xlons,jcol,whoami_all &
            ,z ,qes ,he ,hes ,t ,q ,po,z1 ,psur,zo,qeso,heo,heso,tn,qo,us,vs ,omeg,xz     &
            ,h_sfc_flux,le_sfc_flux,tsur, dx,stochastic_sig,zws,ztexec,zqexec, xland      &
            ,kpbl,k22,klcl,kbcon,ktop,aa0,aa1,sig,xaa0,hkb,xmb,pre,edto                   &
            ,zo_cup,dhdt,rho,zuo,zdo,up_massentro,up_massdetro,outt, outq,outqc,outu,outv)
      endif

      !--- for output only
      if( trim(cumulus) == 'deep') then
            AA1_       (:) = AA1       (:)
            AA0_       (:) = AA0       (:)
            AA1_RADPBL_(:) = AA1_RADPBL(:)
            
            if(dicycle == 2) then 
              AA1_ADV_ (:) = Q_adv     (:)
            else
            ! AA1_ADV_ (:) = AA1_ADV    (:)
            ! AA1_ADV_ (:) = wlpool_bcon(:)
              AA1_ADV_ (:) = vshear (:)
            ! AA1_ADV_ (:) = depth_neg_buoy (:)
            ! AA1_ADV_ (:) = cin1 (:)
            endif
            do i=its,itf
               if(ierr(i) == 0) cycle
               kbcon(i)=1
               ktop(i)=1
               klcl(i)=1
               jmin(i)=1
               k22(i)=1
            enddo
      endif

      if(liq_ice_number_conc == 1) then
         call get_liq_ice_number_conc(itf,ktf,its,ite, kts,kte,ierr,ktop &
            ,dtime,rho,outqc,tempco,outnliq,outnice)
      endif
      !
      !
      !--------------------------------------------------------------------------------------------!
      !- section for atmospheric composition
      !--------------------------------------------------------------------------------------------!
      if(use_tracer_transp==1)  then

         !--only for debug
         if(use_gate) then
            if(jl==1) then
               se_chem_update(1,:,:) = se_chem(1,:,:)
            else
               se_chem(1,:,:) = se_chem_update(1,:,:)
            endif
            do i=its,itf
               if(ierr(i) /= 0) cycle
               massi=0.
               do k=kts,ktop(i)
                  dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                  massi  =   massi + se_chem(1,i,k)*dp/g
               enddo
            enddo
         endif
         !--only for debug

         !-1) get mass mixing ratios at the cloud levels

         call cup_env_clev_chem(mtp,se_chem,se_cup_chem,ierr,itf,ktf,its,ite, kts,kte)

         !-2) determine in-cloud tracer mixing ratios
         !
         ! a) chem - updraft
         !- note: here "sc_up_chem" stores the total in-cloud tracer mixing ratio (i.e., including the portion
         !        embedded in the condensates).
         call get_incloud_sc_chem_up(cumulus,FSCAV,mtp,se_chem,se_cup_chem,sc_up_chem,pw_up_chem,tot_pw_up_chem      &
            ,zo_cup,rho,po,po_cup,qrco,tempco,pwo,zuo,up_massentro,up_massdetro               &
            ,vvel2d,vvel1d,start_level,k22,kbcon,ktop,klcl,ierr,xland,itf,ktf,its,ite, kts,kte)

         ! b) chem - downdraft
         call get_incloud_sc_chem_dd(cumulus,FSCAV,mtp,se_chem,se_cup_chem,sc_dn_chem,pw_dn_chem ,pw_up_chem,sc_up_chem &
            ,tot_pw_up_chem,tot_pw_dn_chem                                                       &
            ,zo_cup,rho,po_cup,qrcdo,pwdo,pwevo,edto,zdo,dd_massentro,dd_massdetro ,pwavo,pwo     &
            ,jmin,ierr,itf,ktf,its,ite, kts,kte)
         !
         !-3) determine the vertical transport including mixing, scavenging and evaporation
         !
         !---a) change per unit mass that a model cloud would modify the environment
         do i=its,itf
            if(ierr(i) /= 0) cycle

            !- flux form + source/sink terms + time explicit + FCT
            if(USE_FLUX_FORM == 1 .and. alp1 == 0. ) then

               if(use_fct == 0 ) then
                  do k=kts,ktop(i)
                     dp=100.*(po_cup(i,k)-po_cup(i,k+1))

                     out_chem(:,i,k) =-(zuo(i,k+1)*(sc_up_chem(:,i,k+1)-se_cup_chem(:,i,k+1) ) -                 &
                                        zuo(i,k  )*(sc_up_chem(:,i,k  )-se_cup_chem(:,i,k  ) ))*g/dp             &
                                      +(zdo(i,k+1)*(sc_dn_chem(:,i,k+1)-se_cup_chem(:,i,k+1) ) -                 &
                                        zdo(i,k  )*(sc_dn_chem(:,i,k  )-se_cup_chem(:,i,k  ) ))*g/dp*edto(i)
                  enddo

               else

                  !-- FCT scheme for the subsidence transport: d(M_env*S_env)/dz
                  sub_tend  = 0.
                  trcflx_in = 0.
                  dtime_max = dtime
                  massflx (i,:)=0.

                  do k=kts+1,ktop(i)+1
                     dp           = 100.*(po_cup(i,k)-po_cup(i,k+1))
                     trcflx_in (:,k) =-(zuo(i,k)  -edto(i)*zdo(i,k))*se_cup_chem(:,i,k) !* xmb(i)
                     massflx   (i,k) =-(zuo(i,k)  -edto(i)*zdo(i,k))           !* xmb(i)
                     dtime_max=min(dtime_max,.5*dp)
                  enddo
                  !- if dtime_max<dtime => needs a loop to update from t to t+dtime (check this!)
                  !if( dtime_max < dtime ) stop "dtime_max < dtime in GF scheme"

                  do ispc=1,mtp
                     call fct1d3 (ktop(i),kte,dtime_max,po_cup(i,:),se_chem(ispc,i,:),massflx(i,:),trcflx_in(ispc,:),sub_tend(ispc,:))
                  enddo

                  do k=kts,ktop(i)
                     dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                     out_chem(:,i,k) = -(zuo(i,k+1)*(sc_up_chem(:,i,k+1)) - zuo(i,k)*(sc_up_chem(:,i,k)))*g/dp     &
                                       +(zdo(i,k+1)*(sc_dn_chem(:,i,k+1)) - zdo(i,k)*(sc_dn_chem(:,i,k)))*g/dp*edto(i)

                     !- update with the subsidence term from FCT scheme
                     out_chem(:,i,k) = out_chem(:,i,k) + sub_tend(:,k)

                  enddo
               endif

                !- include evaporation (this term must not be applied to the tracer 'QW')
               if(USE_TRACER_EVAP == 1 .and. cumulus /= 'shallow') then
                  do k=kts,ktop(i)
                     dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                     out_chem(:,i,k) = out_chem(:,i,k)    &
                        - 0.5*edto(i)*(zdo(i,k)*pw_dn_chem(:,i,k)+zdo(i,k+1)*pw_dn_chem(:,i,k+1))*g/dp !&  ! evaporated ( pw_dn < 0 => E_dn > 0)
                                       !*chem_name_mask_evap(:) !-- to avoid the "Dry Mass Violation"
                  enddo
               endif

               !- include scavenging
               if(USE_TRACER_SCAVEN > 0 .and. cumulus /= 'shallow') then
                  do k=kts,ktop(i)
                     dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                     out_chem(:,i,k) = out_chem(:,i,k) &
                        - 0.5*(zuo(i,k)*pw_up_chem(:,i,k)+zuo(i,k+1)*pw_up_chem(:,i,k+1))*g/dp  ! incorporated in rainfall (<0)
                  enddo
               endif
            endif ! IF(USE_FLUX_FORM == 1 .and. alp1 == 0. )

            !- flux form + source/sink terms + time explicit/implicit + upstream
            if(USE_FLUX_FORM == 1 .and. alp1 > 0. ) then

               alp0=1.-alp1
               do k=kts,ktop(i)+1
                  fp(k) = 0.5*(zenv(i,k)+abs(zenv(i,k)))
                  fm(k) = 0.5*(zenv(i,k)-abs(zenv(i,k)))
               enddo

               do k=kts,ktop(i)
                  dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                  beta1 = dtime*g/dp
                  aa(k) =    alp1*beta1*fm(k)
                  bb(k) = 1.+alp1*beta1*(fp(k)-fm(k+1))
                  cc(k) =   -alp1*beta1*fp(k+1)

                  ddtr(:,k) = se_chem(:,i,k) - (zuo(i,k+1)*sc_up_chem(:,i,k+1) - zuo(i,k)*sc_up_chem(:,i,k))*beta1      &
                                             + (zdo(i,k+1)*sc_dn_chem(:,i,k+1) - zdo(i,k)*sc_dn_chem(:,i,k))*beta1*edto(i)

                  !- include evaporation (this term must not be applied to the tracer 'QW')
                  if(USE_TRACER_EVAP == 1 .and. cumulus /= 'shallow') then
                     out_chem(:,i,k) = out_chem(:,i,k)    &
                        - 0.5*edto(i)*(zdo(i,k)*pw_dn_chem(:,i,k)+zdo(i,k+1)*pw_dn_chem(:,i,k+1))*beta1 !&  ! evaporated ( pw_dn < 0 => E_dn > 0)
                                    !*chem_name_mask_evap(:) !-- to avoid the "Dry Mass Violation"
                  endif

                  !- include scavenging
                  if(USE_TRACER_SCAVEN > 0 .and. cumulus /= 'shallow') then
                     out_chem(:,i,k) = out_chem(:,i,k) &
                        - 0.5*(zuo(i,k)*pw_up_chem(:,i,k)+zuo(i,k+1)*pw_up_chem(:,i,k+1))*beta1  ! incorporated in rainfall (<0)
                  endif

                  ddtr(:,k) = ddtr(:,k) + out_chem(:,i,k) + &
                     alp0*beta1*(-fm(k)*se_chem(:,i,max(kts,k-1)) +(fm(k+1)-fp(k))*se_chem(:,i,k) +fp(k+1)*se_chem(:,i,k+1))

               enddo
               do ispc = 1, mtp
                  if(CHEM_NAME_MASK (ispc) == 0 ) cycle
                  call tridiag (ktop(i),aa (kts:ktop(i)),bb (kts:ktop(i)),cc (kts:ktop(i)),ddtr (ispc,kts:ktop(i)))
                  out_chem(ispc,i,kts:ktop(i))=(ddtr(ispc,kts:ktop(i))-se_chem(ispc,i,kts:ktop(i)))/dtime
               enddo
            endif !USE_FLUX_FORM == 1 .and. alp1 > 0.

            !- flux form + source/sink terms + time explicit + upstream with anti-diffusion step (Smolarkiewicz 1983)
            if(USE_FLUX_FORM == 2 .or. USE_FLUX_FORM == 3) then
               if(USE_FLUX_FORM == 2)  lstep = -1 ! upstream + anti-diffusion step
               if(USE_FLUX_FORM == 3)  lstep =  1 ! only upstream
               alp0=1.

               if(ierr(i) /= 0) cycle
               !--- Zenv here have the following reference:  < 0 => downward motion
               zenv(i,:) = -(zuo(i,:)-edto(i)*zdo(i,:))

               do istep=1,lstep,-2

                  if(istep == 1) then
                     ddtr_upd(:,:) = se_chem(:,i,:)
                     do k=kts,ktop(i)+1
                        fp_mtp(:,k) = 0.5*(zenv(i,k)+abs(zenv(i,k)))
                        fm_mtp(:,k) = 0.5*(zenv(i,k)-abs(zenv(i,k)))
                     enddo
                  else
                     ddtr_upd(:,kts:ktop(i)+1) = se_chem(:,i,kts:ktop(i)+1) + out_chem(:,i,kts:ktop(i)+1)*dtime
                     zenv_diff(:,kts) = 0.
                     do k=kts,ktop(i)+1
                        dz =  zo_cup(i,k+1)-zo_cup(i,k)
                        zenv_diff (:,k+1) = 1.08*( dz*abs(zenv(i,k+1)) - dtime * zenv(i,k+1)**2 ) &
                           * (ddtr_upd(:,k+1) - ddtr_upd(:,k))               &
                           /((ddtr_upd(:,k+1) + ddtr_upd(:,k) + 1.e-16)*dz)
                     enddo
                     do k=kts,ktop(i)+1
                        fp_mtp(:,k) = 0.5*(zenv_diff(:,k)+abs(zenv_diff(:,k)))
                        fm_mtp(:,k) = 0.5*(zenv_diff(:,k)-abs(zenv_diff(:,k)))
                     enddo
                  endif

                  do k=kts,ktop(i)
                     dp=-100.*(po_cup(i,k)-po_cup(i,k+1))
                     beta1 = dtime*g/dp
                     ddtr(:,k) = ddtr_upd(:,k) + alp0*beta1*( &
                         (fp_mtp(:,k+1)*ddtr_upd(:,k)           +fm_mtp(:,k+1)*ddtr_upd(:,k+1)) &
                        -(fp_mtp(:,k  )*ddtr_upd(:,max(kts,k-1))+fm_mtp(:,k  )*ddtr_upd(:,k  )) )
                  enddo
                  do ispc = 1, mtp
                     if(CHEM_NAME_MASK (ispc) == 0 ) cycle
                     out_chem(ispc,i,kts:ktop(i))=(ddtr(ispc,kts:ktop(i))-se_chem(ispc,i,kts:ktop(i)))/dtime
                  enddo

               enddo ! anti-diff steps

               do k=kts,ktop(i)
                  dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                  beta1 = g/dp

                  out_chem(:,i,k)  =   out_chem(:,i,k)                                                           &
                     - (zuo(i,k+1)*sc_up_chem(:,i,k+1) - zuo(i,k)*sc_up_chem(:,i,k))*beta1      &
                     + (zdo(i,k+1)*sc_dn_chem(:,i,k+1) - zdo(i,k)*sc_dn_chem(:,i,k))*beta1*edto(i)

                  !- include evaporation (this term must not be applied to the tracer 'QW')
                  if(USE_TRACER_EVAP == 1 .and. cumulus /= 'shallow') then
                     out_chem(:,i,k) = out_chem(:,i,k)    &
                        - 0.5*edto(i)*(zdo(i,k)*pw_dn_chem(:,i,k)+zdo(i,k+1)*pw_dn_chem(:,i,k+1))*beta1 !&  ! evaporated ( pw_dn < 0 => E_dn > 0)
                                    !*chem_name_mask_evap(:) !-- to avoid the "Dry Mass Violation"
                  endif

                  !- include scavenging
                  if(USE_TRACER_SCAVEN > 0 .and. cumulus /= 'shallow') then
                     out_chem(:,i,k) = out_chem(:,i,k) &
                        - 0.5*(zuo(i,k)*pw_up_chem(:,i,k)+zuo(i,k+1)*pw_up_chem(:,i,k+1))*beta1  ! incorporated in rainfall (<0)
                  endif
               enddo
            endif ! USE_FLUX_FORM == 2 .or. USE_FLUX_FORM == 3

            !--- check mass conservation for tracers
            do ispc = 1, mtp
               if(CHEM_NAME_MASK (ispc) == 0 ) cycle
               trash_ (:) = 0.
               trash2_(:) = 0.
               evap_  (:) = 0.
               wetdep_(:) = 0.
               residu_(:) = 0.
               do k=kts,ktop(i)
                  dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                  evap   =  -0.5*(zdo(i,k)*pw_dn_chem(ispc,i,k)+zdo(i,k+1)*pw_dn_chem(ispc,i,k+1))*g/dp*edto(i)
                  wetdep =   0.5*(zuo(i,k)*pw_up_chem(ispc,i,k)+zuo(i,k+1)*pw_up_chem(ispc,i,k+1))*g/dp

                  evap_  (ispc) =   evap_  (ispc) + evap  *dp/g
                  wetdep_(ispc) =   wetdep_(ispc) + wetdep*dp/g
                  residu_(ispc) =   residu_(ispc) + (wetdep - evap)*dp/g

                  !trash_ (ispc) =   trash_ (ispc) + (out_chem (ispc,i,k) - evap + wetdep)*dp/g
                  trash_ (ispc) =   trash_ (ispc) + (out_chem (ispc,i,k)                )*dp/g

                  trash2_(ispc) =   trash2_(ispc) + se_chem(ispc,i,k)*dp/g
               enddo
               if(residu_(ispc) < 0.) then
                  beta1 = g/(po_cup(i,kts)-po_cup(i,ktop(i)+1))
                  do k=kts,ktop(i)
                     out_chem(ispc,i,k)=out_chem(ispc,i,k)+residu_(ispc)*beta1
                  enddo
               endif

               !if(evap_  (ispc) > wetdep_(ispc)) then
                !print*,"budget=",ispc,evap_  (ispc), wetdep_(ispc),trash_ (ispc),trim(CHEM_NAME(ispc))!,trash_ (ispc),trash2_(ispc)
                !call flush(6)
                !endif
               !if(evap_  (ispc) > wetdep_(ispc)) stop " eva<wet "
                !if(abs(trash_(ispc)) >1.e-6 ) then
               !  if (MAPL_AM_I_ROOT())  write(6,*)'=> mass_cons=',trash_(ispc),spacing(trash2_(ispc)),trim(CHEM_NAME(ispc)),trim(cumulus)
               !endif
            enddo

         enddo ! loop 'i'

         if(use_gate) then
            !--only for debug
            do i=its,itf
               if(ierr(i) /= 0) cycle
               massf = 0.
               do k=kts,ktop(i)
                  se_chem_update(ispc_CO,i,k) = se_chem_update(ispc_CO,i,k) + out_chem(ispc_CO,i,k) * dtime
                  !se_chem_update(ispc_CO,i,k) = max(0.,se_chem_update(ispc_CO,i,k))
                  dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                  evap_  (ispc_CO)=  -0.5*(zdo(i,k)*pw_dn_chem(ispc_CO,i,k)+zdo(i,k+1)*pw_dn_chem(ispc_CO,i,k+1))*g/dp*edto(i)
                  wetdep_(ispc_CO)=   0.5*(zuo(i,k)*pw_up_chem(ispc_CO,i,k)+zuo(i,k+1)*pw_up_chem(ispc_CO,i,k+1))*g/dp
                  massf = massf + se_chem_update(1,i,k)*dp/g + (- evap_(ispc_CO) + wetdep_(ispc_CO))*dp/g
               enddo
               if(abs((massf-massi)/(1.e-12+massi))>1.e-6) print*,"mass con=>",(massf-massi)/(1.e-12+massi)
            enddo
19          format(1x,I3,1x,5e14.3)
18          format(1x,I3,1x,4e14.3)
20          format(1x,I3,1x,11e16.6)
         !--only for debug
         endif

      !--------------------------------------------------------------------------------------------!
      endif !- end of section for atmospheric composition
      !--------------------------------------------------------------------------------------------!


      !
      !IF(use_gate) THEN
      !   do i=its,itf
      !  massf = 0.
      !  if(ierr(i) /= 0) cycle
      !  do k=kts,ktop(i)
      !     se_chem_update(1,i,k) = se_chem_update(1,i,k) + outmpql(lsmp,i,k)* dtime
      !     se_chem_update(2,i,k) = se_chem_update(2,i,k) + outmpqi(lsmp,i,k)* dtime
      !     se_chem_update(3,i,k) = se_chem_update(3,i,k) + outmpcf(lsmp,i,k)* dtime
      !     mpql (lsmp,i,k)    = se_chem_update(1,i,k)
      !     mpqi (lsmp,i,k)    = se_chem_update(2,i,k)
      !     mpcf (lsmp,i,k)    = se_chem_update(3,i,k)
      !  enddo
      !   enddo
      !ENDIF


      !- for debug/diag
          if(cumulus == 'deep') then
            do i=its,itf
               !if(ierr(i) /= 0) cycle
               aa0_    (i)  = aa0      (i)
               aa1_    (i)  = aa1      (i)
               aa1_bl_ (i)  = aa1_bl   (i)
               tau_bl_ (i)  = tau_bl   (i)
               tau_ec_ (i)  = tau_ecmwf(i)
            !  if(use_memory == 0) tau_ec_ (i)  = x_add_buoy(i)
            enddo
          endif

      !
      !- begin: for GATE soundings-------------------------------------------
      if(use_gate .or. wrtgrads) then
         if(cumulus == 'deep'   ) then
            cty='1'
            nvarbegin =  0
         endif
         if(cumulus == 'shallow') then
            cty='2'
            nvarbegin =101
         endif
         if(cumulus == 'mid'    ) then
            cty='3'
            nvarbegin =201
         endif
         do i=its,itf
            !if(ierr(i).eq.0) then
            !- 2-d section
            do k=kts,ktf  !max(1,ktop(i))
               nvar=nvarbegin

               if(cumulus == 'deep'   ) &
                  call set_grads_var(jl,k,nvar,zo(i,k),"zo"//cty ,' height','3d')
               !        call set_grads_var(jl,k,nvar,po(i,k),"po"//cty ,' press','3d')

               dp=100.*(po_cup(i,k)-po_cup(i,k+1))
               E_dn  = -0.5*(pwdo(i,k)+pwdo(i,k+1))*g/dp*edto(i)*86400.*xlv/cp*xmb(i) ! pwdo < 0 and E_dn must > 0
               C_up  = dellaqc(i,k)+(zuo(i,k+1)* qrco(i,k+1) - zuo(i,k  )* qrco(i,k  )  )*g/dp &
                  +0.5*(pwo (i,k)+pwo (i,k+1))*g/dp
               C_up = - C_up*86400.*xlv/cp*xmb(i)

               trash =-(zuo(i,k+1)*(qco (i,k+1)-qo_cup(i,k+1) ) -                 &
                  zuo(i,k  )*(qco (i,k  )-qo_cup(i,k  ) ) )*g/dp
               trash2=+(zdo(i,k+1)*(qcdo(i,k+1)-qo_cup(i,k+1) ) -                 &
                  zdo(i,k  )*(qcdo(i,k  )-qo_cup(i,k  ) ) )*g/dp*edto(i)

               trash  = trash *86400.*xlv/cp*xmb(i)
               trash2 = trash2*86400.*xlv/cp*xmb(i)

               env_mf = 0.5* ((zuo(i,k+1)-zdo(i,k+1)*edto(i)) + (zuo(i,k)-zdo(i,k)*edto(i)))
               resten_H = dellah(i,k) - subten_H(i,k)
               resten_Q = dellaQ(i,k) - subten_Q(i,k)
               resten_T =(1./cp)*(resten_H-xlv*resten_Q)
               !trash2 = qco   (i,k  )! zuo(i,k+1)*(qco (i,k+1)-qo_cup(i,k+1) ) !*g/dp
               !trash  = qo_cup(i,k  )! zuo(i,k  )*(qco (i,k  )-qo_cup(i,k  ) ) !*g/dp
               trash2 = zuo(i,k+1)*(qco (i,k+1)-qo_cup(i,k+1) )*1000 !*g/dp
               trash  = zuo(i,k  )*(qco (i,k  )-qo_cup(i,k  ) )*1000  !*g/dp

               call set_grads_var(jl,k,nvar,out_chem(1,i,k)*86400,"outchem"//cty ,' outchem','3d')
               call set_grads_var(jl,k,nvar,sc_up_chem(1,i,k),"scup"//cty ,' sc_chem','3d')
               call set_grads_var(jl,k,nvar,sc_dn_chem(1,i,k),"scdn"//cty ,' sc_chem','3d')
               call set_grads_var(jl,k,nvar,massi,"mi"//cty ,' initial mass','2d')
               call set_grads_var(jl,k,nvar,massf,"mf"//cty ,' final mass','2d')
               call set_grads_var(jl,k,nvar,se_chem(1,i,k),"se"//cty ,' se_chem','3d')
               call set_grads_var(jl,k,nvar,se_cup_chem(1,i,k),"secup"//cty ,' se_cup_chem','3d')
               !-- only for debug
               !call set_grads_var(jl,k,nvar,se_chem_update(1,i,k),"newse"//cty ,' new se_chem','3d')
               if(APPLY_SUB_MP == 1) then
                  kmp=lsmp
                  call set_grads_var(jl,k,nvar,outmpqi(kmp,i,k)*86400*1000,"outqi"//cty ,' outmpqi','3d')
                  call set_grads_var(jl,k,nvar,outmpql(kmp,i,k)*86400*1000,"outql"//cty ,' outmpql','3d')
                  call set_grads_var(jl,k,nvar,outmpcf(kmp,i,k)*86400,"outcf"//cty ,' outmpcf','3d')

                  call set_grads_var(jl,k,nvar,mpqi(kmp,i,k),"mpqi"//cty ,' mpqi','3d')
                  call set_grads_var(jl,k,nvar,mpql(kmp,i,k),"mpql"//cty ,' mpql','3d')
                  call set_grads_var(jl,k,nvar,mpcf(kmp,i,k),"mpcf"//cty ,' mpcf','3d')
               endif
               call set_grads_var(jl,k,nvar,env_mf,"sub"//cty ,' sub','3d')


               if(LIQ_ICE_NUMBER_CONC == 1) then
                  call set_grads_var(jl,k,nvar,outnice(i,k)*86400.,"outnice"//cty ,'out # ice1/day','3d')
                  call set_grads_var(jl,k,nvar,outnliq(i,k)*86400.,"outnliq"//cty ,'out # liq /day','3d')
               endif
               call set_grads_var(jl,k,nvar,zuo(i,k)/xmb(i),"zup"//cty,'norm m flux up ','3d')
               call set_grads_var(jl,k,nvar,zdo(i,k)/xmb(i),"zdn"//cty,'norm m flux dn ','3d')
               call set_grads_var(jl,k,nvar,zenv(i,k),"zenv"//cty,'norm m flux env ','3d')

               call set_grads_var(jl,k,nvar,-edto(i)*xmb(i)*zdo(i,k),"mdn"//cty ,'m flux down (kg/s/m^2)','3d')
               call set_grads_var(jl,k,nvar,up_massentro(i,k),"upent"//cty ,'up_massentr(kg/s/m^2)','3d')
               call set_grads_var(jl,k,nvar,xmb(i)*up_massdetro(i,k),"updet"//cty ,'up_massdetr(kg/s/m^2)','3d')
               call set_grads_var(jl,k,nvar,outt(i,k)*86400.,"outt"//cty ,'outt K/day','3d')

               call set_grads_var(jl,k,nvar,resten_T*86400.,     "rest"//cty ,'residuo T K/day','3d')
               call set_grads_var(jl,k,nvar,resten_H*86400./cp,    "resh"//cty ,'residuo H J/kg/day','3d')
               call set_grads_var(jl,k,nvar,resten_Q*86400.*xlv/cp,"resq"//cty ,'residuo q K/day   ','3d')
               call set_grads_var(jl,k,nvar,subten_T(i,k)*86400.,       "subt"//cty ,'subT K/day','3d')
               call set_grads_var(jl,k,nvar,subten_H(i,k)*86400./cp,    "subh"//cty ,'subH J/kg/day','3d')
               call set_grads_var(jl,k,nvar,subten_Q(i,k)*86400.*xlv/cp,"subq"//cty ,'subq K/day   ','3d')

               call set_grads_var(jl,k,nvar,outq(i,k)*86400.*xlv/cp,"outq"//cty ,'outq K/s','3d')
               call set_grads_var(jl,k,nvar,outqc(i,k)*86400.*xlv/cp,"outqc"//cty ,'outqc K/day','3d')
               call set_grads_var(jl,k,nvar,pre(i)*3600.,"precip"//cty ,'precip mm','2d')
               call set_grads_var(jl,k,nvar,prec_flx(i,k)*3600.,"precflx"//cty ,'prec flx mm','3d')
               call set_grads_var(jl,k,nvar,pwo(i,k),"pwo"//cty ,' xx','3d')
               call set_grads_var(jl,k,nvar,outu(i,k)*86400.,"outu"//cty ,'out_U m/s/day','3d')
               call set_grads_var(jl,k,nvar,outv(i,k)*86400.,"outv"//cty ,'out_V m/s/day','3d')
               call set_grads_var(jl,k,nvar,xmb(i),"xmb"//cty ,'xmb kg/m2/s','2d')
               call set_grads_var(jl,k,nvar,vvel2d(i,k),"W2d"//cty ,'W /m/s','3d')
               call set_grads_var(jl,k,nvar,vvel1d(i),"W1d"//cty ,'W1s /m/s','2d')
               call set_grads_var(jl,k,nvar,us(i,k),"us"//cty ,'U /m/s','3d')
               call set_grads_var(jl,k,nvar,outu(i,k)*86400./(1.e-16+xmb(i)),"delu"//cty ,'dellu','3d')
               call set_grads_var(jl,k,nvar,evap_bcb(i,k)*1000.,"evcb"//cty ,'g/kg','3d')

               call set_grads_var(jl,k,nvar,tot_pw_up_chem(1,i),"pwup"//cty ,'pwup','2d')
               call set_grads_var(jl,k,nvar,tot_pw_dn_chem(1,i),"pwdn"//cty ,'pwdn','2d')
               !----
               !----
               call set_grads_var(jl,k,nvar,xmb(i)*dellah(i,k)*86400./cp,"delh"//cty ,'dellah K/day','3d')
               call set_grads_var(jl,k,nvar,xmb(i)*dellaq(i,k)*86400.*xlv/cp, "dellq"//cty ,'dellaq K/day','3d')
               call set_grads_var(jl,k,nvar,xmb(i)*dellaqc(i,k)*86400.*xlv/cp,"dellqc"//cty ,'dellaqc K/day','3d')
               call set_grads_var(jl,k,nvar,xmb(i),"xmb"//cty,'m flux up (kg/s/m^2)','2d')
               call set_grads_var(jl,k,nvar,aa1(i),"aa1"//cty,'AA1 J/kg3)','2d')
               call set_grads_var(jl,k,nvar,float(ierr(i)),"ierr"//cty ,'ierr #','2d')
               call set_grads_var(jl,k,nvar,xmb(i)*dd_massentro(i,k),"ddent"//cty ,'dd_massentr(kg/s/m^2)','3d')
               call set_grads_var(jl,k,nvar,xmb(i)*dd_massdetro(i,k),"dddet"//cty ,'dd_massdetr(kg/s/m^2)','3d')
               !!      go to 333
               call set_grads_var(jl,k,nvar,hc(i,k),"hc"//cty ,' hc','3d')
               call set_grads_var(jl,k,nvar,hco(i,k),"hco"//cty ,' hco','3d')
               call set_grads_var(jl,k,nvar,dby(i,k),"dby"//cty ,' dbuo','3d')
               !call set_grads_var(jl,k,nvar,QCUP(i,k),"qcup"//cty ,'C_UP','3d')
               call set_grads_var(jl,k,nvar,t_cup(i,k)-273.15,"te"//cty ,' K','3d')
               call set_grads_var(jl,k,nvar,1000.*q_cup(i,k),"qe"//cty ,' kg kg-1','3d')
               call set_grads_var(jl,k,nvar,he_cup(i,k),"he"//cty ,' he','3d')
               call set_grads_var(jl,k,nvar,HKB(i),"hkb"//cty ,' H','2d')
               call set_grads_var(jl,k,nvar,HKB(i),"hkb"//cty ,' H','2d')
               call set_grads_var(jl,k,nvar,1000.*zqexec(i),"qex"//cty ,' qex','2d')
               call set_grads_var(jl,k,nvar,z_cup(i,max(1,k22  (i))),"zs"//cty ,' m','2d')
               call set_grads_var(jl,k,nvar,z_cup(i,max(1,kbcon(i))),"zbcon"//cty ,' m','2d')
               call set_grads_var(jl,k,nvar,z_cup(i,max(1,ktop (i))),"ztop"//cty ,' m','2d')
               call set_grads_var(jl,k,nvar,z_cup(i,max(1,klcl (i))),"zlcl"//cty ,' m','2d')
               call set_grads_var(jl,k,nvar,z_cup(i,max(1,jmin (i))),"zjmin"//cty ,' m','2d')
               call set_grads_var(jl,k,nvar,zws(i),"ws"//cty ,' m/s','2d')
               call set_grads_var(jl,k,nvar,clfrac(i,k),"clfrac"//cty ,'shcf #','3d')
               call set_grads_var(jl,k,nvar,entr_rate_2d(i,k),"entr"//cty ,' m-1','3d')
               call set_grads_var(jl,k,nvar,cd(i,k),"detr"//cty ,' m-1','3d')
               call set_grads_var(jl,k,nvar,pwdo(i,k),"pwd"//cty ,' xx','3d')
               call set_grads_var(jl,k,nvar,edto(i),"edt"//cty ,'edt kg/m2/s','2d')
               call set_grads_var(jl,k,nvar,E_DN,"EVAP"//cty ,' xx','3d')
               call set_grads_var(jl,k,nvar,C_UP,"CUP"//cty ,' xx','3d')
               !       call set_grads_var(jl,k,nvar,trash,"TUP"//cty ,' xx','3d')
               !       call set_grads_var(jl,k,nvar,trash2,"TDN"//cty ,' xx','3d')
               call set_grads_var(jl,k,nvar,trash,"F1"//cty ,' F1','3d')
               call set_grads_var(jl,k,nvar,trash2,"F2"//cty ,' F2','3d')

               call set_grads_var(jl,k,nvar,p_liq_ice(i,k),"pli"//cty ,'#','3d')
               call set_grads_var(jl,k,nvar,melting_layer(i,k),"cpli"//cty ,'#','3d')
               call set_grads_var(jl,k,nvar,t(i,k),"t"//cty ,'temp K','3d')
               call set_grads_var(jl,k,nvar,tn(i,k),"tn"//cty ,'temp K','3d')
               call set_grads_var(jl,k,nvar,1000.*q(i,k),"q"//cty ,'q g/kg','3d')
               call set_grads_var(jl,k,nvar,1000.*qo(i,k),"qn"//cty ,'q g/kg','3d')
               call set_grads_var(jl,k,nvar,1000.*qrco(i,k),"qrc"//cty ,'q g/kg','3d')
               call set_grads_var(jl,k,nvar,1000.*(q(i,k)+outq(i,k)*dtime),"qnc"//cty ,'q upd conv g/kg','3d')
               call set_grads_var(jl,k,nvar,1000.*(qo(i,k)+outq(i,k)*dtime),"qnall"//cty ,'q upd all g/kg','3d')
               call set_grads_var(jl,k,nvar,1000.*qrr(i,k),"qrr"//cty ,'qrr g/kg','3d')
               call set_grads_var(jl,k,nvar,1000.*qco(i,k),"qc"//cty ,'qc g/kg','3d')
               call set_grads_var(jl,k,nvar,1000.*qo_cup(i,k),"qcup"//cty ,'qc g/kg','3d')
               call set_grads_var(jl,k,nvar,1000.*qeso_cup(i,k),"qescup"//cty ,'qc g/kg','3d')

                 !~ call set_grads_var(jl,k,nvar,aa0(i),"a0"//cty,'aa0','2d')
                 !~ call set_grads_var(jl,k,nvar,aa1_fa(i),"aa1fa"//cty,'aa1fa','2d')
                 !~ call set_grads_var(jl,k,nvar,aa1_bl(i),"aa1bl"//cty,'aa1bl','2d')
                 !~ call set_grads_var(jl,k,nvar,aa0_bl(i),"aa0bl"//cty,'aa0bl','2d')
                 !~ call set_grads_var(jl,k,nvar,aa1(i),"a1"//cty,'aa1','2d')
                 !~ call set_grads_var(jl,k,nvar,aa1(i)/(1.e-6+tau_ecmwf(i)),"mb13"//cty,'aa0','2d')
                 !~ call set_grads_var(jl,k,nvar,xaa0(i),"xa0"//cty,'xaa0','2d')
                 !~ call set_grads_var(jl,k,nvar,(XAA0(I)-AA1(I))/MBDT(I),"xk"//cty,'xk','2d')
333         continue
            enddo
            if(wrtgrads .and. .not. use_gate) then
               call wrt_bin_ctl(1,kte,po(1,1:kte),cumulus)
            endif
         enddo
      endif
   !- end  : for GATE soundings-------------------------------------------
   !
   !
   !-------------------------- not in use ------------------------------------------------------!
   !--- get cloud fraction
   !
   ! do i=its,itf
   !    clfrac(i,:)=0.
   !    if(ierr(i) /= 0) cycle
   !    dummy1(kts:ktf) = xmb(i)* zuo(i,kts:ktf)
   !    dummy2(kts:ktf) = 100.*po_cup(i,kts:ktf)
   !    call get_cloud_fraction(ktf,kts,ktf                                                   &
   !     ,dummy2(kts:ktf),zo_cup(i,kts:ktf),tn_cup(i,kts:ktf),qo_cup(i,kts:ktf) &
   !     ,qco (i,kts:ktf),  qrco(i,kts:ktf),  dummy1(kts:ktf),clfrac(i,kts:ktf) )
   ! enddo
   !--------------------------------------------------------------------------------------------!
   !

   end subroutine CUP_gf
   !------------------------------------------------------------------------------------
   subroutine cup_dd_edt(cumulus,ierr,us,vs,z,ktop,kbcon,edt,p,pwav, &
      pw,ccn,pwev,edtmax,edtmin,maxens2,edtc,psum2,psumh,    &
      ccnclean,rho,aeroevap,itf,ktf,ipr,jpr,its,ite, kts,kte,vshear )

      implicit none

      character *(*),        intent (in )  :: cumulus
      integer  ,intent (in  )             ::                           &
         ipr,jpr,aeroevap,itf,ktf,its,ite, kts,kte
      integer, intent (in   )              :: maxens2
      !
      ! ierr error value, maybe modified in this routine
      !
      real,    dimension (its:ite,kts:kte)                              &
         ,intent (in   )                   ::                           &
         rho,us,vs,z,p,pw
      real,    dimension (its:ite)                                      &
         ,intent (in   )                   ::                           &
         pwav,pwev,ccn,psum2,psumh,edtmax,edtmin
      real                                                              &
         ,intent (in   )                   ::                           &
         ccnclean
      integer, dimension (its:ite)                                      &
         ,intent (in   )                   ::                           &
         ktop,kbcon
      integer, dimension (its:ite)                                      &
         ,intent (inout)                   ::                           &
         ierr
      real,    dimension (its:ite,1:maxens2)                            &
         ,intent (out  )                   ::                           &
         edtc
      real,    dimension (its:ite)                                      &
         ,intent (out  )                   ::                           &
         edt,vshear
      !
      !  local variables in this routine
      !
      integer i,k,kk
      real    einc,pef,pefb,prezk,zkbc
      real,    dimension (its:ite)         ::   vws,sdp
      real :: pefc,aeroadd,rhoc,dp,prop_c
      real, parameter ::  alpha3 = 1.9 ,beta3  = -1.13

      !
      !--- DETERMINE DOWNDRAFT STRENGTH IN TERMS OF WINDSHEAR
      !
      ! */ calculate an average wind shear over the depth of the cloud
      !
      edt   =0.
      vws   =0.
      sdp   =0.
      vshear=0.
      edtc  =0.

      if(cumulus=='shallow') return

      do i=its,itf
         if(ierr(i) /= 0)cycle
         do kk = kbcon(i),ktop(i)
            dp = p(i,kk) - p(i,kk+1)
            vws(i) = vws(i) + (abs((us(i,kk+1)-us(i,kk))/(z(i,kk+1)-z(i,kk)))  + &
                               abs((vs(i,kk+1)-vs(i,kk))/(z(i,kk+1)-z(i,kk)))) * dp
            sdp(i) = sdp(i) + dp
         enddo
         vshear(i) = 1.e3 * vws(i) / sdp(i)
      end do

      do i=its,itf
         if(ierr(i) /= 0)cycle
         pef = (1.591-0.639*vshear(i)+0.0953*(vshear(i)**2) -0.00496*(vshear(i)**3))

         !print*,"shear=",vshear(i),pef,1-max(min(pef,0.9),0.1)
         pef = min(pef,0.9)
         pef = max(pef,0.1)
         edt(i) = 1.-pef

         !
         !--- cloud base precip efficiency
         !
         if(use_rebcb == 0) then 
           zkbc  = z(i,kbcon(i))*3.281e-3
           prezk = 0.02
           if(zkbc > 3.0) prezk=0.96729352+zkbc*(-0.70034167+zkbc* &
                        (0.162179896+zkbc*(- 1.2569798e-2+zkbc*(4.2772e-4-zkbc*5.44e-6))))
           if(zkbc > 25.) prezk=2.4
           pefb = 1./(1.+prezk)
           pefb = min(pefb,0.9)
           pefb = max(pefb,0.1)
           edt(i) = 1.-0.5*(pefb+pef)
         endif


         if(aeroevap.gt.1)then
            aeroadd=(ccnclean**beta3)*((psumh(i))**(alpha3-1)) !*1.e6
            !if(i.eq.ipr)write(0,*)'edt',ccnclean,psumh(i),aeroadd
            !prop_c=.9/aeroadd
            prop_c=.5*(pefb+pef)/aeroadd
            aeroadd=(ccn(i)**beta3)*((psum2(i))**(alpha3-1)) !*1.e6
            !if(i.eq.ipr)write(0,*)'edt',ccn(i),psum2(i),aeroadd,prop_c
            aeroadd=prop_c*aeroadd
            pefc=aeroadd
            if(pefc.gt.0.9)pefc=0.9
            if(pefc.lt.0.1)pefc=0.1
            EDT(I)=1.-pefc
            if(aeroevap.eq.2)EDT(I)=1.-.25*(pefb+pef+2.*pefc)
         endif

         !--- edt here is 1-precipeff!
         if(ZERO_DIFF==1) then
            edtc(i,1)=0.8*edt(i)
         else
            edtc(i,1)=    edt(i)
         endif
      enddo
      do i=its,itf
         if(ierr(i) /= 0)cycle
         edtc(i,1) = -edtc(i,1)*pwav(i)/pwev(i)
         edtc(i,1) = min(edtmax(i), edtc(i,1))
         edtc(i,1) = max(edtmin(i), edtc(i,1))
      enddo

   end subroutine cup_dd_edt
   !------------------------------------------------------------------------------------
   subroutine cup_dd_moisture(cumulus,ierrc,zd,hcd,hes_cup,qcd,qes_cup,                &
      pwd,q_cup,z_cup,dd_massentr,dd_massdetr,jmin,ierr,                       &
      gamma_cup,pwev,bu,qrcd, q,he,t_cup,iloop,t_wetbulb,q_wetbulb,qco,pwavo,  &
      itf,ktf,its,ite, kts,kte                              )

      implicit none

      character(len=*), intent(in) :: cumulus
      integer         , intent(in) :: itf,ktf,its,ite, kts,kte
      ! q       = environmental q on model levels
      ! q_cup   = environmental q on model cloud levels
      ! qes_cup = saturation q on model cloud levels
      ! hes_cup = saturation h on model cloud levels
      ! hcd = h in model cloud
      ! bu = buoancy term
      ! zd = normalized downdraft mass flux
      ! gamma_cup = gamma on model cloud levels
      ! mentr_rate = entrainment rate
      ! qcd  = cloud q (including liquid water) after entrainment
      ! qrch = saturation q in cloud
      ! pwd  = evaporate at that level
      ! pwev = total normalized integrated evaoprate (I2)
      ! entr = entrainment rate
      ! cdd  = detrainment function
      !
      real,    dimension (its:ite) ,intent (in  )          ::           &
         t_wetbulb,q_wetbulb, pwavo
      real,    dimension (its:ite,kts:kte) ,intent (in   ) ::           &
         zd,t_cup,hes_cup,hcd,qes_cup,q_cup,z_cup,                      &
         dd_massentr,dd_massdetr,gamma_cup,q,he,qco
      integer  ,intent (in   )                             ::           &
         iloop
      integer, dimension (its:ite) ,intent (in   )         ::           &
         jmin
      integer, dimension (its:ite) ,intent (inout)         ::           &
         ierr
      real,    dimension (its:ite,kts:kte) ,intent (out  ) ::           &
         qcd,qrcd,pwd
      real,    dimension (its:ite) ,intent (out  )         ::           &
         pwev,bu
      character*128 :: ierrc(its:ite)
      !
      !  local variables in this routine
      !
      integer   ::     i,k
      real      ::     dh,dz,dq_eva,denom,fix_evap
      !
      bu  =0.  !-- buoyancy
      qcd =0.  !-- in-downdradt water vapor mixing ratio
      qrcd=0.  !-- saturation water vapor mixing ratio
      pwev=0.  !-- column integrated rain evaporation (normalized)
      pwd =0.  !-- rain evaporation at layer k

      if(cumulus == 'shallow') return
      !
      do i=its,itf
         if(ierr(i) /= 0) cycle

         !-- boundary condition in jmin ('level of free sinking')
         k=jmin(i)
         dz=z_cup(i,k+1)-z_cup(i,k)

         qcd(i,k)=q_cup(i,k)

         if(use_wetbulb==1) then
            !--option 1
            !qcd(i,k)=q_wetbulb(i)
            !--option 2
            qcd(i,k)=0.5*(q_wetbulb(i)+qco(i,k)) ! mixture 50% env air + updraft
         endif

         dh=hcd(i,k)-hes_cup(i,k)

         if(dh.lt.0)then
            qrcd(i,k)=(qes_cup(i,k)+(1./xlv)*(gamma_cup(i,k) /(1.+gamma_cup(i,k)))*dh)
         else
            qrcd(i,k)=qes_cup(i,k)
         endif

         pwd (i,k) = zd(i,k)*min(0.,qcd(i,k)-qrcd(i,k))
         qcd (i,k) = qrcd(i,k)
         pwev(i)   = pwev(i)+pwd(i,k)
         bu  (i)   = dz*dh

         do k=jmin(i)-1,kts,-1

            dz=z_cup(i,k+1)-z_cup(i,k)

            !-- downward transport + mixing
            denom = (zd(i,k+1)-0.5*dd_massdetr(i,k)+dd_massentr(i,k))
            if( denom == 0.0 )then
               qcd(i,k)= qcd(i,k+1)
            else
               qcd(i,k)=(qcd(i,k+1)*zd(i,k+1) -0.5*dd_massdetr(i,k)*qcd(i,k+1)+ &
                  dd_massentr(i,k)*q  (i,k)    )/ denom
            endif
            !
            !--- to be negatively buoyant, hcd should be smaller than hes!
            !--- ideally, dh should be negative till dd hits ground, but that is not always
            !--- the case
            !
            dh    = hcd(i,k)-hes_cup(i,k)
            bu(i) = bu(i)+dz*dh
            qrcd(i,k)=qes_cup(i,k)+(1./xlv)*(gamma_cup(i,k) /(1.+gamma_cup(i,k)))*dh

            !-- rain water evaporation amount at layer k
            dq_eva=qcd(i,k)-qrcd(i,k)

            if(dq_eva.gt.0.)then
               dq_eva=0.
               qrcd(i,k)=qcd(i,k)
            endif
            !-- amount of the evaporated rain water
            pwd(i,k)=zd(i,k)*dq_eva  ! kg[water vapor]/kg[air]

            !-- source term for in-downdraft water vapor mixing ratio
            qcd(i,k)=qrcd(i,k)     ! => equiv to qcd = qcd - dq_eva !( -dq_eva >0 => source term for qcd)

            !-- total evaporated rain water
            pwev(i)=pwev(i)+pwd(i,k)

           !-- for GEOS diagnostic
           ! evap(i,k) = - edt * xmb * zd * dq_eva = - edt * xmb * pwd (i,k)
            ! downdfrat temp = (hcd(i,k)-qcd(i,k)*xlv-g*z_cup(i,k))/cp - 273.15

         enddo

         if(pwev(i).ge.0.and.iloop.eq.1)then
            ierr(i)=70
            ierrc(i)="problem with buoy in cup_dd_moisture"
         endif
         if(bu(i).ge.0.and.iloop.eq.1)then
            ierr(i)=73
            ierrc(i)="problem2 with buoy in cup_dd_moisture"
         endif

         if(ZERO_DIFF==0 .and. EVAP_FIX==1) then
            if(abs(pwev(i)) > pwavo(i) .and. ierr(i) == 0)then
               fix_evap = pwavo(i)/(1.e-16+abs(pwev(i)))
               pwev(i)  = 0.

               do k=jmin(i),kts,-1
                  pwd(i,k) = pwd (i,k)*fix_evap
                  pwev(i)  = pwev(i) + pwd(i,k)
                  dq_eva   = pwd (i,k)/(1.e-16+zd(i,k))
                  qcd(i,k) = qrcd(i,k) + dq_eva
               enddo
               if(pwev(i) .ge. 0.)then
                  ierr(i)=70
                  ierrc(i)="problem with buoy in cup_dd_moisture"
               endif
            endif
         endif
      enddo!--- end loop over i

   end subroutine cup_dd_moisture

   !------------------------------------------------------------------------------------

   subroutine cup_env(z,qes,he,hes,t,q,p,z1,psur,ierr,itest,                     &
      itf,ktf,its,ite, kts,kte             )

      implicit none

      integer ,intent (in)                ::                         &
         itf,ktf,its,ite, kts,kte
      !
      ! ierr error value, maybe modified in this routine
      ! q           = environmental mixing ratio
      ! qes         = environmental saturation mixing ratio
      ! t           = environmental temp
      ! tv          = environmental virtual temp
      ! p           = environmental pressure
      ! z           = environmental heights
      ! he          = environmental moist static energy
      ! hes         = environmental saturation moist static energy
      ! psur        = surface pressure
      ! z1          = terrain elevation
      !
      !
      real,    dimension (its:ite,kts:kte)                              &
         ,intent (in   )                   ::                           &
         p,t,q
      real,    dimension (its:ite,kts:kte)                              &
         ,intent (out  )                   ::                           &
         he,hes,qes
      real,    dimension (its:ite,kts:kte)                              &
         ,intent (inout)                   ::                           &
         z
      real,    dimension (its:ite)                                      &
         ,intent (in   )                   ::                           &
         psur,z1
      integer, dimension (its:ite)                                      &
         ,intent (in)                   ::                           &
         ierr
      integer                                                           &
         ,intent (in   )                   ::                           &
         itest
      !
      !  local variables in this routine
      !

      integer                              ::                           &
         i,k,iph
      !      real, dimension (1:2) :: AE,BE,HT
      real, dimension (its:ite,kts:kte) :: tv
      real :: e,tvbar,pqsat
      !      real, external :: satvap
      !      real :: satvap

      he  =0.0
      hes =0.0
      qes =0.0

      if(SATUR_CALC == 0) then
         do k=kts,ktf
            do i=its,itf
               if(ierr(i).eq.0)then

                  e=satvap(t(i,k))
                  qes(i,k)=0.622*e/max(1.e-8,(p(i,k)-e))

                  if(QES(I,K).le.1.e-08  ) QES(I,K)=1.e-08
                  if(QES(I,K).gt.max_qsat) QES(I,K)=max_qsat
                  if(QES(I,K).lt.Q(I,K)  ) QES(I,K)=Q(I,K)
                  !       IF(Q(I,K).GT.QES(I,K))Q(I,K)=QES(I,K)
                  TV(I,K)=T(I,K)+.608*Q(I,K)*T(I,K)
               endif
            enddo
         enddo

      else

         !--- better formulation for the mixed phase regime
         do k=kts,ktf
            do i=its,itf
               if(ierr(i).eq.0)then

                  pqsat=satur_spec_hum(t(i,k),p(i,k))
                  qes(i,k)=pqsat
                  !print*,"qes=",k,p(i,k),1000*qes(i,k),1000*pqsat

                  qes(i,k) = min(max_qsat, max(1.e-08,qes(i,k)))
                  qes(i,k) = max(qes(i,k), q(i,k))
                  tv (i,k) = t(i,k)+.608*q(i,k)*t(i,k)
               endif
            enddo
         enddo
      endif

      !
      !--- z's are calculated with changed h's and q's and t's
      !--- if itest=2
      !
      if(itest.eq.1 .or. itest.eq.0)then
         do i=its,itf
            if(ierr(i).eq.0)then
               Z(I,1)=max(0.,Z1(I))-(ALOG(P(I,1))- &
                  ALOG(PSUR(I)))*287.*TV(I,1)/g
            endif
         enddo

         ! --- calculate heights
         do K=kts+1,ktf
            do i=its,itf
               if(ierr(i).eq.0)then
                  TVBAR=.5*TV(I,K)+.5*TV(I,K-1)
                  Z(I,K)=Z(I,K-1)-(ALOG(P(I,K))- &
                     ALOG(P(I,K-1)))*287.*TVBAR/g
               endif
            enddo
         enddo
      else if(itest.eq.2)then
         do k=kts,ktf
            do i=its,itf
               if(ierr(i).eq.0)then
                  z(i,k)=(he(i,k)-1004.*t(i,k)-2.5e6*q(i,k))/g
                  z(i,k)=max(1.e-3,z(i,k))
               endif
            enddo
         enddo
      else if(itest.eq.-1)then
      endif
      !
      !--- calculate moist static energy - HE
      !    saturated moist static energy - HES
      !
      do k=kts,ktf
         do i=its,itf
            if(ierr(i) /= 0) cycle
            if(itest.le.0) he (i,k)=g*z(i,k)+cp*t(i,k)+xlv*q(i,k)

            hes(i,k)=g*z(i,k)+cp*t(i,k)+xlv*qes(i,k)

            if(he(i,k).ge.hes(i,k)) he(i,k)=hes(i,k)
         enddo
      enddo

   end subroutine cup_env
   !------------------------------------------------------------------------------------
   subroutine cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,he_cup, &
      us, vs,u_cup,v_cup,                                   &
      hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur,tsur,        &
      ierr,z1,itf,ktf,its,ite, kts,kte                      )
      implicit none
      integer ,intent (in   )              :: itf,ktf, its,ite, kts,kte
      !
      ! ierr error value, maybe modified in this routine
      ! q           = environmental mixing ratio
      ! q_cup       = environmental mixing ratio on cloud levels
      ! qes         = environmental saturation mixing ratio
      ! qes_cup     = environmental saturation mixing ratio on cloud levels
      ! t           = environmental temp
      ! t_cup       = environmental temp on cloud levels
      ! p           = environmental pressure
      ! p_cup       = environmental pressure on cloud levels
      ! z           = environmental heights
      ! z_cup       = environmental heights on cloud levels
      ! he          = environmental moist static energy
      ! he_cup      = environmental moist static energy on cloud levels
      ! hes         = environmental saturation moist static energy
      ! hes_cup     = environmental saturation moist static energy on cloud levels
      ! gamma_cup   = gamma on cloud levels
      ! psur        = surface pressure
      ! z1          = terrain elevation
      !
      !
      real,    dimension (its:ite,kts:kte)                              &
         ,intent (in   )                   ::                           &
         qes,q,he,hes,z,p,t,us, vs
      real,    dimension (its:ite,kts:kte)                              &
         ,intent (out  )                   ::                           &
         qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,u_cup,v_cup
      real,    dimension (its:ite)                                      &
         ,intent (in   )                   ::                           &
         psur,z1,tsur
      integer, dimension (its:ite)                                      &
         ,intent (in)                   ::                              &
         ierr
      !
      !  local variables in this routine
      !
      integer                              :: i,k
      real                                 :: p1,p2,ct1,ct2,rho
      integer, save                        :: irun = 0

      qes_cup  =0.
      q_cup    =0.
      hes_cup  =0.
      he_cup   =0.
      z_cup    =0.
      p_cup    =0.
      t_cup    =0.
      gamma_cup=0.
      u_cup    =0.
      v_cup    =0.

      if( clev_grid == 2 ) then
         !--original formulation
         do k=kts+1,ktf
            do i=its,itf
               if(ierr(i) /= 0) cycle
               qes_cup(i,k)=.5*(qes(i,k-1)+qes(i,k))
               q_cup  (i,k)=.5*(q  (i,k-1)+q  (i,k))
               hes_cup(i,k)=.5*(hes(i,k-1)+hes(i,k))
               he_cup (i,k)=.5*(he (i,k-1)+he (i,k))
               if(he_cup(i,k).gt.hes_cup(i,k))he_cup(i,k)=hes_cup(i,k)
               z_cup  (i,k)=.5*(z(i,k-1)+z(i,k))
               p_cup  (i,k)=.5*(p(i,k-1)+p(i,k))
               t_cup  (i,k)=.5*(t(i,k-1)+t(i,k))
               gamma_cup(i,k)=(xlv/cp)*(xlv/(rv*t_cup(i,k) &
                  *t_cup(i,k)))*qes_cup(i,k)
               u_cup  (i,k)=.5*(us(i,k-1)+us(i,k))
               v_cup  (i,k)=.5*(vs(i,k-1)+vs(i,k))

            enddo
         enddo
         do i=its,itf
            if(ierr(i) /= 0) cycle
            qes_cup(i,1)=qes(i,1)
            q_cup(i,1)=q(i,1)
            !hes_cup(i,1)=hes(i,1)
            !he_cup(i,1)=he(i,1)
            hes_cup(i,1)=g*z1(i)+cp*t(i,1)+xlv*qes(i,1)
            he_cup (i,1)=g*z1(i)+cp*t(i,1)+xlv*q  (i,1)
            !z_cup(i,1)=.5*(z(i,1)+z1(i))
            !p_cup(i,1)=.5*(p(i,1)+psur(i))
            z_cup(i,1)=z1(i)
            p_cup(i,1)=psur(i)
            t_cup(i,1)=t(i,1)
            gamma_cup(i,1)=xlv/cp*(xlv/(rv*t_cup(i,1) &
               *t_cup(i,1)))*qes_cup(i,1)
            u_cup(i,1)=us(i,1)
            v_cup(i,1)=vs(i,1)
         enddo

       !do k=kts,ktf
       ! i=1
       !        print*,"air_dens=",k,z_cup(i,k),p_cup(i,k),(p_cup(i,k)-p_cup(i,k+1))/(z_cup(i,k+1)-z_cup(i,k))/g
       !enddo

      elseif( clev_grid == 0) then
         !--- weigthed mean
         do i=its,itf
            if(ierr(i) /= 0) cycle
            p_cup  (i,1)=psur(i)
            z_cup  (i,1)=z1(i)
            do k=kts,ktf-1
               p_cup (i,k+1) = 2.0*p(i,k) - p_cup(i,k)
               z_cup (i,k+1) = 2.0*z(i,k) - z_cup(i,k)
            enddo

            ! ----------- p,T          k+1
            !p1
            ! ----------- p_cup,T_cup  k+1
            !p2
            ! ----------- p,T          k
            !
            ! ----------- p_cup,T_cup  k

            do k=kts,ktf-1
               p1=abs((p    (i,k+1) - p_cup(i,k+1))/(p(i,k+1)-p(i,k)))
               p2=abs((p_cup(i,k+1) - p    (i,k  ))/(p(i,k+1)-p(i,k)))

               t_cup  (i,k+1) = p1*t  (i,k) + p2*t  (i,k+1)

               u_cup  (i,k+1) = p1*us (i,k) + p2*us (i,k+1)
               v_cup  (i,k+1) = p1*vs (i,k) + p2*vs (i,k+1)
               q_cup  (i,k+1) = p1*q  (i,k) + p2*q  (i,k+1)
               he_cup (i,k+1) = p1*he (i,k) + p2*he (i,k+1)

               qes_cup(i,k+1) = p1*qes(i,k) + p2*qes(i,k+1)
               hes_cup(i,k+1) = p1*hes(i,k) + p2*hes(i,k+1)

               if(he_cup(i,k+1).gt.hes_cup(i,k+1))he_cup(i,k+1)=hes_cup(i,k+1)

               gamma_cup(i,k+1)=(xlv/cp)*(xlv/(rv*t_cup(i,k+1) &
                  *t_cup(i,k+1)))*qes_cup(i,k+1)

            enddo
            !--- surface level from X(kts) and X_cup(kts+1) determine X_cup(kts)
            k=kts
            p1=abs(p    (i,k  )-p_cup(i,k))
            p2=abs(p_cup(i,k+1)-p_cup(i,k))

            ct1=(p1+p2)/p2
            ct2=p1/p2

            t_cup  (i,k) = ct1*t  (i,k) - ct2*t_cup(i,k+1)
            q_cup  (i,k) = ct1*q  (i,k) - ct2*q_cup(i,k+1)

            u_cup  (i,k) = ct1*us (i,k) - ct2*u_cup(i,k+1)
            v_cup  (i,k) = ct1*vs (i,k) - ct2*v_cup(i,k+1)
            qes_cup(i,k) = ct1*qes(i,k) - ct2*qes_cup(i,k+1)

            hes_cup(i,k)=g*z_cup(i,k)+cp*t_cup(i,k)+xlv*qes_cup(i,k)
            he_cup (i,k)=g*z_cup(i,k)+cp*t_cup(i,k)+xlv*q_cup  (i,k)

            if(he_cup(i,k).gt.hes_cup(i,k))he_cup(i,k)=hes_cup(i,k)

            gamma_cup(i,k)=xlv/cp*(xlv/(rv*t_cup(i,k)*t_cup(i,k)))*qes_cup(i,k)
         enddo

      elseif(clev_grid == 1) then
         !--- based on Tiedke (1989)
         do i=its,itf
            if(ierr(i) /= 0) cycle
            do k=ktf, kts+1,-1

               qes_cup(i,k) = qes(i,k)
               q_cup  (i,k) = q  (i,k)
               p_cup  (i,k) = 0.5*(p(i,k-1)+p(i,k))
               z_cup  (i,k) = 0.5*(z(i,k-1)+z(i,k))
               t_cup  (i,k) = (max(cp*t(i,k-1)+g*z(i,k-1),cp*t(i,k)+g*z(i,k)) - g*z_cup(i,k))/cp

               if(qes(i,k) < max_qsat) &
                  call get_interp(qes_cup(i,k),t_cup(i,k),p_cup(i,k),qes_cup(i,k),t_cup(i,k))

               q_cup  (i,k) = min(q(i,k),qes(i,k)) + qes_cup(i,k) - qes(i,k)
               q_cup  (i,k) = max(q_cup  (i,k) ,0.0)

            enddo
            !---level kts
            qes_cup(i,1)= qes (i,1)
            q_cup  (i,1)= q   (i,1)
            z_cup  (i,1)= z1  (i)
            p_cup  (i,1)= psur(i)

            t_cup  (i,1)= (cp*t(i,1)+g*z(i,1) - g*z_cup(i,1))/cp

            hes_cup(i,1)=g*z_cup(i,1)+cp*t_cup(i,1)+xlv*qes_cup(i,1)
            he_cup (i,1)=g*z_cup(i,1)+cp*t_cup(i,1)+xlv*q_cup  (i,1)

            gamma_cup(i,1)=xlv/cp*(xlv/(rv*t_cup(i,1)*t_cup(i,1)))*qes_cup(i,1)
            u_cup(i,1)=us(i,1)
            v_cup(i,1)=vs(i,1)

            do k=ktf,kts+1,-1
               p1=max(cp*t_cup(i,k)+g*z_cup(i,k), cp*t_cup(i,k-1)+g*z_cup(i,k-1))
               t_cup(i,k) = (p1-g*z_cup(i,k))/cp

               hes_cup(i,k)=cp*t_cup(i,k)+xlv*qes_cup(i,k)+g*z_cup  (i,k)
               he_cup (i,k)=cp*t_cup(i,k)+xlv*q_cup  (i,k)+g*z_cup  (i,k)
               if(he_cup(i,k).gt.hes_cup(i,k))he_cup(i,k)=hes_cup(i,k)

               gamma_cup(i,k)=(xlv/cp)*(xlv/(rv*t_cup(i,k)*t_cup(i,k)))*qes_cup(i,k)
               u_cup    (i,k)=us(i,k)
               v_cup    (i,k)=vs(i,k)
            enddo
         enddo
      else
         stop "cup_env_clev"
      endif

      return
      !IF( MAPL_AM_I_ROOT() .and. irun == 0) then
      irun = 1
      do i=its,itf
         if(ierr(i) == 0) then
            do k=kts,kte-1
               rho=100*(p_cup(i,k)-p_cup(i,k+1))/(z_cup(i,k+1)-z_cup(i,k))/g ! air dens by hidrostatic balance (kg/m3)
               write(23,101) i,k,z_cup(i,k),p_cup(i,k),t_cup(i,k),q_cup(i,k)*1000.,he_cup(i,k),u_cup(i,k),v_cup(i,k),rho

               rho=100*(p(i,k)-p(i,k+1))/(z(i,k+1)-z(i,k))/g
               write(25,101) i,k,z    (i,k),p    (i,k),t    (i,k),q    (i,k)*1000.,he   (i,k),us   (i,k),vs   (i,k),rho

101            format(2i3,8F15.5)
            enddo
            exit
         !            goto 400
         endif
      enddo
     !ENDIF
   !400 continue

   end subroutine cup_env_clev
   !------------------------------------------------------------------------------------

   subroutine cup_forcing_ens_3d_mid(aa0,aa1,xaa0,mbdt,dtime,ierr,&
      po_cup,ktop,k22,kbcon,kpbl,ichoice, maxens,maxens3, &
      itf,ktf,its,ite, kts,kte, &
      tau_ecmwf,aa1_bl,xf_dicycle, &
      dhdt,xff_mid,zws,hc,hco,he_cup,heo_cup,wlpool,xf_coldpool)

      implicit none
      ! aa0     = cloud work function without forcing effects
      ! aa1     = cloud work function with forcing effects
      ! xaa0    = cloud work function with cloud effects
      ! mbdt    = arbitrary numerical parameter
      ! dtime   = dt over which forcing is applied
      ! kbcon   = LFC of parcel from k22
      ! k22     = updraft originating level
      ! ichoice = flag if only want one closure
      ! name    = deep,mid or shallow convection flag
      !
      integer  ,intent (in   )    :: itf,ktf,its,ite, kts,kte,maxens,maxens3
      integer, dimension (its:ite)          ,intent (in) ::      &
         k22,kbcon,ktop,kpbl
      real,    dimension (its:ite,kts:kte)  ,intent (in) ::      &
         po_cup,dhdt,hc,hco,he_cup,heo_cup
      real,    dimension (its:ite)          ,intent (in) ::      &
         xaa0
      real,    dimension (its:ite)          ,intent (in) ::      &
         aa1,zws,mbdt,  aa0
      real                                  ,intent (in) :: dtime
      integer                               ,intent (in) :: ichoice
      real,    dimension (its:ite)          ,intent (in) :: aa1_bl,tau_ecmwf,wlpool
      integer, dimension (its:ite)          ,intent (inout) :: ierr
      real,    dimension (its:ite)          ,intent (inout) :: xf_dicycle,xf_coldpool
      real,    dimension (its:ite,1:maxens3),intent (out)   :: xff_mid
      !
      !  local variables in this routine
      !
      real,    dimension (1:maxens) :: xk
      integer                       :: i,k
      real                          :: xff_dicycle, trash, blqe,xff_ens1,mf_ens1
      do i=its,itf
         !-initialization
         xff_mid     (i,:)= 0.
         xf_dicycle  (i)  = 0.

         if(ierr(i) /= 0)cycle

         !- Think about this:
         !xff0= (AA1(I)-AA0(I))/DTIME
         !if(xff0.lt.0.) xff_dicycle = 0.

         XK(1)=(XAA0(I)-(AA1(I)))/MBDT(i)

         if(xk(1).le.0.and.xk(1).gt.-0.1*mbdt(i)) xk(1)=-0.1*mbdt(i)
         if(xk(1).gt.0.and.xk(1).lt.1.e-2       ) xk(1)=1.e-2

         !- closure 3 for mid
         if(xk(1) < 0.) xff_mid(i,3)=max(0., -(AA1(i)/tau_ecmwf(i))/xk(1))
      enddo

      do i=its,itf
         if(ierr(i) /= 0) cycle
         !- Boundary layer quasi-equilibrium (Raymond 1995)
         if(k22(i).lt.kpbl(i)+1)then
            blqe=0.
            do k=kts,kbcon(i) !- orig formulation
               !do k=kts,kpbl(i)
               blqe = blqe+100.*dhdt(i,k)*(po_cup(i,k)-po_cup(i,k+1))/g
            enddo
            !trash = max((hc (i,kbcon(i))-he_cup (i,kbcon(i))),1.e1)!- orig formulation
            trash = max((hco(i,kbcon(i))-heo_cup(i,kbcon(i))),1.e1)
            xff_mid(i,2) = max(0.,blqe/trash)
         endif

         !- W* closure (Grant,2001)
         xff_mid(i,1)=0.03*zws(i)
      enddo
      
   end subroutine cup_forcing_ens_3d_mid
   !------------------------------------------------------------------------------------

   subroutine cup_minimi(ARRAY,KS,KEND,KT,ierr,              &
      itf,ktf,its,ite, kts,kte                     )

      implicit none
      !
      !  on input
      !

      ! only local dimensions are need as of now in this routine

      integer                                                           &
         ,intent (in   )                   ::                           &
         itf,ktf,                                    &
         its,ite, kts,kte
      ! array input array
      ! x output array with return values
      ! kt output array of levels
      ! ks,kend  check-range
      real,    dimension (its:ite,kts:kte)                              &
         ,intent (in   )                   ::                           &
         array
      integer, dimension (its:ite)                                      &
         ,intent (in   )                   ::                           &
         ierr,ks,kend
      integer, dimension (its:ite)                                      &
         ,intent (out  )                   ::                           &
         kt
      real,    dimension (its:ite)         ::                           &
         x
      integer                              ::                           &
         i,k,kstop

      do i=its,itf
         KT(I)=KS(I)
         if(ierr(i).eq.0)then
            X(I)=ARRAY(I,KS(I))
            KSTOP=MAX(KS(I)+1,KEND(I))
            !
            do K=KS(I)+1,KSTOP
               if(ARRAY(I,K).lt.X(I)) then
                  X(I)=ARRAY(I,K)
                  KT(I)=K
               endif
            enddo
         endif
      enddo

   end subroutine cup_MINIMI
   !------------------------------------------------------------------------------------
   subroutine cup_up_aa0(aa0,z_cup,zu,dby,GAMMA_CUP,t_cup,   &
      k22,klcl,kbcon,ktop,ierr,                      &
      itf,ktf,its,ite, kts,kte,integ_interval        )
      ! aa0 cloud work function
      ! gamma_cup = gamma on model cloud levels
      ! t_cup = temperature (Kelvin) on model cloud levels
      ! dby = buoancy term
      ! zu= normalized updraft mass flux
      ! z = heights of model levels
      ! ierr error value, maybe modified in this routine
      !
      implicit none
      !
      ! on input
      integer,intent (in   )              ::  itf,ktf,its,ite, kts,kte
      real,    dimension (its:ite,kts:kte) ,intent (in   )  ::  &
         z_cup,zu,gamma_cup,t_cup,dby
      integer, dimension (its:ite)         ,intent (in   )  ::  &
         k22,klcl,kbcon,ktop
      character(len=*),optional, intent(in) :: integ_interval
      ! input and output
      integer, dimension (its:ite) ,intent (in   )  :: ierr
      real,    dimension (its:ite) ,intent (out  )  :: aa0
      !
      !  local variables in this routine
      integer                             ::  i,k
      real                                ::  dz,da,aa_2,aa_1
      integer, dimension (its:ite)        ::  kbeg,kend
      !
      !
      !  initialize array to zero.
      aa0(:)=0.
      !  set domain of integration
      if(present(integ_interval)) then
         if(integ_interval == 'BL') then
            kbeg(:) = kts
            kend(:) = kbcon(:)-1
         elseif(integ_interval == 'CIN') then
            kbeg(:) = KTS  ! k22(:) !klcl (:) ! kts
            kend(:) = kbcon(:) ! kbcon(:)-1
         else
            stop "unknown range in cup_up_aa0"
         endif
      else
         kbeg(:) = kbcon(:)
         kend(:) = ktop (:)
      endif

      do i=its,itf
         if(ierr(i) /= 0)cycle
         do k= kbeg(i),kend(i)

            dz=z_cup(i,k+1)-z_cup(i,k)
            aa_1=zu(i,k  )*(g/(cp*t_cup(i,k  )))*dby(i,k  )/(1.+gamma_cup(i,k  ))
            aa_2=zu(i,k+1)*(g/(cp*t_cup(i,k+1)))*dby(i,k+1)/(1.+gamma_cup(i,k+1))
            da=0.5*(aa_1+aa_2)*dz

            aa0(i)=aa0(i)+da

           !aa0(i)=aa0(i)+max(0.,da)
         enddo
      enddo

   end subroutine cup_up_aa0
   !------------------------------------------------------------------------------------

   subroutine cup_up_moisture(name,start_level,klcl,ierr,ierrc,z_cup,qc,qrc,pw,pwav,hc,tempc,xland,&
      po,p_cup,kbcon,ktop,cd,dby,clw_all,                   &
      t_cup,q,gamma_cup,zu,qes_cup,k22,qe_cup,           &
      zqexec,use_excess,ccn,rho,                            &
      up_massentr,up_massdetr,psum,psumh,c1d,x_add_buoy,  &
      vvel2d,vvel1d,zws,entr_rate_2d,                         &
      itest,itf,ktf,ipr,jpr,its,ite, kts,kte       )

      implicit none
      real, parameter :: bdispm = 0.366       !berry--size dispersion (maritime)
      real, parameter :: bdispc = 0.146       !berry--size dispersion (continental)
      real, parameter :: T_BF   = 268.16 , T_ice_BF = 235.16
      real, parameter :: rk = 3 ! or 2
      real, parameter :: xexp = 2.
      !
      !  on input
      integer  ,intent (in   ) ::  use_excess,itest,itf,ktf    &
         ,its,ite,ipr,jpr, kts,kte
      ! cd= detrainment function
      ! q = environmental q on model levels
      ! qe_cup = environmental q on model cloud levels
      ! qes_cup = saturation q on model cloud levels
      ! dby = buoancy term
      ! cd= detrainment function
      ! zu = normalized updraft mass flux
      ! gamma_cup = gamma on model cloud levels
      !
      character *(*)                    ,intent (in) ::  name
      integer, dimension (its:ite)      ,intent (in) ::  kbcon,ktop,k22,klcl,start_level
      real,  dimension (its:ite,kts:kte),intent (in) ::  t_cup,p_cup,rho,q,zu,gamma_cup       &
         ,qe_cup,hc,po,up_massentr,up_massdetr &
         ,dby,qes_cup,z_cup,cd,c1d

      real,  dimension (its:ite)        ,intent (in) ::  zqexec,xland,x_add_buoy
      real,  dimension (its:ite)        ,intent (in) ::  zws,ccn
      real,  dimension (its:ite,kts:kte),intent (in) ::  entr_rate_2d
      real,  dimension (its:ite,kts:kte),intent (in) ::  vvel2d
      real,  dimension (its:ite        ),intent (in) ::  vvel1d
      !
      ! input and output
      !
      ! ierr error value, maybe modified in this routine
      integer, dimension (its:ite)  ,intent (inout)                   ::  ierr
      ! qc = cloud q (including liquid water) after entrainment
      ! qrch = saturation q in cloud
      ! qrc = liquid water content in cloud after rainout
      ! pw = condensate that will fall out at that level
      ! pwav = totan normalized integrated condensate (I1)
      ! c0 = conversion rate (cloud to rain)

      real,   dimension (its:ite,kts:kte),intent (out)   :: qc,qrc,pw,clw_all,tempc
      real,   dimension (its:ite)        ,intent (out)   :: pwav,psum,psumh
      character*128                      ,intent (inout) :: ierrc(its:ite)
      !
      !  local variables in this routine
      !
      integer                              ::                           &
         iounit,iprop,i,k,k1,k2,n,nsteps
      real                                 ::                           &
         dp,rhoc,dh,qrch,dz,radius,berryc0,q1,berryc
      real :: qaver,denom,aux,cx0,qrci,step,cbf,qrc_crit_BF,min_liq,qavail
      real delt,tem1,qrc_0,cup
      !real,   dimension (its:ite,kts:kte) :: qc2

      !--- no precip for small clouds
      !if(name.eq.'shallow')  c0 = c0_shal
      !if(name.eq.'mid'    )  c0 = c0_mid
      !if(name.eq.'deep'   )  c0 = c0_deep
      do i=its,itf
         pwav (i)=0.
         psum (i)=0.
         psumh(i)=0.
      enddo
      do k=kts,ktf
         do i=its,itf
            pw      (i,k)=0.
            clw_all (i,k)=0.
            tempc   (i,k)=t_cup (i,k)
            qrc     (i,k)=0.          !--- liq/ice water
            qc      (i,k)=qe_cup(i,k) !--- total water: liq/ice = vapor water
         !qc2     (i,k)=qe_cup(i,k) !--- total water: liq/ice = vapor water
         enddo
      enddo

      !--- get boundary condition for qc
      do i=its,itf
         if(ierr(i)  /= 0) cycle
         call get_cloud_bc(name,kts,kte,ktf,xland(i),po(i,kts:kte),qe_cup (i,kts:kte),qaver,k22(i))
         qc  (i,kts:start_level(i)) = qaver + zqexec(i) +     1.* x_add_buoy(i)/xlv
        !qc  (i,kts:start_level(i)) = qaver + zqexec(i) +     0.67* x_add_buoy(i)/xlv
         qrc (i,kts:start_level(i)) = 0.
         !qc  (i,kts:start_level(i)) = qaver + zqexec(i) + 0.5*x_add_buoy(i)/xlv
         !qc2 (i,kts:start_level(i)) = qaver + zqexec(i) + 0.5*x_add_buoy(i)/xlv
      enddo

      !--- option to produce linear fluxes in the sub-cloud layer.
      if(name == 'shallow' .and. use_linear_subcl_mf == 1) then
         do i=its,itf
            if(ierr(i) /= 0) cycle
            call get_delmix(name,kts,kte,ktf,xland(i),start_level(i),po(i,kts:kte) &
               ,qe_cup(i,kts:kte), qc(i,kts:kte))
         enddo
      endif
      do i=its,itf
         if (ierr(i) /= 0) cycle

         do k=start_level(i) + 1,ktop(i) + 1

            DZ=Z_cup(i,K)-Z_cup(i,K-1)
            !
            !--- saturation  in cloud, this is what is allowed to be in it
            !
            QRCH = qes_cup(I,K)+(1./XLV)*(GAMMA_cup(i,k)/(1.+GAMMA_cup(i,k)))*DBY(I,K)

            !-    1. steady state plume equation, for what could
            !-       be in cloud without condensation
            denom =  (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
            if(denom > 0.) then

               qc (i,k)=  ( qc (i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)* qc(i,k-1) +   &
                  up_massentr(i,k-1)* q (i,k-1)     &
                  )/ denom

               if(k==start_level(i)+1) qc(i,k)= qc(i,k) + (zqexec(i) + 0.5*x_add_buoy(i)/xlv) &
                                              * up_massentr(i,k-1)/denom
                                              !--- assuming no liq/ice water in the environment
                                              qrc(i,k)=  ( qrc(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)* qrc(i,k-1)   &
                                              )/ denom

            else
               qc (i,k)=    qc (i,k-1)
               qrc(i,k)=    qrc(i,k-1)
            endif

            !            qc2(i,k)= ( (1.-0.5*entr_rate_2d(i,k-1)*dz)*qc2(i,k-1)     &
            !                              + entr_rate_2d(i,k-1)*dz *q  (i,k-1) ) / &
            !                        (1.+0.5*entr_rate_2d(i,k-1)*dz)

            !-- updraft temp
            tempc(i,k) = (1./cp)*(hc(i,k)-g*z_cup(i,k)-xlv*QRCH)

            !--- total condensed water before rainout
            clw_all(i,k)= max(0.,qc(i,k)-qrch)

            qrc   (i,k) = min(clw_all(i,k),qrc(i,k))

            !--- production term => condensation/diffusional growth
            cup         = max(0.,qc(i,k)-qrch-qrc(i,k))/dz

            if(c0 < 1.e-6)then
               qrc (i,k) = clw_all(i,k)
               qc  (i,k) = qrc(i,k)+min(qc(i,k),qrch)
               pwav(i)   = 0.
               psum(i)   = psum(i)+clw_all(i,k)*zu(i,k) *dz
               cycle
            endif

            if (autoconv == 1 ) then
               min_liq  = qrc_crit * ( xland(i)*1. + (1.-xland(i))*0.7 )
               if(name.eq.'mid') min_liq = min_liq*0.5 
               
               cx0     = (c1d(i,k)+c0)*DZ
               qrc(i,k)= clw_all(i,k)/(1.+cx0)
               pw (i,k)= cx0*max(0.,qrc(i,k) - min_liq)! units kg[rain]/kg[air]
               !pw (i,k)= cx0*qrc(i,k)    ! units kg[rain]/kg[air]

               !--- convert pw to normalized pw
               pw (i,k)=pw(i,k)*zu(i,k)

            elseif (autoconv == 5 ) then
               !  c0_deep     = 1.5e-3; c0_mid     = 1.5e-3 ; qrc_crit        = 1.e-4 !(kg/kg)

               min_liq  = qrc_crit * ( xland(i)*0.4 + (1.-xland(i))*1. )

               if(clw_all(i,k) <= min_liq) then !=> more heating at upper levels, more detrained ice

                  qrc(i,k)= clw_all(i,k)
                  pw(i,k) = 0.
               else

                  cx0     = (c1d(i,k)+c0)*(1.+ 0.33*fract_liq_f(tempc(i,k)))
                  !cx0     = (c1d(i,k)+c0)*(1.+ 2.*fract_liq_f(tempc(i,k)))
                  !--- v0
                  qrc(i,k)= qrc(i,k)*exp(-cx0*dz) + (cup/cx0)*(1.-exp(-cx0*dz))
                  qrc(i,k)= max(qrc(i,k),min_liq)
                  pw (i,k)= max(0.,clw_all(i,k)-qrc(i,k)) ! units kg[rain]/kg[air]
                  qrc(i,k)= clw_all(i,k)-pw(i,k)
                  !--- v1
                  !  qrc_0   = qrc(i,k)
                  !  qrc(i,k)= (qrc_0-min_liq)*exp(-cx0*dz) + (cup/cx0)*(1.-exp(-cx0*dz))+min_liq
                  !  qrc(i,k)= max(qrc(i,k),min_liq)
                  !  pw (i,k)= max(0.,clw_all(i,k)-qrc(i,k)) ! units kg[rain]/kg[air]
                  !  qrc(i,k)= clw_all(i,k)-pw(i,k)

                  !  qrc(i,k)= (clw_all(i,k)-min_liq)*exp(-cx0*dz)+min_liq
                  !  pw (i,k)= clw_all(i,k)-qrc(i,k) ! units kg[rain]/kg[air]
                  !--- v3
                  !  qrc(i,k)= (clw_all(i,k)-min_liq) / (1.+cx0*dz)+min_liq
                  !  pw (i,k)= cx0*dz*(qrc(i,k)-min_liq) ! units kg[rain]/kg[air]
                  !  print*,"BG=",k,real(cx0*1.e+3,4),real(pw(i,k),4),real(qrc(i,k),4)&
                  !              ,real(clw_all(i,k)-pw(i,k)-qrc(i,k),4) !==> must be zero

                  !--- convert pw to normalized pw
                  pw (i,k)= pw(i,k)*zu(i,k)
               endif

            elseif (autoconv == 6 ) then
               min_liq  = 0.5* qrc_crit * (xland(i)*1.5+(1.-xland(i))*2.5)

               if(clw_all(i,k) <= min_liq) then
                  qrc(i,k)= clw_all(i,k)
                  pw(i,k) = 0.
               else
                  cx0     = (c1d(i,k)+c0)*dz
                  qrc(i,k)= (clw_all(i,k))*exp(-cx0)
                  pw (i,k)= clw_all(i,k) - qrc(i,k)
                  pw (i,k)= pw(i,k)*zu(i,k)
               endif
                !
                !print*,"6mass=",pw(i,k)/(1.e-12+zu(i,k))+qrc(i,k),clw_all(i,k)
            elseif (autoconv == 7 ) then
               min_liq  = 0.5* qrc_crit * (xland(i)*1.5+(1.-xland(i))*2.5)

               if(clw_all(i,k) <= min_liq) then
                  qrc(i,k)= clw_all(i,k)
                  pw(i,k) = 0.
               else
                  cx0     = c1d(i,k)+c0
                  qrc_0   = qrc(i,k)
                  qrc(i,k)= qrc_0*exp(-cx0*dz) + (cup/cx0)*(1.-exp(-cx0*dz))

                  pw (i,k)= max(clw_all(i,k) - qrc(i,k),0.)
                  qrc(i,k)= clw_all(i,k) - pw (i,k)
                  pw (i,k)= pw(i,k)*zu(i,k)
               endif
                !
                !print*,"6mass=",pw(i,k)/(1.e-12+zu(i,k))+qrc(i,k),clw_all(i,k)

            elseif (autoconv == 3 ) then
               min_liq  =qrc_crit ! * (xland(i)*1.5+(1.-xland(i))*2.5)

               if(clw_all(i,k) <= min_liq) then
                  qrc(i,k)= clw_all(i,k)
                  pw(i,k) = 0.
               else
                  DELT=-5.
                  if(t_cup(i,k) > 273.16 + DELT) then
                     aux = 1.
                  else
                     aux = 1. * exp(0.07* (t_cup(i,k) - (273.16 + DELT)))
                  endif
                  cx0     = aux*c0
                  !                      cx0     = max(cx0,c0)
                  !                      cx0     = max(cx0,0.25*c0)
                  cx0     = max(cx0,0.50*c0)
                  qrc_0   = qrc(i,k)
                  qrc(i,k)= qrc_0*exp(-cx0*dz) + (cup/cx0)*(1.-exp(-cx0*dz))
                  qrc(i,k)= min(clw_all(i,k), qrc(i,k))
                  pw (i,k)= clw_all(i,k) - qrc(i,k)
                  pw (i,k)= pw(i,k)*zu(i,k)
                 !if(pw(i,k)<0.) stop " pw<0 autoc 3"
               endif

            elseif (autoconv == 4 ) then

               min_liq  = (xland(i)*0.3+(1.-xland(i))*0.5)*1.e-3

               if(clw_all(i,k) > min_liq) then

                  tem1 = fract_liq_f(tempc(i,k))
                  cbf  = 1.
                  if(tempc(i,k) < T_BF) cbf=1.+0.5*sqrt(min(max(T_BF-tempc(i,k),0.),T_BF-T_ice_BF))
                  !qrc_crit_BF = 0.5e-3/cbf
                  qrc_crit_BF = 3.e-4/cbf
                  cx0 = c0*cbf*(tem1*1.3+(1.-tem1))/(0.75*min(15.,max(vvel2d(i,k),1.)))

                  !---solution 1 by Runge-Kutta
                  !step = cx0*dz
                  !do n=int(rk),1,-1
                  !  aux     = qrc(i,k)/qrc_crit_BF
                  !  pw (i,k)= auto_rk(n,step,aux,xexp,qrc(i,k))
                  !  qrc(i,k)= max(clw_all(i,k) - pw(i,k), min_liq)
                  !enddo
                  !---

                  !---solution 2 by Runge-Kutta
                  !qrc_0 = qrc(i,k)
                  !step  = cx0*dz
                  !do n = int(rk),1,-1
                  !aux      = qrc(i,k)/qrc_crit_BF
                  !pw (i,k) =-step*qrc(i,k)*(1.0-exp(-aux**xexp))/float(n) + cup*dz/float(n)
                  !pw (i,k) = max(-qrc_0, pw(i,k))
                  !qrc(i,k) = qrc_0 + pw(i,k)
                  !enddo
                  !---

                  !---analytical solution
                  qrc_0   = qrc(i,k)
                  cx0     = cx0* (1.- exp(- (qrc_0/qrc_crit_BF)**2))
                  qrc(i,k)= qrc_0*exp(-cx0*dz) + (cup/cx0)*(1.-exp(-cx0*dz))

                  pw (i,k)= max(clw_all(i,k) - qrc(i,k),0.)
                    !--- convert PW to normalized PW
                  pw (i,k) = pw(i,k)*zu(i,k)

                   !if(pw(i,k)<-1.e-12) stop " pw<0 autoc 4"
               else
                  pw (i,k) = 0.0
                  qrc(i,k) = clw_all(i,k)
               endif
            endif
            !- total water (vapor + condensed) in updraft after the rainout
            qc(i,k)=qrc(i,k)+min(qc(i,k),qrch)

            !--- integrated normalized condensates
            pwav(i)=pwav(i)+pw(i,k)
            psum(i)=psum(i)+clw_all(i,k)*zu(i,k) *dz

         enddo

         if(ZERO_DIFF==0) then
            if(pwav(i) < 0.) then
               ierr(i)=66
               ierrc(i)="pwav negative"
            endif
         endif

      enddo

      !--- get back water vapor qc
      do i=its,itf
         if (ierr(i)  /= 0) cycle
         do k=kts,ktop(i)+1
            qc(i,k)= qc(i,k)-qrc(i,k)
           !if(qc(i,k) < 0.)stop " qc negative"
         enddo
      enddo

   end subroutine cup_up_moisture

   !------------------------------------------------------------------------------------
   subroutine cup_up_moisture_light(name,start_level,klcl,ierr,ierrc,z_cup,qc,qrc,pw,pwav,hc,tempc,xland &
      ,po,p_cup,kbcon,ktop,cd,dby,clw_all,t_cup,q,gamma_cup,zu  &
      ,qes_cup,k22,qe_cup,zqexec,use_excess,rho                 &
      ,up_massentr,up_massdetr,psum,psumh,c1d,x_add_buoy        &
      ,itest,itf,ktf,ipr,jpr,its,ite, kts,kte                   )

      implicit none
      !
      !  on input
      integer  ,intent (in   ) ::  use_excess,itest,itf,ktf    &
         ,its,ite,ipr,jpr, kts,kte
      ! cd= detrainment function
      ! q = environmental q on model levels
      ! qe_cup = environmental q on model cloud levels
      ! qes_cup = saturation q on model cloud levels
      ! dby = buoancy term
      ! cd= detrainment function
      ! zu = normalized updraft mass flux
      ! gamma_cup = gamma on model cloud levels
      !
      character *(*)                    ,intent (in) ::  name
      integer, dimension (its:ite)      ,intent (in) ::  kbcon,ktop,k22,klcl,start_level
      real,  dimension (its:ite,kts:kte),intent (in) ::  t_cup,p_cup,rho,q,zu,gamma_cup       &
         ,qe_cup,hc,po,up_massentr,up_massdetr &
         ,dby,qes_cup,z_cup,cd,c1d

      real,  dimension (its:ite)        ,intent (in) ::  zqexec,xland,x_add_buoy
      !
      ! input and output
      !
      ! ierr error value, maybe modified in this routine
      integer, dimension (its:ite)  ,intent (inout)                   ::  ierr
      ! qc = cloud q (including liquid water) after entrainment
      ! qrch = saturation q in cloud
      ! qrc = liquid water content in cloud after rainout
      ! pw = condensate that will fall out at that level
      ! pwav = totan normalized integrated condensate (I1)
      ! c0 = conversion rate (cloud to rain)

      real,   dimension (its:ite,kts:kte),intent (out)   :: qc,qrc,pw,clw_all,tempc
      real,   dimension (its:ite)        ,intent (out)   :: pwav,psum,psumh
      character*128                      ,intent (inout) :: ierrc(its:ite)
      !
      !  local variables in this routine
      !
      integer                              ::                           &
         iounit,iprop,i,k,k1,k2,n,nsteps
      real                                 ::                           &
         dp,rhoc,dh,qrch,dz,radius,berryc0,q1,berryc
      real :: qaver,denom,aux,cx0,qrci,step,cbf,qrc_crit_BF,min_liq,qavail,delt_hc_glac
      real delt,tem1

      !!--- no precip for small clouds
      !if(name.eq.'shallow')  c0 = c0_shal
      !if(name.eq.'mid'    )  c0 = c0_mid
      !if(name.eq.'deep'   )  c0 = c0_deep
      do i=its,itf
         pwav (i)=0.
         psum (i)=0.
         psumh(i)=0.
      enddo
      do k=kts,ktf
         do i=its,itf
            pw      (i,k)=0.
            qrc     (i,k)=0.
            clw_all (i,k)=0.
            tempc   (i,k)=t_cup (i,k)
            qc      (i,k)=qe_cup(i,k)
         enddo
      enddo

      !--- get boundary condition for qc
      do i=its,itf
         if(ierr(i)  /= 0) cycle
         call get_cloud_bc(name,kts,kte,ktf,xland(i),po(i,kts:kte),qe_cup (i,kts:kte),qaver,k22(i))
         qc (i,kts:start_level(i)) = qaver + zqexec(i) + 0.5*x_add_buoy(i)/xlv
      enddo

      do i=its,itf
         if (ierr(i)  /= 0) cycle

         do k=start_level(i)+1,ktop(i) + 1

            DZ=Z_cup(i,K)-Z_cup(i,K-1)
            !
            !--- saturation  in cloud, this is what is allowed to be in it
            !
            QRCH = qes_cup(I,K)+(1./XLV)*(GAMMA_cup(i,k)/(1.+GAMMA_cup(i,k)))*DBY(I,K)

            !-    1. steady state plume equation, for what could
            !-       be in cloud without condensation
            denom =  (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
            if(denom > 0.) then
               qc (i,k)=  ( qc (i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)* qc(i,k-1) +   &
                  up_massentr(i,k-1)* q (i,k-1)     &
                  )/ denom
               if(k==start_level(i)+1) qc(i,k)= qc(i,k) +  (zqexec(i) + 0.5*x_add_buoy(i)/xlv)&
                                               * up_massentr(i,k-1)/denom
            else
               qc (i,k)=    qc (i,k-1)
            endif

            !--- total condensed water before rainout
            clw_all(i,k)=max(0.,QC(I,K)-QRCH)
            !--- updraft temp
            tempc(i,k) = (1./cp)*(hc(i,k) - g*z_cup(i,k)-xlv*QRCH)

            !--add glaciation effect on the MSE
            if(MELT_GLAC) then
               delt_hc_glac = clw_all(i,k)*(1.- fract_liq_f(tempc(i,k)))*xlf

               tempc(i,k) = tempc(i,k)+(1./cp)*delt_hc_glac
            endif

            cx0     = (c1d(i,k)+c0)*DZ
            if(c0 < 1.e-6) cx0 = 0.

            qrc(i,k)= clw_all(i,k)/(1.+cx0)
            pw (i,k)= cx0*max(0.,qrc(i,k) -qrc_crit)! units kg[rain]/kg[air]
            !--- convert pw to normalized pw
            pw (i,k)= pw(i,k)*zu(i,k)

            !- total water (vapor + condensed) in updraft after the rainout
            qc(i,k)=qrc(i,k)+min(qc(i,k),qrch)

         enddo
      enddo

      !- get back water vapor qc
      do i=its,itf
         if (ierr(i)  /= 0) cycle
         do k=kts,ktop(i)+1
            qc(i,k)= qc(i,k)-qrc(i,k)
         enddo
      enddo

   end subroutine cup_up_moisture_light

   !------------------------------------------------------------------------------------

   real function satvap(temp2)
      implicit none
      real :: temp2, temp, toot, toto, eilog, tsot,  &
         ewlog, ewlog2, ewlog3, ewlog4
      temp = temp2-273.155
      if (temp.lt.-20.) then   !!!! ice saturation
         toot = 273.16 / temp2
         toto = 1 / toot
         eilog = -9.09718 * (toot - 1) - 3.56654 * (log(toot) / &
            log(10.)) + .876793 * (1 - toto) + (log(6.1071) / log(10.))
         satvap = 10 ** eilog
      else
         tsot = 373.16 / temp2
         ewlog = -7.90298 * (tsot - 1) + 5.02808 * &
            (log(tsot) / log(10.))
         ewlog2 = ewlog - 1.3816e-07 * &
            (10 ** (11.344 * (1 - (1 / tsot))) - 1)
         ewlog3 = ewlog2 + .0081328 * &
            (10 ** (-3.49149 * (tsot - 1)) - 1)
         ewlog4 = ewlog3 + (log(1013.246) / log(10.))
         satvap = 10 ** ewlog4
      endif

   end function
   !------------------------------------------------------------------------------------

   subroutine cup_up_aa1bl(version,aa1_bl,aa1_fa,aa1,t,tn,q,qo,dtime,po_cup,  &
      z_cup,zu,dby,GAMMA_CUP,t_cup,rho,                  &
      klcl,kpbl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte,&
      xland,ztexec,xlons,xlats, h_sfc_flux,le_sfc_flux,  &
      tau_bl,tau_ecmwf,t_star,cumulus ,tn_bl,qo_bl       )

      implicit none
      character*(*), intent(in)                         :: cumulus
      integer      , intent (in   )                     :: &
         itf,ktf,its,ite, kts,kte,version
      ! aa0 cloud work function
      ! gamma_cup = gamma on model cloud levels
      ! t_cup = temperature (Kelvin) on model cloud levels
      ! dby = buoancy term
      ! zu= normalized updraft mass flux
      ! z = heights of model levels
      ! ierr error value, maybe modified in this routine
      !
      real,    dimension (its:ite,kts:kte) ,intent (in  ) ::    &
         z_cup,zu,gamma_cup,t_cup,dby,t,tn,q,qo,po_cup,rho,tn_bl,qo_bl
      integer, dimension (its:ite)         ,intent (in  ) ::    &
         klcl,kbcon,ktop,kpbl
      real                                 ,intent(in   ) :: dtime, t_star

      real,    dimension (its:ite)         ,intent (in  ) ::    &
         xland,ztexec,xlons,xlats, h_sfc_flux,le_sfc_flux       &
         ,aa1,tau_bl,tau_ecmwf
      !
      ! input and output
      integer, dimension (its:ite),intent (inout) ::  ierr
      real,    dimension (its:ite),intent (out  ) ::  aa1_bl,aa1_fa
      !  local variables in this routine
      !
      integer                              ::    i,k,iprloc
      real                                 ::    dz,da,aa_1,aa_2,tcup,da_bl,a1_bl

      !
      !
      aa1_bl (:)=0.
      if(version == 0 ) then
         do i=its,itf
            if(ierr(i) /= 0 ) cycle
            !***       do k=kts,kbcon(i)
            do k=kts,kpbl(i)
               dz = g * (z_cup (i,k+1)-z_cup (i,k))
               da = dz*(tn(i,k)*(1.+0.608*qo(i,k))-t(i,k)*(1.+0.608*q(i,k)))/dtime
               aa1_bl(i)=aa1_bl(i)+da ! Units : J K / (kg seg)
            enddo
         enddo
      elseif(version==1) then
         do i=its,itf
            if(ierr(i) /= 0 ) cycle
            do k=kts,kpbl(i)
               dz = (z_cup (i,k+1)-z_cup (i,k))
               aa_1=(g/(cp*t_cup(i,k  )))*dby(i,k  )*zu(i,k  )
               aa_2=(g/(cp*t_cup(i,k+1)))*dby(i,k+1)*zu(i,k+1)
               da=0.5*(aa_1+aa_2)*dz! Units : J / kg
               aa1_bl(i)=aa1_bl(i)+da
            enddo
         enddo
      else
         stop "unknown version option in routine: cup_up_aa1bl"
      endif

      return
      
      aa1_fa (:)=0.
      do i=its,itf
         if(ierr(i) /= 0)cycle
         do k= kbcon(i),ktop(i)

            dz=z_cup(i,k+1)-z_cup(i,k)
            aa_1=(g/(cp*((t_cup(i,k  )))))*dby(i,k  )/(1.+gamma_cup(i,k  ))*zu(i,k)
            aa_2=(g/(cp*((t_cup(i,k+1)))))*dby(i,k+1)/(1.+gamma_cup(i,k+1))*zu(i,k+1)
            da=0.5*(aa_1+aa_2)*dz

            aa1_fa(i)=aa1_fa(i)+da
         enddo
      enddo

   end subroutine cup_up_aa1bl
   !------------------------------------------------------------------------------------
   subroutine get_lateral_massflux(itf,ktf, its,ite, kts,kte,min_entr_rate               &
      ,ierr,ktop,zo_cup,zuo,cd,entr_rate_2d,po_cup          &
      ,up_massentro, up_massdetro ,up_massentr, up_massdetr &
      ,draft,kbcon,k22,kpbl,up_massentru,up_massdetru,lambau)
      implicit none
      character *(*), intent (in) :: draft
      integer, intent(in) :: itf,ktf, its,ite, kts,kte
      real,    intent(in) :: min_entr_rate
      integer, intent(in)   , dimension(its:ite)            :: ierr,ktop,kbcon,k22,kpbl
      real,    intent(in)   , dimension(its:ite), optional  :: lambau
      real,    intent(in)   , dimension(its:ite,kts:kte) :: zo_cup,zuo,po_cup
      real,    intent(inout), dimension(its:ite,kts:kte) :: cd,entr_rate_2d
      real,    intent(  out), dimension(its:ite,kts:kte) :: up_massentro, up_massdetro&
         ,up_massentr,  up_massdetr
      real,    intent(  out), dimension(its:ite,kts:kte),  &
         optional :: up_massentru,up_massdetru
      !
      !-- local vars
      integer :: i,k, turn, ismooth1,ismooth2
      real :: dz, mass1,mass2,dp,rho,zuo_ave
      logical  :: SMOOTH
      integer, parameter :: MASS_U_OPTION = 1
      integer, parameter :: SMOOTH_DEPTH  = 2 ! -- increasing this parameter,
                                              ! -- strongly damps the heat/drying rates, precip ...
      integer ::  incr1=1, incr2=1, nlay, k_ent
      !---
      !---non-zero-diff-APR-08-2020
      SMOOTH = .false.
      if(USE_SMOOTH_PROF == 1)  SMOOTH = .true.
      if(ZERO_DIFF       == 1)  SMOOTH = .false.
      !---non-zero-diff-APR-08-2020

      up_massentro(:,:)=0.
      up_massdetro(:,:)=0.
      if(present(up_massentru) .and. present(up_massdetru))then
         up_massentru(:,:)=0.
         up_massdetru(:,:)=0.
      endif
      nlay = int(kte/90)

      do i=its,itf
         if(ierr(i)/= 0)cycle

         !-will not allow detrainment below the location of the maximum zu
         ! if(draft=='shallow'.or.draft == 'mid') cd(i,1:maxloc(zuo(i,:),1)-2)=0.0

          !-will not allow detrainment below cloud base or in the PBL
         if(draft=='shallow') then
            cd(i,1:max(kbcon(i),kpbl(i))+nlay)=0.0

         else
            cd(i,1:maxloc(zuo(i,:),1)+nlay)=0.0
         endif

         !- mass entrainment and detrainment are defined on model levels
         do k=kts,maxloc(zuo(i,:),1)
            !=> below location of maximum value zu -> change entrainment

            dz=zo_cup(i,k+1)-zo_cup(i,k)
            zuo_ave = 0.5*(zuo(i,k+1)+zuo(i,k))

            up_massdetro(i,k)=cd(i,k)*dz*zuo_ave

            up_massentro(i,k)=zuo(i,k+1)-zuo(i,k)+up_massdetro(i,k)
            up_massentro(i,k)=max(up_massentro(i,k),min_entr_rate*dz*zuo_ave)

            !-- check dd_massdetro in case of dd_massentro has been changed above
            up_massdetro(i,k)=-zuo(i,k+1)+zuo(i,k)+up_massentro(i,k)

            !if(zuo(i,k-1).gt.0.) then
            cd          (i,k)=up_massdetro(i,k)/(dz*zuo_ave)
            entr_rate_2d(i,k)=up_massentro(i,k)/(dz*zuo_ave)
           !endif
           !if(draft=='shallow')print*,"ent1=",k,real(entr_rate_2d(i,k),4)!,real((min(zo_cup(i,k_ent)/zo_cup(i,k-1),1.)))

         enddo

         !--- limit the effective entrainment rate
         k_ent=maxloc(zuo(i,:),1)
         do k=k_ent+1,ktop(i)-1
            entr_rate_2d(i,k)=entr_rate_2d(i,k_ent)*(min(zo_cup(i,k_ent)/zo_cup(i,k),1.))
            entr_rate_2d(i,k)=max(min_entr_rate, entr_rate_2d(i,k))
          !if(draft=='shallow')print*,"ent2=",k,real(entr_rate_2d(i,k),4),real((min(zo_cup(i,k_ent)/zo_cup(i,k),1.)))
         enddo
         entr_rate_2d(i,ktop(i):kte)=0.

         !=================
         if(SMOOTH .and. draft /= 'shallow' ) then
            !---smoothing the transition zone (from maxloc(zu)-1 to maxloc(zu)+1)

            ismooth1 = max(kts+2, maxloc(zuo(i,:),1) - SMOOTH_DEPTH)
            ismooth2 = min(ktf-2, maxloc(zuo(i,:),1) + SMOOTH_DEPTH)
            !if(draft == 'shallow') ismooth1 = max(ismooth1,max(kbcon(i),kpbl(i))+nlay)+1

            do k=ismooth1,ismooth2
               dz=zo_cup(i,k+1)-zo_cup(i,k)

               zuo_ave = 0.5*(zuo(i,k+1)+zuo(i,k))

               up_massentro(i,k)=0.5*(entr_rate_2d(i,k)*dz*zuo_ave+up_massentro(i,k-1))

               up_massdetro(i,k)=zuo(i,k)+up_massentro(i,k)-zuo(i,k+1)

               if(up_massdetro(i,k).lt.0.)then
                  up_massdetro(i,k)=0.
                  up_massentro(i,k)=zuo(i,k+1)-zuo(i,k)
                  entr_rate_2d(i,k)=(up_massentro(i,k))/(dz*zuo_ave)
               endif
               if(zuo_ave > 0.) &
                  cd(i,k)=up_massdetro(i,k)/(dz*zuo_ave)
            enddo

            do k=ismooth1,ismooth2
               dz=zo_cup(i,k+1)-zo_cup(i,k)

               zuo_ave = 0.5*(zuo(i,k+1)+zuo(i,k))

               up_massdetro(i,k)=0.5*(cd(i,k)*dz*zuo_ave+up_massdetro(i,k-1))
               up_massentro(i,k)=zuo(i,k+1)-zuo(i,k)+up_massdetro(i,k)

               if(up_massentro(i,k).lt.0.)then
                  up_massentro(i,k)=0.
                  up_massdetro(i,k)=zuo(i,k)-zuo(i,k+1)
                  cd(i,k)=up_massdetro(i,k)/(dz*zuo_ave)
               endif
               if(zuo_ave > 0.) &
                  entr_rate_2d(i,k)=(up_massentro(i,k))/(dz*zuo_ave)
            enddo
          !-----end of the transition zone
         endif
         !=================

         do k=maxloc(zuo(i,:),1)+incr1 ,ktop(i)
            !=> above location of maximum value zu -> change detrainment
            dz=zo_cup(i,k+1)-zo_cup(i,k)
            zuo_ave = 0.5*(zuo(i,k+1)+zuo(i,k))

            up_massentro(i,k)=entr_rate_2d(i,k)*dz*zuo_ave

            up_massdetro(i,k)=zuo(i,k)+up_massentro(i,k)-zuo(i,k+1)
            up_massdetro(i,k)=max(up_massdetro(i,k),0.0)
            !-- check up_massentro in case of dd_up_massdetro has been changed above
            up_massentro(i,k)=-zuo(i,k)+up_massdetro(i,k)+zuo(i,k+1)

            if(zuo_ave.gt.0.) then
               cd          (i,k)=up_massdetro(i,k)/(dz*zuo_ave)
               entr_rate_2d(i,k)=up_massentro(i,k)/(dz*zuo_ave)
            endif
         enddo

         do k=kts,kte
            up_massentr(i,k)=up_massentro(i,k)
            up_massdetr(i,k)=up_massdetro(i,k)
         enddo
         if(present(up_massentru) .and. present(up_massdetru))then
            if(mass_U_option==1) then
               do k=kts+1,kte
                  !--       for weaker mixing
                  up_massentru(i,k-1)=up_massentro(i,k-1)+lambau(i)*up_massdetro(i,k-1)
                  up_massdetru(i,k-1)=up_massdetro(i,k-1)+lambau(i)*up_massdetro(i,k-1)
               !--       for stronger mixing
               ! up_massentru(i,k-1)=up_massentro(i,k-1)+lambau(i)*up_massentro(i,k-1)
               ! up_massdetru(i,k-1)=up_massdetro(i,k-1)+lambau(i)*up_massentro(i,k-1)
               enddo
            else
               turn=maxloc(zuo(i,:),1)
               do k=kts+1,turn
                  up_massentru(i,k-1)=up_massentro(i,k-1)+lambau(i)*up_massentro(i,k-1)
                  up_massdetru(i,k-1)=up_massdetro(i,k-1)+lambau(i)*up_massentro(i,k-1)
               enddo
               do k=turn+1,kte
                  up_massentru(i,k-1)=up_massentro(i,k-1)+lambau(i)*up_massdetro(i,k-1)
                  up_massdetru(i,k-1)=up_massdetro(i,k-1)+lambau(i)*up_massdetro(i,k-1)
               enddo
            endif
         endif
         do k=ktop(i)+1,kte
            cd          (i,k)=0.
            entr_rate_2d(i,k)=0.
         enddo

      enddo ! i
      !---- check mass conservation
      do i=its,itf
         if(ierr(i)/= 0)cycle
         do k=kts+1,kte

            dz =      zo_cup(i,k)-zo_cup(i,k-1)
            dp = 100*(po_cup(i,k)-po_cup(i,k-1))
            rho= -dp/dz/g
            mass1= (zuo(i,k)-zuo(i,k-1)) - up_massentro(i,k-1)+up_massdetro(i,k-1)
            !print*,"masscons=",mass1!,-rho*g*(zuo(i,k)-zuo(i,k-1))/dp, (zuo(i,k)-zuo(i,k-1))/dz,( up_massentro(i,k-1)-up_massdetro(i,k-1))/dz,rho
            mass2= (zuo(i,k)-zuo(i,k-1)) - up_massentru(i,k-1)+up_massdetru(i,k-1)
         enddo
      enddo
   end subroutine get_lateral_massflux

   !------------------------------------------------------------------------------------
   subroutine get_lateral_massflux_down(cumulus,itf,ktf, its,ite, kts,kte              &
      ,ierr,jmin,zo_cup,zdo,xzd,zd,cdd,mentrd_rate_2d      &
      ,dd_massentro,dd_massdetro ,dd_massentr, dd_massdetr &
      ,draft,mentrd_rate,dd_massentru,dd_massdetru,lambau)

      implicit none
      character *(*), intent (in) :: draft,cumulus
      integer, intent(in):: itf,ktf, its,ite, kts,kte
      real,    intent(in)   , dimension(its:ite)         :: mentrd_rate
      integer, intent(in)   , dimension(its:ite)         :: ierr,jmin
      real,    intent(in)   , dimension(its:ite        ) :: lambau
      real,    intent(in)   , dimension(its:ite,kts:kte) :: zo_cup,zdo
      real,    intent(inout), dimension(its:ite,kts:kte) :: cdd,mentrd_rate_2d,xzd,zd
      real,    intent(  out), dimension(its:ite,kts:kte) :: dd_massentro, dd_massdetro&
         ,dd_massentr,  dd_massdetr
      real,    intent(  out), dimension(its:ite,kts:kte), optional &
         :: dd_massentru, dd_massdetru
      integer ::i,ki
      real :: dzo

      cdd          = 0.
      dd_massentr  = 0.
      dd_massdetr  = 0.
      dd_massentro = 0.
      dd_massdetro = 0.
      if(present(dd_massentru).and.present(dd_massdetru))then
         dd_massentru = 0.
         dd_massdetru = 0.
      endif
      if(cumulus == 'shallow') return

      do i=its,itf
         if(ierr(i) /= 0) cycle

         mentrd_rate_2d(i,1:jmin(i))   =mentrd_rate(i)
         cdd           (i,1:jmin(i)-1) =mentrd_rate(i)
         mentrd_rate_2d(i,1)=0.

         do ki=jmin(i)   ,maxloc(zdo(i,:),1),-1

            !=> from jmin to maximum value zd -> change entrainment
            dzo=zo_cup(i,ki+1)-zo_cup(i,ki)
            dd_massdetro(i,ki)=cdd(i,ki)*dzo*zdo(i,ki+1)
            !XXX
            dd_massentro(i,ki)= zdo(i,ki)-zdo(i,ki+1)+dd_massdetro(i,ki)
            dd_massentro(i,ki)= MAX(0.,dd_massentro(i,ki))
            !-- check dd_massdetro in case of dd_massentro has been changed above
            dd_massdetro(i,ki)= dd_massentro(i,ki)-zdo(i,ki)+zdo(i,ki+1)

              !~ if(dd_massentro(i,ki).lt.0.)then
              !~ dd_massentro(i,ki)=0.
              !~ dd_massdetro(i,ki)=zdo(i,ki+1)-zdo(i,ki)
              !~ if(zdo(i,ki+1) > 0.0)&
              !~ cdd(i,ki)=dd_massdetro(i,ki)/(dzo*zdo(i,ki+1))
              !~ endif
              !~ if(zdo(i,ki+1) > 0.0)&
              !~ mentrd_rate_2d(i,ki)=dd_massentro(i,ki)/(dzo*zdo(i,ki+1))
         enddo

         do ki=maxloc(zdo(i,:),1)-1,kts,-1
            !=> from maximum value zd to surface -> change detrainment
            dzo=zo_cup(i,ki+1)-zo_cup(i,ki)
            dd_massentro(i,ki)=mentrd_rate_2d(i,ki)*dzo*zdo(i,ki+1)
            !XXX
            dd_massdetro(i,ki) = zdo(i,ki+1)+dd_massentro(i,ki)-zdo(i,ki)
            dd_massdetro(i,ki) = MAX(0.0,dd_massdetro(i,ki))
            !-- check dd_massentro in case of dd_massdetro has been changed above
            dd_massentro(i,ki) = dd_massdetro(i,ki)+zdo(i,ki)-zdo(i,ki+1)


             !~ if(dd_massdetro(i,ki).lt.0.)then
             !~ dd_massdetro(i,ki)=0.
             !~ dd_massentro(i,ki)=zdo(i,ki)-zdo(i,ki+1)
             !~ if(zdo(i,ki+1) > 0.0)&
             !~ mentrd_rate_2d(i,ki)=dd_massentro(i,ki)/(dzo*zdo(i,ki+1))
             !~ endif
             !~ if(zdo(i,ki+1) > 0.0)&
             !~ cdd(i,ki)= dd_massdetro(i,ki)/(dzo*zdo(i,ki+1))
         enddo

         do ki=jmin(i),kts,-1
            xzd(i,ki)= zdo(i,ki)
            zd (i,ki)= zdo(i,ki)
            dd_massentr(i,ki)=dd_massentro(i,ki)
            dd_massdetr(i,ki)=dd_massdetro(i,ki)
         enddo
         if(present(dd_massentru).and.present(dd_massdetru))then
            do ki=jmin(i),kts,-1
               dd_massentru(i,ki)=dd_massentro(i,ki)+lambau(i)*dd_massdetro(i,ki)
               dd_massdetru(i,ki)=dd_massdetro(i,ki)+lambau(i)*dd_massdetro(i,ki)
            enddo
         endif
      enddo

   end subroutine get_lateral_massflux_down
   !------------------------------------------------------------------------------------

   subroutine get_zi_gf(j,its,ite,kts,kte,istart,iend,ktf,ierr,kzi,pbl,tkeg, &
      rcpg,z,ztop,tkmin)

      implicit none
      integer,intent(in):: j,its,ite,kts,kte, ktf,istart,iend
      integer :: kzimax,ktke_max,i,k
      real tkmin,tke_tmp
      real,    dimension(its:ite,kts:kte) :: tkeg,rcpg,z
      real,    dimension(its:ite)          :: ztop,pbl
      integer, dimension(its:ite)          :: kzi,ierr

      real, parameter :: rcpmin=1.e-5 , pblhmax=3000.

      do i=istart,iend
         kzi(i)  = 1
         !    if(ierr(i).eq.0)then
         !         tke_tmp = 0.
         ktke_max= 1
         kzimax  =ktf-1
         !---  max level for kzi
         do K=kts,ktf
            if(z(i,k).ge. pblhmax+ztop(i)) then
               kzimax = min(k,ktf-1)
               !if(j==8 .and. i==10) print*,"1",z(i,k), pblhmax,ztop(i),kzimax
               exit
            endif
         enddo
         !---

         !         !level of max tke  below kzimax and w/out clouds
         !         do  k=kts,kzimax
         !           !print*,k,tkeg(i,k), tke_tmp,ktke_max,kzimax
         !           if(rcpg(i,k) .lt. rcpmin) then
         !             if( tkeg(i,k) .ge. tke_tmp) then
         !               tke_tmp = tkeg(i,k)
         !               cycle
         !             else
         !               ktke_max= max(1,k-1)
         !               exit
         !             endif
         !           endif
         !         enddo
         !    201         continue
         !             print*,ktke_max

         do k=ktke_max,kzimax+1

            !           if(tkeg(i,k) .gt. 1.1*tkmin .and. rcpg(i,k) .lt. rcpmin )  then
            if(tkeg(i,k) .gt. 1.1*tkmin )  then
               kzi(i) = k
               !if(j==8 .and. i==10) print*,"I",k,rcpg(i,k),tkeg(i,k),kzi(i),z(i,k)-ztop(i)
               cycle

            else
               kzi(i) = max(1,k-1)
               !if(j==8 .and. i==10) print*,"F",k,rcpg(i,k),tkeg(i,k),kzi(i),z(i,k)-ztop(i)
               exit
            endif


         enddo
         kzi(i) = max(1     ,kzi(i))
         !print*,"1",kzi(i),i
         kzi(i) = min(kzimax,kzi(i))
         !print*,"2",kzi(i),icall flush(6)
         pbl(i) = max( z(i,kzi(i))-ztop(i), z(i,1)-ztop(i) )
      enddo

   end subroutine get_zi_gf
   !------------------------------------------------------------------------------------
   subroutine get_zu_zd_pdf(cumulus, draft,ierr,kb,kt,zu,kts,kte,ktf,kpbli,k22,kbcon,klcl,po_cup,psur&
      ,xland,random)

      implicit none
      integer, intent(in   ) :: kts,kte,ktf,kpbli,k22,kbcon,kt,kb,klcl
      integer, intent(inout) :: ierr
      real   , intent(in   ) :: po_cup(kts:kte),psur,xland,random
      real   , intent(inout) :: zu(kts:kte)
      character*(*), intent(in) ::draft,cumulus
      !- local var
      integer :: kk,add,i,nrec=0,k,kb_adj,kpbli_adj,level_max_zu,ktarget
      real :: zumax,ztop_adj,a2,beta, alpha,kratio,tunning,FZU,krmax,dzudk,hei_updf,hei_down
      real :: zuh(kts:kte),zul(kts:kte),  pmaxzu ! pressure height of max zu for deep
      real,   parameter :: px =45./120. ! px sets the pressure level of max zu. its range is from 1 to 120.
      real,   parameter :: px2=45./120. ! px sets the pressure level of max zu. its range is from 1 to 120.
      integer:: itest                   ! 5=gamma+beta, 4=gamma, 1=beta

      integer:: minzu,maxzul,maxzuh,kstart
      logical :: do_smooth

      !-------- gama pdf
      real, parameter :: beta_deep=1.25,g_beta_deep=0.8974707
      integer :: k1
      real :: lev_start,g_alpha2,g_a,y1,x1,g_b,a,b,alpha2,csum,zubeg,wgty,dp_layer,slope
      real, dimension(30) :: x_alpha,g_alpha
      data  (x_alpha(k),k=1,30)/                                    &
         3.699999,3.699999,3.699999,3.699999,          &
         3.024999,2.559999,2.249999,2.028571,1.862500, &
         1.733333,1.630000,1.545454,1.475000,1.415385, &
         1.364286,1.320000,1.281250,1.247059,1.216667, &
         1.189474,1.165000,1.142857,1.122727,1.104348, &
         1.087500,1.075000,1.075000,1.075000,1.075000, &
         1.075000/
      data (g_alpha(k),k=1,30)/                                         &
         4.1706450,4.1706450,4.1706450,4.1706450,          &
         2.0469250,1.3878370,1.1330030,1.012418,0.9494680, &
         0.9153771,0.8972442,0.8885444,0.8856795,0.8865333,&
         0.8897996,0.8946404,0.9005030,0.9070138,0.9139161,&
         0.9210315,0.9282347,0.9354376,0.9425780,0.9496124,&
         0.9565111,0.9619183,0.9619183,0.9619183,0.9619183,&
         0.9619183/
      !-------- gama pdf
      DO_SMOOTH = .false.
      if(USE_SMOOTH_PROF == 1)  DO_SMOOTH = .true.

      !-- fill zu with zeros
      itest=-999
      zu =0.0
      zuh=0.0
      zul=0.0
      if(zero_diff==1) then
         if(draft == "deep_up" .and. xland >  0.90) itest=11 !ocean
         if(draft == "deep_up" .and. xland <= 0.90) itest=12 !land
         if(draft == "mid_up" ) itest= 5
      else

         if(draft == "deep_up"                    ) itest=cum_zuform(deep)  !ocean/land
         if(draft == "mid_up"                     ) itest=cum_zuform(mid )
      endif

      !---------------------------------------------------------
      if(itest==5 .and. draft == "mid_up" ) then

         !--- part 1 GAMMA format
         csum =0.
         zubeg=0.
         lev_start=min(.9,.1+csum*.013)
         kb_adj=max(kb,2)
         kb_adj=min(kb_adj,kt-1)
         if(kb_adj==kt) stop "kb_adj==kt"

         tunning = 0.30
         alpha2  = (tunning*(beta_deep -2.)+1.)/(1.-tunning)

         do k=27,3,-1
            if(x_alpha(k) >= alpha2)exit
         enddo
         k1=k+1
         if(x_alpha(k1) .ne. x_alpha(k1-1))then
            a=x_alpha(k1)-x_alpha(k1-1)
            b=x_alpha(k1-1)*(k1) -(k1-1)*x_alpha(k1)
            x1= (alpha2-b)/a
            y1=a*x1+b
            g_a=g_alpha(k1)-g_alpha(k1-1)
            g_b=g_alpha(k1-1)*k1 - (k1-1)*g_alpha(k1)
            g_alpha2=g_a*x1+g_b
         else
            g_alpha2=g_alpha(k1)
         endif

         fzu = gammaBrams(alpha2 + beta_deep)/(g_alpha2*g_beta_deep)
         fzu=0.01*fzu
         do k=kb_adj,min(kte,kt)
            kratio= (po_cup(k)-po_cup(kb_adj))/(po_cup(kt)-po_cup(kb_adj))
            zu(k) = zubeg+FZU*kratio**(alpha2-1.0) * (1.0-kratio)**(beta_deep-1.0)

         enddo
         !- normalize ZU
         zu(kts:min(kte,kt+1))= zu(kts:min(kte,kt+1))/ (1.e-12+maxval(zu(kts:min(kte,kt+1))))

         !--- part 2: BETA format
         pmaxzu=psur-px*(psur-po_cup(kt))
         kb_adj=minloc(abs(po_cup(kts:kt)-pmaxzu),1)
         kb_adj=max(kb,kb_adj)
         kb_adj=min(kb_adj,kt)
         !beta=4.  !=> must be larger than 1
                   !=> higher makes the profile sharper
                   !=> around the maximum zu
         !- 2nd approach for beta and alpha parameters
         !- the tunning parameter must be between 0.5 (low  level max zu)
         !-                                   and 1.5 (high level max zu)
         !tunning= 1.0
         tunning = 0.6
         !
         beta    = 2.0/tunning
         alpha   = tunning*beta
         !
         !- this alpha constrains the location of the maximun ZU to be at
         !- "kb_adj" vertical level
         alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))
         !
         ! imposing zu(ktop) = 0
         do k=klcl-1,min(kte,kt)
            kratio= float(k)/float(kt+1)
            zuh(k) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
         enddo
         !- normalize ZU
         zuh(kts:min(kte,kt+1))= zuh(kts:min(kte,kt+1))/ (1.e-12+maxval(zuh(kts:min(kte,kt+1))))

         !--- part 3: BETA format from sfc to max zuh, then GAMMA format
         do k=kts,max(kts,maxloc(zuh(:),1)-2)
            zu(k)=zuh(k)
         enddo
         do k=max(kts,maxloc(zuh(:),1)-1),min(maxloc(zuh(:),1)+1,kt)
            zu(k)=0.5*(zu(k)+zuh(k))
         enddo

         !-- special treatment below k22/klcl
         do k=klcl,kts+1,-1
            zu(k)=zu(k+1)*0.5
         enddo
         !-- smooth section
         if(do_smooth) then
            !--from surface
            zul(kts+1)=zu(kts+1)*0.25
            do k=kts+2,maxloc(zu,1)
               zul(k)=(zu(k-1)+zu(k))*0.5
            enddo
            do k=kts+1,maxloc(zu,1)
               zu(k)=(zul(k)+zu(k))*0.5
            enddo

            !--from ktop
            zul(kt-1)=zu(kt-1)*0.1
             !print*,"ZUMD=",kt,zu(kt),zul(kt)
            do k=kt-2,max( kt-min(maxloc(zu,1),5), kts ),-1
               zul(k)=(zul(k+1)+zu(k))*0.5
            enddo
            wgty=0.
            do k=kt,max( kt-min(maxloc(zu,1),5), kts ),-1
               wgty=wgty+1./(float(min(maxloc(zu,1),5))+1)
               zu(k)=zul(k)*(1.-wgty)+ zu(k)*wgty
               !print*,"zuMD=",k,zu(k),zul(k),(zul(k)+zu(k))*0.5,min(maxloc(zu,1),5),wgty
            enddo
         endif
         zu(kts)=0.
      !---------------------------------------------------------
      elseif(itest==20                         ) then       !--- land/ocean

         hei_updf=(1.-xland)*hei_updf_LAND+xland*hei_updf_OCEAN
         !- add a randomic perturbation
         hei_updf = hei_updf + random

         !- for gate soundings
         !hei_updf = max(0.1, min(1.,float(JL)/100.))
         !beta =1.0+float(JL)/100. * 5.

         !--hei_updf parameter goes from 0 to 1 = rainfall decreases with hei_updf
         pmaxzu  =  (psur-100.) * (1.- 0.5*hei_updf) + 0.6*( po_cup(kt) ) * 0.5*hei_updf

         !- beta parameter: must be larger than 1, higher makes the profile sharper around the maximum zu
         !beta    = max(1.1, 2.1 - 0.5*hei_updf)
         beta=2.2
         
         kb_adj=minloc(abs(po_cup(kts:kt)-pmaxzu),1)
         kb_adj=max(kb,kb_adj)
         kb_adj=min(kb_adj,kt)

         !- this alpha constrains the location of the maximun ZU to be at "kb_adj" vertical level
         alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))

         !
         do k=klcl-1,min(kte,kt)
            kratio= float(k)/float(kt+1)
            zu(k) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
         enddo

         !-- special treatment below k22/klcl
         do k=klcl,kts+1,-1
            zu(k)=zu(k+1)*0.5
         enddo

        !if(use_linear_subcl_mf == 1) then
        !    zu(kts)=0.
        !    kstart=kbcon
        !    slope=(zu(kstart)-zu(kts))/(po_cup(kstart)-po_cup(kts) + 1.e-6)
        !    do k=kstart-1,kts+1,-1
        !       zu(k) = zu(kstart)-slope*(po_cup(kstart)-po_cup(k))
        !      !print*,"k=",zu(kstart),zu(k),zu(kts)
        !    enddo
        !    go to 333
        ! endif
  
         !-- smooth section
         if(do_smooth) then
            !--from surface
            zul(kts+1)=zu(kts+1)*0.25
            do k=kts+2,maxloc(zu,1)
               zul(k)=zu(k-1)*0.8+zu(k)*0.2
            enddo
            do k=kts+1,maxloc(zu,1)
               zu(k)=(zul(k)+zu(k))*0.5
            enddo

            !--from ktop
            zul(kt)=zu(kt)*0.1
            do k=kt-1,max( kt-min(maxloc(zu,1),5),kts) ,-1
               zul(k)=(zul(k+1)+zu(k))*0.5
            enddo
            wgty=0.0
            do k=kt,max( kt-min(maxloc(zu,1),5), kts ),-1
               wgty=wgty+1./(float(min(maxloc(zu,1),5))+1)
               zu(k)=zul(k)*(1.-wgty)+ zu(k)*wgty
            enddo
         endif
         zu(kts)=0.
!333 continue         

      !---------------------------------------------------------
      elseif(itest==10) then

         if(xland  < 0.90 ) then !- over land
            hei_updf= hei_updf_LAND
         else
            hei_updf= hei_updf_OCEAN
         endif

         !- for gate soundings
         !if(gate) hei_updf = max(0.1, min(1.,float(JL)/100.))
         !print*,"JL=",jl,hei_updf

         pmaxzu=780.
         kb_adj=minloc(abs(po_cup(kts:kt)-pmaxzu),1)!;print*,"1=",kb_adj,po_cup(kb_adj)
         kb_adj=max(kb,kb_adj)
         kb_adj=min(kb_adj,kt)
         beta    = 5.0

         !- this alpha constrains the location of the maximun ZU to be at "kb_adj" vertical level
         alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))
         do k=klcl-1,min(kte,kt)
            kratio= float(k)/float(kt+1)
            zul(k) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
             !print*,"1",k,zul(k),kb_adj,pmaxzu
         enddo
         zul(kts:min(kte,kt))= zul(kts:min(kte,kt))/ (1.e-9+maxval(zul(kts:min(kte,kt)),1))

         !-----------
         pmaxzu=po_cup(kt)+300.
         kb_adj=minloc(abs(po_cup(kts:kt)-pmaxzu),1)
         kb_adj=max(kb,kb_adj)
         kb_adj=min(kb_adj,kt)
         beta    = 1.5

         !- this alpha constrains the location of the maximun ZU to be at "kb_adj" vertical level
         alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))
         do k=klcl-1,min(kte,kt)
            kratio= float(k)/float(kt+1)
            zuh(k) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
         enddo
         zuh(kts:min(kte,kt))= zuh(kts:min(kte,kt))/ (1.e-9+maxval(zuh(kts:min(kte,kt)),1))

         !increasing contribuition of zuh => more heating at upper levels/less precip
         zu(:)=(1.-hei_updf)*zul(:) + hei_updf*zuh(:)
        
         !-- special treatment below k22/klcl
         do k=klcl,kts+1,-1
            zu(k)=zu(k+1)*0.5
         enddo
         !-- smooth section
         if(do_smooth) then
            !--from surface
            zul(kts+1)=zu(kts+1)*0.25
            do k=kts+2,maxloc(zu,1)
               zul(k)= zu(k-1)*0.8+zu(k)*0.2
            enddo
            do k=kts+1,maxloc(zu,1)
               zu(k)=(zul(k)+zu(k))*0.5
            enddo

            !--from ktop
            zul(kt)=zu(kt)*0.1
            do k=kt-1,max( kt-min(maxloc(zu,1),5), kts ),-1
               zul(k)=(zul(k+1)+zu(k))*0.5
            enddo

            wgty=0.
            do k=kt,max( kt-min(maxloc(zu,1),5), kts ),-1
               wgty=wgty+1./(float(min(maxloc(zu,1),5))+1)
               zu(k)=zul(k)*(1.-wgty)+ zu(k)*wgty
               !print*,"zu=",k,zu(k),zul(k),(zul(k)+zu(k))*0.5,min(maxloc(zu,1),5),wgty
            enddo
         endif
         zu(kts)=0.
      !---------------------------------------------------------
      elseif(draft == "shallow_up") then
         kb_adj   =kts     ! level where mass flux starts
         kpbli_adj=kpbli
         if(kpbli_adj < kb_adj .or. kpbli_adj >= kt ) then
            kpbli_adj = kb_adj + 1
         endif

         !- location of the maximum Zu: dp_layer mbar above PBL height
         !dp_layer     = 10. !mbar
         !level_max_zu = minloc(abs(po_cup(kts:kt+1)-(po_cup(kpbli_adj)-dp_layer)),1)
         !

         k1           = max(kbcon,kpbli_adj)
         !- location of the maximum Zu: dp_layer mbar above k1 height
         hei_updf     =(1.-xland)*hei_updf_LAND+xland*hei_updf_OCEAN

         !hei_updf = (float(JL)-20)/40. ; print*,"JL=",jl,hei_updf

         dp_layer     = hei_updf*(po_cup(k1)-po_cup(kt))

         level_max_zu = minloc(abs(po_cup(kts:kt+1)-(po_cup(k1)-dp_layer)),1)
         level_max_zu = min(level_max_zu,kt -1)
         level_max_zu = max(level_max_zu,kts+1)

         krmax        = float(level_max_zu)/float(kt+1)
         krmax        = min(krmax,0.99)

         beta= beta_sh!smaller => sharper detrainment layer
         !beta= ((1.-xland)*0.43 +xland)*beta_sh

         !beta= 3.0!smaller => sharper detrainment layer
         !beta = 1.+4.*(float(JL))/40.

         !- this alpha imposes the maximum zu at kpbli
         alpha=1.+krmax*(beta-1.)/(1.-krmax)
         !alpha=min(6.,alpha)

         !- to check if dZu/dk = 0 at k=kpbli_adj
         !kratio=krmax
         !dzudk=(alpha-1.)*(kratio**(alpha-2.)) * (1.-kratio)**(beta-1.) - &
         !          (kratio**(alpha-1.))*((1.-kratio)**(beta-2.))*(beta-1.)

         !- Beta PDF
         do k=kts+1,min(kte,kt)
            kratio= float(k)/float(kt+1)
            zu(k)=kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
         enddo
         zu(kts)=0.
         !
         !-- special treatment below kbcon - linear Zu
         if(use_linear_subcl_mf == 1) then
            kstart=kbcon
            slope=(zu(kstart)-zu(kts))/(po_cup(kstart)-po_cup(kts) + 1.e-6)
            do k=kstart-1,kts+1,-1
               zu(k) = zu(kstart)-slope*(po_cup(kstart)-po_cup(k))
              !print*,"k=",zu(kstart),zu(k),zu(kts)
            enddo
         endif
         !-- special treatment below kclcl
         !do k=(klcl-1),kts+1,-1
         !  zu(k)=zu(k+1)*0.5
         !enddo
         !
         !-- smooth section
         !IF( .not. do_smooth) then
         ! zul(kts+1)=zu(kts+1)*0.1
         ! do k=kts+2,maxloc(zu,1)
         !    zul(k)=(zu(k-1)+zu(k))*0.5
         ! enddo
         ! do k=kts+1,maxloc(zu,1)
         !    zu(k)=(zul(k)+zu(k))*0.5
         ! enddo
         !ENDIF
         !zu(kts)=0.

      !---------------------------------------------------------
      elseif(draft == "DOWN" ) then
         if(cumulus == 'shallow') return
         if(cumulus == 'mid' ) beta =2.5
         if(cumulus == 'deep') beta =2.5

         hei_down=(1.-xland)*hei_down_LAND+xland*hei_down_OCEAN

         !---non-zero-diff-APR-08-2020
         if(zero_diff==1) hei_down= 0.5
         !---non-zero-diff-APR-08-2020

         pmaxzu= hei_down * po_cup(kt) + (1.-hei_down)*psur
         kb_adj=minloc(abs(po_cup(kts:kt)-pmaxzu),1)

         !- this alpha constrains the location of the maximun ZU to be at "kb_adj" vertical level
         alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))

         do k=kts+1,min(kt+1,ktf)
            kratio= float(k)/float(kt+1)
            zu(k)= kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
         enddo
         !-- smooth section
         if(do_smooth) then
            zul(kts+1)=zu(kts+1)*0.1
            wgty=0.
            do k=kts+2,maxloc(zu,1)
               wgty=wgty+1./(float(max(2,maxloc(zu,1)))-1.)
               wgty=0.5
               !print*,"zD1=",k,zu(k),zul(k-1),wgty,zul(k-1)*(1.-wgty)+ zu(k)*wgty
               zul(k)=zul(k-1)*(1.-wgty)+ zu(k)*wgty
            enddo
            wgty=0.
            do k=kts+1,maxloc(zu,1)
               wgty=wgty+1./(float(max(2,maxloc(zu,1)))-1.)
               wgty=0.5
               !print*,"zD2=",k,zu(k),zul(k),wgty,zul(k)*(1.-wgty)+ zu(k)*wgty
               zu(k)=zul(k)*(1.-wgty)+ zu(k)*wgty
            enddo
         endif
         zu(kts)=0.

      endif

      !---------------------------------------------------------

      if( maxval(zu(kts:min(kte,kt+1)),1) <= 0.0) then
         zu=0.0
         ierr=51 !ierr(i)=51
      else
         !- normalize ZU
         zu(kts:min(kte,kt+1))= zu(kts:min(kte,kt+1))/ (1.e-9+maxval(zu(kts:min(kte,kt+1)),1))
      endif

   end subroutine get_zu_zd_pdf
    !------------------------------------------------------------------------------------

   subroutine get_zu_zd_pdf_orig(draft,ierr,kb,kt,zs,zuf,ztop,zu,kts,kte,ktf)

      implicit none
      integer, intent(in) ::kb,kt,kts,kte,ktf
      real, intent(in) :: Zs,Zuf,Ztop
      real, intent(inout) :: zu(kts:kte)
      integer, intent(inout) :: ierr
      character*(*), intent(in) ::draft

      !- local var
      integer :: add,i,nrec=0,k,kb_adj
      real ::zumax,ztop_adj
      real ::beta, alpha,kratio,tunning

      !- kb cannot be at 1st level
      kb_adj=max(kb,2)

      !-- fill zu with zeros
      zu=0.0

      if(draft == "UP" .or. draft == "up" ) then
         if(kt<=kb_adj) then
            !stop "ktop must be larger than kbcon"
            ierr=99
            return
         endif
         !beta=4.  !=> must larger than 1
                   !=> higher makes the profile sharper
                   !=> around the maximum zu
         add=0     !=> additional levels above kbcon, where
                   !=> the maximum zu will resides
         kb_adj=kb_adj+add
         kb_adj=max(10,kb_adj)

         !- this alpha constrains the location of the maximun ZU to be at
         !- "kb_adj" vertical level
         !alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))

         !- 2nd approach for beta and alpha parameters
         !- the tunning parameter must be between 0.5 (low  level max zu)
         !-                                   and 1.5 (high level max zu)
         tunning = 0.6
         beta    = 2.0/tunning
         alpha   = tunning*beta

          !- Beta PDF
         do k=kts,min(kte,kt+1)
            kratio= float(k)/float(kt+1)

            zu(k) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
         enddo

      elseif(draft == "DOWN" .or. draft == "DOWNM") then
         add=0    !=> additional levels above kbcon, where
                  !=> the maximum zu will resides
         beta=4.  !=> must larger than 1
                  !=> higher makes the profile sharper
                  !=> around the maximum zu
         alpha= 0.25*beta

         !- for downdrafts kb = jmin(i)-levadj
         kb_adj=kb_adj+add

         !- 2nd approach for beta and alpha parameters
         !- the tunning parameter must be between 0.5 (low  level max zu)
         !-                                   and 1.5 (high level max zu)
         tunning = 1.
         beta    = 2.0/tunning
         alpha   = tunning*beta

         !- Beta PDF
         do k=kts,min(kte,kt)
            kratio= float(k)/float(kt)

            zu(k+1) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
         enddo

      elseif(draft == "shallow" .or. draft == "SHALLOW") then

         alpha= 3.
         beta = 2.*alpha
         kb_adj=1 ! level where mass flux starts

         !- Beta PDF
         do k=kts+kb_adj-1,min(kte,kt+1)
            kratio=float(k+1-kb_adj)/float(kt+1)  !-kb_adj+1)

            zu(k)=kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
         enddo

      else
         print*, "unknown type of flow" ,draft
         stop "routine get_zu_zd"

      endif

      !- normalize ZU
      zu(kts:min(kte,kt+1))= zu(kts:min(kte,kt+1))/ maxval(zu(kts:min(kte,kt+1)))

      !--- Sanity checks
      if(beta <= 1) stop "beta must be larger than 1"

      if(minval(zu(:)) < 0.0 ) then
         print*," zu has negative values for ", draft
         stop   " zu < zero"
      endif
      if(maxval(zu(:)) > 1.0 ) then
         print*," zu has values greater than 1 for ", draft
         stop   " zu  >  one"
      endif

      return

   !OPEN(19,FILE= 'zu.gra', FORM='unformatted',ACCESS='direct'&
   !       ,STATUS='unknown',RECL=4)
   ! DO k = kts,kte
   !    nrec=nrec+1
   !     WRITE(19,REC=nrec) zu(k)
   ! END DO
   !close (19)

   end subroutine get_zu_zd_pdf_orig
   !------------------------------------------------------------------------------------
   subroutine cup_up_cape(aa0,z,zu,dby,GAMMA_CUP,t_cup,                  &
      k22,kbcon,ktop,ierr,tempco,qco,qrco, qo_cup,   &
      itf,ktf,its,ite, kts,kte                       )

      implicit none
      integer ,intent (in   )                   ::        &
         itf,ktf, its,ite, kts,kte

      ! aa0 = dummy array for CAPE (total cape)
      ! gamma_cup = gamma on model cloud levels
      ! t_cup = temperature (Kelvin) on model cloud levels
      ! dby = buoancy term
      ! zu= normalized updraft mass flux
      ! z = heights of model levels
      ! ierr = error value, maybe modified in this routine
      ! tempco = in-cloud temperature (Kelvin) on model cloud levels
      ! qco    = in-cloud water vapor mixing ratio on model cloud levels
      ! qo_cup = environ water vapor mixing ratio on model cloud levels
      ! qrco   = in-cloud liquid water mixing ratio on model cloud levels

      real,    dimension (its:ite,kts:kte) ,intent (in   )    ::       &
         z,zu,gamma_cup,t_cup,dby,tempco,qco,qrco, qo_cup
      integer, dimension (its:ite)         ,intent (in   )    ::       &
         k22,kbcon,ktop
      !
      ! input and output
      integer, dimension (its:ite)         ,intent (inout)    ::       &
         ierr
      real,    dimension (its:ite)         ,intent (out  )    ::       &
         aa0
      !
      !  local variables in this routine
      integer       ::   i,k
      real          ::   dz,daa0
      !
      aa0(:)=0.
      do i=its,itf
         if(ierr(i) == 0) then
            do k=kbcon(i),ktop(i)
               dz=z(i,k)-z(i,max(1,k-1))
               daa0=g*dz*(   (tempco(i,k)*(1.+0.608*qco   (i,k))  - t_cup(i,k)*(1.+0.608*qo_cup(i,k)))&
                  / (t_cup (i,k)*(1.+0.608*qo_cup(i,k))) &
                  )
               aa0(i)=aa0(i)+max(0.,daa0)
              !~ print*,"cape",k,AA0(I),tempco(i,k),t_cup(i,k), qrco  (i,k)
            enddo
         endif
      enddo
   end subroutine cup_up_cape
   !------------------------------------------------------------------------------------

   subroutine FLIPZ(flip,mzp)
      implicit none
      integer, intent(in) :: mzp
      integer, dimension(mzp), intent(inout) :: flip
      integer :: m,k
      m=mzp
      do k=1,mzp
         flip(k)=m
         m=m-1
      enddo
   end subroutine FLIPZ
   !------------------------------------------------------------------------------------

   subroutine set_index_loops( ims,ime, jms,jme, kms,kme,    &
      its,ite, jts,jte, kts,kte,    &
      mxp,myp,mzp                   )

      implicit none
      integer, intent(in)         :: mxp,myp,mzp
      integer, intent(inout)      :: ims,ime, jms,jme, kms,kme,&
         its,ite, jts,jte, kts,kte


      ims=1
      ime=mxp
      jms=1
      jme=myp
      kms=1
      kme=mzp
      its=1
      ite=mxp
      jts=1
      jte=myp
      kts=1
      kte=mzp

   end subroutine set_index_loops
   !------------------------------------------------------------------------------------

   subroutine get_vars(LM,mxp,myp,Q,T,PLE,ZLE,ZLO,PLO,PK,TH)
      implicit none
      integer, intent(in) :: LM,mxp,myp
      real, intent(in) , dimension(mxp,myp,0:LM) :: PLE
      real, intent(in) , dimension(mxp,myp,LM)   :: T,Q

      real, intent(out), dimension(mxp,myp,0:LM) :: ZLE
      real, intent(out), dimension(mxp,myp,LM)   :: ZLO,PLO,PK,TH

      !-local var
      real, parameter:: MAPL_GRAV   = 9.80665                ! m^2/s
      real, parameter:: MAPL_AIRMW  = 28.965                 ! kg/Kmole
      real, parameter:: MAPL_H2OMW  = 18.015                 ! kg/Kmole
      real, parameter:: MAPL_RUNIV  = 8314.47                ! J/(Kmole K)
      real, parameter:: MAPL_RDRY   = MAPL_RUNIV/MAPL_AIRMW  ! J/(kg K)
      real, parameter:: MAPL_CPDRY  = 3.5*MAPL_RDRY          ! J/(kg K)
      real, parameter:: MAPL_KAPPA  = MAPL_RDRY/MAPL_CPDRY   ! (2.0/7.0)
      real, parameter:: MAPL_EPSILON= MAPL_H2OMW/MAPL_AIRMW  ! --
      real, parameter:: MAPL_RGAS   = MAPL_RDRY              ! J/(kg K) (DEPRECATED)
      real, parameter:: MAPL_CP     = MAPL_RGAS/MAPL_KAPPA   ! J/(kg K) (DEPRECATED)
      real, parameter:: MAPL_VIREPS = 1.0/MAPL_EPSILON-1.0   !          (DEPRECATED)
      integer :: L
      real, dimension(mxp,myp,0:LM) :: CNV_PLE,PKE

      CNV_PLE  = PLE*.01
      PLO      = 0.5*(CNV_PLE(:,:,0:LM-1) +  CNV_PLE(:,:,1:LM  ) )
      PKE      = (CNV_PLE/1000.)**(MAPL_RGAS/MAPL_CP)
      PK       = (PLO/1000.)**(MAPL_RGAS/MAPL_CP)
      TH       = T/PK

      ZLE(:,:,LM) = 0.
      do L=LM,1,-1
         ZLE(:,:,L-1) = TH (:,:,L) * (1.+MAPL_VIREPS*Q(:,:,L))
         ZLO(:,:,L  ) = ZLE(:,:,L) + (MAPL_CP/MAPL_GRAV)*( PKE(:,:,L)-PK (:,:,L  ) ) * ZLE(:,:,L-1)
         ZLE(:,:,L-1) = ZLO(:,:,L) + (MAPL_CP/MAPL_GRAV)*( PK (:,:,L)-PKE(:,:,L-1) ) * ZLE(:,:,L-1)
      end do

        !.. IF( MAPL_AM_I_ROOT()) then
             !.. print*,"1get-vars =============================================================="
             !.. do L=LM,1,-1
                 !.. print*,"PLE/PLO",L,PLO(1,1,L),PLE(1,1,L)
             !.. end do
             !.. print*,"PLE/PLO",0,PLO(1,1,1),PLE(1,1,0)
             !.. print*,"2get-vars =============================================================="
             !.. call flush(6)
        !.. ENDIF

   end subroutine get_vars

   !------------------------------------------------------------------------------------
   subroutine get_cloud_bc(cumulus,kts,kte,ktf,xland,po,array,x_aver,k22,add,Tpert)
      implicit none
      character *(*)   ,intent (in) :: cumulus
      integer,intent(in)            :: kts,kte,ktf,k22
      real   ,intent(in)            :: array(kts:kte),po(kts:kte),xland
      real   ,optional ,intent(in)  :: add
      real   ,optional ,intent(in)  :: Tpert(kts:kte)
      real   ,intent(out)           :: x_aver
      integer                       :: i,local_order_aver,order_aver, i_beg,i_end,ic
      real,    parameter            :: frac_ave_layer_ocean= 0.3
      real                          :: count,dp,dp_layer,effec_frac,x_ave_layer

      !-- dimensions of the average:
      !-- a) to pick the value at k22 level, instead of an average between
      !--    (k22-order_aver, ..., k22-1, k22) set order_aver=kts
      !-- b) to average between kts and k22 => set order_aver = k22
      !order_aver = 4    !=> bc_meth 0: average between k22, k22-1, k22-2 ...
                        !=> bc_meth 1: average between ... k22+1,k22, k22-1 ...
      !-- order_aver = kts !=> average between k22, k22-1 and k22-2

      if(bc_meth == 0) then

         order_aver = 3
         local_order_aver=min(k22,order_aver)

         x_aver=0.
         do i = kts,local_order_aver
            x_aver = x_aver + array(k22-i+1)
         enddo
         x_aver = x_aver/float(local_order_aver)

      elseif(bc_meth == 1) then
         effec_frac  = (1.-xland) +xland*frac_ave_layer_ocean
         x_ave_layer = ave_layer*effec_frac


         i_beg = minloc(abs(po(kts:ktf)-(po(k22)+0.5*x_ave_layer)),1)
         i_end = minloc(abs(po(kts:ktf)-(po(k22)-0.5*x_ave_layer)),1)
         i_beg = min(ktf,max(i_beg,kts))
         i_end = min(ktf,max(i_end,kts))

         if(i_beg >= i_end) then
            x_aver   = array(k22)
            dp_layer = 0.
            ic       = i_beg

         else
            dp_layer = 1.e-06
            x_aver   = 0.
            ic       = 0
            do i = i_beg,ktf
               dp = -(po(i+1)-po(i))
               if(dp_layer + dp <= x_ave_layer)  then
                  dp_layer =  dp_layer  + dp
                  x_aver   =  x_aver    + array(i)*dp

               else
                  dp       =  x_ave_layer - dp_layer
                  dp_layer =  dp_layer    + dp
                  x_aver   =  x_aver      + array(i)*dp
                  exit
               endif
            enddo
            x_aver = x_aver/dp_layer
            ic  = max(i_beg,i)
         endif
         !print*,"xaver1=",real(x_aver,4),real(dp_layer,4)

         !-- this perturbation is included only for MSE
         if(present(Tpert)) x_aver = x_aver + cp*maxval(Tpert(i_beg:ic))  ! version 2 - maxval in the layer

      endif
      if(present(add)) x_aver = x_aver + add

   end subroutine get_cloud_bc

   !------------------------------------------------------------------------------------
   subroutine get_lcl(t0,pp0,r0,tlcl,plcl,dzlcl)
      implicit none
      real,intent(in ) :: t0,pp0,r0
      real,intent(out) :: tlcl,plcl,dzlcl

      real, parameter :: &
         cpg      = 102.45             &
         ,   rgas     = 287.               &
         ,   cp       = 1004.              &
         ,   p00      = 1.e5               &
         ,   g        = 9.80               &
         ,   rocp     = rgas / cp          &
         ,   p00i     = 1. / p00           &
         ,   cpor     = cp / rgas          &
         ,   cpi      = 1. / cp            &
         ,   p00k     = 26.870941          &  !  = p00 ** rocp
         ,   p00ki    = 1. / p00k

      integer :: nitt,ip
      real :: p0k,pi0i,ttth0,ttd,dz,pki,pppi,ti,rvs,e
      !real, external :: td,satvap

      !================
      !-simpler,cheaper method
      ttd=td(pp0,r0)
      tlcl=ttd-(0.001296*ttd+0.1963)*(t0-ttd)
      plcl=pp0*(tlcl/t0)**cpor
      dzlcl=127*(t0-ttd)
      if(dzlcl.le.0.)dzlcl=-999.
      !print*,"1st meth",tlcl,plcl,dzlcl;call flush(6)
      return
   !      !================
   !      !-2nd method
   !      dzlcl=-999.
   !      ip=0
   !11 continue
   !
   !   plcl=pp0
   !   tlcl=t0
   !   p0k=pp0**rocp
   !   pi0i=p0k/p00k*cp
   !   ttth0=t0*p00k/p0k
   !   ttd=td(pp0,r0)
   !   dz=cpg*(t0-ttd)
   !   if(dz.le.0.)then
   !      dzlcl=-999.
   !      return
   !   endif
   !   do nitt=1,50
   !      pki=pi0i-g*dz/(ttth0*(1.+.61*r0))
   !      pppi=(pki/cp)**cpor*p00
   !      ti=ttth0*pki/cp
   !      e=100.*satvap(ti)
   !      rvs= ( 0.622*e )/ max(1.e-8,(pppi-e))
   !      !print*,'1',nitt,rvs,r0,ttd,ti,dz
   !      if(abs(rvs-r0).lt..00003)goto 110
   !      ttd=td(pppi,r0)
   !      dz=dz+cp/g*(ti-ttd)
   !        !print*,'2',nitt,rvs-r0,ttd,ti,dz
   !   enddo
   !   print*, 'no converge in LCL:',t0,pp0,r0
   !   ip=ip+1
   !   if(ip==1)go to 11
   !   return
   !
   !110 continue
   !    !- solution for LCL
   !    plcl=pppi
   !    tlcl=ti
   !    dzlcl=dz !displacement
   !    !print*,"2nd meth",tlcl,plcl,dz
   end subroutine get_lcl
   !------------------------------------------------------------------------------------

   real function td(p,rs)
      implicit none
      real :: rr,rs,es,esln,p
      rr=rs+1e-8
      es=p*rr/(.622+rr)
      esln=log(es)
      td=(35.86*esln-4947.2325)/(esln-23.6837)
      return
   end function td
   !------------------------------------------------------------------------------------
   subroutine get_inversion_layers(cumulus,ierr,psur,po_cup,to_cup,zo_cup,k_inv_layers,&
      dtempdz,itf,ktf,its,ite, kts,kte)

      implicit none
      integer,                              intent (in )  :: itf,ktf,its,ite,kts,kte
      character *(*),                       intent (in )  :: cumulus
      integer, dimension (its:ite),         intent (inout):: ierr
      real,    dimension (its:ite),         intent (in )  :: psur
      real,    dimension (its:ite,kts:kte), intent (in )  :: po_cup,to_cup,zo_cup
      real,    dimension (its:ite,kts:kte), intent (out)  :: dtempdz
      integer, dimension (its:ite,kts:kte), intent (out)  :: k_inv_layers
      real :: dzm,delp, first_deriv(kts:kte),sec_deriv(kts:kte),distance(kts:kte)
      integer:: i,k,ilev,kk,k1,ix,k800,k550,ist
      integer, parameter :: extralayer = 0 !- makes plume top higher
      integer, dimension (its:ite,kts:kte)  :: local_k_inv_layers
      !
      !-initialize k_inv_layers as 1 (non-existent layer)_
      k_inv_layers= 1 !integer
      dtempdz     = 0.0
      first_deriv = 0.0
      sec_deriv   = 0.0
      distance    = 0.0
      local_k_inv_layers=1
      ist=3

      do i = its,itf
         if(ierr(i) /= 0) cycle
         !- displacement from local surface pressure level
         delp=1000.-psur(i)

         !- 2nd method
         ! DO k = kts+1,ktf-2
           !dtempdz(i,k)=   ( deriv3(zo_cup(i,k), zo_cup(i,kts:ktf), to_cup(i,kts:ktf), ktf-kts+1, 1,ierr(i)))
           !!! sec_deriv(k)=abs( deriv3(zo_cup(i,k), zo_cup(i,kts:ktf), to_cup(i,kts:ktf), ktf-kts+1, 2))
           !print*,"2=",k,z_cup(i,k),dtempdz(i,k),
         ! ENDDO
         ! if(ierr(i) /= 0) cycle

         !-1st method
         !-  get the 1st derivative
         do k = kts+ist,ktf-ist
            first_deriv(k)  = (to_cup(i,k+1)-to_cup(i,k-1))/(zo_cup(i,k+1)-zo_cup(i,k-1))
         enddo
         first_deriv(kts      :kts+ist-1)  =first_deriv(kts+ist)
         first_deriv(ktf-ist+1:kte      )  =first_deriv(ktf-ist)

         dtempdz  (i,:)  = first_deriv(:)

         !-  get the abs of the 2nd derivative
         do k = kts+ist+1,ktf-ist-1
            sec_deriv(k)= abs((first_deriv(k+1)-first_deriv(k-1))/(zo_cup(i,k+1)-zo_cup(i,k-1)))
         enddo
         sec_deriv(kts    :kts+ist)=sec_deriv(kts+ist+1)
         sec_deriv(ktf-ist:kte    )=sec_deriv(ktf-ist-1)

         ix=1
         do kk=kts+ist+2,ktf-ist-2
            if(sec_deriv(kk) < sec_deriv(kk+1) .and. sec_deriv(kk) < sec_deriv(kk-1)) then
               local_k_inv_layers(i,ix)=kk
               ix  =ix+1
            endif
         enddo

         !- 2nd criteria
         do k=kts+ist+2,ktf-ist-2
            kk=local_k_inv_layers(i,k)
            if(kk == 1) cycle
            if( dtempdz(i,kk) < dtempdz(i,kk-1) .and. dtempdz(i,kk) < dtempdz(i,kk+1) ) then ! the layer is not a local maximum
               local_k_inv_layers(i,k) = 1
            endif
         enddo

      enddo


      !- find the locations of inversions around 800 and 550 hPa
      do i = its,itf
         !----------------
         !k_inv_layers(i,mid)=1
         !----------------
         if(ierr(i) /= 0) cycle
         !- displacement from local surface pressure level
         delp=1000.-psur(i)
         !----------------
         !k_inv_layers(i,mid)=21
         !cycle
         !----------------
         if( trim(cumulus)=='shallow') then
            !- now find the closest layers of 800 and 550 hPa.
            !- this is for shallow convection k800
            do k=kts,ktf
               distance(k)=abs(po_cup(i,local_k_inv_layers(i,k))-(750.-delp))
            enddo
            k800=minloc(abs(distance(kts:ktf)),1)

            if( k800 <= kts .or. k800 >= ktf - 4) then
               k_inv_layers(i,shal)= ktf
               !ierr(i)=8
            else
               !-save k800 in the k_inv_layers array
               k_inv_layers(i,shal)=local_k_inv_layers(i,k800) +extralayer
            endif
            !if(  k_inv_layers(i,shal) <= kts .or. k_inv_layers(i,shal) >= ktf-4) then
         !print*,"SHAL_k_inv_layers=",k_inv_layers(i,shal),ierr(i)
         !ierr(i)=11
            !endif

         elseif( trim(cumulus)=='mid') then
            !- this is for mid/congestus convection k500
            do k=kts,ktf
               distance(k)=abs(po_cup(i,local_k_inv_layers(i,k))-(550.-delp))
            enddo
            k550=minloc(abs(distance(kts:ktf)),1)

            if( k550 <= kts .or. k550 >= ktf - 4) then
               k_inv_layers(i,mid) = 1
               ierr(i)=8
            else
               !-save k550 in the k_inv_layers array
               k_inv_layers(i,mid )=local_k_inv_layers(i,k550) +extralayer
            endif
            if(  k_inv_layers(i,mid) <= kts .or. k_inv_layers(i,mid) >= ktf-4) then
               !print*,"MID_k_inv_layers=",k_inv_layers(i,MID),ierr(i)
               ierr(i)=12
            endif
         else
            k_inv_layers(i,:)=1
            ierr(i)=88
         endif

      enddo

   contains
      real function deriv3(xx, xi, yi, ni, m,ierr)
         !====================================================================
         ! Evaluate first- or second-order derivatives
         ! using three-point Lagrange interpolation
         ! written by: Alex Godunov (October 2009)
         !--------------------------------------------------------------------
         ! input ...
         ! xx    - the abscissa at which the interpolation is to be evaluated
         ! xi()  - the arrays of data abscissas
         ! yi()  - the arrays of data ordinates
         ! ni - size of the arrays xi() and yi()
         ! m  - order of a derivative (1 or 2)
         ! output ...
         ! deriv3  - interpolated value
         !============================================================================*/

         implicit none
         integer, parameter :: n=3
         real   , intent(in):: xx
         integer, intent(in):: ni, m
         real   , intent(in) :: xi(ni), yi(ni)
         real:: x(n), f(n)
         integer i, j, k, ix
         integer, intent(inout) :: ierr

         ! exit if too high-order derivative was needed,
         if (m > 2) then
            deriv3 = 0.0
            return
         end if

         ! if x is ouside the xi(1)-xi(ni) interval set deriv3=0.0
         if (xx < xi(1) .or. xx > xi(ni)) then
            deriv3 = 0.0
            ierr=8
            !stop "problem with 2nd derivative-deriv3 routine"
            return
         endif

         ! a binary (bisectional) search to find i so that xi(i-1) < x < xi(i)
         i = 1
         j = ni
         do while (j > i+1)
            k = (i+j)/2
            if (xx < xi(k)) then
               j = k
            else
               i = k
            endif
         enddo

         ! shift i that will correspond to n-th order of interpolation
         ! the search point will be in the middle in x_i, x_i+1, x_i+2 ...
         i = i + 1 - n/2

         ! check boundaries: if i is ouside of the range [1, ... n] -> shift i
         if (i < 1) i=1
         if (i + n > ni) i=ni-n+1

         !  old output to test i
         !  write(*,100) xx, i
         !  100 format (f10.5, I5)

         ! just wanted to use index i
         ix = i

         ! initialization of f(n) and x(n)
         do i=1,n
            f(i) = yi(ix+i-1)
            x(i) = xi(ix+i-1)
         end do

         ! calculate the first-order derivative using Lagrange interpolation
         if (m == 1) then
            deriv3 =          (2.0*xx - (x(2)+x(3)))*f(1)/((x(1)-x(2))*(x(1)-x(3)))
            deriv3 = deriv3 + (2.0*xx - (x(1)+x(3)))*f(2)/((x(2)-x(1))*(x(2)-x(3)))
            deriv3 = deriv3 + (2.0*xx - (x(1)+x(2)))*f(3)/((x(3)-x(1))*(x(3)-x(2)))
         ! calculate the second-order derivative using Lagrange interpolation
         else
            deriv3 =          2.0*f(1)/((x(1)-x(2))*(x(1)-x(3)))
            deriv3 = deriv3 + 2.0*f(2)/((x(2)-x(1))*(x(2)-x(3)))
            deriv3 = deriv3 + 2.0*f(3)/((x(3)-x(1))*(x(3)-x(2)))
         endif
      end function deriv3
   end subroutine get_inversion_layers
   !------------------------------------------------------------------------------------

   subroutine alloc_grads_arr(n,mzp,task,jl)
      implicit none
      integer, intent(in)    :: n,mzp,task
      integer, intent(inout) :: jl
      integer :: nvar

      if(task == 1) then
         jl = n
         allocate (cupout(nvar_grads))
         do nvar=1,nvar_grads
            allocate(cupout(nvar)%varp(n,mzp))
            allocate(cupout(nvar)%varn(3))
            cupout(nvar)%varp(:,:)=0.0
            cupout(nvar)%varn(:)  ="xxxx"
         enddo
      else
         do nvar=1,nvar_grads
            deallocate(cupout(nvar)%varp)
            deallocate(cupout(nvar)%varn)
         enddo
         deallocate(cupout)
      endif

   end subroutine alloc_grads_arr
   !------------------------------------------------------------------------------------

   subroutine set_grads_var(i,k,nvar,f,name1,name2,name3)
      implicit none
      integer, intent(in)    :: i,k
      integer, intent(inout) :: nvar
      real, intent(in) :: f
      character*(*), intent(in) :: name1,name2,name3

      cupout(nvar)%varp(i,k)= f
      cupout(nvar)%varn(1)=name1
      cupout(nvar)%varn(2)=name2
      cupout(nvar)%varn(3)=name3
      nvar=nvar+1
      if(nvar>nvar_grads) stop 'nvar>nvar_grads'

   end subroutine set_grads_var
   !------------------------------------------------------------------------------------

   subroutine wrt_bin_ctl(n,mzp,p2d,cumulus)
      implicit none
      integer, intent(in):: n,mzp
      character*(*), intent(in) :: cumulus
      real, dimension(mzp),intent(in):: p2d
      integer:: nvartotal,klevgrads(200),jk,int_byte_size,nvar,maxklevgrads
      real   :: real_byte_size
      real, parameter :: undef=-9.99e33
      integer :: nrec=0
      integer :: recSize

      maxklevgrads=min(60,mzp)
      runname='15geos5_'//cumulus
      runlabel=runname

      print*,"writing grads control file:',trim(runname)//'.ctl",ntimes
      call flush(6)
      !
      !number of variables to be written
      nvartotal=0
      do nvar=1,nvar_grads
         if(cupout(nvar)%varn(1) .ne. "xxxx") nvartotal=nvartotal+1
         if(cupout(nvar)%varn(3)  ==  "3d"  ) klevgrads(nvar)=maxklevgrads
         if(cupout(nvar)%varn(3)  ==  "2d"  ) klevgrads(nvar)=1
      enddo

      !- binary file
      inquire (iolength=int_byte_size) real_byte_size  ! inquire by output list

      print*, 'opening grads file:',trim(runname)//'.gra'
      recSize=size(cupout(nvar)%varp,1)*real_byte_size
      if(ntimes == 1) then
         open(19,file= trim(runname)//'.gra',form='unformatted',&
            access='direct',status='replace',recl=recSize)
      else
         open(19,file= trim(runname)//'.gra',form='unformatted',&
            access='direct',status='old', recl=recSize)
      endif

      do nvar=1,nvar_grads
         if(cupout(nvar)%varn(1) .ne. "xxxx") then
            do jk=1,klevgrads(nvar)
               nrec=nrec+1
               !write(19)          real((cupout(nvar)%varp(:,jk)),4)
               write(19,rec=nrec)  real((cupout(nvar)%varp(:,jk)),4)
            enddo
         endif
      enddo
      close (19)
      !-setting vertical dimension '0' for 2d var
      where(klevgrads==1)klevgrads=0
      !- ctl file
      open(20,file=trim(runname)//'.ctl',status='unknown')
      write(20,2001) '^'//trim(runname)//'.gra'
      write(20,2002) 'undef -9.99e33'
      write(20,2002) 'options sequential byteswapped' ! zrev'
      write(20,2002) 'title '//trim(runlabel)
      write(20,2003) 1,0.,1. ! units m/km
      write(20,2004) n,1.,1.
      write(20,2005) maxklevgrads,(p2d(jk),jk=1,maxklevgrads)
      write(20,2006) ntimes,'00:00Z24JAN1999','10mn'
      write(20,2007) nvartotal
      do nvar=1,nvar_grads
         if(cupout(nvar)%varn(1) .ne. "xxxx") then
            !
            write(20,2008) cupout(nvar)%varn(1)(1:len_trim(cupout(nvar)%varn(1))),klevgrads(nvar)&
               ,cupout(nvar)%varn(2)(1:len_trim(cupout(nvar)%varn(2)))
         endif
      enddo
      write(20,2002) 'endvars'
      close(20)

2001  format('dset ',a)
2002  format(a)
2003  format('xdef ',i4,' linear ',2f15.3)
2004  format('ydef ',i4,' linear ',2f15.3)
2005  format('zdef ',i4,' levels ',60f8.3)
2006  format('tdef ',i4,' linear ',2a15)
2007  format('vars ',i4)
2008  format(a10,i4,' 99 ',a40)!'[',a8,']')
2055  format(60f7.0)
133   format (1x,F7.0)

   end subroutine wrt_bin_ctl
   !------------------------------------------------------------------------------------

   subroutine writetxt(mzp,t,ple,th1,pk,q1,u1,zle,zlo      &
      ,DYNF_Q ,DYNF_T ,DYNF_PLE, DYNF_UA  &
      !
      ,theta,pp,rv,up,zm3d,zt3d,vp,omega&
      ,sflux_r,sflux_t,topt,xland,sfc_press,dx,kpbl,temp2m,dt_moist&
      )

      implicit none
      integer, intent(in)    :: mzp
      real, dimension(mzp)   :: t ,th1 ,pk ,q1 ,u1 ,zlo &
         ,DYNF_Q ,DYNF_T , DYNF_UA

      real, dimension(mzp)   :: theta ,pp ,rv ,up ,zm3d ,zt3d,vp,omega
      real, dimension(0:mzp) :: ple ,zle ,DYNF_PLE
      real :: sflux_r,sflux_t,topt,xland,sfc_press,dx,temp2m,dt_moist

      integer :: k,kpbl

      write(8,*) "================================================"
      write(7,*) "================================================"
      write(7,*) kpbl,sflux_r,sflux_t,topt,xland,sfc_press,dx,temp2m,dt_moist
      do k=1,mzp
         write(8,10) k,ple(k),t(k),th1(k),pk(k),1000.*q1(k),u1(k),zle(k),zlo(k),86400.*DYNF_Q(k)
         write(7,11) k,theta(k),pp(k),1000.*rv(k),up(k),zm3d(k),zt3d(k),vp(k),omega(k)
      enddo
      call flush(7)
      call flush(8)
10    format(1x,i4,9F11.3)
11    format(1x,i4,8F11.3)
   end subroutine writetxt
   !------------------------------------------------------------------------------------

   subroutine get_cloud_fraction(  mzp, kts, ktf, &
      PPABS, PZZ, PT, PRV, QCO, QRCO, PMFLX, PCLDFR ,PRC, PRI )
      !!    PURPOSE
      !!    -------
      !!**  Routine to diagnose cloud fraction and liquid and ice condensate mixing ratios
      !!**  METHOD
      !!    ------
      !!    Based on the large-scale fields of temperature, water vapor, and possibly
      !!    liquid and solid condensate, the conserved quantities r_t and h_l are constructed
      !!    and then fractional cloudiness, liquid and solid condensate is diagnosed.
      !!
      !!    The total variance is parameterized as the sum of  stratiform/turbulent variance
      !!    and a convective variance.
      !!    The turbulent variance is parameterized as a function of first-order moments, and
      !!    the convective variance is modelled as a function of the convective mass flux (units kg/s m^2)
      !!    as provided by a mass flux convection scheme.
      !!
      !!    Nota: if the host model does not use prognostic values for liquid and solid condensate
      !!    or does not provide a convective mass flux, put all these values to zero.
      !!    Also, it is supposed that vertical model levels are numbered from
      !!    1 to MZP, where 1 is the first model level above the surface
      !!
      !!    ------------------
      !!    REFERENCE
      !!    ---------
      !!      Chaboureau J.P. and P. Bechtold (J. Atmos. Sci. 2002)
      !!      Chaboureau J.P. and P. Bechtold (JGR/AGU 2005)
      !!
      !!    AUTHOR
      !!    ------
      !!      P. BECHTOLD       * Laboratoire d'Aerologie *
      !!
      !!    MODIFICATIONS
      !!    -------------
      !!      Original    13/06/2001
      !!      modified    20/03/2002 : add convective Sigma_s and improve turbulent
      !!                               length-scale in boundary-layer and near tropopause
      !!      adapted     09/12/2016 : adapted to GEOS-5 by Saulo Freitas
      !-------------------------------------------------------------------------------
      !*       0.    DECLARATIONS
      !              ------------
      implicit none
      !
      !-------------------------------------------------------------------------------
      !
      !*       1.    Set the fundamental thermodynamical constants
      !              these have the same values (not names) as in ARPEGE IFS
      !              -------------------------------------------------------
      real, parameter :: XP00   = 1.e5        ! reference pressure
      real, parameter :: XPI    = 3.141592654 ! Pi
      real, parameter ::  XG    = 9.80665     ! gravity constant
      real, parameter :: XMD    = 28.9644e-3  ! molecular weight of dry air
      real, parameter :: XMV    = 18.0153e-3  ! molecular weight of water vapor
      real, parameter :: XRD    = 287.05967   ! gaz constant for dry air
      real, parameter :: XRV    = 461.524993  ! gaz constant for water vapor
      real, parameter :: XCPD   = 1004.708845 ! specific heat of dry air
      real, parameter :: XCPV   = 1846.1      ! specific heat of water vapor
      real, parameter :: XRHOLW = 1000.       ! density of liquid water
      real, parameter :: XCL    = 4218.       ! specific heat of liquid water
      real, parameter :: XCI    = 2106.       ! specific heat of ice
      real, parameter :: XTT    = 273.16      ! triple point temperature
      real, parameter :: XLVTT  = 2.5008e6    ! latent heat of vaporisation at XTT
      real, parameter :: XLSTT  = 2.8345e6    ! latent heat of sublimation at XTT
      real, parameter :: XLMTT  = 0.3337e6    ! latent heat of melting at XTT
      real, parameter :: XESTT  = 611.14      ! saturation pressure at XTT
      real, parameter :: XALPW  = 60.22416    ! constants in saturation pressure over liquid water
      real, parameter :: XBETAW = 6822.459384
      real, parameter :: XGAMW  = 5.13948
      real, parameter :: XALPI  = 32.62116    ! constants in saturation pressure over ice
      real, parameter :: XBETAI = 6295.421
      real, parameter :: XGAMI  = 0.56313
      logical, parameter :: LUSERI = .true. ! logical switch to compute both
                                                ! liquid and solid condensate (LUSERI=.TRUE.)
                                                ! or only liquid condensate (LUSERI=.FALSE.)
      !
      !*       0.1   Declarations of dummy arguments :
      !
      !
      integer,              intent(in)   :: mzp     ! vertical dimension
      integer,              intent(in)   :: kts     ! vertical  computations start at
      !                                             ! KTS that is at least 1
      integer,              intent(in)   :: ktf     ! vertical computations can be
                                                    ! limited to MZP + 1 - KTF
                                                    ! default=1
      real, dimension(mzp), intent(in)    :: PPABS  ! pressure (Pa)
      real, dimension(mzp), intent(in)    :: PZZ    ! height of model levels (m)
      real, dimension(mzp), intent(in)    :: PT     ! grid scale T  (K)
      real, dimension(mzp), intent(in)    :: PRV    ! grid scale water vapor mixing ratio (kg/kg)
      real, dimension(mzp), intent(in)    :: PMFLX  ! convective mass flux (kg/(s m^2))
      real, dimension(mzp), intent(in)    :: QRCO   ! sub-grid scale liq water mixing ratio (kg/kg)
      real, dimension(mzp), intent(in)    :: QCO    ! in-cloud water mixing ratio (kg/kg)
      real, dimension(mzp), intent(inout),optional :: PRC    ! grid scale r_c mixing ratio (kg/kg)
      real, dimension(mzp), intent(inout),optional :: PRI    ! grid scale r_i (kg/kg)
      real, dimension(mzp), intent(out)   :: PCLDFR ! fractional cloudiness (between 0 and 1)
      !
      !
      !*       0.2   Declarations of local variables :
      !
      integer  ::  JKT, JKP, JKM,K     ! loop index
      real, dimension(mzp) :: ZTLK, ZRT       ! work arrays for T_l, r_t
      real, dimension(mzp) :: ZL              ! length-scale
      integer   :: ITPL    ! top levels of tropopause/highest inversion
      real      :: ZTMIN   ! min Temp. related to ITPL
      real, dimension(mzp) :: LOC_PRC,LOC_PRI
      !
      real :: ZTEMP, ZLV, ZLS, ZTL, ZPV, ZQSL, ZPIV, ZQSI, ZFRAC, ZCOND, ZCPD ! thermodynamics
      real :: ZLL, DZZ, ZZZ ! length scales
      real :: ZAH, ZA, ZB, ZSBAR, ZQ1, ZSIGMA, ZDRW, ZDTL ! related to computation of Sig_s
      real :: ZSIG_CONV,  ZSIGMA_NOCONV,  ZQ1_NOCONV      ! convective part of Sig_s
      !
      !*       0.3  Definition of constants :
      !
      !-------------------------------------------------------------------------------
      !
      real :: ZL0     = 600.        ! tropospheric length scale
                                    ! changed to 600 m instead of 900 m to give a consistent
                                    ! value (linear increase) in general 500 m deep oceanic
                                    ! mixed layer - but could be put back to 900 m if wished
      real :: ZCSIGMA = 0.2         ! constant in sigma_s parameterization
      real :: ZCSIG_CONV = 0.30e-2  ! scaling factor for ZSIG_CONV as function of mass flux
      !
      !
      logical :: ONLY_CONVECTIVE_CLOUD_FRACTION=.true. ! set .false. for the total cloud fraction
      !-------------------------------------------------------------------------------
      !RETURN
      !
      if(PRESENT(PRC)) then
         LOC_PRC(:)=PRC(:)
      else
         LOC_PRC(:)=0.0
      endif
      if(PRESENT(PRI)) then
         LOC_PRI(:)=PRI(:)
      else
         LOC_PRI(:)=0.0
      endif

      PCLDFR(:) = 0. ! Initialize values
      !

      JKT = MZP+1-KTS
      !-will limit the model vertical column to 60 hPa
      do K=KTF,KTS,-1
         if(PPABS(k) > 60.*100.) then
            JKT = k
            !PRINT*,"JKT=",K,MZP+1-KTS ;CALL FLUSH(6)
            exit
         endif
      enddo

      do K=KTS,JKT
         ZTEMP  = PT(k)
          !latent heat of vaporisation/sublimation
         ZLV    = XLVTT + ( XCPV - XCL ) * ( ZTEMP - XTT )
         ZLS    = XLSTT + ( XCPV - XCI ) * ( ZTEMP - XTT )

         !store temperature at saturation and total water mixing ratio
         ZRT(k)   = PRV(k) + LOC_PRC(k) + LOC_PRI(k)
         ZCPD     = XCPD  + XCPV*PRV(k) + XCL*LOC_PRC(k) + XCI*LOC_PRI(k)
         ZTLK(k)  = ZTEMP - ZLV*LOC_PRC(k)/ZCPD - ZLS*LOC_PRI(k)/ZCPD
      end do

      !-------------------------------------------------------------------------------
      ! Determine tropopause/inversion  height from minimum temperature

      ITPL  = KTS+1
      ZTMIN = 400.
      do k = KTS+1,JKT-1
         if ( PT(k) < ZTMIN ) then
            ZTMIN = PT(k)
            ITPL  = K
         endif
      end do

      ! Set the mixing length scale - used for computing the "turbulent part" of Sigma_s

      ZL(:) = 20.
      do k = KTS+1,JKT

         ! free troposphere
         ZL(k) = ZL0
         JKP   = ITPL
         ZZZ   = PZZ(k) -  PZZ(KTS)
            ! approximate length for boundary-layer : linear increase
         if ( ZL0 > ZZZ )  ZL(k) = ZZZ
            ! gradual decrease of length-scale near and above tropopause/top inversion
         if ( ZZZ > 0.9*(PZZ(JKP)-PZZ(KTS)) ) &
            ZL(k) = .6 * ZL(K-1)
      end do
      !-------------------------------------------------------------------------------

      do k=KTS+1,JKT-1
         JKP=k+1
         JKM=k-1
         ZTEMP  = PT(k)

         !latent heat of vaporisation/sublimation
         ZLV    = XLVTT + ( XCPV - XCL ) * ( ZTEMP - XTT )
         ZLS    = XLSTT + ( XCPV - XCI ) * ( ZTEMP - XTT )

         ZCPD   = XCPD + XCPV*PRV(k) + XCL*LOC_PRC(k) + XCI*LOC_PRI(k)
         !temperature at saturation
         ZTL    = ZTEMP - ZLV*LOC_PRC(k)/ZCPD - ZLS*LOC_PRI(k)/ZCPD

         !saturated water vapor mixing ratio over liquid water
         ZPV    = MIN(EXP( XALPW - XBETAW / ZTL - XGAMW * LOG( ZTL ) ),0.99*PPABS(k))
         ZQSL   = XRD / XRV * ZPV / ( PPABS(k) - ZPV )

          !saturated water vapor mixing ratio over ice
         ZPIV   = MIN(EXP( XALPI - XBETAI / ZTL - XGAMI * LOG( ZTL ) ),0.99*PPABS(k))
         ZQSI   = XRD / XRV * ZPIV / ( PPABS(k) - ZPIV )

         !interpolate between liquid and solid as function of temperature
         ! glaciation interval is specified here to 20 K
         ZFRAC = ( ZTL  - 250.16 ) / ( XTT - 250.16 )  ! liquid/solid fraction
         ZFRAC = MAX( 0., MIN(1., ZFRAC ) )

         if(.not. LUSERI) ZFRAC=1.
         ZQSL = ( 1. - ZFRAC ) * ZQSI + ZFRAC * ZQSL
         ZLV  = ( 1. - ZFRAC ) * ZLS  + ZFRAC * ZLV

         !coefficients a and b
         ZAH  = ZLV * ZQSL / ( XRV * ZTL**2 ) * (XRV * ZQSL / XRD + 1.)
         !orig  ZAH  = ZLV * ZQSL / ( XRV * ZTL**2 )

         ZA   = 1. / ( 1. + ZLV/ZCPD * ZAH )
         ZB   = ZAH * ZA

         !- parameterize Sigma_s with first_order closure
         DZZ    =  PZZ (JKP)  - PZZ(JKM)
         ZDRW   =  ZRT (JKP)  - ZRT(JKM)
         ZDTL   =  ZTLK(JKP) - ZTLK(JKM) + XG/ZCPD * DZZ
         ZLL    =  ZL(k)

         !- standard deviation due to convection
         ZSIG_CONV = ZCSIG_CONV * PMFLX(k) / ZA

         !- turb + conv
         ZSIGMA = SQRT( MAX( 1.e-25, ZCSIGMA*ZCSIGMA* ZLL*ZLL/(DZZ*DZZ) * ( &
            ZA*ZA*ZDRW*ZDRW - 2.*ZA*ZB*ZDRW*ZDTL   &
            + ZB*ZB*ZDTL*ZDTL                      ) &
            + ZSIG_CONV * ZSIG_CONV ) )

         !- zsigma should be of order 4.e-4 in lowest 5 km of atmosphere
         ZSIGMA = MAX( ZSIGMA, 1.e-10 )

         !- normalized saturation deficit
         ZSBAR = ZA * ( ZRT (k) - ZQSL )
         !- "Q1" parameter
         ZQ1   = ZSBAR / ZSIGMA

         !- total cloud fraction
         PCLDFR(k) = MAX( 0., MIN(1.,0.5+0.36*ATAN(1.55*ZQ1)) )

         if(ONLY_CONVECTIVE_CLOUD_FRACTION) then
            !- get cloud fraction associated with ONLY the sub-grid scale convective part
            !- this sigma does not include the sub-grid scale convective part
            ZSIGMA_NOCONV = SQRT( MAX( 1.e-25, ZCSIGMA*ZCSIGMA* ZLL*ZLL/(DZZ*DZZ) * ( &
               ZA*ZA*ZDRW*ZDRW - 2.*ZA*ZB*ZDRW*ZDTL   &
               + ZB*ZB*ZDTL*ZDTL  )))
            !- zsigma should be of order 4.e-4 in lowest 5 km of atmosphere
            ZSIGMA_NOCONV = MAX( ZSIGMA_NOCONV, 1.e-10 )
            ZQ1_NOCONV = ZSBAR / ZSIGMA_NOCONV

            !- cloud fraction associated with ONLY convective part ("total-turb")
            PCLDFR(k) = 0.36*(ATAN(1.55*ZQ1)-ATAN(1.55*ZQ1_NOCONV))

            PCLDFR(k) = MAX( 0., MIN(1.,PCLDFR(k)) )

         endif
         !- newer formulation, see GMD 2015
         !PCLDFR(k) = MAX( 0., MIN(1.,0.5+0.34*ATAN(1.85*ZQ1+2.33)) )
         !- this is area fraction of cloud cores
         !PCLDFR(k) = MAX( 0., MIN(1.,0.292/ZQ1**2) )

         cycle
         !total condensate diagnostic (not being used)
         if (ZQ1 > 0. .and. ZQ1 <= 2. ) then
            !orig   ZCOND =     EXP(-1.)+.66*ZQ1+.086*ZQ1*ZQ1
            ZCOND = MIN(EXP(-1.)+.66*ZQ1+.086*ZQ1**2, 2.) ! We use the MIN function for continuity
         else if (ZQ1 > 2.) then
            ZCOND = ZQ1
         else
            ZCOND = EXP( 1.2*ZQ1-1. )
         end if
         ZCOND = ZCOND * ZSIGMA

         if ( zcond < 1.e-12) then
            zcond = 0.
            pcldfr(k) = 0.
         end if
         if ( pcldfr(k) == 0.) then
            zcond = 0.
         end if

         LOC_PRC(k) = ZFRAC * ZCOND ! liquid condensate
         if (LUSERI) then
            LOC_PRI(k) = (1.-ZFRAC) * ZCOND   ! solid condensate
         end if

      !---
      ! compute s'rl'/Sigs^2
      ! used in w'rl'= w's' * s'rl'/Sigs^2
      !  PSIGRC(k) = PCLDFR(k)   ! Gaussian relation
      !
      ! s r_c/ sig_s^2
      !    PSIGRC(JI,JJ,JK) = PCLDFR(JI,JJ,JK)  ! use simple Gaussian relation
      !
      !    multiply PSRCS by the lambda3 coefficient
      !
      !      PSIGRC(JI,JJ,JK) = 2.*PCLDFR(JI,JJ,JK) * MIN( 3. , MAX(1.,1.-ZQ1) )
      ! in the 3D case lambda_3 = 1.
      !      INQ1 = MIN( MAX(-22,FLOOR(2*ZQ1) ), 10)
      !      ZINC = 2.*ZQ1 - INQ1
      !
      !      PSIGRC(k) =  MIN(1.,(1.-ZINC)*ZSRC_1D(INQ1)+ZINC*ZSRC_1D(INQ1+1))
      !
      !      PSIGRC(k) = PSIGRC(k)* MIN( 3. , MAX(1.,1.-ZQ1) )
      !---
      end do
     !
   end subroutine get_cloud_fraction
   !------------------------------------------------------------------------------------
   subroutine cup_cloud_limits(name,ierrc,ierr,cap_inc,cap_max_in                                   &
      ,heo_cup,heso_cup,qo_cup,qeso_cup,po,po_cup,z_cup,heo,hkbo,qo,qeso    &
      ,entr_rate_2d,hcot,k22,kbmax,klcl,kbcon,ktop,depth_neg_buoy,frh,Tpert &
      ,start_level_,use_excess,zqexec,ztexec, x_add_buoy,xland              &
      ,itf,ktf,its,ite, kts,kte)

      implicit none
      character *(*), intent (in)         ::     name

      integer ,intent (in   )             ::    &
         itf,ktf,its,ite, kts,kte,use_excess

      real,    dimension (its:ite,kts:kte)  ,intent (in   )     ::    &
         heo_cup,heso_cup,po_cup,z_cup,heo,qo_cup,qeso_cup,po,qo,qeso,Tpert
      real,    dimension (its:ite)          ,intent (in   )     ::    &
         cap_max_in,cap_inc,xland
      real,    dimension (its:ite)          ,intent (in   )     ::    &
         zqexec,ztexec,x_add_buoy
      integer, dimension (its:ite)          ,intent (in   )     ::    &
         kbmax,start_level_
      integer, dimension (its:ite)          ,intent (inout)     ::    &
         kbcon,ierr,ktop,klcl,k22
      character*128                    ,intent (inout) :: ierrc(its:ite)
      real, dimension (its:ite)        ,intent (inout) :: hkbo,depth_neg_buoy,frh
      real, dimension (its:ite,kts:kte),intent (in)    :: entr_rate_2d
      real, dimension (its:ite,kts:kte),intent (inout) :: hcot

      !  local variables in this routine

      real, parameter              :: frh_crit_O=0.7
      real, parameter              :: frh_crit_L=0.7  !--- test 0.5
      real                         :: delz_oversh !--- height of cloud overshoot is 10% higher than the LNB.
                                                  !--- Typically it can 2 - 2.5km higher, but it depends on
                                                  !--- the severity of the thunderstorm.

      real,    dimension (its:ite) ::   cap_max
      integer                      ::   i,k,k1,k2,kfinalzu
      real                         ::   plus,hetest,dz,dbythresh,denom &
         ,dzh,del_cap_max,fx,x_add,Z_overshoot,frh_crit
      real   , dimension (kts:kte) ::   dby
      integer, dimension (its:ite) ::   start_level

      delz_oversh = OVERSHOOT
      hcot = 0.0
      dby  = 0.0
      start_level = 0
      cap_max(:)  = cap_max_in(:)

      do i=its,itf
         if(ierr(i) /= 0) cycle
         if(ZERO_DIFF==1) then
            start_level(i) = klcl(i)
         else
            start_level(i) = start_level_(i)
         endif

         do k=kts,start_level(i)
            hcot(i,k) = hkbo(i) ! assumed no entraiment between these layers
         enddo
      enddo
      !
      !--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
      !
      loop0: do i=its,itf
         !-default value
         kbcon         (i)=kbmax(i)+3
         depth_neg_buoy(i)=0.
         frh           (i)=0.
         if(ierr(i) /= 0) cycle


         loop1:  do while(ierr(i) == 0)

            kbcon(i)=start_level(i)
            do k=start_level(i)+1,KBMAX(i)+3
               dz=z_cup(i,k)-z_cup(i,k-1)
               hcot(i,k)= ( (1.-0.5*entr_rate_2d(i,k-1)*dz)*hcot(i,k-1)     &
                                  + entr_rate_2d(i,k-1)*dz *heo (i,k-1) ) / &
                          (1.+0.5*entr_rate_2d(i,k-1)*dz)
               if(k==start_level(i)+1) then
                  x_add    = (xlv*zqexec(i)+cp*ztexec(i)) + x_add_buoy(i)
                  hcot(i,k)= hcot(i,k) +  x_add
               endif
            enddo

            loop2:      do while (hcot(i,kbcon(i)) < HESO_cup(i,kbcon(i)))
               kbcon(i)=kbcon(i)+1
               if(kbcon(i).gt.kbmax(i)+2) then
                  ierr(i)=3
                  ierrc(i)="could not find reasonable kbcon in cup_kbcon : above kbmax+2 "
                  exit loop2
               endif
                !print*,"kbcon=",kbcon(i);call flush(6)
            enddo loop2

            if(ierr(i) /= 0) cycle loop0

            !---     cloud base pressure and max moist static energy pressure
            !---     i.e., the depth (in mb) of the layer of negative buoyancy
            depth_neg_buoy(i) = - (po_cup(i,kbcon(i))-po_cup(i,start_level(i)))

            if(MOIST_TRIGGER == 1) then
               frh(i)=0.
               dzh = 0
               do k=k22(i),kbcon(i)
                  dz     = z_cup(i,k)-z_cup(i,max(k-1,kts))
                  frh(i) = frh(i) + dz*(qo(i,k)/qeso(i,k))
                  dzh    = dzh + dz
                  !print*,"frh=", k,dz,qo(i,k)/qeso(i,k)
               enddo
               frh(i) = frh(i)/(dzh+1.e-16)
               frh_crit =frh_crit_O*xland(i) + frh_crit_L*(1.-xland(i))

               !fx     = 4.*(frh(i) - frh_crit)* abs(frh(i) - frh_crit) !-quadratic
               fx     = ((2./0.78)*exp(-(frh(i) - frh_crit)**2)*(frh(i) - frh_crit)) !- exponential
               fx     = max(-1.,min(1.,fx))

               del_cap_max = fx* cap_inc(i)
               cap_max(i)  = min(max(cap_max_in(i) + del_cap_max, 10.),150.)
               !print*,"frh=", frh(i),kbcon(i),del_cap_max, cap_max(i)!,  cap_max_in(i)
            endif

            !- test if the air parcel has enough energy to reach the positive buoyant region
            if(cap_max(i) >= depth_neg_buoy(i)) cycle loop0


            !--- use this for just one search (original k22)
            !            if(cap_max(i) < depth_neg_buoy(i)) then
            !                    ierr(i)=3
            !                    ierrc(i)="could not find reasonable kbcon in cup_cloud_limits"
            !            endif
            !            cycle loop0
            !---

            !- if am here -> kbcon not found for air parcels from k22 level
            k22(i)=k22(i)+1
            !--- increase capmax
            !if(USE_MEMORY == 2000) cap_max(i)=cap_max(i)+cap_inc(i)

            !- get new hkbo
            x_add = (xlv*zqexec(i)+cp*ztexec(i)) +  x_add_buoy(i)
            call get_cloud_bc(name,kts,kte,ktf,xland(i),po(i,kts:kte),heo_cup (i,kts:kte),hkbo (i),k22(i),x_add,Tpert(i,kts:kte))
            !
            start_level(i)=start_level(i)+1
            !
            hcot(i,start_level(i))=hkbo (i)
         enddo loop1
         !--- last check for kbcon
         if(kbcon(i) == kts) then
            ierr(i)=33
            ierrc(i)="could not find reasonable kbcon in cup_kbcon = kts"
         endif
      enddo loop0

      !
      !      if(NAME /= 'shallow') return
      !
      !--- DETERMINE THE LEVEL OF NEUTRAL BUOYANCY - KTOP
      !
      do i=its,itf
         ktop(i) = ktf-1
         if(ierr(i) /= 0) cycle
         !~ dby(:)=0.0

         start_level(i)=kbcon(i)

         do k=start_level(i)+1,ktf-1
            dz=z_cup(i,k)-z_cup(i,k-1)
            denom = 1.+0.5*entr_rate_2d(i,k-1)*dz
            if(denom == 0.) then
               hcot(i,k)=hcot(i,k-1)
            else
               hcot(i,k)=( (1.-0.5*entr_rate_2d(i,k-1)*dz)*hcot(i,k-1)    &
                  +entr_rate_2d(i,k-1)*dz *heo (i,k-1) )/ denom
            endif
         enddo
         do k=start_level(i)+1,ktf-1

            if(hcot(i,k) < heso_cup(i,k) )then
               ktop(i)  =  k - 1
               exit
            endif
         enddo
         if(ktop(i).le.kbcon(i)+1) ierr(i)=41

         !----------------
         if(OVERSHOOT > 1.e-6 .and. ierr(i) == 0) then
            Z_overshoot = (1. + delz_oversh) * z_cup(i,ktop(i))
            do k=ktop(i),ktf-2
               if(Z_overshoot < z_cup(i,k)) then
                  ktop(i) = min(k-1, ktf-2)
                  exit
               endif
            enddo
         endif
      enddo
   end subroutine cup_cloud_limits
   !------------------------------------------------------------------------------------
   subroutine get_buoyancy(itf,ktf, its,ite, kts,kte,ierr,klcl,kbcon,ktop &
      ,hc,he_cup,hes_cup,dby,z_cup)

      implicit none
      integer, intent (in   )   :: itf,ktf,its,ite, kts,kte
      integer, dimension (its:ite)           ,intent (in)     ::      &
         ierr,klcl,kbcon,ktop
      real,    dimension (its:ite,kts:kte)   ,intent (in   )  ::      &
         hc,he_cup,hes_cup,z_cup
      real,    dimension (its:ite,kts:kte)   ,intent (out  )  ::      &
         dby
      integer :: i,k

      do i=its,itf
         dby (i,:)=0.
         if(ierr(i) /= 0) cycle
         do k=kts,klcl(i)
            dby (i,k)=hc (i,k)-he_cup (i,k)
         enddo
         do  k=klcl(i)+1,ktop(i)+1
            dby (i,k)=hc (i,k)-hes_cup (i,k)
         enddo
      enddo
   end subroutine get_buoyancy

   !------------------------------------------------------------------------------------
   subroutine cup_up_vvel(vvel2d,vvel1d,zws,entr_rate_2d,cd ,z,z_cup,zu,dby,GAMMA_CUP,t_cup &
      ,tempco,qco,qrco,qo,start_level,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte &
      ,wlpool,wlpool_bcon,task)

      implicit none
      real, parameter :: ctea=1./3. ,cteb=2., visc=2000., eps=0.622
      integer,intent (in   )              ::  itf,ktf,its,ite, kts,kte,task
      real,    dimension (its:ite,kts:kte) ,intent (in   )  ::  &
         z,z_cup,zu,gamma_cup,t_cup,dby,entr_rate_2d,cd,tempco,qco,qrco,qo

      integer, dimension (its:ite)         ,intent (in   )  ::  &
         klcl,kbcon,ktop,start_level
      real,    dimension (its:ite)         ,intent (in   )  ::  &
         zws
      real,    dimension (its:ite)        ,intent (inout)   ::  wlpool

      ! input and output
      integer, dimension (its:ite)        ,intent (inout) ::  ierr
      real,    dimension (its:ite,kts:kte),intent (out  ) ::  vvel2d
      real,    dimension (its:ite        ),intent (out  ) ::  vvel1d
      real,    dimension (its:ite)        ,intent (inout) ::  wlpool_bcon

      !
      !  local variables in this routine
      integer                             ::  i,k,k1,nvs
      real                                ::  dz,BU,dw2,dw1,kx,dz1m,Tv,Tve,vs,ftun1,ftun2,ke
      real   , parameter :: f=2., C_d=0.506, gam=0.5, beta=1.875 !,ftun1=0.5, ftun2=0.8
      logical, parameter :: smooth=.true.
      integer, parameter :: n_smooth=1

      ftun1=0.25
      ftun2=1.
      if(ZERO_DIFF==1) then
         ftun1=1.
         ftun2=0.5
      endif
      
      if(task == 1) then 

       do i=its,itf
         !-- initialize arrays to zero.
         vvel1d(i  ) = 0.0
         vvel2d(i,:) = 0.0

         if(ierr(i) /= 0) cycle
         vvel2d(i,kts:kbcon(i))= max(1.,max(wlpool_bcon(i)**2,zws(i)**2))

         loop0:  do k= kbcon(i),ktop(i)

            dz=z_cup(i,k+1)-z_cup(i,k)

            Tve= 0.5* ( t_cup (i,k  )*(1.+(qo (i,k  )/eps)/(1.+qo (i,k  ))) + &
               t_cup (i,k+1)*(1.+(qo (i,k+1)/eps)/(1.+qo (i,k+1))))

            Tv = 0.5* ( tempco(i,k  )*(1.+(qco(i,k  )/eps)/(1.+qco(i,k  ))) + &
               tempco(i,k+1)*(1.+(qco(i,k+1)/eps)/(1.+qco(i,k+1)) ))

            BU = g*( (Tv-Tve)/Tve -  ftun2*0.50*(qrco(i,k+1)+qrco(i,k) ))

            dw1 = 2./(f*(1.+gam)) * BU * dz
            if(ZERO_DIFF==1) then
               kx  =               max(entr_rate_2d(i,k),cd(i,k))*dz
            else
               kx  = (1.+beta*C_d)*max(entr_rate_2d(i,k),cd(i,k))*dz*ftun1
            endif

            dw2 =  (vvel2d(i,k)) -2.*kx * (vvel2d(i,k))

            vvel2d(i,k+1)=(dw1+dw2)/(1.+kx)

            if( vvel2d(i,k+1)< 0.) then
               vvel2d(i,k+1) = 0.5* vvel2d(i,k)
            endif

         enddo loop0
       enddo
       if(smooth) then
         if(ZERO_DIFF==1) then
            do i=its,itf
               if(ierr(i) /= 0)cycle
               do k=kts,ktop(i)-2
                  nvs=0
                  vs =0.
                  do k1 = max(k-n_smooth,kts),min(k+n_smooth,ktf)
                     nvs = nvs + 1
                     vs =  vs + vvel2d(i,k1)
                  enddo
                  vvel2d(i,k) = vs/(1.e-16+float(nvs))
               enddo
            enddo
       else
            do i=its,itf
               if(ierr(i) /= 0)cycle
               do k=kts,ktop(i)+1
                  vs =0.
                  dz1m= 0.
                  do k1 = max(k-n_smooth,kts),min(k+n_smooth,ktf)
                     dz   = z_cup(i,k1+1)-z_cup(i,k1)
                     vs   =  vs + dz*vvel2d(i,k1)
                     dz1m = dz1m + dz
                  enddo
                  vvel2d(i,k) = vs/(1.e-16+dz1m)
               !if(k>ktop(i)-3)print*,"v2=",k,ktop(i),sqrt(vvel2d(i,k)),sqrt(vvel2d(i,ktop(i)))
               enddo
            enddo
         endif
       endif

      !-- convert to vertical velocity
       do i=its,itf
         if(ierr(i) /= 0)cycle
         vvel2d(i,:)= sqrt(max(0.1,vvel2d(i,:)))

         if(maxval(vvel2d(i,:)) < 1.0) then
            ierr(i)=54
         !  print*,"ierr=54",maxval(vvel2d(i,:))
         endif


         !-- sanity check
         where(vvel2d(i,:) < 1. ) vvel2d(i,:) = 1.
         where(vvel2d(i,:) > 20.) vvel2d(i,:) = 20.
         if(ZERO_DIFF==0)         vvel2d(i,ktop(i)+1:kte) = 0.1

         !-- get the column average vert velocity
         do k= kbcon(i),ktop(i)
            dz=z_cup(i,k+1)-z_cup(i,k)
            vvel1d(i)=vvel1d(i)+vvel2d(i,k)*dz
         !print*,"w=",k,z_cup(i,k),vvel2d(i,k)
         enddo
         vvel1d(i)=vvel1d(i)/(z_cup(i,ktop(i)+1)-z_cup(i,kbcon(i))+1.e-16)
         vvel1d(i)=max(1.,vvel1d(i))
       enddo
      !-------------------------------------------------------------------
      else
      !-------------------------------------------------------------------
      do i=its,itf
         if(ierr(i) /= 0) cycle
         ke= wlpool(i)**2

         loop1:  do k= start_level(i),kbcon(i)

            dz=z_cup(i,k+1)-z_cup(i,k)

            Tve= 0.5* ( t_cup (i,k  )*(1.+(qo (i,k  )/eps)/(1.+qo (i,k  ))) + &
                        t_cup (i,k+1)*(1.+(qo (i,k+1)/eps)/(1.+qo (i,k+1))))

            Tv = 0.5* ( tempco(i,k  )*(1.+(qco(i,k  )/eps)/(1.+qco(i,k  ))) + &
                        tempco(i,k+1)*(1.+(qco(i,k+1)/eps)/(1.+qco(i,k+1)) ))

            BU = g*( (Tv-Tve)/Tve -  ftun2*0.50*(qrco(i,k+1)+qrco(i,k) ))

            dw1 = 2./(f*(1.+gam)) * BU * dz
            kx  = (1.+beta*C_d)*max(entr_rate_2d(i,k),cd(i,k))*dz*ftun1
            
            dw2 =  ke -2.*kx *ke

            ke=max(0.,(dw1+dw2)/(1.+kx))
            
!            vvel2d(i,k)=sqrt(ke)
         enddo loop1
         wlpool_bcon(i) = sqrt (ke)
         !print*,"wlpool=",wlpool(i),sqrt (ke)
       enddo
      !-------------------------------------------------------------------
      endif

   end subroutine cup_up_vvel
   !------------------------------------------------------------------------------------
   subroutine cup_output_ens_3d(cumulus,xff_shal,xff_mid,xf_ens,ierr,dellat,dellaq,dellaqc,  &
      outtem,outq,outqc,zu,pre,pw,xmb,ktop,                     &
      nx,nx2,ierr2,ierr3,pr_ens, maxens3,ensdim,sig,xland1,     &
      ichoice,ipr,jpr,itf,ktf,its,ite, kts,kte,                 &
      xf_dicycle,outu,outv,dellu,dellv,dtime,po_cup,kbcon,      &
      dellabuoy,outbuoy, dellampqi,outmpqi,dellampql,outmpql,   &
      dellampcf,outmpcf ,nmp, rh_dicycle_fct,xf_coldpool,wlpool_bcon)
      implicit none
      !
      !  on input
      !

      ! only local dimensions are need as of now in this routine

      integer   ,intent (in   )            ::                           &
         ichoice,ipr,jpr,itf,ktf,                                       &
         its,ite, kts,kte
      integer, intent (in   )              ::                           &
         ensdim,nx,nx2,maxens3,nmp
      ! xf_ens = ensemble mass fluxes
      ! pr_ens = precipitation ensembles
      ! dellat = change of temperature per unit mass flux of cloud ensemble
      ! dellaq = change of q per unit mass flux of cloud ensemble
      ! dellaqc = change of qc per unit mass flux of cloud ensemble
      ! outtem = output temp tendency (per s)
      ! outq   = output q tendency (per s)
      ! outqc  = output qc tendency (per s)
      ! pre    = output precip
      ! xmb    = total base mass flux
      ! xfac1  = correction factor
      ! pw = pw -epsilon*pd (ensemble dependent)
      ! ierr error value, maybe modified in this routine
      !
      character *(*), intent (in)          ::    cumulus
      real,    dimension (its:ite,1:ensdim)                             &
         ,intent (inout)                   ::                           &
         xf_ens,pr_ens
      real,    dimension (its:ite,kts:kte)                              &
         ,intent (out  )                   ::                           &
         outtem,outq,outqc,outu,outv,outbuoy

      real,    dimension (nmp,its:ite,kts:kte)                          &
         ,intent (out  )                   ::                           &
         outmpqi,outmpql,outmpcf
      real,    dimension (its:ite,kts:kte)                              &
         ,intent (in  )                   ::                            &
         zu,po_cup
      real,   dimension (its:ite)                                       &
         ,intent (in  )                   ::                            &
         sig, rh_dicycle_fct,wlpool_bcon
      real,   dimension (its:ite,maxens3)                               &
         ,intent (in  )                   ::                            &
         xff_mid
      real,    dimension (its:ite)                                      &
         ,intent (out  )                   ::                           &
         pre,xmb
      real,    dimension (its:ite)                                      &
         ,intent (inout  )                 ::                           &
         xland1
      real,    dimension (its:ite,kts:kte)                              &
         ,intent (inout   )                   ::                           &
         dellat,dellaqc,dellaq,pw,dellu,dellv,dellabuoy
      real,    dimension (nmp,its:ite,kts:kte)                          &
         ,intent (in   )                   ::                           &
         dellampqi,dellampql,dellampcf

      integer, dimension (its:ite)                                      &
         ,intent (in   )                   ::                           &
         ktop,kbcon
      integer, dimension (its:ite)                                      &
         ,intent (inout)                   ::                           &
         ierr,ierr2,ierr3
      real,    intent(inout), dimension (its:ite) :: xf_dicycle,xf_coldpool
      real,    intent(in) :: dtime
      real,   dimension (its:ite,shall_closures)                       &
         ,intent (in  )                   ::                           &
         xff_shal
      !
      !  local variables in this routine
      !

      integer                              ::  i,k,n,ncount,zmax,kk,kqmx,ktmx
      real                                 ::  outtes,ddtes,dtt,dtq,dtqc&
         ,dtpw,prerate,fixouts,dp, xfixQ,xfixT
      real                                 ::  dtts,dtqs,fsum, rcount
      real,    dimension (its:ite)         ::  xmb_ave,xmbmax
      real,    dimension (kts:kte,8)       ::  tend2d
      real,    dimension (8)               ::  tend1d
      real,    dimension (its:ite,8)       ::  check_cons_I,check_cons_F
      !
      do k=kts,ktf
         do i=its,itf
            outtem   (i,k)=0.
            outq     (i,k)=0.
            outqc    (i,k)=0.
            outu     (i,k)=0.
            outv     (i,k)=0.
            outbuoy  (i,k)=0.
         enddo
      enddo
      do i=its,itf
         pre(i)  =0.
         xmb(i)  =0.
         xmb_ave(i) =0.
      enddo

      do i=its,itf
         if(ierr(i).eq.0)then
            do n=1,maxens3
               if(pr_ens(i,n).le.0.)then
                  xf_ens(i,n)=0.
               endif
            enddo
         endif
      enddo
      !
      !--- calculate ensemble average mass fluxes
      !
      if(cumulus == 'deep') then
         do i=its,itf
            if(ierr(i).eq.0)then
               k=0
               xmb_ave(i)=0.
               do n=1,maxens3
                  k=k+1
                  xmb_ave(i)=xmb_ave(i)+xf_ens(i,n)
               enddo
               !- 'ensemble' average mass flux
               xmb_ave(i)=xmb_ave(i)/float(k)
            endif
         enddo

      !- mid (congestus type) convection
      elseif(cumulus=='mid') then
         if(ichoice .le. 3) then
            do i=its,itf
               if(ierr(i) /= 0) cycle
               if(ichoice == 0) then
                  xmb_ave(i)=0.3333*(xff_mid(i,1)+xff_mid(i,2)+xff_mid(i,3))
               else
                  xmb_ave(i)= xff_mid(i,ichoice)
               endif
            enddo
         else
            stop 'For mid ichoice must be 0,1,2,3'
         endif

      !- shallow  convection
      elseif(cumulus=='shallow') then
         do i=its,itf
            if(ierr(i) /= 0) cycle

            if(ichoice > 0) then
               xmb_ave(i)=xff_shal(i,ichoice)
            else
               fsum=0.
               xmb_ave(i)=0.
               do k=1,shall_closures
                  !- heat engine closure is not working properly
                  !- turning it off for now.
                  if(k.ge.4 .and. k.le.6) cycle                  
                  xmb_ave(i)=xmb_ave(i)+xff_shal(i,k)
                  fsum=fsum+1.
               enddo
               !- ensemble average of mass flux
               xmb_ave(i)=xmb_ave(i)/fsum
            endif
         enddo
      endif
      !- apply the mean tropospheric RH control on diurnal cycle (Tian GRL 2022)
       if(cumulus == 'deep' .and. rh_dicycle == 1) then
         do i=its,itf
           if(ierr(i) /= 0) cycle
           xf_dicycle(i) =  xf_dicycle(i) * rh_dicycle_fct(i)
         enddo
      endif


!if(cumulus == 'deep' ) then 
!  do i=its,itf
!!!      if(ierr(i) /= 0) cycle
!      print*,'xmbs',xmb_ave(i),xf_dicycle(i),xf_coldpool(i),wlpool_bcon(i)
!      call flush(6)
!  enddo
!endif


      !- set the updraft mass flux, do not allow negative values and apply the diurnal cycle closure
      do i=its,itf
         if(ierr(i) /= 0) cycle
         !- mass flux of updradt at cloud base
         xmb(i) = xmb_ave(i)

         !- add kinetic energy at the gust front of the cold pools
         xmb(i) = xmb(i) + xf_coldpool(i)

         !- diurnal cycle closure
         xmb(i) = xmb(i) - xf_dicycle(i)
         if(xmb(i) .le. 0.)then
            ierr(i)=13
            xmb (i)=0.
         endif
      enddo
      !-apply the scale-dependence Arakawa's approach
      do i=its,itf
         if(ierr(i) /= 0) cycle
         !- scale dependence
         xmb(i)=sig(i)*xmb(i)

         !- apply the adjust factor for tunning
         !xmb(i) = FADJ_MASSFLX * xmb(i)

         if(xmb(i) == 0. ) ierr(i)=14
         if(xmb(i) > 100.) ierr(i)=15
      enddo

      !--- sanity check for mass flux
      !
      do i=its,itf
         if(ierr(i) /= 0) cycle
         xmbmax(i)=100.*(po_cup(i,kbcon(i))-po_cup(i,kbcon(i)+1))/(g*dtime)
         xmb(i) = min(xmb(i),xmbmax(i))
      enddo

      !--- check outtem and and outq for high values
      !--- criteria: if abs (dT/dt or dQ/dt) > 100 K/day => fix xmb
      if( MAX_TQ_TEND < -1.e-6) then
         do i=its,itf
            if(ierr(i) /= 0) cycle
            fixouts=xmb(i) *86400.*max(maxval(abs(dellat(i,kts:ktop(i)))),&
                              (xlv/cp)*maxval(abs(dellaq(i,kts:ktop(i)))) )

            if(fixouts > abs(MAX_TQ_TEND)) then ! K/day
               fixouts=abs(MAX_TQ_TEND)/(fixouts)
               xmb   (i)  = xmb   (i)  *fixouts
               xf_ens(i,:)= xf_ens(i,:)*fixouts
            endif
         enddo
      endif
       !--- criteria: if abs (dT/dt or dQ/dt) > 100 K/day => fix dT/dt, dQ/dt and xmb
      if( MAX_TQ_TEND > 1.e-6) then
        do i=its,itf

           if(ierr(i) /= 0) cycle
           tend1d=0.
           do k=kts,ktop(i)
              dp         = (po_cup(i,k)-po_cup(i,k+1))
              tend1d(1)  = tend1d(1)  +  dp*xmb(i) * 86400.*(dellat(i,k))
 
              if(xmb(i) * 86400.*abs(dellat(i,k)) > MAX_TQ_TEND ) & 
                 dellat(i,k)=MAX_TQ_TEND/(xmb(i)*86400)*sign(1., dellat(i,k))

              tend1d(2)  = tend1d(2)  +  dp*xmb(i) * 86400.*(dellat(i,k))         
           enddo
      
           do k=kts,ktop(i)
              dp         = (po_cup(i,k)-po_cup(i,k+1))
              tend1d(3)  = tend1d(3)  +  dp*xmb(i) * 86400.*(dellaq(i,k))*(xlv/cp)

              if(xmb(i) * 86400.*abs(dellaq(i,k))*(xlv/cp)  > MAX_TQ_TEND ) &    
                 dellaq(i,k)=MAX_TQ_TEND/(xmb(i)*86400*(xlv/cp))*sign(1., dellaq(i,k))
            
              tend1d(4)  = tend1d(4)  +  dp*xmb(i) * 86400.*(dellaq(i,k))*(xlv/cp)
           enddo
           xfixT = tend1d(1)/(1.e-6+tend1d(2))
           xfixQ = tend1d(3)/(1.e-6+tend1d(4))
         
           xmb(i) = xmb(i)/ max(1.,max(xfixQ,xfixT)) 
          !   print*,"tend",
        enddo
      endif 
      !
      !-- now do feedback
      !
      do i=its,itf
         if(ierr(i) /= 0) cycle
         do k=kts,ktop(i)
            pre      (i)  = pre(i) + pw(i,k)*xmb(i)

            outtem   (i,k)= dellat     (i,k)*xmb(i)
            outq     (i,k)= dellaq     (i,k)*xmb(i)
            outqc    (i,k)= dellaqc    (i,k)*xmb(i)
            outu     (i,k)= dellu      (i,k)*xmb(i)
            outv     (i,k)= dellv      (i,k)*xmb(i)
            outbuoy  (i,k)= dellabuoy  (i,k)*xmb(i)
         enddo
         xf_ens (i,:)= sig(i)*xf_ens(i,:)


         if(APPLY_SUB_MP == 1) then
            do k=kts,ktop(i)
               outmpqi(:,i,k)= dellampqi(:,i,k)*xmb(i)
               outmpql(:,i,k)= dellampql(:,i,k)*xmb(i)
               outmpcf(:,i,k)= dellampcf(:,i,k)*xmb(i)
            enddo
            outmpqi(:,i,ktop(i):ktf)=0.
            outmpql(:,i,ktop(i):ktf)=0.
            outmpcf(:,i,ktop(i):ktf)=0.
         endif
      enddo
      !
      !--  smooth the tendencies (future work: include outbuoy, outmpc* and tracers)
      !
      if(USE_SMOOTH_TEND < 0) then

         do i=its,itf
            if(ierr(i) /= 0) cycle
            tend2d=0.

            !--- get the initial integrals
            rcount = 1.e-6
            tend1d=0.
            do k=kts,ktop(i)
               dp         = (po_cup(i,k)-po_cup(i,k+1))
               rcount     = rcount + dp
               tend1d(1)  = tend1d(1)  +  dp*outtem (i,k)
               tend1d(2)  = tend1d(2)  +  dp*outq   (i,k)
               tend1d(3)  = tend1d(3)  +  dp*outqc  (i,k)
               tend1d(4)  = tend1d(4)  +  dp*outu   (i,k)
               tend1d(5)  = tend1d(5)  +  dp*outv   (i,k)

            enddo
            check_cons_I(i,1:5) = tend1d(1:5)/rcount

            !--- make the smoothness procedure
            do k=kts,ktop(i)
               rcount = 1.e-6
               tend1d=0.

               do kk= max(kts,k-abs(USE_SMOOTH_TEND)),min(ktop(i),k+abs(USE_SMOOTH_TEND))

                  dp=(po_cup(i,kk)-po_cup(i,kk+1))
                  rcount = rcount + dp

                  tend1d(1)  = tend1d(1)  +  dp*outtem (i,kk)
                  tend1d(2)  = tend1d(2)  +  dp*outq   (i,kk)
                  tend1d(3)  = tend1d(3)  +  dp*outqc  (i,kk)
                  tend1d(4)  = tend1d(4)  +  dp*outu   (i,kk)
                  tend1d(5)  = tend1d(5)  +  dp*outv   (i,kk)
               enddo
               tend2d(k,1:5)  = tend1d(1:5) /rcount
            enddo
            !--- get the final/smoother tendencies
            do k=kts,ktop(i)
               outtem (i,k) = tend2d(k,1)
               outq   (i,k) = tend2d(k,2)
               outqc  (i,k) = tend2d(k,3)
               outu   (i,k) = tend2d(k,4)
               outv   (i,k) = tend2d(k,5)
            enddo

            cycle ! the rest is only for checking preservation

            !--- check the final integrals
            rcount = 1.e-6
            tend1d=0.
            do k=kts,ktop(i)
               dp         = (po_cup(i,k)-po_cup(i,k+1))
               rcount     = rcount + dp
               tend1d(1)  = tend1d(1)  +  dp*outtem (i,k)
               tend1d(2)  = tend1d(2)  +  dp*outq   (i,k)
               tend1d(3)  = tend1d(3)  +  dp*outqc  (i,k)
               tend1d(4)  = tend1d(4)  +  dp*outu   (i,k)
               tend1d(5)  = tend1d(5)  +  dp*outv   (i,k)
            enddo
            !--- get the ratio between initial and final integrals.
            check_cons_F(i,1:5) = tend1d(1:5)/rcount

            !--- apply correction to preserve the integrals
            do kk=1,5
               if(abs(check_cons_F(i,kk))>0.) then
                  check_cons_F(i,kk) = abs(check_cons_I(i,kk)/check_cons_F(i,kk))
               else
                  check_cons_F(i,kk) = 1.
               endif
            enddo

            do k=kts,ktop(i)
               outtem (i,k) = outtem (i,k) * check_cons_F(i,1)
               outq   (i,k) = outq   (i,k) * check_cons_F(i,2)
               outqc  (i,k) = outqc  (i,k) * check_cons_F(i,3)
               outu   (i,k) = outu   (i,k) * check_cons_F(i,4)
               outv   (i,k) = outv   (i,k) * check_cons_F(i,5)
            enddo

           !    print*,"check=",real( (check_cons_F(i,3)+check_cons_F(i,2))/ &
           !              (check_cons_I(i,3)+check_cons_I(i,2)),4)!&
           !               ,real( check_cons_F(i,2)/check_cons_I(i,2),4)&
           !             ,real( check_cons_F(i,1)/check_cons_I(i,1),4)
            !-----
         enddo

      endif ! USE_SMOOTH_TEND < 0

   end subroutine cup_output_ens_3d
   !------------------------------------------------------------------------------------

   subroutine cup_forcing_ens_3d(itf,ktf,its,ite, kts,kte,ens4,ensdim,ichoice,maxens,maxens2,maxens3&
      ,ierr,ierr2,ierr3,k22,kbcon,ktop      &
      ,xland,aa0,aa1,xaa0,mbdt,dtime        &
      ,xf_ens,mconv,qo                      &
      ,p_cup,omeg,zd,zu,pr_ens,edt          &
      ,tau_ecmwf,aa1_bl,xf_dicycle,xk_x     &
      ,alpha_adv,Q_adv,aa1_radpbl,aa1_adv,wlpool,xf_coldpool)

      implicit none

      integer                                                           &
         ,intent (in   )                   ::                           &
         itf,ktf,its,ite, kts,kte,ens4
      integer, intent (in   )              ::                           &
         ensdim,maxens,maxens2,maxens3
      !
      ! ierr error value, maybe modified in this routine
      ! pr_ens = precipitation ensemble
      ! xf_ens = mass flux ensembles
      ! massfln = downdraft mass flux ensembles used in next timestep
      ! omeg = omega from large scale model
      ! mconv = moisture convergence from large scale model
      ! zd      = downdraft normalized mass flux
      ! zu      = updraft normalized mass flux
      ! aa0     = cloud work function without forcing effects
      ! aa1     = cloud work function with forcing effects
      ! xaa0    = cloud work function with cloud effects
      ! edt     = epsilon
      ! dir     = "storm motion"
      ! mbdt    = arbitrary numerical parameter
      ! dtime   = dt over which forcing is applied
      ! iact_gr_old = flag to tell where convection was active
      ! kbcon       = LFC of parcel from k22
      ! k22         = updraft originating level
      ! ichoice       = flag if only want one closure (usually set to zero!)
      ! name        = deep or shallow convection flag
      !
      real,    dimension (its:ite,1:ensdim)                     &
         ,intent (inout)                   ::                           &
         pr_ens
      real,    dimension (its:ite,1:ensdim)                     &
         ,intent (out  )                   ::                           &
         xf_ens
      real,    dimension (its:ite,kts:kte)                              &
         ,intent (in   )                   ::                           &
         zd,zu,p_cup,qo
      real,    dimension (its:ite,kts:kte,1:ens4)                       &
         ,intent (in   )                   ::                           &
         omeg
      real,    dimension (its:ite)                                      &
         ,intent (in   )                   ::                           &
         xaa0
      real,    dimension (its:ite)                                      &
         ,intent (in   )                   ::                           &
         aa1,edt,xland
      real,    dimension (its:ite)                                      &
         ,intent (inout)                   ::                           &
         mconv
      real,    dimension (its:ite)                                      &
         ,intent (in   )                   ::                           &
         aa0
      real,    dimension (its:ite)                                      &
         ,intent (in   )                   ::                           &
         mbdt
      real                                                              &
         ,intent (in   )                   ::                           &
         dtime
      integer, dimension (its:ite)                                      &
         ,intent (in   )                   ::                           &
         k22,kbcon,ktop
      integer, dimension (its:ite)                                      &
         ,intent (inout)                   ::                           &
         ierr,ierr2,ierr3
      integer                                                           &
         ,intent (in   )                   ::                           &
         ichoice
      real,    intent(in)   , dimension (its:ite) :: aa1_bl,tau_ecmwf &
         ,alpha_adv,Q_adv,aa1_radpbl,aa1_adv,wlpool
      real,    intent(inout), dimension (its:ite) :: xf_dicycle,xk_x,xf_coldpool
      !- local var
      real  :: xff_dicycle
      !
      !  local variables in this routine
      !

      real,    dimension (1:maxens3)       ::                           &
         xff_ens3
      real,    dimension (its:ite)        ::                           &
         xk
      integer                              ::                           &
         i,k,nall,n,ne,nens,nens3
      real                                 ::                           &
         a1,a_ave,xff0,xomg

      real :: betajb,ke
      integer :: kk
      real, dimension (its:ite) :: ens_adj!,xmbmax
      !
      ens_adj(:)=1.

      !--- LARGE SCALE FORCING
      !
      do i=its,itf
         xf_ens(i,1:16)= 0.
         if(ierr(i) /=  0)cycle

         xff0 = (AA1(I)-AA0(I))/DTIME
         !-- default
         xff_ens3(1) = max(0.,(AA1(I)   -AA0(I))/dtime)

         xff_ens3(2) = xff_ens3(1)
         xff_ens3(3) = xff_ens3(1)
         xff_ens3(16)= xff_ens3(1)
         !
         !--- more like Brown (1979), or Frank-Cohen (199?)
         !--- omeg is in Pa/s
         xomg=0.
         kk=0
         xff_ens3(4)=0.
         do k=max(kts,kbcon(i)-1),kbcon(i)+1
            !-  betajb=(zu(i,k)-edt(i)*zd(i,k))
            betajb=1.
            !if(betajb .gt. 0.)then
            xomg=xomg-omeg(i,k,1)/g/betajb
            kk=kk+1
            !endif
         enddo
         if(kk.gt.0)xff_ens3(4)=xomg/float(kk) ! kg[air]/m^3 * m/s
         xff_ens3(4) = max(0.0, xff_ens3(4))
         xff_ens3(5) = xff_ens3(4)
         xff_ens3(6) = xff_ens3(4)
         xff_ens3(14)= xff_ens3(4)
         !
         !--- more like Krishnamurti et al.;
         !
         !mconv(i) = 0.
         !do k=k22(i),ktop(i)
         !    mconv(i)=mconv(i)+omeg(i,k,1)*(qo(i,k+1)-qo(i,k))/g
         !enddo
         !- 2nd option (assuming that omeg(ktop)*q(ktop)<< omeg(kbcon)*q(kbcon))
         mconv(i)  = -omeg(i,kbcon(i),1)*qo(i,kbcon(i))/g ! (kg[air]/m^3)*m/s*kg[water]/kg[air]

         mconv(i)  = max(0., mconv(i))
         xff_ens3(7) = mconv(i)
         xff_ens3(8) = xff_ens3(7)
         xff_ens3(9) = xff_ens3(7)
         xff_ens3(15)= xff_ens3(7)
         !
         !---- more like  Betchold et al (2014). Note that AA1 already includes the forcings tendencies
         xff_ens3(10)= AA1(i)/tau_ecmwf(i)

         xff_ens3(11)= xff_ens3(10)
         xff_ens3(12)= xff_ens3(10)
         xff_ens3(13)= xff_ens3(10)

         !
         if(ichoice == 0)then
            if(xff0 < 0.)then
               xff_ens3( 1)=0.
               xff_ens3( 2)=0.
               xff_ens3( 3)=0.
               xff_ens3(16)=0.

               xff_ens3(10)=0.
               xff_ens3(11)=0.
               xff_ens3(12)=0.
               xff_ens3(13)=0.
            endif
         endif

         xk(i)=(XAA0(I)-(AA1(I)))/mbdt(i)
         if(xk(i).le.0. .and. xk(i).gt.-0.1*mbdt(i)) xk(i)=-0.1*mbdt(i)
         if(xk(i).gt.0. .and. xk(i).lt.1.e-2       ) xk(i)=1.e-2
         !
         !---  over water, enfor!e small cap for some of the closures
         !
         if(xland(i).lt.0.1)then
            if(ierr2(i).gt.0.or.ierr3(i).gt.0)then
               xff_ens3(1:16) = ens_adj(i)*xff_ens3(1:16)
            endif
         endif
         !
         !--- special treatment for stability closures
         !
         if(xk(i).lt.0.)then
            if(xff_ens3( 1).gt.0.)xf_ens(i, 1)=max(0.,-xff_ens3( 1)/xk(i))
            if(xff_ens3( 2).gt.0.)xf_ens(i, 2)=max(0.,-xff_ens3( 2)/xk(i))
            if(xff_ens3( 3).gt.0.)xf_ens(i, 3)=max(0.,-xff_ens3( 3)/xk(i))
            if(xff_ens3(16).gt.0.)xf_ens(i,16)=max(0.,-xff_ens3(16)/xk(i))
         else
            xff_ens3(1 )=0.
            xff_ens3(2 )=0.
            xff_ens3(3 )=0.
            xff_ens3(16)=0.
         endif

         xf_ens(i,4) =max(0.,xff_ens3(4) )
         xf_ens(i,5) =max(0.,xff_ens3(5) )
         xf_ens(i,6) =max(0.,xff_ens3(6) )
         xf_ens(i,14)=max(0.,xff_ens3(14))

         a1=max(1.e-3,pr_ens(i,7) )
         xf_ens(i,7) =max(0.,xff_ens3(7)/a1)
         a1=max(1.e-3,pr_ens(i,8) )
         xf_ens(i,8) =max(0.,xff_ens3(8)/a1)
         a1=max(1.e-3,pr_ens(i,9) )
         xf_ens(i,9) =max(0.,xff_ens3(9)/a1)
         a1=max(1.e-3,pr_ens(i,15))
         xf_ens(i,15)=max(0.,xff_ens3(15)/a1)
         if(xk(i).lt.0.)then
            xf_ens(i,10)= max(0.,-xff_ens3(10)/xk(i))
            xf_ens(i,11)= max(0.,-xff_ens3(11)/xk(i))
            xf_ens(i,12)= max(0.,-xff_ens3(12)/xk(i))
            xf_ens(i,13)= max(0.,-xff_ens3(13)/xk(i))
         else
            xf_ens(i,10)= 0.
            xf_ens(i,11)= 0.
            xf_ens(i,12)= 0.
            xf_ens(i,13)= 0.
         endif


         if(ichoice.ge.1)then
            xf_ens(i,1:16) =xf_ens(i,ichoice)
         endif

         !---special combination for 'ensemble closure':
         !---over the land, only applies closures 1 and 10.
         !if(ichoice == 0 .and. xland(i) < 0.1)then
         !  xf_ens(i,1:16) =0.5*(xf_ens(i,10)+xf_ens(i,1))
         !endif

         !---over the land, only applies closure 10. (only for GEOS-5)
         !if(zero_diff == 0 .and. ichoice == 0) then
         !   xf_ens(i,1:16)=(1.-xland(i))*xf_ens(i,10)+xland(i)*xf_ens(i,1:16)
         !endif

      !------------------------------------
      enddo
      !-
      !- diurnal cycle mass flux closure
      !-
      if(DICYCLE==1 .or. DICYCLE==2 )then

         do i=its,itf
            xf_dicycle(i) = 0.
!           if(ierr(i) /=  0 .or. p_cup(i,kbcon(i))< 950. )cycle
            if(ierr(i) /=  0)cycle

            !--- Bechtold et al (2014)
            !xff_dicycle  = (AA1(i)-AA1_BL(i))/tau_ecmwf(i)

            !--- Bechtold et al (2014) + Becker et al (2021)
            xff_dicycle  = (1.- alpha_adv(i))* AA1(i) +  alpha_adv(i)*AA1_RADPBL(i) &
                              + alpha_adv(i)*Q_adv(i) - AA1_BL(i)

            !xff_dicycle  = Q_adv(i) 

            xff_dicycle  = xff_dicycle /tau_ecmwf(i)

            if(xk(i).lt.0) xf_dicycle(i)= max(0.,-xff_dicycle/xk(i))
            xf_dicycle(i)= xf_ens(i,10)-xf_dicycle(i)
!----------
!            if(xk(i).lt.0) then 
!                xf_dicycle(i)= max(0.,-xff_dicycle/xk(i))
!                xf_dicycle(i)= xf_ens(i,10)-xf_dicycle(i)
!            else
!                xf_dicycle(i)= 0.
!            endif
!            xf_dicycle(i)= xf_ens(i,10)-xf_dicycle(i)
!----------
         enddo

      elseif( DICYCLE==3) then
        do i=its,itf
            xf_dicycle(i) = 0.

            if(ierr(i) /=  0)cycle

            xff_dicycle  = (1.- alpha_adv(i))* AA1(i)                               & 
!                        +      alpha_adv(i) *(AA1_RADPBL(i) + AA1_ADV(i) - AA0(i)) &
                         +      alpha_adv(i) *(AA1_RADPBL(i) + AA1_ADV(i)         ) &
                         - AA1_BL(i)
!--------tmp           
!           xff_dicycle  =  AA1_ADV(i)
!--------tmp


            xff_dicycle  = xff_dicycle /tau_ecmwf(i)

            if(xk(i).lt.0) xf_dicycle(i)= max(0.,-xff_dicycle/xk(i))
            xf_dicycle(i)= xf_ens(i,10)-xf_dicycle(i)
        enddo
      
      elseif( DICYCLE==4) then
         do i=its,itf
            xf_dicycle(i) = 0.
            if(ierr(i) /=  0)cycle
            !the signal "-" is to convert from Pa/s to kg/m2/s
            if(xk_x(i) > 0.) xf_dicycle(i)= max(0., -AA1_BL(I))/xk_x(i)

            xf_ens    (i,:) = xf_dicycle(i)
            xf_dicycle(i)   = 0.0
         enddo

      else
         xf_dicycle(:)=0.0

      endif
      !------------------------------------
      !-
      !- add the kinetic energy at the gust front at the
      !- mass flux closure
      !-
      if(add_coldpool_clos == 4 )then
         do i=its,itf
            if(ierr(i) /= 0 .or. xk(i) >= 0) cycle
            xf_coldpool(i) = -(0.5*wlpool(i)**2/tau_ecmwf(i)) /xk(i)
         enddo
      endif 
   end subroutine cup_forcing_ens_3d

   !------------------------------------------------------------------------------------
   subroutine get_partition_liq_ice(ierr,tn,z1,zo_cup,po_cup, p_liq_ice,melting_layer         &
      ,itf,ktf,its,ite, kts,kte, cumulus          )
      implicit none
      character *(*), intent (in)                          :: cumulus
      integer  ,intent (in   )                             :: itf,ktf, its,ite, kts,kte
      integer  ,intent (in   ), dimension(its:ite)         :: ierr
      real     ,intent (in   ), dimension(its:ite)         :: z1
      real     ,intent (in   ), dimension(its:ite,kts:kte) :: tn,po_cup,zo_cup
      real     ,intent (inout), dimension(its:ite,kts:kte) :: p_liq_ice,melting_layer
      integer :: i,k
      real    :: dp, height
      real, dimension(its:ite) :: norm
      real, parameter ::  T1=276.16, Z_meltlayer1=4000.,Z_meltlayer2=6000.,delT=3.
      p_liq_ice    (:,:) = 1.
      melting_layer(:,:) = 0.
      !-- get function of T for partition of total condensate into liq and ice phases.
      if(MELT_GLAC .and. cumulus == 'deep') then
         do k=kts,ktf
            do i=its,itf
               if(ierr(i) /= 0) cycle
               p_liq_ice(i,k) = fract_liq_f(tn(i,k))
            enddo
         enddo
         !        go to 650
         !
         !-- define the melting layer (the layer will be between T_0+1 < TEMP < T_1
         !-- definition em terms of temperatura
         do k=kts,ktf
            do i=its,itf
               if(ierr(i) /= 0) cycle
               if    (tn(i,k) <= T_0-delT) then
                  melting_layer(i,k) = 0.

               elseif(  tn(i,k) < T_0+delT .and. tn(i,k) > T_0-delT) then
                  melting_layer(i,k) =  ((tn(i,k)-(T_0-delt))/(2.*delT))**2

               else
                  melting_layer(i,k) = 1.
               endif
               melting_layer(i,k) = melting_layer(i,k)*(1.-melting_layer(i,k))
            enddo
         enddo
         !go to 655
         !650 continue
         !        !-- definition em terms of height above local terrain
         !        DO k=kts,ktf
         !          DO i=its,itf
         !             if(ierr(i) /= 0) cycle
         !             height= zo_cup(i,k)+z1(i)
         !             if   (height > Z_meltlayer2 ) then
         !                melting_layer(i,k) = 0.
         !
         !             elseif(height > Z_meltlayer1  .and. height < Z_meltlayer2 ) then
         !
         !                melting_layer(i,k) =  ((height - Z_meltlayer1)/(Z_meltlayer2-Z_meltlayer1))**2.
         !
         !
         !             else
         !                melting_layer(i,k) = 1.
         !             endif
         !             melting_layer(i,k) = melting_layer(i,k)*(1.-melting_layer(i,k))
         !          ENDDO
         !        ENDDO
         !
         !
         !         655 continue
         !-normalize vertical integral of melting_layer to 1
         norm(:)=0.
         do k=kts,ktf-1
            do i=its,itf
               if(ierr(i) /= 0) cycle
               dp = 100.*(po_cup(i,k)-po_cup(i,k+1))
               norm(i) = norm(i) + melting_layer(i,k)*dp/g
            enddo
         enddo
         do i=its,itf
            if(ierr(i) /= 0) cycle
            melting_layer(i,:)=melting_layer(i,:)/(norm(i)+1.e-6)*(100*(po_cup(i,kts)-po_cup(i,ktf))/g)
          !print*,"i2=",i,maxval(melting_layer(i,:)),minval(melting_layer(i,:)),norm(i)
         enddo
         !--check
      !       norm(:)=0.
      !        DO k=kts,ktf-1
      !          DO i=its,itf
      !             dp = 100.*(po_cup(i,k)-po_cup(i,k+1))
      !             norm(i) = norm(i) + melting_layer(i,k)*dp/g/(100*(po_cup(i,kts)-po_cup(i,ktf))/g)
      !             !print*,"n=",i,k,norm(i)
      !          ENDDO
      !        ENDDO

      !~ ELSE
         !~ p_liq_ice    (:,:) = 1.
         !~ melting_layer(:,:) = 0.
      endif
   end  subroutine get_partition_liq_ice

   !------------------------------------------------------------------------------------
   subroutine get_melting_profile(ierr,tn_cup,po_cup, p_liq_ice,melting_layer,qrco    &
      ,pwo,edto,pwdo,melting                               &
      ,itf,ktf,its,ite, kts,kte, cumulus                   )
      implicit none
      character *(*), intent (in)                          :: cumulus
      integer  ,intent (in   )                             :: itf,ktf, its,ite, kts,kte
      integer  ,intent (in   ), dimension(its:ite)         :: ierr
      real     ,intent (in   ), dimension(its:ite)         :: edto
      real     ,intent (in   ), dimension(its:ite,kts:kte) :: tn_cup,po_cup,qrco,pwo &
         ,pwdo,p_liq_ice,melting_layer
      real     ,intent (inout), dimension(its:ite,kts:kte) :: melting
      integer :: i,k
      real    :: dp
      real, dimension(its:ite)         :: norm,total_pwo_solid_phase
      real, dimension(its:ite,kts:kte) :: pwo_solid_phase,pwo_eff

      if(MELT_GLAC .and. cumulus == 'deep') then

         norm                  = 0.0
         pwo_solid_phase       = 0.0
         pwo_eff               = 0.0
         melting               = 0.0
         !-- set melting mixing ratio to zero for columns that do not have deep convection
         do i=its,itf
            if(ierr(i) > 0) melting(i,:) = 0.
         enddo

         !-- now, get it for columns where deep convection is activated
         total_pwo_solid_phase(:)=0.

         do k=kts,ktf-1
            do i=its,itf
               if(ierr(i) /= 0) cycle
               dp = 100.*(po_cup(i,k)-po_cup(i,k+1))

               !-- effective precip (after evaporation by downdraft)
               !-- pwdo is not defined yet
               !pwo_eff(i,k) = 0.5*(pwo(i,k)+pwo(i,k+1) + edto(i)*(pwdo(i,k)+pwdo(i,k+1)))
               pwo_eff(i,k) = 0.5*(pwo(i,k)+pwo(i,k+1))

               !-- precipitation at solid phase(ice/snow)
               pwo_solid_phase(i,k) = (1.-p_liq_ice(i,k))*pwo_eff(i,k)

               !-- integrated precip at solid phase(ice/snow)
               total_pwo_solid_phase(i) = total_pwo_solid_phase(i)+pwo_solid_phase(i,k)*dp/g
            enddo
         enddo

         do k=kts,ktf
            do i=its,itf
               if(ierr(i) /= 0) cycle
               !-- melting profile (kg/kg)
               melting(i,k) = melting_layer(i,k)*(total_pwo_solid_phase(i)/(100*(po_cup(i,kts)-po_cup(i,ktf))/g))
               !print*,"mel=",k,melting(i,k),pwo_solid_phase(i,k),po_cup(i,k)
            enddo
         enddo

      !-- check conservation of total solid phase precip
      !       norm(:)=0.
      !        DO k=kts,ktf-1
      !          DO i=its,itf
      !             dp = 100.*(po_cup(i,k)-po_cup(i,k+1))
      !             norm(i) = norm(i) + melting(i,k)*dp/g
      !          ENDDO
      !        ENDDO
      !
      !       DO i=its,itf
      !         print*,"cons=",i,norm(i),total_pwo_solid_phase(i)
      !        ENDDO
      !--

      else
         !-- no melting allowed in this run
         melting     (:,:) = 0.
      endif
   end  subroutine get_melting_profile
   !------------------------------------------------------------------------------------
   subroutine  ke_to_heating(itf,ktf,its,ite, kts,kte,ktop,ierr &
      ,po_cup,us,vs,dellu,dellv,dellat)

      implicit none
      integer                              ,intent (in   ) :: itf,ktf,its,ite, kts,kte
      integer, dimension (its:ite)         ,intent (in   ) :: ierr,ktop
      real   , dimension (its:ite,kts:kte) ,intent (in   ) :: po_cup,us,vs,dellu,dellv
      real   , dimension (its:ite,kts:kte) ,intent (inout) :: dellat

      real :: dts,fp,dp,fpi
      integer ::i,k

      ! since kinetic energy is being dissipated, add heating accordingly (from ECMWF)
      !
      do i=its,itf
         if(ierr(i) /= 0) cycle
         dts=0.
         fpi=0.
         do k=kts,ktop(i)
            dp=(po_cup(i,k)-po_cup(i,k+1))*100.
            !total KE dissiptaion estimate
            dts = dts - (dellu(i,k)*us(i,k)+dellv(i,k)*vs(i,k))*dp/g
            !
            ! fpi needed for calcualtion of conversion to pot. energyintegrated
            fpi = fpi + sqrt(dellu(i,k)*dellu(i,k) + dellv(i,k)*dellv(i,k))*dp
         enddo
         if(fpi.gt.0.)then
            do k=kts,ktop(i)
               fp= sqrt((dellu(i,k)*dellu(i,k)+dellv(i,k)*dellv(i,k)))/fpi
               dellat(i,k) = dellat(i,k) + fp*dts*g/cp
            enddo
         endif
      enddo

   end subroutine  ke_to_heating
   !------------------------------------------------------------------------------------
   function auto_rk(n,step,aux,xexp,qrc1) result(PW)
      integer, intent(in) :: n
      real   , intent(in) :: step,aux,qrc1,xexp
      real                :: PW

      PW=step*qrc1*(1.0-exp(-aux**xexp))/float(n)

   end function auto_rk

   !------------------------------------------------------------------------------------
   subroutine get_incloud_sc_chem_up(cumulus,fscav,mtp,se,se_cup,sc_up,pw_up,tot_pw_up_chem&
      ,z_cup,rho,po,po_cup  &
      ,qrco,tempco,pwo,zuo,up_massentro,up_massdetro,vvel2d,vvel1d  &
      ,start_level,k22,kbcon,ktop,klcl,ierr,xland,itf,ktf,its,ite, kts,kte)
      implicit none
      !-inputs
      integer                               ,intent (in)  :: itf,ktf, its,ite, kts,kte
      integer                               ,intent (in)  :: mtp
      integer, dimension (its:ite)          ,intent (in)  :: ierr ,kbcon,ktop,k22,klcl,start_level
      character *(*)                        ,intent (in)  :: cumulus
      real, dimension (mtp)                 ,intent (in)  :: FSCAV
      real, dimension (mtp ,its:ite,kts:kte),intent (in)  :: se,se_cup
      real, dimension (its:ite,kts:kte)     ,intent (in)  :: z_cup,rho,po_cup,qrco,tempco,pwo,zuo &
         ,up_massentro,up_massdetro,po


      real,    dimension (its:ite,kts:kte)  ,intent (in)  :: vvel2d
      real,    dimension (its:ite        )  ,intent (in)  :: vvel1d,xland

      !-outputs
      real, dimension (mtp ,its:ite,kts:kte),intent (out) :: sc_up,pw_up
      real, dimension (mtp ,its:ite        ),intent (out) :: tot_pw_up_chem

      !-locals
      real, parameter :: scav_eff = 0.6  ! for smoke : Chuang et al. (1992) J. Atmos. Sci.
      real   , dimension (mtp ,its:ite) ::  sc_b
      real   , dimension (mtp) :: conc_mxr
      real :: x_add,dz,XZZ,XZD,XZE,denom,henry_coef,w_upd,fliq,dp
      integer :: i,k,ispc
      real, parameter :: cte_w_upd = 10. ! m/s
      !    real, parameter :: kc = 5.e-3  ! s-1
      real, parameter :: kc = 2.e-3  ! s-1        !!! autoconversion parameter in GF is lower than what is used in GOCART
      real, dimension (mtp ,its:ite,kts:kte) ::  factor_temp

      !--initialization
      sc_up          = se_cup
      pw_up          = 0.0
      tot_pw_up_chem = 0.0

      if(USE_TRACER_SCAVEN==2 .and. cumulus /= 'shallow') then
         factor_temp = 1.
         do i=its,itf
            if(ierr(i) /= 0) cycle
            do ispc = 1,mtp
               ! - if tracer is type "carbon" then set coefficient to 0 for hydrophobic
               if( TRIM(CHEM_name (ispc)(1:len_trim('OCphobic') )) == 'OCphobic') factor_temp(ispc,:,:) = 0.0

               ! - suppress scavenging most aerosols at cold T except BCn1 (hydrophobic), dust, and HNO3
               if( TRIM(CHEM_name (ispc)(1:len_trim('BCphobic') )) == 'BCphobic') then
                  where(tempco < 258.) factor_temp(ispc,:,:) = 0.0
               endif

               if( TRIM(CHEM_name (ispc)) == 'sulfur'   .or. &

                  TRIM(CHEM_name (ispc)(1:len_trim('ss') )) == 'ss'  .or. & ! 'seasalt'

                  TRIM(CHEM_name (ispc)) == 'SO2'      .or. &
                  TRIM(CHEM_name (ispc)) == 'SO4'      .or. &

                  TRIM(CHEM_name (ispc)) == 'nitrate'  .or. &
                  TRIM(CHEM_name (ispc)) == 'bromine'  .or. &
                  TRIM(CHEM_name (ispc)) == 'NH3'      .or. &
                  TRIM(CHEM_name (ispc)) == 'NH4a'          ) then

                  where(tempco < 258.) factor_temp(ispc,:,:) = 0.0
               endif

            enddo
         enddo
      endif

      do i=its,itf
         if(ierr(i) /= 0) cycle
         !start_level(i) = klcl(i)
         !start_level(i) = k22(i)

         do ispc=1,mtp
            call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),se_cup(ispc,i,kts:kte),sc_b(ispc,i),k22(i))
         enddo
         do k=kts,start_level(i)
            sc_up   (:,i,k) = sc_b  (:,i)
             !sc_up   (:,i,k) = se_cup(:,i,k)
         enddo
      enddo

      do i=its,itf
         if(ierr(i) /= 0) cycle
         loopk:      do k=start_level(i)+1,ktop(i)+1

            !-- entr,detr, mass flux ...
            XZZ=             zuo(i,k-1)
            XZD=0.5*up_massdetro(i,k-1)
            XZE=    up_massentro(i,k-1)
            denom =  (XZZ-XZD+XZE)

            !-- transport + mixing
            if(denom > 0.) then
               sc_up(:,i,k) = (sc_up(:,i,k-1)*XZZ - sc_up(:,i,k-1)*XZD + se(:,i,k-1)*XZE) / denom
            else
               sc_up(:,i,k) = sc_up(:,i,k-1)
            endif

            !-- scavenging section
            if(USE_TRACER_SCAVEN==0 .or. cumulus == 'shallow') cycle loopk
            dz=z_cup(i,k)-z_cup(i,k-1)

            !-- in-cloud vert velocity for scavenging formulation 2
            !           w_upd = cte_w_upd
            !           w_upd = vvel1d(i)
            w_upd = vvel2d(i,k)

            do ispc = 1,mtp
               if(fscav(ispc) > 1.e-6) then ! aerosol scavenging

                  !--formulation 1 as in GOCART with RAS conv_par
                  if(USE_TRACER_SCAVEN==1) &
                     pw_up(ispc,i,k) = max(0.,sc_up(ispc,i,k)*(1.-exp(- FSCAV(ispc) * (dz/1000.))))

                  !--formulation 2 as in GOCART
                  if(USE_TRACER_SCAVEN==2) &
                     pw_up(ispc,i,k) = max(0.,sc_up(ispc,i,k)*(1.-exp(- CHEM_ADJ_AUTOC(ispc) * kc * (dz/w_upd)))*factor_temp(ispc,i,k))

                  !--formulation 3 - orignal GF conv_par
                  if(USE_TRACER_SCAVEN==3) then
                     !--- cloud liquid water tracer concentration
                     conc_mxr(ispc)  =  scav_eff* sc_up(ispc,i,k) !unit [kg(aq)/kg(air)]  for aerosol/smoke
                     !---   aqueous-phase concentration in rain water
                     pw_up(ispc,i,k) = conc_mxr(ispc)*pwo(i,k)/(1.e-8+qrco(i,k))
                  endif

                  !---(in cloud) total mixing ratio in gas and aqueous phases
                  sc_up(ispc,i,k) = sc_up(ispc,i,k) - pw_up(ispc,i,k)

                   !
               elseif(Hcts(ispc)%hstar>1.e-6) then ! tracer gas phase scavenging

                  !--- equilibrium tracer concentration - Henry's law
                  henry_coef=henry(ispc,tempco(i,k),rho(i,k))

                  if(USE_TRACER_SCAVEN==3) then
                     !--- cloud liquid water tracer concentration
                     conc_mxr(ispc) = (henry_coef*qrco(i,k) /(1.+henry_coef*qrco(i,k)) )* sc_up(ispc,i,k)
                     !
                     !---   aqueous-phase concentration in rain water
                     pw_up(ispc,i,k) = conc_mxr(ispc)*pwo(i,k)/(1.e-8+qrco(i,k))

                  else

                     !-- this the 'alpha' parameter in Eq 8 of Mari et al (2000 JGR) = X_aq/X_total
                     fliq = henry_coef*qrco(i,k) /(1.+henry_coef*qrco(i,k))

                     !---   aqueous-phase concentration in rain water
                     pw_up(ispc,i,k) = max(0.,sc_up(ispc,i,k)*(1.-exp(-fliq* CHEM_ADJ_AUTOC(ispc) *kc*dz/w_upd)))!*factor_temp(ispc,i,k))

                  endif

                  !---(in cloud) total mixing ratio in gas and aqueous phases
                  sc_up(ispc,i,k) = sc_up(ispc,i,k) - pw_up(ispc,i,k)

                   !
                   !---(in cloud)  mixing ratio in aqueous phase
                   !sc_up_aq(ispc,i,k) = conc_mxr(ispc) !if using set to zero at the begin.
               endif
            enddo
            !
            !-- total aerosol/gas in the rain water
            dp=100.*(po_cup(i,k)-po_cup(i,k+1))

            tot_pw_up_chem(:,i) = tot_pw_up_chem(:,i) + pw_up(:,i,k)*dp/g
         enddo loopk
      !
      !----- get back the in-cloud updraft gas-phase mixing ratio : sc_up(ispc,k)
      !          do k=start_level(i)+1,ktop(i)+1
      !            do ispc = 1,mtp
      !             sc_up(ispc,i,k) = sc_up(ispc,i,k) - sc_up_aq(ispc,i,k)
      !            enddo
      !          enddo
      enddo
   end subroutine get_incloud_sc_chem_up
   !---------------------------------------------------------------------------------------------------
   function henry(ispc,temp,rhoair) result(henry_coef)
      !--- calculate Henry's constant for solubility of gases into cloud water
      !--- inputs : ak0(ispc), dak(ispc),  hstar(ispc), dhr(ispc)
      implicit none
      integer, intent(in) :: ispc
      real   , intent(in) :: temp,rhoair
      real :: henry_coef
      real :: fct ,tcorr, corrh

      !--- define some constants!
      real, parameter:: rgas_ =  8.205e-2 ! atm M^-1 K^-1 ! 8.314 gas constant [J/(mol*K)]
      real, parameter:: avogad=  6.022e23! Avogadro constant [1/mol]
      real, parameter:: rhoH2O=  999.9668! density of water [kg/m3]
      real, parameter:: temp0 =    298.15! standard temperature [K]
      real, parameter:: temp0i= 1./298.15! inverse of standard temperature [K]
      real, parameter:: MWH2O =   18.02! molecular mass of water [kg/kmol]
      real, parameter:: MWAIR =   28.97! effective molecular mass of air [kg/kmol]
      real, parameter:: conv3 = avogad / 1.0e6!  [mol(g)/m3(air)]  to [molec(g)/cm3(air)]
      real, parameter:: conv4 = 100.   !  [m]   to [cm]
      real, parameter:: conv5 = 1000.     !  [m^3]      to [l]
      real, parameter:: conv7 = 1/conv5   !  [l]    to [m^3]
      real, parameter:: conv6 = 1. / 101325.  !  [Pa]      to [atm]
      real, parameter:: hplus = 1.175e-4     !  for cloud water. pH is asuumed to be 3.93: pH=3.93 =>hplus=10**(-pH)

      ! aqueous-phase concentrations XXXa [mol/m3(air)]!
      ! gas-phase concentrations XXXg [mol/m3(air)]!
      ! Henry constants XXXh for scavenging [mol/(l*atm)]!
      ! converted to [(mol(aq)/m3(aq))/(mol(g)/m3(air))], i.e. dimensionless!
      ! in equilibrium XXXa = XXXh * LWC * XXXg!
      tcorr = 1./temp - temp0i

      !-P. Colarco corrected the expression below
      !fct   = conv7 * rgas_ * temp ! - for henry_coef in units 1/m3
      fct   =         rgas_ * temp ! - for henry_coef dimensioless


      !-taking into account the acid dissociation constant
      ! ak=ak0*exp(dak*(1/t-1/298))
      corrh=1.+Hcts(ispc)%ak0 * exp(Hcts(ispc)%dak * tcorr)/hplus

      !-- for concentration in mol[specie]/mol[air] - Eq 5 in 'Compilation of Henry's law constants (version 4.0) for
      !-- water as solvent, R. Sander, ACP 2015'.
      henry_coef =  Hcts(ispc)%hstar* exp(Hcts(ispc)%dhr*tcorr) * fct * corrh

   end function henry
   !---------------------------------------------------------------------------------------------------
   subroutine get_incloud_sc_chem_dd(cumulus,FSCAV,mtp,se,se_cup,sc_dn,pw_dn ,pw_up,sc_up                          &
      ,tot_pw_up_chem,tot_pw_dn_chem                                                 &
      ,z_cup,rho,po_cup,qrcdo,pwdo,pwevo,edto,zdo,dd_massentro,dd_massdetro,pwavo,pwo &
      ,jmin,ierr,itf,ktf,its,ite, kts,kte)
      implicit none
      !-inputs
      integer                               ,intent (in)  :: itf,ktf, its,ite, kts,kte
      integer                               ,intent (in)  :: mtp
      integer, dimension (its:ite)          ,intent (in)  :: ierr ,jmin
      character *(*)                        ,intent (in)  :: cumulus
      real, dimension (mtp ,its:ite,kts:kte),intent (in)  :: se,se_cup,pw_up,sc_up
      real, dimension (mtp)                 ,intent (in)  :: FSCAV
      real, dimension (its:ite)             ,intent (in)  :: edto,pwavo,pwevo
      real, dimension (its:ite,kts:kte)     ,intent (in)  :: z_cup,rho,po_cup&
         ,qrcdo,pwdo,zdo,dd_massentro,dd_massdetro,pwo
      real, dimension (mtp ,its:ite        ),intent (in)  :: tot_pw_up_chem

      !-outputs
      real, dimension (mtp ,its:ite,kts:kte),intent (out) :: sc_dn,pw_dn
      real, dimension (mtp ,its:ite        ),intent (out) :: tot_pw_dn_chem

      !-locals
      real   , dimension (mtp) :: conc_mxr
      real :: x_add,dz,XZZ,XZD,XZE,denom, evaporate,pwdper,x1,frac_evap,dp,xkk
      integer :: i,k,ispc

      sc_dn          = 0.0
      pw_dn          = 0.0
      tot_pw_dn_chem = 0.0
      if(cumulus == 'shallow') return

      do i=its,itf
         if(ierr(i) /= 0) cycle

         !--- fration of the total rain that was evaporated
         frac_evap = - pwevo(i)/(1.e-16+pwavo(i))

         !--- scalar concentration in-cloud - downdraft

         !--- at k=jmim
         k=jmin(i)
         pwdper = pwdo(i,k)/(1.e-16+pwevo(i)) *frac_evap  ! > 0
         if(USE_TRACER_EVAP == 0 ) pwdper = 0.0

         dp= 100.*(po_cup(i,k)-po_cup(i,k+1))

         do ispc=1,mtp
            !--downdrafts will be initiate with a mixture of 50% environmental and in-cloud concentrations
            sc_dn(ispc,i,k) = se_cup(ispc,i,k)
                !sc_dn(ispc,i,k) = 0.9*se_cup(ispc,i,k)+0.1*sc_up(ispc,i,k)

            pw_dn(ispc,i,k) = - pwdper * tot_pw_up_chem(ispc,i)*g/dp
            sc_dn(ispc,i,k) = sc_dn(ispc,i,k) - pw_dn(ispc,i,k)
            tot_pw_dn_chem(ispc,i) = tot_pw_dn_chem(ispc,i) + pw_dn(ispc,i,k)*dp/g
         enddo
         !
         !--- calculate downdraft mass terms
         do k=jmin(i)-1,kts,-1
            XZZ=              zdo(i,k+1)
            XZD= 0.5*dd_massdetro(i,k  )
            XZE=     dd_massentro(i,k  )

            denom =  (XZZ-XZD+XZE)
            !-- transport + mixing
            if(denom > 0.) then
               sc_dn(:,i,k) = (sc_dn(:,i,k+1)*XZZ - sc_dn(:,i,k+1)*XZD + se(:,i,k)*XZE) / denom
            else
               sc_dn(:,i,k) = sc_dn(:,i,k+1)
            endif
            !
            !-- evaporation term
            if(USE_TRACER_EVAP == 0 )cycle

            dp= 100.*(po_cup(i,k)-po_cup(i,k+1))

            !-- fraction of evaporated precip per layer
            pwdper   = pwdo(i,k)/(1.e-16+pwevo(i))! > 0

            !-- fraction of the total precip that was actually evaporated at layer k
            pwdper   = pwdper * frac_evap

            !-- sanity check
            pwdper   = min(1.,max(pwdper,0.))

            do ispc=1,mtp
               !-- amount evaporated by the downdraft from the precipitation
               pw_dn(ispc,i,k) = - pwdper * tot_pw_up_chem (ispc,i)*g/dp ! < 0. => source term for the downdraft tracer concentration

                !if(ispc==1) print*,"pw=",pwdper,tot_pw_up_chem (ispc,i),pwevo(i)/pwavo(i),pwdo(i,k)/(1.e-16+pwo(i,k))

               !-- final tracer in the downdraft
               sc_dn(ispc,i,k) = sc_dn(ispc,i,k) - pw_dn(ispc,i,k) ! observe that -pw_dn is > 0.

               !-- total evaporated tracer
               tot_pw_dn_chem(ispc,i) = tot_pw_dn_chem(ispc,i) + pw_dn(ispc,i,k)*dp/g

                 !print*,"to=",k,tot_pw_dn_chem(ispc,i),pwdo(i,k)/(1.e-16+pwevo(i)),frac_evap,tot_pw_dn_chem(ispc,i)/tot_pw_up_chem (ispc,i)

            enddo
         enddo
         !
      enddo
   end subroutine get_incloud_sc_chem_dd
   !---------------------------------------------------------------------------------------------------
   subroutine interface_aerchem(mtp,itrcr,aer_chem_mech, cnames,qnames, fscav, fscav_int)
      implicit none
      integer, intent(in) :: mtp,ITRCR
      character(len=*),                intent(in)   :: AER_CHEM_MECH
      character(len=*),dimension(mtp) ,intent(in)   :: CNAMES,QNAMES
      real            ,dimension(ITRCR) ,intent(in) :: FSCAV

      real,    dimension(mtp) , intent(out)   :: FSCAV_INT
      !-local vars
      integer :: ispc, len_ACM, len_spc, irun=0
      character(len=100) :: TMP_AER_NAME
      ispc_CO             = 1
      CHEM_NAME_MASK      = 1
      CHEM_NAME_MASK_EVAP = 1
      CHEM_ADJ_AUTOC      = 1.0

      !- GOCART + PCHEM section
      if(AER_CHEM_MECH=='GOCART') then
         len_ACM=len(TRIM(AER_CHEM_MECH))
         do ispc=1,mtp
            FSCAV_INT(ispc)   = fscav(ispc)
            len_spc           = len(TRIM(cnames(ispc)))
            CHEM_name (ispc)  = TRIM(qnames(ispc))!(len_ACM+1:len_spc))
            if(TRIM(CHEM_name (ispc)) == 'CO') ispc_CO=ispc

            !-the tracers below are not being transported :
            if(trim(chem_name(ispc)) == "NCPI"    ) CHEM_NAME_MASK     (ispc) = 0
            if(trim(chem_name(ispc)) == "NCPL"    ) CHEM_NAME_MASK     (ispc) = 0
            if(trim(chem_name(ispc)) == "QGRAUPEL") CHEM_NAME_MASK     (ispc) = 0
            if(trim(chem_name(ispc)) == "QRAIN"   ) CHEM_NAME_MASK     (ispc) = 0
            if(trim(chem_name(ispc)) == "QSNOW"   ) CHEM_NAME_MASK     (ispc) = 0
            !if(trim(chem_name(ispc)) == "QW"      ) CHEM_NAME_MASK     (ispc) = 0
            !if(trim(chem_name(ispc)) == "QW"      ) CHEM_NAME_MASK_EVAP(ispc) = 0

            if(trim(chem_name(ispc)(1:3)) == "du0") CHEM_ADJ_AUTOC     (ispc) = 1.0
            if(trim(chem_name(ispc)(1:3)) == "ss0") CHEM_ADJ_AUTOC     (ispc) = 1.0
            if(trim(chem_name(ispc)(1:3)) == "OCp") CHEM_ADJ_AUTOC     (ispc) = 1.0
            if(trim(chem_name(ispc)(1:3)) == "BCp") CHEM_ADJ_AUTOC     (ispc) = 1.0
            if(trim(chem_name(ispc)(1:3)) == "SO4") CHEM_ADJ_AUTOC     (ispc) = 1.0
         enddo
         !-----------------------------------------------------------------------------------------
         !---temporary section to fill Henrys cts for N2O and CH4 of PCHEM chemical mechanism
         do ispc=1,mtp
            if (TRIM(CHEM_name (ispc)) == 'N2O' .or. &
               TRIM(CHEM_name (ispc)) == 'CH4'      &
               )then

               call get_HenrysLawCts(TRIM(CHEM_name (ispc)), &
                  Hcts(ispc)%hstar,Hcts(ispc)%dhr,Hcts(ispc)%ak0,Hcts(ispc)%dak)
            endif
         enddo
        !***   all tracers not having wet deposition _must_ have all Hcts reset to zero here.
        !.. WHERE( Hcts(:)%hstar == -1.) Hcts(:)%hstar = 0.
        !.. WHERE( Hcts(:)%dhr   == -1.) Hcts(:)%dhr   = 0.
        !.. WHERE( Hcts(:)%ak0   == -1.) Hcts(:)%ak0   = 0.
        !.. WHERE( Hcts(:)%dak   == -1.) Hcts(:)%dak   = 0.
      !-----------------------------------------------------------------------------------------
      endif
      return
      !IF( MAPL_AM_I_ROOT() .and. irun == 0)THEN
      irun = 1
      print*,"=========================================================================";call flush(6)
      print*," the following tracers will be transport by GF scheme                    ";call flush(6)

      write(10,*)"================= table of tracers in the GF conv transport ==================="
      write(10,*)" SPC,  CHEM_NAME,  FSCAV          - the four Henrys law cts  -   Transport Flag - kc adjust"
      do ispc=1,mtp
         write(10,121) ispc,trim(chem_name(ispc)), FSCAV_INT(ispc),Hcts(ispc)%hstar &
            ,Hcts(ispc)%dhr,Hcts(ispc)%ak0,Hcts(ispc)%dak, CHEM_NAME_MASK(ispc),CHEM_ADJ_AUTOC (ispc)
         if( CHEM_NAME_MASK     (ispc) == 1) then
            print*,"GF is doing transp and wet removal of: ",trim(chem_name(ispc))
            call flush(6)
         endif
      enddo
      print*,"=========================================================================";call flush(6)
121   format(1x,I4,A10,5F14.5,I4,F14.5)
      !ENDIF
   end subroutine interface_aerchem

   !-----------------------------------------------------------------------------------------

   subroutine fct1d3 (ktop,n,dt,z,tracr,massflx,trflx_in,del_out)

      ! --- modify a 1-D array of tracer fluxes for the purpose of maintaining
      ! --- monotonicity (including positive-definiteness) in the tracer field
      ! --- during tracer transport.

      ! --- the underlying transport equation is   (d tracr/dt) = - (d trflx/dz)
      ! --- where  dz = |z(k+1)-z(k)| (k=1,...,n) and  trflx = massflx * tracr

      ! --- note: tracr is carried in grid cells while z and fluxes are carried on
      ! --- interfaces. interface variables at index k are at grid location k-1/2.
      ! --- sign convention: mass fluxes are considered positive in +k direction.

      ! --- massflx and trflx_in  must be provided independently to allow the
      ! --- algorithm to generate an auxiliary low-order (diffusive) tracer flux
      ! --- as a stepping stone toward the final product trflx_out.

      implicit none
      integer,intent(in) :: n,ktop                        ! number of grid cells
      real   ,intent(in) :: dt                            ! transport time step
      real   ,intent(in) :: z(n+0)                        ! location of cell interfaces
      real   ,intent(in) :: tracr(n)                      ! the transported variable
      real   ,intent(in) :: massflx  (n+0)                ! mass flux across interfaces
      real   ,intent(in) :: trflx_in (n+0)                ! original tracer flux
      real   ,intent(out):: del_out  (n+0)                ! modified tracr flux
      real               :: trflx_out(n+0)                ! modified tracr flux
      integer k,km1,kp1
      logical :: NaN, error=.false., vrbos=.false.
      real dtovdz(n),trmax(n),trmin(n),flx_lo(n+0),antifx(n+0),clipped(n+0),  &
         soln_hi(n),totlin(n),totlout(n),soln_lo(n),clipin(n),clipout(n),arg
      real,parameter :: epsil=1.e-22           ! prevent division by zero
      real,parameter :: damp=1.                ! damper of antidff flux (1=no damping)

      logical, parameter :: hi_order = .false.

      NaN(arg) = .not. (arg.ge.0. .or. arg.le.0.)        ! NaN detector
      soln_lo(:)=0.
      antifx (:)=0.
      clipout(:)=0.
      flx_lo (:)=0.

      do k=1,ktop
         dtovdz(k)=.01*dt/abs(z(k+1)-z(k))                ! time step / grid spacing
      !     if (z(k).eq.z(k+1)) error=.true.
      end do
      if (vrbos .or. error) print '(a/(8es10.3))','(fct1d) dtovdz =',dtovdz(1:ktop)

      do k=2,ktop
         if (massflx(k) > 0.) then
            flx_lo(k)=massflx(k)*tracr(k-1)              ! low-order flux, upstream
         else
            flx_lo(k)=massflx(k)*tracr(k)                ! low-order flux, upstream
         endif
         antifx(k)=trflx_in(k)-flx_lo(k)                ! antidiffusive flux
      end do
      flx_lo(  1)   =trflx_in(  1)
      flx_lo(ktop+1)=trflx_in(ktop+1)
      antifx(  1)   =0.
      antifx(ktop+1)=0.
      ! --- clip low-ord fluxes to make sure they don't violate positive-definiteness
      do k=1,ktop
         totlout(k)=max(0.,flx_lo(k+1))-min(0.,flx_lo(k  ))         ! total flux out
         clipout(k)=min(1.,tracr(k)/max(epsil,totlout(k))/ (1.0001*dtovdz(k)))
      end do

      do k=2,ktop
         if (massflx(k).ge.0.)  then
            flx_lo(k)=flx_lo(k)*clipout(k-1)
         else
            flx_lo(k)=flx_lo(k)*clipout(k)
         endif
      end do
      if (massflx(1)     .lt.0.) flx_lo(1)     =flx_lo(1)     *clipout(1)
      if (massflx(ktop+1).gt.0.) flx_lo(ktop+1)=flx_lo(ktop+1)*clipout(ktop)

      ! --- a positive-definite low-order (diffusive) solution can now be  constructed

      do k=1,ktop
         soln_lo  (k)=tracr(k)-(flx_lo(k+1)-flx_lo(k))*dtovdz(k)        ! low-ord solutn
         del_out  (k)=-g*(flx_lo(k+1)-flx_lo(k))*dtovdz(k)/dt
      end do

      if(.not. hi_order) return

      soln_hi  (:)=0.
      clipin   (:)=0.
      trmin    (:)=0.
      trmax    (:)=0.
      clipped  (:)=0.
      trflx_out(:)=0.


      do k=1,ktop
         km1=max(1,k-1)
         kp1=min(n,k+1)
         trmax(k)=       max(soln_lo(km1),soln_lo(k),soln_lo(kp1),        &
            tracr  (km1),tracr  (k),tracr  (kp1))        ! upper bound
         trmin(k)=max(0.,min(soln_lo(km1),soln_lo(k),soln_lo(kp1),        &
            tracr  (km1),tracr  (k),tracr  (kp1)))       ! lower bound
      end do

      do k=1,ktop
         totlin (k)=max(0.,antifx(k  ))-min(0.,antifx(k+1))                ! total flux in
         totlout(k)=max(0.,antifx(k+1))-min(0.,antifx(k  ))                ! total flux out

         clipin (k)=min(damp,(trmax(k)-soln_lo(k))/max(epsil,totlin (k)) / (1.0001*dtovdz(k)))
         clipout(k)=min(damp,(soln_lo(k)-trmin(k))/max(epsil,totlout(k)) / (1.0001*dtovdz(k)))

         if (NaN(clipin (k))) print *,'(fct1d) error: clipin is NaN,  k=',k
         if (NaN(clipout(k))) print *,'(fct1d) error: clipout is NaN,  k=',k

         if (clipin(k).lt.0.) then
            print 100,'(fct1d) error: clipin < 0 at k =',k,                        &
               'clipin',clipin(k),'trmax',trmax(k),'soln_lo',soln_lo(k),        &
               'totlin',totlin(k),'dt/dz',dtovdz(k)
            error=.true.
         end if
         if (clipout(k).lt.0.) then
            print 100,'(fct1d) error: clipout < 0 at k =',k,                        &
               'clipout',clipout(k),'trmin',trmin(k),'soln_lo',soln_lo(k),        &
               'totlout',totlout(k),'dt/dz',dtovdz(k)
            error=.true.
         end if
100      format (a,i3/(4(a10,"=",es9.2)))
      end do

      do k=2,ktop
         if (antifx(k).gt.0.)  then
            clipped(k)=antifx(k)*min(clipout(k-1),clipin(k))
         else
            clipped(k)=antifx(k)*min(clipout(k),clipin(k-1))
         end if
         trflx_out(k)=flx_lo(k)+clipped(k)
         if (NaN(trflx_out(k)))  then
            print *,'(fct1d) error: trflx_out is NaN,  k=',k
            error=.true.
         end if
      end do

      trflx_out(     1)=trflx_in(     1)
      trflx_out(ktop+1)=trflx_in(ktop+1)
      do k=1,ktop
         soln_hi(k)=tracr(k)-(trflx_out(k+1)-trflx_out(k))*dtovdz(k)
         del_out(k) =     -g*(trflx_out(k+1)-trflx_out(k))*dtovdz(k)/dt
        !write(32,*)'3',k,soln_lo(k),soln_hi(k)
      end do

      if (vrbos .or. error) then
         do k=2,ktop
            write(32,99)k,                   &
               'tracr(k)', tracr(k),            &
               'flx_in(k)', trflx_in(k),        &
               'flx_in(k+1)', trflx_in(k+1),    &
               'flx_lo(k)', flx_lo(k),          &
               'flx_lo(k+1)', flx_lo(k+1),      &
               'soln_lo(k)', soln_lo(k),        &
               'trmin(k)', trmin(k),            &
               'trmax(k)', trmax(k),            &
               'totlin(k)', totlin(k),          &
               'totlout(k)', totlout(k),        &
               'clipin(k-1)', clipin(k-1),      &
               'clipin(k)', clipin(k),          &
               'clipout(k-1)', clipout(k-1),    &
               'clipout(k)', clipout(k),        &
               'antifx(k)', antifx(k),          &
               'antifx(k+1)', antifx(k+1),      &
               'clipped(k)', clipped(k),        &
               'clipped(k+1)', clipped(k+1),    &
               'flx_out(k)', trflx_out(k),      &
               'flx_out(k+1)', trflx_out(k+1),  &
               'dt/dz(k)', dtovdz(k),           &
               'final', tracr(k)-(trflx_out(k+1)-trflx_out(k))*dtovdz(k)
99          format ('(trc1d)   k =',i4/(3(a13,'=',es13.6)))
         end do
         if (error) stop '(fct1d error)'
      end if
   end subroutine fct1d3
   !---------------------------------------------------------------------------------------------------
   subroutine tridiag (m,a,b,c,f)
      !-- this routine solves the problem: aa*f(k-1,t+1) + bb*f(k,t+1) + cc*f(k+1,t+1) = dd
      !-- an updated "f" at time t+1 is the output
      implicit none
      integer, intent(in) :: m
      real, dimension(m), intent(inout) :: a,b,c
      real, dimension(m), intent(inout) :: f
      !--locals
      real, dimension(m) :: q
      integer :: k
      real :: p

      c(m)=0.
      q(1)=-c(1)/b(1)
      f(1)= f(1)/b(1)
      do k=2,m
         p  = 1./( b(k)+a(k)*q(k-1) )
         q(k) = -c(k)*p
         f(k) = p*(f(k) - a(k)*f(k-1))
      enddo
      do k=m-1,1,-1
         f(k) = f(k) +q(k)*f(k+1)
      enddo

   end subroutine tridiag
   !---------------------------------------------------------------------------------------------------
   subroutine bidiag (m,b,c,f)
      !-- this routine solves the problem:  bb*f(k,t+1) + cc*f(k+1,t+1) = dd
      !-- an updated "f" at time t+1 is the output
      implicit none
      integer, intent(in) :: m
      real, dimension(m), intent(inout) :: b,c
      real, dimension(m), intent(inout) :: f
      !--locals
      real, dimension(m) :: q
      integer :: k
      real :: p

      c(m)=0.
      q(1)=-c(1)/b(1)
      f(1)= f(1)/b(1)
      do k=2,m
         p  = 1./b(k)
         q(k) = -c(k)*p
         f(k) =  f(k)*p
      enddo
      do k=m-1,1,-1
         f(k) = f(k) +q(k)*f(k+1)
      enddo

   end subroutine bidiag
   !---------------------------------------------------------------------------------------------------
   subroutine cup_env_clev_chem(mtp,se_chem,se_cup_chem,ierr,itf,ktf,its,ite, kts,kte)

      implicit none
      !-inputs
      integer ,intent (in   )                   :: itf,ktf, its,ite, kts,kte
      integer ,intent (in   )                   :: mtp
      integer, dimension (its:ite),intent (in)  :: ierr
      real, dimension (mtp,its:ite,kts:kte),intent (in)  ::   se_chem
      !-outputs
      real, dimension (mtp,its:ite,kts:kte),intent (out) ::   se_cup_chem
      !-locals
      integer ::  i,k
      integer, parameter ::  clev_option = 2 !- use option 2

      !
      if(clev_option == 1) then
         !-- original version
         do i=its,itf
            if(ierr(i) /= 0) cycle
            do k=kts+1,ktf
               se_cup_chem(1:mtp,i,k)=0.5*(se_chem(1:mtp,i,k-1)+se_chem(1:mtp,i,k))
            enddo
            se_cup_chem(1:mtp,i,kts)=se_chem(1:mtp,i,kts)
            se_cup_chem(1:mtp,i,kte)=se_chem(1:mtp,i,ktf)
         enddo
      else
         !-- version 2: se_cup (k+1/2) = se(k) => smoother profiles
         do i=its,itf
            if(ierr(i) /= 0) cycle
            do k=kts,ktf
               se_cup_chem(1:mtp,i,k)=se_chem(1:mtp,i,k)
            enddo
         enddo
      endif

   end subroutine cup_env_clev_chem
   !---------------------------------------------------------------------------------------------------

   subroutine rain_evap_below_cloudbase(cumulus,itf,ktf, its,ite, kts,kte,ierr,kbcon,ktop,xmb,psur,xland&
      ,qo_cup,t_cup,po_cup,qes_cup,pwavo,edto,pwevo,pwo,pwdo&
      ,pre,prec_flx,evap_flx,outt,outq,outbuoy,evap_bcb)

      implicit none
      real, parameter :: alpha1=5.44e-4 & !1/sec
         ,alpha2=5.09e-3 & !unitless
         ,alpha3=0.5777  & !unitless
         ,c_conv=0.05      !conv fraction area, unitless

      character*(*)                      ,intent(in)    :: cumulus
      integer                            ,intent(in)    :: itf,ktf, its,ite, kts,kte
      integer, dimension(its:ite)        ,intent(in)    :: ierr,kbcon,ktop
      real,    dimension(its:ite)        ,intent(in)    :: psur,xland,pwavo,edto,pwevo,xmb
      real,    dimension(its:ite,kts:kte),intent(in)    :: po_cup,qo_cup,qes_cup,pwo,pwdo,t_cup
      real,    dimension(its:ite)        ,intent(inout) :: pre
      real,    dimension(its:ite,kts:kte),intent(inout) :: outt,outq,outbuoy,prec_flx,evap_flx

      real,    dimension(its:ite,kts:kte),intent(out)   :: evap_bcb

      !-- locals
      integer :: i,k
      real    :: RH_cr , del_t,del_q,dp,q_deficit, pqsat, temp_pre
      real    :: RH_cr_OCEAN,RH_cr_LAND
      real,    dimension(its:ite) :: tot_evap_bcb,eff_c_conv

      if(cumulus == 'shallow') then
         RH_cr_OCEAN   = 1.
         RH_cr_LAND    = 1.
         eff_c_conv(:) = min(0.2,max(xmb(:),c_conv))
      else
         RH_cr_OCEAN   = 0.95 !test 0.90
         RH_cr_LAND    = 0.90 
         eff_c_conv(:) = c_conv
      endif

      prec_flx     = 0.0
      evap_flx     = 0.0
      tot_evap_bcb = 0.0
      if(c0 < 1.e-6 ) return

      do i=its,itf

         if(ierr(i) /= 0) cycle

         !-- critical rel humidity  - check this, if the value is too small, not evapo will take place.
         RH_cr=RH_cr_OCEAN*xland(i)+RH_cr_LAND*(1.0-xland(i))

         !if(xland(i)  < 0.90 ) then !- over land
         !  RH_cr = RH_cr_LAND
         !else
         !  RH_cr = RH_cr_OCEAN
         !endif

         do k=ktop(i),kts,-1

            dp = 100.*(po_cup(i,k)-po_cup(i,k+1))

            !p_liq_ice(i,k) = fract_liq_f(tempco(i,k))

            !---rainfall evaporation below cloud base
            if(k <= kbcon(i)) then
               q_deficit = max(0.,(RH_cr*qes_cup(i,k) -qo_cup(i,k)))
               !pqsat=satur_spec_hum(t_cup(i,k),po_cup(i,k))

               !--units here: kg[water]/kg[air}/sec
               evap_bcb(i,k) = eff_c_conv(i) * alpha1 * q_deficit * &
                  ( sqrt(po_cup(i,k)/psur(i))/alpha2 * prec_flx(i,k+1)/eff_c_conv(i) )**alpha3

               !--units here: kg[water]/kg[air}/sec * kg[air]/m3 * m = kg[water]/m2/sec
               evap_bcb(i,k)= evap_bcb(i,k)*dp/g

            else

               evap_bcb(i,k)=0.0

            endif

            !-- before anything check if the evaporation already consumed all precipitation
            temp_pre = pre(i) - evap_bcb(i,k)
            if (temp_pre < 0.)  evap_bcb(i,k) =  pre(i)

            !-- get the net precitation flux after the local evaporation and downdraft
            prec_flx(i,k) = prec_flx(i,k+1) - evap_bcb(i,k) + xmb(i)*(pwo(i,k) + edto(i)*pwdo(i,k))
            prec_flx(i,k) = max(0.,prec_flx(i,k))

            evap_flx(i,k) = evap_flx(i,k+1) + evap_bcb(i,k) - xmb(i)*edto(i)*pwdo(i,k)
            evap_flx(i,k) = max(0.,evap_flx(i,k))

            tot_evap_bcb(i) = tot_evap_bcb(i)+evap_bcb(i,k)

            !-- feedback
            del_q =  evap_bcb(i,k)*g/dp          ! > 0., units: kg[water]/kg[air}/sec
            del_t = -evap_bcb(i,k)*g/dp*(xlv/cp) ! < 0., units: K/sec

            outq   (i,k) = outq   (i,k) + del_q
            outt   (i,k) = outt   (i,k) + del_t
            !--- comment out 17nov
            !outbuoy(i,k) = outbuoy(i,k) + cp*del_t+xlv*del_q

            pre(i) = pre(i) - evap_bcb(i,k)


            !--for future use (rain and snow precipitation fluxes)
            !prec_flx_rain(k) = prec_flx(i,k)*(1.-p_liq_ice(k))
            !prec_flx_snow(k) = prec_flx(i,k)*    p_liq_ice(k)


         enddo

         if(pre(i)<0.) then
            print*,"prec evap neg for cumulus=",pre(i),trim(cumulus)
            call flush(6)
            !stop '@subroutine rain_evap_below_cloudbase'
         endif

      enddo

   end subroutine rain_evap_below_cloudbase
   !---------------------------------------------------------------------------------------------------

   subroutine get_precip_fluxes(cumulus,klcl,kbcon,ktop,k22,ierr,xland,pre,xmb  &
      ,pwo,pwavo,edto,pwevo,pwdo,t_cup,tempco                   &
      ,prec_flx,evap_flx                                        &
      ,itf,ktf,its,ite, kts,kte)

      implicit none
      character *(*)            , intent (in) :: cumulus
      integer                    ,intent (in) :: itf,ktf,its,ite,kts,kte
      integer, dimension(its:ite),intent (in) :: kbcon,ktop,k22,klcl,ierr
      real,    dimension(its:ite),intent (in) :: xland,pwavo,pwevo,edto,pre,xmb
      real,    dimension(its:ite,kts:kte),intent (in)  :: pwo,pwdo,t_cup,tempco
      real,    dimension(its:ite,kts:kte),intent (out) :: prec_flx,evap_flx !-- units kg[water]/m2/s

      !-- locals
      integer :: i,k
      prec_flx = 0.0
      evap_flx = 0.0
      if(c0 < 1.e-6 ) return

      do i=its,itf
         if (ierr(i)  /= 0) cycle

         do k=ktop(i),kts,-1

            !--- precipitation flux (at 'cup' levels), units: kg[water]/m2/s
            prec_flx(i,k) = prec_flx(i,k+1) + xmb(i)*(pwo(i,k) + edto(i)*pwdo(i,k))
            prec_flx(i,k) = max(0.,prec_flx(i,k))

            !--- evaporation flux (at 'cup' levels), units: kg[water]/m2/s
            evap_flx(i,k) = evap_flx(i,k+1) - xmb(i)*edto(i)*pwdo(i,k)
            evap_flx(i,k) = max(0.,evap_flx(i,k))


            !
            !--for future use (rain and snow precipitation fluxes)
           !p_liq_ice(i,k) = fract_liq_f(tempco(i,k))
            !prec_flx_rain(k) = prec_flx(i,k)*(1.-p_liq_ice(k))
            !prec_flx_snow(k) = prec_flx(i,k)*    p_liq_ice(k)
         enddo

           !if(prec_flx   (i,kts) .ne. pre(i)) then
             !print*,"error=",100.*(prec_flx   (i,kts) - pre(i))/(1.e-16+pre(i)),pre(i),prec_flx   (i,kts)
             !STOP 'problem with water balance'
           !endif
      enddo
   end   subroutine get_precip_fluxes
   !------------------------------------------------------------------------------------
   real function satur_spec_hum(pt,press) result(pqsat)
      implicit none
      real   , intent(in ) :: pt,press ! Kelvin, hPa
      !real   , intent(out) :: pqsat  !saturation specific humidity kg/kg

      !---locals
      real :: zew,zqs,zcor,foealfcu,foeewmcu
      real, parameter ::                 &
         RD=287.06                          &
         ,RV=461.52               &
         ,RTT=273.16              &
         ,RETV=RV/RD-1.0             &
         ,R2ES=611.21*RD/RV             &
         ,R3LES=17.502            &
         ,R3IES=22.587            &
         ,R4LES=32.19                &
         ,R4IES=-0.7              &
         ,RTWAT=RTT               &
         ,RTICE=RTT-23.              &
         ,RTICECU=RTT-23.            &
         ,RTWAT_RTICE_R=1./(RTWAT-RTICE)     &
         ,RTWAT_RTICECU_R=1./(RTWAT-RTICECU)

      foealfcu = min(1.0,((max(rticecu,min(rtwat,pt))-rticecu)*rtwat_rticecu_r)**2)
      foeewmcu = r2es *(foealfcu *exp(r3les*(pt-rtt)/(pt-r4les))+&
         (1.0-foealfcu)*exp(r3ies*(pt-rtt)/(pt-r4ies)))

      zew  = foeewmcu
      zqs  = zew/(100.*press)
      if(1.0-retv*zqs > 0. )then
         zcor = 1.0/(1.0-retv*zqs)  ! divide by zero
         pqsat= zqs*zcor
      else
         pqsat= max_qsat
      endif

   end function satur_spec_hum
   !---------------------------------------------------------------------------------------------------
   subroutine get_jmin(cumulus,itf,ktf,its,ite, kts,kte,ierr,kdet,ktop,kbcon,jmin,ierrc  &
      ,beta,depth_min,heso_cup,zo_cup,melting_layer)

      implicit none
      character *(*)             ,intent (in)    :: cumulus
      real                       ,intent (in)    :: depth_min
      integer                    ,intent (in)    :: itf,ktf,its,ite,kts,kte
      integer, dimension(its:ite),intent (in)    :: ktop,kbcon
      real,    dimension(its:ite,kts:kte),intent (in)  ::  heso_cup,zo_cup,melting_layer

      integer, dimension(its:ite),intent (inout) :: ierr,jmin,kdet
      real                       ,intent (out)   :: beta
      character*128,              intent (out)   :: ierrc(its:ite)


      !-- locals
      integer :: i,k,jmini,ki
      real    :: dh,dz
      real,    dimension(its:ite,kts:kte)  ::  hcdo
      logical :: keep_going

      if(cumulus == 'deep') beta=0.05
      if(cumulus == 'mid' ) beta=0.02

      if(cumulus == 'shallow'  ) then
         beta    = 0.02
         jmin(:) = 0
         return
      endif

      do i=its,itf
         if(ierr(i) /= 0) cycle

         if(cumulus == 'deep' .and. melt_glac) jmin(i)=max(jmin(i),maxloc(melting_layer(i,:),1))

         !
         !--- check whether it would have buoyancy, if there where
         !--- no entrainment/detrainment
         !
         jmini = jmin(i)
         keep_going = .true.
         do while ( keep_going )
            keep_going = .false.
            if ( jmini - 1 .lt. kdet(i)   ) kdet(i) = jmini-1
            if ( jmini     .ge. ktop(i)-1 ) jmini = ktop(i) - 2
            ki = jmini
            hcdo(i,ki)=heso_cup(i,ki)
            dz=zo_cup(i,ki+1)-zo_cup(i,ki)
            dh=0.
            do k=ki-1,1,-1
               hcdo(i,k)=heso_cup(i,jmini)
               dz=zo_cup(i,k+1)-zo_cup(i,k)
               dh=dh+dz*(hcdo(i,k)-heso_cup(i,k))
               if(dh.gt.0.)then
                  jmini=jmini-1
                  if ( jmini .gt. 5 ) then
                     keep_going = .true.
                  else
                     ierr(i)=9
                     ierrc(i) = "could not find jmini9"
                     exit
                  endif
               endif
            enddo
         enddo
         jmin(i) = jmini
         if ( jmini .le. 5 ) then
            ierr(i)=4
            ierrc(i) = "could not find jmini4"
         endif
      enddo
      !
      ! - must have at least depth_min m between cloud convective base and cloud top.
      !
      do i=its,itf
         if(ierr(i) /= 0) cycle
         if ( jmin(i) - 1 .lt. kdet(i)) kdet(i) = jmin(i)-1
         if (-zo_cup(i,kbcon(i))+zo_cup(i,ktop(i)).lt.depth_min)then
            ierr(i)=6
            ierrc(i)="cloud depth very shallow"
         endif
      enddo

   end subroutine get_jmin
   !------------------------------------------------------------------------------------
   subroutine precip_cwv_factor(itf,ktf,its,ite,kts,kte,ierr,t,po,qo,po_cup,cumulus,p_cwv_ave)
      implicit none
      character *(*), intent (in)                        :: cumulus
      integer  ,intent (in )                             :: itf,ktf, its,ite, kts,kte
      integer  ,intent (in ), dimension(its:ite)         :: ierr
      real     ,intent (in ), dimension(its:ite,kts:kte) :: t,po,qo,po_cup
      real     ,intent (out), dimension(its:ite)         :: p_cwv_ave

      !--locals
      integer :: i,k
      real    :: dp, trash
      real    ,dimension(its:ite) :: w_col,w_ccrit,t_troposph
      real, parameter :: fpkup=0.8  !-- 90% of precip occurs above 80% of critical w

      p_cwv_ave = 0.0
      if(cumulus /= 'deep') return
      !
      !-- get the pickup of ensemble ave prec, following Neelin et al 2009.
      !
      do i=its,itf
         w_col     (i)= 0.
         w_ccrit   (i)= 0.
         t_troposph(i)= 0.
         if(ierr(i) /= 0) cycle
         trash=0.
         loopN:    do k=kts,ktf
            if(po(i,k) .lt. 200.) exit loopN

            dp=100.*(po_cup(i,k)-po_cup(i,k+1))
            trash=trash+dp/g

            w_col     (i)= w_col     (i) + qo(i,k)*dp/g ! unit mm
            t_troposph(i)= t_troposph(i) + t (i,k)*dp/g
         enddo loopN
         !--average temperature
         t_troposph(i) =   t_troposph(i)/(1.e-8+trash)! unit K
         !
         !--- wcrit given by Neelin et al 2009.
         w_ccrit(i) = max(0.,56.2 + 2.1*(t_troposph(i)-268.)) ! unit mm
         !
         !--- pickup (normalized by the factor 'a')
         !-- <p>=a[(w-w_c)/w_c]**beta, a=0.15, beta=0.23
         !
         p_cwv_ave(i) = (max(0.,w_col(i)-fpkup*w_ccrit(i))/(1.e-8+fpkup*w_ccrit(i)))**0.23
         p_cwv_ave(i) = max(0., min(1., p_cwv_ave(i)))

           !print*,"NEE=",i,w_col(i),t_troposph(i),w_ccrit(i),p_cwv_ave    (i)
           !print*,"=================================================="
      enddo
   end subroutine precip_cwv_factor
   !------------------------------------------------------------------------------------
   subroutine get_wetbulb(jmin,qo_cup,t_cup,po_cup ,q_wetbulb,t_wetbulb)

      implicit none
      integer ,intent (in   ) :: jmin
      real    ,intent (in   ) :: qo_cup,t_cup,po_cup
      real    ,intent (inout) :: q_wetbulb,t_wetbulb

      !---locals
      real ::  zqp, zcond, zcond1, zcor, zqsat
      real :: psp, pt , pq
      real :: z3es,   z4es, z5alcp, zaldcp
      real :: ptare, evap
      real :: foedelta,foeewmcu,foealfcu,foedemcu,foeldcpmcu

      real, parameter :: &
          RD=287.06                             &
         ,RV=461.52                             &
         ,RCPD=1004.71                          &
         ,RTT=273.16                            &
         ,RHOH2O=1000.                          &
         ,RLVTT=2.5008e+6                       &
         ,RLSTT=2.8345e+6                       &
         ,RETV=RV/RD-1.0                        &
         ,RLMLT=RLSTT-RLVTT                     &
         ,RCPV=4.*RV                            &
         ,R2ES=611.21*RD/RV                     &
         ,R3LES=17.502                          &
         ,R3IES=22.587                          &
         ,R4LES=32.19                           &
         ,R4IES=-0.7                            &
         ,R5LES=R3LES*(RTT-R4LES)               &
         ,R5IES=R3IES*(RTT-R4IES)               &
         ,R5ALVCP=R5LES*RLVTT/RCPD              &
         ,R5ALSCP=R5IES*RLSTT/RCPD              &
         ,RALVDCP=RLVTT/RCPD                    &
         ,RALSDCP=RLSTT/RCPD                    &
         ,RALFDCP=RLMLT/RCPD                    &
         ,RTWAT=RTT                             &
         ,RTBER=RTT-5.                          &
         ,RTBERCU=RTT-5.0                       &
         ,RTICE=RTT-23.                         &
         ,RTICECU=RTT-23.                       &
         ,RTWAT_RTICE_R=1./(RTWAT-RTICE)        &
         ,RTWAT_RTICECU_R=1./(RTWAT-RTICECU)    &
         ,RVTMP2=RCPV/RCPD-1.                   &
         ,ZQMAX=0.5

      !-- for testing
      !              PSP                   TEMP        Q                     ZCOND1
      ! input   85090.0000000000        289.140030372766     1.105078557441815E-002
      ! output  85090.0000000000        287.230570412846     1.181792062536557E-002 -2.761256206705639E-005
      ! PT  = 289.140030372766
      ! PQ  = 1.105078557441815E-002
      ! PSP = 85090.
      !----------------------

      !-- environmental values
      PT  = t_cup       ! K
      PQ  = qo_cup      ! kg/kg
      PSP = po_cup*100. ! hPa

      if (PT > RTT) then
         Z3ES=R3LES
         Z4ES=R4LES
         Z5ALCP=R5ALVCP
         ZALDCP=RALVDCP
      else
         Z3ES=R3IES
         Z4ES=R4IES
         Z5ALCP=R5ALSCP
         ZALDCP=RALSDCP
      endif

      !--- get wet bulb thermo properties --------------------------
      PTARE = PT
      ZQP    =1.0/PSP

      FOEALFCU = MIN(1.0,((MAX(RTICECU,MIN(RTWAT,PTARE))-RTICECU)*RTWAT_RTICECU_R)**2)
      FOEEWMCU = R2ES *(FOEALFCU *EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))+&
         (1.0-FOEALFCU)*EXP(R3IES*(PTARE-RTT)/(PTARE-R4IES)))
      ZQSAT=FOEEWMCU*ZQP

      ZQSAT=MIN(max_qsat,ZQSAT)
      ZCOR=1.0/(1.0-RETV  *ZQSAT)
      ZQSAT=ZQSAT*ZCOR

      FOEDEMCU =  FOEALFCU *R5ALVCP*(1.0/(PTARE-R4LES)**2)+&
         (1.0-FOEALFCU)*R5ALSCP*(1.0/(PTARE-R4IES)**2)

      ZCOND=(PQ-ZQSAT)/(1.0+ZQSAT*ZCOR*FOEDEMCU)

      ZCOND=MIN(ZCOND,0.0)

      FOELDCPMCU= FOEALFCU*RALVDCP+(1.0-FOEALFCU)*RALSDCP
      PT=PT+FOELDCPMCU*ZCOND

      PQ=PQ-ZCOND

      !--update PTARE
      PTARE = PT

      FOEALFCU = MIN(1.0,((MAX(RTICECU,MIN(RTWAT,PTARE))-RTICECU)*RTWAT_RTICECU_R)**2)
      FOEEWMCU = R2ES *(FOEALFCU *EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))+&
         (1.0-FOEALFCU)*EXP(R3IES*(PTARE-RTT)/(PTARE-R4IES)))
      ZQSAT=FOEEWMCU*ZQP

      ZQSAT=MIN(0.5,ZQSAT)
      ZCOR=1.0/(1.0-RETV  *ZQSAT)
      ZQSAT=ZQSAT*ZCOR

      FOEDEMCU =  FOEALFCU *R5ALVCP*(1.0/(PTARE-R4LES)**2)+&
         (1.0-FOEALFCU)*R5ALSCP*(1.0/(PTARE-R4IES)**2)
      ZCOND1=(PQ-ZQSAT)/(1.0+ZQSAT*ZCOR*FOEDEMCU)

      if(ZCOND == 0.0)ZCOND1=MIN(ZCOND1,0.0)
      FOELDCPMCU= FOEALFCU*RALVDCP+(1.0-FOEALFCU)*RALSDCP
      PT=PT+FOELDCPMCU*ZCOND1
      PQ=PQ-ZCOND1

      !-- set output --------------------------
      q_wetbulb =  PQ
      t_wetbulb =  PT
      evap      = -ZCOND1 != q_wetbulb-qo_cup, source for water vapor
   end subroutine get_wetbulb
   !------------------------------------------------------------------------------------
   subroutine cup_forcing_ens_3d_shal(itf,ktf,its,ite,kts,kte,dtime,ichoice     &
      ,ierrc,ierr,klcl,kpbl,kbcon,k22,ktop       &
      ,xmb,tsur,cape,h_sfc_flux,le_sfc_flux,zws  &
      ,po, hco, heo_cup,po_cup,t_cup,dhdt,rho    &
      ,xff_shal2d,xf_dicycle,tke_pbl,wlpool,xf_coldpool)

      implicit none
      integer                               ,intent (in) :: itf,ktf,its,ite, kts,kte,ichoice
      integer ,dimension (its:ite)          ,intent (in) :: klcl,kpbl,kbcon,k22,ktop
      real                                  ,intent (in) :: dtime
      real    ,dimension (its:ite)          ,intent (in) :: tsur,cape,h_sfc_flux,le_sfc_flux &
                                                           ,zws,tke_pbl,wlpool
      real    ,dimension (its:ite,kts:kte)  ,intent (in) :: po,hco,heo_cup,po_cup,t_cup,dhdt,rho
      integer ,dimension (its:ite)          ,intent (inout):: ierr
      character*128,dimension (its:ite)     ,intent (inout):: ierrc
      real    ,dimension (its:ite)               ,intent (inout)  :: xmb,xf_dicycle,xf_coldpool
      real    ,dimension (its:ite,shall_closures),intent (out)  :: xff_shal2d

      !---local vars
      real   ,dimension (its:ite)    :: xmbmax
      integer :: i,k,kbase
      real    :: blqe,trash,tcold,fin,fsum,efic,thot,dp
      real    ,dimension (shall_closures)  :: xff_shal
      !-- tuning numbers for the TKE-based closure for shallow convection   
      real,parameter :: k1 = 1.2, cloud_area = 0.15

      do i=its,itf
         xmb       (i)     = 0.
         xf_dicycle(i)     = 0.
         if(ierr(i) /= 0 ) cycle

         xmbmax(i)=100.*(po(i,kbcon(i))-po(i,kbcon(i)+1))/(g*dtime)

         !- limiting the mass flux at cloud base
         xmbmax(i)=min(xmbmaxshal,xmbmax(i))
         
         !- cloud base
         kbase=kbcon(i)
         !kbase=klcl(i)

         !--- closure from Grant (2001): ichoice = 1
         xff_shal(1)=.030*zws(i)*rho(i,kpbl(i))        
         xff_shal(2)=xff_shal(1)
         xff_shal(3)=xff_shal(1)
         
         !--- closure from the heat-engine principle : ichoice = 4
         !- Renno and Ingersoll(1996), Souza et al (1999)
         !- get the averaged environment temperature between cloud base
         !- and cloud top
         tcold=0.
         do k=kbase,ktop(i)
            dp   = po_cup(i,k)-po_cup(i,k+1)
            tcold= tcold + t_cup(i,k)*dp
         enddo
         tcold=tcold/(po_cup(i,kbase)-po_cup(i,ktop(i)+1))

         !-surface temperature
         thot=tsur(i)  ! + ztexec(i)
         !- thermodynamic eficiency
         !efic = max(0.05, (thot-tcold)/thot )
         efic = max(0.0, (thot-tcold)/thot )

         !- total heat flux from surface
         fin = max(0.0, h_sfc_flux(i)+le_sfc_flux(i))

         !--- mass flux at cloud base
         !if(cape(i) > 0.0 .and. h_sfc_flux(i) >0.0 ) then
         if(cape(i) > 0.0  ) then
            xff_shal(4) = efic * fin / cape(i)
         else
            xff_shal(4) = 0.0
         endif
         xff_shal(5)=xff_shal(4)
         xff_shal(6)=xff_shal(4)

         !--- closure from boundary layer QE (Raymond 1995): ichoice = 7
         blqe=0.
         trash=0.
         if(k22(i).lt.kpbl(i)+1)then
            do k=kts,kbase
               blqe=blqe+100.*dhdt(i,k)*(po_cup(i,k)-po_cup(i,k+1))/g
            enddo
            trash = max((hco(i,kbase)-heo_cup(i,kbase)),1.e1)
            xff_shal(7)=max(0.,blqe/trash)
         else
            xff_shal(7)=0.0
         endif
         xff_shal(8)= xff_shal(7)
         xff_shal(9)= xff_shal(7)

         !--- new closure based on the PBL TKE mean (Zheng et al, 2020 GRL): ichoice = 10
         !-- shallow cumulus active area is for now keept by 0.15 (Zheng 2021 p. commun.)
         !-- k1 is 'slope' of the curve between Wb x (TKE_PBL)**0.5
         !--        and varies between 1.2 (from lidar) to 1.6 (from WRF and SAM models)
         xff_shal(10) = cloud_area * rho(i,kbase) * k1 * sqrt(tke_pbl(i))
         xff_shal(11) = xff_shal(10)
         xff_shal(12) = xff_shal(10)

         !--- store all closures for later.
         xff_shal2d(i,:) = xff_shal(:)

      enddo
   end subroutine cup_forcing_ens_3d_shal

   !------------------------------------------------------------------------------------
   subroutine cup_up_lightning(itf,ktf,its,ite, kts,kte, ierr, kbcon,ktop,xland,cape &
      ,zo,zo_cup,t_cup,t,tempco,qrco,po_cup,rho,prec_flx     &
      ,lightn_dens)

      !=====================================================================================
      !- Lightning parameterization based on:
      !- "A Lightning Parameterization for the ECMWF Integrated Forecasting System"
      !-  P. Lopez, 2016 MWR
      !
      !- Coded/adapted to the GF scheme by Saulo Freitas (10-Aug-2019)
      !=====================================================================================
      implicit none
      integer                            ,intent(in)  :: itf,ktf, its,ite, kts,kte
      integer, dimension(its:ite)        ,intent(in)  :: ierr,kbcon,ktop
      real,    dimension(its:ite)        ,intent(in)  :: cape,xland
      real,    dimension(its:ite,kts:kte),intent(in)  :: po_cup,zo_cup,t_cup,t,tempco,zo &
                                                        ,qrco,rho,prec_flx

      real,    dimension(its:ite)        ,intent(out) :: lightn_dens ! lightning flash density
                                                                     ! rate (units: 1/km2/day)

      !-- locals
      real, parameter :: V_graup     = 3.0  ! m/s
      real, parameter :: V_snow      = 0.5  ! m/s
      real, parameter :: beta_land   = 0.70 ! 1
      real, parameter :: beta_ocean  = 0.45 ! 1
      real, parameter :: alpha       = 37.5 ! 1
      real, parameter :: t_initial   =  0.0 + 273.15 ! K
      real, parameter :: t_final     = -25. + 273.15 ! K

      integer :: i, k, k_initial, k_final
      real    :: Q_R, z_base,beta,prec_flx_fr,dz
      real,    dimension(kts:kte) :: p_liq_ice, q_graup,q_snow

      do i=its,itf
         lightn_dens(i) = 0.0
         if(ierr(i) /= 0) cycle

         beta= xland(i)*beta_ocean + (1.-xland(i))*beta_land

         q_graup(:) = 0.
         q_snow (:) = 0.

         do k=kts,ktop(i)

            p_liq_ice(k) = fract_liq_f(tempco(i,k))

            prec_flx_fr=   p_liq_ice(k)*prec_flx(i,k)/rho(i,k)

            q_graup(k) =      beta *prec_flx_fr/V_graup ! - graupel mixing ratio (kg/kg)
            q_snow (k) =  (1.-beta)*prec_flx_fr/V_snow  ! - snow    mixing ratio (kg/kg)

         enddo

         k_initial = minloc(abs(tempco(i,kbcon(i):ktop(i))-t_initial),1)+kbcon(i)-1
         k_final   = minloc(abs(tempco(i,kbcon(i):ktop(i))-t_final  ),1)+kbcon(i)-1

         Q_R = 0.0
         do k = k_initial, k_final
            dz  = zo(i,k)-zo(i,k-1)
            Q_R = Q_R + dz*rho(i,k)*(q_graup(k)*(qrco(i,k)+q_snow(k)))
           !print*,"qr=",q_r,tempco(i,k)-273.15,k,tempco(i,k)-t_initial
         enddo

         z_base = zo_cup(i,kbcon(i))/1000. ! km

         !---
         !--- lightning flash density (units: number of flashes/km2/day) - equation 5
         !--- (to compare with Lopez 2016's results, convert to per year: lightn_dens*365)
         !
         lightn_dens(i) = alpha * Q_R *sqrt (max(0.,cape(i))) * min(z_base,1.8)**2
        !
      enddo
   end subroutine cup_up_lightning

   !------------------------------------------------------------------------------------
   subroutine cup_up_rain(cumulus,klcl,kbcon,ktop,k22,ierr,xland        &
      ,zo_cup,qco,qrco,pwo,pwavo,po,p_cup,t_cup,tempco  &
      ,zuo,up_massentr,up_massdetr,vvel2d,rho           &
      ,qrr                                              &
      ,itf,ktf,its,ite, kts,kte)

      implicit none
      character *(*)            , intent (in) :: cumulus
      integer                    ,intent (in) :: itf,ktf,its,ite,kts,kte
      integer, dimension(its:ite),intent (in) :: kbcon,ktop,k22,klcl,ierr
      real,    dimension(its:ite),intent (in) :: xland,pwavo
      real,    dimension(its:ite,kts:kte),intent (in)  ::   &
         zo_cup,qco,qrco,pwo,po,p_cup,t_cup,zuo       &
         ,up_massentr,up_massdetr,vvel2d,tempco,rho

      !--for future use (rain water mixing ratio)
      real,    dimension(its:ite,kts:kte),intent (out) :: qrr      !-- units kg[water]/kg[air]

      !-- locals
      integer :: i,k
      real :: tmp
      integer, dimension(its:ite) :: start_level
      real :: dz,z1,zrnew,zc,zd,zint,z2,zrold,denom,fall_fact,wup, exp1,R_vr
      real,    dimension(kts:kte) :: prec_flx_rain,prec_flx_snow
      real,    dimension(kts:kte) :: pw,p_liq_ice ! - rainfall source
      real,    parameter :: rho1000mb = 1.2 , rhow = 1000., N_r = 0.1 ! cm^-3, rainfall drops number concen
      real,    parameter :: exp_KR    = 1./5. &! Kuo & Raymond 1983
         ,exp_SB    = 2./3.  ! Seifert & Beheng 2006 eq 27

      qrr = 0.
      if(c0 < 1.e-6 ) return

      !--- rain water mixing ratio
      do i=its,itf
         if (ierr(i)  /= 0) cycle

         do k=ktop(i),kts,-1

            p_liq_ice(k) = fract_liq_f(tempco(i,k))

            !--- transport + mixing
            denom = zuo(i,k+1)-.5*up_massdetr(i,k)+up_massentr(i,k)
            if(denom > 0.) then
               qrr(i,k) = (qrr(i,k+1)*zuo(i,k+1)-.5*up_massdetr(i,k)* qrr(i,k+1))/ denom
            else
               qrr(i,k) =  qrr(i,k+1)
            endif

            !--- rain source
            pw(k)= pwo(i,k)/(1.e-16+zuo(i,k))

            !-- rainfall sedimentation
            !-- Kuo & Raymond 1983
            !-- fallout of rain (21.18 * qrr^0.2 have m/s as units with qrr in kg/kg)
            !---                                half velocity for ice phase
            fall_fact = 21.18 *( p_liq_ice(k) + 0.5*(1.-p_liq_ice(k) ))
            exp1      = exp_KR
            !
            !-- Seifert & Beheng 2006 eq 27, units= m/s
            ! fall_fact = 159.*sqrt(rho1000mb/rho(i,k))*( p_liq_ice(k) + 0.5*(1.-p_liq_ice(k) ))
            ! exp1      = exp_SB

            !-- Kogan 2013
            R_vr = (4.*MAPL_PI*rhow/(3.*rho(i,k)))**(-1./3.) * (qrr(i,k) + pw(k))**(1./3.) * &
               (N_r)**(-1./3.)! cm, N_r is the rainfall drops number concentration, and not CCN or N_c

            R_vr = max(40., R_vr * 1.e-2 * 1.e+6 )   ! micrometer
            fall_fact = 1.e-2*(2.4*R_vr-62.0)        ! m/s
            exp1      = 0.

            wup      = min(15.,max(2.,vvel2d(i,k)))
            z2       = fall_fact/wup

            !--exact solution
            if(qrr(i,k) + pw(k)>0.) then
               !-- this is the sedimentation speed divided by W_up
               zd= z2* (qrr(i,k) + pw(k))**exp1
               zint= EXP(-zd)
               zrold=qrr(i,k)
               zc=pw(k)
               zrnew=zrold*zint+zc/zd*(1.-zint)
               zrnew=max(0.0,min(qrr(i,k) + pw(k),zrnew))
            else
               zrnew=0.
            endif
            qrr(i,k)=zrnew

          !--solution with rk3
          !z1= qrr(i,k)
          !z1= qrr(i,k) +  1./3.*(pw(k) - z2*(z1)**exp1 * z1)
          !z1= qrr(i,k) +  0.500*(pw(k) - z2*(z1)**exp1 * z1)
          !z1= qrr(i,k) +  1.000*(pw(k) - z2*(z1)**exp1 * z1)
          !qrr(i,k)=z1
          !---
          !--- solution 2
          !qrr(i,k)= qrr(i,k) + pw(k) - (fall_fact*(qrr(i,k) + pw(k))**exp1) / wup * qrr(i,k)

          !--- solution 3
          !tmp = 0.5*(qrr(i,k) + pw(k) + qrr(i,k+1) + pw(k+1))
          !qrr(i,k)= qrr(i,k) + pw(k) - z2*(tmp)**exp1 * 0.5*(qrr(i,k)+qrr(i,k+1))

          !print*,"rr3=",k,zo_cup(i,k), pwo(i,k)*1000.,qrr(i,k)*1000.,fall_fact*(qrr(i,k) + pw(k))**exp1
         enddo
      enddo

   end subroutine cup_up_rain

   !------------------------------------------------------------------------------------
   subroutine get_condensation(q_old,t_old,po_cup,q_new,t_new)

      !-- calculate condensation and adjust t and q accordingly
      implicit none
      real    ,intent (in   ) :: po_cup,q_old,t_old ! before condensation
      real    ,intent (inout) ::        q_new,t_new ! after  condensation

      !---locals
      real ::  zqp, zcond, zcond1, zcor, zqsat,zi,zl,zf
      real :: psp, pt , pq
      real :: z3es,   z4es, z5alcp, zaldcp
      real :: ptare, cond
      real :: foeewmcu,foealfcu,foedemcu,foeldcpmcu

      real, parameter :: &
          RD=287.06                             &
         ,RV=461.52                             &
         ,RCPD=1004.71                          &
         ,RTT=273.16                            &
         ,RHOH2O=1000.                          &
         ,RLVTT=2.5008e+6                       &
         ,RLSTT=2.8345e+6                       &
         ,RETV=RV/RD-1.0                        &
         ,RLMLT=RLSTT-RLVTT                     &
         ,RCPV=4.*RV                            &
         ,R2ES=611.21*RD/RV                     &
         ,R3LES=17.502                          &
         ,R3IES=22.587                          &
         ,R4LES=32.19                           &
         ,R4IES=-0.7                            &
         ,R5LES=R3LES*(RTT-R4LES)               &
         ,R5IES=R3IES*(RTT-R4IES)               &
         ,R5ALVCP=R5LES*RLVTT/RCPD              &
         ,R5ALSCP=R5IES*RLSTT/RCPD              &
         ,RALVDCP=RLVTT/RCPD                    &
         ,RALSDCP=RLSTT/RCPD                    &
         ,RALFDCP=RLMLT/RCPD                    &
         ,RTWAT=RTT                             &
         ,RTBER=RTT-5.                          &
         ,RTBERCU=RTT-5.0                       &
         ,RTICE=RTT-23.                         &
         ,RTICECU=RTT-23.                       &
         ,RTWAT_RTICE_R=1./(RTWAT-RTICE)        &
         ,RTWAT_RTICECU_R=1./(RTWAT-RTICECU)    &
         ,RVTMP2=RCPV/RCPD-1.                   &
         ,ZQMAX=0.5

      !----------------------
      !-- for testing
      !----------------------
      !                          PSP (hPA)            TEMP(K)              Q(kg/kg)                ZCOND1(kg/kg)
      ! input          1   98020.0000000000        295.163188640152     1.745679200989956E-002
      ! output         1   98020.0000000000        295.513490789916     1.731605637618801E-002 -9.779916453243843E-007
      !----------------------
      ! input       157   85090.0000000000        288.089188935407     1.399404805052166E-002
      ! output      157   85090.0000000000        289.294751760460     1.350970717999820E-002 -1.146268822756454E-005
      !----------------------
      ! PT  = 288.089188935407
      ! PQ  = 1.399404805052166E-002
      ! PSP = 85090
      !----------------------

      !-- initial values
      PT  = t_old       ! K
      PQ  = q_old       ! kg/kg
      PSP = po_cup*100. ! hPa


      !--- get condensation in moist ascent --------------------------
      PTARE = PT
      ZQP    =1.0/PSP

      ZL=1.0/(PT-R4LES)
      ZI=1.0/(PT-R4IES)

      FOEALFCU = MIN(1.0,((MAX(RTICECU,MIN(RTWAT,PTARE))-RTICECU)*RTWAT_RTICECU_R)**2)
      ZQSAT=R2ES *(FOEALFCU *EXP(R3LES*(PTARE-RTT)*ZL)+&
         (1.0-FOEALFCU)*EXP(R3IES*(PTARE-RTT)*ZI))

      ZQSAT=ZQSAT*ZQP
      ZQSAT=MIN(0.5,ZQSAT)
      ZCOR=1.0-RETV*ZQSAT

      ZF=FOEALFCU*R5ALVCP*ZL**2 + (1.0-FOEALFCU)*R5ALSCP*ZI**2
      ZCOND=(PQ*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*ZF)

      if(ZCOND > 0.0)then
         FOELDCPMCU= FOEALFCU*RALVDCP+(1.0-FOEALFCU)*RALSDCP
         PT=PT+FOELDCPMCU*ZCOND
         PTARE = PT
         PQ=PQ-ZCOND

         ZL=1.0/(PT-R4LES)
         ZI=1.0/(PT-R4IES)


         FOEALFCU = MIN(1.0,((MAX(RTICECU,MIN(RTWAT,PTARE))-RTICECU)*RTWAT_RTICECU_R)**2)
         ZQSAT=R2ES *(FOEALFCU* EXP(R3LES*(PT-RTT)*ZL)+&
            (1.0-FOEALFCU)*EXP(R3IES*(PT-RTT)*ZI))

         ZQSAT=ZQSAT*ZQP
         ZQSAT=ZQSAT-0.5*(ABS(0.5-ZQSAT)-(0.5-ZQSAT))


         ZCOR=1.0-RETV*ZQSAT
         ZF=FOEALFCU*R5ALVCP*ZL**2 + (1.0-FOEALFCU)*R5ALSCP*ZI**2

         ZCOND1=(PQ*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*ZF)
         if(ZCOND ==  0.0)ZCOND1=0.0
         FOELDCPMCU= FOEALFCU*RALVDCP+(1.0-FOEALFCU)*RALSDCP
         PT=PT+FOELDCPMCU*ZCOND1
         PQ=PQ-ZCOND1
      endif

      !-- FINAL --------------------------
      q_new =  PQ
      t_new =  PT
      cond  = -ZCOND1 != q_old-qnew, source for the liquid water
   end subroutine get_condensation

   !------------------------------------------------------------------------------------
   subroutine get_interp(q_old,t_old,po_cup,q_new,t_new)
      implicit none
      real    ,intent (in   ) :: po_cup ! original
      real    ,intent (inout) :: q_old,t_old,q_new,t_new ! extrapolated

      !---locals
      real ::  zqp, zcond1, zcor, zqsat
      real ::  psp, pt , pq, ptare
      real ::  FOEALFCU, FOEEWMCU,FOEDEMCU,FOELDCPMCU
      real, parameter :: &
          RD=287.06                             &
         ,RV=461.52                             &
         ,RCPD=1004.71                          &
         ,RTT=273.16                            &
         ,RHOH2O=1000.                          &
         ,RLVTT=2.5008e+6                       &
         ,RLSTT=2.8345e+6                       &
         ,RETV=RV/RD-1.0                        &
         ,RLMLT=RLSTT-RLVTT                     &
         ,RCPV=4.*RV                            &
         ,R2ES=611.21*RD/RV                     &
         ,R3LES=17.502                          &
         ,R3IES=22.587                          &
         ,R4LES=32.19                           &
         ,R4IES=-0.7                            &
         ,R5LES=R3LES*(RTT-R4LES)               &
         ,R5IES=R3IES*(RTT-R4IES)               &
         ,R5ALVCP=R5LES*RLVTT/RCPD              &
         ,R5ALSCP=R5IES*RLSTT/RCPD              &
         ,RALVDCP=RLVTT/RCPD                    &
         ,RALSDCP=RLSTT/RCPD                    &
         ,RALFDCP=RLMLT/RCPD                    &
         ,RTWAT=RTT                             &
         ,RTBER=RTT-5.                          &
         ,RTBERCU=RTT-5.0                       &
         ,RTICE=RTT-23.                         &
         ,RTICECU=RTT-23.                       &
         ,RTWAT_RTICE_R=1./(RTWAT-RTICE)        &
         ,RTWAT_RTICECU_R=1./(RTWAT-RTICECU)    &
         ,RVTMP2=RCPV/RCPD-1.                   &
         ,ZQMAX=0.5

      integer :: i

      PT  = t_old       ! K
      PQ  = q_old       ! kg/kg
      PSP = po_cup*100. ! hPa

      !-- for testing
      !              PSP                   TEMP        Q                     ZCOND1
      ! input    27940.0000000000        236.604976804749       3.220181796223121E-004
      ! output   27940.0000000000        236.361132108860       4.084506812610067E-004
      !  PT  = 236.604976804749      ! K
      !  PQ  = 3.220181796223121E-004       ! kg/kg
      !  PSP = 27940. ! hPa
      !----------------------
      !print*,"1",PSP,PT,PQ

      ZQP   =1.0/PSP
      do i=1,2
         PTARE = PT

         FOEALFCU = MIN(1.0,((MAX(RTICECU,MIN(RTWAT,PTARE))-RTICECU)*RTWAT_RTICECU_R)**2)
         FOEEWMCU = R2ES *(FOEALFCU *EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))+&
            (1.0-FOEALFCU)*EXP(R3IES*(PTARE-RTT)/(PTARE-R4IES)))
         ZQSAT=FOEEWMCU*ZQP

         !    if(1.0-RETV  *ZQSAT == 0.) then
         !
         !      print*,"ZQSAT=",ZQP,FOEEWMCU,q_old,t_old,po_cup,q_new,t_new
         !3.5491847E-02   46.36052      0.5000000       249.8219
         !  0.2817549      0.5000000       249.8219
         !      call flush(6)
         !      stop 3333
         !    endif

         ZCOR=1.0/(1.0-RETV  *ZQSAT)
         ZQSAT=ZQSAT*ZCOR

         FOEDEMCU =  FOEALFCU *R5ALVCP*(1.0/(PTARE-R4LES)**2)+&
            (1.0-FOEALFCU)*R5ALSCP*(1.0/(PTARE-R4IES)**2)


         ZCOND1=(PQ-ZQSAT)/(1.0+ZQSAT*ZCOR*FOEDEMCU)

         FOELDCPMCU= FOEALFCU*RALVDCP+(1.0-FOEALFCU)*RALSDCP
         PT=PT+FOELDCPMCU*ZCOND1
         PQ=PQ-ZCOND1
      enddo
      !-- FINAL --------------------------
      q_new =  PQ
      t_new =  PT
      !print*,"2",PSP,PT,PQ
      !print*,"E",100*(PT-236.361132108860)/236.361132108860,100*(PQ-4.084506812610067E-004)/4.084506812610067E-004
   end subroutine get_interp
   !------------------------------------------------------------------------------------
   real function fract_liq_f(temp2) ! temp2 in Kelvin, fraction between 0 and 1.
      implicit none
      real,intent(in)  :: temp2 ! K
      real             :: temp,ptc
      real, parameter  :: max_temp = 46. !Celsius
      select case(FRAC_MODIS)

         case (1)
            temp = temp2-273.16 !Celsius
            temp = min(max_temp,max(-max_temp,temp))
            ptc  = 7.6725 + 1.0118*temp + 0.1422*temp**2 + &
               0.0106*temp**3 + 3.39e-4 * temp**4    + &
               3.95e-6 * temp**5
            fract_liq_f = 1./(1.+exp(-ptc))

         !WMP skew ice fraction for deep convective clouds
         !       fract_liq_f = fract_liq_f**4
         !WMP
         case default
            fract_liq_f =  min(1., (max(0.,(temp2-t_ice))/(t_0-t_ice))**2)

      end select

   end function
   !------------------------------------------------------------------------------------

   subroutine prepare_temp_pertubations(kts,kte,ktf,its,ite,itf,jts,jte,jtf,dt,xland,topo,zm &
      ,temp,rqvften,rthblten,rthften,Tpert_h,Tpert_v&
      ,AA1_CIN,AA1_BL)!<<<<<<<<<<
      implicit none
      integer ,intent(in) :: kts,kte,ktf,its,ite,itf,jts,jte,jtf
      real    ,intent(in) :: dt

      real    ,dimension(its:ite,jts:jte)        , intent(in) :: xland   !flag < 1 para land
                                                                         !flag  =1 para water
      real    ,dimension(its:ite,jts:jte)        , intent(in) :: topo

      real    ,dimension(kts:kte,its:ite,jts:jte), intent(in) ::        &
         zm       &
         ,rthften  &
         ,rqvften  &
         ,rthblten &
         ,temp

      real    ,dimension(kts:kte,its:ite,jts:jte), intent(out)::Tpert_h,Tpert_v


      real   ,dimension(its:ite,jts:jte)  ,intent(inout)  ::AA1_CIN,AA1_BL!<<<<<<<<<<

      !-- local vars
      integer :: i,j,k,kr,i1,i2,i3,j1,j2,j3
      real    :: aveh_qmax,aveh_qmin,avev_qmax,avev_qmin,coef_h,coef_v,dz1,dz2,dz3
      real, dimension(kts:kte)                  :: avev_adv_q,avev_t
      real, dimension(kts:kte,its:ite,jts:jte)  :: ave_T,ave_Q_adv
      integer, dimension(its:ite,jts:jte)       :: is_ocean

      !--avoid the trigger over the land areas
      is_ocean(:,:) = 0
      do j = jts,jtf
         do i= its,itf
            if(xland(i,j) > 0.999 .and. topo(i,j) < 10.) is_ocean(i,j) = 1
         enddo
      enddo

      !-reset/initialization
      Tpert_h=0.
      Tpert_v=0.
      ave_T=0.
      ave_Q_adv=0.
      avev_adv_q=0.
      avev_adv_q =0.
      avev_t=0.

      !-- calculate 9-point average of moisture advection and temperature using halo (Horizontal)
      !-- in the future these lines must be moved to the dynamics, because there is not "halo" in physics.
      do j = jts,jtf
         j1 = max (j-1,jts)
         j2 = j
         j3 = min(j+1,jtf)
         do i= its,itf
            i1 = max (i-1,its)
            i2 = i
            i3 = min(i+1,itf)
            if(is_ocean(i,j) == 0) cycle  ! do it only over ocean regions

            do k = kts, ktf
               kr = k

               !--think about using temp _OR_ temp_new(kr,i,j)= temp(kr,i,j) + (rthblten(kr,i,j)+rthften(kr,i,j))*dt
                  ave_T (k,i,j) = ( temp(kr,i1,j1) + temp(kr,i1,j2) + temp(kr,i1,j3) + & ! row 1
                  temp(kr,i2,j1) + temp(kr,i2,j2) + temp(kr,i2,j3) + & ! row 2
                  temp(kr,i3,j1) + temp(kr,i3,j2) + temp(kr,i3,j3)   & ! row 3
                  ) / 9.
               !-- advection forcing of q
                  ave_Q_adv (k,i,j) = ( rqvften(kr,i1,j1) + rqvften(kr,i1,j2) + rqvften(kr,i1,j3) + & ! row 1
                  rqvften(kr,i2,j1) + rqvften(kr,i2,j2) + rqvften(kr,i2,j3) + & ! row 2
                  rqvften(kr,i3,j1) + rqvften(kr,i3,j2) + rqvften(kr,i3,j3)   & ! row 3
                  ) / 9.

            enddo
         enddo
      enddo

      !--search for max/min moisture advection in 9-point range, calculate horizontal T-perturbation (tpert_h)

      do j = jts,jtf
         j1 = max (j-1,jts)
         j2 = j
         j3 = min(j+1,jtf)
         do i= its,itf
            i1 = max (i-1,its)
            i2 = i
            i3 = min(i+1,itf)

            if(is_ocean(i,j) == 0) cycle

            do k = kts, ktf
               kr = k
               aveh_qmax = maxval( ave_Q_adv(k,i1:i3,j1:j3) )
               aveh_qmin = minval( ave_Q_adv(k,i1:i3,j1:j3) )

               if(aveh_qmax > aveh_qmin )then
                  coef_h = (ave_Q_adv(k,i,j) - aveh_qmin)/(aveh_qmax-aveh_qmin)
               else
                  coef_h = 0.
               endif
               coef_h=min(coef_h,1.0)
               coef_h=max(coef_h,0.0)

               !--think about using temp _OR_ temp_new(kr,i,j)= temp(kr,i,j) + (rthblten(kr,i,j)+rthften(kr,i,j))*dt
               Tpert_h(k,i,j)=coef_h *(temp(kr,i,j)-ave_T(k,i,j))

            enddo
         enddo
      enddo


      !--search for max/min moisture advection in 3-point vertical range, calculate vertical T-perturbation (tpert_v)
      do j = jts,jtf
         do i= its,itf
            if(is_ocean(i,j) == 0) cycle

            do k = kts+1,ktf-2
               kr=k
               dz1 = zm(kr  ,i,j) -  zm(kr-1,i,j)
               dz2 = zm(kr+1,i,j) -  zm(kr  ,i,j)
               dz3 = zm(kr+2,i,j) -  zm(kr+1,i,j)

               !--think about using temp _OR_ temp_ne(kr,i,j)w=temp(kr,i,j) + (rthblten(kr,i,j)+rthften(kr,i,j))*dt
               avev_t     (k) = ( dz1*temp   (kr-1,i,j) + dz2*temp   (kr,i,j)+ dz3*temp   (kr+1,i,j) )/ (dz1+dz2+dz3)
               avev_adv_q (k) = ( dz1*rqvften(kr-1,i,j) + dz2*rqvften(kr,i,j)+ dz3*rqvften(kr+1,i,j) )/ (dz1+dz2+dz3)

            enddo
            avev_t     (kts)       =  avev_t     (kts+1)
            avev_adv_q (kts)       =  avev_adv_q (kts+1)
            avev_t     (ktf-1:ktf) =  avev_t     (ktf-2)
            avev_adv_q (ktf-1:ktf) =  avev_adv_q (ktf-2)


            do k = kts+1,ktf-2
               kr=k
               avev_qmax = maxval( avev_adv_q (k-1:k+1) )
               avev_qmin = minval( avev_adv_q (k-1:k+1) )
               if(avev_qmax > avev_qmin) then
                  coef_v = (avev_adv_q(k)-avev_qmin)/(avev_qmax-avev_qmin)
               else
                  coef_v = 0.
               endif

               !--think about using temp _OR_ temp_new(kr,i,j)=temp(kr,i,j) + (rthblten(kr,i,j)+rthften(kr,i,j))*dt
               Tpert_v(k,i,j)=coef_v*( temp(kr,i,j)-avev_t(k) )

            enddo
            Tpert_v(kts,i,j)=  Tpert_v(kts+1,i,j)
            Tpert_v(ktf,i,j)=  Tpert_v(ktf-1,i,j)

         enddo
      enddo

      !--avoid the trigger over the land areas
      do j = jts,jtf
         do i= its,itf
            if(is_ocean(i,j) == 1) cycle  ! ocean areas
            do k=kts,ktf
               Tpert_v(k,i,j)= 0.! Tpert_v(k,i,j) * xland(i,j)
               Tpert_h(k,i,j)= 0.! Tpert_h(k,i,j) * xland(i,j)
            enddo

             !print*,"TperH=",i,j,whoami_all,maxval(Tpert_h(:,i,j)),minval(Tpert_h(:,i,j))
             !print*,"TperV=",i,j,whoami_all,maxval(Tpert_V(:,i,j)),minval(Tpert_V(:,i,j))
         enddo
      enddo

      !-----
      return
      !-----

      !---- check balance
      do j = jts,jtf
         do i= its,itf
            dz2=0.
            AA1_BL (i,j)=0.
            AA1_CIN(i,j)=0.

            do k = kts,13 !-- 13 ~ 900 hPa
               dz1 = zm(k+1,i,j) -  zm(k,i,j)
               dz2 = dz2 + dz1
               AA1_BL  (i,j)=AA1_BL  (i,j)+dz1*Tpert_v(k,i,j)
               AA1_cin (i,j)=AA1_cin (i,j)+dz1*Tpert_h(k,i,j)
            enddo
            AA1_BL  (i,j)=AA1_BL  (i,j)/(dz2+1.e-16)
            AA1_cin (i,j)=AA1_cin (i,j)/(dz2+1.e-16)
         enddo
      enddo


       ! New trigger function
       ! Na vertical, fa\E7a isto apenas dentro da camada de 60 hPa entorno do k22.
       ! IF(trigger.eq.2) then
       !         DTLCL=amax1(tpart_h1D(KLCL)+tpart_v1D(KLCL),0.0)
       !  ENDIF

   end subroutine prepare_temp_pertubations
   !------------------------------------------------------------------------------------
   subroutine SOUND(part,cumulus,int_time,dtime,ens4,itf,ktf,its,ite, kts,kte,xlats,xlons,jcol,whoami_all  &
      ,z ,qes ,he ,hes ,t ,q ,po,z1 ,psur,zo,qeso,heo,heso,tn,qo,us,vs ,omeg,xz    &
      ,h_sfc_flux,le_sfc_flux,tsur, dx,stochastic_sig,zws,ztexec,zqexec, xland     &
      ,kpbl,k22,klcl,kbcon,ktop,aa0,aa1,sig,xaa0,hkb,xmb,pre,edto                  &
      ,zo_cup,dhdt,rho,zuo,zdo,up_massentro,up_massdetro,outt, outq,outqc,outu,outv)

      implicit none
      character*(*) ,intent(in)    :: cumulus
      integer, intent(in) ::ens4, itf,ktf,its,ite, kts,kte,jcol,whoami_all,part
      real   , intent(in) ::int_time,dtime
      real,    dimension (its:ite,kts:kte) :: z ,qes ,he ,hes ,t ,q ,po,zo,qeso,heo,heso,tn,qo &
         ,us,vs ,xz,zo_cup,dhdt,rho &
         ,zuo,zdo,up_massentro,up_massdetro,outt, outq,outqc,outu,outv
      real,    dimension (its:ite,kts:kte,1:ens4) :: omeg

      integer, dimension (its:ite) ::kpbl,k22,klcl,kbcon,ktop
      real,    dimension (its:ite) ::h_sfc_flux,le_sfc_flux,tsur, dx,stochastic_sig,zws,ztexec &
         ,zqexec,xlats,xlons,xland,z1,psur
      real,    dimension (its:ite) ::aa0,aa1,xaa0,hkb,xmb,pre,edto,sig

      !---locals
      integer :: i,k,X_kte,X_i,X_jcol,X_k,X_WHOAMI_ALL
      real :: X_time
      character(200) :: lixo
      real,    dimension (its:ite) :: X_stochastic_sig,X_xland
      real, parameter :: LATSND=-10., LONSND= 301., DELTX=0.2
      !      real, parameter :: LATSND= -8.72, LONSND= 186.6, DELTX=0.2

      if(trim(rundata)=="NONE") then

         if( mod(int_time,3600.) < dtime ) then
            open(15, file="dataLXXX.dat_"//trim(cumulus),status='unknown',position="APPEND")

            if(part == 1) then

               do i=its,itf

                  if( xlats(i) > LATSND-DELTX .and. xlats(i) < LATSND+DELTX ) then
                     if( xlons(i) > LONSND-DELTX .and. xlons(i) < LONSND+DELTX ) then

                        print*,"==============================================="
                        print*,"00>", i,jcol,xlats(i),xlons(i),whoami_all,int_time/3600.
                        call flush(6)

                        write(15,*) "====begin====="
                        write(15,*) "i,jcol,xlats(i),xlons(i),int_time/3600."
                        write(15,*)  i,jcol,xlats(i),xlons(i),int_time/3600.

                        write(15,*) "kte,z1(i),psur(i),tsur(i),xland(i)"
                        write(15,*)  kte,z1(i),psur(i),tsur(i),xland(i)

                        write(15,*) "h_sfc_flux(i),le_sfc_flux(i),ztexec(i),zqexec(i)"
                        write(15,*)  h_sfc_flux(i),le_sfc_flux(i),ztexec(i),zqexec(i)

                        write(15,*) "stochastic_sig(i), dx(i),zws(i),kpbl(i)"
                        write(15,*)  stochastic_sig(i), dx(i),zws(i),kpbl(i)

                        write(15,*) "=>k zo po t tn-t q qo-q us vs qes he hes qeso-qes heo-he heso-hes dhdt omeg"
                        do k = kts,kte
                           write(15,100)   k,zo(i,k), po(i,k) ,t  (i,k) ,tn(i,k)-t(i,k), q(i,k) ,qo(i,k)-q(i,k)  &
                              ,us (i,k), vs(i,k) ,qes(i,k) ,he(i,k) ,hes(i,k) &
                              ,qeso(i,k)- qes(i,k) &
                              ,heo (i,k)- he (i,k) &
                              ,heso(i,k)- hes(i,k) &
                              ,dhdt(i,k), omeg(i,k,1:ens4)
                        enddo


                     endif
                  endif
               enddo

            else

               do i=its,itf
                  if(xlats(i)  > LATSND-DELTX .and. xlats(i) < LATSND+DELTX ) then
                     if(xlons(i)  > LONSND-DELTX .and. xlons(i) < LONSND+DELTX ) then

                        write(15,*) "====outputs======="
                        write(15,*) "L=",i,jcol,xlats(i),xlons(i),whoami_all
                        write(15,*) "A=",aa0(i),aa1(i),xaa0(i),sig(i)
                        write(15,*) "K=",k22(i),klcl(i),kpbl(i),kbcon(i),ktop(i)
                        write(15,*) "Z=",zo_cup(i,k22(i))-z1(i),zo_cup(i,klcl(i))-z1(i),zo_cup(i,kpbl(i))-z1(i)&
                           ,zo_cup(i,kbcon(i))-z1(i),zo_cup(i,ktop(i))-z1(i)
                        write(15,*) "H=",hkb(i)/cp,edto(i)
                        write(15,*) "T=",maxval(outt(i,1:ktop(i)))*86400.,maxval(outq(i,1:ktop(I)))*86400.*1000.,&
                           minval(outt(i,1:ktop(i)))*86400.,minval(outq(i,1:ktop(I)))*86400.*1000.
                        write(15,*) "P=",xmb(i)*1000.,'g/m2/s',3600*pre(i),'mm/h'
                        if(xmb(i)>0.0) then
                           write(15,*) "=> k zo po zuo,zdo,up_massentro,up_massdetro,outt, outq,outqc,outu,outv"
                           do k = kts,kte
                              write(15,101) k,zo(i,k), po(i,k) &
                                 ,zuo(i,k),zdo(i,k),up_massentro(i,k),up_massdetro(i,k),outt(i,k)*86400. &
                                 ,outq(i,k)*86400.*1000.,outqc(i,k)*86400.*1000.,outu(i,k)*86400.,outv(i,k)*86400.

                           enddo
                        endif
                        write(15,*) "=====end=========="
                     endif
                  endif
               enddo
            endif
            close(15)
         endif
      else

         if(part == 1) then

            open(15, file=trim(rundata),status='old')
            i=1
            read(15,*) lixo
            read(15,*) lixo
            read(15,*) X_i,X_jcol,xlats(i),xlons(i),X_TIME
            read(15,*) lixo
            read(15,*) X_kte,z1(i),psur(i),tsur(i),X_xland(i)
                  !-- check
            if(X_kte .ne. kte) stop " X_kte .ne. kte "
            read(15,*) lixo
            read(15,*) h_sfc_flux(i),le_sfc_flux(i),ztexec(i),zqexec(i)
            read(15,*) lixo
            read(15,*) X_stochastic_sig(i), dx(i),zws(i),kpbl(i)
            read(15,*) lixo
            do k = kts,kte
               read(15,100)  X_k,  zo(i,k), po(i,k) ,t  (i,k) ,tn  (i,k) ,q  (i,k) ,qo  (i,k), us  (i,k) ,vs(i,k) &
                  , qes(i,k) ,he(i,k) ,hes(i,k) ,qeso(i,k) ,heo(i,k) ,heso(i,k), dhdt(i,k) &
                  ,omeg(i,k,1:ens4)
            enddo
            close(15)
            !---settings
            tn    (i,:) = t  (i,:) + tn   (i,:) ! input is delta(T)
            qo    (i,:) = q  (i,:) + qo   (i,:) ! input is delta(Q)
            qeso  (i,:) = qes(i,:) + qeso (i,:) ! input is delta(Q)
            heo   (i,:) = he (i,:) + heo  (i,:) ! input is delta(H)
            heso  (i,:) = hes(i,:) + heso (i,:) ! input is delta(HO)
            xz    (i,:) = zo (i,:)
            z     (i,:) = zo (i,:)
            rho   (i,:) = 1.e2*po(i,:)/(rgas*t(i,:))

         else

            do i=its,itf
               if(xlats(i)  > LATSND-DELTX .and. xlats(i) < LATSND+DELTX ) then
                  if(xlons(i)  > LONSND-DELTX .and. xlons(i) < LONSND+DELTX ) then

                     print*, "====outputs======="
                     print*, "A=",aa0(i),aa1(i),xaa0(i),sig(i)
                     print*, "K=",k22(i),klcl(i),kpbl(i),kbcon(i),ktop(i)
                     print*, "Z=",zo_cup(i,k22(i))-z1(i),zo_cup(i,klcl(i))-z1(i),zo_cup(i,kpbl(i))-z1(i)&
                        ,zo_cup(i,kbcon(i))-z1(i),zo_cup(i,ktop(i))-z1(i)
                     print*, "H=",hkb(i)/cp,edto(i)
                     print*, "T=",maxval(outt(i,1:ktop(i)))*86400.,maxval(outq(i,1:ktop(I)))*86400.*1000.,&
                        minval(outt(i,1:ktop(i)))*86400.,minval(outq(i,1:ktop(I)))*86400.*1000.
                     print*, "P=",xmb(i)*1000.,'g/m2/s',3600*pre(i),'mm/h'
                     if(xmb(i)>0.0) then
                        print*, "=> k zo po zuo,zdo,up_massentro,up_massdetro,outt, outq,outqc,outu,outv"
                        do k = kts,kte
                           write(*,101) k,zo(i,k), po(i,k) &
                              ,zuo(i,k),zdo(i,k),up_massentro(i,k),up_massdetro(i,k),outt(i,k)*86400. &
                              ,outq(i,k)*86400.*1000.,outqc(i,k)*86400.*1000.,outu(i,k)*86400.,outv(i,k)*86400.
                        enddo
                     endif
                  endif
               endif
            enddo
         endif
      endif
100   format(1x,i4,16e16.8)
101   format(1x,i4,11e16.8)

   end subroutine SOUND
   !------------------------------------------------------------------------------------
   subroutine cloud_dissipation(cumulus,itf,ktf, its,ite, kts,kte,ierr,kbcon,ktop,dtime,xmb,xland&
      ,qo_cup,qeso_cup,po_cup,outt,outq,outqc,zuo,vvel2d,rho_hydr       &
      ,qrco,sig,tempco,qco,tn_cup,heso_cup,zo)

      implicit none
      character*(*)                      ,intent(in)    :: cumulus
      integer                            ,intent(in)    :: itf,ktf, its,ite, kts,kte
      real                               ,intent(in)    :: dtime
      integer, dimension(its:ite)        ,intent(in)    :: ierr,kbcon,ktop
      real,    dimension(its:ite)        ,intent(in)    :: xmb,xland,sig
      real,    dimension(its:ite,kts:kte),intent(in)    :: po_cup,qo_cup,qeso_cup,zuo,vvel2d,rho_hydr&
         ,tempco,qco,tn_cup,heso_cup,zo
      real,    dimension(its:ite,kts:kte),intent(inout) :: outt,outq,outqc,qrco

      !-- locals
      integer :: i,k
      real    ::  del_t,del_q,dp,frh
      real    :: qrc_diss,fractional_area,outqc_diss,outq_mix,outt_diss,outt_mix,tempx,qvx
      real, parameter :: cloud_lifetime= 1800.

      integer, parameter :: versionx = 2
      do i=its,itf

         if(ierr(i) /= 0) cycle

         do k=ktop(i),kbcon(i),-1

            !--- cloud liq/ice remained in the convection plume
            qrc_diss = max(0., qrco(i,k) - outqc(i,k) * dtime)

            !dp  = 100.*(po_cup(i,k)-po_cup(i,k+1))

            !--- get relative humidity
            frh = 0. !min(qo_cup(i,k)/qeso_cup(i,k),1.)

            !--- estimation of the fractional area
            fractional_area = (xmb(i)/sig(i)) * zuo(i,k) / (rho_hydr(i,k)*vvel2d(i,k))

            !--- source of enviroment moistening/cooling due to the 'remained' cloud dissipation into it.
            outqc_diss = ( qrc_diss * (1.-frh) ) / cloud_lifetime

            if(versionx==1 .or. COUPL_MPHYSICS .eqv. .false.) then

               outt_diss  = -outqc_diss*(xlv/cp) !--- cooling

               !--- source of enviroment moistening/warming due to the 'remained' in-cloud water vapor mixing into it.
               !  qvx   = qco   (i,k)
               !  tempx = tempco(i,k)
               qvx   = qeso_cup(i,k)
               tempx =(heso_cup(i,k)-g*zo(i,k)-xlv*qeso_cup(i,k))/cp

               outq_mix = ( qvx   - qo_cup(i,k) ) / cloud_lifetime

               outt_mix = ( tempx - tn_cup(i,k) ) / cloud_lifetime

               !-- feedback
               del_q = (outqc_diss + outq_mix) * use_cloud_dissipation * fractional_area ! units: kg[water]/kg[air}/sec
               del_t = (outt_diss  + outt_mix) * use_cloud_dissipation * fractional_area ! units: K/sec

               outq (i,k) = outq (i,k) + del_q
               outt (i,k) = outt (i,k) + del_t

            else

               outqc(i,k) = outqc(i,k) + outqc_diss* fractional_area * use_cloud_dissipation

            endif

             !print*,"diss2=",k,real(outqc_diss*86400.*1000),real(sqrt(1.-sig(i)),4),real( fractional_area*100.,4)

            qrco    (i,k) = max(0., qrco(i,k) - outqc_diss * use_cloud_dissipation * fractional_area *dtime)
           !if(qrco (i,k) <0.) print*,"qrc<0",trim(cumulus),qrco(i,k)

         enddo
      enddo

   end subroutine cloud_dissipation
   !------------------------------------------------------------------------------------
   subroutine GF_convpar_init (mynum)
      implicit none
      integer,intent(in) :: mynum
      integer :: nlunit = 4
      character (len = 64) :: fn_nml = 'GF_ConvPar_nml'
      logical :: exists
      namelist /GF_NML/ icumulus_gf, closure_choice,use_scale_dep,dicycle                                  &
         ,use_tracer_transp, use_tracer_scaven,use_flux_form, use_tracer_evap,downdraft,use_fct  &
         ,use_rebcb, vert_discr, satur_calc, clev_grid,apply_sub_mp, alp1         &
         ,sgs_w_timescale, lightning_diag,autoconv, bc_meth,overshoot,use_wetbulb &
         ,c1,c0_deep, qrc_crit,lambau_deep,lambau_shdn,c0_mid                     &
         ,cum_max_edt_land  ,cum_max_edt_ocean, cum_hei_down_land                 &
         ,cum_hei_down_ocean,cum_hei_updf_land, cum_hei_updf_ocean                &
         ,cum_entr_rate ,tau_deep,tau_mid                                         &
         ,zero_diff ,use_momentum_transp ,moist_trigger,frac_modis                &
         ,cum_use_excess,cum_ave_layer,adv_trigger,use_smooth_prof, evap_fix      &
         ,use_cloud_dissipation,use_smooth_tend,use_gustiness, use_random_num     &
         ,dcape_threshold,beta_sh,c0_shal,use_linear_subcl_mf,liq_ice_number_conc &
         ,alpha_adv_tuning,cap_maxs,sig_factor,cum_fadj_massflx,lcl_trigger       &
         ,rh_dicycle,cum_t_star, convection_tracer, tau_ocea_cp, tau_land_cp      &
         ,use_memory, add_coldpool_prop ,mx_buoy1, mx_buoy2,max_tq_tend,cum_zuform&
         ,add_coldpool_clos,add_coldpool_diff

      inquire (file = trim (fn_nml), exist = exists)
      if (.not. exists) then
         write (6, *) 'GF_convpar_nml :: namelist file: ', trim (fn_nml), ' does not exist'
         stop 31415
      else
         open (nlunit,file=fn_nml,status='old',form='formatted')
         read (nlunit,nml=GF_NML)
         close(nlunit)
      endif
      if(mynum==1) then
         !- print the namelist
         print*,"           "
         print*,"------------- GF ConvPar namelist -------------"
         print*,"!---- the main controls"
         print*, 'icumulus_gf        ' , icumulus_gf
         print*, 'cum_entr           ' , real(cum_entr_rate      ,4)
         print*, 'closure_choice     ' , closure_choice
         print*, 'use_scale_dep      ' , use_scale_dep
         print*, 'sig_factor         ' , real(sig_factor         ,4)
         print*, 'dicycle            ' , dicycle
         print*, 't_star             ' , real(cum_t_star         ,4)
         print*, 'rh_dicycle         ' , rh_dicycle
         print*, 'alpha_adv_tuning   ' , real(alpha_adv_tuning   ,4)
         print*, 'cap_maxs           ' , real(cap_maxs           ,4)
         print*, 'moist_trigger      ' , moist_trigger
         print*, 'adv_trigger        ' , adv_trigger
         print*, 'lcl_trigger        ' , lcl_trigger
         print*, 'dcape_threshold    ' , real(dcape_threshold    ,4)
         print*, 'tau_deep,tau_mid   ' , real(tau_deep,4),real(tau_mid,4)
         print*, 'sgs_w_timescale    ' , sgs_w_timescale
         print*, 'convection_tracer  ' , convection_tracer
         print*, 'add_coldpool_prop  ' , add_coldpool_prop
         print*, 'add_coldpool_clos  ' , add_coldpool_clos
         print*, 'add_coldpool_diff  ' , add_coldpool_diff
         print*, 'tau_ocea_cp        ' , tau_ocea_cp
         print*, 'tau_land_cp        ' , tau_land_cp
         print*, 'mx_buoy1 - kJ/kg   ' , mx_buoy1*1.e-3
         print*, 'mx_buoy2 - kJ/kg   ' , mx_buoy2*1.e-3
         print*, 'use_memory         ' , use_memory
         
         print*,'!--- controls rainfall evaporation'
         print*, 'use_rebcb          ' , use_rebcb
         print*, 'downdraft          ' , downdraft
         print*, 'max_edt_land       ' , real(cum_max_edt_land   ,4)
         print*, 'max_edt_ocean      ' , real(cum_max_edt_ocean  ,4)

         print*,'!---- boundary condition specification'
         print*, 'bc_meth            ' , bc_meth
         print*, 'cum_use_excess     ' , cum_use_excess
         print*, 'cum_ave_layer      ' , real(cum_ave_layer      ,4)

         print*,'!---- for mass flux profiles - (deep ,shallow ,congestus)'
         print*, 'cum_zuform         ' , cum_zuform
         print*, 'hei_down_land      ' , real(cum_hei_down_land  ,4)
         print*, 'hei_down_ocean     ' , real(cum_hei_down_ocean ,4)
         print*, 'hei_updf_land      ' , real(cum_hei_updf_land  ,4)
         print*, 'hei_updf_ocean     ' , real(cum_hei_updf_ocean ,4)
         print*, 'beta_sh            ' , real(beta_sh            ,4)
         print*, 'use_linear_subcl_mf' , use_linear_subcl_mf
         print*, 'use_smooth_prof    ' , use_smooth_prof
         print*, 'use_smooth_tend    ' , use_smooth_tend
         print*, 'use_random_num     ' , use_random_num

         print*,'!---- the cloud microphysics'
         print*, 'autoconv           ' , autoconv
         print*, 'c0_deep            ' , real(c0_deep            ,4)
         print*, 'c0_mid             ' , real(c0_mid             ,4)
         print*, 'c0_shal            ' , real(c0_shal            ,4)
         print*, 'c1                 ' , real(c1                 ,4)
         print*, 'qrc_crit           ' , real(qrc_crit           ,4)

         print*, '!--- for momentum transport'
         print*, 'use_momentum_trans ' , use_momentum_transp
         print*, 'lambau_deep        ' , real(lambau_deep        ,4)
         print*, 'lambau_shdn        ' , real(lambau_shdn        ,4)

         print*, '!--- for tracer transport'
         print*, 'use_tracer_transp  ' , use_tracer_transp
         print*, 'use_tracer_scaven  ' , use_tracer_scaven
         print*, 'use_flux_form      ' , use_flux_form
         print*, 'use_fct            ' , use_fct
         print*, 'use_tracer_evap    ' , use_tracer_evap
         print*, 'apply_sub_mp       ' , apply_sub_mp
         print*, 'alp1               ' , real(alp1               ,4)

         print*,'!---- couplings w/ other parameterizations'
         print*, 'lightning_diag     ' , lightning_diag
         print*, 'overshoot          ' , real(overshoot          ,4)
         print*, 'liq_ice_number_conc' , liq_ice_number_conc
         print*, 'use_gustiness      ' , use_gustiness

         print*, '!----misc controls'
         print*, 'zero_diff          ' , zero_diff
         print*, 'frac_modis         ' , frac_modis
         print*, 'evap_fix           ' , evap_fix
         print*, 'use_cloud_dissipat ' , real(use_cloud_dissipation,4)
         print*, 'use_wetbulb        ' , use_wetbulb
         print*, 'clev_grid          ' , clev_grid
         print*, 'vert_discr         ' , vert_discr
         print*, 'satur_calc         ' , satur_calc
         print*, 'max_tq_tend        ' , real(max_tq_tend,4)
         print*, 'cum_fadj_massflx   ' , real(cum_fadj_massflx     ,4)
         print*,"========================================================================"
         call flush(6)
      endif
   end subroutine GF_convpar_init
   !------------------------------------------------------------------------------------
   subroutine get_liq_ice_number_conc(itf,ktf,its,ite, kts,kte,ierr,ktop    &
      ,dtime,rho,outqc,tempco,outnliq,outnice)

      implicit none
      integer,   intent (in )  :: itf,ktf,its,ite,kts,kte
      real,      intent (in )  :: dtime

      integer, dimension (its:ite)         ,intent (in )  :: ierr,ktop
      real,    dimension (its:ite,kts:kte) ,intent (in )  :: outqc,tempco,rho
      real,    dimension (its:ite,kts:kte) ,intent (out)  :: outnliq,outnice

      integer :: i,k
      real    :: fr,tqliq,tqice,dtinv

      real,   dimension (its:ite,kts:kte) :: nwfa   ! in the future set this as NCPL
      real,   dimension (its:ite,kts:kte) :: nifa   ! in the future set this as NCPI


      nwfa(:,:) =  99.e7  ! in the future set this as NCPL
      nifa(:,:) =  0.     ! in the future set this as NCPI
      dtinv    = 1./dtime
      do i=its,itf
         if(ierr(i) /= 0) cycle

         do k=kts,ktop(i)+1

            fr    = fract_liq_f(tempco(i,k))
            tqliq = dtime * outqc(i,k)* rho(i,k) * fr
            tqice = dtime * outqc(i,k)* rho(i,k) * (1.-fr)

            outnice(i,k) = max(0.0,  make_IceNumber    (tqice, tempco(i,k))/rho(i,k))
            outnliq(i,k) = max(0.0,  make_DropletNumber(tqliq, nwfa  (i,k))/rho(i,k))

         enddo
         !-- convert in tendencies
         outnice = outnice * dtinv ! unit [1/s]
         outnliq = outnliq * dtinv ! unit [1/s]
      !--- for update
      ! nwfa =nwfa + outnliq*dtime
      ! nifa =nifa + outnice*dtime

      enddo

   end subroutine get_liq_ice_number_conc
   !+---+-----------------------------------------------------------------+
   !+---+-----------------------------------------------------------------+
   !----- module_mp_thompson_make_number_concentrations
   !- Developed by H. Barnes @ NOAA/OAR/ESRL/GSL Earth Prediction Advancement Division
   !-----------------------------------------------------------------------
   !      Q_ice              is cloud ice mixing ratio, units of kg/m3
   !      Q_cloud            is cloud water mixing ratio, units of kg/m3
   !      Q_rain             is rain mixing ratio, units of kg/m3
   !      temp               is air temperature in Kelvin
   !      make_IceNumber     is cloud droplet number mixing ratio, units of number per m3
   !      make_DropletNumber is rain number mixing ratio, units of number per kg of m3
   !      make_RainNumber    is rain number mixing ratio, units of number per kg of m3
   !      qnwfa              is number of water-friendly aerosols in number per kg

   !+---+-----------------------------------------------------------------+
   !+---+-----------------------------------------------------------------+
   elemental real function make_IceNumber (Q_ice, temp)

      implicit none
      real, parameter:: ice_density = 890.0
      real, parameter:: pi = 3.1415926536
      real, intent(in):: q_ice, temp
      integer :: idx_rei
      real :: corr, reice, deice
      double precision :: lambda

      !+---+-----------------------------------------------------------------+
      !..Table of lookup values of radiative effective radius of ice crystals
      !.. as a function of Temperature from -94C to 0C.  Taken from WRF RRTMG
      !.. radiation code where it is attributed to Jon Egill Kristjansson
      !.. and coauthors.
      !+---+-----------------------------------------------------------------+

      real, dimension(95), parameter:: retab = (/                       &
         5.92779, 6.26422, 6.61973, 6.99539, 7.39234,                   &
         7.81177, 8.25496, 8.72323, 9.21800, 9.74075, 10.2930,          &
         10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319,          &
         15.0351, 15.8799, 16.7674, 17.6986, 18.6744, 19.6955,          &
         20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125,          &
         27.7895, 28.6450, 29.4167, 30.1088, 30.7306, 31.2943,          &
         31.8151, 32.3077, 32.7870, 33.2657, 33.7540, 34.2601,          &
         34.7892, 35.3442, 35.9255, 36.5316, 37.1602, 37.8078,          &
         38.4720, 39.1508, 39.8442, 40.5552, 41.2912, 42.0635,          &
         42.8876, 43.7863, 44.7853, 45.9170, 47.2165, 48.7221,          &
         50.4710, 52.4980, 54.8315, 57.4898, 60.4785, 63.7898,          &
         65.5604, 71.2885, 75.4113, 79.7368, 84.2351, 88.8833,          &
         93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424,          &
         124.954, 130.630, 136.457, 142.446, 148.608, 154.956,          &
         161.503, 168.262, 175.248, 182.473, 189.952, 197.699,          &
         205.728, 214.055, 222.694, 231.661, 240.971, 250.639 /)

      if (Q_ice == 0) then
         make_IceNumber = 0
         return
      end if

      !+---+-----------------------------------------------------------------+
      !..From the model 3D temperature field, subtract 179K for which
      !.. index value of retab as a start.  Value of corr is for
      !.. interpolating between neighboring values in the table.
      !+---+-----------------------------------------------------------------+

      idx_rei = int(temp-179.)
      idx_rei = min(max(idx_rei,1),94)
      corr = temp - int(temp)
      reice = retab(idx_rei)*(1.-corr) + retab(idx_rei+1)*corr
      deice = 2.*reice * 1.e-6

      !+---+-----------------------------------------------------------------+
      !..Now we have the final radiative effective size of ice (as function
      !.. of temperature only).  This size represents 3rd moment divided by
      !.. second moment of the ice size distribution, so we can compute a
      !.. number concentration from the mean size and mass mixing ratio.
      !.. The mean (radiative effective) diameter is 3./Slope for an inverse
      !.. exponential size distribution.  So, starting with slope, work
      !.. backwords to get number concentration.
      !+---+-----------------------------------------------------------------+

      lambda = 3.0 / deice
      make_IceNumber = Q_ice * lambda*lambda*lambda / (PI*Ice_density)

      !+---+-----------------------------------------------------------------+
      !..Example1: Common ice size coming from Thompson scheme is about 30 microns.
      !.. An example ice mixing ratio could be 0.001 g/kg for a temperature of -50C.
      !.. Remember to convert both into MKS units.  This gives N_ice=357652 per kg.
      !..Example2: Lower in atmosphere at T=-10C matching ~162 microns in retab,
      !.. and assuming we have 0.1 g/kg mixing ratio, then N_ice=28122 per kg,
      !.. which is 28 crystals per liter of air if the air density is 1.0.
      !+---+-----------------------------------------------------------------+

      return
   end function make_IceNumber

   !+---+-----------------------------------------------------------------+
   !+---+-----------------------------------------------------------------+

   elemental real function make_DropletNumber (Q_cloud, qnwfa)

      implicit none
      real, intent(in):: q_cloud, qnwfa
      real, parameter:: am_r = pi*1000./6.
      real, dimension(15), parameter:: g_ratio = (/24,60,120,210,336,   &
         504,720,990,1320,1716,2184,2730,3360,4080,4896/)
      double precision:: lambda, qnc
      real:: q_nwfa, x1, xDc
      integer:: nu_c

      if (Q_cloud == 0) then
         make_DropletNumber = 0
         return
      end if

      !+---+

      q_nwfa = MAX(99.e6, MIN(qnwfa,5.e10))
      nu_c = MAX(2, MIN(NINT(2.5e10/q_nwfa), 15))

      x1 = MAX(1., MIN(q_nwfa*1.e-9, 10.)) - 1.
      xDc = (30. - x1*20./9.) * 1.e-6

      lambda = (4.0D0 + nu_c) / xDc
      qnc = Q_cloud / g_ratio(nu_c) * lambda*lambda*lambda / am_r
      make_DropletNumber = SNGL(qnc)

      return
   end function make_DropletNumber

   !+---+-----------------------------------------------------------------+
   !+---+-----------------------------------------------------------------+

   elemental real function make_RainNumber (Q_rain, temp)

      implicit none

      real, intent(in):: q_rain, temp
      double precision:: lambda, n0, qnr
      real, parameter:: am_r = pi*1000./6.

      if (Q_rain == 0) then
         make_RainNumber = 0
         return
      end if

      !+---+-----------------------------------------------------------------+
      !.. Not thrilled with it, but set Y-intercept parameter to Marshal-Palmer value
      !.. that basically assumes melting snow becomes typical rain. However, for
      !.. -2C < T < 0C, make linear increase in exponent to attempt to keep
      !.. supercooled collision-coalescence (warm-rain) similar to drizzle rather
      !.. than bigger rain drops.  While this could also exist at T>0C, it is
      !.. more difficult to assume it directly from having mass and not number.
      !+---+-----------------------------------------------------------------+

      N0 = 8.e6

      if (temp .le. 271.15) then
         N0 = 8.e8
      elseif (temp .gt. 271.15 .and. temp.lt.273.15) then
         N0 = 8. * 10**(279.15-temp)
      endif

      lambda = SQRT(SQRT(N0*am_r*6.0/Q_rain))
      qnr = Q_rain / 6.0 * lambda*lambda*lambda / am_r
      make_RainNumber = SNGL(qnr)

      return
   end function make_RainNumber
   !+---+-----------------------------------------------------------------+
   !DSM {
   pure function intfuncgamma(x, y) result(z)
      real :: z
      real, intent(in) :: x, y

      z = x**(y-1.0) * exp(-x)
   end function intfuncgamma

   function gammaBrams(a) result(g)
      real :: g
      real, intent(in) :: a

      real, parameter :: small = 1.0e-4
      integer, parameter :: points = 100000

      real :: infty, dx, p, sp(2, points), x
      integer :: i
      logical :: correction

      x = a

      correction = .false.
      ! value with x<1 gives \infty, so we use
      ! \Gamma(x+1) = x\Gamma(x)
      ! to avoid the problem
      if ( x < 1.0 ) then
         correction = .true.
         x = x + 1
      end if

      ! find a "reasonable" infinity...
      ! we compute this integral indeed
      ! \int_0^M dt t^{x-1} e^{-t}
      ! where M is such that M^{x-1} e^{-M} ≤ \epsilon
      infty = 1.0e4
      do while ( intfuncgamma(infty, x) > small )
         infty = infty * 10.0
      end do

      ! using simpson
      dx = infty/real(points)
      sp = 0.0
      forall(i=1:points/2-1) sp(1, 2*i) = intfuncgamma(2.0*(i)*dx, x)
      forall(i=1:points/2) sp(2, 2*i - 1) = intfuncgamma((2.0*(i)-1.0)*dx, x)
      g = (intfuncgamma(0.0, x) + 2.0*sum(sp(1,:)) + 4.0*sum(sp(2,:)) + &
         intfuncgamma(infty, x))*dx/3.0

      if ( correction ) g = g/a

   end function gammaBrams
   !DSM}
   !----------------------------------------------------------------------
   subroutine gen_random(its,ite,use_random_num,random)
      implicit none
      integer, intent(in)  :: its,ite
      real,    intent(in)  :: use_random_num
      real,    intent(out) :: random(its:ite)

      !-local vars
      integer   :: i
      integer(8) :: iran, ranseed = 0

      call system_clock(ranseed)
      ranseed=mod(ranseed,2147483646)+1 !seed between 1 and 2^31-2
      iran = -ranseed

      !-- ran1 produces numbers between [ 0,1]
      !-- random        will be between [-1,1]
      !-- with use_random_num the interval will be [-use_random_num,+use_random_num]
      do i=its,ite
         random(i) = use_random_num * 2.0*(0.5-real(RAN1(IRAN),4))
        !print*,"ran=",i,random(i)
      enddo

      if(maxval(abs(random)) > use_random_num) stop "random > use_random_num"

   end subroutine gen_random
   !----------------------------------------------------------------------

   real(8) function ran1(idum)

      ! This is contributed code standardized by Yong Wang
      ! Random number generator taken from Press et al.
      !
      ! Returns numbers in the range 0-->1
      !
      ! Their description...
      ! "Minimal" random number generator of Park and Miller with Bays-Durham
      ! shuffle and added safeguards. Returns a uniform deviate between 0.0 and 1.0
      ! (exclusive of the endpoint values). Call with idum a negative integer to
      ! initialize; thereafter, do not alter idum between successive calls in a
      ! sequence. RNMX should approximate the largest floating value that is less
      ! than 1.

      !use shr_kind_mod,    only: r8 => shr_kind_r8, i8 => shr_kind_i8
      implicit none
      integer(8), parameter:: ntab = 32,iq = 127773,ia = 16807,ir = 2836, &
         im = 2147483647,ndiv = 1+(im-1)/ntab

      real(8), parameter:: am = 1.0/im,eps = 1.2e-7,rnmx = 1.0-eps

      integer(8), intent(inout):: idum

      integer(8):: iy
      integer(8), dimension(ntab):: iv
      !save iv,iy
      data iv /ntab*0/, iy /0/
      integer(8):: j,k

      !
      if (idum.le.0.or.iy.eq.0) then
         ! initialize
         idum = max(-idum,1)
         do j = ntab+8,1,-1
            k = idum/iq
            idum = ia*(idum-k*iq)-ir*k
            if (idum.lt.0) idum = idum+im
            if (j.le.ntab) iv(j) = idum
         end do
         iy = iv(1)
      end if
      !
      k = idum/iq
      ! compute idum = mod(ia*idum,im) without overflows by schrage's method
      idum = ia*(idum-k*iq)-ir*k
      if (idum.lt.0) idum = idum+im
      ! j will be in the range 1-->ntab
      j = 1+iy/ndiv
      ! output previously stored value and refill the shuffle table
      iy = iv(j)
      iv(j) = idum
      ran1 = min(am*iy,rnmx)

   end function ran1
   !----------------------------------------------------------------------

   subroutine get_delmix(cumulus,kts,kte,ktf,xland,subcl_level,po,ain,aout)
      implicit none
      character *(*)   ,intent (in) :: cumulus
      integer,intent(in)            :: kts,kte,ktf,subcl_level
      real   ,intent(in)            :: ain(kts:kte),po(kts:kte),xland
      real   ,intent(inout)         :: aout(kts:kte)

      !-- local var
      real :: x1,x2,dp,del,qc
      integer :: k

      !-
      qc = aout(kts)

      x2=0.
      x1=0.
      do k = kts,subcl_level
         dp = po(k+1)-po(k)
         x2 = x2 + dp
         x1 = x1 + dp*ain(k)
      enddo
      del = abs(qc-x1/(x2+1.e-12))
      aout(kts:subcl_level) =  ain(kts:subcl_level) + del

   end subroutine get_delmix
   !----------------------------------------------------------------------

   subroutine get_Qadv(cumulus,itf,ktf,its,ite,kts,kte,ierr,dt,q,qo,qo_adv,po,po_cup &
      ,qeso, Q_adv,col_sat_adv,alpha_adv,tau_bl,zo_cup,kbcon,ktop)
      implicit none
      character *(*), intent (in)                      :: cumulus
      integer ,intent (in)                             :: itf,ktf, its,ite, kts,kte
      real    ,intent (in)                             :: dt
      integer ,intent (in) ,dimension(its:ite)         :: ierr,kbcon,ktop

      real ,intent (in)    ,dimension(its:ite,kts:kte) :: q,qo,qo_adv,po_cup,qeso,po,zo_cup
      real ,intent (in)    ,dimension(its:ite)         :: tau_bl

      real ,intent (inout) ,dimension(its:ite)         :: Q_adv,col_sat_adv,alpha_adv

      !--locals
      integer :: i,k
      real :: dp, layer,H_cloud,dz
      real ,parameter :: ptop = 60.
      !-- get the advective moisture tendency scaled with the relative humidity
      !--  Q_adv = integral( q/q*  DQv/Dt_adv dp), see Eq 1 Becker et al(2021 QJRMS)
      !-- units here are "J m^-3" _or_  "J kg^-1"

      do i=its,itf
         col_sat_adv(i) = 0.   !check if it needs be inout, perhavps only local var

         if(ierr(i) /= 0) cycle

         alpha_adv  (i) = alpha_adv_tuning
         layer = 0.

         loopN: do k=kts,ktf
            if(po(i,k) < ptop) exit loopN

            !dp=100.*(po_cup(i,k+1)-po_cup(i,k)) ! dp < 0.
            dz=      zo_cup(i,k+1)-zo_cup(i,k)  ! dz > 0

            !-- integral over dp
            !Q_adv(i) = Q_adv(i) + dp*(qo(i,k)/qeso(i,k))*(qo_adv(i,k)-q(i,k))/dt

            !-- integral over dz
            Q_adv(i) = Q_adv(i) + dz*(qo(i,k)/qeso(i,k))*(qo_adv(i,k)-q(i,k))/dt

            col_sat_adv(i) = col_sat_adv(i) + dz* qo(i,k)/qeso(i,k)

            layer = layer + dz

         enddo loopN
         !-- get the column-average saturation fraction
         col_sat_adv(i) = col_sat_adv(i)/(1.e-8+layer)

         !--check if the col-ave saturation fraction is over the threshold
         if(col_sat_adv(i) > col_sat_adv_threshold) then

            alpha_adv(i) = 0.0
            cycle

         endif

         !-- check if cloud top _OR_cloud layer   !<<<< check this
         H_cloud =  zo_cup(i,ktop(i))- zo_cup(i,kbcon(i))

         !-- convert Q_adv to units as in Eq (1) => J m^-3
         !Q_adv(i) = - Q_adv(i) * tau_bl(i) * xlv / (g * H_cloud)

         !-- convert Q_adv to units as in cloud work function => J kg-1
         Q_adv(i) =  Q_adv(i) * tau_bl(i) * xlv / (H_cloud)

          !if(abs(q_adv(i))>1.) print*,"Qadv=",i,q_adv(i),Q_adv_dz(i)call flush(6)
      enddo

   end subroutine get_Qadv
!----------------------------------------------------------------------
 subroutine rh_controls(itf,ktf,its,ite,kts,kte,ierr,t,po,qo,qeso,po_cup &
                       ,cumulus,rh_entr_factor, rh_dicycle_fct           &
                       ,entr_rate_input, entr_rate,xlons,dt)
           
      implicit none
      character *(*), intent (in)                         :: cumulus
      integer  ,intent (in )                              :: itf,ktf, its,ite, kts,kte
      integer  ,intent (in )  ,dimension(its:ite)         :: ierr
      real     ,intent (in )                              :: entr_rate_input,dt
      real     ,intent (in )  ,dimension(its:ite)         :: xlons
      real     ,intent (in )  ,dimension(its:ite,kts:kte) :: t,po,qo,po_cup,qeso
      real     ,intent (inout),dimension(its:ite)         :: entr_rate 
      real     ,intent (inout),dimension(its:ite)         :: rh_entr_factor,rh_dicycle_fct

      !--locals
      integer :: i,k  
      real*8  :: y,x
      real    :: dpg, trash, dayhr, p_start = 1000.
      real    ,dimension(its:ite) :: frh,dayhrr
      real    ,parameter :: ref_local_time = 8., ftun3=0.25
      logical ,parameter :: free_troposphere = .true. 
      
      if(moist_trigger /= 2 .and. rh_dicycle == 0)return  
      !
      !
      !-- ave rh from 1000 -> 450 hPa, following Tian et al 2022 GRL.
      ! OR
      !-- ave rh from 800 -> 450 hPa accounts only for the ´free troposphere'
      if(free_troposphere) p_start = 800. 

      do i=its,itf        
         frh(i) = 0.
         trash  = 0.
  
         loopN:    do k=kts,ktf
            if( po(i,k) .gt. p_start .and. po(i,k) .lt. 450.) cycle loopN
            dpg=100.*(po_cup(i,k)-po_cup(i,k+1))/g
            trash=trash+dpg
            frh(i)= frh(i) + (qo(i,k)/qeso(i,k))*dpg 
         enddo loopN
         
         !--average relative humidity
         frh(i) =   100.*frh(i)/(1.e-8+trash) ! no unit
         frh(i) = max(1., min(100., frh(i)))
         !
         !--- this is for the moist_trigger = 2
          x = dble(frh(i))          
          y = 9.192833D0 - 0.2529055D0*x + 0.002344832d0*x**2 &
            - 0.000007230408D0*x**3
          rh_entr_factor(i) = real(y,4)
      
          
          !--- local time
          dayhr  = (time_in / 3600. + float(itime1_in/100) &
                 + float(mod(itime1_in,100)) / 60.) 
          dayhrr(i) = mod(dayhr+xlons(i)/15.+24., 24.)
          
           !print*,"FRH=",i,frh(i),rh_dicycle_fct(i),dayhrr,time_in/3600.,xlons(i)
           !print*,"LONS=",i,dayhrr,time_in/3600.,xlons(i)
           !print*,"=================================================="
      enddo
      if(moist_trigger == 2) then 
            entr_rate(:) = entr_rate_input *rh_entr_factor(:)
            !print*,"rh-entr-fac",minval(rh_entr_factor),maxval(rh_entr_factor)
      endif

      if(rh_dicycle == 1) then 
       do i=its,itf
         if(abs ( dayhrr(i) - ref_local_time) < 1. .or. time_in < dt+1. ) &  
           !--- ftun3 controls the domain of the dicycle closure 
           !    ftun3 => 1 the closure is insensitive to the mean tropospheric RH
           !    ftun3 => 0 it depends on the RH, with mininum = ftun3
           !rh_dicycle_fct(i) = ftun3 +(1. - ftun3)*&
           !                 (1.-(atan((frh(i)-60.)/10.)+atan(50.))/3.1016)/0.9523154
           rh_dicycle_fct(i) = ftun3 +(1. - ftun3)*&
                            (1.-(atan((frh(i)-55.)/10.)+atan(55.))/3.1016)
           !print*,"fac=",xlons(i),frh(i), rh_dicycle_fct(i);call flush(6)
       enddo
      endif
      !-- code to test the atan function
      !  do i = 1 ,100 !relative humidity
      !      y = 0.25 +0.75*(1.-(atan((float(i)-60.)/10.)+atan(50.))/3.1016)/0.9523154
      !      print*,i,y
      !  enddo
      
      !print*,"FRH",maxval(frh),minval(frh),maxval(rh_dicycle_fct),minval(rh_dicycle_fct)
      !call flush(6)
   end subroutine rh_controls
   !------------------------------------------------------------------------------------
   real function coldpool_start(CNV_TR)
      implicit none
      real,intent(in)  :: CNV_TR
      real             :: f1
      real,parameter   :: width = 100. !orig 100.

      f1= min (mx_buoy2,CNV_TR)
      !--- f1 > mx_buoy1 => coldpool_start ---> 1 
      coldpool_start =  (1.35+atan( (f1-mx_buoy1)/width))/2.8
      coldpool_start =  max(0.00,min(coldpool_start,1.00))
     !coldpool_start =  max(0.05,min(coldpool_start,0.95))
   end function
   !------------------------------------------------------------------------------------

end  module ConvPar_GF_GEOS5

