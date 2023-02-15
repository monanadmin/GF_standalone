module modConvParGF
   !! ## GF Convective parametrization
   !!
   !! ![](https://i.ibb.co/LNqGy3S/logo-Monan-Color-75x75.png)
   !! ## MONAN
   !!
   !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
   !!
   !! E-mail: <mailto:saulo.freitas@inpe.br>
   !!
   !! Date: 2014
   !!
   !! #####Version: 0.1.0
   !!
   !! ---
   !! **Full description**:
   !!
   !! This convective parameterization is build to attempt                      
   !! a smooth transition to cloud resolving scales as proposed                 
   !! by Arakawa et al (2011, ACP). The scheme is  described                    
   !! in the paper Grell and Freitas (ACP, 2014).                               
   !!
   !! Implemented in GEOS5 GCM by Saulo Freitas (July 2016)                     
   !! Use the following references for this implementation:                     
   !! Freitas et al (2018, JAMES/AGU, https://doi.org/10.1029/2017MS001251)     
   !! Freitas et al (2021, GMD/EGU,   https://doi.org/10.5194/gmd-14-5393-2021) 
   !! Please, contact Saulo Freitas (saulo.r.de.freitas@gmail.com) for comments 
   !! questions, bugs, etc.                                                     
   !!
   !! ** History**:
   !!
   !! -Adapted for BRAMS 6.0 by Saulo Freitas (November 2021)                    
   !! -Refactoring by Luiz Flavio Rodrigues at 20 December 2021 (Monday)         
   !! Keywords using ; are separeted, some loops receives exit instead goto,    
   !! The identation was fixed and all keywords are lowercase                   
   !! -Refactoring by GCC (INPE) at 20 January 2023 using fprettify and manual
   !! changes according MONAN rules code patterns DTN 01
   !!
   !! --- 
   !! ** Licence **:
   !!
   !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
   !!
   !! This program is free software: you can redistribute it and/or modify
   !! it under the terms of the GNU General Public License as published by
   !! the  Free  Software  Foundation, either version 3 of the License, or
   !! (at your option) any later version.
   !!
   !! This program is distributed in the hope that it  will be useful, but
   !! ** WITHOUT  ANY  WARRANTY **;  without  even  the   implied   warranty  of
   !! **MERCHANTABILITY** or **FITNESS FOR A  PARTICULAR PURPOSE**.  See  the, GNU
   !! GNU General Public License for more details.
   !!
   !! You should have received a copy  of the GNU General  Public  License
   !! along with this program.  If not, see [GNU Public License](https://www.gnu.org/licenses/gpl-3.0.html).
   !!
   use modGate, only: cupout, rundata, p_nvar_grads, jl, p_use_gate, ppres, ptemp, pq, pu &
                  ,   runlabel, runname, pv, pvervel, pgeo, zqr, zadvq, zadvt
   use modHenrysLawCts, only: getHenrysLawCts
   use modConstants, only: c_rgas, c_cp, c_rm, c_p00, c_tcrit, c_grav, c_cpor, c_alvl, c_pi &    
                        ,  c_akmin, c_tkmin_ms, c_ccnclean, c_t00, c_t_ice, c_xlf, c_max_qsat &
                        ,  p_xmbmaxshal, p_mintracer, c_smaller_qv, c_t01, c_t100, c_temp0i &
                        ,  c_rgas_atm, c_hplus, c_r2es, c_r3les, c_r4ies, c_r4les, c_retv &
                        ,  c_rticecu, c_rtwat_rticecu_r, c_r3ies, c_r5alscp, c_r5alvcp, c_ralsdcp &
                        ,  c_ralvdcp, c_rtice, c_rtwat_rtice_r, kind_ib, kind_rb

   implicit none

   private
   public p_maxiens, icumulus_gf, closure_choice, p_deep, p_shal, p_mid &
      , use_scale_dep, dicycle, tau_deep, tau_mid, hcts &
      , use_tracer_transp, use_tracer_scaven, use_memory, convection_tracer &
      , use_flux_form, use_tracer_evap, downdraft, use_fct &
      , use_rebcb, vert_discr, satur_calc, clev_grid, apply_sub_mp, alp1 &
      , sgs_w_timescale, lightning_diag, tau_ocea_cp, tau_land_cp &
      , autoconv, bc_meth, overshoot, use_wetbulb &
      , c1, c0_deep, qrc_crit, lambau_deep, lambau_shdn, c0_mid &
      , cum_max_edt_land, cum_max_edt_ocean, cum_hei_down_land &
      , cum_hei_down_ocean, cum_hei_updf_land, cum_hei_updf_ocean &
      , use_momentum_transp, cum_entr_rate &
      , p_nmp, p_lsmp, p_cnmp, moist_trigger, frac_modis, max_tq_tend &
      , cum_fadj_massflx, cum_use_excess, cum_ave_layer, adv_trigger &
      , use_smooth_prof, output_sound, use_cloud_dissipation &
      , use_smooth_tend, GFConvparInit, beta_sh, c0_shal &
      , use_linear_subcl_mf, cap_maxs, liq_ice_number_conc, alpha_adv_tuning &
      , sig_factor, lcl_trigger, rh_dicycle, add_coldpool_prop &
      , add_coldpool_clos, mx_buoy1, mx_buoy2, cum_t_star, cum_zuform &
      , add_coldpool_diff, initModConvParGF, modConvParGF_initialized

   public modConvParGFDriver, MakeDropletNumber, MakeIceNumber, fractLiqF &
      , use_gustiness, use_random_num, dcape_threshold, ColdPoolStart
   !
   !=================================================
   ! module parameters
   !=================================================
   character(len=*), parameter :: p_source_name = "modCobvParGF.F90"
   !! Name of source file
   integer, parameter :: p_maxiens = 3
   !! plume spectral size
   integer, parameter :: p_deep = 1
   !! plume spectral size
   integer, parameter :: p_shal = 2
   !! plume spectral size
   integer, parameter :: p_mid = 3
   !! plume spectral size
   character(len=10), parameter, dimension(p_maxiens)  :: p_cumulus_type = (/ &
                                                          'deep      ' &
                                                        , 'shallow   ' &
                                                        , 'mid       ' &
                                                        /)
   !! Cumulus type
   integer, parameter  :: p_nmp = 2
   !! number of microphysics schemes in the host model
   integer, parameter  :: p_lsmp = 1
   !! number of microphysics schemes in the host model
   integer, parameter  :: p_cnmp = 2
   !! number of microphysics schemes in the host model
   
   !-- General internal controls for the diverse options in GF
   logical, parameter :: p_entr_new = .true.  
   !! new entr formulation
   logical, parameter :: p_coupl_mphysics = .true.  
   !! coupling with cloud microphysics (do not change to false)
   logical, parameter :: p_melt_glac = .true.  
   !! turn ON/OFF ice phase/melting
   logical, parameter :: p_feed_3D_model = .true.  
   !! set "false" to not feedback the AGCM with the
   !! heating/drying/transport conv tendencies
   integer, parameter :: p_aeroevap = 1             
   !! rainfall evaporation (1) orig  - (2) mix orig+new - (3) new
   integer, parameter :: p_maxens = 1
   ! ensemble one on cap_max
   integer, parameter :: p_maxens2 = 1
   ! ensemble two on precip efficiency
   integer, parameter :: p_maxens3 = 16
   !ensemble three done in cup_forcing_ens16 for G3d
   integer, parameter :: p_ensdim = p_maxens*p_maxens2*p_maxens3
   !!
   integer, parameter :: p_ens4 = 1
   !!
   real, parameter :: p_pgcon = 0.0
   !! proportionality constant to estimate pressure
   !! gradient of updraft (Zhang and Wu, 2003, JAS) => REAL, PARAMETER ::    pgcon=-0.55

   integer, parameter :: p_max_n_spec = 200
   !!
   integer, parameter :: p_shall_closures = 12
   !!
   integer, parameter :: p_on = 1
   !! ON integer paremeters
   integer, parameter :: p_off = 0 
   !! OFF integer paremeters
   !=================================================
   ! End of module parameters
   !=================================================

   !=================================================
   ! namelist variables
   !=================================================
   real :: CUM_ENTR_RATE(p_maxiens)
   !! -- gross entraiment rate: deep, shallow, congestus
   real :: ALP1
   !! 0/0.5/1: apply subsidence transport of LS/anvil cloud fraction using
   !!         time implicit discretization
   real ::  TAU_DEEP
   !! deep      convective timescale
   real ::  TAU_MID
   !! congestus convective timescale
   real :: MAX_TQ_TEND
   !! max T,Q tendency allowed (100 K/day)
   real :: OVERSHOOT
   !!  values = 0 or 1 simulates an additional capping layer above the cloud top
   real ::  C0_DEEP
   !! default= 3.e-3   conversion rate (cloud to rain, m-1) - for deep      plume
   real ::  C0_MID
   !! default= 2.e-3   conversion rate (cloud to rain, m-1) - for congestus plume
   real ::  C0_SHAL
   !! default= 0.e-3   conversion rate (cloud to rain, m-1) - for shallow   plume
   real ::  QRC_CRIT
   !! default= 2.e-4   kg/kg
   real ::  C1 
   !! default= 1.e-3   conversion rate (cloud to rain, m-1) - for the 'C1d' detrainment approach
   real ::  LAMBAU_DEEP
   !! default= 2.0 lambda parameter for deep/congestus convection momentum transp
   real ::  LAMBAU_SHDN
   !! default= 2.0 lambda parameter for shallow/downdraft convection momentum transp
   real :: CUM_HEI_DOWN_LAND(p_maxiens)
   !! [0.2,0.8] height of the max Z Downdraft , default = 0.50
   real :: CUM_HEI_DOWN_OCEAN(p_maxiens)
   !! [0.2,0.8] height of the max Z Downdraft , default = 0.50
   real :: CUM_HEI_UPDF_LAND(p_maxiens)
   !! [0.2,0.8] height of the max Z Updraft   , default = 0.35
   real :: CUM_HEI_UPDF_OCEAN(p_maxiens)
   !! [0.2,0.8] height of the max Z Updraft   , default = 0.35
   real :: CUM_MAX_EDT_LAND(p_maxiens)
   !! maximum evap fraction allowed over the land  ,default= 0.9
   real :: CUM_MAX_EDT_OCEAN(p_maxiens)
   !! maximum evap fraction allowed over the ocean ,default= 0.9
   real :: CUM_FADJ_MASSFLX(p_maxiens)
   !! multiplicative factor for tunning the mass flux at cloud base
   != default = 1.0
   real :: CUM_T_STAR(p_maxiens)
   !! scale temperature for the diurnal cycle closure
   real :: CUM_AVE_LAYER(p_maxiens)
   !! layer depth for average the properties
   !! of source air parcels (mbar) = (/15., -99., -99./)
   !! scale temperature for the diurnal cycle closure
   real :: DCAPE_THRESHOLD
   !! CAPE time rate threshold for ADV_TRIGGER = 1 (J kg^-1 hr^-1)
   !! typical range is [-200,200] J/kg/hr, Wu et all (2007) recomends ~ 70 J/kg/hr
   !! 55 J/kg/hr is indicated for the Amazon basin (Song&Zhang 2017)
   real :: TAU_OCEA_CP
   !! default 7200 s:  cold pool timescale over the oceans
   real :: TAU_LAND_CP
   !! default 7200 s:  cold pool timescale over land 
   real :: MX_BUOY1
   !! 250.5 J/kg: Lower range of buoyancy-excess 
   real :: MX_BUOY2
   != 20004.0 J/kg: Upper range of buoyancy-excess (temp exc=10 K, q deficit=4 g/kg => mx_buoy ~ 20 kJ/kg)
   real :: USE_CLOUD_DISSIPATION
   !! values 0. or 1. to acccount for the cloud dissipation at the decayment phase (implies on extra moistening and cooling)
   real :: USE_RANDOM_NUM
   !! stochastic pertubation for the height of maximum Zu
   real :: BETA_SH
   !! default 2.2: Beta factor of the Beta PDF (only for shallow plume)
   real :: ALPHA_ADV_TUNING
   !! default = 0.8 : tuning parameter for the Becker et al (2021) closure
   real :: SIG_FACTOR
   !! default 2.2 : exponential factor for the sigma determination
   real :: CAP_MAXS
   !! defalut= 50 hPa : max distance (hPa) the air parcel is allowed to go up looking for the LFC
   integer :: ICUMULUS_GF(p_maxiens)
   !! define each plume will be activated (1 true, 0 false): deep, shallow, congestus
   integer :: CLOSURE_CHOICE(p_maxiens)
   !! -- choice for the closures:
   !! --  deep   : 0 ensemble (all)          , 1 GR, 4 ll omega, 7 moist conv, 10 PB
   !! --  shallow: 0 ensemble (all)          , 1 Wstar, 4 heat-engine, 7 BLQE, 10 TKE-based
   !! --  mid    : 0 ensemble (Wstar/BLQE/PB), 1 Wstar, 2 BLQE, 3 PB, 4 PB_BL
   integer :: CUM_ZUFORM(p_maxiens)
   !! -- zu updraft format : deep, shallow, congestus
   !! for deep: 10 or 20, congestus: 20
   integer :: USE_TRACER_TRANSP
   != 0/1     - default 1
   integer :: USE_TRACER_SCAVEN
   !! 0/1/2/3 - default 2
   integer :: USE_FLUX_FORM
   !! 1/2/3   - default 1
   integer :: USE_FCT 
   !! 0/1     - default 1 (only for USE_FLUX_FORM     = 2)
   integer :: USE_TRACER_EVAP 
   !! 0/1     - default 1 (only for USE_TRACER_SCAVEN > 0)
   integer :: CONVECTION_TRACER
   !! 0/1:  turn ON/OFF the "convection" tracer
   integer :: USE_MEMORY
   !! -1/0/1/2 .../10    !-
   integer :: ADD_COLDPOOL_PROP
   !! -1,0,1,2,3 add coldpool propagation
   integer :: ADD_COLDPOOL_CLOS
   !! add the kinetic energy at leading of the gust front
   integer :: ADD_COLDPOOL_DIFF
   !! add vert/horizontal diffusion to the cold pool propaga
   integer :: USE_SCALE_DEP
   !! 0/1:  scale dependence flag, default = 1
   integer :: DICYCLE
   !! 0/1/2:  diurnal cycle closure, default = 1
   !! 2 uses Qadv closure (Becker et al 2021)
   integer :: RH_DICYCLE
   !! controls of RH on the diurnal cycle (see Tian et al 2022 GRL)
   integer :: CLEV_GRID
   !! 0/1/2: interpolation method to define environ state at the
   !! cloud levels (at face layer), default = 0
   !! CLEV_GRID = 0 default method
   !! CLEV_GRID = 1 interpolation method based on Tiedtke (1989)
   !! CLEV_GRID = 2 for GATE soundings only
   integer :: USE_REBCB
   !! 0/1: turn ON/OFF rainfall evap below cloud base, default = 0
   integer :: VERT_DISCR
   !! 0/1: 1=new vert discretization, default = 0
   integer :: SATUR_CALC
   !! 0/1: 1=new saturation specific humidity calculation, default = 0
   integer :: SGS_W_TIMESCALE
   !! 0/1: vertical velocity for tau_ecmwf, default = 0
   integer :: LIGHTNING_DIAG
   !! 0/1: do LIGHTNING_DIAGgnostics based on Lopez (2016, MWR)
   integer :: APPLY_SUB_MP
   !! 0/1: subsidence transport applied the to grid-scale/anvil ice/liq mix
   !!     ratio and cloud fraction
   integer :: USE_WETBULB
   !! 0/1
   integer :: BC_METH
   !! boundary condition determination for the plumes
   !! 0: simple arithmetic mean around the source level
   !! 1: mass weighted mean around the source level
   integer :: AUTOCONV
   !! 1, 3 or 4 autoconversion formulation: (1) Kessler,
   !! (3) Kessler with temp dependence, (4) Sundvisqt
   integer :: USE_MOMENTUM_TRANSP
   !! 0/1:  turn ON/OFF conv transp of momentum
   integer :: DOWNDRAFT
   !! 0/1:  turn ON/OFF downdrafts, default = 1
   integer :: USE_SMOOTH_PROF
   !! 1 makes the normalized mass flux, entr and detraiment profiles smoother
   integer :: USE_SMOOTH_TEND
   !! 0 => OFF, > 0 produces smoother tendencies (e.g.: for 1=> makes average between k-1,k,k+1)
   !! deep, shallow, congestus
   integer :: CUM_USE_EXCESS(p_maxiens)
   !! use T,Q excess sub-grid scale variability
   integer :: MOIST_TRIGGER 
   !! relative humidity effects on the cap_max trigger function
   integer :: FRAC_MODIS
   !! use fraction liq/ice content derived from MODIS/CALIPO sensors
   integer :: ADV_TRIGGER
   !! dcape trigger
   integer :: LCL_TRIGGER
   !! greater than zero, activates the LCL trigger which requires the lcl height
   !! be lower than the pbl height, only for shallow convection
   integer :: USE_GUSTINESS
   !! not in use
   integer :: USE_LINEAR_SUBCL_MF
   !! only for shallow plume
   integer :: LIQ_ICE_NUMBER_CONC
   !! include drop/ice number mixing ratio convective tendencies
   !=================================================
   ! End of namelist variables
   !=================================================

!=================================================
   ! module internal variables  -
   !=================================================
   real :: hei_down_land     
   !! [0.2,0.8] height of the max Z Downdraft , default = 0.50
   real :: hei_down_ocean    
   !! [0.2,0.8] height of the max Z Downdraft , default = 0.50
   real :: hei_updf_land     
   !! [0.2,0.8] height of the max Z Updraft   , default = 0.35
   real :: hei_updf_ocean    
   !! [0.2,0.8] height of the max Z Updraft   , default = 0.35
   real :: max_edt_land      
   !! default= 0.9 - maximum evap fraction allowed over the land
   real :: max_edt_ocean     
   !! default= 0.9 - maximum evap fraction allowed over the ocean
   real :: fadj_massflx     
   !! default= 1.0 - multiplicative factor for the mass flux at cloud base
   real :: t_star         
   !! Scale Temperature for the DC closure
   real :: ave_layer         
   !! layer depth for average the properties of source air parcels (mbar)
   real :: c0                
   !! autoconversion constant
   real :: col_sat_adv_threshold
   !! suppress Qadv closure for col_sat_adv > col_sat_adv_threshold
   real ::  chem_adj_autoc(p_max_n_spec)
   !!
   real :: time_in
   !!
   real :: int_time
   !!
   integer :: ispc_co
   !!
   integer :: whoami_all
   !!
   integer :: jcol
   !!
   integer :: itime1_in
   !!
   integer :: nrec
   !!
   integer :: ntimes
   !!
   integer ::  use_excess        
   !! default= 1   - use T,Q excess sub-grid scale variability
   integer :: output_sound
   !! outputs a "GEOS" vertical profile for the GF stand alone model
   integer :: ind_chem(p_max_n_spec)
   !!
   integer :: chem_name_mask(p_max_n_spec)
   !!
   integer :: chem_name_mask_evap(p_max_n_spec)
   !!
   logical :: use_c1d
   !! turn ON/OFF the 'c1d' detrainment approach, don't change this.
   logical :: first_guess_w
   !! use it to calculate a 1st guess of the updraft vert velocity
   logical :: wrtgrads
   !!
   character(len=100)  ::  chem_name(p_max_n_spec)
   !!
   type t_hcts_vars
      real :: hstar
      !!
      real :: dhr
      !!
      real :: ak0
      !!
      real :: dak
      !!
   end type t_hcts_vars
   type(t_hcts_vars), allocatable :: hcts(:)
   !!
   logical :: modConvParGF_initialized
   !=================================================
   ! End of module internal variables  -
   !=================================================

contains

   !-----------------------------------------------------------------------
   subroutine modConvParGFDriver(mxp, myp, mzp, mtp, nmp, time, itime1 &
                           , ims, ime, jms, jme, kms, kme &
                           , its, ite, jts, jte, kts, kte &
                           , flip &
                           , fscav &
                           , mynum &
                           , dt &
                           , dx2d &
                           , stochastic_sig &
                           , zm &
                           , zt &
                           , dm &
                           , lons &
                           , lats &
                           , aot500 &
                           , temp2m &
                           , sflux_r &
                           , sflux_t &
                           , qexcp &
                           , hexcp &
                           , wlpool &
                           , topt &
                           , xland &
                           , sfc_press &
                           , kpbl &
                           , tke_pbl &
                           , col_sat &
                           , u_wind &
                           , v_wind &
                           , w_wind &
                           , temp &
                           , press &
                           , rvap &
                           , mp_ice &
                           , mp_liq &
                           , mp_cf &
                           , curr_rvap &
                           , tracer &!-note: uses GEOS-5 data structure
                           !---- forcings---
                           , buoy_exc &
                           , rthften &! gsf_t
                           , rqvften &! gsf_q
                           , rth_advten &!advf_t
                           , rthblten &!sgsf_t
                           , rqvblten &!sgsf_q
                           !---- output ----
                           , conprr &
                           , lightn_dens &
                           , rh_dicycle_fct &
                           , rthcuten &
                           , rqvcuten &
                           , rqccuten &
                           , rnlcuten &
                           , rnicuten &
                           , rucuten &
                           , rvcuten &
                           , sub_mpqi &
                           , sub_mpql &
                           , sub_mpcf &
                           , rbuoycuten &
                           , rchemcuten &
                           , revsu_gf &
                           , prfil_gf &
                           !
                           , do_this_column &
                           , ierr4d &
                           , jmin4d &
                           , klcl4d &
                           , k224d &
                           , kbcon4d &
                           , ktop4d &
                           , kstabi4d &
                           , kstabm4d &
                           , cprr4d &
                           , xmb4d &
                           , edt4d &
                           , pwav4d &
                           , sigma4d &
                           , pcup5d &
                           , up_massentr5d &
                           , up_massdetr5d &
                           , dd_massentr5d &
                           , dd_massdetr5d &
                           , zup5d &
                           , zdn5d &
                           , prup5d &
                           , prdn5d &
                           , clwup5d &
                           , tup5d &
                           , conv_cld_fr5d &
                           !-- for debug/diagnostic
                           , aa0, aa1, aa1_adv, a1_radpbl, aa1_bl, aa2, aa3, aa1_cin, tau_bl, tau_ec &
                           , var2d, var3d_agf, var3d_bgf, var3d_cgf, var3d_dgf)
      !! ## Driver for modConvParGF
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Driver for modConvParGF
      !!
      !! ** History**:
      !!
      !! - Itenizado_as_alterações_ao_longo_do_tempo (genérica)
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !!
      implicit none
      !Parameters:
      character(len=*), parameter :: procedureName = 'modConvParGFDriver' ! Subroutine name
   
      !Variables (input, output, inout)
            !------------------------------------------------------------------------
      integer, intent(in) :: ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, mzp, mxp, myp, mtp, nmp, mynum &
                           , itime1
      integer, intent(in) :: flip(:)
      integer, intent(in) :: kpbl(its:ite, jts:jte)

      real, intent(in) :: dt
      real, intent(in) :: time
      real, intent(in) :: zm(kts:kte, its:ite, jts:jte)
      real, intent(in) :: zt(kts:kte, its:ite, jts:jte)
      real, intent(in) :: u_wind(kts:kte, its:ite, jts:jte)
      real, intent(in) :: v_wind(kts:kte, its:ite, jts:jte)
      real, intent(in) :: w_wind(kts:kte, its:ite, jts:jte)
      real, intent(in) :: rvap(kts:kte, its:ite, jts:jte)
      real, intent(in) :: temp(kts:kte, its:ite, jts:jte)
      real, intent(in) :: press(kts:kte, its:ite, jts:jte)
      real, intent(in) :: dm(kts:kte, its:ite, jts:jte)
      real, intent(in) :: curr_rvap(kts:kte, its:ite, jts:jte)
      real, intent(in) :: buoy_exc(kts:kte, its:ite, jts:jte)
      real, intent(in) :: qexcp(kts:kte, its:ite, jts:jte)
      real, intent(in) :: hexcp(kts:kte, its:ite, jts:jte)
      real, intent(in) :: rthften(kts:kte, its:ite, jts:jte)
      real, intent(in) :: rqvften(kts:kte, its:ite, jts:jte)
      real, intent(in) :: rth_advten(kts:kte, its:ite, jts:jte)
      real, intent(in) :: rthblten(kts:kte, its:ite, jts:jte)
      real, intent(in) :: rqvblten(kts:kte, its:ite, jts:jte)
      real, intent(in) :: mp_ice(nmp, kts:kte, its:ite, jts:jte)
      real, intent(in) :: mp_liq(nmp, kts:kte, its:ite, jts:jte)
      real, intent(in) :: mp_cf(nmp, kts:kte, its:ite, jts:jte)
      real, intent(in) :: tracer(its:ite, jts:jte, kts:kte, mtp)
      !! tracer has different data structure   (i,j,k,ispc) *********
      real, intent(in) :: topt(its:ite, jts:jte)
      real, intent(in) :: aot500(its:ite, jts:jte)
      real, intent(in) :: temp2m(its:ite, jts:jte)
      real, intent(in) :: sfc_press(its:ite, jts:jte)
      real, intent(in) :: sflux_r(its:ite, jts:jte)
      real, intent(in) :: sflux_t(its:ite, jts:jte)
      real, intent(in) :: xland(its:ite, jts:jte)
      real, intent(in) :: lons(its:ite, jts:jte)
      real, intent(in) :: lats(its:ite, jts:jte)
      real, intent(in) :: stochastic_sig(its:ite, jts:jte)
      real, intent(in) :: tke_pbl(its:ite, jts:jte)
      real, intent(in) :: col_sat(its:ite, jts:jte)

      integer, intent(inout) :: do_this_column(its:ite, jts:jte)
      
      !- for convective transport and cloud/radiation (OUT)
      integer, intent(inout) :: ierr4d(mxp, myp, p_maxiens)
      integer, intent(inout) :: jmin4d(mxp, myp, p_maxiens)
      integer, intent(inout) :: klcl4d(mxp, myp, p_maxiens)
      integer, intent(inout) :: k224d(mxp, myp, p_maxiens)
      integer, intent(inout) :: kbcon4d(mxp, myp, p_maxiens)
      integer, intent(inout) :: ktop4d(mxp, myp, p_maxiens)
      integer, intent(inout) :: kstabi4d(mxp, myp, p_maxiens)
      integer, intent(inout) :: kstabm4d(mxp, myp, p_maxiens)

      real, intent(inout) :: fscav(:)
      real, intent(inout) :: rh_dicycle_fct(its:ite, jts:jte)

      !input but communicate to another subroutine
      real, intent(inout) :: dx2d(its:ite, jts:jte)
      real, intent(inout) :: wlpool(its:ite, jts:jte)
      real, intent(inout) :: pcup5d(mxp, myp, mzp, p_maxiens)
      real, intent(inout) :: up_massentr5d(mxp, myp, mzp, p_maxiens)
      real, intent(inout) :: up_massdetr5d(mxp, myp, mzp, p_maxiens)
      real, intent(inout) :: dd_massentr5d(mxp, myp, mzp, p_maxiens)
      real, intent(inout) :: dd_massdetr5d(mxp, myp, mzp, p_maxiens)
      real, intent(inout) :: zup5d(mxp, myp, mzp, p_maxiens)
      real, intent(inout) :: zdn5d(mxp, myp, mzp, p_maxiens)
      real, intent(inout) :: prup5d(mxp, myp, mzp, p_maxiens)
      real, intent(inout) :: prdn5d(mxp, myp, mzp, p_maxiens)
      real, intent(inout) :: clwup5d(mxp, myp, mzp, p_maxiens)
      real, intent(inout) :: tup5d(mxp, myp, mzp, p_maxiens)
      real, intent(inout) :: conv_cld_fr5d(mxp, myp, mzp, p_maxiens)
      real, intent(inout) :: cprr4d(mxp, myp, p_maxiens)
      real, intent(inout) :: xmb4d(mxp, myp, p_maxiens)
      real, intent(inout) :: edt4d(mxp, myp, p_maxiens)
      real, intent(inout) :: pwav4d(mxp, myp, p_maxiens)
      real, intent(inout) :: sigma4d(mxp, myp, p_maxiens)

      !--for debug
      real, intent(inout) :: aa0(mxp, myp)
      real, intent(inout) :: aa1(mxp, myp)
      real, intent(inout) :: aa1_adv(mxp, myp)
      real, intent(inout) :: a1_radpbl(mxp, myp)
      real, intent(inout) :: aa2(mxp, myp)
      real, intent(inout) :: aa3(mxp, myp)
      real, intent(inout) :: aa1_bl(mxp, myp)
      real, intent(inout) :: aa1_cin(mxp, myp)
      real, intent(inout) :: tau_bl(mxp, myp)
      real, intent(inout) :: tau_ec(mxp, myp)
      real, intent(inout) :: var2d(mxp, myp)

      real, intent(out) :: conprr(its:ite, jts:jte)
      real, intent(out) :: lightn_dens(its:ite, jts:jte)
      real, intent(out) :: rthcuten(kts:kte, its:ite, jts:jte)
      real, intent(out) :: rqvcuten(kts:kte, its:ite, jts:jte)
      real, intent(out) :: rqccuten(kts:kte, its:ite, jts:jte)
      real, intent(out) :: rnlcuten(kts:kte, its:ite, jts:jte)
      real, intent(out) :: rnicuten(kts:kte, its:ite, jts:jte)
      real, intent(out) :: rucuten(kts:kte, its:ite, jts:jte)
      real, intent(out) :: rvcuten(kts:kte, its:ite, jts:jte)
      real, intent(out) :: rbuoycuten(kts:kte, its:ite, jts:jte)
      real, intent(out) :: revsu_gf(kts:kte, its:ite, jts:jte)
      real, intent(out) :: prfil_gf(kts:kte, its:ite, jts:jte)
      real, intent(out) :: var3d_agf(kts:kte, its:ite, jts:jte)
      real, intent(out) :: var3d_bgf(kts:kte, its:ite, jts:jte)
      real, intent(out) :: var3d_cgf(kts:kte, its:ite, jts:jte)
      real, intent(out) :: var3d_dgf(kts:kte, its:ite, jts:jte)
      real, intent(out) :: sub_mpqi(nmp, kts:kte, its:ite, jts:jte)
      real, intent(out) :: sub_mpql(nmp, kts:kte, its:ite, jts:jte)
      real, intent(out) :: sub_mpcf(nmp, kts:kte, its:ite, jts:jte)
      real, intent(out) :: rchemcuten(mtp, kts:kte, its:ite, jts:jte)
      !! rchemcuten uses the GF data structure (ispc,k,i,j) *********


      !----------------------------------------------------------------------
      ! local variabels

      ! basic environmental input includes
      ! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
      ! convection for this call only and at that particular gridpoint
      !
      real, dimension(kts:kte, its:ite, jts:jte)  :: tpert_h, tpert_v

      real, dimension(its:ite, jts:jte) ::  rtgt

      real, dimension(its:ite, kts:kte) :: &
         zo, temp_old, qv_old, PO, US, VS, rhoi, phil &
         , temp_new_dp, qv_new_dp, temp_new_sh, qv_new_sh, z2d &
         , tkeg, rcpg, dhdt, temp_new_md, qv_new_md &
         , temp_new_bl, qv_new_bl, dm2d, temp_tendqv, qv_curr &
         , buoy_exc2d, revsu_gf_2d, prfil_gf_2d, var3d_Agf_2d &
         , var3d_Bgf_2d, temp_new, qv_new, Tpert_2d &
         , temp_new_adv, qv_new_adv

      real, dimension(its:ite, kts:kte, p_maxiens) :: &
         outt, outq, outqc, outu, outv, outbuoy, outnliq, outnice

      real, dimension(mtp, its:ite, kts:kte)         :: se_chem
      real, dimension(mtp, its:ite, kts:kte, p_maxiens) :: out_chem

      real, dimension(nmp, its:ite, kts:kte)         :: mpqi, mpql, mpcf
      real, dimension(nmp, its:ite, kts:kte, p_maxiens) :: outmpqi, outmpql, outmpcf

      real, dimension(its:ite)   :: ter11, xlandi, pbl, zws, ccn, psur &
                                    , ztexec, zqexec, h_sfc_flux, le_sfc_flux, tsur &
                                    , xlons, xlats, fixout_qv, cum_ztexec, cum_zqexec

      real, dimension(its:ite, kts:kte, 1:p_ens4)      ::  omeg

      real, dimension(kts:kte) :: min_tend, distance
      integer, dimension(its:ite) :: kpbli, last_ierr

      integer :: i, j, k, kr, n, itf, jtf, ktf, ispc, zmax, status

      real :: dp, dq, exner, dtdt, pten, pqen, paph, zrho, pahfs, pqhfl, zkhvfl, pgeoh
      real :: fixouts, dt_inv

      real, dimension(mxp, myp, -1:5) :: dummy_precip
      integer :: imemory, irun, jlx, kk, kss, plume, ii_plume

      !----------------------------------------------------------------------
      !-do not change this
      itf = ite
      ktf = kte - 1
      jtf = jte
      int_time = int_time + dt
      whoami_all = mynum
      time_in = time
      itime1_in = itime1
      !----------------------------------------------------------------------
      if (abs(C1) > 0.) use_c1d = .true.

      !-- big loop over j dimension
      do j = jts, jtf
         jcol = J

         ! Initialize the variables for each J
         call initDriverJ(its, itf, aot500(:,j), xland(:,j), sfc_press(:,j), temp2m(:,j), topt(:,j), kpbl(:,j) &
                       , lons(:,j), lats(:,j), rtgt(:,j), ztexec, zqexec, last_ierr, fixout_qv, conprr(:,j) &
                       , lightn_dens(:,j), var2d(:,j), revsu_gf_2d, prfil_gf_2d, var3d_agf_2d, var3d_bgf_2d &
                       , Tpert_2d, temp_tendqv, outt, outu, outv, outq, outqc, outnice, outnliq, outbuoy, outmpqi &
                       , outmpql, outmpcf, omeg, out_chem, ccn, xlandi, psur, tsur, ter11, kpbli, xlons, xlats)

         ! Initialize the variables for each J for all k
         call initDriverK(its, itf, kts, ktf, kte, mtp, dt, rtgt(:,j), topt(:,j), zt(:,:,j), press(:,:,j), temp(:,:,j) &
                        , rvap(:,:,j), curr_rvap(:,:,j), u_wind(:,:,j), v_wind(:,:,j), w_wind(:,:,j), dm(:,:,j) &
                        , buoy_exc(:,:,j) , rth_advten(:,:,j), rqvften(:,:,j), flip, mp_ice(:,:,:,j), mp_liq(:,:,:,j) &
                        , mp_cf(:,:,:,j), tracer(:,j,:,:), zo, po, temp_old, qv_old, qv_curr, rhoi, tkeg, rcpg, us, vs &
                        , dm2d, buoy_exc2d, temp_new_adv, qv_new_adv, mpqi, mpql, mpcf, se_chem, omeg)

         !- pbl  (i) = depth of pbl layer (m)
         !- kpbli(i) = index of zo(i,k)
         !call get_zi_gf(j,its,ite,kts,kte,its,itf,ktf,ierrs,kpbli,pbl,&
         !             tkeg,rcpg,zo,ter11,tkmin)
         do i = its, itf
            pbl(i) = zo(i, kpbli(i)) - topt(i, j)
            !print*,"PBL=",kpbli(i),zo(i,kpbli(i)),topt(i,j),pbl  (i)
         end do

         !- begin: for GATE soundings-------------------------------------------
         !- this section is intended for model developments only and must
         !- not be used for normal runs.
         if (p_use_gate) then
            if (CLEV_GRID == 0) stop "use_gate requires CLEV_GRID 1 or 2"
            if (USE_TRACER_TRANSP == 1) then
               ispc_co = 1
               if (.not. allocated(hcts)) allocate (hcts(mtp))
               chem_name_mask(:) = 1
               !--- dummy initization FSCAV
               do i = 1, mtp
                  !FSCAV(i) = 0.1  !km^-1

                  fscav(i) = 1.e-5  !km^-1
                  hcts(i)%hstar = 0.0 !8.300e+4! 2.4E+3 !59.
                  hcts(i)%dhr = 0.0 !7400.   !5000.  !4200.
                  hcts(i)%ak0 = 0.0
                  hcts(i)%dak = 0.0
                  ! H2O2      0.00000      8.300e+4    7400.00000       0.00000       0.00000
                  ! HNO3      0.00000      2.100e+5    8700.00000       0.00000       0.00000
                  ! NH3       0.00000      59.00000    4200.00000       0.00000       0.00000
                  ! SO2       0.00000      2.400e+3    5000.00000       0.00000       0.00000
               end do
               do i = its, itf
                  se_chem(1:mtp, i, kts:kpbli(i) - 1) = 1.+1.e-6
                  do k = kpbli(i), kte
                     se_chem(1:mtp, i, k) = 1.*exp(-(max(0., 0.9*float(k - kpbli(i)))/float(kpbli(i)))) + 1.e-6
                  end do
                  do k = kts + 1, kte - 1
                     se_chem(1:mtp, i, k) = 1./3.*(se_chem(1:mtp, i, k) + se_chem(1:mtp, i, k - 1) + se_chem(1:mtp, i, k + 1))
                  end do
               end do
            end if

            !--- only for GATE soundingg
            if (trim(rundata) == "GATE.dat") then
               jlx = jl
               !jlx= 1 ! to run with only one soundings
               !jlx= 42 ! to run with only one soundings

               do i = its, itf
                  do k = kts, kte
                     po(i, k) = 0.5*(ppres(jlx, k) + ppres(jlx, min(kte, k + 1)))
                     temp_old(i, k) = ptemp(jlx, k) + 273.15
                     qv_old(i, k) = pq(jlx, k)/1000.
                     us(i, k) = pu(jlx, k)
                     vs(i, k) = pv(jlx, k)
                     omeg(i, k, :) = pvervel(jlx, k)
                     phil(i, k) = pgeo(jlx, k)*c_grav   !geo
                     rhoi(i, k) = 1.e2*po(i, k)/(c_rgas*temp_old(i, k))
                  end do

                  do k = kts, kte
                     mpql(:, i, k) = 0.
                     mpql(:, i, k) = 0.
                     mpcf(:, i, k) = 0.
                     if (po(i, k) > 900. .or. po(i, k) < 300.) cycle
                     pqen = exp((-3.e-5*(po(i, k) - 550.)**2))
                     pten = min(1., (max(0., (temp_old(i, k) - c_t_ice))/(c_t00 - c_t_ice))**2)
                     mpql(:, i, k) = 3.*pqen*pten
                     mpqi(:, i, k) = 3.*pqen*(1.-pten)
                     mpcf(:, i, k) = (mpqi(:, i, k) + mpql(:, i, k))*100.
                  end do

                  do k = kts, kte
                     zo(i, k) = 0.5*(phil(i, k) + phil(i, min(kte, k + 1)))/c_grav    !meters
                  end do
                  ter11(i) = phil(i, 1)/c_grav  ! phil is given in g*h.
                  psur(i) = ppres(jlx, 1)
                  tsur(i) = temp2m(i, j) !temp_old(i,1)
                  kpbli(i) = 5
                  pbl(i) = zo(i, kpbli(i))
                  zws(i) = 1.0 ! wstar
                  do k = kts, ktf
                     temp_new(i, k) = temp_old(i, k) + dt*(zadvt(jlx, k) + zqr(jlx, k))/86400.
                     qv_new(i, k) = qv_old(i, k) + dt*zadvq(jlx, k)

                     temp_new_dp(i, k) = temp_old(i, k) + dt*(zadvt(jlx, k) + zqr(jlx, k))/86400.
                     qv_new_dp(i, k) = qv_old(i, k) + dt*zadvq(jlx, k)

                     temp_new_md(i, k) = temp_new_dp(i, k)
                     qv_new_md(i, k) = qv_new_dp(i, k)
                     temp_new_bl(i, k) = temp_new_dp(i, k)
                     qv_new_bl(i, k) = qv_new_dp(i, k)
                     temp_new_adv(i, k) = temp_old(i, k) + dt*zadvt(jlx, k)/86400.
                     qv_new_adv(i, k) = qv_old(i, k) + dt*zadvq(jlx, k)
                  end do
               end do
            end if
         end if !- end:   for GATE soundings-------------------------------------------
         !
         !- get execess T and Q for source air parcels
         do i = its, itf
            pten = temp_old(i, 1)
            pqen = qv_old(i, 1)
            paph = 100.*psur(i)
            zrho = paph/(287.04*(temp_old(i, 1)*(1.+0.608*qv_old(i, 1))))
            !- sensible and latent sfc fluxes for the heat-engine closure
            h_sfc_flux(i) = zrho*real(c_cp)*sflux_t(i, j)!W/m^2
            le_sfc_flux(i) = zrho*real(c_alvl)*sflux_r(i, j)!W/m^2
            !
            !- local le and h fluxes for W*
            pahfs = -sflux_t(i, j)*zrho*1004.64  !W/m^2
            pqhfl = -sflux_r(i, j)                !kg/m^2/s
            !- buoyancy flux (h+le)
            zkhvfl = (pahfs/1004.64 + 0.608*pten*pqhfl)/zrho ! K m s-1
            !- depth of 1st model layer
            !- (zo(1)-top is ~ 1/2 of the depth of 1st model layer, => mult by 2)
            pgeoh = 2.*(zo(i, 1) - topt(i, j))*c_grav ! m+2 s-2
            !-convective-scale velocity w*
            !- in the future, change 0.001 by ustar^3
            zws(i) = max(0., 0.001 - 1.5*0.41*zkhvfl*pgeoh/pten) ! m+3 s-3

            if (zws(i) > tiny(pgeoh)) then
               !-convective-scale velocity w*
               zws(i) = 1.2*zws(i)**.3333
               !- temperature excess
               ztexec(i) = max(0., -1.5*pahfs/(zrho*zws(i)*1004.64)) ! K
               !print*,"exce1=",pahfs,zrho,ztexec(i),zws(i),pgeoh,zo(i,1),topt(i,j)
               !call flush(6)
               !- moisture  excess
               zqexec(i) = max(0., -1.5*pqhfl/(zrho*zws(i)))        !kg kg-1
            end if   ! zws > 0
            !
            !- zws for shallow convection closure (Grant 2001)
            !- depth of the pbl
            pgeoh = pbl(i)*c_grav
            !-convective-scale velocity W* (m/s)
            zws(i) = max(0., 0.001 - 1.5*0.41*zkhvfl*pgeoh/pten)
            zws(i) = 1.2*zws(i)**.3333
         end do
         !
         !------ CALL CUMULUS PARAMETERIZATION
         !

         do ii_plume = 1, p_maxiens

            if (ii_plume == 1) then
               plume = p_shal
               c0 = C0_SHAL
            end if
            if (ii_plume == 2) then
               plume = p_deep
               c0 = C0_DEEP
            end if
            if (ii_plume == 3) then
               plume = p_mid
               c0 = C0_MID
            end if

            if (ICUMULUS_GF(plume) /= p_on) cycle

            hei_down_land = CUM_HEI_DOWN_LAND(plume)
            hei_down_ocean = CUM_HEI_DOWN_OCEAN(plume)
            hei_updf_land = CUM_HEI_UPDF_LAND(plume)
            hei_updf_ocean = CUM_HEI_UPDF_OCEAN(plume)
            max_edt_land = CUM_MAX_EDT_LAND(plume)
            max_edt_ocean = CUM_MAX_EDT_OCEAN(plume)
            fadj_massflx = CUM_FADJ_MASSFLX(plume)
            use_excess = CUM_USE_EXCESS(plume)
            ave_layer = CUM_AVE_LAYER(plume)
            t_star = CUM_T_STAR(plume)
            !print*,"plume=",plume,shal,mid,deep

            !-- set minimum/max for excess of T and Q
            if (use_excess == 0) then
               cum_ztexec(:) = 0.
               cum_zqexec(:) = 0.
            elseif (use_excess == 1) then
               cum_ztexec(:) = ztexec(:)
               cum_zqexec(:) = zqexec(:)
            elseif (use_excess == 2) then
               do i = its, itf
                  cum_zqexec(i) = min(5.e-4, max(1.e-4, zqexec(i)))! kg kg^-1
                  cum_ztexec(i) = min(0.5, max(0.2, ztexec(i)))! Kelvin
               end do
            else
               do i = its, itf
                  if (xlandi(i) > 0.98) then ! ocean
                     cum_zqexec(i) = min(8.e-4, max(5.e-4, zqexec(i)))! kg kg^-1
                     cum_ztexec(i) = min(1., max(0.5, ztexec(i)))! Kelvin
                  else                      ! land
                     cum_ztexec(i) = ztexec(i)
                     cum_zqexec(i) = zqexec(i)
                  end if
               end do
            end if
            !
            !--- replace 'q' and 't' excess in case of use of the cold pool scheme
            !
            if (CONVECTION_TRACER == 1 .and. plume == p_deep) then
               if (USE_GUSTINESS == 1) then
                  k = 2 ! surface in brams
                  do i = its, itf
                     cum_ztexec(i) = (hexcp(k, i, j) - qexcp(k, i, j)*real(c_alvl))/real(c_cp)
                     cum_zqexec(i) = qexcp(k, i, j)
                  end do
               else
                  cum_ztexec(:) = 0.
                  cum_zqexec(:) = 0.
               end if
            end if
            !
            !
            !-- shallow convection
            !
            if (plume == p_shal) then
               do i = its, itf
                  do k = kts, ktf
                     kr = k!+1 <<<<
                     if (p_use_gate) then
                        dhdt(i, k) = real(c_cp)*(temp_new_dp(i, k) - temp_old(i, k)) + real(c_alvl)*(qv_new_dp(i, k) - qv_old(i, k))
                        temp_new(i, k) = temp_new_dp(i, k)
                        qv_new(i, k) = qv_new_dp(i, k)
                     else

                        temp_new(i, k) = temp_old(i, k) + (rthblten(kr, i, j) + rthften(kr, i, j))*dt
                        qv_new(i, k) = qv_old(i, k) + (rqvblten(kr, i, j) + rqvften(kr, i, j))*dt
                        qv_new(i, k) = max(c_smaller_qv, qv_new(i, k))

                        !- only pbl forcing changes moist static energy
                        dhdt(i, k) = real(c_cp)*(rthblten(kr, i, j)) + real(c_alvl)*(rqvblten(kr, i, j))

                        !- all forcings change moist static energy
                        dhdt(i, k) = dhdt(i, k) + real(c_cp)*rthften(kr, i, j) + real(c_alvl)*rqvften(kr, i, j)

                     end if
                  end do
               end do
            end if
            !
            !--- deep convection
            if (plume == p_deep) then

               if (p_use_gate) then
                  do k = kts, ktf
                     do i = its, itf
                        temp_new(i, k) = temp_new_dp(i, k)
                        qv_new(i, k) = qv_new_dp(i, k)
                     end do
                  end do
               else
                  do k = kts, ktf
                     do i = its, itf
                        kr = k!+1 <<<<
                        temp_new(i, k) = temp_old(i, k) + (rthblten(kr, i, j) + rthften(kr, i, j))*dt
                        qv_new(i, k) = qv_old(i, k) + (rqvblten(kr, i, j) + rqvften(kr, i, j))*dt

                        temp_new_bl(i, k) = temp_old(i, k) + (rthblten(kr, i, j))*dt
                        qv_new_bl(i, k) = qv_old(i, k) + (rqvblten(kr, i, j))*dt
                     end do
                  end do
               end if
            end if
            !
            !--- mid/congestus type convection
            if (plume == p_mid) then

               if (p_use_gate) then
                  do k = kts, ktf
                     do i = its, itf
                        temp_new(i, k) = temp_new_dp(i, k)
                        qv_new(i, k) = qv_new_dp(i, k)
                     end do
                  end do
               else
                  do i = its, itf
                     do k = kts, ktf
                        kr = k!+1 <<<<

                        temp_new(i, k) = temp_old(i, k) + (rthblten(kr, i, j) + rthften(kr, i, j))*dt
                        qv_new(i, k) = qv_old(i, k) + (rqvblten(kr, i, j) + rqvften(kr, i, j))*dt
                        qv_new(i, k) = max(c_smaller_qv, qv_new(i, k))

                        !- only pbl forcing changes moist static energy
                        dhdt(i, k) = real(c_cp)*(rthblten(kr, i, j)) + real(c_alvl)*(rqvblten(kr, i, j))

                        !- all forcings change moist static energy
                        dhdt(i, k) = dhdt(i, k) + real(c_cp)*rthften(kr, i, j) + real(c_alvl)*rqvften(kr, i, j)

                        !- temp/water vapor modified only by bl processes
                        temp_new_bl(i, k) = temp_old(i, k) + (rthblten(kr, i, j))*dt
                        qv_new_bl(i, k) = qv_old(i, k) + (rqvblten(kr, i, j))*dt

                     end do
                  end do
               end if
            end if
            !

            call cupGf(its, ite, kts, kte, itf, ktf, mtp, nmp, fscav &
                        , p_cumulus_type(plume) &
                        , CLOSURE_CHOICE(plume) &
                        , CUM_ENTR_RATE(plume) &
                        , CUM_USE_EXCESS(plume) &
                        !- input data
                        , dx2d(:, j) &
                        , stochastic_sig(:, j) &
                        , col_sat(:, j) &
                        , tke_pbl(:, j) &
                        , rh_dicycle_fct(:, j) &
                        , wlpool(:, j) &
                        , dt &
                        , kpbli &
                        , cum_ztexec &
                        , cum_zqexec &
                        , ccn &
                        , rhoi &
                        , omeg &
                        , temp_old &
                        , qv_old &
                        , ter11 &
                        , h_sfc_flux &
                        , le_sfc_flux &
                        , xlons &
                        , xlats &
                        , xlandi &
                        , temp_new &
                        , qv_new &
                        , temp_new_bl &
                        , qv_new_bl &
                        , temp_new_adv &
                        , qv_new_adv &
                        , zo &
                        , po &
                        , tsur &
                        , psur &
                        , us &
                        , vs &
                        , dm2d &
                        , se_chem &
                        , zws &
                        , dhdt &
                        , buoy_exc2d &
                        , mpqi &
                        , mpql &
                        , mpcf &
                        , last_ierr(:) &
                        !output data
                        , outt(:, :, plume) &
                        , outq(:, :, plume) &
                        , outqc(:, :, plume) &
                        , outu(:, :, plume) &
                        , outv(:, :, plume) &
                        , outnliq(:, :, plume) &
                        , outnice(:, :, plume) &
                        , outbuoy(:, :, plume) &
                        , outmpqi(:, :, :, plume) &
                        , outmpql(:, :, :, plume) &
                        , outmpcf(:, :, :, plume) &
                        , out_chem(:, :, :, plume) &
                        !- for convective transport
                        , ierr4d(:, j, plume) &
                        , jmin4d(:, j, plume) &
                        , klcl4d(:, j, plume) &
                        , k224d(:, j, plume) &
                        , kbcon4d(:, j, plume) &
                        , ktop4d(:, j, plume) &
                        , kstabi4d(:, j, plume) &
                        , kstabm4d(:, j, plume) &
                        , cprr4d(:, j, plume) &
                        , xmb4d(:, j, plume) &
                        , edt4d(:, j, plume) &
                        , pwav4d(:, j, plume) &
                        , sigma4d(:, j, plume) &
                        , pcup5d(:, j, :, plume) &
                        , up_massentr5d(:, j, :, plume) &
                        , up_massdetr5d(:, j, :, plume) &
                        , dd_massentr5d(:, j, :, plume) &
                        , dd_massdetr5d(:, j, :, plume) &
                        , zup5d(:, j, :, plume) &
                        , zdn5d(:, j, :, plume) &
                        , prup5d(:, j, :, plume) &
                        , prdn5d(:, j, :, plume) &
                        , clwup5d(:, j, :, plume) &
                        , tup5d(:, j, :, plume) &
                        , conv_cld_fr5d(:, j, :, plume) &
                        !-- for debug/diag
                        , aa0(:, j), aa1(:, j), aa1_adv(:, j), a1_radpbl(:, j), aa2(:, j), aa3(:, j) &
                        , aa1_bl(:, j), aa1_cin(:, j), tau_bl(:, j), tau_ec(:, j) &
                        !-- for diag
                        , lightn_dens(:, j) &
                        , var2d(:, j) &
                        , revsu_gf_2d &
                        , prfil_gf_2d &
                        , var3d_agf_2d &
                        , var3d_bgf_2d &
                        , Tpert_2d &
                        )

         end do !- plume

         !--- reset ierr4d to value different of zero in case the correspondent
         !--- plume (shalllow, congestus, deep) was not actually used
         do n = 1, p_maxiens
            if (ICUMULUS_GF(n) == p_off) ierr4d(:, j, n) = -99
         end do

         do i = its, itf
            do_this_column(i, j) = 0
            loop1: do n = 1, p_maxiens
               if (ierr4d(i, j, n) == 0) then
                  do_this_column(i, j) = 1
                  exit loop1
               end if
            end do loop1
         end do
         !----------- check for negative water vapor mix ratio
         do i = its, itf
            if (do_this_column(i, j) == 0) cycle
            do k = kts, ktf
               temp_tendqv(i, k) = outq(i, k, p_shal) + outq(i, k, p_deep) + outq(i, k, p_mid)
            end do

            do k = kts, ktf
               distance(k) = qv_curr(i, k) + temp_tendqv(i, k)*dt
            end do

            if (minval(distance(kts:ktf)) < 0.) then
               zmax = minloc(distance(kts:ktf), 1)

               if (abs(temp_tendqv(i, zmax)*dt) < p_mintracer) then
                  fixout_qv(i) = 0.999999
                  !fixout_qv(i)= 0.
               else
                  fixout_qv(i) = ((c_smaller_qv - qv_curr(i, zmax)))/(temp_tendqv(i, zmax)*dt)
               end if
               fixout_qv(i) = max(0., min(fixout_qv(i), 1.))
            end if
         end do
         !------------ feedback
         !-- deep convection
         do i = its, itf
            if (do_this_column(i, j) == 0) cycle
            cprr4d(i, j, p_deep) = cprr4d(i, j, p_deep)*fixout_qv(i)
            cprr4d(i, j, p_mid) = cprr4d(i, j, p_mid)*fixout_qv(i)
            cprr4d(i, j, p_shal) = cprr4d(i, j, p_shal)*fixout_qv(i)
            conprr(i, j) = (cprr4d(i, j, p_deep) + cprr4d(i, j, p_mid) + cprr4d(i, j, p_shal))
            conprr(i, j) = max(0., conprr(i, j))
         end do

         !-- deep + shallow + mid convection
         do i = its, itf
            if (do_this_column(i, j) == 0) cycle
            do k = kts, kte
               kr = k!+1
               !- feedback the tendencies from convection
               rthcuten(kr, i, j) = (outt(i, k, p_shal) + outt(i, k, p_deep) + outt(i, k, p_mid))*fixout_qv(i)

               rqvcuten(kr, i, j) = (outq(i, k, p_shal) + outq(i, k, p_deep) + outq(i, k, p_mid))*fixout_qv(i)

               rqccuten(kr, i, j) = (outqc(i, k, p_shal) + outqc(i, k, p_deep) + outqc(i, k, p_mid))*fixout_qv(i)

               revsu_gf(kr, i, j) = revsu_gf_2d(i, k)*fixout_qv(i) !-- already contains deep and mid amounts.

               !---these arrays are only for the deep plume mode
               prfil_gf(kr, i, j) = prfil_gf_2d(i, k)*fixout_qv(i) !-- ice/liq prec flux of the deep plume
               !VAR3d_aGF(kr,i,j)= var3d_gf_2d(i,k)              !-- vertical velocity of the deep plume
               var3d_agf(kr, i, j) = outt(i, k, p_mid)*fixout_qv(i)   !--
               var3d_bgf(kr, i, j) = outq(i, k, p_mid)*fixout_qv(i)   !--

               if (ICUMULUS_GF(p_shal) == p_off) then
                  var3d_cgf(kr, i, j) = outqc(i, k, p_deep)*fixout_qv(i)  !--
                  var3d_dgf(kr, i, j) = outqc(i, k, p_mid)*fixout_qv(i)  !--
               else
                  var3d_cgf(kr, i, j) = outt(i, k, p_shal)*fixout_qv(i)   !--
                  var3d_dgf(kr, i, j) = outq(i, k, p_shal)*fixout_qv(i)   !--
               end if

            end do
         end do
         if (USE_MOMENTUM_TRANSP > 0) then
            do i = its, itf
               if (do_this_column(i, j) == 0) cycle
               do k = kts, kte
                  kr = k!+1
                  rucuten(kr, i, j) = (outu(i, k, p_deep) + outu(i, k, p_mid) + outu(i, k, p_shal))*fixout_qv(i)
                  rvcuten(kr, i, j) = (outv(i, k, p_deep) + outv(i, k, p_mid) + outv(i, k, p_shal))*fixout_qv(i)
               end do
            end do
         end if

         if (APPLY_SUB_MP == 1) then
            do i = its, itf
               if (do_this_column(i, j) == 0) cycle
               do k = kts, kte
                  kr = k!+1
                  sub_mpql(:, kr, i, j) = (outmpql(:, i, k, p_deep) + outmpql(:, i, k, p_mid) + outmpql(:, i, k, p_shal)) &
                                        * fixout_qv(i)
                  sub_mpqi(:, kr, i, j) = (outmpqi(:, i, k, p_deep) + outmpqi(:, i, k, p_mid) + outmpqi(:, i, k, p_shal)) &
                                        * fixout_qv(i)
                  sub_mpcf(:, kr, i, j) = (outmpcf(:, i, k, p_deep) + outmpcf(:, i, k, p_mid) + outmpcf(:, i, k, p_shal)) &
                                        * fixout_qv(i)
               end do
            end do
         end if

         if (LIQ_ICE_NUMBER_CONC == 1) then
            do i = its, itf
               if (do_this_column(i, j) == 0) cycle
               do k = kts, kte
                  kr = k!+1
                  rnicuten(kr, i, j) = (outnice(i, k, p_shal) + outnice(i, k, p_deep) + outnice(i, k, p_mid))*fixout_qv(i)
                  rnlcuten(kr, i, j) = (outnliq(i, k, p_shal) + outnliq(i, k, p_deep) + outnliq(i, k, p_mid))*fixout_qv(i)
               end do
            end do
         end if

         if (USE_TRACER_TRANSP == 1) then
            do i = its, itf
               if (do_this_column(i, j) == 0) cycle
               do k = kts, kte
                  kr = k!+1
                 rchemcuten(:, kr, i, j) = (out_chem(:, i, k, p_deep) + out_chem(:, i, k, p_mid) + out_chem(:, i, k, p_shal)) &
                                         * fixout_qv(i)
               end do
            end do

            !- constrain positivity for tracers
            do i = its, itf
               if (do_this_column(i, j) == 0) cycle

               do ispc = 1, mtp
                  if (chem_name_mask(ispc) == 0) cycle

                  do k = kts, ktf
                     distance(k) = se_chem(ispc, i, k) + rchemcuten(ispc, k, i, j)*dt
                  end do

                  !-- fixer for mass of tracer
                  if (minval(distance(kts:ktf)) < 0.) then
                     zmax = minloc(distance(kts:ktf), 1)

                     if (abs(rchemcuten(ispc, zmax, i, j)*dt) < p_mintracer) then
                        fixouts = 0.999999
                        !fixouts= 0.
                     else
                        fixouts = ((p_mintracer - se_chem(ispc, i, zmax)))/(rchemcuten(ispc, zmax, i, j)*dt)
                     end if
                     if (fixouts > 1. .or. fixouts < 0.) fixouts = 0.

                     rchemcuten(ispc, kts:ktf, i, j) = fixouts*rchemcuten(ispc, kts:ktf, i, j)
                  end if
               end do
            end do
         end if

         if (CONVECTION_TRACER == 1) then
            do i = its, itf
               if (do_this_column(i, j) == 0) cycle
               do k = kts, kte
                  kr = k!+1
                  rbuoycuten(kr, i, j) = (outbuoy(i, k, p_deep) + outbuoy(i, k, p_mid) + outbuoy(i, k, p_shal))*fixout_qv(i)
                  !print*,"RBUOYCUTEN", RBUOYCUTEN (kr,i,j),outbuoy(i,k,deep),&
                  !             outbuoy(i,k,mid),outbuoy(i,k,shal)
               end do
            end do

            !----- for output only
            !if(use_gustiness==1 .or. use_gustiness ==2 ) then
            !print*,"H-T",cp*1.1*maxval(sflux_t(:,j))&
            !            ,xlv*1.1*maxval(sflux_r(:,j)),maxval(ztexec),maxval(zqexec)
            !sflux_t(:,j) = ztexec(:)
            !sflux_r(:,j) = zqexec(:)
            !endif
         end if

         !-----memory
         !AA3(:,j)=cprr4d(:,j,deep) *fixout_qv(:)
         !AA2(:,j)=cprr4d(:,j,mid)  *fixout_qv(:)
      end do

   end subroutine modConvParGFDriver

   ! ---------------------------------------------------------------------------------------------------
   subroutine cupGF(its, ite, kts, kte, itf, ktf, mtp, nmp, fscav &
                     , cumulus &
                     , ichoice &
                     , entr_rate_input &
                     , use_excess &
                     !input data
                     , dx &
                     , stochastic_sig &
                     , col_sat &
                     , tke_pbl &
                     , rh_dicycle_fct &
                     , wlpool &
                     , dtime &
                     , kpbl &
                     , ztexec &
                     , zqexec &
                     , ccn &
                     , rho &
                     , omeg &
                     , t &
                     , q &
                     , z1 &
                     , h_sfc_flux &
                     , le_sfc_flux &
                     , xlons &
                     , xlats &
                     , xland &
                     , tn &
                     , qo &
                     , tn_bl &
                     , qo_bl &
                     , tn_adv &
                     , qo_adv &
                     , zo &
                     , po &
                     , tsur &
                     , psur &
                     , us &
                     , vs &
                     , dm2d &
                     , se_chem &
                     , zws &
                     , dhdt &
                     , buoy_exc &
                     , mpqi &
                     , mpql &
                     , mpcf &
                     , last_ierr &
                     !output data
                     , outt &
                     , outq &
                     , outqc &
                     , outu &
                     , outv &
                     , outnliq &
                     , outnice &
                     , outbuoy &
                     , outmpqi &
                     , outmpql &
                     , outmpcf &
                     , out_chem &
                     !- for convective transport
                     , ierr &
                     , jmin &
                     , klcl &
                     , k22 &
                     , kbcon &
                     , ktop &
                     , kstabi &
                     , kstabm &
                     , pre &
                     , xmb &
                     , edto &
                     , pwavo &
                     , sig &
                     , po_cup &
                     , up_massentro &
                     , up_massdetro &
                     , dd_massentro &
                     , dd_massdetro &
                     , zuo &
                     , zdo &
                     , pwo &
                     , pwdo &
                     , qrco &
                     , tup &
                     , clfrac &
                     !- for convective transport-end
                     !- for debug/diag
                     , aaa0_, aa1_, aa1_adv_, aa1_radpbl_, aa2_, aa3_, aa1_bl_, aa1_cin_, tau_bl_, tau_ec_ &
                     , lightn_dens &
                     , var2d &
                     , revsu_gf &
                     , prfil_gf &
                     , var3d_agf &
                     , var3d_bgf &
                     , tpert &
                     )
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: procedureName = 'cupGf' ! Subroutine name
      logical, parameter:: p_use_inv_layers = .true.
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts &
                           , kte, ichoice, use_excess  &
                           , mtp, nmp

      real, intent(in) :: col_sat(:)
      real, intent(in) :: stochastic_sig(:)
      real, intent(in) :: tke_pbl(:)
      !- atmos composition arrays
      real, intent(in) :: fscav(:)
      real, intent(in) :: dtime
      real, intent(in) :: entr_rate_input
      
      character(len=*), intent(in) :: cumulus

      integer, intent(inout) :: kpbl(:)
      integer, intent(inout) :: last_ierr(:)
      integer, intent(inout) :: ierr(:)
      integer, intent(inout) :: jmin(:)
      integer, intent(inout) :: klcl(:)
      integer, intent(inout) :: k22(:)
      integer, intent(inout) :: kbcon(:)
      integer, intent(inout) :: ktop(:)
      integer, intent(inout) :: kstabi(:)
      integer, intent(inout) :: kstabm(:)

      real, intent(inout) :: outu(:, :)
      real, intent(inout) :: outv(:, :)
      real, intent(inout) :: outt(:, :)
      !! output temp tendency (per s)
      real, intent(inout) :: outq(:, :)
      !! output q tendency (per s)
      real, intent(inout) :: outqc(:, :)
      !! output qc tendency (per s)
      real, intent(inout) :: outbuoy(:, :)
      real, intent(inout) :: revsu_gf(:, :)
      real, intent(inout) :: prfil_gf(:, :)
      real, intent(inout) :: var3d_agf(:, :)
      real, intent(inout) :: var3d_bgf(:, :)
      real, intent(inout) :: outnliq(:, :)
      real, intent(inout) :: outnice(:, :)
      ! basic environmental input includes
      ! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
      ! convection for this call only and at that particular gridpoint
      real, intent(inout) :: buoy_exc(:, :)
      real, intent(inout) :: ccn(:)
      real, intent(inout) :: clfrac(:, :)
      real, intent(inout) :: dd_massdetro(:, :)
      real, intent(inout) :: dd_massentro(:, :)
      real, intent(inout) :: dhdt(:, :)
      real, intent(inout) :: dm2d(:, :)
      real, intent(inout) :: dx(:)
      real, intent(inout) :: edto(:)
      real, intent(inout) :: h_sfc_flux(:)
      real, intent(inout) :: le_sfc_flux(:)
      real, intent(inout) :: mpcf(:, :, :)
      real, intent(inout) :: mpqi(:, :, :)
      real, intent(inout) :: mpql(:, :, :)
      real, intent(inout) :: omeg(:, :, :)
      real, intent(inout) :: out_chem(: , :, :)
      real, intent(inout) :: outmpcf(:, :, :)
      real, intent(inout) :: outmpqi(:, :, :)
      real, intent(inout) :: outmpql(:, :, :)
      real, intent(inout) :: po_cup(:, :)
      real, intent(inout) :: po(:, :)
      real, intent(inout) :: psur(:)
      real, intent(inout) :: pwavo(:)
      real, intent(inout) :: pwdo(:, :)
      real, intent(inout) :: pwo(:, :)
      real, intent(inout) :: q(:, :)
      !! environmental mixing ratio
      real, intent(inout) :: qo_adv(:, :)
      real, intent(inout) :: qo_bl(:, :)
      real, intent(inout) :: qo(:, :)
      real, intent(inout) :: qrco(:, :)
      real, intent(inout) :: rh_dicycle_fct(:)
      real, intent(inout) :: rho(:, :)
      real, intent(inout) :: se_chem(:, :, :)
      real, intent(inout) :: t(:, :)
      !! environmental temp
      real, intent(inout) :: tn_adv(:, :)
      real, intent(inout) :: tn_bl(:, :)
      real, intent(inout) :: tn(:, :)
      real, intent(inout) :: tpert(:, :)
      real, intent(inout) :: tsur(:)
      real, intent(inout) :: tup(:, :)
      real, intent(inout) :: up_massdetro(:, :)
      real, intent(inout) :: up_massentro(:, :)
      real, intent(inout) :: us(:, :)
      real, intent(inout) :: vs(:, :)
      real, intent(inout) :: wlpool(:)
      real, intent(inout) :: xland(:)
      real, intent(inout) :: xlats(:)
      real, intent(inout) :: xlons(:)
      real, intent(inout) :: xmb(:)
      real, intent(inout) :: z1(:)
      real, intent(inout) :: zdo(:, :)
      real, intent(inout) :: zqexec(:)
      real, intent(inout) :: ztexec(:)
      real, intent(inout) :: zuo(:, :)
      real, intent(inout) :: zws(:)
      !-- debug/diag
      real, intent(inout) :: aaa0_(:)
      real, intent(inout) :: aa1_(:)
      real, intent(inout) :: aa1_adv_(:)
      real, intent(inout) :: aa1_bl_(:)
      real, intent(inout) :: aa1_cin_(:)
      real, intent(inout) :: aa1_radpbl_(:)
      real, intent(inout) :: aa2_(:)
      real, intent(inout) :: aa3_(:)
      real, intent(inout) :: tau_bl_(:)
      real, intent(inout) :: tau_ec_(:)

      real, intent(out) :: lightn_dens(:)
      real, intent(out) :: pre(:)
      real, intent(out) :: sig(:)
      real, intent(out) :: var2d(:)

      !Local variables:
      ! local ensemble dependent variables in this routine
      real, dimension(its:ite, 1:p_maxens2) :: edtc
      real, dimension(its:ite, 1:p_ensdim) :: xf_ens, pr_ens
      !
      !*******the following are your basic environmental
      !          variables. They carry a "_cup" if they are
      !          on model cloud levels (staggered). They carry
      !          an "o"-ending (z becomes zo), if they are the forced
      !          variables. They are preceded by x (z becomes xz)
      !          to indicate modification by some typ of cloud
      !
      ! z           = heights of model levels
      ! qes         = environmental saturation mixing ratio
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

      integer, dimension(its:ite) :: kzdown, kdet, kb, ierr2 &
                                   , ierr3, kbmax, ierr_dummy &
                                   , start_level
      integer, dimension(its:ite, kts:kte) :: k_inv_layers

      integer :: iloop, nall, iedt, nens, nens3, ki, i, k, KK, iresult, nvar, nvarbegin
      integer :: jprnt, k1, k2, kbegzu, kdefi, kfinalzu, kstart, jmini, imid, k_free_trop
      integer :: iversion, step, fase, start_k22, ispc, kmp, istep, lstep, status
      integer :: X_kte, X_k, X_i, X_jcol

      integer :: i_wb = 0,ipr = 0, jpr = 0, bl = 1, fa = 2

      real, dimension(its:ite) :: edtx, hkb, hkbo, xhkb, qkb, pwevo &
                                , bu, bud, cap_max, xland1, vshear &
                                , cap_max_increment, psum, psumh, sigd &
                                , mconv, rescale_entrain, entr_rate &
                                , mentrd_rate, aa0_bl, aa1_bl, tau_bl &
                                , tau_ecmwf, wmean, aa1_fa, aa1_tmp, hkbo_x &
                                , aa2, aa3, cin0, cin1, edtmax, edtmin &
                                , aa1_lift, aa_tmp, aa_ini, aa_adv &
                                , daa_adv_dt, wlpool_bcon, xaa0_x, xk_x &
                                , xf_dicycle, mbdt, xf_coldpool, vvel1d &
                                , x_add_buoy, lambau_dn, lambau_dp &
                                , q_wetbulb, t_wetbulb, col_sat_adv &
                                , q_adv, alpha_adv, aa1_radpbl, aa1_adv &
                                , p_cwv_ave, cape, depth_neg_buoy, frh_bcon &
                                , check_sig, random, rh_entr_factor
      real, dimension(its:ite) :: aa0
      !! cloud work function without forcing effects
      real, dimension(its:ite) :: aa1
      !! cloud work function with forcing effects
      real, dimension(its:ite) :: xaa0
      !! cloud work function with cloud effects
      real, dimension(its:ite) :: edt
      !! epsilon
      real, dimension(kts:kte)   :: dummy1, dummy2
      real, dimension(its:ite, kts:kte) :: entr_rate_2d, mentrd_rate_2d &
                                         , he, hes, qes, z, heo, heso, qeso, zo, zu, zd &
                                         , xhe, xhes, xqes, xz, xt, xq &
                                         , qes_cup, q_cup, he_cup, hes_cup, z_cup, p_cup, gamma_cup, t_cup &
                                         , qeso_cup, qo_cup, heo_cup, heso_cup, zo_cup, gammao_cup, tn_cup &
                                         , xqes_cup, xq_cup, xhe_cup, xhes_cup, xz_cup &
                                         , xt_cup, hcot, evap_bcb &
                                         , dby, hc, clw_all &
                                         , dbyo, qco, qrcdo, hcdo, qcdo, dbydo, hco &
                                         , xdby, xzu, xzd, xhc, cupclw, pwo_eff &
                                         , dsubq, dsubh, dellabuoy, u_cup, v_cup, uc &
                                         , vc, ucd, vcd, dellu, dellv, up_massentr &
                                         , up_massdetr, dd_massentr, dd_massdetr &
                                         , subten_h, subten_q, subten_t
      real, dimension(its:ite, kts:kte) :: tn_x, qo_x, qeso_x, heo_x, heso_x, zo_cup_x &
                                         , qeso_cup_x, qo_cup_x, heo_cup_x, heso_cup_x &
                                         , po_cup_x, gammao_CUP_x, tn_cup_x, hco_x, DBYo_x &
                                         , u_cup_x, v_cup_x, xhe_x, xhes_x, xt_x, xq_x &
                                         , xqes_x, xqes_cup_x, xq_cup_x, xhe_cup_x &
                                         , xhes_cup_x, gamma_cup_x, xt_cup_x, dtempdz &
                                         , vvel2d, tempco, tempcdo, p_liq_ice, melting_layer &
                                         , melting, c1d, up_massentru, up_massdetru &
                                         , dd_massentru, dd_massdetru, prec_flx, evap_flx &
                                         , qrr, massflx, zenv, rho_hydr, dtdt, dqdt !alpha_h, alpha_q
      real, dimension(its:ite, kts:kte) :: cd
      !! detrainment function for updraft
      real, dimension(its:ite, kts:kte) :: cdd
      !! detrainment function for downdraft
      real, dimension(its:ite, kts:kte) :: dellah
      !! 
      real, dimension(its:ite, kts:kte) :: dellaq
      !! change of q per unit mass flux of cloud ensemble
      real, dimension(its:ite, kts:kte) :: dellat
      !! change of temperature per unit mass flux of cloud ensemble
      real, dimension(its:ite, kts:kte) :: dellaqc
      !! change of qc per unit mass flux of cloud ensemble
      real, dimension(its:ite, 1:p_maxens3) ::  xff_mid
      real, dimension(its:ite, p_shall_closures) :: xff_shal
      real, dimension(mtp, its:ite, kts:kte) :: se_cup_chem, sc_up_chem, sc_dn_chem &
                                              , pw_up_chem, pw_dn_chem
      real, dimension(mtp, its:ite) ::  tot_pw_up_chem, tot_pw_dn_chem
      real, dimension(mtp, kts:kte) :: trcflx_in, sub_tend, ddtr, ddtr_upd, zenv_diff &
                                     , fp_mtp, fm_mtp
      real, dimension(mtp) :: min_tend_chem, dummy_chem, delS_up, delS_dn, env_sub &
                            , outchem1, evap_, wetdep_, trash_, trash2_, residu_
      real, dimension(nmp, its:ite, kts:kte) :: dellampqi, dellampql, dellampcf
      real, dimension(kts:kte) ::  aa, bb, cc, ddu, ddv, ddh, ddq, fp, fm
      real, dimension(nmp, kts:kte) ::  dd
      real, dimension(its:ite) :: x_dx, x_stochastic_sig
      real, dimension(kts:kte, 8) ::  tend2d
      real, dimension(8) ::  tend1d
      real, dimension(its:ite, 8) ::  check_cons_i, check_cons_f
      
            !-- only for debug (atmos composition)
      real, allocatable, dimension(:, :, :), save    ::   se_chem_update

      real ::  massi, massf, dtime_max, evap, wetdep, umean
      real :: day, dz, dzo, radius, entrd_rate, zcutdown, depth_min, zkbmax &
            , z_detr, zktop, massfld, dh, trash, frh, xlamdd, radiusd, frhd &
            , effec_entrain, detdo1, detdo2, entdo, dp, subin, detdo, entup &
            , detup, subdown, entdoj, entupk, detupk, totmas, min_entr_rate
      real :: tot_time_hr, beta, env_mf, env_mf_p, env_mf_m, ts, denom, denom_u
      real :: dsubh_aver, dellah_aver, x_add, cap_max_inc, tlll, plll, rlll &
            , tlcl, plcl, dzlcl, zlll, c_up, e_dn, g_rain, trash2, pgc, bl2dp &
            , trash3, ke_mx, s1, s2, q1, q2, rz_env, factor, cwv, entr_threshold &
            , resten_h, resten_q, resten_t, min_deep_top, min_shall_top
      real :: alp0, beta1, beta2, dp_p, dp_m, delt1, delt2, delt_tvv, rcount

      logical :: keep_going

      character(len=128) :: ierrc(its:ite)
      character(len=2) :: cty
      character(len=200) :: lixo

      !----------------------------------------------------------------------
      !--only for debug
      if (p_use_gate) then
         if (.not. allocated(se_chem_update)) allocate (se_chem_update(3, its:ite, kts:kte))
         if (jl == 1) then
            !    se_chem_update(1,:,:) = mpql (lsmp,:,:)
            !    se_chem_update(2,:,:) = mpqi (lsmp,:,:)
            !    se_chem_update(3,:,:) = mpcf (lsmp,:,:)
         else
            !    mpql (lsmp,:,:)= se_chem_update(1,:,:)
            !    mpqi (lsmp,:,:)= se_chem_update(2,:,:)
            !    mpcf (lsmp,:,:)= se_chem_update(3,:,:)
         end if
      end if

      call initCupGF(its, itf, ite, ktf, max_edt_ocean, max_edt_land, xland, zo, cumulus, random, depth_min, zkbmax, z_detr &
                  ,  zcutdown, kstabm, kbmax, ierr2 , ierr3, lambau_dp, lambau_dn, edtmin, edtmax, xland1, cap_max &
                  ,  aa0, aa1, aa2, aa3, aa1_bl, aa1_fa, aa0_bl, q_adv, aa1_radpbl, aa1_adv, alpha_adv, cin1, xk_x &
                  ,  edt , edto, tau_bl, q_wetbulb, t_wetbulb, tau_ecmwf, xf_dicycle, x_add_buoy, xf_coldpool &
                  ,  wlpool_bcon, cap_max_increment, mbdt, z, xz, hcdo, cupclw, qrcdo, hcot, c1d, xf_ens, pr_ens &
                  , evap_bcb, ierrc)

      ! --- environmental conditions, FIRST HEIGHTS
      ! --- calculate moist static energy, heights, qes
      call cupEnv(z, qes, he, hes, t, q, po, z1, psur, ierr, -1, itf, ktf, its, ite, kts, kte)
      call cupEnv(zo, qeso, heo, heso, tn, qo, po, z1, psur, ierr, -1, itf, ktf, its, ite, kts, kte)

      ! --- outputs a model sounding for the stand-alone code (part 1)
      if (output_sound == 1) then
         call sound(1, cumulus, int_time, dtime, p_ens4, itf, ktf, its, ite, kts, kte, xlats, xlons, jcol, whoami_all &
                    , z, qes, he, hes, t, q, po, z1, psur, zo, qeso, heo, heso, tn, qo, us, vs, omeg, xz &
                    , h_sfc_flux, le_sfc_flux, tsur, dx, stochastic_sig, zws, ztexec, zqexec, xland &
                    , kpbl, k22, klcl, kbcon, ktop, aa0, aa1, sig, xaa0, hkb, xmb, pre, edto &
                    , zo_cup, dhdt, rho, zuo, zdo, up_massentro, up_massdetro, outt, outq, outqc, outu, outv)
      end if

      ! --- environmental values on cloud levels
      call cupEnvCLev(t, qes, q, he, hes, z, po, qes_cup, q_cup, he_cup, us, vs, u_cup, v_cup, hes_cup, z_cup, p_cup &
                    , gamma_cup, t_cup, psur, tsur, ierr, z1, itf, ktf, its, ite, kts, kte)

      call cupEnvCLev(tn, qeso, qo, heo, heso, zo, po, qeso_cup, qo_cup, heo_cup, us, vs, u_cup, v_cup, heso_cup, zo_cup &
                    , po_cup, gammao_cup, tn_cup, psur, tsur, ierr, z1, itf, ktf, its, ite, kts, kte)

      ! --- get air density at full layer (model levels) by hydrostatic balance (kg/m3)
      do i = its, itf
         rho_hydr(i, :) = 0.0
         if (ierr(i) /= 0) cycle
         do k = kts, ktf
            rho_hydr(i, k) = 100.*(po_cup(i, k) - po_cup(i, k + 1))/(zo_cup(i, k + 1) - zo_cup(i, k))/c_grav
            !print*,"rhohidr=",k,rho_hydr(i,k),po_cup(i,k+1),zo_cup(i,k+1)
         end do
      end do

      ! --- partition between liq/ice cloud contents
      call getPartitionLiqIce(ierr, tn, z1, zo_cup, po_cup, p_liq_ice, melting_layer, itf, ktf, its, ite, kts, kte, cumulus)

      do i = its, itf
         if (ierr(i) /= 0) cycle
         do k = kts, ktf
            if (zo_cup(i, k) .gt. zkbmax + z1(i)) then
               kbmax(i) = k
               exit
            end if
         end do
         !--- level where detrainment for downdraft starts
         do k = kts, ktf
            if (zo_cup(i, k) .gt. z_detr + z1(i)) then
               kdet(i) = k
               exit
            end if
         end do
      end do

      ! --- determine level with highest moist static energy content - k22
      if (trim(cumulus) == 'shallow') then
         start_k22 = 1
      else
         start_k22 = 2
      end if
      k22(:) = kts
      do i = its, itf
         if (ierr(i) /= 0) cycle
         k22(i) = maxloc(heo_cup(i, start_k22:kbmax(i) + 1), 1) + start_k22 - 1
         k22(i) = max(k22(i), start_k22)
         if (trim(cumulus) == 'shallow') then
            k22(i) = min(2, k22(i))

            if (K22(i) .gt. kbmax(i)) then
               ierr(i) = 2
               ierrc(i) = "could not find k22"
            end if
         else
            if (k22(i) > kbmax(i)) then
               !- let's try k22=start_k22 for the cases k22>kbmax
               k22(i) = start_k22
               cycle
            end if

         end if
      end do

      !-- get the pickup of ensemble ave prec, following Neelin et al 2009.
      call precipCwvFactor(itf, ktf, its, ite, kts, kte, ierr, tn, po, qo, po_cup, cumulus, p_cwv_ave)

      !------- determine LCL for the air parcels around K22
      do i = its, itf
         klcl(i) = k22(i) ! default value
         if (ierr(i) == 0) then
            !tlll, rlll,plll - temp, water vapor and pressure of the source air parcel
            x_add = max(0., zqexec(i))
            call getCloudBc(cumulus, kts, kte, ktf, xland(i), po(i, kts:kte), q_cup(i, kts:kte), rlll, k22(i), x_add)
            x_add = max(0., ztexec(i))
            call getCloudBc(cumulus, kts, kte, ktf, xland(i), po(i, kts:kte), t_cup(i, kts:kte), tlll, k22(i), x_add)
            call getCloudBc(cumulus, kts, kte, ktf, xland(i), po(i, kts:kte), p_cup(i, kts:kte), plll, k22(i))
            !-get LCL
            call getLcl(tlll, 100.*plll, rlll, tlcl, plcl, dzlcl)

            if (dzlcl >= 0.) then ! LCL found (if dzlcl<0 => not found)
               call getCloudBc(cumulus, kts, kte, ktf, xland(i), po(i, kts:kte), z_cup(i, kts:kte), zlll, k22(i))
               do k = kts, ktf
                  if (z_cup(i, k) .gt. zlll + dzlcl) then
                     klcl(i) = max(k, k22(i))
                     exit
                  end if
               end do
               klcl(i) = min(klcl(i), ktf - 4)
            end if
         end if
         !write(12,111)'MDlcl',tlcl,plcl,dzlcl,klcl(i),ierr(i)
         !111      format(1x,A5,3F10.2,2i4)
      end do

      !-- check if LCL height is below PBL height to allow shallow convection
      if (LCL_TRIGGER > 0 .and. trim(cumulus) == 'shallow') then
         do i = its, itf
            if (ierr(i) /= 0) cycle
            if (klcl(i) > max(1, kpbl(i) - LCL_TRIGGER)) then
               ierr(i) = 21
               ierrc(i) = 'for shallow convection:  LCL height < PBL height'
            end if
         end do
         !print*,"LCL",maxval(klcl),minval(klcl),maxval(kpbl),minval(kpbl)
      end if

      !--- define entrainment/detrainment profiles for updrafts
      !- initial entrainment/detrainment
      entr_rate(:) = entr_rate_input
      min_entr_rate = entr_rate_input*0.1

      !-- cold pool parameterization and convective memory
      if (CONVECTION_TRACER == 1 .and. trim(cumulus) == 'deep') then
         if (USE_MEMORY >= 0) then
            do i = its, itf
               if (ierr(i) /= 0) cycle
               !x_add_buoy(i) = min(mx_buoy2, maxval(buoy_exc(i,kts:klcl(i))))

               call getCloudBc(trim(cumulus), kts, kte, ktf, xland(i), po(i, kts:kte), buoy_exc(i, kts:kte), x_add_buoy(i), kts)
               ! buoy_exc (i,kts:kte),x_add_buoy (i),klcl(i))

            end do
            !print*,"BU=",maxval(x_add_buoy),minval(x_add_buoy)
         end if

         !-- avoid extra-buoyancy where rained before
         if (USE_MEMORY == 1 .or. USE_MEMORY == 12) then
            where (aa2_ > 10./3600.)
               x_add_buoy = 0.0
               wlpool = 0.0
            end where
         end if
         !-- avoid extra-buoyancy where rained before
         if (USE_MEMORY == 4 .or. USE_MEMORY == 14) then
            do i = its, itf
               if (ierr(i) /= 0) cycle
               if (aa2_(i) > 1.e-6 .and. x_add_buoy(i) < 1000. .and. x_add_buoy(i) > 250.) then
                  x_add_buoy(i) = 0.0
                  wlpool(i) = 0.0
                  ierr(i) = 100
               end if
            end do
            !where(AA2_ > 1.e-6. .and. x_add_buoy < 1000.)
            !   x_add_buoy = 0.0
            !   wlpool     = 0.0
            !   ierr       = 100
            !end where
         end if

         if (USE_MEMORY == 2 .or. USE_MEMORY == 12 .or. USE_MEMORY == 14) then
            !- initial entrainment/detrainment
            entr_rate(:) = entr_rate_input       !* 2.0
            min_entr_rate = entr_rate_input*0.1 !* 2.0

            do i = its, itf !-- reduce entr rate, where cold pools exist
               if (ierr(i) /= 0) cycle
               !entr_rate(i) = max(0.1, 1.-ColdPoolStart(x_add_buoy(i))) * entr_rate(i)
               !entr_rate(i) = max(0.5, 1.-ColdPoolStart(x_add_buoy(i))) * entr_rate(i)
               entr_rate(i) = max(0.7, 1.-ColdPoolStart(x_add_buoy(i)))*entr_rate(i)
               !entr_rate(i) = max(0.8, 1.-ColdPoolStart(x_add_buoy(i))) * entr_rate(i)
            end do
            !print*,"ENT",1000.*maxval(entr_rate),1000.*minval(entr_rate)&
            !            ,ColdPoolStart(maxval((x_add_buoy(:)))),ColdPoolStart(minval((x_add_buoy(:))))
         end if

         if (USE_MEMORY == 3 .or. ADD_COLDPOOL_CLOS >= 1) then ! increase capmax
            do i = its, itf
               if (ierr(i) /= 0) cycle
               cap_max(i) = cap_max(i) + ColdPoolStart(x_add_buoy(i))*35.
            end do
         end if
         if (ADD_COLDPOOL_CLOS == 3) then ! increase x_add_buoy
            do i = its, itf
               if (ierr(i) /= 0) cycle
               x_add_buoy(i) = x_add_buoy(i) + 0.5*wlpool(i)**2
            end do
         end if

         !--- temporary for output
         if (USE_GUSTINESS == 0) aa3_(:) = x_add_buoy(:)
         if (USE_GUSTINESS == 4) then
            aa3_(:) = x_add_buoy(:); x_add_buoy(:) = 0.0
         end if

         !-- using ztexc and zqexc as perturbation:
         if (USE_GUSTINESS == 1 .or. USE_GUSTINESS == 2) then
            aa3_(:) = real(c_cp)*ztexec(:) + real(c_alvl)*zqexec(:)
            x_add_buoy(:) = 0.
         end if

      end if
      !
      !--- determine the entrainment dependent on environmental moist (here relative humidity)
      !--- also the controls of RH on the diurnal cycle (see Tian et al 2022 GRL)
      if (trim(cumulus) == 'deep') call rhControls(itf, ktf, its, ite, kts, kte, ierr, tn, po, qo, qeso, po_cup, cumulus &
                                                 , rh_entr_factor, rh_dicycle_fct, entr_rate_input, entr_rate, xlons, dtime)

      !--- determine the vertical entrainment/detrainment rates, the level of convective cloud base -kbcon-
      !--- and the scale dependence factor (sig).
      do i = its, itf
         entr_rate_2d(i, :) = entr_rate(i)
         cd(i, :) = entr_rate(i)
         if (ierr(i) /= 0) cycle

         if (trim(cumulus) /= 'shallow') then
            do k = kts, ktf
               frh = min(qo_cup(i, k)/qeso_cup(i, k), 1.)
               !-------------------------------------------
               if (p_entr_new) then
                  !- v 2
                  if (k >= klcl(i)) then
                     !entr_rate_2d(i,k)=entr_rate(i)*(1.3-frh)*(qeso_cup(i,k)/qeso_cup(i,klcl(i)))**3
                     entr_rate_2d(i, k) = entr_rate(i)*(1.3 - frh)*(qeso_cup(i, k)/qeso_cup(i, klcl(i)))**1.25
                  else
                     entr_rate_2d(i, k) = entr_rate(i)*(1.3 - frh)
                  end if
                  cd(i, k) = 0.75e-4*(1.6 - frh)
                  entr_rate_2d(i, k) = max(entr_rate_2d(i, k), min_entr_rate)
               else
                  !- v 1
                  entr_rate_2d(i, k) = max(entr_rate(i)*(1.3 - frh)*max(min(1., (qeso_cup(i, k) &
                                     / qeso_cup(i, klcl(i)))**1.25), 0.1), 1.e-5)
                  if (trim(cumulus) == 'deep') cd(i, k) = 1.e-2*entr_rate(i)
                  if (trim(cumulus) == 'mid') cd(i, k) = 0.75*entr_rate_2d(i, k)
               end if
            end do
         else
            do k = kts, ktf
               frh = min(qo_cup(i, k)/qeso_cup(i, k), 1.)
               !entr_rate_2d(i,k)=entr_rate(i)*(1.3-frh)*max(min(1.,(qeso_cup(i,max(k,klcl(i)))&
               !                                                    /qeso_cup(i,klcl(i)))**3) ,0.1)
               entr_rate_2d(i, k) = entr_rate(i)*(1.3 - frh)*max(min(1., (qeso_cup(i, max(k, klcl(i))) &
                                  / qeso_cup(i, klcl(i)))**1), 0.1)

               ! entr_rate_2d(i,k)=entr_rate(i)*(1.3-frh)*(min(z(i,klcl(i))/z(i,k),1.))
               ! entr_rate_2d(i,k) = max(entr_rate_2d(i,k),min_entr_rate)
               !print*,"ent=",k,real(z(i,k),4),real(min(z(i,klcl(i))/z(i,k),1.),4),real(entr_rate_2d(i,k)*1000.,4)

               cd(i, k) = 0.75*entr_rate_2d(i, k)!+0.5e-3
            end do
         end if
      end do

      !--- start_level
      start_level(:) = klcl(:)
      !start_level(:)=  KTS

      !--- determine the moist static energy of air parcels at source level
      do i = its, itf
         if (ierr(i) /= 0) cycle
         x_add = (real(c_alvl)*zqexec(i) + real(c_cp)*ztexec(i)) + x_add_buoy(i)
         call getCloudBc(cumulus, kts, kte, ktf, xland(i), po(i, kts:kte), he_cup(i, kts:kte), hkb(i), k22(i), x_add &
                       , tpert(i, kts:kte))
         call getCloudBc(cumulus, kts, kte, ktf, xland(i), po(i, kts:kte), heo_cup(i, kts:kte), hkbo(i), k22(i), x_add &
                       , tpert(i, kts:kte))
         !print*,"xi=",i, xlv*zqexec(i)+cp*ztexec(i) , x_add_buoy(i), hkbo(i)
      end do
      !
      !--- determine the level of convective cloud base  - kbcon
      call cupCloudLimits(cumulus, ierrc, ierr, cap_max_increment, cap_max, heo_cup, heso_cup, qo_cup &
                            , qeso_cup, po, po_cup, zo_cup, heo, hkbo, qo, qeso, entr_rate_2d, hcot, k22, kbmax &
                            , klcl, kbcon, ktop, depth_neg_buoy, frh_bcon, tpert, start_level &
                            , use_excess, zqexec, ztexec, x_add_buoy, xland, itf, ktf, its, ite, kts, kte)

      !--- scale dependence factor (sig), version new
      if (USE_SCALE_DEP == 0 .or. trim(cumulus) == 'shallow') then
         sig(:) = 1.
      else
         do i = its, itf
            sig(i) = 0.
            if (ierr(i) /= 0) cycle
            !--original
            !sig(i) = 1.0-0.9839*exp(-0.09835.  *(dx(i)/1000.))
            !-- for similar curve as in IFS/EC, use sig_factor = 0.22
            sig(i) = 1.0 - exp(-SIG_FACTOR*(dx(i)/1000.))
            !print*,"sig=",sig(i),dx(i), sig_factor
            if (stochastic_sig(i) /= 1.0) then
               sig(i) = sig(i)**(stochastic_sig(i)*max(0.9, 0.9*sig(i)))
            end if
            sig(i) = max(0.001, min(sig(i), 1.))
         end do
         !print*,'sig',maxval(sig),minval(sig),maxval(dx),minval(dx)
      end if

      !--- define entrainment/detrainment profiles for downdrafts
      if (p_entr_new) then
         mentrd_rate(:) = entr_rate(:)*0.3
      else
         mentrd_rate(:) = entr_rate(:)
      end if
      do i = its, itf
         cdd(i, kts:kte) = mentrd_rate(i)
      end do
      !- scale dependence factor
      sigd(:) = 1.
      if (DOWNDRAFT == 0) sigd(:) = 0.0

      !--- update hkb/hkbo in case of k22 is redefined in 'cup_kbon'
      do i = its, itf
         if (ierr(i) /= 0) cycle
         x_add = (real(c_alvl)*zqexec(i) + real(c_cp)*ztexec(i)) + x_add_buoy(i)
         call getCloudBc(cumulus, kts, kte, ktf, xland(i), po(i, kts:kte), he_cup(i, kts:kte), hkb(i), k22(i), x_add &
                         , tpert(i, kts:kte))
         call getCloudBc(cumulus, kts, kte, ktf, xland(i), po(i, kts:kte), heo_cup(i, kts:kte), hkbo(i), k22(i), x_add &
                         , tpert(i, kts:kte))
      end do

      !--- increase detrainment in stable layers
      call cupMinimi(HEso_cup, Kbcon, kstabm, kstabi, ierr, itf, ktf, its, ite, kts, kte)

      !--- option for using the inversion layers as a barrier for the convection development
      if (trim(cumulus) == 'mid') then
         if (p_use_inv_layers) then
            !--- get inversion layers
            call getInversionLayers(cumulus, ierr, psur, po_cup, tn_cup, zo_cup, k_inv_layers, dtempdz, itf, ktf, its, ite &
                                  , kts, kte)
            do i = its, itf
               if (ierr(i) /= 0) cycle
               ktop(i) = min(ktop(i), k_inv_layers(i, p_mid))
               !print*,"ktop=",ktop(i),k_inv_layers(i,mid)
            end do
         end if

         !-- check if ktop is above 450hPa layer for mid convection
         do i = its, itf
            if (ierr(i) /= 0) cycle
            !print*,"sta=",Kbcon(i),kstabm(i),kstabi(i),p_cup(i,ktop(i)),z_cup(i,kstabi(i))
            if (po_cup(i, ktop(i)) < 450.) then
               ierr(i) = 25
               ierrc(i) = 'mid convection with cloud top above 450 hPa (~ 7km asl)'
            end if
         end do

         !-- check if ktop is below 750hPa layer for mid convection
         do i = its, itf
            if (ierr(i) /= 0) cycle
            if (po_cup(i, ktop(i)) > 750.) then
               ierr(i) = 55
               ierrc(i) = 'ktop too low for mid'
            end if
         end do
      end if

      if (trim(cumulus) == 'shallow') then
         if (p_use_inv_layers) then
            call getInversionLayers(cumulus, ierr, psur, po_cup, tn_cup, zo_cup, k_inv_layers, dtempdz, itf, ktf, its, ite, kts &
                                  , kte)
            do i = its, itf
               if (ierr(i) /= 0) cycle
               ktop(i) = min(ktop(i), k_inv_layers(i, p_shal))
            end do
         end if

         !--- Check if ktop is above 700hPa layer for shallow convection
         do i = its, itf
            if (ierr(i) /= 0) cycle
            min_shall_top = 700.
            !if(icumulus_gf(mid) == 0) min_shall_top=500.
            if (po_cup(i, ktop(i)) < min_shall_top) then
               ierr(i) = 26
               ierrc(i) = 'shallow convection wit h cloud top above min_shall_top hPa'
            end if
         end do
      end if

      do i = its, itf
         if (ktop(i) <= kbcon(i)) then
            ierr(i) = 5
            ierrc(i) = 'ktop too small'
         end if
      end do

      if (trim(cumulus) == 'deep') then
         min_deep_top = 500.
         if (ICUMULUS_GF(p_mid) == 0) min_deep_top = 750.
         do i = its, itf
            if (ierr(i) /= 0) cycle
            if (po_cup(i, ktop(i)) > min_deep_top) then
               ierr(i) = 55
               ierrc(i) = 'ktop too low for deep'
            end if
         end do
      end if

      !-- avoid double-counting with shallow scheme (deep and mid)
      do i = its, itf
         if (ierr(i) /= 0) cycle
         if (last_ierr(i) == 0) then
            !--- if 'mid' => last was 'shallow'
            ! if(cumulus == 'mid' .and. po_cup(i,ktop(i)) > 700.) then
            !   ierr(i)=27
            !   ierrc(i)='avoiding double-counting shallow and mid'
            ! endif
            !--- if 'mid' => last was 'shallow'
            if (trim(cumulus) == 'mid') then
               ierr(i) = 27
               ierrc(i) = 'avoiding double-counting deep and mid'
            end if
         end if
      end do

      !--- determine the normalized mass flux profile for updraft
      do i = its, itf
         zuo(i, :) = 0.
         if (ierr(i) /= 0) cycle
         call getZuZdPdf(trim(cumulus), trim(cumulus)//"_up", ierr(i), k22(i), ktop(i), zuo(i, kts:kte), kts, kte, ktf &
                     ,   kpbl(i), k22(i), kbcon(i), klcl(i), po_cup(i, kts:kte), psur(i), xland(i), random(i))
      end do

      do i = its, itf
         if (ierr(i) /= 0) cycle
         xzu(i, :) = zuo(i, :)
         zu(i, :) = zuo(i, :)
      end do

      ! calculate mass entrainment and detrainment
      call getLateralMassFlux(itf, ktf, its, ite, kts, kte, min_entr_rate, ierr, ktop, zo_cup, zuo, cd, entr_rate_2d, po_cup &
                                , up_massentro, up_massdetro, up_massentr, up_massdetr, cumulus, kbcon, k22, kpbl, up_massentru &
                                , up_massdetru, lambau_dp)
      uc = 0.
      vc = 0.
      hc = 0.
      hco = 0.
      do i = its, itf
         if (ierr(i) .eq. 0) then
            do k = kts, start_level(i)
               hc(i, k) = hkb(i)
               hco(i, k) = hkbo(i)
               !-get uc and vc as average between layers below k22
               call getCloudBc(cumulus, kts, kte, ktf, xland(i), po(i, kts:kte), u_cup(i, kts:kte), uc(i, k), k22(i))
               call getCloudBc(cumulus, kts, kte, ktf, xland(i), po(i, kts:kte), v_cup(i, kts:kte), vc(i, k), k22(i))
            end do
         end if
      end do

      !--- 1st guess for moist static energy and dbyo (not including ice phase)
      do i = its, itf
         if (ierr(i) /= 0) cycle
         do k = start_level(i) + 1, ktop(i) + 1  ! mass cons option
            denom = (zu(i, k - 1) - .5*up_massdetro(i, k - 1) + up_massentro(i, k - 1))
            if (denom > 0.0) then
               hco(i, k) = (hco(i, k - 1)*zuo(i, k - 1) - .5*up_massdetro(i, k - 1)*hco(i, k - 1) + up_massentro(i, k - 1) &
                         * heo(i, k - 1))/denom
               if (k == start_level(i) + 1) then
                  x_add = (real(c_alvl)*zqexec(i) + real(c_cp)*ztexec(i)) + x_add_buoy(i)
                  hco(i, k) = hco(i, k) + x_add*up_massentro(i, k - 1)/denom
               end if
            else
               hco(i, k) = hco(i, k - 1)
            end if
         end do
         do k = ktop(i) + 2, ktf
            hco(i, k) = heso_cup(i, k)!=heo_cup(i,k)
         end do
      end do

      !--- Get buoyancy of updrafts
      call getBuoyancy(itf, ktf, its, ite, kts, kte, ierr, klcl, kbcon, ktop, hco, heo_cup, heso_cup, dbyo, zo_cup)

      !--- get "c1d" profile ----------------------------------------
      if (trim(cumulus) == 'deep' .and. use_c1d) then
         do i = its, itf
            if (ierr(i) /= 0) cycle
            c1d(i, kbcon(i) + 1:ktop(i) - 1) = abs(C1)
         end do
      end if

      if (first_guess_w .or. AUTOCONV == 4) then
         call cupUpMoistureLight(cumulus, start_level, klcl, ierr, ierrc, zo_cup, qco, qrco, pwo, pwavo, hco, tempco, xland &
                                    , po, p_cup, kbcon, ktop, cd, dbyo, clw_all, t_cup, qo, gammao_cup, zuo &
                                    , qeso_cup, k22, qo_cup, zqexec, use_excess, rho, up_massentr, up_massdetr &
                                    , psum, psumh, c1d, x_add_buoy, 1, itf, ktf, ipr, jpr, its, ite, kts, kte)

         call cupUpVVel(vvel2d, vvel1d, zws, entr_rate_2d, cd, zo, zo_cup, zuo, dbyo, gammao_cup, tn_cup, tempco, qco, qrco, qo  &
                      , start_level, klcl, kbcon, ktop, ierr, itf, ktf, its, ite, kts, kte, wlpool, wlpool_bcon, 1)
      end if

      !--- calculate moisture properties of updraft
      call cupUpMoisture(cumulus, start_level, klcl, ierr, ierrc, zo_cup, qco, qrco, pwo, pwavo, hco, tempco, xland &
                           , po, p_cup, kbcon, ktop, cd, dbyo, clw_all, t_cup, qo, gammao_cup, zuo, qeso_cup &
                           , k22, qo_cup, zqexec, use_excess, ccn, rho, up_massentr, up_massdetr, psum &
                           , psumh, c1d, x_add_buoy, vvel2d, vvel1d, zws, entr_rate_2d &
                           , 1, itf, ktf, ipr, jpr, its, ite, kts, kte)

      do i = its, itf
         if (ierr(i) /= 0) cycle
         cupclw(i, kts:ktop(i) + 1) = qrco(i, kts:ktop(i) + 1)
      end do

      !--- get melting profile
      call getMeltingProfile(ierr, tn_cup, po_cup, p_liq_ice, melting_layer, qrco, pwo, edto, pwdo, melting, itf, ktf, its, ite &
                           , kts, kte, cumulus)

      !--- updraft moist static energy + momentum budget
      !--- option to produce linear fluxes in the sub-cloud layer.
      if (trim(cumulus) == 'shallow' .and. USE_LINEAR_SUBCL_MF == 1) then
         do i = its, itf
            if (ierr(i) /= 0) cycle
            call getDelmix(cumulus, kts, kte, ktf, xland(i), start_level(i), po(i, kts:kte), he_cup(i, kts:kte), hc(i, kts:kte))
            call getDelmix(cumulus, kts, kte, ktf, xland(i), start_level(i), po(i, kts:kte) , heo_cup(i, kts:kte), hco(i, kts:kte))
         end do
      end if

      do i = its, itf
         if (ierr(i) /= 0) cycle

         do k = start_level(i) + 1, ktop(i) + 1  ! mass cons option
            denom = (zu(i, k - 1) - .5*up_massdetr(i, k - 1) + up_massentr(i, k - 1))
            denom_u = (zu(i, k - 1) - .5*up_massdetru(i, k - 1) + up_massentru(i, k - 1))
            if (denom > 0.0 .and. denom_u > 0.0) then
               hc(i, k) = (hc(i, k - 1)*zu(i, k - 1) - .5*up_massdetr(i, k - 1)*hc(i, k - 1) + &
                           up_massentr(i, k - 1)*he(i, k - 1))/denom
               hco(i, k) = (hco(i, k - 1)*zuo(i, k - 1) - .5*up_massdetro(i, k - 1)*hco(i, k - 1) + &
                            up_massentro(i, k - 1)*heo(i, k - 1))/denom
               if (k == start_level(i) + 1) then
                  x_add = (real(c_alvl)*zqexec(i) + real(c_cp)*ztexec(i)) + x_add_buoy(i)
                  hco(i, k) = hco(i, k) + x_add*up_massentro(i, k - 1)/denom
                  hc(i, k) = hc(i, k) + x_add*up_massentr(i, k - 1)/denom
               end if
               !assuming zuo=zu,up_massdetro=up_massdetr, ...
               !(zuo(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))
               uc(i, k) = (uc(i, k - 1)*zu(i, k - 1) - .5*up_massdetru(i, k - 1)*uc(i, k - 1) + &
                           up_massentru(i, k - 1)*us(i, k - 1) &
                           - p_pgcon*.5*(zu(i, k) + zu(i, k - 1))*(u_cup(i, k) - u_cup(i, k - 1)))/denom_u

               vc(i, k) = (vc(i, k - 1)*zu(i, k - 1) - .5*up_massdetru(i, k - 1)*vc(i, k - 1) + &
                           up_massentru(i, k - 1)*vs(i, k - 1) &
                           - p_pgcon*.5*(zu(i, k) + zu(i, k - 1))*(v_cup(i, k) - v_cup(i, k - 1)))/denom_u
            else
               hc(i, k) = hc(i, k - 1)
               hco(i, k) = hco(i, k - 1)
               uc(i, k) = uc(i, k - 1)
               vc(i, k) = vc(i, k - 1)
            end if
            !---meltglac-------------------------------------------------
            !- includes glaciation effects on HC,HCO
            !                    ------ ice content --------
            !print*,"H=",hc (i,k),(1.-p_liq_ice(i,k))*qrco(i,k)*xlf,hc (i,k)+(1.-p_liq_ice(i,k))*qrco(i,k)*xlf
            hc(i, k) = hc(i, k) + (1.-p_liq_ice(i, k))*qrco(i, k)*c_xlf
            hco(i, k) = hco(i, k) + (1.-p_liq_ice(i, k))*qrco(i, k)*c_xlf
         end do

         do k = ktop(i) + 2, ktf
            hc(i, k) = hes_cup(i, k)!= he_cup(i,k)
            uc(i, k) = u_cup(i, k)
            vc(i, k) = v_cup(i, k)
            hco(i, k) = heso_cup(i, k)!=heo_cup(i,k)
            zu(i, k) = 0.
            zuo(i, k) = 0.
         end do
      end do

      !--- Get buoyancy of updrafts
      call getBuoyancy(itf, ktf, its, ite, kts, kte, ierr, klcl, kbcon, ktop, hc, he_cup, hes_cup, dby, z_cup)
      call getBuoyancy(itf, ktf, its, ite, kts, kte, ierr, klcl, kbcon, ktop, hco, heo_cup, heso_cup, dbyo, zo_cup)

      if (CONVECTION_TRACER == 1 .and. SGS_W_TIMESCALE == 1 .and. trim(cumulus) == 'deep') then
         !--- compute vertical velocity
         !
         !call cup_up_vvel(vvel2d,vvel1d,zws,entr_rate_2d,cd,zo,zo_cup,zuo,dbyo,GAMMAo_CUP,tn_cup &
         !                ,tempco,qco,qrco,qo,start_level,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte&
         !                ,wlpool,wlpool_bcon,2)
         wlpool_bcon(:) = wlpool(:)
         !
         !--- trigger function based on KE > CIN
         if (ADD_COLDPOOL_CLOS == 2) then
            call cupUpAa0(cin1, zo_cup, zuo, dbyo, gammao_cup, tn_cup, k22, klcl, kbcon, ktop, ierr, itf, ktf, its &
                          , ite, kts, kte, 'CIN')
            do i = its, itf
               if (ierr(i) /= 0) cycle
               ke_mx = 0.5*max(wlpool_bcon(i)**2, zws(i)**2) + 1.e-6
               if (ke_mx < abs(min(cin1(i), 0.))) ierr(i) = 500
            end do
         end if
      end if

      if (.not. first_guess_w) then
         !--- calculate in-cloud/updraft air temperature for vertical velocity
         do i = its, itf
            if (ierr(i) == 0) then
               do k = kts, ktf
                  tempco(i, k) = (1./real(c_cp))*(hco(i, k) - c_grav*zo_cup(i, k) - real(c_alvl)*qco(i, k))
               end do
               tempco(i, kte) = tn_cup(i, kte)
            else
               tempco(i, :) = tn_cup(i, :)
            end if
         end do

         !--- vertical velocity
         call cupUpVVel(vvel2d, vvel1d, zws, entr_rate_2d, cd, zo, zo_cup, zuo, dbyo, gammao_cup, tn_cup &
                , tempco, qco, qrco, qo, start_level, klcl, kbcon, ktop, ierr, itf, ktf, its, ite, kts, kte, wlpool, wlpool_bcon, 1)
      end if

      !---- new rain
      !
      !--- calculate rain mixing ratio in updrafts
      !
      !       call cup_up_rain(cumulus,klcl,kbcon,ktop,k22,ierr,xland         &
      !                       ,zo_cup,qco,qrco,pwo,pwavo,po,p_cup,t_cup,tempco&
      !                       ,zuo,up_massentr,up_massdetr,vvel2d,rho         &
      !                       ,qrr                                            &
      !                       ,itf,ktf,its,ite, kts,kte)
      !--- DOWNDRAFT section
      !
      do i = its, itf
         kzdown(i) = 0
         if (ierr(i) .eq. 0) then
            zktop = (zo_cup(i, ktop(i)) - z1(i))*.6
            zktop = min(zktop + z1(i), zcutdown + z1(i))
            do k = kts, ktf
               if (zo_cup(i, k) .gt. zktop) then
                  kzdown(i) = k
                  go to 37
               end if
            end do
         end if
37       continue
      end do

      !--- downdraft originating level - jmin
      call cupMinimi(heso_cup, k22, kzdown, jmin, ierr, itf, ktf, its, ite, kts, kte)

      call getJmin(cumulus, itf, ktf, its, ite, kts, kte, ierr, kdet, ktop, kbcon, jmin, ierrc, beta, depth_min, heso_cup, zo_cup &
                 , melting_layer)

      !--- this calls routine to get downdrafts normalized mass flux
      do i = its, itf
         zd(i, :) = 0.
         if (ierr(i) /= 0) cycle
         call getZuZdPdf(trim(cumulus), "DOWN", ierr(i), kdet(i), jmin(i), zdo(i, :), kts, kte, ktf, kpbl(i), k22(i), kbcon(i) &
                       , klcl(i), po_cup(i, kts:kte), psur(i), xland(i), random(i))
      end do

      !---  calls routine to get lateral mass fluxes associated with downdrafts
      call getLateralMassFluxDown(trim(cumulus), itf, ktf, its, ite, kts, kte, ierr, jmin, zo_cup, zdo, xzd, zd, cdd &
                                , mentrd_rate_2d , dd_massentro, dd_massdetro, dd_massentr, dd_massdetr, cumulus, mentrd_rate &
                                , dd_massentru, dd_massdetru, lambau_dn)

      !---  calls routine to get wet bulb temperature and moisture at jmin
      if (USE_WETBULB == 1 .and. trim(cumulus) /= 'shallow') then
         do i = its, itf
            if (ierr(i) /= 0) cycle
            k = jmin(i)
            call getWetbulb(jmin(i), qo_cup(i, k), t_cup(i, k), po_cup(i, k), q_wetbulb(i), t_wetbulb(i))
            !print*,"wb       =",jmin,qo_cup(i,k),t_cup(i,k),q_wetbulb(i),t_wetbulb(i)
            !print*,"evap/cool=",q_wetbulb(i)-qo_cup(i,k),t_wetbulb(i)-t_cup(i,k)
         end do
      end if

      !--- downdraft moist static energy + moisture budget
      do i = its, itf
         hcdo(i, :) = heso_cup(i, :)
         ucd(i, :) = u_cup(i, :)
         vcd(i, :) = v_cup(i, :)
         dbydo(i, :) = 0.
      end do

      do i = its, itf
         bud(i) = 0.
         if (ierr(i) /= 0 .or. trim(cumulus) == 'shallow') cycle
         i_wb = 0
         !--for future test
         if (USE_WETBULB == 1) then
            !--option 1
            !hcdo(i,jmin(i))=cp*t_wetbulb(i)+xlv*q_wetbulb(i)+zo_cup(i,jmin(i))*g
            !--option 2
            hcdo(i, jmin(i)) = 0.5*(real(c_cp)*t_wetbulb(i) + real(c_alvl)*q_wetbulb(i) + zo_cup(i, jmin(i))*c_grav + hc(i, jmin(i)))
            i_wb = 1
         end if

         dbydo(i, jmin(i)) = hcdo(i, jmin(i)) - heso_cup(i, jmin(i))
         bud(i) = dbydo(i, jmin(i))*(zo_cup(i, jmin(i) + 1) - zo_cup(i, jmin(i)))

         do ki = jmin(i) - i_wb, kts, -1!do ki=jmin(i)-1,1,-1
            denom = zdo(i, ki + 1) - 0.5*dd_massdetro(i, ki) + dd_massentro(i, ki)
            denom_u = zdo(i, ki + 1) - 0.5*dd_massdetru(i, ki) + dd_massentru(i, ki)
            !-tmp fix for denominator being zero
            if (denom > 0.0 .and. denom_u > 0.0) then
               dzo = zo_cup(i, ki + 1) - zo_cup(i, ki)

               ucd(i, ki) = (ucd(i, ki + 1)*zdo(i, ki + 1) - .5*dd_massdetru(i, ki)*ucd(i, ki + 1) + dd_massentru(i, ki)*us(i, ki) &
                          - p_pgcon*zdo(i, ki + 1)*(us(i, ki + 1) - us(i, ki)))/denom_u
               vcd(i, ki) = (vcd(i, ki + 1)*zdo(i, ki + 1) - .5*dd_massdetru(i, ki)*vcd(i, ki + 1) + dd_massentru(i, ki)*vs(i, ki) &
                          - p_pgcon*zdo(i, ki + 1)*(vs(i, ki + 1) - vs(i, ki)))/denom_u

               hcdo(i, ki) = (hcdo(i, ki + 1)*zdo(i, ki + 1) - .5*dd_massdetro(i, ki)*hcdo(i, ki + 1) + dd_massentro(i, ki) &
                           * heo(i, ki))/denom

               dbydo(i, ki) = hcdo(i, ki) - heso_cup(i, ki)
               !if(i.eq.ipr)write(0,*)'ki,bud = ',ki,bud(i),hcdo(i,ki)
               bud(i) = bud(i) + dbydo(i, ki)*dzo
            else
               ucd(i, ki) = ucd(i, ki + 1)
               vcd(i, ki) = vcd(i, ki + 1)
               hcdo(i, ki) = hcdo(i, ki + 1)
            end if
         end do
         if (bud(i) .gt. 0) then
            ierr(i) = 7
            ierrc(i) = 'downdraft is not negatively buoyant '
         end if
      end do

      !--- calculate moisture properties of downdraft
      call cupDdMoisture(cumulus, ierrc, zdo, hcdo, heso_cup, qcdo, qeso_cup, pwdo, qo_cup, zo_cup, dd_massentro, dd_massdetro &
                     ,   jmin, ierr, gammao_cup, pwevo, bu, qrcdo, qo, heo, tn_cup, 1, t_wetbulb, q_wetbulb, qco, pwavo, itf, ktf &
                     ,   its, ite, kts, kte)
                           !--test    pwevo,bu,qrcdo,qo,heo,t_cup,1,t_wetbulb,q_wetbulb,qco,pwavo,       &

      !--- calculate workfunctions for updrafts
      call cupUpAa0(aa0, z_cup, zu, dby, GAMMA_CUP, t_cup, k22, klcl, kbcon, ktop, ierr, itf, ktf, its, ite, kts, kte)
      call cupUpAa0(aa1, zo_cup, zuo, dbyo, gammao_cup, tn_cup, k22, klcl, kbcon, ktop, ierr, itf, ktf, its, ite, kts, kte)

      do i = its, itf
         if (ierr(i) /= 0) cycle
         if (aa1(i) .eq. 0.) then
            ierr(i) = 17
            ierrc(i) = "cloud work function zero"
         end if
      end do

      !--- Implements Becker et al (2021) closure, part 1
      if ((DICYCLE == 2 .or. DICYCLE == 3) .and. trim(cumulus) == 'deep') then
         do ki = 1, 2
            if (DICYCLE == 2 .and. ki == 2) cycle
            if (ki == 1) then
               !-- get the cloud work function for updrafts associated only with RAD + PBL
               tn_x = t + tn - tn_adv
               qo_x = q + qo - qo_adv

               !-- to check => aa1_radpbl=aa1
          !! tn_x = tn
          !! qo_x = qo
            end if
            if (ki == 2) then
               !-- get the cloud work function for updrafts associated only with Qv-advection
               !tn_x = tn_adv  ! orig
               tn_x = t        ! v2
               qo_x = qo_adv
            end if
            ierr_dummy = ierr

            call cupEnv(zo, qeso_x, heo_x, heso_x, tn_x, qo_x, po, z1, psur, ierr_dummy, -1, itf, ktf, its, ite, kts, kte)
            call cupEnvCLev(tn_x, qeso_x, qo_x, heo_x, heso_x, zo, po, qeso_cup_x, qo_cup_x, heo_cup_x, us, vs &
                              , u_cup_x, v_cup_x, heso_cup_x, zo_cup_x, po_cup_x, gammao_CUP_x, tn_cup_x, psur, tsur &
                              , ierr_dummy, z1, itf, ktf, its, ite, kts, kte)

            !--- get MSE
            do i = its, itf
               if (ierr_dummy(i) /= 0) cycle
               x_add = (real(c_alvl)*zqexec(i) + real(c_cp)*ztexec(i)) + x_add_buoy(i)
               call getCloudBc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),heo_cup_x(i,kts:kte),hkbo_x(i),k22(i),x_add &
                               , tpert(i,kts:kte))
               hco_x(i, kts:start_level(i)) = hkbo_x(i)

               do k = start_level(i) + 1, ktop(i) + 1  ! mass cons option
                  denom = (zu(i, k - 1) - .5*up_massdetr(i, k - 1) + up_massentr(i, k - 1))
                  if (denom > 0.0) then

                     hco_x(i, k) = (hco_x(i, k - 1)*zuo(i, k - 1) - .5*up_massdetro(i, k - 1)*hco_x(i, k - 1) &
                                 + up_massentro(i, k - 1)*heo_x(i, k - 1))/denom
                     if (k == start_level(i) + 1) then
                        x_add = (real(c_alvl)*zqexec(i) + real(c_cp)*ztexec(i)) + x_add_buoy(i)
                        hco_x(i, k) = hco_x(i, k) + x_add*up_massentro(i, k - 1)/denom
                     end if
                  else
                     hco_x(i, k) = hco_x(i, k - 1)
                  end if
                  !- includes glaciation effects on HCO_X
                  hco_x(i, k) = hco_x(i, k) + (1.-p_liq_ice(i, k))*qrco(i, k)*c_xlf
               end do
               hco_x(i, ktop(i) + 2:ktf) = heso_cup_x(i, ktop(i) + 2:ktf)
            end do

            call getBuoyancy(itf, ktf, its, ite, kts, kte, ierr_dummy, klcl, kbcon, ktop, hco_x, heo_cup_x, heso_cup_x, dbyo_x &
                           , zo_cup_x)

            if (ki == 1) then  ! RAD+PBL only
               call cupUpAa0(aa1_radpbl, zo_cup_x, zuo, dbyo_x, gammao_CUP_x, tn_cup_x, k22, klcl, kbcon, ktop, ierr_dummy, itf &
                           , ktf, its, ite, kts, kte)
               !-- get AA1_ADV
               !aa1_adv = aa1 + aa0 - aa1_radpbl
            end if

            if (ki == 2) & ! ADV of Qv only
               call cupUpAa0(aa1_adv, zo_cup_x, zuo, dbyo_x, gammao_CUP_x, tn_cup_x, k22, klcl, kbcon, ktop, ierr_dummy, itf, ktf &
                           , its, ite, kts, kte)
            !Observe that :
            !aa1 ~ aa0 + (aa1_radpbl-aa0) + (aa1_adv-aa0)
         end do ! ki
      end if
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
      do i = its, itf
         if (ierr(i) == 0) then
            do k = kts, ktf
               tempcdo(i, k) = (1./real(c_cp))*(hcdo(i, k) - c_grav*zo_cup(i, k) - real(c_alvl)*qcdo(i, k))
            end do
         else
            tempcdo(i, :) = tn_cup(i, :)
         end if
      end do

      !--- diurnal cycle section
      !--- Bechtold et al 2008 time-scale of cape removal
      if (trim(cumulus) == 'deep') then
         tau_ecmwf(:) = TAU_DEEP; wmean(:) = 3. !  mean vertical velocity m/s
      else
         tau_ecmwf(:) = TAU_MID; wmean(:) = 3.
      end if
      !--- we shall let all scale dependence on the sig parameter
      !tau_ecmwf(:)= tau_ecmwf(:) * (1. + 1.66 * (dx(:)/(125*1000.)))! dx must be in meters
      if (SGS_W_TIMESCALE == 1 .and. trim(cumulus) == 'deep') then
         do i = its, itf
            if (ierr(i) /= 0) cycle
            !- mean vertical velocity based on integration of vertical veloc equation
            wmean(i) = min(max(vvel1d(i), 3.), 20.)

            !- time-scale cape removal from Bechtold et al. 2008
            tau_ecmwf(i) = (zo_cup(i, ktop(i)) - zo_cup(i, kbcon(i)))/wmean(i)
            !tau_ecmwf(i)= min(10800., max(720.,tau_ecmwf(i)))
            tau_ecmwf(i) = min(10800., max(1000., tau_ecmwf(i)))
         end do
!====
!  do i=its,itf
!      if(ierr(i) /= 0) cycle
!      print*,'tauec',wlpool_bcon(i),vvel1d(i), wmean(i) ,tau_ecmwf(i)
!      call flush(6)
!  enddo
!====
      end if

      !--- Implements the Bechtold et al (2014) and Becker et al (2021) closures
      do i = its, itf
         if (ierr(i) /= 0) cycle
         !- over water
         !   umean= 2.0+sqrt(0.5*(US(i,1)**2+VS(i,1)**2+US(i,kbcon(i))**2+VS(i,kbcon(i))**2))
         !   tau_bl(i) = (zo_cup(i,kbcon(i))- z1(i)) /umean
         !- over land
         !   tau_bl(i) = tau_ecmwf(i)
         !-----------
         umean = 2.0 + sqrt(0.5*(us(i, 1)**2 + vs(i, 1)**2 + us(i, kbcon(i))**2 + vs(i, kbcon(i))**2))
         !--                    - over land -            -          over ocean       s   -
         tau_bl(i) = (1.-xland(i))*tau_ecmwf(i) + xland(i)*(zo_cup(i, kbcon(i)) - z1(i))/umean
      end do

      if (DICYCLE <= 3 .and. trim(cumulus) == 'deep') then

         !-- calculate "pcape" or equivalent cloud work function from the BL forcing only
         iversion = 0
         call cupUpAa1Bl(iversion, aa1_bl, aa1_fa, aa1, t, tn, q, qo, dtime, po_cup, zo_cup, zuo, dbyo, gammao_cup, tn_cup &
                        ,rho, klcl, kpbl, kbcon, ktop, ierr, itf, ktf, its, ite, kts, kte, xland, ztexec, xlons, xlats &
                        , h_sfc_flux, le_sfc_flux, tau_bl, tau_ecmwf, t_star, cumulus, tn_bl, qo_bl)

         do i = its, itf
            if (ierr(i) /= 0) cycle
            aa1_bl(i) = (aa1_bl(i)/t_star)*tau_bl(i) ! units J/kg
            !aa1_bl(i) = (aa1_bl(i)/T_star) * tau_bl(i) - cin1(i)
            aa1_bl(i) = min(2000., abs(aa1_bl(i)))*sign(1., aa1_bl(i))
         end do

         !--- Adds Becker et al (2021) closure, part 2
         if (DICYCLE == 2) then
            call getQadv(cumulus, itf, ktf, its, ite, kts, kte, ierr, dtime, q, qo, qo_adv, po, po_cup, qeso, q_adv, col_sat_adv &
                       , alpha_adv, tau_bl, zo_cup, kbcon, ktop)
         end if

         if (DICYCLE == 3) then
            do i = its, itf
               if (ierr(i) /= 0) cycle
               aa1_adv(i) = (aa1_adv(i) - aa0(i))*tau_bl(i)/dtime
            end do
         end if

         !--- Implements the Zhang(2002) closure
      elseif (DICYCLE == 4 .and. trim(cumulus) == 'deep') then
         !- T and Q profiles modified only by RAD+ADV tendencies
         do i = its, itf
            if (ierr(i) /= 0) cycle
            tn_x(i, kts:ktf) = tn(i, kts:ktf) - tn_bl(i, kts:ktf) + t(i, kts:ktf)
            qo_x(i, kts:ktf) = qo(i, kts:ktf) - qo_bl(i, kts:ktf) + q(i, kts:ktf)
         end do
         !--- calculate moist static energy, heights, qes, ... only by free troposphere tendencies
         call cupEnv(zo, qeso_x, heo_x, heso_x, tn_x, qo_x, po, z1, psur, ierr, -1, itf, ktf, its, ite, kts, kte)
         !--- environmental values on cloud levels only by FT tendencies
         call cupEnvCLev(tn_x, qeso_x, qo_x, heo_x, heso_x, zo, po, qeso_cup_x, qo_cup_x, heo_cup_x, us, vs, u_cup, v_cup, &
                           heso_cup_x, zo_cup, po_cup, gammao_CUP_x, tn_cup_x, psur, tsur, ierr, z1, itf, ktf, its, ite, kts, kte)
         !--- this is (DT_ve/Dt)_adv+rad
         do i = its, itf
            if (ierr(i) /= 0) cycle
            aa3(i) = 0.
            do k = max(kbcon(i), kts + 1), ktop(i)
               dp = -(log(100.*po(i, k)) - log(100.*po(i, k - 1))) !no units
               aa3(i) = aa3(i) - (tn_cup_x(i, k)*(1.+0.608*qo_cup_x(i, k)) - t_cup(i, k)*(1.+0.608*q_cup(i, k)))*dp/dtime !units = K/s
               !print*,"tve=",k,aa3(i),tn_cup_x(i,k)*(1.+0.608*qo_cup_x(i,k)),&
               !               t_cup (i,k)*(1.+0.608*q_cup   (i,k)),dp
            end do
         end do
         do i = its, itf
            if (ierr(i) /= 0) cycle
            !- this is (DCAPE_env/Dt)_adv+rad
            !aa1_bl(i) = -aa3(i)
            !- Zhang threshold:  65 J/kg/hour => 65/(Rd *3600)= 63 10^-6 K/s
            aa1_bl(i) = aa3(i) - (63.e-6)!*1.5
            !print*,"dcape_env=",aa3(i),aa1_bl(i)
            if (xland(i) > 0.90) aa1_bl(i) = 1.4*aa1_bl(i) !- over water
         end do
         !--- this is (DT_ve/Dt)_cu
         do i = its, itf
            dtdt(i, :) = 0.
            dqdt(i, :) = 0.
            if (ierr(i) /= 0) cycle
            do k = max(kbcon(i), kts + 1), ktop(i)
               dp = 100.*(po_cup(i, k + 1) - po_cup(i, k))
               rz_env = 0.5*(zuo(i, k + 1) + zuo(i, k) - (zdo(i, k + 1) + zdo(i, k))*edto(i))
               s2 = real(c_cp)*tn_cup_x(i, k + 1) + c_grav*zo_cup(i, k + 1)
               s1 = real(c_cp)*tn_cup_x(i, k) + c_grav*zo_cup(i, k)
               q2 = qo_cup_x(i, k + 1)
               q1 = qo_cup_x(i, k)

               dqdt(i, k) = -rz_env*(q2 - q1)*c_grav/dp
               dtdt(i, k) = -(1./real(c_cp))*rz_env*(s2 - s1)*c_grav/dp

               dqdt(i, k) = dqdt(i, k) + (up_massdetro(i, k)*0.5*(qco(i, k + 1) + qco(i, k) - (q2 + q1)) + edto(i) &
                          * dd_massdetro(i, k)*0.5*(qcdo(i, k + 1) + qcdo(i, k) - (q2 + q1)))*c_grav/dp

               dtdt(i, k) = dtdt(i, k) + (up_massdetro(i, k)*0.5*(tempco(i, k + 1) + tempco(i, k) - (tn_cup_x(i, k + 1) &
                          + tn_cup_x(i, k))) + edto(i)*dd_massdetro(i, k)*0.5*(tempcdo(i, k + 1) + tempcdo(i, k) &
                          - (tn_cup_x(i, k + 1) + tn_cup_x(i, k))))*c_grav/dp
               !print*,"dtdt=",k, dtdt(i,k),zuo(i,k+1),zdo(i,k+1),dqdt(i,k)
            end do
            xk_x(i) = 0.
            do k = max(kbcon(i), kts + 1), ktop(i)
               dp = -(log(100.*po_cup(i, k + 1)) - log(100.*po_cup(i, k)))      ! no units here
               xk_x(i) = xk_x(i) + ((1.+0.608*qo_x(i, k))*dtdt(i, k) + 0.608*tn_x(i, k)*dqdt(i, k))*dp !  units=K m/Pa s2
               !=> aa3/xk_x will have units of kg/m2/s for the mass flux at cloud base.
               !print*,"xk_x=",k, xk_x(i),dtdt(i,k),dqdt(i,k)
            end do
         end do
      end if

      !--- Trigger function based on Xie et al (2019)
      if (ADV_TRIGGER == 1 .and. trim(cumulus) == 'deep') then
         daa_adv_dt = 0.
         do step = 1, 2
            !--- calculate moist static energy, heights, qes, ... only by ADV tendencies
            if (step == 1) then
               tn_x = t
               qo_x = q
            else
               tn_x = tn_adv
               qo_x = qo_adv
            end if
            call cupEnv(zo, qeso_x, heo_x, heso_x, tn_x, qo_x, po, z1, psur, ierr, -1, itf, ktf, its, ite, kts, kte)
            call cupEnvCLev(tn_x, qeso_x, qo_x, heo_x, heso_x, zo, po, qeso_cup_x, qo_cup_x, heo_cup_x, us, vs &
                              , u_cup_x, v_cup_x, heso_cup_x, zo_cup_x, po_cup_x, gammao_CUP_x, tn_cup_x, psur, tsur &
                              , ierr, z1, itf, ktf, its, ite, kts, kte)

            !--- get MSE
            do i = its, itf
               if (ierr(i) /= 0) cycle
               call getCloudBc(cumulus, kts, kte, ktf, xland(i), po(i, kts:kte), heo_cup_x(i, kts:kte), hkbo_x(i), k22(i))
               hco_x(i, kts:start_level(i)) = hkbo_x(i)

               do k = start_level(i) + 1, ktop(i) + 1

                  denom = (zuo(i, k - 1) - .5*up_massdetro(i, k - 1) + up_massentro(i, k - 1))
                  if (denom > 0.0) then
                     hco_x(i, k) = (hco_x(i, k - 1)*zuo(i, k - 1) - .5*up_massdetro(i, k - 1)*hco_x(i, k - 1) &
                                 + up_massentro(i, k - 1)*heo_x(i, k - 1))/denom
                  else
                     hco_x(i, k) = hco_x(i, k - 1)
                  end if
               end do
               hco_x(i, ktop(i) + 2:ktf) = heso_cup_x(i, ktop(i) + 2:ktf)
            end do
            call getBuoyancy(itf, ktf, its, ite, kts, kte, ierr, klcl, kbcon, ktop, hco_x, heo_cup_x, heso_cup_x, dbyo_x, zo_cup_x)
            !--- get cloud work function
            aa_tmp = 0.
            call cupUpAa0(aa_tmp, zo_cup_x, zuo, dbyo_x, gammao_CUP_x, tn_cup_x, k22, klcl, kbcon, ktop, ierr, itf, ktf, its, ite &
                        , kts, kte)

            if (step == 1) aa_ini = aa_tmp ! cloud work function initial
            if (step == 2) aa_adv = aa_tmp ! cloud work function modified by advection tendencies
         end do
         !
         do i = its, itf
            if (ierr(i) /= 0) cycle

            daa_adv_dt(i) = (aa_adv(i) - aa_ini(i))/dtime

            !print*,"daa_adv_dt J. kg-1 hr-1=",daa_adv_dt(i)*3600.
            ! call flush(6)

            if (daa_adv_dt(i) > DCAPE_THRESHOLD/3600. .and. aa_ini(i) > 0.) cycle !
            ierr(i) = 90
            ierrc(i) = "dcape trigger not satisfied"

         end do
         !--- only for output
         aaa0_(:) = daa_adv_dt(:)*3600. ! J/kg/hour

      end if

      !--- determine downdraft strength in terms of windshear
      call cupDdEdt(cumulus, ierr, us, vs, zo, ktop, kbcon, edt, po, pwavo,pwo, ccn, pwevo, edtmax, edtmin, p_maxens2 &
                  , edtc, psum, psumh, rho, p_aeroevap, itf, ktf, ipr, jpr, its, ite, kts, kte, vshear)

      do iedt = 1, p_maxens2
         do i = its, itf
            if (ierr(i) .eq. 0) then
               edto(i) = sigd(i)*edtc(i, iedt)
               edt(i) = edto(i)
            end if
         end do

         !--- get the environmental mass flux
         do i = its, itf
            zenv(i, :) = 0.0
            if (ierr(i) /= 0) cycle
            zenv(i, :) = zuo(i, :) - edto(i)*zdo(i, :)
         end do

         !--- check mass conservation
         do i = its, itf
            if (ierr(i) /= 0) cycle
            do k = kts, ktop(i)
               ! these three are only used at or near mass detrainment and/or entrainment levels
               entupk = 0.
               detupk = 0.
               entdoj = 0.
               ! detrainment and entrainment for downdrafts
               detdo = edto(i)*dd_massdetro(i, k)
               entdo = edto(i)*dd_massentro(i, k)
               ! entrainment/detrainment for updraft
               entup = up_massentro(i, k)
               detup = up_massdetro(i, k)
               ! subsidence by downdrafts only
               subin = -zdo(i, k + 1)*edto(i)
               subdown = -zdo(i, k)*edto(i)
               if (k .eq. ktop(i)) then
                  detupk = zuo(i, ktop(i))
                  subin = 0.
                  subdown = 0.
                  detdo = 0.
                  entdo = 0.
                  entup = 0.
                  detup = 0.
               end if
               totmas = subin - subdown + detup - entup - entdo + &
                        detdo - entupk - entdoj + detupk + zuo(i, k + 1) - zuo(i, k)
               if (abs(totmas) .gt. 1.e-6) then
                  write (6, *) '**mass cons: k,ktop,zo(ktop),totmas,subin,subdown,detup,entup,detdo,entdo,entupk,detupk'
                  write (6, 123) 'mass*1.e+6', k, ktop(i), zo(i, ktop(i)), totmas*1.e+6, subin*1.e+6, subdown*1.e+6, detup*1.e+6 &
                                , entup*1.e+6, detdo*1.e+6, entdo*1.e+6, entupk*1.e+6, detupk*1.e+6
123               format(1x, A11, 2i5, 10e12.5)
                  ! call error_fatal ( 'totmas .gt.1.e-6' )
               end if
            end do   ! k
         end do

         !--- change per unit mass that a model cloud would modify the environment
         !--- 1. in bottom layer
         dellu = 0.
         dellv = 0.
         dellah = 0.
         dellat = 0.
         dellaq = 0.
         dellaqc = 0.
         dellabuoy = 0.
         subten_h = 0.
         subten_q = 0.
         subten_t = 0.

         if (VERT_DISCR == 0) then
            do i = its, itf
               if (ierr(i) /= 0) cycle
               do k = kts, ktop(i)
                  dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))
                  dellu(i, k) = -(zuo(i, k + 1)*(uc(i, k + 1) - u_cup(i, k + 1)) - zuo(i, k)*(uc(i, k) - u_cup(i, k)))*c_grav/dp &
                              + (zdo(i, k + 1)*(ucd(i, k + 1) - u_cup(i, k + 1)) - zdo(i, k)*(ucd(i, k) - u_cup(i, k)))*c_grav/dp &
                              * edto(i)

                  dellv(i, k) = -(zuo(i, k + 1)*(vc(i, k + 1) - v_cup(i, k + 1)) - zuo(i, k)*(vc(i, k) - v_cup(i, k)))*c_grav/dp &
                              + (zdo(i, k + 1)*(vcd(i, k + 1) - v_cup(i, k + 1)) - zdo(i, k)*(vcd(i, k) - v_cup(i, k)))*c_grav/dp &
                              * edto(i)
               end do   ! k
            end do

            do i = its, itf
               trash = 0.0
               trash2 = 0.0
               if (ierr(i) .eq. 0) then
                  do k = kts, ktop(i)
                     dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))
                     dellah(i, k) = -(zuo(i, k + 1)*(hco(i, k + 1) - heo_cup(i, k + 1)) - zuo(i, k)*(hco(i, k) - heo_cup(i, k))) &
                                  * c_grav/dp + (zdo(i, k + 1)*(hcdo(i, k + 1) - heo_cup(i, k + 1)) - zdo(i, k)*(hcdo(i, k) &
                                  - heo_cup(i, k)))*c_grav/dp*edto(i)

                     !---meltglac-------------------------------------------------
                     dellah(i, k) = dellah(i, k) + c_xlf*((1.-p_liq_ice(i, k))*0.5*(qrco(i, k + 1) + qrco(i, k)) - melting(i, k)) &
                                  * c_grav/dp

                     !-- for output only
                     subten_h(i, k) = -(zuo(i, k + 1)*(-heo_cup(i, k + 1)) - zuo(i, k)*(-heo_cup(i, k)))*c_grav/dp &
                                    + (zdo(i, k + 1)*(-heo_cup(i, k + 1)) - zdo(i, k)*(-heo_cup(i, k)))*c_grav/dp*edto(i)

                     !- check H conservation
                     trash2 = trash2 + (dellah(i, k))*dp/c_grav

                     !-- take out cloud liquid/ice water for detrainment
                     detup = up_massdetro(i, k)
                     if (trim(cumulus) == 'mid' .or. trim(cumulus) == 'shallow') then
                        dellaqc(i, k) = detup*0.5*(qrco(i, k + 1) + qrco(i, k))*c_grav/dp
                     elseif (trim(cumulus) == 'deep') then
                        if (.not. use_c1d) then
                           dellaqc(i, k) = detup*0.5*(qrco(i, k + 1) + qrco(i, k))*c_grav/dp
                        elseif (C1 > 0.0) then
                           if (k == ktop(i)) then
                              dellaqc(i, k) = detup*0.5*(qrco(i, k + 1) + qrco(i, k))*c_grav/dp
                           else
                              dz = zo_cup(i, k + 1) - zo_cup(i, k)
                              dellaqc(i, k) = zuo(i, k)*c1d(i, k)*qrco(i, k)*dz/dp*c_grav
                           end if
                        else
                           if (k == ktop(i)) then
                              dellaqc(i, k) = detup*0.5*(qrco(i, k + 1) + qrco(i, k))*c_grav/dp
                           else
                              dz = zo_cup(i, k + 1) - zo_cup(i, k)
                              dellaqc(i, k) = (zuo(i, k)*c1d(i, k)*qrco(i, k)*dz/dp*c_grav + detup*0.5*(qrco(i, k + 1) &
                                            + qrco(i, k))*c_grav/dp)*0.5
                           end if
                        end if
                     end if

                     g_rain = 0.5*(pwo(i, k) + pwo(i, k + 1))*c_grav/dp
                     e_dn = -0.5*(pwdo(i, k) + pwdo(i, k + 1))*c_grav/dp*edto(i) ! pwdo < 0 and E_dn must > 0

                     !-- condensation source term = detrained + flux divergence of
                     !-- cloud liquid/ice water (qrco) + converted to rain
                     c_up = dellaqc(i, k) + (zuo(i, k + 1)*qrco(i, k + 1) - zuo(i, k)*qrco(i, k))*c_grav/dp + g_rain

                     !-- water vapor budget
                     !-- = flux divergence z*(Q_c - Q_env)_up_and_down &
                     !--   - condensation term + evaporation
                     dellaq(i, k) = -(zuo(i, k + 1)*(qco(i, k + 1) - qo_cup(i, k + 1)) - zuo(i, k)*(qco(i, k) - qo_cup(i, k))) &
                                  * c_grav/dp + (zdo(i, k + 1)*(qcdo(i, k + 1) - qo_cup(i, k + 1)) - zdo(i, k)*(qcdo(i, k) &
                                  - qo_cup(i, k)))*c_grav/dp*edto(i) - c_up + e_dn

                     !-- for output only
                     subten_q(i, k) = -(zuo(i, k + 1)*(-qo_cup(i, k + 1)) - zuo(i, k)*(-qo_cup(i, k)))*c_grav/dp &
                                    + (zdo(i, k + 1)*(-qo_cup(i, k + 1)) - zdo(i, k)*(-qo_cup(i, k)))*c_grav/dp*edto(i)

                     !- check water conservation liq+condensed (including rainfall)
                     trash = trash + (dellaq(i, k) + dellaqc(i, k) + g_rain - e_dn)*dp/c_grav

                     dellabuoy(i, k) = edto(i)*dd_massdetro(i, k)*0.5*(dbydo(i, k + 1) + dbydo(i, k))*c_grav/dp

                     !write(3,*)'=>H= ',k,real(trash2,4),real(dellah(i,k),4)
                     !write(4,*)'=>W= ',k,real(trash,4),real(dellaq(i,k),4)
                  end do   ! k
                  !--- test only with double precision:
                  !write(0,*)'=>H/W-FINAL= ',real(trash2,4),real(trash,4),k22(i),kbcon(i),ktop(i)
                  !if(abs(trash)>1.e-6 .or. abs(trash2) > 1.e-6) then
                  !    write(0,*)'=> not water mass or H cons for deep= ',i,trash,trash2
                  !    !stop 33
                  !endif
               end if
            end do
         elseif (VERT_DISCR == 1) then
            !---- convective transport of momentum
            if (ALP1 == 0.) then !-- fully time explicit
               do i = its, itf
                  if (ierr(i) /= 0) cycle
                  do k = kts, ktop(i)
                     dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))

                     dellu(i, k) = -(zuo(i, k + 1)*(uc(i, k + 1) - u_cup(i, k + 1)) - zuo(i, k)*(uc(i, k) - u_cup(i, k))) & 
                                 * c_grav/dp + (zdo(i, k + 1)*(ucd(i, k + 1) - u_cup(i, k + 1)) - zdo(i, k)*(ucd(i, k)  &
                                 - u_cup(i, k)))*c_grav/dp*edto(i)

                     dellv(i, k) = -(zuo(i, k + 1)*(vc(i, k + 1) - v_cup(i, k + 1)) - zuo(i, k)*(vc(i, k) - v_cup(i, k))) &
                                 * c_grav/dp + (zdo(i, k + 1)*(vcd(i, k + 1) - v_cup(i, k + 1)) - zdo(i, k)*(vcd(i, k) & 
                                 - v_cup(i, k)))*c_grav/dp*edto(i)
                  end do   ! k
               end do
            elseif (ALP1 > 0.) then              !-- time alp0*explict + ALP1*implicit + upstream
               alp0 = 1.-ALP1
               do i = its, itf
                  if (ierr(i) /= 0) cycle
                  do k = kts, ktop(i) + 1
                     fp(k) = 0.5*(zenv(i, k) + abs(zenv(i, k)))
                     fm(k) = 0.5*(zenv(i, k) - abs(zenv(i, k)))
                  end do

                  do k = kts, ktop(i)
                     dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))

                     beta1 = dtime*c_grav/dp
                     aa(k) = ALP1*beta1*fm(k)
                     bb(k) = 1.+ALP1*beta1*(fp(k) - fm(k + 1))
                     cc(k) = -ALP1*beta1*fp(k + 1)

                     ddu(k) = us(i, k) - (zuo(i, k + 1)*uc(i, k + 1) - zuo(i, k)*uc(i, k))*beta1 + (zdo(i, k + 1)*ucd(i, k + 1) &
                            - zdo(i, k)*ucd(i, k))*beta1*edto(i)

                     ddu(k) = ddu(k) + alp0*beta1*(-fm(k)*us(i, max(kts, k - 1)) + (fm(k + 1) - fp(k))*us(i, k) + fp(k + 1) &
                            * us(i, k + 1))

                     ddv(k) = vs(i, k) - (zuo(i, k + 1)*vc(i, k + 1) - zuo(i, k)*vc(i, k))*beta1 + (zdo(i, k + 1)*vcd(i, k + 1) &
                            - zdo(i, k)*vcd(i, k))*beta1*edto(i)

                     ddv(k) = ddv(k) + alp0*beta1*(-fm(k)*vs(i, max(kts, k - 1)) + (fm(k + 1) - fp(k))*vs(i, k)  + fp(k + 1) &
                            * vs(i, k + 1))

                  end do
                  call tridiag(ktop(i), aa(kts:ktop(i)), bb(kts:ktop(i)), cc(kts:ktop(i)), ddu(kts:ktop(i)))
                  dellu(i, kts:ktop(i)) = (ddu(kts:ktop(i)) - us(i, kts:ktop(i)))/dtime

                  call tridiag(ktop(i), aa(kts:ktop(i)), bb(kts:ktop(i)), cc(kts:ktop(i)), ddv(kts:ktop(i)))
                  dellv(i, kts:ktop(i)) = (ddv(kts:ktop(i)) - vs(i, kts:ktop(i)))/dtime
               end do
            end if

            !--- convective transport of MSE and Q/Qc
            !if(USE_FLUX_FORM == 1) then
            do i = its, itf
               if (ierr(i) /= 0) cycle

               !--- moist static energy : flux form + source/sink terms + time explicit
               !
               !   if(use_fct == 0 .or. adjustl(cumulus) == 'shallow') then
               if (USE_FCT == 0) then
                  do k = kts, ktop(i)
                     dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))
                     dellah(i, k) = -(zuo(i, k + 1)*(hco(i, k + 1) - heo_cup(i, k + 1)) - zuo(i, k)*(hco(i, k) - heo_cup(i, k))) &
                                  * c_grav/dp + (zdo(i, k + 1)*(hcdo(i, k + 1) - heo_cup(i, k + 1)) - zdo(i, k)*(hcdo(i, k) &
                                  - heo_cup(i, k)))*c_grav/dp*edto(i)

                     dellah(i, k) = dellah(i, k) + c_xlf*((1.-p_liq_ice(i, k))* 0.5*(qrco(i, k + 1) + qrco(i, k)) - melting(i, k)) &
                                  * c_grav/dp

                     !--- for output only
                     subten_h(i, k) = -(zuo(i, k + 1)*(-heo_cup(i, k + 1)) - zuo(i, k)*(-heo_cup(i, k)))*c_grav/dp &
                                    + (zdo(i, k + 1)*(-heo_cup(i, k + 1)) - zdo(i, k)*(-heo_cup(i, k)))*c_grav/dp*edto(i)
                  end do   ! k
               else

                  !-- FCT scheme for the subsidence transport: d(M_env*S_env)/dz
                  sub_tend(1, :) = 0. ! dummy array
                  trcflx_in(1, :) = 0. ! dummy array
                  massflx(i, :) = 0.
                  dtime_max = dtime

                  do k = kts, ktop(i)
                     dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))
                     trcflx_in(1, k) = -(zuo(i, k) - edto(i)*zdo(i, k))*heo_cup(i, k) !* xmb(i)
                     massflx(i, k) = -(zuo(i, k) - edto(i)*zdo(i, k))        !* xmb(i)
                     dtime_max = min(dtime_max, .5*dp)
                  end do
                  call fct1d3(ktop(i), kte, dtime_max, po_cup(i, :), heo(i, :), massflx(i, :), trcflx_in(1, :), sub_tend(1, :))

                  do k = kts, ktop(i)
                     dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))
                     dellah(i, k) = -(zuo(i, k + 1)*hco(i, k + 1) - zuo(i, k)*hco(i, k))*c_grav/dp + (zdo(i, k + 1)*hcdo(i, k + 1) &
                                  - zdo(i, k)*hcdo(i, k))*c_grav/dp*edto(i)

                     dellah(i, k) = dellah(i, k) + c_xlf*((1.-p_liq_ice(i, k))* 0.5*(qrco(i, k + 1) + qrco(i, k)) - melting(i, k)) &
                                  * c_grav/dp
                     !- update with subsidence term from the FCT scheme
                     dellah(i, k) = dellah(i, k) + sub_tend(1, k)
                     !--- for output only
                     subten_h(i, k) = sub_tend(1, k)
                  end do   ! k
               end if
            end do

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
            do i = its, itf
               if (ierr(i) /= 0) cycle
               do k = kts, ktop(i)
                  dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))

                  !-- take out cloud liquid/ice water for detrainment
                  detup = up_massdetro(i, k)
                  if (trim(cumulus) == 'mid' .or. trim(cumulus) == 'shallow') then

                     dellaqc(i, k) = detup*0.5*(qrco(i, k + 1) + qrco(i, k))*c_grav/dp

                  elseif (trim(cumulus) == 'deep') then
                     if (.not. use_c1d) then
                        dellaqc(i, k) = detup*0.5*(qrco(i, k + 1) + qrco(i, k))*c_grav/dp
                     elseif (C1 > 0.0) then
                        if (k == ktop(i)) then
                           dellaqc(i, k) = detup*0.5*(qrco(i, k + 1) + qrco(i, k))*c_grav/dp
                        else
                           dz = zo_cup(i, k + 1) - zo_cup(i, k)
                           dellaqc(i, k) = zuo(i, k)*c1d(i, k)*qrco(i, k)*dz/dp*c_grav
                        end if
                     else
                        if (k == ktop(i)) then
                           dellaqc(i, k) = detup*0.5*(qrco(i, k + 1) + qrco(i, k))*c_grav/dp
                        else
                           dz = zo_cup(i, k + 1) - zo_cup(i, k)
                           dellaqc(i, k) = (zuo(i, k)*c1d(i, k)*qrco(i, k)*dz/dp*c_grav + detup*0.5*(qrco(i, k + 1) + qrco(i, k)) &
                                         * c_grav/dp)*0.5
                        end if
                     end if
                  end if

                  g_rain = 0.5*(pwo(i, k) + pwo(i, k + 1))*c_grav/dp
                  e_dn = -0.5*(pwdo(i, k) + pwdo(i, k + 1))*c_grav/dp*edto(i) ! pwdo < 0 and E_dn must > 0

                  !-- condensation source term = detrained + flux divergence of
                  !-- cloud liquid/ice water (qrco) + converted to rain
                  c_up = dellaqc(i, k) + (zuo(i, k + 1)*qrco(i, k + 1) - zuo(i, k)*qrco(i, k))*c_grav/dp + g_rain

                  !-- water vapor budget
                  !-- = flux divergence z*(Q_c - Q_env)_up_and_down  - condensation term + evaporation
                  dellaq(i, k) = -(zuo(i, k + 1)*qco(i, k + 1) - zuo(i, k)*qco(i, k))*c_grav/dp  + (zdo(i, k + 1)*qcdo(i, k + 1) &
                               - zdo(i, k)*qcdo(i, k))*c_grav/dp*edto(i) - c_up + e_dn

                  !--- source of cold pools
                  dellabuoy(i, k) = edto(i)*dd_massdetro(i, k)*0.5*(dbydo(i, k + 1) + dbydo(i, k))*c_grav/dp
               end do
               !        if(use_fct == 0 .or. adjustl(cumulus) == 'shallow') then
               if (USE_FCT == 0) then
                  do k = kts, ktop(i)
                     dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))
                     sub_tend(1, k) = -(zuo(i, k + 1)*(-qo_cup(i, k + 1)) - zuo(i, k)*(-qo_cup(i, k)))*c_grav/dp &
                                    + (zdo(i, k + 1)*(-qo_cup(i, k + 1)) - zdo(i, k)*(-qo_cup(i, k)))*c_grav/dp*edto(i)
                  end do
               else
                  !-- FCT scheme for the subsidence transport: d(M_env*S_env)/dz
                  sub_tend(1, :) = 0. ! dummy array
                  trcflx_in(1, :) = 0. ! dummy array
                  massflx(i, :) = 0.
                  dtime_max = dtime
                  do k = kts, ktop(i)
                     dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))
                     trcflx_in(1, k) = -(zuo(i, k) - edto(i)*zdo(i, k))*qo_cup(i, k) !* xmb(i)
                     massflx(i, k) = -(zuo(i, k) - edto(i)*zdo(i, k))             !* xmb(i)
                     dtime_max = min(dtime_max, .5*dp)
                  end do
                  call fct1d3(ktop(i), kte, dtime_max, po_cup(i, :), qo(i, :), massflx(i, :), trcflx_in(1, :), sub_tend(1, :))
               end if

               !--- add the contribuition from the environ subsidence
               dellaq(i, kts:ktop(i)) = dellaq(i, kts:ktop(i)) + sub_tend(1, kts:ktop(i))

               !--- for output only
               subten_q(i, kts:ktop(i)) = sub_tend(1, kts:ktop(i))

               !- check H and water conservation liq+condensed (including rainfall)
               trash = 0.
               trash2 = 0.0
               do k = kts, ktop(i)
                  dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))
                  g_rain = 0.5*(pwo(i, k) + pwo(i, k + 1))*c_grav/dp
                  e_dn = -0.5*(pwdo(i, k) + pwdo(i, k + 1))*c_grav/dp*edto(i)
                  trash = trash + (dellaq(i, k) + dellaqc(i, k) + g_rain - e_dn)*dp/c_grav
                  trash2 = trash2 + dellah(i, k)*c_grav/dp + c_xlf*((1.-p_liq_ice(i, k))*0.5*(qrco(i, k + 1) + qrco(i, k)) &
                         - melting(i, k))*c_grav/dp
               end do   ! k
            end do
         end if ! vertical discretization formulation

         !--- apply environmental subsidence on grid-scale ice and liq water contents, and cloud fraction (Upwind scheme)
         if (APPLY_SUB_MP == 1) then
            dellampqi = 0.
            dellampql = 0.
            dellampcf = 0.

            do i = its, itf
               if (ierr(i) /= 0) cycle
               do k = kts, ktop(i)
                  dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))

                  !--- apply environmental subsidence on grid-scale/anvil ice and liq water contents (Upwind scheme)
                  !
                  env_mf = -0.5*(zenv(i, k + 1) + zenv(i, k))
                  env_mf_m = min(env_mf, 0.)*c_grav/dp
                  env_mf_p = max(env_mf, 0.)*c_grav/dp

                  dellampqi(:, i, k) = -(env_mf_m*(mpqi(:, i, k + 1) - mpqi(:, i, k)) + env_mf_p*(mpqi(:, i, k) &
                                     - mpqi(:, i, max(k - 1, kts))))
                  dellampql(:, i, k) = -(env_mf_m*(mpql(:, i, k + 1) - mpql(:, i, k)) + env_mf_p*(mpql(:, i, k) &
                                     - mpql(:, i, max(k - 1, kts))))

                  !--- apply environmental subsidence on grid-scale/anvil cloud fraction
                  dellampcf(:, i, k) = -(env_mf_m*(mpcf(:, i, k + 1) - mpcf(:, i, k)) + env_mf_p*(mpcf(:, i, k) &
                                     - mpcf(:, i, max(k - 1, kts))))
               end do

               !--- apply environmental subsidence on grid-scale and anvil cloud fraction using time implicit/explict method
               if (ALP1 > 0.) then
                  alp0 = 1.0 - ALP1
                  do k = kts, ktop(i)
                     dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))
                     env_mf = -0.5*(zenv(i, k + 1) + zenv(i, k))
                     env_mf_m = min(env_mf, 0.)*c_grav/dp
                     env_mf_p = max(env_mf, 0.)*c_grav/dp

                     beta1 = -env_mf_m
                     beta2 = -env_mf_p

                     aa(k) = ALP1*beta2             ! coef of f(k-1,t+1),
                     bb(k) = 1.+ALP1*beta1 - ALP1*beta2  ! coef of f(k  ,t+1),
                     cc(k) = -ALP1*beta1             ! coef of f(k+1,t+1),

                     !-- this is the rhs of the discretization
                     dd(:, k) = (1.-alp0*beta1 + alp0*beta2)*mpcf(:, i, k) + alp0*beta1*mpcf(:, i, k + 1) - alp0*beta2  &
                              * mpcf(:, i, max(kts, k - 1)) ! coef of  f(k-1,t),
                  end do
                  do kmp = 1, nmp
                     !-- this routine solves the problem: aa*f(k-1,t+1) + bb*f(k,t+1) + cc*f(k+1,t+1) = dd
                     call tridiag(ktop(i), aa(kts:ktop(i)), bb(kts:ktop(i)), cc(kts:ktop(i)), dd(kmp, kts:ktop(i)))

                     dellampcf(kmp, i, kts:ktop(i)) = dd(kmp, kts:ktop(i)) - mpcf(kmp, i, kts:ktop(i))
                  end do
               end if
            end do
         end if

         !--- make the smoothness procedure
         if (USE_SMOOTH_TEND > 0) then
            do i = its, itf
               if (ierr(i) /= 0) cycle
               tend2d = 0.

               do k = kts, ktop(i)
                  rcount = 1.e-8
                  tend1d = 0.
                  do kk = max(kts, k - USE_SMOOTH_TEND), min(ktop(i), k + USE_SMOOTH_TEND)
                     dp = (po_cup(i, kk) - po_cup(i, kk + 1))
                     rcount = rcount + dp
                     tend1d(1) = tend1d(1) + dp*dellah(i, kk)
                     tend1d(2) = tend1d(2) + dp*dellaq(i, kk)
                     tend1d(3) = tend1d(3) + dp*dellaqc(i, kk)
                     tend1d(4) = tend1d(4) + dp*dellu(i, kk)
                     tend1d(5) = tend1d(5) + dp*dellv(i, kk)
                  end do
                  tend2d(k, 1:5) = tend1d(1:5)/rcount
               end do
               !--- get the final/smoother tendencies
               do k = kts, ktop(i)
                  dellah(i, k) = tend2d(k, 1)
                  dellaq(i, k) = tend2d(k, 2)
                  dellaqc(i, k) = tend2d(k, 3)
                  dellu(i, k) = tend2d(k, 4)
                  dellv(i, k) = tend2d(k, 5)
               end do
            end do
         end if ! USE_SMOOTH_TEND == 1

         !--- using dellas, calculate changed environmental profiles
         do k = kts, ktf
            do i = its, itf
               dellat(i, k) = 0.
               if (ierr(i) /= 0) cycle
               !
               xhe(i, k) = (dellah(i, k))*mbdt(i) + heo(i, k)
               xq(i, k) = (dellaq(i, k) + dellaqc(i, k))*mbdt(i) + qo(i, k)
               if (xq(i, k) .le. 0.) xq(i, k) = 1.e-08

               !- do not feed dellat with dellaqc if the detrainment of liquid water
               !- will be used as a source for cloud microphysics
               if (p_coupl_mphysics) then
                  dellat(i, k) = (1./real(c_cp))*(dellah(i, k) - real(c_alvl)*dellaq(i, k))
               else
                  !---meltglac-------------------------------------------------
                  dellat(i, k) = (1./real(c_cp))*(dellah(i, k) - real(c_alvl)*(dellaq(i, k) + dellaqc(i, k)) &
                               * (1.+(c_xlf/real(c_alvl))*(1.-p_liq_ice(i, k))))
                  !DELLAT (I,K)=(1./cp)*( DELLAH(I,K)  -xlv*(DELLAQ(I,K) + DELLAQC(i,k)))

                  !-adding dellaqc to dellaq:
                  dellaq(i, k) = dellaq(i, k) + dellaqc(i, k)
                  dellaqc(i, k) = 0.0
               end if
               !---meltglac-------------------------------------------------
               xt(i, k) = ((1./real(c_cp))*dellah(i, k) - (real(c_alvl)/real(c_cp))*(dellaq(i, k) + dellaqc(i, k) &
                        * (1.+(c_xlf/real(c_alvl))*(1.-p_liq_ice(i, k)))))*mbdt(i) + tn(i, k)
               !XT(I,K)=((1./cp)*DELLAH(i,k)-(xlv/cp)*(DELLAQ(i,k)+DELLAQC(i,k)))*MBDT(i)+TN(I,K)

               !--- temp tendency due to the environmental subsidence
               subten_t(i, k) = (1./real(c_cp))*(subten_h(i, k) - real(c_alvl)*subten_q(i, k))
            end do
         end do
         do i = its, itf
            if (ierr(i) /= 0) cycle
            !XHKB(I)=(dsubh(i,k22(i))+DELLAH(I,K22(i)))*MBDT+HKBO(I)
            xhe(i, ktf) = heo(i, ktf)
            xq(i, ktf) = qo(i, ktf)
            xt(i, ktf) = tn(i, ktf)
            if (xq(i, ktf) .le. 0.) xq(i, ktf) = 1.e-08
         end do
         !- new way for defining XHKB
         do i = its, itf
            if (ierr(i) /= 0) cycle
            !XHKB(I)= DELLAH(I,K22(i))*MBDT+HKBO(I)
            !-note that HKBO already contains the contribuition from
            !-ztexec and zqexec
            call getCloudBc(cumulus, kts, kte, ktf, xland(i), po(i, kts:kte), dellah(i, kts:kte), dellah_aver, k22(i))
            xhkb(i) = dellah_aver*mbdt(i) + hkbo(i)
         end do

         !--- calculate moist static energy, heights, qes
         call cupEnv(xz, xqes, xhe, xhes, xt, xq, po, z1, psur, ierr, -1, itf, ktf, its, ite, kts, kte)

         !--- environmental values on cloud levels
         call cupEnvCLev(xt, xqes, xq, xhe, xhes, xz, po, xqes_cup, xq_cup, xhe_cup, us, vs, u_cup, v_cup, xhes_cup, xz_cup &
                       , po_cup, gamma_cup, xt_cup, psur, tsur, ierr, z1, itf, ktf, its, ite, kts, kte)
         !
         !--- static control
         !
         !--- moist static energy inside cloud
         !
         do i = its, itf
            xhc(i, :) = 0.
            if (ierr(i) /= 0) cycle
            do k = kts, start_level(i) !k22(i)
               xhc(i, k) = xhkb(i)
            end do
         end do
         !
         !--- option to produce linear fluxes in the sub-cloud layer.
         if (trim(cumulus) == 'shallow' .and. USE_LINEAR_SUBCL_MF == 1) then
            do i = its, itf
               if (ierr(i) /= 0) cycle
               call getDelmix(cumulus, kts, kte, ktf, xland(i), start_level(i), po(i, kts:kte), xhe_cup(i, kts:kte) &
                            , xhc(i, kts:kte))
            end do
         end if
         do i = its, itf
            if (ierr(i) /= 0) cycle
            do k = start_level(i) + 1, ktop(i) + 1  ! mass cons option
               denom = (xzu(i, k - 1) - .5*up_massdetro(i, k - 1) + up_massentro(i, k - 1))
               if (denom == 0.0) then
                  xhc(i, k) = xhc(i, k - 1)
               else
                  xhc(i, k) = (xhc(i, k - 1)*xzu(i, k - 1) - .5*up_massdetro(i, k - 1)*xhc(i, k - 1) + up_massentro(i, k - 1) &
                            * xhe(i, k - 1))/denom
                  if (k == start_level(i) + 1) then
                     x_add = (real(c_alvl)*zqexec(i) + real(c_cp)*ztexec(i)) + x_add_buoy(i)
                     xhc(i, k) = xhc(i, k) + x_add*up_massentro(i, k - 1)/denom
                  end if
               end if
               !
               !- include glaciation effects on XHC
               !                                   ------ ice content --------
               xhc(i, k) = xhc(i, k) + c_xlf*(1.-p_liq_ice(i, k))*qrco(i, k)
            end do
            do k = ktop(i) + 2, ktf
               xhc(i, k) = xhes_cup(i, k)
               xzu(i, k) = 0.
            end do
         end do
         call getBuoyancy(itf, ktf, its, ite, kts, kte, ierr, klcl, kbcon, ktop, xhc, xhe_cup, xhes_cup, xdby, xz_cup)
         !
         !--- workfunctions for updraft
         !
         call cupUpAa0(xaa0, xz_cup, xzu, xdby, GAMMA_CUP, xt_cup, k22, klcl, kbcon, ktop, ierr, itf, ktf, its, ite, kts, kte)

         do nens = 1, p_maxens
            do i = its, itf
               if (ierr(i) /= 0) cycle
               !~ xaa0_ens(i,nens)=xaa0(i)
               do k = kts, ktop(i)
                  do nens3 = 1, p_maxens3
                     if (nens3 .eq. 7) then
                        !--- b=0
                        pr_ens(i, nens3) = pr_ens(i, nens3) + pwo(i, k) + edto(i)*pwdo(i, k)
                        !--- b=beta
                     else if (nens3 .eq. 8) then
                        pr_ens(i, nens3) = pr_ens(i, nens3) + pwo(i, k) + edto(i)*pwdo(i, k)
                        !--- b=beta/2
                     else if (nens3 .eq. 9) then
                        pr_ens(i, nens3) = pr_ens(i, nens3) + pwo(i, k) + edto(i)*pwdo(i, k)
                     else
                        pr_ens(i, nens3) = pr_ens(i, nens3) + pwo(i, k) + edto(i)*pwdo(i, k)
                     end if
                  end do
               end do
               if (pr_ens(i, 7) .lt. 1.e-6 .and. C0_MID > 0. .and. trim(cumulus) /= 'shallow') then
                  ierr(i) = 18
                  ierrc(i) = "total normalized condensate too small"
                  do nens3 = 1, p_maxens3
                     pr_ens(i, nens3) = 0.
                  end do
               end if
               do nens3 = 1, p_maxens3
                  if (pr_ens(i, nens3) < 1.e-5) pr_ens(i, nens3) = 0.
               end do
            end do
         end do
         !
         !--- LARGE SCALE FORCING
         !
         !do i=its,itf
         !   ierr2(i)=ierr(i)
         !   ierr3(i)=ierr(i)
         !enddo
         !
         !--- calculate cloud base mass flux
         !
         if (trim(cumulus) == 'deep') &
            call cupForcingEns3d(itf, ktf, its, ite, kts, kte, p_ens4, p_ensdim, ichoice, p_maxens, p_maxens2, p_maxens3 &
                                    , ierr, ierr2, ierr3, k22, kbcon, ktop, xland1, aa0, aa1, xaa0, mbdt, dtime &
                                    , xf_ens, mconv, qo, po_cup, omeg, zdo, zuo, pr_ens, edto, tau_ecmwf, aa1_bl, xf_dicycle &
                                    , xk_x, alpha_adv, q_adv, aa1_radpbl, aa1_adv, wlpool_bcon, xf_coldpool)

         if (trim(cumulus) == 'mid') &
            call cupForcingEns3dMid(aa0, aa1, xaa0, mbdt, dtime, ierr, po_cup, ktop, k22, kbcon, kpbl, ichoice, p_maxens &
                                  , p_maxens3, itf, ktf, its, ite, kts, kte, tau_ecmwf, aa1_bl, xf_dicycle, dhdt, xff_mid, zws, hc &
                                  , hco, he_cup, heo_cup, wlpool_bcon, xf_coldpool)

         if (trim(cumulus) == 'shallow') then
            call cupUpCape(cape, z, zu, dby, gamma_cup, t_cup, k22, kbcon, ktop, ierr, tempco, qco, qrco, qo_cup, itf, ktf, its &
                         , ite, kts, kte)

            call cupForcingEns3dShal(itf, ktf, its, ite, kts, kte, dtime, ichoice, ierrc, ierr, klcl, kpbl, kbcon, k22, ktop &
                                         , xmb, tsur, cape, h_sfc_flux, le_sfc_flux, zws, po, hco, heo_cup, po_cup, t_cup, dhdt &
                                         , rho, xff_shal, xf_dicycle, tke_pbl, wlpool_bcon, xf_coldpool)
         end if
         !
         !
         !--- get the net precipitation at surface
         !
         do i = its, itf
            if (ierr(i) == 0) then
               pwo_eff(i, :) = pwo(i, :) + edto(i)*pwdo(i, :)
            else
               pwo_eff(i, :) = 0.
            end if
         end do

      end do

      !--- Include kinetic energy dissipation converted to heating
      call keToHeating(itf, ktf, its, ite, kts, kte, ktop, ierr, po_cup, us, vs, dellu, dellv, dellat)

      !--- feedback
      call cupOutputEns3d(cumulus, xff_shal, xff_mid, xf_ens, ierr, dellat, dellaq, dellaqc, outt, outq, outqc, zuo, pre, pwo_eff &
                        , xmb, ktop, p_maxens2, p_maxens, ierr2, ierr3, pr_ens, p_maxens3, p_ensdim, sig, xland1, ichoice, ipr &
                        , jpr, itf, ktf, its, ite, kts, kte, xf_dicycle, outu, outv, dellu, dellv, dtime, po_cup, kbcon, dellabuoy &
                        , outbuoy, dellampqi, outmpqi, dellampql, outmpql, dellampcf, outmpcf, nmp, rh_dicycle_fct, xf_coldpool &
                        , wlpool_bcon)

      !--- get the net precipitation flux (after downdraft evaporation)
      call getPrecipFluxes(cumulus, klcl, kbcon, ktop, k22, ierr, xland, pre, xmb, pwo, pwavo, edto, pwevo, pwdo, t_cup, tempco &
                        ,  prec_flx, evap_flx, itf, ktf, its, ite, kts, kte)

      !--- rainfall evap below cloud base
      if (USE_REBCB == 1) &
         call rainEvapBelowCloudBase(cumulus, itf, ktf, its, ite, kts, kte, ierr, kbcon, ktop, xmb, psur, xland, qo_cup, t_cup &
                                   , po_cup, qes_cup, pwavo, edto, pwevo, pwo, pwdo, pre, prec_flx, evap_flx, outt, outq &
                                   , outbuoy, evap_bcb)

      !--- includes effects of the remained cloud dissipation into the enviroment
      if (USE_CLOUD_DISSIPATION >= 0.) &
         call cloudDissipation(cumulus, itf, ktf, its, ite, kts, kte, ierr, kbcon, ktop, dtime, xmb, xland, qo_cup, qeso_cup &
                             , po_cup, outt, outq, outqc, zuo, vvel2d, rho_hydr, qrco, sig, tempco, qco, tn_cup, heso_cup, zo)

      !--- get the total (deep+congestus) evaporation flux for output (units kg/kg/s)
      do i = its, itf
         if (ierr(i) /= 0) cycle
         do k = kts, ktop(i)
            dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))
            !--- add congestus and deep plumes, and convert to kg/kg/s
            revsu_gf(i, k) = revsu_gf(i, k) + evap_flx(i, k)*c_grav/dp
         end do
      end do

      !--- get lightning flashes density (parameterization from Lopez 2016, MWR)
      if (LIGHTNING_DIAG == 1 .and. trim(cumulus) == 'deep') then
         call cupUpCape(cape, z, zu, dby, gamma_cup, t_cup, k22, kbcon, ktop, ierr, tempco, qco, qrco, qo_cup, itf, ktf, its, ite &
                      , kts, kte)

         call cupUpLightning(itf, ktf, its, ite, kts, kte, ierr, kbcon, ktop, xland, cape, zo, zo_cup, t_cup, t, tempco, qrco &
                           , po_cup, rho, prec_flx, lightn_dens)
      end if

      !--- for outputs (only deep plume)
      if (trim(cumulus) == 'deep') then
         do i = its, itf
            if (ierr(i) /= 0) cycle
            var2d(i) = p_cwv_ave(i)
            do k = kts, ktop(i) + 1
               prfil_gf(i, k) = prec_flx(i, k)
               var3d_agf(i, k) = vvel2d(i, k)
            end do
         end do
      end if

      !--- for tracer convective transport / outputs
      do i = its, itf
         if (ierr(i) /= 0) cycle
         do k = kts, ktf
            !clwup5d     (i,k) = qrco (i,k) !ice/liquid water
            !tup         (i,k) = (1./cp)*(hco(i,k)-g*zo_cup(i,k)-xlv*qco(i,k))!in-updraft temp
            tup(i, k) = tempco(i, k) !in-updraft temp
         end do
         tup(i, kte) = t_cup(i, kte)
      end do

      !--- convert mass fluxes, etc...
      do i = its, itf
         if (ierr(i) /= 0) cycle
         pwavo(i) = xmb(i)*pwavo(i)
         pwevo(i) = xmb(i)*pwevo(i)
         zuo(i, :) = xmb(i)*zuo(i, :)
         zdo(i, :) = xmb(i)*zdo(i, :)
         pwo(i, :) = xmb(i)*pwo(i, :)
         pwdo(i, :) = xmb(i)*pwdo(i, :)
         up_massentro(i, :) = xmb(i)*up_massentro(i, :)
         up_massdetro(i, :) = xmb(i)*up_massdetro(i, :)
         dd_massentro(i, :) = xmb(i)*dd_massentro(i, :)
         dd_massdetro(i, :) = xmb(i)*dd_massdetro(i, :)
         zenv(i, :) = xmb(i)*zenv(i, :)
      end do

      !--for output only.
      do i = its, itf
         subten_q(i, :) = xmb(i)*subten_q(i, :)
         subten_h(i, :) = xmb(i)*subten_h(i, :)
         subten_t(i, :) = xmb(i)*subten_t(i, :)
      end do

      !--- outputs a model sounding for the stand-alone code (part 2)
      if (output_sound == 1) then
         call sound(2, cumulus, int_time, dtime, p_ens4, itf, ktf, its, ite, kts, kte, xlats, xlons, jcol, whoami_all &
                    , z, qes, he, hes, t, q, po, z1, psur, zo, qeso, heo, heso, tn, qo, us, vs, omeg, xz &
                    , h_sfc_flux, le_sfc_flux, tsur, dx, stochastic_sig, zws, ztexec, zqexec, xland &
                    , kpbl, k22, klcl, kbcon, ktop, aa0, aa1, sig, xaa0, hkb, xmb, pre, edto &
                    , zo_cup, dhdt, rho, zuo, zdo, up_massentro, up_massdetro, outt, outq, outqc, outu, outv)
      end if

      !--- for output only
      if (trim(cumulus) == 'deep') then
         aa1_(:) = aa1(:)
         aaa0_(:) = aa0(:)
         aa1_radpbl_(:) = aa1_radpbl(:)

         if (DICYCLE == 2) then
            aa1_adv_(:) = q_adv(:)
         else
            ! AA1_ADV_ (:) = AA1_ADV    (:)
            ! AA1_ADV_ (:) = wlpool_bcon(:)
            aa1_adv_(:) = vshear(:)
            ! AA1_ADV_ (:) = depth_neg_buoy (:)
            ! AA1_ADV_ (:) = cin1 (:)
         end if
         do i = its, itf
            if (ierr(i) == 0) cycle
            kbcon(i) = 1
            ktop(i) = 1
            klcl(i) = 1
            jmin(i) = 1
            k22(i) = 1
         end do
      end if

      if (LIQ_ICE_NUMBER_CONC == 1) then
         call getLiqIceNumberConc(itf, ktf, its, ite, kts, kte, ierr, ktop, dtime, rho, outqc, tempco, outnliq, outnice)
      end if

      !- section for atmospheric composition
      if (USE_TRACER_TRANSP == 1) then

         !--only for debug
         if (p_use_gate) then
            if (jl == 1) then
               se_chem_update(1, :, :) = se_chem(1, :, :)
            else
               se_chem(1, :, :) = se_chem_update(1, :, :)
            end if
            do i = its, itf
               if (ierr(i) /= 0) cycle
               massi = 0.
               do k = kts, ktop(i)
                  dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))
                  massi = massi + se_chem(1, i, k)*dp/c_grav
               end do
            end do
         end if
         !--only for debug

         !-1) get mass mixing ratios at the cloud levels
         call cupEnvClevChem(mtp, se_chem, se_cup_chem, ierr, itf, ktf, its, ite, kts, kte)

         !-2) determine in-cloud tracer mixing ratios
         !
         ! a) chem - updraft
         !- note: here "sc_up_chem" stores the total in-cloud tracer mixing ratio (i.e., including the portion
         !        embedded in the condensates).
         call getInCloudScChemUp(cumulus, fscav, mtp, se_chem, se_cup_chem, sc_up_chem, pw_up_chem, tot_pw_up_chem, zo_cup, rho &
                               , po, po_cup, qrco, tempco, pwo, zuo, up_massentro, up_massdetro, vvel2d, vvel1d, start_level, k22 &
                               , kbcon, ktop, klcl, ierr, xland, itf, ktf, its, ite, kts, kte)

         ! b) chem - downdraft
         call getInCloudScChemDd(cumulus, fscav, mtp, se_chem, se_cup_chem, sc_dn_chem, pw_dn_chem, pw_up_chem, sc_up_chem &
                               , tot_pw_up_chem, tot_pw_dn_chem, zo_cup, rho, po_cup, qrcdo, pwdo, pwevo, edto, zdo, dd_massentro  &
                               , dd_massdetro, pwavo, pwo, jmin, ierr, itf, ktf, its, ite, kts, kte)
         !
         !-3) determine the vertical transport including mixing, scavenging and evaporation
         !
         !---a) change per unit mass that a model cloud would modify the environment
         do i = its, itf
            if (ierr(i) /= 0) cycle

            !- flux form + source/sink terms + time explicit + FCT
            if (USE_FLUX_FORM == 1 .and. ALP1 == 0.) then

               if (USE_FCT == 0) then
                  do k = kts, ktop(i)
                     dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))

                     out_chem(:, i, k) = -(zuo(i, k + 1)*(sc_up_chem(:, i, k + 1) - se_cup_chem(:, i, k + 1)) - zuo(i, k) &
                                       * (sc_up_chem(:, i, k) - se_cup_chem(:, i, k)))*c_grav/dp + (zdo(i, k + 1) &
                                       * (sc_dn_chem(:, i, k + 1) - se_cup_chem(:, i, k + 1)) - zdo(i, k)*(sc_dn_chem(:, i, k) &
                                       - se_cup_chem(:, i, k)))*c_grav/dp*edto(i)
                  end do

               else

                  !-- FCT scheme for the subsidence transport: d(M_env*S_env)/dz
                  sub_tend = 0.
                  trcflx_in = 0.
                  dtime_max = dtime
                  massflx(i, :) = 0.

                  do k = kts + 1, ktop(i) + 1
                     dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))
                     trcflx_in(:, k) = -(zuo(i, k) - edto(i)*zdo(i, k))*se_cup_chem(:, i, k) !* xmb(i)
                     massflx(i, k) = -(zuo(i, k) - edto(i)*zdo(i, k))           !* xmb(i)
                     dtime_max = min(dtime_max, .5*dp)
                  end do
                  !- if dtime_max<dtime => needs a loop to update from t to t+dtime (check this!)
                  !if( dtime_max < dtime ) stop "dtime_max < dtime in GF scheme"

                  do ispc = 1, mtp
                     call fct1d3(ktop(i), kte, dtime_max, po_cup(i, :), se_chem(ispc, i, :), massflx(i, :), trcflx_in(ispc, :) &
                              ,  sub_tend(ispc, :))
                  end do

                  do k = kts, ktop(i)
                     dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))
                     out_chem(:, i, k) = -(zuo(i, k + 1)*(sc_up_chem(:, i, k + 1)) - zuo(i, k)*(sc_up_chem(:, i, k)))*c_grav/dp &
                                       + (zdo(i, k + 1)*(sc_dn_chem(:, i, k + 1)) - zdo(i, k)*(sc_dn_chem(:, i, k))) *c_grav/dp &
                                       * edto(i)

                     !- update with the subsidence term from FCT scheme
                     out_chem(:, i, k) = out_chem(:, i, k) + sub_tend(:, k)

                  end do
               end if

               !- include evaporation (this term must not be applied to the tracer 'QW')
               if (USE_TRACER_EVAP == 1 .and. trim(cumulus) /= 'shallow') then
                  do k = kts, ktop(i)
                     dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))
                     out_chem(:,i,k) = out_chem(:,i,k) - 0.5*edto(i)*(zdo(i,k)*pw_dn_chem(:,i,k)+zdo(i,k+1) * pw_dn_chem(:,i,k+1)) &
                                     * c_grav/dp !&  ! evaporated ( pw_dn < 0 => E_dn > 0)
                     !*chem_name_mask_evap(:) !-- to avoid the "Dry Mass Violation"
                  end do
               end if

               !- include scavenging
               if (USE_TRACER_SCAVEN > 0 .and. trim(cumulus) /= 'shallow') then
                  do k = kts, ktop(i)
                     dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))
                     out_chem(:, i, k) = out_chem(:, i, k) - 0.5*(zuo(i, k)*pw_up_chem(:, i, k) + zuo(i, k + 1) &
                                       * pw_up_chem(:, i, k + 1))*c_grav/dp  ! incorporated in rainfall (<0)
                  end do
               end if
            end if ! IF(USE_FLUX_FORM == 1 .and. ALP1 == 0. )

            !- flux form + source/sink terms + time explicit/implicit + upstream
            if (USE_FLUX_FORM == 1 .and. ALP1 > 0.) then

               alp0 = 1.-ALP1
               do k = kts, ktop(i) + 1
                  fp(k) = 0.5*(zenv(i, k) + abs(zenv(i, k)))
                  fm(k) = 0.5*(zenv(i, k) - abs(zenv(i, k)))
               end do

               do k = kts, ktop(i)
                  dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))
                  beta1 = dtime*c_grav/dp
                  aa(k) = ALP1*beta1*fm(k)
                  bb(k) = 1.+ALP1*beta1*(fp(k) - fm(k + 1))
                  cc(k) = -ALP1*beta1*fp(k + 1)

                  ddtr(:, k) = se_chem(:, i, k) - (zuo(i, k + 1)*sc_up_chem(:, i, k + 1) - zuo(i, k)*sc_up_chem(:, i, k))*beta1 &
                             + (zdo(i, k + 1)*sc_dn_chem(:, i, k + 1) - zdo(i, k)*sc_dn_chem(:, i, k))*beta1*edto(i)

                  !- include evaporation (this term must not be applied to the tracer 'QW')
                  if (USE_TRACER_EVAP == 1 .and. trim(cumulus) /= 'shallow') then
                     out_chem(:,i,k) = out_chem(:,i,k) - 0.5*edto(i)*(zdo(i,k)*pw_dn_chem(:,i,k)+zdo(i,k+1) * pw_dn_chem(:,i,k+1)) &
                                     * beta1 !&  ! evaporated ( pw_dn < 0 => E_dn > 0)
                     !*chem_name_mask_evap(:) !-- to avoid the "Dry Mass Violation"
                  end if

                  !- include scavenging
                  if (USE_TRACER_SCAVEN > 0 .and. trim(cumulus) /= 'shallow') then
                     out_chem(:, i, k) = out_chem(:, i, k) - 0.5*(zuo(i, k)*pw_up_chem(:, i, k) + zuo(i, k + 1) &
                                       * pw_up_chem(:, i, k + 1))*beta1  ! incorporated in rainfall (<0)
                  end if

                  ddtr(:, k) = ddtr(:, k) + out_chem(:, i, k) + alp0*beta1*(-fm(k)*se_chem(:, i, max(kts, k - 1)) + (fm(k + 1) &
                             - fp(k))*se_chem(:, i, k) + fp(k + 1)*se_chem(:, i, k + 1))

               end do
               do ispc = 1, mtp
                  if (chem_name_mask(ispc) == 0) cycle
                  call tridiag(ktop(i), aa(kts:ktop(i)), bb(kts:ktop(i)), cc(kts:ktop(i)), ddtr(ispc, kts:ktop(i)))
                  out_chem(ispc, i, kts:ktop(i)) = (ddtr(ispc, kts:ktop(i)) - se_chem(ispc, i, kts:ktop(i)))/dtime
               end do
            end if !USE_FLUX_FORM == 1 .and. ALP1 > 0.

            !- flux form + source/sink terms + time explicit + upstream with anti-diffusion step (Smolarkiewicz 1983)
            if (USE_FLUX_FORM == 2 .or. USE_FLUX_FORM == 3) then
               if (USE_FLUX_FORM == 2) lstep = -1 ! upstream + anti-diffusion step
               if (USE_FLUX_FORM == 3) lstep = 1 ! only upstream
               alp0 = 1.

               if (ierr(i) /= 0) cycle
               !--- Zenv here have the following reference:  < 0 => downward motion
               zenv(i, :) = -(zuo(i, :) - edto(i)*zdo(i, :))

               do istep = 1, lstep, -2

                  if (istep == 1) then
                     ddtr_upd(:, :) = se_chem(:, i, :)
                     do k = kts, ktop(i) + 1
                        fp_mtp(:, k) = 0.5*(zenv(i, k) + abs(zenv(i, k)))
                        fm_mtp(:, k) = 0.5*(zenv(i, k) - abs(zenv(i, k)))
                     end do
                  else
                     ddtr_upd(:, kts:ktop(i) + 1) = se_chem(:, i, kts:ktop(i) + 1) + out_chem(:, i, kts:ktop(i) + 1)*dtime
                     zenv_diff(:, kts) = 0.
                     do k = kts, ktop(i) + 1
                        dz = zo_cup(i, k + 1) - zo_cup(i, k)
                        zenv_diff(:, k + 1) = 1.08*(dz*abs(zenv(i, k + 1)) - dtime*zenv(i, k + 1)**2) &
                                              *(ddtr_upd(:, k + 1) - ddtr_upd(:, k)) &
                                              /((ddtr_upd(:, k + 1) + ddtr_upd(:, k) + 1.e-16)*dz)
                     end do
                     do k = kts, ktop(i) + 1
                        fp_mtp(:, k) = 0.5*(zenv_diff(:, k) + abs(zenv_diff(:, k)))
                        fm_mtp(:, k) = 0.5*(zenv_diff(:, k) - abs(zenv_diff(:, k)))
                     end do
                  end if

                  do k = kts, ktop(i)
                     dp = -100.*(po_cup(i, k) - po_cup(i, k + 1))
                     beta1 = dtime*c_grav/dp
                     ddtr(:, k) = ddtr_upd(:, k) + alp0*beta1*((fp_mtp(:, k + 1)*ddtr_upd(:, k) + fm_mtp(:, k + 1) &
                                * ddtr_upd(:, k + 1)) - (fp_mtp(:, k)*ddtr_upd(:, max(kts, k - 1)) + fm_mtp(:, k)*ddtr_upd(:, k)))
                  end do
                  do ispc = 1, mtp
                     if (chem_name_mask(ispc) == 0) cycle
                     out_chem(ispc, i, kts:ktop(i)) = (ddtr(ispc, kts:ktop(i)) - se_chem(ispc, i, kts:ktop(i)))/dtime
                  end do

               end do ! anti-diff steps

               do k = kts, ktop(i)
                  dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))
                  beta1 = c_grav/dp

                  out_chem(:, i, k) = out_chem(:, i, k) - (zuo(i, k + 1)*sc_up_chem(:, i, k + 1) - zuo(i, k)*sc_up_chem(:, i, k)) &
                                    * beta1 + (zdo(i, k + 1)*sc_dn_chem(:, i, k + 1) - zdo(i, k)*sc_dn_chem(:, i, k))*beta1*edto(i)

                  !- include evaporation (this term must not be applied to the tracer 'QW')
                  if (USE_TRACER_EVAP == 1 .and. trim(cumulus) /= 'shallow') then
                     out_chem(:,i,k) = out_chem(:,i,k) - 0.5*edto(i)*(zdo(i,k)*pw_dn_chem(:,i,k)+zdo(i,k+1) * pw_dn_chem(:,i,k+1)) &
                                     * beta1 !&  ! evaporated ( pw_dn < 0 => E_dn > 0)
                     !*chem_name_mask_evap(:) !-- to avoid the "Dry Mass Violation"
                  end if

                  !- include scavenging
                  if (USE_TRACER_SCAVEN > 0 .and. trim(cumulus) /= 'shallow') then
                     out_chem(:, i, k) = out_chem(:, i, k) - 0.5*(zuo(i, k)*pw_up_chem(:, i, k) + zuo(i, k + 1) &
                                       * pw_up_chem(:, i, k + 1))*beta1  ! incorporated in rainfall (<0)
                  end if
               end do
            end if ! USE_FLUX_FORM == 2 .or. USE_FLUX_FORM == 3

            !--- check mass conservation for tracers
            do ispc = 1, mtp
               if (chem_name_mask(ispc) == 0) cycle
               trash_(:) = 0.
               trash2_(:) = 0.
               evap_(:) = 0.
               wetdep_(:) = 0.
               residu_(:) = 0.
               do k = kts, ktop(i)
                  dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))
                  evap = -0.5*(zdo(i, k)*pw_dn_chem(ispc, i, k) + zdo(i, k + 1)*pw_dn_chem(ispc, i, k + 1))*c_grav/dp*edto(i)
                  wetdep = 0.5*(zuo(i, k)*pw_up_chem(ispc, i, k) + zuo(i, k + 1)*pw_up_chem(ispc, i, k + 1))*c_grav/dp

                  evap_(ispc) = evap_(ispc) + evap*dp/c_grav
                  wetdep_(ispc) = wetdep_(ispc) + wetdep*dp/c_grav
                  residu_(ispc) = residu_(ispc) + (wetdep - evap)*dp/c_grav

                  !trash_ (ispc) =   trash_ (ispc) + (out_chem (ispc,i,k) - evap + wetdep)*dp/g
                  trash_(ispc) = trash_(ispc) + (out_chem(ispc, i, k))*dp/c_grav

                  trash2_(ispc) = trash2_(ispc) + se_chem(ispc, i, k)*dp/c_grav
               end do
               if (residu_(ispc) < 0.) then
                  beta1 = c_grav/(po_cup(i, kts) - po_cup(i, ktop(i) + 1))
                  do k = kts, ktop(i)
                     out_chem(ispc, i, k) = out_chem(ispc, i, k) + residu_(ispc)*beta1
                  end do
               end if

               !if(evap_  (ispc) > wetdep_(ispc)) then
               !print*,"budget=",ispc,evap_  (ispc), wetdep_(ispc),trash_ (ispc),trim(CHEM_NAME(ispc))!,trash_ (ispc),trash2_(ispc)
               !call flush(6)
               !endif
               !if(evap_  (ispc) > wetdep_(ispc)) stop " eva<wet "
               !if(abs(trash_(ispc)) >1.e-6 ) then
               !  if (MAPL_AM_I_ROOT())  write(6,*)'=> mass_cons=',trash_(ispc),spacing(trash2_(ispc)),trim(CHEM_NAME(ispc)),trim(cumulus)
               !endif
            end do

         end do ! loop 'i'

         if (p_use_gate) then
            !--only for debug
            do i = its, itf
               if (ierr(i) /= 0) cycle
               massf = 0.
               do k = kts, ktop(i)
                  se_chem_update(ispc_co, i, k) = se_chem_update(ispc_co, i, k) + out_chem(ispc_co, i, k)*dtime
                  !se_chem_update(ispc_CO,i,k) = max(0.,se_chem_update(ispc_CO,i,k))
                  dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))
                  evap_(ispc_co) = -0.5*(zdo(i, k)*pw_dn_chem(ispc_co, i, k) + zdo(i, k + 1)*pw_dn_chem(ispc_co, i, k + 1)) &
                                 * c_grav/dp*edto(i)
                  wetdep_(ispc_co) = 0.5*(zuo(i, k)*pw_up_chem(ispc_co, i, k) + zuo(i, k + 1)*pw_up_chem(ispc_co, i, k + 1)) &
                                   * c_grav/dp
                  massf = massf + se_chem_update(1, i, k)*dp/c_grav + (-evap_(ispc_co) + wetdep_(ispc_co))*dp/c_grav
               end do
               if (abs((massf - massi)/(1.e-12 + massi)) > 1.e-6) print *, "mass con=>", (massf - massi)/(1.e-12 + massi)
            end do
19          format(1x, I3, 1x, 5e14.3)
18          format(1x, I3, 1x, 4e14.3)
20          format(1x, I3, 1x, 11e16.6)
            !--only for debug
         end if
      end if !- end of section for atmospheric composition
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
      if (trim(cumulus) == 'deep') then
         do i = its, itf
            !if(ierr(i) /= 0) cycle
            aaa0_(i) = aa0(i)
            aa1_(i) = aa1(i)
            aa1_bl_(i) = aa1_bl(i)
            tau_bl_(i) = tau_bl(i)
            tau_ec_(i) = tau_ecmwf(i)
            !  if(USE_MEMORY == 0) tau_ec_ (i)  = x_add_buoy(i)
         end do
      end if

      !- begin: for GATE soundings-------------------------------------------
      if (p_use_gate .or. wrtgrads) then
         if (trim(cumulus) == 'deep') then
            cty = '1'
            nvarbegin = 0
         end if
         if (trim(cumulus) == 'shallow') then
            cty = '2'
            nvarbegin = 101
         end if
         if (trim(cumulus) == 'mid') then
            cty = '3'
            nvarbegin = 201
         end if
         do i = its, itf
            !if(ierr(i).eq.0) then
            !- 2-d section
            do k = kts, ktf  !max(1,ktop(i))
               nvar = nvarbegin

               if (trim(cumulus) == 'deep') &
                  call setGradsVar(jl, k, nvar, zo(i, k), "zo"//cty, ' height', '3d')
               !        call set_grads_var(jl,k,nvar,po(i,k),"po"//cty ,' press','3d')

               dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))
               e_dn = -0.5*(pwdo(i, k) + pwdo(i, k + 1))*c_grav/dp*edto(i)*86400.*real(c_alvl)/real(c_cp)*xmb(i) ! pwdo < 0 and E_dn must > 0
               c_up = dellaqc(i, k) + (zuo(i, k + 1)*qrco(i, k + 1) - zuo(i, k)*qrco(i, k))*c_grav/dp + 0.5*(pwo(i, k) &
                    + pwo(i, k + 1))*c_grav/dp
               c_up = -c_up*86400.*real(c_alvl)/real(c_cp)*xmb(i)

               trash = -(zuo(i, k + 1)*(qco(i, k + 1) - qo_cup(i, k + 1)) - zuo(i, k)*(qco(i, k) - qo_cup(i, k)))*c_grav/dp
               trash2 = +(zdo(i, k + 1)*(qcdo(i, k + 1) - qo_cup(i, k + 1)) - zdo(i, k)*(qcdo(i, k) - qo_cup(i, k)))*c_grav/dp &
                      * edto(i)

               trash = trash*86400.*real(c_alvl)/real(c_cp)*xmb(i)
               trash2 = trash2*86400.*real(c_alvl)/real(c_cp)*xmb(i)

               env_mf = 0.5*((zuo(i, k + 1) - zdo(i, k + 1)*edto(i)) + (zuo(i, k) - zdo(i, k)*edto(i)))
               resten_h = dellah(i, k) - subten_h(i, k)
               resten_q = dellaq(i, k) - subten_q(i, k)
               resten_t = (1./real(c_cp))*(resten_h - real(c_alvl)*resten_q)
               !trash2 = qco   (i,k  )! zuo(i,k+1)*(qco (i,k+1)-qo_cup(i,k+1) ) !*g/dp
               !trash  = qo_cup(i,k  )! zuo(i,k  )*(qco (i,k  )-qo_cup(i,k  ) ) !*g/dp
               trash2 = zuo(i, k + 1)*(qco(i, k + 1) - qo_cup(i, k + 1))*1000 !*g/dp
               trash = zuo(i, k)*(qco(i, k) - qo_cup(i, k))*1000  !*g/dp

               call setGradsVar(jl, k, nvar, out_chem(1, i, k)*86400, "outchem"//cty, ' outchem', '3d')
               call setGradsVar(jl, k, nvar, sc_up_chem(1, i, k), "scup"//cty, ' sc_chem', '3d')
               call setGradsVar(jl, k, nvar, sc_dn_chem(1, i, k), "scdn"//cty, ' sc_chem', '3d')
               call setGradsVar(jl, k, nvar, massi, "mi"//cty, ' initial mass', '2d')
               call setGradsVar(jl, k, nvar, massf, "mf"//cty, ' final mass', '2d')
               call setGradsVar(jl, k, nvar, se_chem(1, i, k), "se"//cty, ' se_chem', '3d')
               call setGradsVar(jl, k, nvar, se_cup_chem(1, i, k), "secup"//cty, ' se_cup_chem', '3d')
               !-- only for debug
               !call set_grads_var(jl,k,nvar,se_chem_update(1,i,k),"newse"//cty ,' new se_chem','3d')
               if (APPLY_SUB_MP == 1) then
                  kmp = p_lsmp
                  call setGradsVar(jl, k, nvar, outmpqi(kmp, i, k)*86400*1000, "outqi"//cty, ' outmpqi', '3d')
                  call setGradsVar(jl, k, nvar, outmpql(kmp, i, k)*86400*1000, "outql"//cty, ' outmpql', '3d')
                  call setGradsVar(jl, k, nvar, outmpcf(kmp, i, k)*86400, "outcf"//cty, ' outmpcf', '3d')
                  call setGradsVar(jl, k, nvar, mpqi(kmp, i, k), "mpqi"//cty, ' mpqi', '3d')
                  call setGradsVar(jl, k, nvar, mpql(kmp, i, k), "mpql"//cty, ' mpql', '3d')
                  call setGradsVar(jl, k, nvar, mpcf(kmp, i, k), "mpcf"//cty, ' mpcf', '3d')
               end if
               call setGradsVar(jl, k, nvar, env_mf, "sub"//cty, ' sub', '3d')
               if (LIQ_ICE_NUMBER_CONC == 1) then
                  call setGradsVar(jl, k, nvar, outnice(i, k)*86400., "outnice"//cty, 'out # ice1/day', '3d')
                  call setGradsVar(jl, k, nvar, outnliq(i, k)*86400., "outnliq"//cty, 'out # liq /day', '3d')
               end if
               call setGradsVar(jl, k, nvar, zuo(i, k)/xmb(i), "zup"//cty, 'norm m flux up ', '3d')
               call setGradsVar(jl, k, nvar, zdo(i, k)/xmb(i), "zdn"//cty, 'norm m flux dn ', '3d')
               call setGradsVar(jl, k, nvar, zenv(i, k), "zenv"//cty, 'norm m flux env ', '3d')
               call setGradsVar(jl, k, nvar, -edto(i)*xmb(i)*zdo(i, k), "mdn"//cty, 'm flux down (kg/s/m^2)', '3d')
               call setGradsVar(jl, k, nvar, up_massentro(i, k), "upent"//cty, 'up_massentr(kg/s/m^2)', '3d')
               call setGradsVar(jl, k, nvar, xmb(i)*up_massdetro(i, k), "updet"//cty, 'up_massdetr(kg/s/m^2)', '3d')
               call setGradsVar(jl, k, nvar, outt(i, k)*86400., "outt"//cty, 'outt K/day', '3d')
               call setGradsVar(jl, k, nvar, resten_t*86400., "rest"//cty, 'residuo T K/day', '3d')
               call setGradsVar(jl, k, nvar, resten_h*86400./real(c_cp), "resh"//cty, 'residuo H J/kg/day', '3d')
               call setGradsVar(jl, k, nvar, resten_q*86400.*real(c_alvl)/real(c_cp), "resq"//cty, 'residuo q K/day   ', '3d')
               call setGradsVar(jl, k, nvar, subten_t(i, k)*86400., "subt"//cty, 'subT K/day', '3d')
               call setGradsVar(jl, k, nvar, subten_h(i, k)*86400./real(c_cp), "subh"//cty, 'subH J/kg/day', '3d')
               call setGradsVar(jl, k, nvar, subten_q(i, k)*86400.*real(c_alvl)/real(c_cp), "subq"//cty, 'subq K/day   ', '3d')
               call setGradsVar(jl, k, nvar, outq(i, k)*86400.*real(c_alvl)/real(c_cp), "outq"//cty, 'outq K/s', '3d')
               call setGradsVar(jl, k, nvar, outqc(i, k)*86400.*real(c_alvl)/real(c_cp), "outqc"//cty, 'outqc K/day', '3d')
               call setGradsVar(jl, k, nvar, pre(i)*3600., "precip"//cty, 'precip mm', '2d')
               call setGradsVar(jl, k, nvar, prec_flx(i, k)*3600., "precflx"//cty, 'prec flx mm', '3d')
               call setGradsVar(jl, k, nvar, pwo(i, k), "pwo"//cty, ' xx', '3d')
               call setGradsVar(jl, k, nvar, outu(i, k)*86400., "outu"//cty, 'out_U m/s/day', '3d')
               call setGradsVar(jl, k, nvar, outv(i, k)*86400., "outv"//cty, 'out_V m/s/day', '3d')
               call setGradsVar(jl, k, nvar, xmb(i), "xmb"//cty, 'xmb kg/m2/s', '2d')
               call setGradsVar(jl, k, nvar, vvel2d(i, k), "W2d"//cty, 'W /m/s', '3d')
               call setGradsVar(jl, k, nvar, vvel1d(i), "W1d"//cty, 'W1s /m/s', '2d')
               call setGradsVar(jl, k, nvar, us(i, k), "us"//cty, 'U /m/s', '3d')
               call setGradsVar(jl, k, nvar, outu(i, k)*86400./(1.e-16 + xmb(i)), "delu"//cty, 'dellu', '3d')
               call setGradsVar(jl, k, nvar, evap_bcb(i, k)*1000., "evcb"//cty, 'g/kg', '3d')

               call setGradsVar(jl, k, nvar, tot_pw_up_chem(1, i), "pwup"//cty, 'pwup', '2d')
               call setGradsVar(jl, k, nvar, tot_pw_dn_chem(1, i), "pwdn"//cty, 'pwdn', '2d')
               !----
               !----
               call setGradsVar(jl, k, nvar, xmb(i)*dellah(i, k)*86400./real(c_cp), "delh"//cty, 'dellah K/day', '3d')
               call setGradsVar(jl, k, nvar, xmb(i)*dellaq(i, k)*86400.*real(c_alvl)/real(c_cp), "dellq"//cty, 'dellaq K/day', '3d')
               call setGradsVar(jl, k, nvar, xmb(i)*dellaqc(i, k)*86400.*real(c_alvl)/real(c_cp), "dellqc"//cty, 'dellaqc K/day' &
                              , '3d')
               call setGradsVar(jl, k, nvar, xmb(i), "xmb"//cty, 'm flux up (kg/s/m^2)', '2d')
               call setGradsVar(jl, k, nvar, aa1(i), "aa1"//cty, 'AA1 J/kg3)', '2d')
               call setGradsVar(jl, k, nvar, float(ierr(i)), "ierr"//cty, 'ierr #', '2d')
               call setGradsVar(jl, k, nvar, xmb(i)*dd_massentro(i, k), "ddent"//cty, 'dd_massentr(kg/s/m^2)', '3d')
               call setGradsVar(jl, k, nvar, xmb(i)*dd_massdetro(i, k), "dddet"//cty, 'dd_massdetr(kg/s/m^2)', '3d')
               !!     go to 333
               call setGradsVar(jl, k, nvar, hc(i, k), "hc"//cty, ' hc', '3d')
               call setGradsVar(jl, k, nvar, hco(i, k), "hco"//cty, ' hco', '3d')
               call setGradsVar(jl, k, nvar, dby(i, k), "dby"//cty, ' dbuo', '3d')
               !call set_grads_var(jl,k,nvar,QCUP(i,k),"qcup"//cty ,'C_UP','3d')
               call setGradsVar(jl, k, nvar, t_cup(i, k) - 273.15, "te"//cty, ' K', '3d')
               call setGradsVar(jl, k, nvar, 1000.*q_cup(i, k), "qe"//cty, ' kg kg-1', '3d')
               call setGradsVar(jl, k, nvar, he_cup(i, k), "he"//cty, ' he', '3d')
               call setGradsVar(jl, k, nvar, HKB(i), "hkb"//cty, ' H', '2d')
               call setGradsVar(jl, k, nvar, HKB(i), "hkb"//cty, ' H', '2d')
               call setGradsVar(jl, k, nvar, 1000.*zqexec(i), "qex"//cty, ' qex', '2d')
               call setGradsVar(jl, k, nvar, z_cup(i, max(1, k22(i))), "zs"//cty, ' m', '2d')
               call setGradsVar(jl, k, nvar, z_cup(i, max(1, kbcon(i))), "zbcon"//cty, ' m', '2d')
               call setGradsVar(jl, k, nvar, z_cup(i, max(1, ktop(i))), "ztop"//cty, ' m', '2d')
               call setGradsVar(jl, k, nvar, z_cup(i, max(1, klcl(i))), "zlcl"//cty, ' m', '2d')
               call setGradsVar(jl, k, nvar, z_cup(i, max(1, jmin(i))), "zjmin"//cty, ' m', '2d')
               call setGradsVar(jl, k, nvar, zws(i), "ws"//cty, ' m/s', '2d')
               call setGradsVar(jl, k, nvar, clfrac(i, k), "clfrac"//cty, 'shcf #', '3d')
               call setGradsVar(jl, k, nvar, entr_rate_2d(i, k), "entr"//cty, ' m-1', '3d')
               call setGradsVar(jl, k, nvar, cd(i, k), "detr"//cty, ' m-1', '3d')
               call setGradsVar(jl, k, nvar, pwdo(i, k), "pwd"//cty, ' xx', '3d')
               call setGradsVar(jl, k, nvar, edto(i), "edt"//cty, 'edt kg/m2/s', '2d')
               call setGradsVar(jl, k, nvar, e_dn, "EVAP"//cty, ' xx', '3d')
               call setGradsVar(jl, k, nvar, c_up, "CUP"//cty, ' xx', '3d')
               !       call set_grads_var(jl,k,nvar,trash,"TUP"//cty ,' xx','3d')
               !       call set_grads_var(jl,k,nvar,trash2,"TDN"//cty ,' xx','3d')
               call setGradsVar(jl, k, nvar, trash, "F1"//cty, ' F1', '3d')
               call setGradsVar(jl, k, nvar, trash2, "F2"//cty, ' F2', '3d')
               call setGradsVar(jl, k, nvar, p_liq_ice(i, k), "pli"//cty, '#', '3d')
               call setGradsVar(jl, k, nvar, melting_layer(i, k), "cpli"//cty, '#', '3d')
               call setGradsVar(jl, k, nvar, t(i, k), "t"//cty, 'temp K', '3d')
               call setGradsVar(jl, k, nvar, tn(i, k), "tn"//cty, 'temp K', '3d')
               call setGradsVar(jl, k, nvar, 1000.*q(i, k), "q"//cty, 'q g/kg', '3d')
               call setGradsVar(jl, k, nvar, 1000.*qo(i, k), "qn"//cty, 'q g/kg', '3d')
               call setGradsVar(jl, k, nvar, 1000.*qrco(i, k), "qrc"//cty, 'q g/kg', '3d')
               call setGradsVar(jl, k, nvar, 1000.*(q(i, k) + outq(i, k)*dtime), "qnc"//cty, 'q upd conv g/kg', '3d')
               call setGradsVar(jl, k, nvar, 1000.*(qo(i, k) + outq(i, k)*dtime), "qnall"//cty, 'q upd all g/kg', '3d')
               call setGradsVar(jl, k, nvar, 1000.*qrr(i, k), "qrr"//cty, 'qrr g/kg', '3d')
               call setGradsVar(jl, k, nvar, 1000.*qco(i, k), "qc"//cty, 'qc g/kg', '3d')
               call setGradsVar(jl, k, nvar, 1000.*qo_cup(i, k), "qcup"//cty, 'qc g/kg', '3d')
               call setGradsVar(jl, k, nvar, 1000.*qeso_cup(i, k), "qescup"//cty, 'qc g/kg', '3d')

               !~ call set_grads_var(jl,k,nvar,aa0(i),"a0"//cty,'aa0','2d')
               !~ call set_grads_var(jl,k,nvar,aa1_fa(i),"aa1fa"//cty,'aa1fa','2d')
               !~ call set_grads_var(jl,k,nvar,aa1_bl(i),"aa1bl"//cty,'aa1bl','2d')
               !~ call set_grads_var(jl,k,nvar,aa0_bl(i),"aa0bl"//cty,'aa0bl','2d')
               !~ call set_grads_var(jl,k,nvar,aa1(i),"a1"//cty,'aa1','2d')
               !~ call set_grads_var(jl,k,nvar,aa1(i)/(1.e-6+tau_ecmwf(i)),"mb13"//cty,'aa0','2d')
               !~ call set_grads_var(jl,k,nvar,xaa0(i),"xa0"//cty,'xaa0','2d')
               !~ call set_grads_var(jl,k,nvar,(XAA0(I)-AA1(I))/MBDT(I),"xk"//cty,'xk','2d')
333            continue
            end do
            if (wrtgrads .and. .not. p_use_gate) then
               call wrtBinCtl(1, kte, po(1, 1:kte), cumulus)
            end if
         end do
      end if
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

   end subroutine cupGf

   
   subroutine initCupGF(its, itf, ite, ktf, max_edt_ocean, max_edt_land, xland, zo, cumulus, random, depth_min, zkbmax, z_detr &
                     ,  zcutdown, kstabm, kbmax, ierr2 , ierr3, lambau_dp, lambau_dn, edtmin, edtmax, xland1, cap_max &
                     ,  aa0, aa1, aa2, aa3, aa1_bl, aa1_fa, aa0_bl, q_adv, aa1_radpbl, aa1_adv, alpha_adv, cin1, xk_x &
                     ,  edt , edto, tau_bl, q_wetbulb, t_wetbulb, tau_ecmwf, xf_dicycle, x_add_buoy, xf_coldpool &
                     ,  wlpool_bcon, cap_max_increment, mbdt, z, xz, hcdo, cupclw, qrcdo, hcot, c1d, xf_ens, pr_ens &
                     , evap_bcb, ierrc)
      !! ## Initialize variables of cupGF
      !!
      !! Author: Rodrigues, L. F. [LFR]
      !!
      !! E-mail: <mailto:luizfrodrigues@protonmail.com>
      !!
      !! Date: 10Fevereiro2023 14:09
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! INitializa variables of cupGF
      !!
      !! ** History**:
      !!
      !! --- 
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!
   
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'intCupGF' 
      !! subroutine name
   
      !Variables (input, output, inout)
      integer, intent(in) :: its, itf, ite, ktf

      character(len=*), intent(in) :: cumulus
      
      real, intent(in) :: max_edt_ocean
      real, intent(in) :: max_edt_land
      real, intent(in) :: xland(:)
      real, intent(in) :: zo(:, :)

      real, intent(out) :: random(:)
      real, intent(out) :: depth_min
      real, intent(out) :: zkbmax
      real, intent(out) :: z_detr
      real, intent(out) :: zcutdown

      integer, intent(out) :: kstabm        (:)
      integer, intent(out) :: kbmax         (:)
      integer, intent(out) :: ierr2         (:)
      integer, intent(out) :: ierr3         (:)
        
      real, intent(out) :: lambau_dp        (:)
      real, intent(out) :: lambau_dn        (:)
      real, intent(out) :: edtmin           (:)
      real, intent(out) :: edtmax           (:)
      real, intent(out) :: xland1           (:)
      real, intent(out) :: cap_max          (:)
      real, intent(out) :: aa0              (:)
      real, intent(out) :: aa1              (:)
      real, intent(out) :: aa2              (:)
      real, intent(out) :: aa3              (:)
      real, intent(out) :: aa1_bl           (:)
      real, intent(out) :: aa1_fa           (:)
      real, intent(out) :: aa0_bl           (:)
      real, intent(out) :: q_adv            (:)
      real, intent(out) :: aa1_radpbl       (:)
      real, intent(out) :: aa1_adv          (:)
      real, intent(out) :: alpha_adv        (:)
      real, intent(out) :: cin1             (:)
      real, intent(out) :: xk_x             (:)
      real, intent(out) :: edt              (:)
      real, intent(out) :: edto             (:)
      real, intent(out) :: tau_bl           (:)
      real, intent(out) :: q_wetbulb        (:)
      real, intent(out) :: t_wetbulb        (:)
      real, intent(out) :: tau_ecmwf        (:)
      real, intent(out) :: xf_dicycle       (:)
      real, intent(out) :: x_add_buoy       (:)
      real, intent(out) :: xf_coldpool      (:)
      real, intent(out) :: wlpool_bcon      (:)
      real, intent(out) :: cap_max_increment(:)
      real, intent(out) :: mbdt             (:)
      real, intent(out) :: z                (:,:)
      real, intent(out) :: xz               (:,:)
      real, intent(out) :: hcdo             (:,:)
      real, intent(out) :: cupclw           (:,:)
      real, intent(out) :: qrcdo            (:,:)
      real, intent(out) :: hcot             (:,:)
      real, intent(out) :: c1d              (:,:)
      real, intent(out) :: xf_ens           (:,:)
      real, intent(out) :: pr_ens           (:,:)
      real, intent(out) :: evap_bcb         (:,:)
   
      character(len=128), intent(out) :: ierrc(:)

      !Local variables:
      real :: cap_max_inc
      integer :: i_cnt
   
      random = 0.0

      if (trim(cumulus) == 'deep') then
         !--- maximum depth (mb) of capping inversion (larger cap = no convection)
         if (MOIST_TRIGGER == 0) then
            cap_max_inc = 20.
         else
            cap_max_inc = 90.
         endif
         !--- lambda_U parameter for momentum transport
         lambau_dp(:) = LAMBAU_DEEP
         lambau_dn(:) = LAMBAU_SHDN
         do i_cnt = its, itf
            if (xland(i_cnt) > 0.99) then !- over water
               edtmin(i_cnt) = 0.1
               edtmax(i_cnt) = max_edt_ocean  !
            else!- over land
               edtmin(i_cnt) = 0.1
               edtmax(i_cnt) = max_edt_land  !
            end if
         end do
         ! --- minimum depth (m), clouds must have
         depth_min = 1000.
         ! --- max height(m) above ground where updraft air can originate
         zkbmax = 4000.
         ! --- depth(m) over which downdraft detrains all its mass
         z_detr = 1000.
         ! ---  create a real random number in the interval [-use_random_num, +use_random_num]
         if (USE_RANDOM_NUM > 1.e-6) random = genRandom(its, ite, USE_RANDOM_NUM)
      elseif(trim(cumulus) == 'mid') then
         if (MOIST_TRIGGER == 0) then
            cap_max_inc = 10.
         else
            cap_max_inc = 90.
         endif
         lambau_dp(:) = LAMBAU_SHDN
         lambau_dn(:) = LAMBAU_SHDN
         do i_cnt = its, itf
            if (xland(i_cnt) > 0.99) then !- over water
               edtmin(i_cnt) = 0.1
               edtmax(i_cnt) = max_edt_ocean  !
            else!- over land
               edtmin(i_cnt) = 0.1
               edtmax(i_cnt) = max_edt_land  !
            end if
         end do
         if (C0_MID < 1.e-8) edtmin(:) = 0.0
         depth_min = 500.
         zkbmax = 3000.
         z_detr = 300.
      elseif(trim(cumulus) == 'shallow') then
         if (MOIST_TRIGGER == 0) then
            cap_max_inc = 25.
         else
            cap_max_inc = 10.
         endif
         lambau_dp(:) = LAMBAU_SHDN
         lambau_dn(:) = LAMBAU_SHDN
         ! --- max/min allowed value for epsilon (ratio downdraft base mass flux/updraft
         !     base mass flux
         ! --  note : to make the evaporation stronger => increase "edtmin"
         edtmin(:) = 0.0
         edtmax(:) = 0.0
         depth_min = 500.
         zkbmax = 3000.
         z_detr = 300.
      else 
         ! Aqui precisa entrar uma mnesgame de erro - escolha inválida
      end if

      if (p_pgcon .ne. 0.) then
         lambau_dp(:) = 0.
         lambau_dn(:) = 0.
      end if

      do i_cnt = its, itf
         kbmax(i_cnt)       = 1
         kstabm(i_cnt)      = ktf - 1
         ierr2(i_cnt)       = 0
         ierr3(i_cnt)       = 0
         xland1(i_cnt)      = xland(i_cnt) ! 1.
         cap_max(i_cnt)     = CAP_MAXS
         ierrc(i_cnt)       = "ierrtxt"
         aa0(i_cnt)         = 0.0
         aa1(i_cnt)         = 0.0
         aa2(i_cnt)         = 0.0
         aa3(i_cnt)         = 0.0
         aa1_bl(i_cnt)      = 0.0
         aa1_fa(i_cnt)      = 0.0
         aa0_bl(i_cnt)      = 0.0
         q_adv(i_cnt)       = 0.0
         aa1_radpbl(i_cnt)  = 0.0
         aa1_adv(i_cnt)     = 0.0
         alpha_adv(i_cnt)   = 0.0
         cin1(i_cnt)        = 0.0
         xk_x(i_cnt)        = 0.0
         edt(i_cnt)         = 0.0
         edto(i_cnt)        = 0.0
         tau_bl(i_cnt)      = 0.0
         q_wetbulb(i_cnt)   = 0.0
         t_wetbulb(i_cnt)   = 0.0
         tau_ecmwf(i_cnt)   = 0.0
         xf_dicycle(i_cnt)  = 0.0
         x_add_buoy(i_cnt)  = 0.0
         xf_coldpool(i_cnt) = 0.0
         wlpool_bcon(i_cnt) = 0.0
         z(i_cnt, :)        = zo(i_cnt, :)
         xz(i_cnt, :)       = zo(i_cnt, :)
         hcdo(i_cnt, :)     = 0.0
         cupclw(i_cnt, :)   = 0.0
         qrcdo(i_cnt, :)    = 0.0
         hcot(i_cnt, :)     = 0.0
         c1d(i_cnt, :)      = 0.0
         xf_ens(i_cnt, :)   = 0.0
         pr_ens(i_cnt, :)   = 0.0
         evap_bcb(i_cnt, :) = 0.0
         cap_max_increment(i_cnt) = cap_max_inc
      end do

      !--- height(m) above which no downdrafts are allowed to originate
      zcutdown = 3000.

      !--- mbdt ~ xmb * timescale
      do i_cnt = its, itf
         mbdt(i_cnt) = 0.1!*dtime*xmb_nm1(i)
         !mbdt(i)= 100.*(p_cup(i,kbcon(i))-p_cup(i,kbcon(i)+1))/(g*dtime)
         !mbdt(i)= 0.1*mbdt(i)
      end do
   
   end subroutine initCupGF


   ! ---------------------------------------------------------------------------------------------------
   function genRandom(its, ite, use_random_num) result(random)
      !! ## Generate a random array
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! —
      !! **Full description**:
      !!
      !! Generate a random array
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !!
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'genRandom' ! Subroutine Name
   
      !Variables (input, output, inout)
      integer, intent(in)  :: its
      integer, intent(in)  :: ite
      real, intent(in) :: use_random_num
      
   
      !Local variables:
      real :: random(its:ite)
      integer   :: i
      integer(8) :: iran, ranseed = 0   
      
      !Code:
      call system_clock(ranseed)
      ranseed = mod(ranseed, 2147483646) + 1 !seed between 1 and 2^31-2
      iran = -ranseed

      !-- Ran1 produces numbers between [ 0,1]
      !-- random        will be between [-1,1]
      !-- with use_random_num the interval will be [-use_random_num,+use_random_num]
      do i = its, ite
         random(i) = use_random_num*2.0*(0.5 - real(Ran1(IRAN), 4))
         !print*,"ran=",i,random(i)
      end do

      if (maxval(abs(random)) > use_random_num) stop "random > use_random_num"
   
   end function genRandom  

   ! ------------------------------------------------------------------------------------
   subroutine cupDdEdt(cumulus, ierr, us, vs, z, ktop, kbcon, edt, p, pwav, pw, ccn, pwev, edtmax, edtmin, maxens2, edtc, psum2 &
                     , psumh, rho, aeroevap, itf, ktf, ipr, jpr, its, ite, kts, kte, vshear)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'cupDdEdt' ! Subroutine name

      real, parameter :: p_alpha3 = 1.9, p_beta3 = -1.13
   
      !Variables (input, output, inout)
      integer, intent(in) :: ipr, jpr, aeroevap, itf, ktf, its, ite, kts, kte, maxens2

      integer, intent(in) :: ktop(:)
      integer, intent(in) :: kbcon(:)

      real, intent(in) :: rho(:, :)
      real, intent(in) :: us(:, :)
      real, intent(in) :: vs(:, :)
      real, intent(in) :: z(:, :)
      real, intent(in) :: p(:, :)
      real, intent(in) :: pw(:, :)
      real, intent(in) :: pwav(:)
      real, intent(in) :: pwev(:)
      real, intent(in) :: ccn(:)
      real, intent(in) :: psum2(:)
      real, intent(in) :: psumh(:)
      real, intent(in) :: edtmax(:)
      real, intent(in) :: edtmin(:)

      character(len=*), intent(in)  :: cumulus

      integer, intent(inout) :: ierr(:)

      real, intent(out) :: edtc(:, :)
      real, intent(out) :: edt(:)
      real, intent(out) :: vshear(:)
      
      !Local variables:
      integer :: i, k, kk
      real :: einc, pef, pefb, prezk, zkbc
      real, dimension(its:ite) :: vws, sdp
      real :: pefc, aeroadd, rhoc, dp, prop_c

      ! determine downdraft strength in terms of windshear
      ! calculate an average wind shear over the depth of the cloud
      edt = 0.
      vws = 0.
      sdp = 0.
      vshear = 0.
      edtc = 0.

      if (trim(cumulus) == 'shallow') return

      do i = its, itf
         if (ierr(i) /= 0) cycle
         do kk = kbcon(i), ktop(i)
            dp = p(i, kk) - p(i, kk + 1)
            vws(i) = vws(i) + (abs((us(i, kk + 1) - us(i, kk))/(z(i, kk + 1) - z(i, kk))) + abs((vs(i, kk + 1) - vs(i, kk)) &
                   / (z(i, kk + 1) - z(i, kk))))*dp
            sdp(i) = sdp(i) + dp
         end do
         vshear(i) = 1.e3*vws(i)/sdp(i)
      end do

      do i = its, itf
         if (ierr(i) /= 0) cycle
         pef = (1.591 - 0.639*vshear(i) + 0.0953*(vshear(i)**2) - 0.00496*(vshear(i)**3))

         !print*,"shear=",vshear(i),pef,1-max(min(pef,0.9),0.1)
         pef = min(pef, 0.9)
         pef = max(pef, 0.1)
         edt(i) = 1.-pef

         !--- cloud base precip efficiency
         if (USE_REBCB == 0) then
            zkbc = z(i, kbcon(i))*3.281e-3
            prezk = 0.02
            if (zkbc > 3.0) prezk = 0.96729352 + zkbc*(-0.70034167 + zkbc * (0.162179896 + zkbc*(-1.2569798e-2 + zkbc &
                                  * (4.2772e-4 - zkbc*5.44e-6))))
            if (zkbc > 25.) prezk = 2.4
            pefb = 1./(1.+prezk)
            pefb = min(pefb, 0.9)
            pefb = max(pefb, 0.1)
            edt(i) = 1.-0.5*(pefb + pef)
         end if

         if (aeroevap .gt. 1) then
            aeroadd = (c_ccnclean**p_beta3)*((psumh(i))**(p_alpha3 - 1)) !*1.e6
            !if(i.eq.ipr)write(0,*)'edt',ccnclean,psumh(i),aeroadd
            !prop_c=.9/aeroadd
            prop_c = .5*(pefb + pef)/aeroadd
            aeroadd = (ccn(i)**p_beta3)*((psum2(i))**(p_alpha3 - 1)) !*1.e6
            !if(i.eq.ipr)write(0,*)'edt',ccn(i),psum2(i),aeroadd,prop_c
            aeroadd = prop_c*aeroadd
            pefc = aeroadd
            if (pefc .gt. 0.9) pefc = 0.9
            if (pefc .lt. 0.1) pefc = 0.1
            EDT(I) = 1.-pefc
            if (aeroevap .eq. 2) EDT(I) = 1.-.25*(pefb + pef + 2.*pefc)
         end if

      end do
      do i = its, itf
         if (ierr(i) /= 0) cycle
         edtc(i, 1) = -edt(i)*pwav(i)/pwev(i)
         edtc(i, 1) = min(edtmax(i), edtc(i, 1))
         edtc(i, 1) = max(edtmin(i), edtc(i, 1))
      end do

   end subroutine cupDdEdt

   !------------------------------------------------------------------------------------
   subroutine cupDdMoisture(cumulus, ierrc, zd, hcd, hes_cup, qcd, qes_cup, pwd, q_cup, z_cup, dd_massentr &
                          , dd_massdetr, jmin, ierr, gamma_cup, pwev, bu, qrcd, q_env, he, t_cup, iloop, t_wetbulb &
                          , q_wetbulb, qco, pwavo, itf, ktf, its, ite, kts, kte)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'cupDdMoisture' ! Subroutine Name
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte   

      integer, intent(in) :: iloop
      integer, intent(in) :: jmin(:)

      ! mentr_rate = entrainment rate
      ! qrch = saturation q in cloud
      ! pwev = total normalized integrated evaoprate (I2)
      ! entr = entrainment rate
      ! cdd  = detrainment function
      !
      real, intent(in) :: t_wetbulb(:)
      real, intent(in) :: q_wetbulb(:)
      real, intent(in) :: pwavo(:)
      real, intent(in) :: zd(:, :)
      !! normalized downdraft mass flux
      real, intent(in) :: t_cup(:, :)
      real, intent(in) :: hes_cup(:, :)
      !! saturation h on model cloud levels
      real, intent(in) :: hcd(:, :)
      !! h in model cloud
      real, intent(in) :: qes_cup(:, :)
      !! saturation q on model cloud levels
      real, intent(in) :: q_cup(:, :)
      !! environmental q on model cloud levels
      real, intent(in) :: z_cup(:, :)
      real, intent(in) :: dd_massentr(:, :)
      real, intent(in) :: dd_massdetr(:, :)
      real, intent(in) :: gamma_cup(:, :)
      !! gamma on model cloud levels
      real, intent(in) :: q_env(:, :)
      !! environmental q on model levels
      real, intent(in) :: he(:, :)
      real, intent(in) :: qco(:, :)

      character(len=*), intent(in) :: cumulus

      integer, intent(inout) :: ierr(:)

      real, intent(out) :: qcd(:, :)
      !! cloud q (including liquid water) after entrainment
      !! in-downdradt water vapor mixing ratio
      real, intent(out) :: qrcd(:, :)
      !! saturation water vapor mixing ratio
      real, intent(out) :: pwd(:, :)
      !! evaporate at that level
      real, intent(out) :: pwev(:)
      !! column integrated rain evaporation (normalized)
      real, intent(out) :: bu(:)
      !! buoancy term

      !Local variables:
      character*128 :: ierrc(its:ite)
      integer :: i, k
      real :: dh, dz, dq_eva, denom, fix_evap
      !
      bu = 0.  
      qcd = 0. 
      qrcd = 0.
      pwev = 0.
      pwd = 0. 

      if (trim(cumulus) == 'shallow') return
      !
      do i = its, itf
         if (ierr(i) /= 0) cycle

         !-- boundary condition in jmin ('level of free sinking')
         k = jmin(i)
         dz = z_cup(i, k + 1) - z_cup(i, k)
         qcd(i, k) = q_cup(i, k)
         if (USE_WETBULB == 1) then
            !--option 1
            !qcd(i,k)=q_wetbulb(i)
            !--option 2
            qcd(i, k) = 0.5*(q_wetbulb(i) + qco(i, k)) ! mixture 50% env air + updraft
         end if
         dh = hcd(i, k) - hes_cup(i, k)
         if (dh .lt. 0) then
            qrcd(i, k) = (qes_cup(i, k) + (1./real(c_alvl))*(gamma_cup(i, k)/(1.+gamma_cup(i, k)))*dh)
         else
            qrcd(i, k) = qes_cup(i, k)
         end if
         pwd(i, k) = zd(i, k)*min(0., qcd(i, k) - qrcd(i, k))
         qcd(i, k) = qrcd(i, k)
         pwev(i) = pwev(i) + pwd(i, k)
         bu(i) = dz*dh
         do k = jmin(i) - 1, kts, -1
            dz = z_cup(i, k + 1) - z_cup(i, k)
            !-- downward transport + mixing
            denom = (zd(i, k + 1) - 0.5*dd_massdetr(i, k) + dd_massentr(i, k))
            if (denom == 0.0) then
               qcd(i, k) = qcd(i, k + 1)
            else
               qcd(i, k) = (qcd(i, k + 1)*zd(i, k + 1) - 0.5*dd_massdetr(i, k)*qcd(i, k + 1) + dd_massentr(i, k)*q_env(i, k))/denom
            end if

            !--- to be negatively buoyant, hcd should be smaller than hes!
            !--- ideally, dh should be negative till dd hits ground, but that is not always
            !--- the case
            dh = hcd(i, k) - hes_cup(i, k)
            bu(i) = bu(i) + dz*dh
            qrcd(i, k) = qes_cup(i, k) + (1./real(c_alvl))*(gamma_cup(i, k)/(1.+gamma_cup(i, k)))*dh

            !-- rain water evaporation amount at layer k
            dq_eva = qcd(i, k) - qrcd(i, k)

            if (dq_eva .gt. 0.) then
               dq_eva = 0.
               qrcd(i, k) = qcd(i, k)
            end if
            !-- amount of the evaporated rain water
            pwd(i, k) = zd(i, k)*dq_eva  ! kg[water vapor]/kg[air]

            !-- source term for in-downdraft water vapor mixing ratio
            qcd(i, k) = qrcd(i, k)     ! => equiv to qcd = qcd - dq_eva !( -dq_eva >0 => source term for qcd)

            !-- total evaporated rain water
            pwev(i) = pwev(i) + pwd(i, k)

            !-- for GEOS diagnostic
            ! evap(i,k) = - edt * xmb * zd * dq_eva = - edt * xmb * pwd (i,k)
            ! downdfrat temp = (hcd(i,k)-qcd(i,k)*xlv-g*z_cup(i,k))/cp - 273.15

         end do

         if (pwev(i) .ge. 0 .and. iloop .eq. 1) then
            ierr(i) = 70
            ierrc(i) = "problem with buoy in cup_dd_moisture"
         end if
         if (bu(i) .ge. 0 .and. iloop .eq. 1) then
            ierr(i) = 73
            ierrc(i) = "problem2 with buoy in cup_dd_moisture"
         end if

         !-- fix evap, in case of not conservation
         if (abs(pwev(i)) > pwavo(i) .and. ierr(i) == 0) then
            fix_evap = pwavo(i)/(1.e-16 + abs(pwev(i)))
            pwev(i) = 0.

            do k = jmin(i), kts, -1
               pwd(i, k) = pwd(i, k)*fix_evap
               pwev(i) = pwev(i) + pwd(i, k)
               dq_eva = pwd(i, k)/(1.e-16 + zd(i, k))
               qcd(i, k) = qrcd(i, k) + dq_eva
            end do
            if (pwev(i) .ge. 0.) then
               ierr(i) = 70
               ierrc(i) = "problem with buoy in cup_dd_moisture"
            end if
         end if
      end do!--- end loop over i

   end subroutine cupDdMoisture

   ! ------------------------------------------------------------------------------------
   subroutine cupEnv(z_heights, qes, he, hes, temp_env, mixratio_env, press_env, z1, psur, ierr, itest, itf, ktf &
                   , its, ite, kts, kte)
      !! ## Determine the thermodynamical variables at full levels
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Determine the thermodynamical variables at full levels:
      !! HE   moist static energy (J/kg)
      !! HES  saturated moist static energy (J/kg)
      !! QES  saturated water vapor mixing ratio (kg/kg)
      !! TV   virtual temperature (K)
      !!
      !! ** History**:
      !! SATUR_CALC = 0 applies the original formulation from Grell (1993)
      !! SATUR_CALC = 1 applies the optional (now default) formulation from Tiedtke et al. (1988)
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'cupEnv' ! Subroutine Name
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte, itest

      integer, intent(in) :: ierr(:)

      real, intent(in) :: press_env(:, :)
      !! environmental pressure
      real, intent(in) :: temp_env(:, :)
      !! environmental temp
      real, intent(in) :: mixratio_env(:, :)
      !! environmental mixing ratio
      real, intent(in) :: psur(:)
      !! surface pressure
      real, intent(in) :: z1(:)
      !! terrain elevation
      real, intent(inout)  :: z_heights(:, :)
      !! environmental heights
      real, intent(out) :: he(:, :)
      !! environmental moist static energy
      real, intent(out) :: hes(:, :)
      !! environmental saturation moist static energy
      real, intent(out) :: qes(:, :)
      !! environmental saturation mixing ratio

      !Local variables
      integer :: i, k, iph
      real, dimension(its:ite, kts:kte) :: tv
      !! environmental virtual temp
      real :: e, tvbar, pqsat
   
      !Code:
      he = 0.0
      hes = 0.0
      qes = 0.0

      if (SATUR_CALC == 0) then
         do k = kts, ktf
            do i = its, itf
               if (ierr(i) /= 0) cycle

                  e = SatVap(temp_env(i, k))
                  qes(i, k) = max(1.e-8, 0.622*e/max(1.e-8, (press_env(i, k) - e)))
                  qes(i, k) = min(qes(i,k),c_max_qsat)
                  qes(i, k) = max(qes(i,k),mixratio_env(i, k))
                  tv(i, k) = temp_env(i, k) + .608*mixratio_env(i, k)*temp_env(i, k)
            end do
         end do
      else
         !--- better formulation for the mixed phase regime
         do k = kts, ktf
            do i = its, itf
                  if (ierr(i) /= 0) cycle
                  qes(i, k) = SaturSpecHum(temp_env(i, k), press_env(i, k))
                  qes(i, k) = min(c_max_qsat, max(1.e-08, qes(i, k)))
                  qes(i, k) = max(qes(i, k), mixratio_env(i, k))
                  tv (i, k) = temp_env(i, k) + .608*mixratio_env(i, k)*temp_env(i, k)
            end do
         end do
      end if

      ! !--- z's are calculated with changed h's and q's and t's
      ! !--- if itest=2
      ! if (itest .eq. 1 .or. itest .eq. 0) then
      !    do i = its, itf
      !       if (ierr(i) .eq. 0) then
      !          z_heights(i, 1) = max(0., z1(i)) - (Alog(press_env(i, 1)) - Alog(psur(i)))*287.*tv(i, 1)/c_grav
      !       end if
      !    end do
      !    ! --- calculate heights
      !    do k = kts + 1, ktf
      !       do i = its, itf
      !          if (ierr(i) .eq. 0) then
      !             tvbar = .5*tv(i, k) + .5*tv(i, k - 1)
      !             z_heights(i, k) = z_heights(i, k - 1) - (Alog(press_env(i, k)) - Alog(press_env(i, k - 1)))*287.*tvbar/c_grav
      !          end if
      !       end do
      !    end do
      ! else if (itest .eq. 2) then
      !    do k = kts, ktf
      !       do i = its, itf
      !          if (ierr(i) .eq. 0) then
      !             z_heights(i, k) = (he(i, k) - 1004.*temp_env(i, k) - 2.5e6*mixratio_env(i, k))/c_grav
      !             z_heights(i, k) = max(1.e-3, z_heights(i, k))
      !          end if
      !       end do
      !    end do
      ! else if (itest .eq. -1) then
      ! end if

      !--- calculate moist static energy - HE
      !    saturated moist static energy - HES
      do k = kts, ktf
         do i = its, itf
            if (ierr(i) /= 0) cycle

            he (i, k) = c_grav*z_heights(i, k) + real(c_cp)*temp_env(i, k) + real(c_alvl)*mixratio_env(i, k)
            hes(i, k) = c_grav*z_heights(i, k) + real(c_cp)*temp_env(i, k) + real(c_alvl)*qes(i, k)
            he (i, k) = min(hes(i, k), he(i, k))

         end do
      end do

   end subroutine cupEnv

   ! ------------------------------------------------------------------------------------
   subroutine cupEnvCLev(temp_env, qes, q, he, hes, z_heights, pres_env, qes_cup, q_cup, he_cup, us, vs, u_cup, v_cup &
                      ,  hes_cup, z_cup, p_cup, gamma_cup, t_cup, psur, tsur, ierr, z1, itf, ktf, its, ite, kts, kte)
      !! ## Determine the thermodynamical variables at half levels
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Determine the thermodynamical variables at half levels. 
      !! The name becomes XXXX_cup
      !!
      !! ** History**:
      !! CLEV_GRID == 0 : optional formulation using weighted mean 
      !! CLEV_GRID == 1 : applies the optional (now default) formulation from Tiedtke et al. (1988)
      !! CLEV_GRID == 2 : original formulation from Grell 1993
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'cupEnvCLev' ! Subroutine Name
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte

      integer, intent(in) :: ierr(:)

      real, intent(in) :: qes(:, :)
      !! environmental saturation mixing ratio
      real, intent(in) :: q(:, :)
      !! environmental mixing ratio
      real, intent(in) :: he(:, :)
      !! environmental moist static energy
      real, intent(in) :: hes(:, :)
      !! environmental saturation moist static energy
      real, intent(in) :: z_heights(:, :)
      !! environmental heights
      real, intent(in) :: pres_env(:, :)
      !! environmental pressure
      real, intent(in) :: temp_env(:, :)
      !! environmental temp
      real, intent(in) :: us(:, :)
      !!
      real, intent(in) :: vs(:, :)
      !!
      real, intent(in) :: psur(:)
      !! surface pressure
      real, intent(in) :: z1(:)
      !! terrain elevation
      real, intent(in) :: tsur(:)
      !! surface temperature

      real, intent(out) :: qes_cup(:, :)
      !! environmental saturation mixing ratio on cloud levels
      real, intent(out) :: q_cup(:, :)
      !! environmental mixing ratio on cloud levels
      real, intent(out) :: he_cup(:, :)
      !! environmental moist static energy on cloud levels
      real, intent(out) :: hes_cup(:, :)
      !! environmental saturation moist static energy on cloud levels
      real, intent(out) :: z_cup(:, :)
      !! environmental heights on cloud levels
      real, intent(out) :: p_cup(:, :)
      !! environmental pressure on cloud levels
      real, intent(out) :: gamma_cup(:, :)
      !! gamma on cloud levels
      real, intent(out) :: t_cup(:, :)
      !! environmental temp on cloud levels
      real, intent(out) :: u_cup(:, :)
      real, intent(out) :: v_cup(:, :)

      !Local variables:
      integer :: i, k
      real  :: p1, p2, ct1, ct2, rho
      integer :: irun = 0

      qes_cup = 0.
      q_cup = 0.
      hes_cup = 0.
      he_cup = 0.
      z_cup = 0.
      p_cup = 0.
      t_cup = 0.
      gamma_cup = 0.
      u_cup = 0.
      v_cup = 0.

      if (CLEV_GRID == 2) then
         !--original formulation
         do k = kts + 1, ktf
            do i = its, itf
               if (ierr(i) /= 0) cycle
               qes_cup(i, k) = .5*(qes(i, k - 1) + qes(i, k))
               q_cup(i, k) = .5*(q(i, k - 1) + q(i, k))
               hes_cup(i, k) = .5*(hes(i, k - 1) + hes(i, k))
               he_cup(i, k) = .5*(he(i, k - 1) + he(i, k))
               if (he_cup(i, k) .gt. hes_cup(i, k)) he_cup(i, k) = hes_cup(i, k)
               z_cup(i, k) = .5*(z_heights(i, k - 1) + z_heights(i, k))
               p_cup(i, k) = .5*(pres_env(i, k - 1) + pres_env(i, k))
               t_cup(i, k) = .5*(temp_env(i, k - 1) + temp_env(i, k))
               gamma_cup(i, k) = (real(c_alvl)/real(c_cp))*(real(c_alvl)/(c_rm*t_cup(i, k) * t_cup(i, k)))*qes_cup(i, k)
               u_cup(i, k) = .5*(us(i, k - 1) + us(i, k))
               v_cup(i, k) = .5*(vs(i, k - 1) + vs(i, k))

            end do
         end do
         do i = its, itf
            if (ierr(i) /= 0) cycle
            qes_cup(i, 1) = qes(i, 1)
            q_cup(i, 1) = q(i, 1)
            !hes_cup(i,1)=hes(i,1)
            !he_cup(i,1)=he(i,1)
            hes_cup(i, 1) = c_grav*z1(i) + real(c_cp)*temp_env(i, 1) + real(c_alvl)*qes(i, 1)
            he_cup(i, 1) = c_grav*z1(i) + real(c_cp)*temp_env(i, 1) + real(c_alvl)*q(i, 1)
            !z_cup(i,1)=.5*(z(i,1)+z1(i))
            !p_cup(i,1)=.5*(p(i,1)+psur(i))
            z_cup(i, 1) = z1(i)
            p_cup(i, 1) = psur(i)
            t_cup(i, 1) = temp_env(i, 1)
            gamma_cup(i, 1) = real(c_alvl)/real(c_cp)*(real(c_alvl)/(c_rm*t_cup(i, 1) &
                                           *t_cup(i, 1)))*qes_cup(i, 1)
            u_cup(i, 1) = us(i, 1)
            v_cup(i, 1) = vs(i, 1)
         end do
         !do k=kts,ktf
         ! i=1
         !        print*,"air_dens=",k,z_cup(i,k),p_cup(i,k),(p_cup(i,k)-p_cup(i,k+1))/(z_cup(i,k+1)-z_cup(i,k))/g
         !enddo
      elseif (CLEV_GRID == 0) then
         !--- weigthed mean
         do i = its, itf
            if (ierr(i) /= 0) cycle
            p_cup(i, 1) = psur(i)
            z_cup(i, 1) = z1(i)
            do k = kts, ktf - 1
               p_cup(i, k + 1) = 2.0*pres_env(i, k) - p_cup(i, k)
               z_cup(i, k + 1) = 2.0*z_heights(i, k) - z_cup(i, k)
            end do
            ! ----------- p,T          k+1
            !p1
            ! ----------- p_cup,T_cup  k+1
            !p2
            ! ----------- p,T          k
            !
            ! ----------- p_cup,T_cup  k
            do k = kts, ktf - 1
               p1 = abs((pres_env(i, k + 1) - p_cup(i, k + 1))/(pres_env(i, k + 1) - pres_env(i, k)))
               p2 = abs((p_cup(i, k + 1) - pres_env(i, k))/(pres_env(i, k + 1) - pres_env(i, k)))

               t_cup(i, k + 1) = p1*temp_env(i, k) + p2*temp_env(i, k + 1)

               u_cup(i, k + 1) = p1*us(i, k) + p2*us(i, k + 1)
               v_cup(i, k + 1) = p1*vs(i, k) + p2*vs(i, k + 1)
               q_cup(i, k + 1) = p1*q(i, k) + p2*q(i, k + 1)
               he_cup(i, k + 1) = p1*he(i, k) + p2*he(i, k + 1)

               qes_cup(i, k + 1) = p1*qes(i, k) + p2*qes(i, k + 1)
               hes_cup(i, k + 1) = p1*hes(i, k) + p2*hes(i, k + 1)

               if (he_cup(i, k + 1) .gt. hes_cup(i, k + 1)) he_cup(i, k + 1) = hes_cup(i, k + 1)

               gamma_cup(i, k + 1) = (real(c_alvl)/real(c_cp))*(real(c_alvl)/(c_rm*t_cup(i, k + 1) * t_cup(i, k + 1))) &
                                   * qes_cup(i, k + 1)
            end do
            !--- surface level from X(kts) and X_cup(kts+1) determine X_cup(kts)
            k = kts
            p1 = abs(pres_env(i, k) - p_cup(i, k))
            p2 = abs(p_cup(i, k + 1) - p_cup(i, k))

            ct1 = (p1 + p2)/p2
            ct2 = p1/p2

            t_cup(i, k) = ct1*temp_env(i, k) - ct2*t_cup(i, k + 1)
            q_cup(i, k) = ct1*q(i, k) - ct2*q_cup(i, k + 1)

            u_cup(i, k) = ct1*us(i, k) - ct2*u_cup(i, k + 1)
            v_cup(i, k) = ct1*vs(i, k) - ct2*v_cup(i, k + 1)
            qes_cup(i, k) = ct1*qes(i, k) - ct2*qes_cup(i, k + 1)

            hes_cup(i, k) = c_grav*z_cup(i, k) + real(c_cp)*t_cup(i, k) + real(c_alvl)*qes_cup(i, k)
            he_cup(i, k) = c_grav*z_cup(i, k) + real(c_cp)*t_cup(i, k) + real(c_alvl)*q_cup(i, k)

            if (he_cup(i, k) .gt. hes_cup(i, k)) he_cup(i, k) = hes_cup(i, k)

            gamma_cup(i, k) = real(c_alvl)/real(c_cp)*(real(c_alvl)/(c_rm*t_cup(i, k)*t_cup(i, k)))*qes_cup(i, k)
         end do
      elseif (CLEV_GRID == 1) then
         !--- based on Tiedke (1989)
         do i = its, itf
            if (ierr(i) /= 0) cycle
            do k = ktf, kts + 1, -1

               qes_cup(i, k) = qes(i, k)
               q_cup(i, k) = q(i, k)
               p_cup(i, k) = 0.5*(pres_env(i, k - 1) + pres_env(i, k))
               z_cup(i, k) = 0.5*(z_heights(i, k - 1) + z_heights(i, k))
               t_cup(i, k) = (max(real(c_cp)*temp_env(i, k - 1) + c_grav*z_heights(i, k - 1), real(c_cp)*temp_env(i, k) + c_grav &
                           * z_heights(i, k)) - c_grav*z_cup(i, k))/real(c_cp)

               if (qes(i, k) < c_max_qsat) &
                  call getInterp(qes_cup(i, k), t_cup(i, k), p_cup(i, k), qes_cup(i, k), t_cup(i, k))

               q_cup(i, k) = min(q(i, k), qes(i, k)) + qes_cup(i, k) - qes(i, k)
               q_cup(i, k) = max(q_cup(i, k), 0.0)
            end do
            !---level kts
            qes_cup(i, 1) = qes(i, 1)
            q_cup(i, 1) = q(i, 1)
            z_cup(i, 1) = z1(i)
            p_cup(i, 1) = psur(i)

            t_cup(i, 1) = (real(c_cp)*temp_env(i, 1) + c_grav*z_heights(i, 1) - c_grav*z_cup(i, 1))/real(c_cp)

            hes_cup(i, 1) = c_grav*z_cup(i, 1) + real(c_cp)*t_cup(i, 1) + real(c_alvl)*qes_cup(i, 1)
            he_cup(i, 1) = c_grav*z_cup(i, 1) + real(c_cp)*t_cup(i, 1) + real(c_alvl)*q_cup(i, 1)

            gamma_cup(i, 1) = real(c_alvl)/real(c_cp)*(real(c_alvl)/(c_rm*t_cup(i, 1)*t_cup(i, 1)))*qes_cup(i, 1)
            u_cup(i, 1) = us(i, 1)
            v_cup(i, 1) = vs(i, 1)

            do k = ktf, kts + 1, -1
               p1 = max(real(c_cp)*t_cup(i, k) + c_grav*z_cup(i, k), real(c_cp)*t_cup(i, k - 1) + c_grav*z_cup(i, k - 1))
               t_cup(i, k) = (p1 - c_grav*z_cup(i, k))/real(c_cp)

               hes_cup(i, k) = real(c_cp)*t_cup(i, k) + real(c_alvl)*qes_cup(i, k) + c_grav*z_cup(i, k)
               he_cup(i, k) = real(c_cp)*t_cup(i, k) + real(c_alvl)*q_cup(i, k) + c_grav*z_cup(i, k)
               if (he_cup(i, k) .gt. hes_cup(i, k)) he_cup(i, k) = hes_cup(i, k)

               gamma_cup(i, k) = (real(c_alvl)/real(c_cp))*(real(c_alvl)/(c_rm*t_cup(i, k)*t_cup(i, k)))*qes_cup(i, k)
               u_cup(i, k) = us(i, k)
               v_cup(i, k) = vs(i, k)
            end do
         end do
      else
         stop "cup_env_clev"
      end if

      return
      !IF( MAPL_AM_I_ROOT() .and. irun == 0) then
      irun = 1
      do i = its, itf
         if (ierr(i) == 0) then
            do k = kts, kte - 1
               rho = 100*(p_cup(i, k) - p_cup(i, k + 1))/(z_cup(i, k + 1) - z_cup(i, k))/c_grav ! air dens by hidrostatic balance (kg/m3)
               write (23, 101) i, k, z_cup(i, k), p_cup(i, k), t_cup(i, k), q_cup(i, k)*1000., he_cup(i, k), u_cup(i, k) &
                             , v_cup(i, k), rho

               rho = 100*(pres_env(i, k) - pres_env(i, k + 1))/(z_heights(i, k + 1) - z_heights(i, k))/c_grav
               write (25, 101) i, k, z_heights(i, k), pres_env(i, k), temp_env(i, k), q(i, k)*1000., he(i, k), us(i, k), vs(i, k) &
                             , rho

101            format(2i3, 8f15.5)
            end do
            exit
            !            goto 400
         end if
      end do

   end subroutine cupEnvCLev

   ! -------------------------------------------------------------------------------------------------------------
   subroutine cupForcingEns3dMid(aa0, aa1, xaa0, mbdt, dtime, ierr, po_cup, ktop, k22, kbcon, kpbl, ichoice &
                               , maxens, maxens3, itf, ktf, its, ite, kts, kte, tau_ecmwf, aa1_bl, xf_dicycle &
                               , dhdt, xff_mid, zws, hc, hco, he_cup, heo_cup, wlpool, xf_coldpool)
      !! ## Determine the mass flux at cloud base for the congestus mode using a set of
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Determine the mass flux at cloud base for the congestus mode using a set of
      !! closures
      !! a) W* closure (Grant,2001)
      !! b) Boundary layer quasi-equilibrium (Raymond 1995)
      !! c) stability removal with a timescale (tau_ecmwf,  Bechtold et al 2008)
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'cupForcingEns3dMid' ! Subroutine Name
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte, maxens, maxens3
      integer, intent(in) :: ichoice
      !! flag if only want one closure

      integer, intent(in) :: k22(:)
      !! updraft originating level
      integer, intent(in) :: kbcon(:)
      !! LFC of parcel from k22
      integer, intent(in) :: ktop(:)
      integer, intent(in) :: kpbl(:)

      real, intent(in) :: po_cup(:, :)
      real, intent(in) :: dhdt(:, :)
      real, intent(in) :: hc(:, :)
      real, intent(in) :: hco(:, :)
      real, intent(in) :: he_cup(:, :)
      real, intent(in) :: heo_cup(:, :)
      real, intent(in) :: aa1_bl(:)
      real, intent(in) :: tau_ecmwf(:)
      real, intent(in) :: wlpool(:)
      real, intent(in) :: xaa0(:)
      !! cloud work function with cloud effects
      real, intent(in) :: aa1(:)
      !! cloud work function with forcing effects
      real, intent(in) :: zws(:)
      real, intent(in) :: mbdt(:)
      !! arbitrary numerical parameter
      real, intent(in) :: aa0(:)
      !! cloud work function without forcing effects
      real, intent(in) :: dtime
      !! dt over which forcing is applied

      integer, intent(inout) :: ierr(:)

      real, intent(inout) :: xf_dicycle(:)
      real, intent(inout) :: xf_coldpool(:)

      real, intent(out)   :: xff_mid(:, :)

      !Local variables:
      real, dimension(1:maxens) :: xk
      integer :: i, k
      real :: xff_dicycle, trash, blqe, xff_ens1, mf_ens1

      do i = its, itf
         !-initialization
         xff_mid(i, :) = 0.
         xf_dicycle(i) = 0.

         if (ierr(i) /= 0) cycle

         xk(1) = (xaa0(i) - (aa1(i)))/mbdt(i)

         if (xk(1) .le. 0 .and. xk(1) .gt. -0.1*mbdt(i)) xk(1) = -0.1*mbdt(i)
         if (xk(1) .gt. 0 .and. xk(1) .lt. 1.e-2) xk(1) = 1.e-2

         !- closure 3 for mid
         if (xk(1) < 0.) xff_mid(i, 3) = max(0., -(aa1(i)/tau_ecmwf(i))/xk(1))
      end do

      do i = its, itf
         if (ierr(i) /= 0) cycle
         !- Boundary layer quasi-equilibrium (Raymond 1995)
         if (k22(i) .lt. kpbl(i) + 1) then
            blqe = 0.
            do k = kts, kbcon(i)
               blqe = blqe + 100.*dhdt(i, k)*(po_cup(i, k) - po_cup(i, k + 1))/c_grav
            end do
            trash = max((hco(i, kbcon(i)) - heo_cup(i, kbcon(i))), 1.e1)
            xff_mid(i, 2) = max(0., blqe/trash)
         end if

         !- W* closure (Grant,2001)
         xff_mid(i, 1) = 0.03*zws(i)
      end do

   end subroutine cupForcingEns3dMid

   ! ------------------------------------------------------------------------------------
   subroutine cupMinimi(array, ks, kend, kt, ierr, itf, ktf, its, ite, kts, kte)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'cupMinimi' ! Subroutine name
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte

      integer, intent(in) :: ierr(:)
      integer, intent(in) :: ks(:)
      !! check-range
      integer, intent(in) :: kend(:)
      !! check-range
      real, intent(in) :: array(:, :)
      !! array input array
      integer, intent(out) :: kt(:)
      !! kt output array of levels

      !Local variables:
      ! only local dimensions are need as of now in this routine
      integer :: i, k, kstop
      real, dimension(its:ite) :: x

       do i = its, itf
         kt(i) = ks(i)
         if (ierr(i) == 0) then
            x(i) = array(i, ks(i))
            kstop = max(ks(i) + 1, kend(i))
            !
            do k = ks(i) + 1, kstop
               if (array(i, k) < x(i)) then
                  x(i) = array(i, k)
                  kt(i) = k
               end if
            end do
         end if
      end do

   end subroutine cupMinimi
      
   ! ------------------------------------------------------------------------------------
   subroutine cupUpAa0(aa0, z_cup, zu, dby, gamma_cup, t_cup,k22, klcl, kbcon, ktop, ierr , itf, ktf, its, ite, kts, kte &
                     , integ_interval)
      !! ## Determine the cloud work function (J/kg) from cloud base to cloud top
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Determine the cloud work function (J/kg) from cloud base to cloud top
      !! It also may be used to calculate two additional CWF
      !! a) the inhibition energy barrier  (from surface to cloud base)
      !! b) in the boundary layer (from surface to cloud base  height)
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'cupUpAa0' ! Subroutine name
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte

      integer, intent(in) :: k22(:)
      integer, intent(in) :: klcl(:)
      integer, intent(in) :: kbcon(:)
      integer, intent(in) :: ktop(:)
      integer, intent(in) :: ierr (:)

      real, intent(in) :: z_cup(:, :)
      real, intent(in) :: zu(:, :)
      !! normalized updraft mass flux
      real, intent(in) :: gamma_cup(:, :)
      !! gamma on model cloud levels
      real, intent(in) :: t_cup(:, :)
      !! temperature (Kelvin) on model cloud levels
      real, intent(in) :: dby(:, :)
      !! buoancy term

      character(len=*), optional, intent(in) :: integ_interval

      real, intent(out)  :: aa0(:)
      !! cloud work function

      !Local variables:
      integer :: i, k
      real :: dz, da, aa_2, aa_1
      integer, dimension(its:ite) ::  kbeg, kend

      !  initialize array to zero.
      aa0(:) = 0.
      !  set domain of integration
      if (present(integ_interval)) then
         if (trim(integ_interval) == 'BL') then
            kbeg(:) = kts
            kend(:) = kbcon(:) - 1
         elseif(trim(integ_interval) == 'CIN') then
            kbeg(:) = kts
            kend(:) = kbcon(:) ! kbcon(:)-1
         else
            stop "unknown range in cup_up_aa0"
         end if
      else
         kbeg(:) = kbcon(:)
         kend(:) = ktop(:)
      end if

      do i = its, itf
         if (ierr(i) /= 0) cycle
         do k = kbeg(i), kend(i)
            dz = z_cup(i, k + 1) - z_cup(i, k)
            aa_1 = zu(i, k)*(c_grav/(real(c_cp)*t_cup(i, k)))*dby(i, k)/(1.+gamma_cup(i, k))
            aa_2 = zu(i, k + 1)*(c_grav/(real(c_cp)*t_cup(i, k + 1)))*dby(i, k + 1)/(1.+gamma_cup(i, k + 1))
            da = 0.5*(aa_1 + aa_2)*dz
            aa0(i) = aa0(i) + da
            !aa0(i)=aa0(i)+max(0.,da)
         end do
      end do

   end subroutine cupUpAa0

   ! ------------------------------------------------------------------------------------
   subroutine cupUpMoisture(name, start_level, klcl, ierr, ierrc, z_cup, qc, qrc, pw, pwav, hc, tempc, xland &
                           ,po, p_cup, kbcon, ktop, cd, dby, clw_all, t_cup, q_env, gamma_cup, zu, qes_cup, k22, qe_cup &
                           ,zqexec, use_excess, ccn, rho, up_massentr, up_massdetr, psum, psumh, c1d, x_add_buoy &
                           ,vvel2d, vvel1d, zws, entr_rate_2d, itest, itf, ktf, ipr, jpr, its, ite, kts, kte)
      !! ## Resolve the 1-d cloud model for the updraft
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Resolve the 1-d cloud model for the updraft
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'cupUpMoisture' ! Subroutine name

      real, parameter :: p_bdispm = 0.366       
      !! berry--size dispersion (maritime)
      real, parameter :: p_dispc = 0.146       
      !! berry--size dispersion (continental)
      real, parameter :: p_t_bf = 268.16, p_t_ice_bf = 235.16
      real, parameter :: p_rk = 3 
      !! Runge-Kutta order: 2 or 3
      real, parameter :: p_xexp = 2.
      !
      real, parameter :: p_frac = 1.
      !
      !Variables (input, output, inout)
      integer, intent(in) :: use_excess, itest, itf, ktf, its, ite, ipr, jpr, kts, kte

      integer, intent(in) :: kbcon(:)
      integer, intent(in) :: ktop(:)
      integer, intent(in) :: k22(:)
      integer, intent(in) :: klcl(:)
      integer, intent(in) :: start_level(:)

      real, intent(in) :: t_cup(:, :)
      real, intent(in) :: p_cup(:, :)
      real, intent(in) :: rho(:, :)
      real, intent(in) :: q_env(:, :)
      !! environmental q on model levels
      real, intent(in) :: zu(:, :)
      !! normalized updraft mass flux
      real, intent(in) :: gamma_cup(:, :)
      !! gamma on model cloud levels
      real, intent(in) :: qe_cup(:, :)
      !! environmental q on model cloud levels
      real, intent(in) :: hc(:, :)
      real, intent(in) :: po(:, :)
      real, intent(in) :: up_massentr(:, :)
      real, intent(in) :: up_massdetr(:, :)
      real, intent(in) :: dby(:, :)
      !! buoancy term
      real, intent(in) :: qes_cup(:, :)
      !! saturation q on model cloud levels
      real, intent(in) :: z_cup(:, :)
      real, intent(in) :: cd(:, :)
      !! detrainment function
      real, intent(in) :: c1d(:, :)
      real, intent(in) ::  entr_rate_2d(:, :)
      real, intent(in) ::  vvel2d(:, :)
      real, intent(in) :: zqexec(:)
      real, intent(in) :: xland(:)
      real, intent(in) :: x_add_buoy(:)
      real, intent(in) :: vvel1d(:)
      real, intent(in) :: zws(:)
      real, intent(in) :: ccn(:)

      character(len=*), intent(in) ::  name

      integer, intent(inout) :: ierr(:)
      !! ierr error value, maybe modified in this routine

      character(len=128), intent(inout) :: ierrc(:)
      !! ierr text message
      
      real, intent(out) :: qc(:, :)
      !! cloud q (including liquid water) after entrainment
      real, intent(out) :: qrc(:, :)
      !! liquid water content in cloud after rainout
      real, intent(out) :: pw(:, :)
      !! condensate that will fall out at that level
      real, intent(out) :: clw_all(:, :)
      real, intent(out) :: tempc(:, :)
      real, intent(out) :: pwav(:)
      !! totan normalized integrated condensate (I1)
      real, intent(out) :: psum(:)
      real, intent(out) :: psumh(:)

      !Local variables:
      integer :: iounit, iprop, i, k, k1, k2, n, nsteps
      real :: dp, rhoc, dh, dz, radius, berryc0, q1, berryc
      real :: qaver, denom, aux, cx0, qrci, step, cbf, qrc_crit_BF, min_liq, qavail
      real :: delt, tem1, qrc_0, cup
      real :: qrch
      !! saturation q in cloud


      !--- initialization
      pwav    = 0.
      psum    = 0.
      psumh   = 0.
      pw      = 0.
      clw_all = 0.
      tempc   = t_cup
      qrc     = 0.     !--- liq/ice water
      qc      = qe_cup !--- total water: liq/ice = vapor water

      !--- get boundary condition for qc
      do i = its, itf
         if (ierr(i) /= 0) cycle
         call getCloudBc(name, kts, kte, ktf, xland(i), po(i, kts:kte), qe_cup(i, kts:kte), qaver, k22(i))
         qc (i, kts:start_level(i)) = qaver + zqexec(i) + p_frac*x_add_buoy(i)/real(c_alvl)
         qrc(i, kts:start_level(i)) = 0.
      end do

      !--- option to produce linear fluxes in the sub-cloud layer.
      if (trim(name) == 'shallow' .and. USE_LINEAR_SUBCL_MF == 1) then
         do i = its, itf
            if (ierr(i) /= 0) cycle
            call getDelmix(name, kts, kte, ktf, xland(i), start_level(i), po(i, kts:kte), qe_cup(i, kts:kte), qc(i, kts:kte))
         end do
      end if
      do i = its, itf
         if (ierr(i) /= 0) cycle

         do k = start_level(i) + 1, ktop(i) + 1

            dz = z_cup(i, k) - z_cup(i, k - 1)
            !
            !--- saturation  in cloud, this is what is allowed to be in it
            !
            qrch = qes_cup(I, K) + (1./real(c_alvl))*(gamma_cup(i, k)/(1.+gamma_cup(i, k)))*dby(i, k)

            !-    1. steady state plume equation, for what could
            !-       be in cloud without condensation
            denom = (zu(i, k - 1) - .5*up_massdetr(i, k - 1) + up_massentr(i, k - 1)) + 1.e-12
            !if (denom > 0.) then

            qc(i, k) = (qc(i, k - 1)*zu(i, k - 1) - .5*up_massdetr(i, k - 1)*qc(i, k - 1) + up_massentr(i, k - 1) &
                       * q_env(i, k - 1))/denom

            if (k == start_level(i) + 1) &
            qc(i, k) = qc(i, k) + (zqexec(i) + p_frac*x_add_buoy(i)/real(c_alvl)) * up_massentr(i, k - 1)/denom
            
            !--- assuming no liq/ice water in the environment
            qrc(i, k) = (qrc(i, k - 1)*zu(i, k - 1) - .5*up_massdetr(i, k - 1)*qrc(i, k - 1))/denom
            
            !else
            !   qc(i, k) = qc(i, k - 1)
            !   qrc(i, k) = qrc(i, k - 1)
            !end if

            !-- updraft temp
            tempc(i, k) = (1./real(c_cp))*(hc(i, k) - c_grav*z_cup(i, k) - real(c_alvl)*qrch)

            !--- total condensed water before rainout
            clw_all(i, k) = max(0., qc(i, k) - qrch)

            qrc(i, k) = min(clw_all(i, k), qrc(i, k))

            !--- production term => condensation/diffusional growth
            cup = max(0., qc(i, k) - qrch - qrc(i, k))/dz

            if (c0 < 1.e-6) then
               qrc(i, k) = clw_all(i, k)
               qc(i, k) = qrc(i, k) + min(qc(i, k), qrch)
               pwav(i) = 0.
               psum(i) = psum(i) + clw_all(i, k)*zu(i, k)*dz
               cycle
            end if

            if (AUTOCONV == 1) then
               min_liq = QRC_CRIT*(xland(i)*1.+(1.-xland(i))*0.7)
               if (name .eq. 'mid') min_liq = min_liq*0.5

               cx0 = (c1d(i, k) + c0)*DZ
               qrc(i, k) = clw_all(i, k)/(1.+cx0)
               pw(i, k) = cx0*max(0., qrc(i, k) - min_liq)! units kg[rain]/kg[air]
               !pw (i,k)= cx0*qrc(i,k)    ! units kg[rain]/kg[air]

               !--- convert pw to normalized pw
               pw(i, k) = pw(i, k)*zu(i, k)

            elseif (AUTOCONV == 5) then
               !  C0_DEEP     = 1.5e-3; C0_MID     = 1.5e-3 ; QRC_CRIT        = 1.e-4 !(kg/kg)

               min_liq = QRC_CRIT*(xland(i)*0.4 + (1.-xland(i))*1.)

               if (clw_all(i, k) <= min_liq) then !=> more heating at upper levels, more detrained ice

                  qrc(i, k) = clw_all(i, k)
                  pw(i, k) = 0.
               else

                  cx0 = (c1d(i, k) + c0)*(1.+0.33*FractLiqF(tempc(i, k)))
                  !cx0     = (c1d(i,k)+c0)*(1.+ 2.*FractLiqF(tempc(i,k)))
                  !--- v0
                  qrc(i, k) = qrc(i, k)*exp(-cx0*dz) + (cup/cx0)*(1.-exp(-cx0*dz))
                  qrc(i, k) = max(qrc(i, k), min_liq)
                  pw(i, k) = max(0., clw_all(i, k) - qrc(i, k)) ! units kg[rain]/kg[air]
                  qrc(i, k) = clw_all(i, k) - pw(i, k)
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
                  pw(i, k) = pw(i, k)*zu(i, k)
               end if

            elseif (AUTOCONV == 6) then
               min_liq = 0.5*QRC_CRIT*(xland(i)*1.5 + (1.-xland(i))*2.5)

               if (clw_all(i, k) <= min_liq) then
                  qrc(i, k) = clw_all(i, k)
                  pw(i, k) = 0.
               else
                  cx0 = (c1d(i, k) + c0)*dz
                  qrc(i, k) = (clw_all(i, k))*exp(-cx0)
                  pw(i, k) = clw_all(i, k) - qrc(i, k)
                  pw(i, k) = pw(i, k)*zu(i, k)
               end if
               !
               !print*,"6mass=",pw(i,k)/(1.e-12+zu(i,k))+qrc(i,k),clw_all(i,k)
            elseif (AUTOCONV == 7) then
               min_liq = 0.5*QRC_CRIT*(xland(i)*1.5 + (1.-xland(i))*2.5)

               if (clw_all(i, k) <= min_liq) then
                  qrc(i, k) = clw_all(i, k)
                  pw(i, k) = 0.
               else
                  cx0 = c1d(i, k) + c0
                  qrc_0 = qrc(i, k)
                  qrc(i, k) = qrc_0*exp(-cx0*dz) + (cup/cx0)*(1.-exp(-cx0*dz))

                  pw(i, k) = max(clw_all(i, k) - qrc(i, k), 0.)
                  qrc(i, k) = clw_all(i, k) - pw(i, k)
                  pw(i, k) = pw(i, k)*zu(i, k)
               end if
               !
               !print*,"6mass=",pw(i,k)/(1.e-12+zu(i,k))+qrc(i,k),clw_all(i,k)

            elseif (AUTOCONV == 3) then
               min_liq = QRC_CRIT ! * (xland(i)*1.5+(1.-xland(i))*2.5)

               if (clw_all(i, k) <= min_liq) then
                  qrc(i, k) = clw_all(i, k)
                  pw(i, k) = 0.
               else
                  DELT = -5.
                  if (t_cup(i, k) > 273.16 + DELT) then
                     aux = 1.
                  else
                     aux = 1.*exp(0.07*(t_cup(i, k) - (273.16 + DELT)))
                  end if
                  cx0 = aux*c0
                  !                      cx0     = max(cx0,c0)
                  !                      cx0     = max(cx0,0.25*c0)
                  cx0 = max(cx0, 0.50*c0)
                  qrc_0 = qrc(i, k)
                  qrc(i, k) = qrc_0*exp(-cx0*dz) + (cup/cx0)*(1.-exp(-cx0*dz))
                  qrc(i, k) = min(clw_all(i, k), qrc(i, k))
                  pw(i, k) = clw_all(i, k) - qrc(i, k)
                  pw(i, k) = pw(i, k)*zu(i, k)
                  !if(pw(i,k)<0.) stop " pw<0 autoc 3"
               end if

            elseif (AUTOCONV == 4) then

               min_liq = (xland(i)*0.3 + (1.-xland(i))*0.5)*1.e-3

               if (clw_all(i, k) > min_liq) then

                  tem1 = FractLiqF(tempc(i, k))
                  cbf = 1.
                  if (tempc(i, k) < p_t_bf) cbf = 1.+0.5*sqrt(min(max(p_t_bf - tempc(i, k), 0.), p_t_bf - p_t_ice_bf))
                  !qrc_crit_BF = 0.5e-3/cbf
                  qrc_crit_BF = 3.e-4/cbf
                  cx0 = c0*cbf*(tem1*1.3 + (1.-tem1))/(0.75*min(15., max(vvel2d(i, k), 1.)))

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
                  qrc_0 = qrc(i, k)
                  cx0 = cx0*(1.-exp(-(qrc_0/qrc_crit_BF)**2))
                  qrc(i, k) = qrc_0*exp(-cx0*dz) + (cup/cx0)*(1.-exp(-cx0*dz))

                  pw(i, k) = max(clw_all(i, k) - qrc(i, k), 0.)
                  !--- convert PW to normalized PW
                  pw(i, k) = pw(i, k)*zu(i, k)

                  !if(pw(i,k)<-1.e-12) stop " pw<0 autoc 4"
               else
                  pw(i, k) = 0.0
                  qrc(i, k) = clw_all(i, k)
               end if
            end if
            !- total water (vapor + condensed) in updraft after the rainout
            qc(i, k) = qrc(i, k) + min(qc(i, k), qrch)

            !--- integrated normalized condensates
            pwav(i) = pwav(i) + pw(i, k)
            psum(i) = psum(i) + clw_all(i, k)*zu(i, k)*dz
         end do
         if (pwav(i) < 0.) then
            ierr(i) = 66
            ierrc(i) = "pwav negative"
         end if
      end do

      !--- get back water vapor qc
      do i = its, itf
         if (ierr(i) /= 0) cycle
         do k = kts, ktop(i) + 1
            qc(i, k) = qc(i, k) - qrc(i, k)
            !if(qc(i,k) < 0.)stop " qc negative"
         end do
      end do

   end subroutine cupUpMoisture

   ! ------------------------------------------------------------------------------------
   subroutine cupUpMoistureLight(name, start_level, klcl, ierr, ierrc, z_cup, qc, qrc, pw, pwav, hc, tempc, xland,  po, p_cup &
                               , kbcon, ktop, cd, dby, clw_all, t_cup, q_env, gamma_cup, zu,  qes_cup, k22, qe_cup, zqexec &
                               , use_excess, rho,  up_massentr, up_massdetr, psum, psumh, c1d, x_add_buoy,  itest, itf, ktf, ipr &
                               , jpr, its, ite, kts, kte)
      !! ## Resolve the 1-d cloud model for the updraft, 1st guess for autoconv=4
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Resolve the 1-d cloud model for the updraft, 1st guess for autoconv=4
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'cupUpMoistureLight' ! Subroutine name
      real, parameter :: p_frac = 1.
      !
      !Variables (input, output, inout)
      integer, intent(in) ::  use_excess, itest, itf, ktf, its, ite, ipr, jpr, kts, kte

      integer, intent(in) :: kbcon(:)
      integer, intent(in) :: ktop(:)
      integer, intent(in) :: k22(:)
      integer, intent(in) :: klcl(:)
      integer, intent(in) :: start_level(:)

      character(len=*), intent(in) ::  name

      real, intent(in) :: t_cup(:, :)
      real, intent(in) :: p_cup(:, :)
      real, intent(in) :: rho(:, :)
      real, intent(in) :: q_env(:, :)
      !! environmental q on model levels
      real, intent(in) :: zu(:, :)
      !! normalized updraft mass flux
      real, intent(in) :: gamma_cup(:, :)
      !! gamma on model cloud levels
      real, intent(in) :: qe_cup(:, :)
      !! environmental q on model cloud levels
      real, intent(in) :: hc(:, :)
      real, intent(in) :: po(:, :)
      real, intent(in) :: up_massentr(:, :)
      real, intent(in) :: up_massdetr(:, :)
      real, intent(in) :: dby(:, :)
      !! buoancy term
      real, intent(in) :: qes_cup(:, :)
      !! saturation q on model cloud levels
      real, intent(in) :: z_cup(:, :)
      real, intent(in) :: cd(:, :)
      !! detrainment function
      real, intent(in) :: c1d(:, :)

      real, intent(in) :: zqexec(:)
      real, intent(in) :: xland(:)
      real, intent(in) :: x_add_buoy(:)

      integer, intent(inout) :: ierr(:)
      !! ierr error value, maybe modified in this routine
      
      character(len=128), intent(inout) :: ierrc(:)

      real, intent(out) :: qc(:, :)
      !! cloud q (including liquid water) after entrainment
      real, intent(out) :: qrc(:, :)
      !! liquid water content in cloud after rainout
      real, intent(out) :: pw(:, :)
      !! condensate that will fall out at that level
      real, intent(out) :: clw_all(:, :)
      real, intent(out) :: tempc(:, :)
      real, intent(out) :: pwav(:)
      !! totan normalized integrated condensate (I1)
      real, intent(out) :: psum(:)
      real, intent(out) :: psumh(:)

      !Local variables:
      integer :: iounit, iprop, i, k, k1, k2, n, nsteps, delt, tem1
      real :: dp, rhoc, dh, dz, radius, berryc0, q1, berryc
      real :: qaver, denom, aux, cx0, qrci, step, cbf, qrc_crit_bf, min_liq, qavail, delt_hc_glac
      real :: qrch
      !! saturation q in cloud
      
      !--- initialization
      pwav    = 0.
      psum    = 0.
      psumh   = 0.
      pw      = 0.
      clw_all = 0.
      tempc   = t_cup
      qrc     = 0.     !--- liq/ice water
      qc      = qe_cup !--- total water: liq/ice = vapor water

      !--- get boundary condition for qc
      do i = its, itf
         if (ierr(i) /= 0) cycle
         call getCloudBc(name, kts, kte, ktf, xland(i), po(i, kts:kte), qe_cup(i, kts:kte), qaver, k22(i))
         qc(i, kts:start_level(i)) = qaver + zqexec(i) + p_frac*x_add_buoy(i)/real(c_alvl)
      end do

      do i = its, itf
         if (ierr(i) /= 0) cycle

         do k = start_level(i) + 1, ktop(i) + 1

            dz = z_cup(i, k) - z_cup(i, k - 1)
            !
            !--- saturation  in cloud, this is what is allowed to be in it
            !
            qrch = qes_cup(i, k) + (1./real(c_alvl))*(gamma_cup(i, k)/(1.+gamma_cup(i, k)))*dby(i, k)

            !-    1. steady state plume equation, for what could
            !-       be in cloud without condensation
            denom = (zu(i, k - 1) - .5*up_massdetr(i, k - 1) + up_massentr(i, k - 1)) + 1.e-12
            !if (denom > 0.) then
            
            qc(i, k) = (qc(i, k - 1)*zu(i, k - 1) - .5*up_massdetr(i, k - 1)*qc(i, k - 1) + up_massentr(i, k - 1) &
                        * q_env(i, k - 1))/denom

            if (k == start_level(i) + 1) &
            qc(i, k) = qc(i, k) + (zqexec(i) + p_frac*x_add_buoy(i)/real(c_alvl)) * up_massentr(i, k - 1)/denom
            
            !else
            !   qc(i, k) = qc(i, k - 1)
            !end if

            !--- total condensed water before rainout
            clw_all(i, k) = max(0., qc(i, k) - qrch)
            !--- updraft temp
            tempc(i, k) = (1./real(c_cp))*(hc(i, k) - c_grav*z_cup(i, k) - real(c_alvl)*qrch)

            !--add glaciation effect on the MSE
            if (p_melt_glac) then
               delt_hc_glac = clw_all(i, k)*(1.-FractLiqF(tempc(i, k)))*c_xlf

               tempc(i, k) = tempc(i, k) + (1./real(c_cp))*delt_hc_glac
            end if

            cx0 = (c1d(i, k) + c0)*dz
            if (c0 < 1.e-6) cx0 = 0.

            qrc(i, k) = clw_all(i, k)/(1.+cx0)
            pw(i, k) = cx0*max(0., qrc(i, k) - QRC_CRIT)! units kg[rain]/kg[air]
            !--- convert pw to normalized pw
            pw(i, k) = pw(i, k)*zu(i, k)

            !- total water (vapor + condensed) in updraft after the rainout
            qc(i, k) = qrc(i, k) + min(qc(i, k), qrch)

         end do
      end do

      !- get back water vapor qc
      do i = its, itf
         if (ierr(i) /= 0) cycle
         do k = kts, ktop(i) + 1
            qc(i, k) = qc(i, k) - qrc(i, k)
         end do
      end do
   end subroutine cupUpMoistureLight

   ! ------------------------------------------------------------------------------------
   function SatVap(temp2) result(r_sat_vap)
      !! ## Determine the saturation water vapor mixing ratio
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! —
      !! **Full description**:
      !!
      !! Determine the saturation water vapor mixing ratio
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !!
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'SatVap' ! Function Name
   
      !Variables (input):
      real, intent(in) :: temp2
   
      !Local variables:
      real :: r_sat_vap !Output

      real :: temp, toot, toto, eilog, tsot, &
              ewlog, ewlog2, ewlog3, ewlog4

      temp = temp2 - c_t01
      if (temp .lt. -20.) then   !!!! ice saturation
         toot = c_t00/temp2
         toto = 1./toot
         eilog = -9.09718*(toot - 1) - 3.56654*(log(toot) / log(10.)) + .876793*(1 - toto) + (log(6.1071)/log(10.))
         r_sat_vap = 10**eilog
      else
         tsot = c_t100/temp2
         ewlog = -7.90298*(tsot - 1) + 5.02808* &
                 (log(tsot)/log(10.))
         ewlog2 = ewlog - 1.3816e-07* &
                  (10**(11.344*(1 - (1/tsot))) - 1)
         ewlog3 = ewlog2 + .0081328* &
                  (10**(-3.49149*(tsot - 1)) - 1)
         ewlog4 = ewlog3 + (log(1013.246)/log(10.))
         r_sat_vap = 10**ewlog4
      end if

   end function SatVap

   ! ------------------------------------------------------------------------------------
   subroutine cupUpAa1Bl(version, aa1_bl, aa1_fa, aa1, t, tn, q, qo, dtime, po_cup, z_cup, zu, dby, gamma_cup, t_cup, rho, klcl &
                       , kpbl, kbcon, ktop, ierr, itf, ktf, its, ite, kts, kte, xland, ztexec, xlons, xlats, h_sfc_flux &
                       , le_sfc_flux , tau_bl, tau_ecmwf, t_star, cumulus, tn_bl, qo_bl)
      !! ## Determine the Pcape parameter from Bechtold et al (2014)
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Determine the Pcape parameter from Bechtold et al (2014) and adapted for cloud work function
      !! formulation be Freitas et al (2021)
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'cupUpAa1Bl' ! Subroutine name
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte, version
   
      integer, intent(in) :: klcl(:)
      integer, intent(in) :: kbcon(:)
      integer, intent(in) :: ktop(:)
      integer, intent(in) :: kpbl(:)

      real, intent(in) :: z_cup(:, :)
      real, intent(in) :: zu(:, :)
      !! normalized updraft mass flux
      real, intent(in) :: gamma_cup(:, :)
      !! gamma on model cloud levels
      real, intent(in) :: t_cup(:, :)
      !! temperature (Kelvin) on model cloud levels
      real, intent(in) :: dby(:, :)
      !! buoancy term
      real, intent(in) :: t(:, :)
      real, intent(in) :: tn(:, :)
      real, intent(in) :: q(:, :)
      real, intent(in) :: qo(:, :)
      real, intent(in) :: po_cup(:, :)
      real, intent(in) :: rho(:, :)
      real, intent(in) :: tn_bl(:, :)
      real, intent(in) :: qo_bl(:, :)

      real, intent(in) :: xland(:)
      real, intent(in) :: ztexec(:)
      real, intent(in) :: xlons(:)
      real, intent(in) :: xlats(:)
      real, intent(in) :: h_sfc_flux(:)
      real, intent(in) :: le_sfc_flux(:)
      real, intent(in) :: aa1(:)
      real, intent(in) :: tau_bl(:)
      real, intent(in) :: tau_ecmwf(:)
      real, intent(in) :: dtime
      real, intent(in) :: t_star

      character(len=*), intent(in) :: cumulus

      integer, intent(inout) ::  ierr(:)
      !! ierr error value, maybe modified in this routine
     
      real, intent(out) :: aa1_bl(:)
      real, intent(out) :: aa1_fa(:)

      !Local variables:
      integer :: i, k, iprloc
      real :: dz, da, aa_1, aa_2, tcup, da_bl, a1_bl

      aa1_bl(:) = 0.
      if (version == 0) then
         do i = its, itf
            if (ierr(i) /= 0) cycle
            !***       do k=kts,kbcon(i)
            do k = kts, kpbl(i)
               dz = c_grav*(z_cup(i, k + 1) - z_cup(i, k))
               da = dz*(tn(i, k)*(1.+0.608*qo(i, k)) - t(i, k)*(1.+0.608*q(i, k)))/dtime
               aa1_bl(i) = aa1_bl(i) + da ! Units : J K / (kg seg)
            end do
         end do
      elseif (version == 1) then
         do i = its, itf
            if (ierr(i) /= 0) cycle
            do k = kts, kpbl(i)
               dz = (z_cup(i, k + 1) - z_cup(i, k))
               aa_1 = (c_grav/(real(c_cp)*t_cup(i, k)))*dby(i, k)*zu(i, k)
               aa_2 = (c_grav/(real(c_cp)*t_cup(i, k + 1)))*dby(i, k + 1)*zu(i, k + 1)
               da = 0.5*(aa_1 + aa_2)*dz! Units : J / kg
               aa1_bl(i) = aa1_bl(i) + da
            end do
         end do
      else
         stop "unknown version option in routine: cup_up_aa1bl"
      end if

      return

      aa1_fa(:) = 0.
      do i = its, itf
         if (ierr(i) /= 0) cycle
         do k = kbcon(i), ktop(i)

            dz = z_cup(i, k + 1) - z_cup(i, k)
            aa_1 = (c_grav/(real(c_cp)*((t_cup(i, k)))))*dby(i, k)/(1.+gamma_cup(i, k))*zu(i, k)
            aa_2 = (c_grav/(real(c_cp)*((t_cup(i, k + 1)))))*dby(i, k + 1)/(1.+gamma_cup(i, k + 1))*zu(i, k + 1)
            da = 0.5*(aa_1 + aa_2)*dz

            aa1_fa(i) = aa1_fa(i) + da
         end do
      end do

   end subroutine cupUpAa1Bl

   !------------------------------------------------------------------------------------
   subroutine getLateralMassFlux(itf, ktf, its, ite, kts, kte, min_entr_rate,  ierr, ktop, zo_cup, zuo, cd, entr_rate_2d, po_cup &
                              ,  up_massentro, up_massdetro, up_massentr, up_massdetr,  draft, kbcon, k22, kpbl, up_massentru &
                              , up_massdetru, lambau)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'getLateralMassFlux' ! Subroutine name

      integer, parameter :: p_mass_u_option = 1
      integer, parameter :: p_smooth_depth = 2 
      !! increasing this parameter,
      !! strongly damps the heat/drying rates, precip ...
   
      !Variables (input, output, inout)
      
      integer, intent(in) :: itf, ktf, its, ite, kts, kte
      integer, intent(in) :: ierr(:)
      integer, intent(in) :: ktop(:)
      integer, intent(in) :: kbcon(:)
      integer, intent(in) :: k22(:)
      integer, intent(in) :: kpbl(:)
      
      real, intent(in) :: min_entr_rate
      real, intent(in) :: zo_cup(:, :)
      real, intent(in) :: zuo(:, :)
      real, intent(in) :: po_cup(:, :)      
      real, intent(in), optional  :: lambau(:)

      character(len=*), intent(in) :: draft

      real, intent(inout) :: cd(:, :)
      real, intent(inout) :: entr_rate_2d(:, :)
      
      real, intent(out) :: up_massentro(:, :)
      real, intent(out) :: up_massdetro(:, :)
      real, intent(out) :: up_massentr(:, :)
      real, intent(out) :: up_massdetr(:, :)
      real, intent(out), optional :: up_massentru(:, :)
      real, intent(out), optional :: up_massdetru(:, :)

      !Local variables:
      integer :: i, k, turn, ismooth1, ismooth2, nlay, k_ent
      real :: dz, mass1, mass2, dp, rho, zuo_ave
      logical  :: smooth

      integer ::  incr1 = 1
      integer ::  incr2 = 1

      smooth = .false.
      if (USE_SMOOTH_PROF == 1) smooth = .true.

      up_massentro(:, :) = 0.
      up_massdetro(:, :) = 0.
      if (present(up_massentru) .and. present(up_massdetru)) then
         up_massentru(:, :) = 0.
         up_massdetru(:, :) = 0.
      end if
      nlay = int(kte/90)

      do i = its, itf
         if (ierr(i) /= 0) cycle

         !-will not allow detrainment below the location of the maximum zu
         ! if(draft=='shallow'.or.draft == 'mid') cd(i,1:maxloc(zuo(i,:),1)-2)=0.0

         !-will not allow detrainment below cloud base or in the PBL
         if (draft == 'shallow') then
            cd(i, 1:max(kbcon(i), kpbl(i)) + nlay) = 0.0

         else
            cd(i, 1:maxloc(zuo(i, :), 1) + nlay) = 0.0
         end if

         !- mass entrainment and detrainment are defined on model levels
         do k = kts, maxloc(zuo(i, :), 1)
            !=> below location of maximum value zu -> change entrainment

            dz = zo_cup(i, k + 1) - zo_cup(i, k)
            zuo_ave = 0.5*(zuo(i, k + 1) + zuo(i, k))

            up_massdetro(i, k) = cd(i, k)*dz*zuo_ave

            up_massentro(i, k) = zuo(i, k + 1) - zuo(i, k) + up_massdetro(i, k)
            up_massentro(i, k) = max(up_massentro(i, k), min_entr_rate*dz*zuo_ave)

            !-- check dd_massdetro in case of dd_massentro has been changed above
            up_massdetro(i, k) = -zuo(i, k + 1) + zuo(i, k) + up_massentro(i, k)

            !if(zuo(i,k-1).gt.0.) then
            cd(i, k) = up_massdetro(i, k)/(dz*zuo_ave)
            entr_rate_2d(i, k) = up_massentro(i, k)/(dz*zuo_ave)
            !endif
            !if(draft=='shallow')print*,"ent1=",k,real(entr_rate_2d(i,k),4)!,real((min(zo_cup(i,k_ent)/zo_cup(i,k-1),1.)))

         end do

         !--- limit the effective entrainment rate
         k_ent = maxloc(zuo(i, :), 1)
         do k = k_ent + 1, ktop(i) - 1
            entr_rate_2d(i, k) = entr_rate_2d(i, k_ent)*(min(zo_cup(i, k_ent)/zo_cup(i, k), 1.))
            entr_rate_2d(i, k) = max(min_entr_rate, entr_rate_2d(i, k))
            !if(draft=='shallow')print*,"ent2=",k,real(entr_rate_2d(i,k),4),real((min(zo_cup(i,k_ent)/zo_cup(i,k),1.)))
         end do
         entr_rate_2d(i, ktop(i):kte) = 0.

         !=================
         if (smooth .and. trim(draft) /= 'shallow') then
            !---smoothing the transition zone (from maxloc(zu)-1 to maxloc(zu)+1)

            ismooth1 = max(kts + 2, maxloc(zuo(i, :), 1) - p_smooth_depth)
            ismooth2 = min(ktf - 2, maxloc(zuo(i, :), 1) + p_smooth_depth)
            !if(draft == 'shallow') ismooth1 = max(ismooth1,max(kbcon(i),kpbl(i))+nlay)+1

            do k = ismooth1, ismooth2
               dz = zo_cup(i, k + 1) - zo_cup(i, k)

               zuo_ave = 0.5*(zuo(i, k + 1) + zuo(i, k))

               up_massentro(i, k) = 0.5*(entr_rate_2d(i, k)*dz*zuo_ave + up_massentro(i, k - 1))

               up_massdetro(i, k) = zuo(i, k) + up_massentro(i, k) - zuo(i, k + 1)

               if (up_massdetro(i, k) .lt. 0.) then
                  up_massdetro(i, k) = 0.
                  up_massentro(i, k) = zuo(i, k + 1) - zuo(i, k)
                  entr_rate_2d(i, k) = (up_massentro(i, k))/(dz*zuo_ave)
               endif
               if(zuo_ave > 0.) &
                  cd(i,k)=up_massdetro(i,k)/(dz*zuo_ave)
            end do

            do k = ismooth1, ismooth2
               dz = zo_cup(i, k + 1) - zo_cup(i, k)

               zuo_ave = 0.5*(zuo(i, k + 1) + zuo(i, k))

               up_massdetro(i, k) = 0.5*(cd(i, k)*dz*zuo_ave + up_massdetro(i, k - 1))
               up_massentro(i, k) = zuo(i, k + 1) - zuo(i, k) + up_massdetro(i, k)

               if (up_massentro(i, k) .lt. 0.) then
                  up_massentro(i, k) = 0.
                  up_massdetro(i, k) = zuo(i, k) - zuo(i, k + 1)
                  cd(i, k) = up_massdetro(i, k)/(dz*zuo_ave)
               end if
               if (zuo_ave > 0.) &
                  entr_rate_2d(i, k) = (up_massentro(i, k))/(dz*zuo_ave)
            end do
            !-----end of the transition zone
         end if
         !=================

         do k = maxloc(zuo(i, :), 1) + incr1, ktop(i)
            !=> above location of maximum value zu -> change detrainment
            dz = zo_cup(i, k + 1) - zo_cup(i, k)
            zuo_ave = 0.5*(zuo(i, k + 1) + zuo(i, k))

            up_massentro(i, k) = entr_rate_2d(i, k)*dz*zuo_ave

            up_massdetro(i, k) = zuo(i, k) + up_massentro(i, k) - zuo(i, k + 1)
            up_massdetro(i, k) = max(up_massdetro(i, k), 0.0)
            !-- check up_massentro in case of dd_up_massdetro has been changed above
            up_massentro(i, k) = -zuo(i, k) + up_massdetro(i, k) + zuo(i, k + 1)

            if (zuo_ave .gt. 0.) then
               cd(i, k) = up_massdetro(i, k)/(dz*zuo_ave)
               entr_rate_2d(i, k) = up_massentro(i, k)/(dz*zuo_ave)
            end if
         end do

         do k = kts, kte
            up_massentr(i, k) = up_massentro(i, k)
            up_massdetr(i, k) = up_massdetro(i, k)
         end do
         if (present(up_massentru) .and. present(up_massdetru)) then
            if (p_mass_u_option == 1) then
               do k = kts + 1, kte
                  !--       for weaker mixing
                  up_massentru(i, k - 1) = up_massentro(i, k - 1) + lambau(i)*up_massdetro(i, k - 1)
                  up_massdetru(i, k - 1) = up_massdetro(i, k - 1) + lambau(i)*up_massdetro(i, k - 1)
                  !--       for stronger mixing
                  ! up_massentru(i,k-1)=up_massentro(i,k-1)+lambau(i)*up_massentro(i,k-1)
                  ! up_massdetru(i,k-1)=up_massdetro(i,k-1)+lambau(i)*up_massentro(i,k-1)
               end do
            else
               turn = maxloc(zuo(i, :), 1)
               do k = kts + 1, turn
                  up_massentru(i, k - 1) = up_massentro(i, k - 1) + lambau(i)*up_massentro(i, k - 1)
                  up_massdetru(i, k - 1) = up_massdetro(i, k - 1) + lambau(i)*up_massentro(i, k - 1)
               end do
               do k = turn + 1, kte
                  up_massentru(i, k - 1) = up_massentro(i, k - 1) + lambau(i)*up_massdetro(i, k - 1)
                  up_massdetru(i, k - 1) = up_massdetro(i, k - 1) + lambau(i)*up_massdetro(i, k - 1)
               end do
            end if
         end if
         do k = ktop(i) + 1, kte
            cd(i, k) = 0.
            entr_rate_2d(i, k) = 0.
         end do

      end do ! i
      !---- check mass conservation
      do i = its, itf
         if (ierr(i) /= 0) cycle
         do k = kts + 1, kte

            dz = zo_cup(i, k) - zo_cup(i, k - 1)
            dp = 100*(po_cup(i, k) - po_cup(i, k - 1))
            rho = -dp/dz/c_grav
            mass1 = (zuo(i, k) - zuo(i, k - 1)) - up_massentro(i, k - 1) + up_massdetro(i, k - 1)
            !print*,"masscons=",mass1!,-rho*g*(zuo(i,k)-zuo(i,k-1))/dp, (zuo(i,k)-zuo(i,k-1))/dz,( up_massentro(i,k-1)-up_massdetro(i,k-1))/dz,rho
            mass2 = (zuo(i, k) - zuo(i, k - 1)) - up_massentru(i, k - 1) + up_massdetru(i, k - 1)
         end do
      end do
   end subroutine getLateralMassFlux

   ! ------------------------------------------------------------------------------------
   subroutine getLateralMassFluxDown(cumulus, itf, ktf, its, ite, kts, kte, ierr, jmin, zo_cup, zdo, xzd, zd, cdd, mentrd_rate_2d &
                                 ,   dd_massentro, dd_massdetro, dd_massentr, dd_massdetr, draft, mentrd_rate, dd_massentru &
                                 ,   dd_massdetru, lambau)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'getLateralMassFluxDown' ! Subroutine name
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte
      
      integer, intent(in) :: ierr(:)
      integer, intent(in) :: jmin(:)
      
      real, intent(in) :: zo_cup(:, :)
      real, intent(in) :: zdo(:, :)
      real, intent(in) :: mentrd_rate(:)
      real, intent(in) :: lambau(:)

      character(len=*), intent(in) :: draft
      character(len=*), intent(in) :: cumulus

      real, intent(inout) :: cdd(:, :)
      real, intent(inout) :: mentrd_rate_2d(:, :)
      real, intent(inout) :: xzd(:, :)
      real, intent(inout) :: zd(:, :)

      real, intent(out) :: dd_massentro(:, :)
      real, intent(out) :: dd_massdetro(:, :)
      real, intent(out) :: dd_massentr(:, :)
      real, intent(out) :: dd_massdetr(:, :)

      real, intent(out), optional :: dd_massentru(:, :)
      real, intent(out), optional :: dd_massdetru(:, :)

      !Local variables:
      integer ::i, ki
      real :: dzo

      cdd = 0.
      dd_massentr = 0.
      dd_massdetr = 0.
      dd_massentro = 0.
      dd_massdetro = 0.
      if (present(dd_massentru) .and. present(dd_massdetru)) then
         dd_massentru = 0.
         dd_massdetru = 0.
      end if
      if (trim(cumulus) == 'shallow') return

      do i = its, itf
         if (ierr(i) /= 0) cycle

         mentrd_rate_2d(i, 1:jmin(i)) = mentrd_rate(i)
         cdd(i, 1:jmin(i) - 1) = mentrd_rate(i)
         mentrd_rate_2d(i, 1) = 0.

         do ki = jmin(i), maxloc(zdo(i, :), 1), -1
            !=> from jmin to maximum value zd -> change entrainment
            dzo = zo_cup(i, ki + 1) - zo_cup(i, ki)
            dd_massdetro(i, ki) = cdd(i, ki)*dzo*zdo(i, ki + 1)
            !XXX
            dd_massentro(i, ki) = zdo(i, ki) - zdo(i, ki + 1) + dd_massdetro(i, ki)
            dd_massentro(i, ki) = max(0., dd_massentro(i, ki))
            !-- check dd_massdetro in case of dd_massentro has been changed above
            dd_massdetro(i, ki) = dd_massentro(i, ki) - zdo(i, ki) + zdo(i, ki + 1)
            !~ if(dd_massentro(i,ki).lt.0.)then
            !~ dd_massentro(i,ki)=0.
            !~ dd_massdetro(i,ki)=zdo(i,ki+1)-zdo(i,ki)
            !~ if(zdo(i,ki+1) > 0.0)&
            !~ cdd(i,ki)=dd_massdetro(i,ki)/(dzo*zdo(i,ki+1))
            !~ endif
            !~ if(zdo(i,ki+1) > 0.0)&
            !~ mentrd_rate_2d(i,ki)=dd_massentro(i,ki)/(dzo*zdo(i,ki+1))
         end do

         do ki = maxloc(zdo(i, :), 1) - 1, kts, -1
            !=> from maximum value zd to surface -> change detrainment
            dzo = zo_cup(i, ki + 1) - zo_cup(i, ki)
            dd_massentro(i, ki) = mentrd_rate_2d(i, ki)*dzo*zdo(i, ki + 1)
            !XXX
            dd_massdetro(i, ki) = zdo(i, ki + 1) + dd_massentro(i, ki) - zdo(i, ki)
            dd_massdetro(i, ki) = max(0.0, dd_massdetro(i, ki))
            !-- check dd_massentro in case of dd_massdetro has been changed above
            dd_massentro(i, ki) = dd_massdetro(i, ki) + zdo(i, ki) - zdo(i, ki + 1)
            !~ if(dd_massdetro(i,ki).lt.0.)then
            !~ dd_massdetro(i,ki)=0.
            !~ dd_massentro(i,ki)=zdo(i,ki)-zdo(i,ki+1)
            !~ if(zdo(i,ki+1) > 0.0)&
            !~ mentrd_rate_2d(i,ki)=dd_massentro(i,ki)/(dzo*zdo(i,ki+1))
            !~ endif
            !~ if(zdo(i,ki+1) > 0.0)&
            !~ cdd(i,ki)= dd_massdetro(i,ki)/(dzo*zdo(i,ki+1))
         end do

         do ki = jmin(i), kts, -1
            xzd(i, ki) = zdo(i, ki)
            zd(i, ki) = zdo(i, ki)
            dd_massentr(i, ki) = dd_massentro(i, ki)
            dd_massdetr(i, ki) = dd_massdetro(i, ki)
         end do
         if (present(dd_massentru) .and. present(dd_massdetru)) then
            do ki = jmin(i), kts, -1
               dd_massentru(i, ki) = dd_massentro(i, ki) + lambau(i)*dd_massdetro(i, ki)
               dd_massdetru(i, ki) = dd_massdetro(i, ki) + lambau(i)*dd_massdetro(i, ki)
            end do
         end if
      end do

   end subroutine getLateralMassFluxDown

   ! -------------------------------------------------------------------------------------
   subroutine getZuZdPdf(cumulus, draft, ierr, kb, kt, zu, kts, kte, ktf, kpbli, k22, kbcon, klcl, po_cup, psur, xland, random)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
     implicit none
     !Parameters:
     character(len=*), parameter :: p_procedure_name = 'getZuZdPdf'

      real, parameter :: p_px = 45./120. 
      !! px sets the pressure level of max zu. its range is from 1 to 120.
      real, parameter :: p_px2 = 45./120. 
      !! px sets the pressure level of max zu. its range is from 1 to 120.
      real, parameter :: p_beta_deep = 1.25, p_g_beta_deep = 0.8974707

      !-------- gama pdf
      real, parameter, dimension(30) :: x_alpha = (/ &
         3.699999, 3.699999, 3.699999, 3.699999, &
         3.024999, 2.559999, 2.249999, 2.028571, 1.862500, &
         1.733333, 1.630000, 1.545454, 1.475000, 1.415385, &
         1.364286, 1.320000, 1.281250, 1.247059, 1.216667, &
         1.189474, 1.165000, 1.142857, 1.122727, 1.104348, &
         1.087500, 1.075000, 1.075000, 1.075000, 1.075000, &
         1.075000/)
      real, parameter, dimension(30) :: g_alpha = (/ &
         4.1706450, 4.1706450, 4.1706450, 4.1706450, &
         2.0469250, 1.3878370, 1.1330030, 1.012418, 0.9494680, &
         0.9153771, 0.8972442, 0.8885444, 0.8856795, 0.8865333, &
         0.8897996, 0.8946404, 0.9005030, 0.9070138, 0.9139161, &
         0.9210315, 0.9282347, 0.9354376, 0.9425780, 0.9496124, &
         0.9565111, 0.9619183, 0.9619183, 0.9619183, 0.9619183, &
         0.9619183/)
      !-------- gama pdf
  
     !Variables (input, output, inout)
      integer, intent(in) :: kts, kte, ktf, kpbli, k22, kbcon, kt, kb, klcl

      character(len=*), intent(in) :: draft
      character(len=*), intent(in) :: cumulus

      real, intent(in) :: po_cup(:)
      real, intent(in) :: psur
      real, intent(in) :: xland
      real, intent(in) :: random

      integer, intent(inout) :: ierr

      real, intent(inout) :: zu(:)
      
     !Local variables:
      integer :: kk, add, i, k, kb_adj, kpbli_adj, level_max_zu, ktarget
      integer:: itest                   
      !! 5=gamma+beta, 4=gamma, 1=beta
      integer:: minzu, maxzul, maxzuh, kstart
      integer :: nrec = 0
      integer :: k1
      real :: zumax, ztop_adj, a2, beta, alpha, kratio, tunning, fzu, krmax, dzudk, hei_updf, hei_down
      real :: lev_start, g_alpha2, g_a, y1, x1, g_b, a, b, alpha2, csum, zubeg, wgty, dp_layer, slope
      real :: zuh(kts:kte), zul(kts:kte)
      real :: pmaxzu      
      !! pressure height of max zu for deep
     
      logical :: do_smooth

      DO_SMOOTH = .false.
      if (USE_SMOOTH_PROF == 1) DO_SMOOTH = .true.

      !-- fill zu with zeros
      itest = -999
      zu = 0.0
      zuh = 0.0
      zul = 0.0

      if (trim(draft) == "deep_up") itest = CUM_ZUFORM(p_deep)  !ocean/land
      if (trim(draft) == "mid_up") itest = CUM_ZUFORM(p_mid)

      !---------------------------------------------------------
      if (itest == 5 .and. trim(draft) == "mid_up") then

         !--- part 1 GAMMA format
         csum = 0.
         zubeg = 0.
         lev_start = min(.9, .1 + csum*.013)
         kb_adj = max(kb, 2)
         kb_adj = min(kb_adj, kt - 1)
         if (kb_adj == kt) stop "kb_adj==kt"

         tunning = 0.30
         alpha2 = (tunning*(p_beta_deep - 2.) + 1.)/(1.-tunning)

         do k = 27, 3, -1
            if (x_alpha(k) >= alpha2) exit
         end do
         k1 = k + 1
         if (x_alpha(k1) .ne. x_alpha(k1 - 1)) then
            a = x_alpha(k1) - x_alpha(k1 - 1)
            b = x_alpha(k1 - 1)*(k1) - (k1 - 1)*x_alpha(k1)
            x1 = (alpha2 - b)/a
            y1 = a*x1 + b
            g_a = g_alpha(k1) - g_alpha(k1 - 1)
            g_b = g_alpha(k1 - 1)*k1 - (k1 - 1)*g_alpha(k1)
            g_alpha2 = g_a*x1 + g_b
         else
            g_alpha2 = g_alpha(k1)
         end if

         fzu = GammaBrams(alpha2 + p_beta_deep)/(g_alpha2*p_g_beta_deep)
         fzu = 0.01*fzu
         do k = kb_adj, min(kte, kt)
            kratio = (po_cup(k) - po_cup(kb_adj))/(po_cup(kt) - po_cup(kb_adj))
            zu(k) = zubeg + fzu*kratio**(alpha2 - 1.0)*(1.0 - kratio)**(p_beta_deep - 1.0)

         end do
         !- normalize ZU
         zu(kts:min(kte, kt + 1)) = zu(kts:min(kte, kt + 1))/(1.e-12 + maxval(zu(kts:min(kte, kt + 1))))

         !--- part 2: BETA format
         pmaxzu = psur - p_px*(psur - po_cup(kt))
         kb_adj = minloc(abs(po_cup(kts:kt) - pmaxzu), 1)
         kb_adj = max(kb, kb_adj)
         kb_adj = min(kb_adj, kt)
         !beta=4.  !=> must be larger than 1
         !=> higher makes the profile sharper
         !=> around the maximum zu
         !- 2nd approach for beta and alpha parameters
         !- the tunning parameter must be between 0.5 (low  level max zu)
         !-                                   and 1.5 (high level max zu)
         !tunning= 1.0
         tunning = 0.6
         !
         beta = 2.0/tunning
         alpha = tunning*beta
         !
         !- this alpha constrains the location of the maximun ZU to be at
         !- "kb_adj" vertical level
         alpha = 1.+(beta - 1.0)*(float(kb_adj)/float(kt + 1))/(1.0 - (float(kb_adj)/float(kt + 1)))
         !
         ! imposing zu(ktop) = 0
         do k = klcl - 1, min(kte, kt)
            kratio = float(k)/float(kt + 1)
            zuh(k) = kratio**(alpha - 1.0)*(1.0 - kratio)**(beta - 1.0)
         end do
         !- normalize ZU
         zuh(kts:min(kte, kt + 1)) = zuh(kts:min(kte, kt + 1))/(1.e-12 + maxval(zuh(kts:min(kte, kt + 1))))

         !--- part 3: BETA format from sfc to max zuh, then GAMMA format
         do k = kts, max(kts, maxloc(zuh(:), 1) - 2)
            zu(k) = zuh(k)
         end do
         do k = max(kts, maxloc(zuh(:), 1) - 1), min(maxloc(zuh(:), 1) + 1, kt)
            zu(k) = 0.5*(zu(k) + zuh(k))
         end do

         !-- special treatment below k22/klcl
         do k = klcl, kts + 1, -1
            zu(k) = zu(k + 1)*0.5
         end do
         !-- smooth section
         if (do_smooth) then
            !--from surface
            zul(kts + 1) = zu(kts + 1)*0.25
            do k = kts + 2, maxloc(zu, 1)
               zul(k) = (zu(k - 1) + zu(k))*0.5
            end do
            do k = kts + 1, maxloc(zu, 1)
               zu(k) = (zul(k) + zu(k))*0.5
            end do

            !--from ktop
            zul(kt - 1) = zu(kt - 1)*0.1
            !print*,"ZUMD=",kt,zu(kt),zul(kt)
            do k = kt - 2, max(kt - min(maxloc(zu, 1), 5), kts), -1
               zul(k) = (zul(k + 1) + zu(k))*0.5
            end do
            wgty = 0.
            do k = kt, max(kt - min(maxloc(zu, 1), 5), kts), -1
               wgty = wgty + 1./(float(min(maxloc(zu, 1), 5)) + 1)
               zu(k) = zul(k)*(1.-wgty) + zu(k)*wgty
               !print*,"zuMD=",k,zu(k),zul(k),(zul(k)+zu(k))*0.5,min(maxloc(zu,1),5),wgty
            end do
         end if
         zu(kts) = 0.
         !---------------------------------------------------------
      elseif (itest == 20) then       !--- land/ocean

         hei_updf = (1.-xland)*hei_updf_land + xland*hei_updf_ocean
         !- add a randomic perturbation
         hei_updf = hei_updf + random

         !- for gate soundings
         !hei_updf = max(0.1, min(1.,float(JL)/100.))
         !beta =1.0+float(JL)/100. * 5.

         !--hei_updf parameter goes from 0 to 1 = rainfall decreases with hei_updf
         pmaxzu = (psur - 100.)*(1.-0.5*hei_updf) + 0.6*(po_cup(kt))*0.5*hei_updf

         !- beta parameter: must be larger than 1, higher makes the profile sharper around the maximum zu
         !beta    = max(1.1, 2.1 - 0.5*hei_updf)
         beta = 2.2

         kb_adj = minloc(abs(po_cup(kts:kt) - pmaxzu), 1)
         kb_adj = max(kb, kb_adj)
         kb_adj = min(kb_adj, kt)

         !- this alpha constrains the location of the maximun ZU to be at "kb_adj" vertical level
         alpha = 1.+(beta - 1.0)*(float(kb_adj)/float(kt + 1))/(1.0 - (float(kb_adj)/float(kt + 1)))

         !
         do k = klcl - 1, min(kte, kt)
            kratio = float(k)/float(kt + 1)
            zu(k) = kratio**(alpha - 1.0)*(1.0 - kratio)**(beta - 1.0)
         end do

         !-- special treatment below k22/klcl
         do k = klcl, kts + 1, -1
            zu(k) = zu(k + 1)*0.5
         end do

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
         if (do_smooth) then
            !--from surface
            zul(kts + 1) = zu(kts + 1)*0.25
            do k = kts + 2, maxloc(zu, 1)
               zul(k) = zu(k - 1)*0.8 + zu(k)*0.2
            end do
            do k = kts + 1, maxloc(zu, 1)
               zu(k) = (zul(k) + zu(k))*0.5
            end do

            !--from ktop
            zul(kt) = zu(kt)*0.1
            do k = kt - 1, max(kt - min(maxloc(zu, 1), 5), kts), -1
               zul(k) = (zul(k + 1) + zu(k))*0.5
            end do
            wgty = 0.0
            do k = kt, max(kt - min(maxloc(zu, 1), 5), kts), -1
               wgty = wgty + 1./(float(min(maxloc(zu, 1), 5)) + 1)
               zu(k) = zul(k)*(1.-wgty) + zu(k)*wgty
            end do
         end if
         zu(kts) = 0.
!333 continue

         !---------------------------------------------------------
      elseif (itest == 10) then

         if (xland < 0.90) then !- over land
            hei_updf = hei_updf_land
         else
            hei_updf = hei_updf_ocean
         end if

         !- for gate soundings
         !if(gate) hei_updf = max(0.1, min(1.,float(JL)/100.))
         !print*,"JL=",jl,hei_updf

         pmaxzu = 780.
         kb_adj = minloc(abs(po_cup(kts:kt) - pmaxzu), 1)!;print*,"1=",kb_adj,po_cup(kb_adj)
         kb_adj = max(kb, kb_adj)
         kb_adj = min(kb_adj, kt)
         beta = 5.0

         !- this alpha constrains the location of the maximun ZU to be at "kb_adj" vertical level
         alpha = 1.+(beta - 1.0)*(float(kb_adj)/float(kt + 1))/(1.0 - (float(kb_adj)/float(kt + 1)))
         do k = klcl - 1, min(kte, kt)
            kratio = float(k)/float(kt + 1)
            zul(k) = kratio**(alpha - 1.0)*(1.0 - kratio)**(beta - 1.0)
            !print*,"1",k,zul(k),kb_adj,pmaxzu
         end do
         zul(kts:min(kte, kt)) = zul(kts:min(kte, kt))/(1.e-9 + maxval(zul(kts:min(kte, kt)), 1))

         !-----------
         pmaxzu = po_cup(kt) + 300.
         kb_adj = minloc(abs(po_cup(kts:kt) - pmaxzu), 1)
         kb_adj = max(kb, kb_adj)
         kb_adj = min(kb_adj, kt)
         beta = 1.5

         !- this alpha constrains the location of the maximun ZU to be at "kb_adj" vertical level
         alpha = 1.+(beta - 1.0)*(float(kb_adj)/float(kt + 1))/(1.0 - (float(kb_adj)/float(kt + 1)))
         do k = klcl - 1, min(kte, kt)
            kratio = float(k)/float(kt + 1)
            zuh(k) = kratio**(alpha - 1.0)*(1.0 - kratio)**(beta - 1.0)
         end do
         zuh(kts:min(kte, kt)) = zuh(kts:min(kte, kt))/(1.e-9 + maxval(zuh(kts:min(kte, kt)), 1))

         !increasing contribuition of zuh => more heating at upper levels/less precip
         zu(:) = (1.-hei_updf)*zul(:) + hei_updf*zuh(:)

         !-- special treatment below k22/klcl
         do k = klcl, kts + 1, -1
            zu(k) = zu(k + 1)*0.5
         end do
         !-- smooth section
         if (do_smooth) then
            !--from surface
            zul(kts + 1) = zu(kts + 1)*0.25
            do k = kts + 2, maxloc(zu, 1)
               zul(k) = zu(k - 1)*0.8 + zu(k)*0.2
            end do
            do k = kts + 1, maxloc(zu, 1)
               zu(k) = (zul(k) + zu(k))*0.5
            end do

            !--from ktop
            zul(kt) = zu(kt)*0.1
            do k = kt - 1, max(kt - min(maxloc(zu, 1), 5), kts), -1
               zul(k) = (zul(k + 1) + zu(k))*0.5
            end do

            wgty = 0.
            do k = kt, max(kt - min(maxloc(zu, 1), 5), kts), -1
               wgty = wgty + 1./(float(min(maxloc(zu, 1), 5)) + 1)
               zu(k) = zul(k)*(1.-wgty) + zu(k)*wgty
               !print*,"zu=",k,zu(k),zul(k),(zul(k)+zu(k))*0.5,min(maxloc(zu,1),5),wgty
            end do
         end if
         zu(kts) = 0.
         !---------------------------------------------------------
      elseif (trim(draft) == "shallow_up") then
         kb_adj = kts     ! level where mass flux starts
         kpbli_adj = kpbli
         if (kpbli_adj < kb_adj .or. kpbli_adj >= kt) then
            kpbli_adj = kb_adj + 1
         end if

         !- location of the maximum Zu: dp_layer mbar above PBL height
         !dp_layer     = 10. !mbar
         !level_max_zu = minloc(abs(po_cup(kts:kt+1)-(po_cup(kpbli_adj)-dp_layer)),1)
         !

         k1 = max(kbcon, kpbli_adj)
         !- location of the maximum Zu: dp_layer mbar above k1 height
         hei_updf = (1.-xland)*hei_updf_land + xland*hei_updf_ocean

         !hei_updf = (float(JL)-20)/40. ; print*,"JL=",jl,hei_updf

         dp_layer = hei_updf*(po_cup(k1) - po_cup(kt))

         level_max_zu = minloc(abs(po_cup(kts:kt + 1) - (po_cup(k1) - dp_layer)), 1)
         level_max_zu = min(level_max_zu, kt - 1)
         level_max_zu = max(level_max_zu, kts + 1)

         krmax = float(level_max_zu)/float(kt + 1)
         krmax = min(krmax, 0.99)

         beta = BETA_SH!smaller => sharper detrainment layer
         !beta= ((1.-xland)*0.43 +xland)*beta_sh

         !beta= 3.0!smaller => sharper detrainment layer
         !beta = 1.+4.*(float(JL))/40.

         !- this alpha imposes the maximum zu at kpbli
         alpha = 1.+krmax*(beta - 1.)/(1.-krmax)
         !alpha=min(6.,alpha)

         !- to check if dZu/dk = 0 at k=kpbli_adj
         !kratio=krmax
         !dzudk=(alpha-1.)*(kratio**(alpha-2.)) * (1.-kratio)**(beta-1.) - &
         !          (kratio**(alpha-1.))*((1.-kratio)**(beta-2.))*(beta-1.)

         !- Beta PDF
         do k = kts + 1, min(kte, kt)
            kratio = float(k)/float(kt + 1)
            zu(k) = kratio**(alpha - 1.0)*(1.0 - kratio)**(beta - 1.0)
         end do
         zu(kts) = 0.
         !
         !-- special treatment below kbcon - linear Zu
         if (USE_LINEAR_SUBCL_MF == 1) then
            kstart = kbcon
            slope = (zu(kstart) - zu(kts))/(po_cup(kstart) - po_cup(kts) + 1.e-6)
            do k = kstart - 1, kts + 1, -1
               zu(k) = zu(kstart) - slope*(po_cup(kstart) - po_cup(k))
               !print*,"k=",zu(kstart),zu(k),zu(kts)
            end do
         end if
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
      elseif (trim(draft) == "DOWN") then
         if (trim(cumulus) == 'shallow') return
         if (trim(cumulus) == 'mid') beta = 2.5
         if (trim(cumulus) == 'deep') beta = 2.5

         hei_down = (1.-xland)*hei_down_land + xland*hei_down_ocean

         pmaxzu = hei_down*po_cup(kt) + (1.-hei_down)*psur
         kb_adj = minloc(abs(po_cup(kts:kt) - pmaxzu), 1)

         !- this alpha constrains the location of the maximun ZU to be at "kb_adj" vertical level
         alpha = 1.+(beta - 1.0)*(float(kb_adj)/float(kt + 1))/(1.0 - (float(kb_adj)/float(kt + 1)))

         do k = kts + 1, min(kt + 1, ktf)
            kratio = float(k)/float(kt + 1)
            zu(k) = kratio**(alpha - 1.0)*(1.0 - kratio)**(beta - 1.0)
         end do
         !-- smooth section
         if (do_smooth) then
            zul(kts + 1) = zu(kts + 1)*0.1
            wgty = 0.
            do k = kts + 2, maxloc(zu, 1)
               wgty = wgty + 1./(float(max(2, maxloc(zu, 1))) - 1.)
               wgty = 0.5
               !print*,"zD1=",k,zu(k),zul(k-1),wgty,zul(k-1)*(1.-wgty)+ zu(k)*wgty
               zul(k) = zul(k - 1)*(1.-wgty) + zu(k)*wgty
            end do
            wgty = 0.
            do k = kts + 1, maxloc(zu, 1)
               wgty = wgty + 1./(float(max(2, maxloc(zu, 1))) - 1.)
               wgty = 0.5
               !print*,"zD2=",k,zu(k),zul(k),wgty,zul(k)*(1.-wgty)+ zu(k)*wgty
               zu(k) = zul(k)*(1.-wgty) + zu(k)*wgty
            end do
         end if
         zu(kts) = 0.

      end if

      if (maxval(zu(kts:min(kte, kt + 1)), 1) <= 0.0) then
         zu = 0.0
         ierr = 51 !ierr(i)=51
      else
         !- normalize ZU
         zu(kts:min(kte, kt + 1)) = zu(kts:min(kte, kt + 1))/(1.e-9 + maxval(zu(kts:min(kte, kt + 1)), 1))
      end if

     end subroutine getZuZdPdf 
     
   ! ------------------------------------------------------------------------------------
   subroutine cupUpCape(aa0, z, zu, dby, gamma_cup, t_cup, k22, kbcon, ktop, ierr, tempco, qco, qrco, qo_cup &
                     ,  itf, ktf, its, ite, kts, kte)
      !! ## Calculate the CAPE for updraft
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Calculate the CAPE for updraft
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'cupUpCape' ! Subroutine name
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte

      integer, intent(in) :: k22(:)
      integer, intent(in) :: kbcon(:)
      integer, intent(in) :: ktop(:)

      real, intent(in) :: z(:, :)
      !! heights of model levels
      real, intent(in) :: zu(:, :)
      !! normalized updraft mass flux
      real, intent(in) :: gamma_cup(:, :)
      ! gamma on model cloud levels
      real, intent(in) :: t_cup(:, :)
      !! temperature (Kelvin) on model cloud levels
      real, intent(in) :: dby(:, :)
      !! buoancy term
      real, intent(in) :: tempco(:, :)
      !! in-cloud temperature (Kelvin) on model cloud levels
      real, intent(in) :: qco(:, :)
      !! in-cloud water vapor mixing ratio on model cloud levels
      real, intent(in) :: qrco(:, :)
      !! in-cloud liquid water mixing ratio on model cloud levels
      real, intent(in) :: qo_cup(:, :)
      !! environ water vapor mixing ratio on model cloud levels

      integer, intent(inout) :: ierr(:)
      !! error value, maybe modified in this routine

      real, intent(out) :: aa0(:)
      !! dummy array for CAPE (total cape)

      !Local variables:
      integer :: i, k
      real :: dz, daa0
      !
      aa0(:) = 0.
      do i = its, itf
         if (ierr(i) == 0) then
            do k = kbcon(i), ktop(i)
               dz = z(i, k) - z(i, max(1, k - 1))
               daa0 = c_grav*dz*((tempco(i, k)*(1.+0.608*qco(i, k)) - t_cup(i, k)*(1.+0.608*qo_cup(i, k))) /(t_cup(i, k) &
                    * (1.+0.608*qo_cup(i, k))) &
                            )
               aa0(i) = aa0(i) + max(0., daa0)
               !~ print*,"cape",k,AA0(I),tempco(i,k),t_cup(i,k), qrco  (i,k)
            end do
         end if
      end do
   end subroutine cupUpCape

   ! ------------------------------------------------------------------------------------
   subroutine getCloudBc(cumulus, kts, kte, ktf, xland, po, array, x_aver, k22, add, tpert)
      !! ## Determine the various boundary conditions
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! General routine to determine the various boundary conditions for the 
      !! updraft of any plume (deep, shallow, congestus)
      !!
      !! ** History**:
      !! BC_METH == 0 use an arithmetic mean
      !! BC_METH == 1 use a weighted mean (this is the default method)
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'getCloudBc' ! Subroutine name

      real, parameter :: p_frac_ave_layer_ocean = 0.3
   
      !Variables (input, output, inout)
      integer, intent(in) :: kts, kte, ktf, k22

      real, intent(in) :: array(:)
      real, intent(in) :: po(:)
      real, intent(in) :: xland
      real, optional, intent(in) :: add
      real, optional, intent(in) :: tpert(:)

      character(len=*), intent(in) :: cumulus

      real, intent(out) :: x_aver

      !Local variables:
      integer  :: i, local_order_aver, order_aver, i_beg, i_end, ic
      real  :: count, dp, dp_layer, effec_frac, x_ave_layer

      !-- dimensions of the average:
      !-- a) to pick the value at k22 level, instead of an average between
      !--    (k22-order_aver, ..., k22-1, k22) set order_aver=kts
      !-- b) to average between kts and k22 => set order_aver = k22
      !order_aver = 4    !=> BC_METH 0: average between k22, k22-1, k22-2 ...
      !=> BC_METH 1: average between ... k22+1,k22, k22-1 ...
      !-- order_aver = kts !=> average between k22, k22-1 and k22-2

      if (BC_METH == 0) then

         order_aver = 3
         local_order_aver = min(k22, order_aver)

         x_aver = 0.
         do i = kts, local_order_aver
            x_aver = x_aver + array(k22 - i + 1)
         end do
         x_aver = x_aver/float(local_order_aver)

      elseif (BC_METH == 1) then
         effec_frac = (1.-xland) + xland*p_frac_ave_layer_ocean
         x_ave_layer = ave_layer*effec_frac

         i_beg = minloc(abs(po(kts:ktf) - (po(k22) + 0.5*x_ave_layer)), 1)
         i_end = minloc(abs(po(kts:ktf) - (po(k22) - 0.5*x_ave_layer)), 1)
         i_beg = min(ktf, max(i_beg, kts))
         i_end = min(ktf, max(i_end, kts))

         if (i_beg >= i_end) then
            x_aver = array(k22)
            dp_layer = 0.
            ic = i_beg

         else
            dp_layer = 1.e-06
            x_aver = 0.
            ic = 0
            do i = i_beg, ktf
               dp = -(po(i + 1) - po(i))
               if (dp_layer + dp <= x_ave_layer) then
                  dp_layer = dp_layer + dp
                  x_aver = x_aver + array(i)*dp

               else
                  dp = x_ave_layer - dp_layer
                  dp_layer = dp_layer + dp
                  x_aver = x_aver + array(i)*dp
                  exit
               end if
            end do
            x_aver = x_aver/dp_layer
            ic = max(i_beg, i)
         end if
         !print*,"xaver1=",real(x_aver,4),real(dp_layer,4)

         !-- this perturbation is included only for MSE
         if (present(tpert)) x_aver = x_aver + real(c_cp)*maxval(tpert(i_beg:ic)) 

      end if
      if (present(add)) x_aver = x_aver + add

   end subroutine getCloudBc 

   !------------------------------------------------------------------------------------
   subroutine getLcl(t0, pp0, r0, tlcl, plcl, dzlcl)
      !! ## Determine the lift of condensation level (LCL) for air parcel
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Determine the lift of condensation level (LCL) for air parcels 
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'getLcl' ! Subroutine name
   
      !Variables (input, output, inout)
      real, intent(in) :: t0
      real, intent(in) :: pp0
      real, intent(in) :: r0

      real, intent(out) :: tlcl
      real, intent(out) :: plcl
      real, intent(out) :: dzlcl
 
      !Local variables:
      integer :: nitt, ip
      real :: p0k, pi0i, ttth0, ttd, dz, pki, pppi, ti, rvs, e

      !Code:
      !real, external :: td,satvap

      !================
      !-simpler,cheaper method
      ttd = Td(pp0, r0)
      tlcl = ttd - (0.001296*ttd + 0.1963)*(t0 - ttd)
      plcl = pp0*(tlcl/t0)**c_cpor
      dzlcl = 127*(t0 - ttd)
      if (dzlcl .le. 0.) dzlcl = -999.
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
   end subroutine getLcl

   !------------------------------------------------------------------------------------
   function Td(ppp, rs) result(r_td)
      !! ## Determine the dew-point temperature (K)
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! —
      !! **Full description**:
      !!
      !! Determine the dew-point temperature (K)
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !!
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'Td' ! Function name
   
      !Variables (input):
      real, intent(in) :: ppp
      real, intent(in) :: rs
   
      !Local variables:
      real :: r_td !output

      real :: rr, es, esln
      rr = rs + 1e-8
      es = ppp*rr/(.622 + rr)
      esln = log(es)
      r_td = (35.86*esln - 4947.2325)/(esln - 23.6837)

   end function Td

   ! ------------------------------------------------------------------------------------
   subroutine getInversionLayers(cumulus, ierr, psur, po_cup, to_cup, zo_cup, k_inv_layers, dtempdz, itf, ktf, its, ite, kts, kte)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'getInversionLayers' ! subroutine name

      integer, parameter :: p_extralayer = 0 
      !! makes plume top higher
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte
      
      real, intent(in) :: psur(:)
      real, intent(in) :: po_cup(:, :)
      real, intent(in) :: to_cup(:, :)
      real, intent(in) :: zo_cup(:, :)
      
      character(len=*), intent(in)  :: cumulus

      integer, intent(inout):: ierr(:)
      
      integer, intent(out)  :: k_inv_layers(:, :)

      real, intent(out)  :: dtempdz(:, :)

      !Local variables:
      integer:: i, k, ilev, kk, k1, ix, k800, k550, ist
      integer :: local_k_inv_layers(its:ite, kts:kte)
      real :: dzm, delp, first_deriv(kts:kte), sec_deriv(kts:kte), distance(kts:kte)

      !-initialize k_inv_layers as 1 (non-existent layer)_
      k_inv_layers = 1 !integer
      dtempdz = 0.0
      first_deriv = 0.0
      sec_deriv = 0.0
      distance = 0.0
      local_k_inv_layers = 1
      ist = 3

      do i = its, itf
         if (ierr(i) /= 0) cycle
         !- displacement from local surface pressure level
         delp = 1000.-psur(i)

         !- 2nd method
         ! DO k = kts+1,ktf-2
         !dtempdz(i,k)=   ( deriv3(zo_cup(i,k), zo_cup(i,kts:ktf), to_cup(i,kts:ktf), ktf-kts+1, 1,ierr(i)))
           !!! sec_deriv(k)=abs( deriv3(zo_cup(i,k), zo_cup(i,kts:ktf), to_cup(i,kts:ktf), ktf-kts+1, 2))
         !print*,"2=",k,z_cup(i,k),dtempdz(i,k),
         ! ENDDO
         ! if(ierr(i) /= 0) cycle

         !-1st method
         !-  get the 1st derivative
         do k = kts + ist, ktf - ist
            first_deriv(k) = (to_cup(i, k + 1) - to_cup(i, k - 1))/(zo_cup(i, k + 1) - zo_cup(i, k - 1))
         end do
         first_deriv(kts:kts + ist - 1) = first_deriv(kts + ist)
         first_deriv(ktf - ist + 1:kte) = first_deriv(ktf - ist)

         dtempdz(i, :) = first_deriv(:)

         !-  get the abs of the 2nd derivative
         do k = kts + ist + 1, ktf - ist - 1
            sec_deriv(k) = abs((first_deriv(k + 1) - first_deriv(k - 1))/(zo_cup(i, k + 1) - zo_cup(i, k - 1)))
         end do
         sec_deriv(kts:kts + ist) = sec_deriv(kts + ist + 1)
         sec_deriv(ktf - ist:kte) = sec_deriv(ktf - ist - 1)

         ix = 1
         do kk = kts + ist + 2, ktf - ist - 2
            if (sec_deriv(kk) < sec_deriv(kk + 1) .and. sec_deriv(kk) < sec_deriv(kk - 1)) then
               local_k_inv_layers(i, ix) = kk
               ix = ix + 1
            end if
         end do

         !- 2nd criteria
         do k = kts + ist + 2, ktf - ist - 2
            kk = local_k_inv_layers(i, k)
            if (kk == 1) cycle
            if (dtempdz(i, kk) < dtempdz(i, kk - 1) .and. dtempdz(i, kk) < dtempdz(i, kk + 1)) then 
               ! the layer is not a local maximum
               local_k_inv_layers(i, k) = 1
            end if
         end do

      end do

      !- find the locations of inversions around 800 and 550 hPa
      do i = its, itf
         !----------------
         !k_inv_layers(i,mid)=1
         !----------------
         if (ierr(i) /= 0) cycle
         !- displacement from local surface pressure level
         delp = 1000.-psur(i)
         !----------------
         !k_inv_layers(i,mid)=21
         !cycle
         !----------------
         if (trim(cumulus) == 'shallow') then
            !- now find the closest layers of 800 and 550 hPa.
            !- this is for shallow convection k800
            do k = kts, ktf
               distance(k) = abs(po_cup(i, local_k_inv_layers(i, k)) - (750.-delp))
            end do
            k800 = minloc(abs(distance(kts:ktf)), 1)

            if (k800 <= kts .or. k800 >= ktf - 4) then
               k_inv_layers(i, p_shal) = ktf
               !ierr(i)=8
            else
               !-save k800 in the k_inv_layers array
               k_inv_layers(i, p_shal) = local_k_inv_layers(i, k800) + p_extralayer
            end if
            !if(  k_inv_layers(i,shal) <= kts .or. k_inv_layers(i,shal) >= ktf-4) then
            !print*,"SHAL_k_inv_layers=",k_inv_layers(i,shal),ierr(i)
            !ierr(i)=11
            !endif

         elseif (trim(cumulus) == 'mid') then
            !- this is for mid/congestus convection k500
            do k = kts, ktf
               distance(k) = abs(po_cup(i, local_k_inv_layers(i, k)) - (550.-delp))
            end do
            k550 = minloc(abs(distance(kts:ktf)), 1)

            if (k550 <= kts .or. k550 >= ktf - 4) then
               k_inv_layers(i, p_mid) = 1
               ierr(i) = 8
            else
               !-save k550 in the k_inv_layers array
               k_inv_layers(i, p_mid) = local_k_inv_layers(i, k550) + p_extralayer
            end if
            if (k_inv_layers(i, p_mid) <= kts .or. k_inv_layers(i, p_mid) >= ktf - 4) then
               !print*,"MID_k_inv_layers=",k_inv_layers(i,MID),ierr(i)
               ierr(i) = 12
            end if
         else
            k_inv_layers(i, :) = 1
            ierr(i) = 88
         end if
      end do

   contains

      function Deriv3(xx, xi, yi, ni, m, ierr) result(deriv3_out)
         !! ## Evaluate first- or second-order derivatives
         !!
         !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
         !!
         !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
         !!
         !! Date: 2014
         !!
         !! #####Version: 0.1.0
         !!
         !! —
         !! **Full description**:
         !!
         !! Evaluate first- or second-order derivatives
         !! using three-point Lagrange interpolation
         !! written by: Alex Godunov (October 2009)
         !! --------------------------------------------------------------------
         !! input ...
         !! xx    - the abscissa at which the interpolation is to be evaluated
         !! xi()  - the arrays of data abscissas
         !! yi()  - the arrays of data ordinates
         !! ni - size of the arrays xi() and yi()
         !! m  - order of a derivative (1 or 2)
         !! output ...
         !! deriv3_out  - interpolated value
         !!
         !!
         !! ** History**:
         !!
         !! - 
         !! ---
         !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
         !!     
         implicit none
         !Parameters:
         character(len=*), parameter :: p_procedure_name = 'Deriv3' ! Function name

         integer, parameter :: p_n = 3
      
         !Variables (input):
         integer, intent(in):: ni, m
         real, intent(in):: xx
         real, intent(in) :: xi(:)
         real, intent(in) :: yi(:)
         
         !Local variables:
         real :: deriv3_out
         
         real:: x(p_n), f(p_n)
         integer i, j, k, ix
         integer, intent(inout) :: ierr

         ! exit if too high-order derivative was needed,
         if (m > 2) then
            deriv3_out = 0.0
            return
         end if

         ! if x is ouside the xi(1)-xi(ni) interval set deriv3=0.0
         if (xx < xi(1) .or. xx > xi(ni)) then
            deriv3_out = 0.0
            ierr = 8
            !stop "problem with 2nd derivative-deriv3 routine"
            return
         end if

         ! a binary (bisectional) search to find i so that xi(i-1) < x < xi(i)
         i = 1
         j = ni
         do while (j > i + 1)
            k = (i + j)/2
            if (xx < xi(k)) then
               j = k
            else
               i = k
            end if
         end do

         ! shift i that will correspond to n-th order of interpolation
         ! the search point will be in the middle in x_i, x_i+1, x_i+2 ...
         i = i + 1 - p_n/2

         ! check boundaries: if i is ouside of the range [1, ... n] -> shift i
         if (i < 1) i = 1
         if (i + p_n > ni) i = ni - p_n + 1

         !  old output to test i
         !  write(*,100) xx, i
         !  100 format (f10.5, I5)

         ! just wanted to use index i
         ix = i

         ! initialization of f(n) and x(n)
         do i = 1, p_n
            f(i) = yi(ix + i - 1)
            x(i) = xi(ix + i - 1)
         end do

         ! calculate the first-order derivative using Lagrange interpolation
         if (m == 1) then
            deriv3_out = (2.0*xx - (x(2) + x(3)))*f(1)/((x(1) - x(2))*(x(1) - x(3)))
            deriv3_out = deriv3_out+ (2.0*xx - (x(1) + x(3)))*f(2)/((x(2) - x(1))*(x(2) - x(3)))
            deriv3_out = deriv3_out+ (2.0*xx - (x(1) + x(2)))*f(3)/((x(3) - x(1))*(x(3) - x(2)))
            ! calculate the second-order derivative using Lagrange interpolation
         else
            deriv3_out = 2.0*f(1)/((x(1) - x(2))*(x(1) - x(3)))
            deriv3_out = deriv3_out+ 2.0*f(2)/((x(2) - x(1))*(x(2) - x(3)))
            deriv3_out = deriv3_out+ 2.0*f(3)/((x(3) - x(1))*(x(3) - x(2)))
         end if
      end function Deriv3
   
   end subroutine getInversionLayers

   ! ---------------------------------------------------------------------------------------------------
   subroutine setGradsVar(i_in, k_in, nvar, f_in, name1, name2, name3)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'setGradsVar' ! Subroutine name
   
      !Variables (input, output, inout)
      integer, intent(in)    :: i_in, k_in
      
      real, intent(in) :: f_in
      
      character(len=*), intent(in) :: name1
      character(len=*), intent(in) :: name2
      character(len=*), intent(in) :: name3

      integer, intent(inout) :: nvar

      cupout(nvar)%varp(i_in, k_in) = f_in
      cupout(nvar)%varn(1) = name1
      cupout(nvar)%varn(2) = name2
      cupout(nvar)%varn(3) = name3
      nvar = nvar + 1
      if (nvar > p_nvar_grads) stop 'nvar>nvar_grads'

   end subroutine setGradsVar

   ! ------------------------------------------------------------------------------------
   subroutine wrtBinCtl(n, mzp, p2d, cumulus)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'wrtBinCtl' ! subroutine name

      real, parameter :: p_undef = -9.99e33
   
      !Variables (input, output, inout)
      integer, intent(in):: n, mzp

      real, intent(in) :: p2d(:)

      character(len=*), intent(in) :: cumulus
      
      !Local variables:
      integer:: nvartotal, klevgrads(200), jk, int_byte_size, nvar, maxklevgrads
      real   :: real_byte_size
      integer :: nrec = 0
      integer :: rec_size

      maxklevgrads = min(60, mzp)
      runname = '15geos5_'//cumulus
      runlabel = runname

      print *, "writing grads control file:',trim(runname)//'.ctl", ntimes
      call flush (6)
      !
      !number of variables to be written
      nvartotal = 0
      do nvar = 1, p_nvar_grads
         if (cupout(nvar)%varn(1) .ne. "xxxx") nvartotal = nvartotal + 1
         if (cupout(nvar)%varn(3) == "3d") klevgrads(nvar) = maxklevgrads
         if (cupout(nvar)%varn(3) == "2d") klevgrads(nvar) = 1
      end do

      !- binary file
      inquire (iolength=int_byte_size) real_byte_size  ! inquire by output list

      print *, 'opening grads file:', trim(runname)//'.gra'
      rec_size = size(cupout(nvar)%varp, 1)*real_byte_size
      if (ntimes == 1) then
         open (19, file=trim(runname)//'.gra', form='unformatted', &
               access='direct', status='replace', recl=rec_size)
      else
         open (19, file=trim(runname)//'.gra', form='unformatted', &
               access='direct', status='old', recl=rec_size)
      end if

      do nvar = 1, p_nvar_grads
         if (cupout(nvar)%varn(1) .ne. "xxxx") then
            do jk = 1, klevgrads(nvar)
               nrec = nrec + 1
               !write(19)          real((cupout(nvar)%varp(:,jk)),4)
               write (19, rec=nrec) real((cupout(nvar)%varp(:, jk)), 4)
            end do
         end if
      end do
      close (19)
      !-setting vertical dimension '0' for 2d var
      where (klevgrads == 1) klevgrads = 0
      !- ctl file
      open (20, file=trim(runname)//'.ctl', status='unknown')
      write (20, 2001) '^'//trim(runname)//'.gra'
      write (20, 2002) 'undef -9.99e33'
      write (20, 2002) 'options sequential byteswapped' ! zrev'
      write (20, 2002) 'title '//trim(runlabel)
      write (20, 2003) 1, 0., 1. ! units m/km
      write (20, 2004) n, 1., 1.
      write (20, 2005) maxklevgrads, (p2d(jk), jk=1, maxklevgrads)
      write (20, 2006) ntimes, '00:00Z24JAN1999', '10mn'
      write (20, 2007) nvartotal
      do nvar = 1, p_nvar_grads
         if (cupout(nvar)%varn(1) .ne. "xxxx") then
            !
            write (20, 2008) cupout(nvar)%varn(1) (1:len_trim(cupout(nvar)%varn(1))), klevgrads(nvar) &
               , cupout(nvar)%varn(2) (1:len_trim(cupout(nvar)%varn(2)))
         end if
      end do
      write (20, 2002) 'endvars'
      close (20)

2001  format('dset ', a)
2002  format(a)
2003  format('xdef ', i4, ' linear ', 2f15.3)
2004  format('ydef ', i4, ' linear ', 2f15.3)
2005  format('zdef ', i4, ' levels ', 60f8.3)
2006  format('tdef ', i4, ' linear ', 2a15)
2007  format('vars ', i4)
2008  format(a10, i4, ' 99 ', a40)!'[',a8,']')
2055  format(60f7.0)
133   format(1x, F7.0)

   end subroutine wrtBinCtl

   ! ------------------------------------------------------------------------------------
   subroutine cupCloudLimits(name, ierrc, ierr, cap_inc, cap_max_in, heo_cup, heso_cup, qo_cup, qeso_cup, po &
                           , po_cup, z_cup, heo, hkbo, qo, qeso, entr_rate_2d, hcot, k22, kbmax, klcl, kbcon &
                           , ktop, depth_neg_buoy, frh, tpert, start_level_, use_excess, zqexec, ztexec &
                           , x_add_buoy, xland, itf, ktf, its, ite, kts, kte)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'cupCloudLimits' ! subroutine name

      real, parameter :: p_frh_crit_O = 0.7
      real, parameter :: p_frh_crit_L = 0.7  !--- test 0.5
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte, use_excess

      integer, intent(in) :: kbmax(:)
      integer, intent(in) :: start_level_(:)

      real, intent(in) :: heo_cup(:, :)
      real, intent(in) :: heso_cup(:, :)
      real, intent(in) :: po_cup(:, :)
      real, intent(in) :: z_cup(:, :)
      real, intent(in) :: heo(:, :)
      real, intent(in) :: qo_cup(:, :)
      real, intent(in) :: qeso_cup(:, :)
      real, intent(in) :: po(:, :)
      real, intent(in) :: qo(:, :)
      real, intent(in) :: qeso(:, :)
      real, intent(in) :: tpert(:, :)
      real, intent(in) :: cap_max_in(:)
      real, intent(in) :: cap_inc(:)
      real, intent(in) :: xland(:)
      real, intent(in) :: zqexec(:)
      real, intent(in) :: ztexec(:)
      real, intent(in) :: x_add_buoy(:)
      real, intent(in) :: entr_rate_2d(:, :)

      character(len=*), intent(in) :: name

      integer, intent(inout) :: kbcon(:)
      integer, intent(inout) :: ierr(:)
      integer, intent(inout) :: ktop(:)
      integer, intent(inout) :: klcl(:)
      integer, intent(inout) :: k22(:)

      real, intent(inout) :: hkbo(:)
      real, intent(inout) :: depth_neg_buoy(:)
      real, intent(inout) :: frh(:)
      real, intent(inout) :: hcot(:, :)

      character(len=128), intent(inout) :: ierrc(:)

      !Local variables:
      integer :: i, k, k1, k2, kfinalzu
      integer, dimension(its:ite) :: start_level
      real   :: delz_oversh 
      !! height of cloud overshoot is 10% higher than the LNB.
      !! Typically it can 2 - 2.5km higher, but it depends on
      !!the severity of the thunderstorm.
      real, dimension(its:ite) :: cap_max
      real :: plus, hetest, dz, dbythresh, denom &
            , dzh, del_cap_max, fx, x_add, Z_overshoot, frh_crit
      real, dimension(kts:kte) ::   dby

      delz_oversh = OVERSHOOT
      hcot = 0.0
      dby = 0.0
      start_level = 0
      cap_max(:) = cap_max_in(:)

      do i = its, itf
         if (ierr(i) /= 0) cycle
         start_level(i) = start_level_(i)
         do k = kts, start_level(i)
            hcot(i, k) = hkbo(i) ! assumed no entraiment between these layers
         end do
      end do

      !--- determine the level of convective cloud base  - kbcon
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
                  x_add    = (c_alvl*zqexec(i)+c_cp*ztexec(i)) + x_add_buoy(i)
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
               frh_crit =p_frh_crit_O*xland(i) + p_frh_crit_L*(1.-xland(i))

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
            x_add = (c_alvl*zqexec(i)+c_cp*ztexec(i)) +  x_add_buoy(i)
            call getCloudBc(name,kts,kte,ktf,xland(i),po(i,kts:kte),heo_cup (i,kts:kte),hkbo (i),k22(i),x_add,Tpert(i,kts:kte))
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


      !--- determine the level of neutral buoyancy - ktop
      do i = its, itf
         ktop(i) = ktf - 1
         if (ierr(i) /= 0) cycle
         !~ dby(:)=0.0

         start_level(i) = kbcon(i)

         do k = start_level(i) + 1, ktf - 1
            dz = z_cup(i, k) - z_cup(i, k - 1)
            denom = 1.+0.5*entr_rate_2d(i, k - 1)*dz
            if (denom == 0.) then
               hcot(i, k) = hcot(i, k - 1)
            else
               hcot(i, k) = ((1.-0.5*entr_rate_2d(i, k - 1)*dz)*hcot(i, k - 1) + entr_rate_2d(i, k - 1)*dz*heo(i, k - 1))/denom
            end if
         end do
         do k = start_level(i) + 1, ktf - 1

            if (hcot(i, k) < heso_cup(i, k)) then
               ktop(i) = k - 1
               exit
            end if
         end do
         if (ktop(i) .le. kbcon(i) + 1) ierr(i) = 41

         !----------------
         if (OVERSHOOT > 1.e-6 .and. ierr(i) == 0) then
            Z_overshoot = (1.+delz_oversh)*z_cup(i, ktop(i))
            do k = ktop(i), ktf - 2
               if (Z_overshoot < z_cup(i, k)) then
                  ktop(i) = min(k - 1, ktf - 2)
                  exit
               end if
            end do
         end if
      end do
   end subroutine cupCloudLimits

   ! ------------------------------------------------------------------------------------
   subroutine getBuoyancy(itf, ktf, its, ite, kts, kte, ierr, klcl, kbcon, ktop, hc, he_cup, hes_cup, dby, z_cup)
      !! ## Detertime the buoyancy of updraft air parcels
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Detertime the buoyancy of updraft air parcels
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'getBuoyancy' ! Nome da subrotina
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte
      integer, intent(in) :: ierr(:)
      integer, intent(in) :: klcl(:)
      integer, intent(in) :: kbcon(:)
      integer, intent(in) :: ktop(:)

      real, intent(in) :: hc(:, :)
      real, intent(in) :: he_cup(:, :)
      real, intent(in) :: hes_cup(:, :)
      real, intent(in) :: z_cup(:, :)

      real, intent(out) :: dby(:, :)

      !Local variables:
      integer :: i, k
   
      do i = its, itf
         dby(i, :) = 0.
         if (ierr(i) /= 0) cycle
         do k = kts, klcl(i)
            dby(i, k) = hc(i, k) - he_cup(i, k)
         end do
         do k = klcl(i) + 1, ktop(i) + 1
            dby(i, k) = hc(i, k) - hes_cup(i, k)
         end do
      end do
      
   end subroutine getBuoyancy

   ! ------------------------------------------------------------------------------------
   subroutine cupUpVVel(vvel2d, vvel1d, zws, entr_rate_2d, cd, z, z_cup, zu, dby, GAMMA_CUP, t_cup,  tempco, qco, qrco, qo &
                      , start_level, klcl, kbcon, ktop, ierr, itf, ktf, its, ite, kts, kte, wlpool, wlpool_bcon, task)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'cupUpVVel'

      integer, parameter :: p_n_smooth = 1
      real, parameter :: p_ctea = 1./3.
      real, parameter :: p_cteb = 2.
      real, parameter :: p_visc = 2000.
      real, parameter :: p_eps = 0.622
      real, parameter :: p_f = 2., p_c_d = 0.506, p_gam = 0.5, p_beta = 1.875 !,ftun1=0.5, ftun2=0.8
      logical, parameter :: p_smooth = .true.

      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte, task

      integer, intent(in) :: klcl(:)
      integer, intent(in) :: kbcon(:)
      integer, intent(in) :: ktop(:)
      integer, intent(in) :: start_level(:)
   
      real, intent(in) :: z(:, :)
      real, intent(in) :: z_cup(:, :)
      real, intent(in) :: zu(:, :)
      real, intent(in) :: gamma_cup(:, :)
      real, intent(in) :: t_cup(:, :)
      real, intent(in) :: dby(:, :)
      real, intent(in) :: entr_rate_2d(:, :)
      real, intent(in) :: cd(:, :)
      real, intent(in) :: tempco(:, :)
      real, intent(in) :: qco(:, :)
      real, intent(in) :: qrco(:, :)
      real, intent(in) :: qo(:, :)
      real, intent(in) :: zws(:)

      integer, intent(inout) :: ierr(:)

      real, intent(inout) :: wlpool(:)
      real, intent(inout) :: wlpool_bcon(:)

      real, intent(out) :: vvel2d(:, :)
      real, intent(out) :: vvel1d(:)

      !Local variables:
      integer :: i, k, k1, nvs
      real :: dz, bu, dw2, dw1, kx, dz1m, tv, tve, vs, ftun1, ftun2, ke

      ftun1 = 0.25
      ftun2 = 1.

      if (task == 1) then
         do i = its, itf
            !-- initialize arrays to zero.
            vvel1d(i) = 0.0
            vvel2d(i, :) = 0.0

            if (ierr(i) /= 0) cycle
            vvel2d(i, kts:kbcon(i)) = max(1., max(wlpool_bcon(i)**2, zws(i)**2))

            loop0: do k = kbcon(i), ktop(i)
               dz = z_cup(i, k + 1) - z_cup(i, k)
               tve = 0.5*(t_cup(i, k)*(1.+(qo(i, k)/p_eps)/(1.+qo(i, k))) &
                        +  t_cup(i, k + 1)*(1.+(qo(i, k + 1)/p_eps)/(1.+qo(i, k + 1))))
               tv = 0.5*(tempco(i, k)*(1.+(qco(i, k)/p_eps)/(1.+qco(i, k))) &
                        + tempco(i, k + 1)*(1.+(qco(i, k + 1)/p_eps)/(1.+qco(i, k + 1))))
               bu = c_grav*((tv - tve)/tve - ftun2*0.50*(qrco(i, k + 1) + qrco(i, k)))
               dw1 = 2./(p_f*(1.+p_gam))*bu*dz
               kx = (1.+p_beta*p_c_d)*max(entr_rate_2d(i, k), cd(i, k))*dz*ftun1
               dw2 = (vvel2d(i, k)) - 2.*kx*(vvel2d(i, k))
               vvel2d(i, k + 1) = (dw1 + dw2)/(1.+kx)

               if (vvel2d(i, k + 1) < 0.) then
                  vvel2d(i, k + 1) = 0.5*vvel2d(i, k)
               end if
            end do loop0
         end do
         if (p_smooth) then
            do i = its, itf
               if (ierr(i) /= 0) cycle
               do k = kts, ktop(i) + 1
                  vs = 0.
                  dz1m = 0.
                  do k1 = max(k - p_n_smooth, kts), min(k + p_n_smooth, ktf)
                     dz = z_cup(i, k1 + 1) - z_cup(i, k1)
                     vs = vs + dz*vvel2d(i, k1)
                     dz1m = dz1m + dz
                  end do
                  vvel2d(i, k) = vs/(1.e-16 + dz1m)
                  !if(k>ktop(i)-3)print*,"v2=",k,ktop(i),sqrt(vvel2d(i,k)),sqrt(vvel2d(i,ktop(i)))
               end do
            end do
         end if

         !-- convert to vertical velocity
         do i = its, itf
            if (ierr(i) /= 0) cycle
            vvel2d(i, :) = sqrt(max(0.1, vvel2d(i, :)))

            if (maxval(vvel2d(i, :)) < 1.0) then
               ierr(i) = 54
               !  print*,"ierr=54",maxval(vvel2d(i,:))
            end if

            !-- sanity check
            where (vvel2d(i, :) < 1.) vvel2d(i, :) = 1.
            where (vvel2d(i, :) > 20.) vvel2d(i, :) = 20.
            vvel2d(i, ktop(i) + 1:kte) = 0.1

            !-- get the column average vert velocity
            do k = kbcon(i), ktop(i)
               dz = z_cup(i, k + 1) - z_cup(i, k)
               vvel1d(i) = vvel1d(i) + vvel2d(i, k)*dz
               !print*,"w=",k,z_cup(i,k),vvel2d(i,k)
            end do
            vvel1d(i) = vvel1d(i)/(z_cup(i, ktop(i) + 1) - z_cup(i, kbcon(i)) + 1.e-16)
            vvel1d(i) = max(1., vvel1d(i))
         end do
      else
         do i = its, itf
            if (ierr(i) /= 0) cycle
            ke = wlpool(i)**2

            loop1: do k = start_level(i), kbcon(i)

               dz = z_cup(i, k + 1) - z_cup(i, k)

               tve = 0.5*(t_cup(i, k)*(1.+(qo(i, k)/p_eps)/(1.+qo(i, k))) + t_cup(i, k + 1)*(1.+(qo(i, k + 1)/p_eps) &
                   / (1.+qo(i, k + 1))))

               tv = 0.5*(tempco(i, k)*(1.+(qco(i, k)/p_eps)/(1.+qco(i, k))) + tempco(i, k + 1)*(1.+(qco(i, k + 1)/p_eps) &
                  / (1.+qco(i, k + 1))))
               bu = c_grav*((tv - tve)/tve - ftun2*0.50*(qrco(i, k + 1) + qrco(i, k)))
               dw1 = 2./(p_f*(1.+p_gam))*bu*dz
               kx = (1.+p_beta*p_c_d)*max(entr_rate_2d(i, k), cd(i, k))*dz*ftun1
               dw2 = ke - 2.*kx*ke
               ke = max(0., (dw1 + dw2)/(1.+kx))

!            vvel2d(i,k)=sqrt(ke)
            end do loop1
            wlpool_bcon(i) = sqrt(ke)
            !print*,"wlpool=",wlpool(i),sqrt (ke)
         end do
      end if

   end subroutine cupUpVVel

   !------------------------------------------------------------------------------------
   subroutine cupOutputEns3d(cumulus, xff_shal, xff_mid, xf_ens, ierr, dellat, dellaq, dellaqc, outtem, outq, outqc &
                           , zu, pre, pw, xmb, ktop, nx, nx2, ierr2, ierr3, pr_ens, maxens3, ensdim, sig, xland1 &
                           , ichoice, ipr, jpr, itf, ktf, its, ite, kts, kte, xf_dicycle, outu, outv, dellu, dellv &
                           , dtime, po_cup, kbcon, dellabuoy, outbuoy, dellampqi, outmpqi, dellampql, outmpql &
                           , dellampcf, outmpcf, nmp, rh_dicycle_fct, xf_coldpool, wlpool_bcon)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'cupOutputEns3d' ! Subroutine name
   
      !Variables (input, output, inout)
      integer, intent(in) :: ichoice, ipr, jpr, itf, ktf, its, ite, kts, kte
      integer, intent(in) :: ensdim, nx, nx2, maxens3, nmp   

      integer, intent(in) :: ktop(:)
      integer, intent(in) :: kbcon(:)

      real, intent(in) :: zu(:, :)
      real, intent(in) :: po_cup(:, :)
      real, intent(in) :: sig(:)
      real, intent(in) :: rh_dicycle_fct(:)
      real, intent(in) :: wlpool_bcon(:)
      real, intent(in) :: xff_mid(:, :)
      real, intent(in) :: dellampqi(: , :, :)
      real, intent(in) :: dellampql(: , :, :)
      real, intent(in) :: dellampcf(: , :, :)
      real, intent(in) :: xff_shal(:, :)
      real, intent(in) :: dtime

      character(len=*), intent(in) :: cumulus

      integer, intent(inout) :: ierr(:)
      !! ierr error value, maybe modified in this routine
      integer, intent(inout) :: ierr2(:)
      integer, intent(inout) :: ierr3(:)

      real, intent(inout) :: xland1(:)
      real, intent(inout) :: dellat(:, :)
      !! change of temperature per unit mass flux of cloud ensemble
      real, intent(inout) :: dellaqc(:, :)
      !! change of qc per unit mass flux of cloud ensemble
      real, intent(inout) :: dellaq(:, :)
      !! change of q per unit mass flux of cloud ensemble
      real, intent(inout) :: pw(:, :)
      !! epsilon*pd (ensemble dependent)
      real, intent(inout) :: dellu(:, :)
      real, intent(inout) :: dellv(:, :)
      real, intent(inout) :: dellabuoy(:, :)
      real, intent(inout) :: xf_ens(:, :)
      !! ensemble mass fluxes
      real, intent(inout) :: pr_ens(:, :)
      !! precipitation ensembles
      real, intent(inout) :: xf_dicycle(:)
      real, intent(inout) :: xf_coldpool(:)

      real, intent(out) :: outtem(:, :)
      !! output temp tendency (per s)
      real, intent(out) :: outq(:, :)
      !! output q tendency (per s)
      real, intent(out) :: outqc(:, :)
      !! output qc tendency (per s)
      real, intent(out) :: outu(:, :)
      real, intent(out) :: outv(:, :)
      real, intent(out) :: outbuoy(:, :)
      real, intent(out) :: outmpqi(: , :, :)
      real, intent(out) :: outmpql(: , :, :)
      real, intent(out) :: outmpcf(: , :, :)
      real, intent(out) :: pre(:)
      !! output precip
      real, intent(out) :: xmb(:)
      !! total base mass flux

      !Local variables:
      integer :: i, k, n, ncount, zmax, kk, kqmx, ktmx
      real :: outtes, ddtes, dtt, dtq, dtqc, dtpw, prerate, fixouts, dp, xfix_q, xfix_t
      real ::  dtts, dtqs, fsum, rcount
      real, dimension(its:ite) :: xmb_ave, xmbmax
      real, dimension(kts:kte, 8) :: tend2d
      real, dimension(8) :: tend1d
      real, dimension(its:ite, 8) :: check_cons_I, check_cons_F
      !
      !-- initialization
      pre     = 0.
      xmb     = 0.
      xmb_ave = 0.
      outtem  = 0.
      outq    = 0.
      outqc   = 0.
      outu    = 0.
      outv    = 0.
      outbuoy = 0.

      do i = its, itf
         if (ierr(i) .eq. 0) then
            do n = 1, maxens3
               if (pr_ens(i, n) .le. 0.) then
                  xf_ens(i, n) = 0.
               end if
            end do
         end if
      end do

      !--- calculate ensemble average mass fluxes
      if (trim(cumulus) == 'deep') then
         do i = its, itf
            if (ierr(i) .eq. 0) then
               k = 0
               xmb_ave(i) = 0.
               do n = 1, maxens3
                  k = k + 1
                  xmb_ave(i) = xmb_ave(i) + xf_ens(i, n)
               end do
               !- 'ensemble' average mass flux
               xmb_ave(i) = xmb_ave(i)/float(k)
            end if
         end do

      !--- mid (congestus type) convection
      elseif (trim(cumulus) == 'mid') then
         if (ichoice .le. 3) then
            do i = its, itf
               if (ierr(i) /= 0) cycle
               if (ichoice == 0) then
                  xmb_ave(i) = 0.3333*(xff_mid(i, 1) + xff_mid(i, 2) + xff_mid(i, 3))
               else
                  xmb_ave(i) = xff_mid(i, ichoice)
               end if
            end do
         else
            stop 'For mid ichoice must be 0,1,2,3'
         end if

      !--- shallow  convection
      elseif (trim(cumulus) == 'shallow') then
         do i = its, itf
            if (ierr(i) /= 0) cycle

            if (ichoice > 0) then
               xmb_ave(i) = xff_shal(i, ichoice)
            else
               fsum = 0.
               xmb_ave(i) = 0.
               do k = 1, p_shall_closures
                  !- heat engine closure is not working properly
                  !- turning it off for now.
                  if (k .ge. 4 .and. k .le. 6) cycle
                  xmb_ave(i) = xmb_ave(i) + xff_shal(i, k)
                  fsum = fsum + 1.
               end do
               !- ensemble average of mass flux
               xmb_ave(i) = xmb_ave(i)/fsum
            end if
         end do
      end if
      !- apply the mean tropospheric RH control on diurnal cycle (Tian GRL 2022)
      if (trim(cumulus) == 'deep' .and. RH_DICYCLE == 1) then
         do i = its, itf
            if (ierr(i) /= 0) cycle
            xf_dicycle(i) = xf_dicycle(i)*rh_dicycle_fct(i)
         end do
      end if

!if(trim(cumulus) == 'deep' ) then
!  do i=its,itf
!!!     if(ierr(i) /= 0) cycle
!      print*,'xmbs',xmb_ave(i),xf_dicycle(i),xf_coldpool(i),wlpool_bcon(i)
!      call flush(6)
!  enddo
!endif

      !- determine the updraft mass flux, do not allow negative values and apply the diurnal cycle closure
      do i = its, itf
         if (ierr(i) /= 0) cycle
         !- mass flux of updradt at cloud base
         xmb(i) = xmb_ave(i)

         !- add kinetic energy at the gust front of the cold pools
         xmb(i) = xmb(i) + xf_coldpool(i)

         !- diurnal cycle closure
         xmb(i) = xmb(i) - xf_dicycle(i)
         if (xmb(i) .le. 0.) then
            ierr(i) = 13
            xmb(i) = 0.
         end if
      end do
      !-apply the scale-dependence Arakawa's approach
      do i = its, itf
         if (ierr(i) /= 0) cycle
         !- scale dependence
         xmb(i) = sig(i)*xmb(i)

         !- apply the adjust factor for tunning
         !xmb(i) = FADJ_MASSFLX * xmb(i)

         if (xmb(i) == 0.) ierr(i) = 14
         if (xmb(i) > 100.) ierr(i) = 15
      end do

      !--- sanity check for mass flux
      !
      do i = its, itf
         if (ierr(i) /= 0) cycle
         xmbmax(i) = 100.*(po_cup(i, kbcon(i)) - po_cup(i, kbcon(i) + 1))/(c_grav*dtime)
         xmb(i) = min(xmb(i), xmbmax(i))
      end do

      !--- check outtem and and outq for unrealistic values
      !
      !--- criteria: if abs (dT/dt or dQ/dt) > MAX_TQ_TEND K/day => fix xmb
      ! if (MAX_TQ_TEND < -1.e-2) then
      !    do i = its, itf
      !       if (ierr(i) /= 0) cycle
      !       fixouts = xmb(i)*86400.*max(maxval(abs(dellat(i, kts:ktop(i)))), (real(c_alvl)/real(c_cp))*maxval(abs &
      !               ( dellaq(i, kts:ktop(i)))))

      !       if (fixouts > abs(MAX_TQ_TEND)) then ! K/day
      !          fixouts = abs(MAX_TQ_TEND)/(fixouts)
      !          xmb(i) = xmb(i)*fixouts
      !          xf_ens(i, :) = xf_ens(i, :)*fixouts
      !       end if
      !    end do
      ! end if
      !--- criteria: if abs (dT/dt or dQ/dt) > 100 K/day => fix dT/dt, dQ/dt and xmb
      if (MAX_TQ_TEND > 1.e-2) then
         do i = its, itf

            if (ierr(i) /= 0) cycle
            tend1d = 0.
            do k = kts, ktop(i)
               dp = (po_cup(i, k) - po_cup(i, k + 1))
               tend1d(1) = tend1d(1) + dp*xmb(i)*86400.*(dellat(i, k))

               if (xmb(i)*86400.*abs(dellat(i, k)) > MAX_TQ_TEND) dellat(i, k) = MAX_TQ_TEND/(xmb(i)*86400)*sign(1., dellat(i, k))

               tend1d(2) = tend1d(2) + dp*xmb(i)*86400.*(dellat(i, k))
            end do

            do k = kts, ktop(i)
               dp = (po_cup(i, k) - po_cup(i, k + 1))
               tend1d(3) = tend1d(3) + dp*xmb(i)*86400.*(dellaq(i, k))*(real(c_alvl)/real(c_cp))

               if (xmb(i)*86400.*abs(dellaq(i, k))*(real(c_alvl)/real(c_cp)) > MAX_TQ_TEND) dellaq(i, k) = MAX_TQ_TEND &
                  / (xmb(i)*86400*(real(c_alvl)/real(c_cp)))*sign(1., dellaq(i, k))

               tend1d(4) = tend1d(4) + dp*xmb(i)*86400.*(dellaq(i, k))*(real(c_alvl)/real(c_cp))
            end do
            xfix_t = tend1d(1)/(1.e-6 + tend1d(2))
            xfix_q = tend1d(3)/(1.e-6 + tend1d(4))

            xmb(i) = xmb(i)/max(1., max(xfix_q, xfix_t))
            !   print*,"tend",
         end do
      end if
      !
      !-- now do feedback
      !
      do i = its, itf
         if (ierr(i) /= 0) cycle
         do k = kts, ktop(i)
            pre(i) = pre(i) + pw(i, k)*xmb(i)

            outtem(i, k) = dellat(i, k)*xmb(i)
            outq(i, k) = dellaq(i, k)*xmb(i)
            outqc(i, k) = dellaqc(i, k)*xmb(i)
            outu(i, k) = dellu(i, k)*xmb(i)
            outv(i, k) = dellv(i, k)*xmb(i)
            outbuoy(i, k) = dellabuoy(i, k)*xmb(i)
         end do
         xf_ens(i, :) = sig(i)*xf_ens(i, :)

         if (APPLY_SUB_MP == 1) then
            do k = kts, ktop(i)
               outmpqi(:, i, k) = dellampqi(:, i, k)*xmb(i)
               outmpql(:, i, k) = dellampql(:, i, k)*xmb(i)
               outmpcf(:, i, k) = dellampcf(:, i, k)*xmb(i)
            end do
            outmpqi(:, i, ktop(i):ktf) = 0.
            outmpql(:, i, ktop(i):ktf) = 0.
            outmpcf(:, i, ktop(i):ktf) = 0.
         end if
      end do
      !
   end subroutine cupOutputEns3d

   !------------------------------------------------------------------------------------
   subroutine cupForcingEns3d(itf, ktf, its, ite, kts, kte, ens4, ensdim, ichoice, maxens, maxens2, maxens3 &
                                 , ierr, ierr2, ierr3, k22, kbcon, ktop, xland, aa0, aa1, xaa0, mbdt, dtime &
                                 , xf_ens, mconv, qo, p_cup, omeg, zd, zu, pr_ens, edt, tau_ecmwf, aa1_bl &
                                 , xf_dicycle, xk_x, alpha_adv, Q_adv, aa1_radpbl, aa1_adv, wlpool, xf_coldpool)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'cupForcingEns3d' ! subroutine name
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte, ens4, ensdim, maxens, maxens2, maxens3

      integer, intent(in) :: k22(:)
      integer, intent(in) :: kbcon(:)
      integer, intent(in) :: ktop(:)
      integer, intent(in) :: ichoice
      !! flag if only want one closure (usually set to zero!)

      ! massfln = downdraft mass flux ensembles used in next timestep
      ! dir     = "storm motion"
      ! iact_gr_old = flag to tell where convection was active
      ! kbcon       = LFC of parcel from k22
      ! k22         = updraft originating level
      ! name        = deep or shallow convection flag
      real, intent(in) :: aa1_bl(:)
      real, intent(in) :: tau_ecmwf(:)
      real, intent(in) :: alpha_adv(:)
      real, intent(in) :: Q_adv(:)
      real, intent(in) :: aa1_radpbl(:)
      real, intent(in) :: aa1_adv(:)
      real, intent(in) :: wlpool(:)
      real, intent(in) :: zd(:, :)
      !! downdraft normalized mass flux
      real, intent(in) :: zu(:, :)
      !! updraft normalized mass flux
      real, intent(in) :: p_cup(:, :)
      real, intent(in) :: qo(:, :)
      real, intent(in) :: omeg(:, :, :)
      !! omega from large scale model
      real, intent(in) :: xaa0(:)
      !! cloud work function with cloud effects
      real, intent(in) :: aa1(:)
      !! cloud work function with forcing effects
      real, intent(in) :: edt(:)
      !! epsilon
      real, intent(in) :: xland(:)
      real, intent(in) :: aa0(:)
      !! cloud work function without forcing effects
      real, intent(in) :: mbdt(:)
      !! arbitrary numerical parameter
      real, intent(in) :: dtime
      !! dt over which forcing is applied

      integer, intent(inout) :: ierr(:)
      !! ierr error value, maybe modified in this routine
      integer, intent(inout) :: ierr2(:)
      integer, intent(inout) :: ierr3(:)

      real, intent(inout) :: pr_ens(:, :)
      !! precipitation ensemble
      real, intent(inout) :: mconv(:)
      !! moisture convergence from large scale model
      real, intent(inout) :: xf_dicycle(:)
      real, intent(inout) :: xk_x(:)
      real, intent(inout) :: xf_coldpool(:)
      
      real, intent(out) :: xf_ens(:, :)
      !! mass flux ensembles

      !Local variables:
      integer :: i, k, nall, n, ne, nens, nens3, kk
      real :: a1, a_ave, xff0, xomg
      real :: betajb, ke
      real  :: xff_dicycle
      real, dimension(1:maxens3) :: xff_ens3
      real, dimension(its:ite) :: xk
      real, dimension(its:ite) :: ens_adj!,xmbmax

      ens_adj(:) = 1.

      ! large scale forcing
      do i = its, itf
         xf_ens(i, 1:16) = 0.
         if (ierr(i) /= 0) cycle

         xff0 = (aa1(I) - aa0(I))/dtime
         !-- default
         xff_ens3(1) = max(0., (aa1(I) - aa0(I))/dtime)

         xff_ens3(2) = xff_ens3(1)
         xff_ens3(3) = xff_ens3(1)
         xff_ens3(16) = xff_ens3(1)
         !
         !--- more like Brown (1979), or Frank-Cohen (199?)
         !--- omeg is in Pa/s
         xomg = 0.
         kk = 0
         xff_ens3(4) = 0.
         do k = max(kts, kbcon(i) - 1), kbcon(i) + 1
            !-  betajb=(zu(i,k)-edt(i)*zd(i,k))
            betajb = 1.
            !if(betajb .gt. 0.)then
            xomg = xomg - omeg(i, k, 1)/c_grav/betajb
            kk = kk + 1
            !endif
         end do
         if (kk .gt. 0) xff_ens3(4) = xomg/float(kk) ! kg[air]/m^3 * m/s
         xff_ens3(4) = max(0.0, xff_ens3(4))
         xff_ens3(5) = xff_ens3(4)
         xff_ens3(6) = xff_ens3(4)
         xff_ens3(14) = xff_ens3(4)
         !
         !--- more like Krishnamurti et al.;
         !
         !mconv(i) = 0.
         !do k=k22(i),ktop(i)
         !    mconv(i)=mconv(i)+omeg(i,k,1)*(qo(i,k+1)-qo(i,k))/g
         !enddo
         !- 2nd option (assuming that omeg(ktop)*q(ktop)<< omeg(kbcon)*q(kbcon))
         mconv(i) = -omeg(i, kbcon(i), 1)*qo(i, kbcon(i))/c_grav ! (kg[air]/m^3)*m/s*kg[water]/kg[air]

         mconv(i) = max(0., mconv(i))
         xff_ens3(7) = mconv(i)
         xff_ens3(8) = xff_ens3(7)
         xff_ens3(9) = xff_ens3(7)
         xff_ens3(15) = xff_ens3(7)
         !
         !---- more like  Betchold et al (2014). Note that AA1 already includes the forcings tendencies
         xff_ens3(10) = aa1(i)/tau_ecmwf(i)

         xff_ens3(11) = xff_ens3(10)
         xff_ens3(12) = xff_ens3(10)
         xff_ens3(13) = xff_ens3(10)

         !
         if (ichoice == 0) then
            if (xff0 < 0.) then
               xff_ens3(1) = 0.
               xff_ens3(2) = 0.
               xff_ens3(3) = 0.
               xff_ens3(16) = 0.

               xff_ens3(10) = 0.
               xff_ens3(11) = 0.
               xff_ens3(12) = 0.
               xff_ens3(13) = 0.
            end if
         end if

         xk(i) = (xaa0(I) - (aa1(I)))/mbdt(i)
         if (xk(i) .le. 0. .and. xk(i) .gt. -0.1*mbdt(i)) xk(i) = -0.1*mbdt(i)
         if (xk(i) .gt. 0. .and. xk(i) .lt. 1.e-2) xk(i) = 1.e-2
         !
         !---  over water, enfor!e small cap for some of the closures
         !
         !if(xland(i).lt.0.1)then
         !   if(ierr2(i).gt.0.or.ierr3(i).gt.0)then
         !      xff_ens3(1:16) = ens_adj(i)*xff_ens3(1:16)
         !   endif
         !endif
         !
         !--- special treatment for stability closures
         !
         if (xk(i) .lt. 0.) then
            if (xff_ens3(1) .gt. 0.) xf_ens(i, 1) = max(0., -xff_ens3(1)/xk(i))
            if (xff_ens3(2) .gt. 0.) xf_ens(i, 2) = max(0., -xff_ens3(2)/xk(i))
            if (xff_ens3(3) .gt. 0.) xf_ens(i, 3) = max(0., -xff_ens3(3)/xk(i))
            if (xff_ens3(16) .gt. 0.) xf_ens(i, 16) = max(0., -xff_ens3(16)/xk(i))
         else
            xff_ens3(1) = 0.
            xff_ens3(2) = 0.
            xff_ens3(3) = 0.
            xff_ens3(16) = 0.
         end if

         xf_ens(i, 4) = max(0., xff_ens3(4))
         xf_ens(i, 5) = max(0., xff_ens3(5))
         xf_ens(i, 6) = max(0., xff_ens3(6))
         xf_ens(i, 14) = max(0., xff_ens3(14))

         a1 = max(1.e-3, pr_ens(i, 7))
         xf_ens(i, 7) = max(0., xff_ens3(7)/a1)
         a1 = max(1.e-3, pr_ens(i, 8))
         xf_ens(i, 8) = max(0., xff_ens3(8)/a1)
         a1 = max(1.e-3, pr_ens(i, 9))
         xf_ens(i, 9) = max(0., xff_ens3(9)/a1)
         a1 = max(1.e-3, pr_ens(i, 15))
         xf_ens(i, 15) = max(0., xff_ens3(15)/a1)
         if (xk(i) .lt. 0.) then
            xf_ens(i, 10) = max(0., -xff_ens3(10)/xk(i))
            xf_ens(i, 11) = max(0., -xff_ens3(11)/xk(i))
            xf_ens(i, 12) = max(0., -xff_ens3(12)/xk(i))
            xf_ens(i, 13) = max(0., -xff_ens3(13)/xk(i))
         else
            xf_ens(i, 10) = 0.
            xf_ens(i, 11) = 0.
            xf_ens(i, 12) = 0.
            xf_ens(i, 13) = 0.
         end if

         if (ichoice .ge. 1) then
            xf_ens(i, 1:16) = xf_ens(i, ichoice)
         end if

         !---special combination for 'ensemble closure':
         !---over the land, only applies closures 1 and 10.
         !if(ichoice == 0 .and. xland(i) < 0.1)then
         !  xf_ens(i,1:16) =0.5*(xf_ens(i,10)+xf_ens(i,1))
         !endif
         !------------------------------------
      end do
      !-
      !- diurnal cycle mass flux closure
      !-
      if (DICYCLE == 1 .or. DICYCLE == 2) then

         do i = its, itf
            xf_dicycle(i) = 0.
!           if(ierr(i) /=  0 .or. p_cup(i,kbcon(i))< 950. )cycle
            if (ierr(i) /= 0) cycle

            !--- Bechtold et al (2014)
            !xff_dicycle  = (AA1(i)-AA1_BL(i))/tau_ecmwf(i)

            !--- Bechtold et al (2014) + Becker et al (2021)
            xff_dicycle = (1.-alpha_adv(i))*aa1(i) + alpha_adv(i)*aa1_radpbl(i) &
                          + alpha_adv(i)*Q_adv(i) - aa1_bl(i)

            !xff_dicycle  = Q_adv(i)

            xff_dicycle = xff_dicycle/tau_ecmwf(i)

            if (xk(i) .lt. 0) xf_dicycle(i) = max(0., -xff_dicycle/xk(i))
            xf_dicycle(i) = xf_ens(i, 10) - xf_dicycle(i)
!----------
!            if(xk(i).lt.0) then
!                xf_dicycle(i)= max(0.,-xff_dicycle/xk(i))
!                xf_dicycle(i)= xf_ens(i,10)-xf_dicycle(i)
!            else
!                xf_dicycle(i)= 0.
!            endif
!            xf_dicycle(i)= xf_ens(i,10)-xf_dicycle(i)
!----------
         end do

      elseif (DICYCLE == 3) then
         do i = its, itf
            xf_dicycle(i) = 0.

            if (ierr(i) /= 0) cycle

            xff_dicycle = (1.-alpha_adv(i))*aa1(i) + alpha_adv(i)*(aa1_radpbl(i) + aa1_adv(i)) - aa1_bl(i)
                          !                        +      alpha_adv(i) *(AA1_RADPBL(i) + AA1_ADV(i) - AA0(i)) &
!--------tmp
!           xff_dicycle  =  AA1_ADV(i)
!--------tmp

            xff_dicycle = xff_dicycle/tau_ecmwf(i)

            if (xk(i) .lt. 0) xf_dicycle(i) = max(0., -xff_dicycle/xk(i))
            xf_dicycle(i) = xf_ens(i, 10) - xf_dicycle(i)
         end do

      elseif (DICYCLE == 4) then
         do i = its, itf
            xf_dicycle(i) = 0.
            if (ierr(i) /= 0) cycle
            !the signal "-" is to convert from Pa/s to kg/m2/s
            if (xk_x(i) > 0.) xf_dicycle(i) = max(0., -aa1_bl(I))/xk_x(i)

            xf_ens(i, :) = xf_dicycle(i)
            xf_dicycle(i) = 0.0
         end do

      else
         xf_dicycle(:) = 0.0

      end if
      !------------------------------------
      !-
      !- add the kinetic energy at the gust front at the
      !- mass flux closure
      !-
      if (ADD_COLDPOOL_CLOS == 4) then
         do i = its, itf
            if (ierr(i) /= 0 .or. xk(i) >= 0) cycle
            xf_coldpool(i) = -(0.5*wlpool(i)**2/tau_ecmwf(i))/xk(i)
         end do
      end if

   end subroutine cupForcingEns3d

   ! ------------------------------------------------------------------------------------
   subroutine getPartitionLiqIce(ierr, tn, z1, zo_cup, po_cup, p_liq_ice, melting_layer, itf, ktf, its, ite, kts, kte, cumulus)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'getPartitionLiqIce' ! subroutine name

      real, parameter ::  p_t1 = 276.16, p_z_meltlayer1 = 4000.
      real, parameter ::  p_z_meltlayer2 = 6000., p_delt = 3.   
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte

      integer, intent(in) :: ierr(:)

      real, intent(in) :: tn(:, :)
      real, intent(in) :: po_cup(:, :)
      real, intent(in) :: zo_cup(:, :)
      real, intent(in) :: z1(:)

      character(len=*), intent(in) :: cumulus

      real, intent(inout) :: p_liq_ice(:, :)
      real, intent(inout) :: melting_layer(:, :)

      !Local variables:
            integer :: i, k
      real :: dp, height
      real, dimension(its:ite) :: norm
   
      p_liq_ice(:, :) = 1.
      melting_layer(:, :) = 0.
      !-- get function of T for partition of total condensate into liq and ice phases.
      if (p_melt_glac .and. trim(cumulus) == 'deep') then
         do k = kts, ktf
            do i = its, itf
               if (ierr(i) /= 0) cycle
               p_liq_ice(i, k) = FractLiqF(tn(i, k))
            end do
         end do

         !-- define the melting layer (the layer will be between T_0+1 < TEMP < T_1
         !-- definition em terms of temperatura
         do k = kts, ktf
            do i = its, itf
               if (ierr(i) /= 0) cycle
               if (tn(i, k) <= c_t00 - p_delt) then
                  melting_layer(i, k) = 0.
               elseif (tn(i, k) < c_t00 + p_delt .and. tn(i, k) > c_t00 - p_delt) then
                  melting_layer(i, k) = ((tn(i, k) - (c_t00 - p_delt))/(2.*p_delt))**2
               else
                  melting_layer(i, k) = 1.
               end if
               melting_layer(i, k) = melting_layer(i, k)*(1.-melting_layer(i, k))
            end do
         end do
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
         norm(:) = 0.
         do k = kts, ktf - 1
            do i = its, itf
               if (ierr(i) /= 0) cycle
               dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))
               norm(i) = norm(i) + melting_layer(i, k)*dp/c_grav
            end do
         end do
         do i = its, itf
            if (ierr(i) /= 0) cycle
            melting_layer(i, :) = melting_layer(i, :)/(norm(i) + 1.e-6)*(100*(po_cup(i, kts) - po_cup(i, ktf))/c_grav)
            !print*,"i2=",i,maxval(melting_layer(i,:)),minval(melting_layer(i,:)),norm(i)
         end do
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
      end if
   end subroutine getPartitionLiqIce

   ! ------------------------------------------------------------------------------------
   subroutine getMeltingProfile(ierr, tn_cup, po_cup, p_liq_ice, melting_layer, qrco, pwo, edto, pwdo, melting, itf, ktf, its &
                              , ite, kts, kte, cumulus)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'getMeltingProfile' ! subroutine name
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte
      integer, intent(in) :: ierr(:)

      real, intent(in) :: edto(:)
      real, intent(in) :: tn_cup(:, :)
      real, intent(in) :: po_cup(:, :)
      real, intent(in) :: qrco(:, :)
      real, intent(in) :: pwo(:, :)
      real, intent(in) :: pwdo(:, :)
      real, intent(in) :: p_liq_ice(:, :)
      real, intent(in) :: melting_layer(:, :)

      character(len=*), intent(in) :: cumulus

      real, intent(inout) :: melting(:, :)

      !Local variables:
      integer :: i, k
      real :: dp
      real, dimension(its:ite)         :: norm, total_pwo_solid_phase
      real, dimension(its:ite, kts:kte) :: pwo_solid_phase, pwo_eff
      
      if (p_melt_glac .and. trim(cumulus) == 'deep') then

         norm = 0.0
         pwo_solid_phase = 0.0
         pwo_eff = 0.0
         melting = 0.0
         !-- set melting mixing ratio to zero for columns that do not have deep convection
         do i = its, itf
            if (ierr(i) > 0) melting(i, :) = 0.
         end do

         !-- now, get it for columns where deep convection is activated
         total_pwo_solid_phase(:) = 0.

         do k = kts, ktf - 1
            do i = its, itf
               if (ierr(i) /= 0) cycle
               dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))
               !-- effective precip (after evaporation by downdraft)
               !-- pwdo is not defined yet
               !pwo_eff(i,k) = 0.5*(pwo(i,k)+pwo(i,k+1) + edto(i)*(pwdo(i,k)+pwdo(i,k+1)))
               pwo_eff(i, k) = 0.5*(pwo(i, k) + pwo(i, k + 1))
               !-- precipitation at solid phase(ice/snow)
               pwo_solid_phase(i, k) = (1.-p_liq_ice(i, k))*pwo_eff(i, k)
               !-- integrated precip at solid phase(ice/snow)
               total_pwo_solid_phase(i) = total_pwo_solid_phase(i) + pwo_solid_phase(i, k)*dp/c_grav
            end do
         end do

         do k = kts, ktf
            do i = its, itf
               if (ierr(i) /= 0) cycle
               !-- melting profile (kg/kg)
               melting(i, k) = melting_layer(i, k)*(total_pwo_solid_phase(i)/(100*(po_cup(i, kts) - po_cup(i, ktf))/c_grav))
               !print*,"mel=",k,melting(i,k),pwo_solid_phase(i,k),po_cup(i,k)
            end do
         end do

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
         melting(:, :) = 0.
      end if
   end subroutine getMeltingProfile

   ! ------------------------------------------------------------------------------------
   subroutine keToHeating(itf, ktf, its, ite, kts, kte, ktop, ierr, po_cup, us, vs, dellu, dellv, dellat)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'keToHeating' ! subroutine name
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte

      integer, intent(in) :: ierr(:)
      integer, intent(in) :: ktop(:)

      real, intent(in) :: po_cup(:, :)
      real, intent(in) :: us(:, :)
      real, intent(in) :: vs(:, :)
      real, intent(in) :: dellu(:, :)
      real, intent(in) :: dellv(:, :)

      real, intent(inout) :: dellat(:, :)

      !Local variables:
      real :: dts, fp, dp, fpi
      integer ::i, k

      ! since kinetic energy is being dissipated, add heating accordingly (from ECMWF)
      do i = its, itf
         if (ierr(i) /= 0) cycle
         dts = 0.
         fpi = 0.
         do k = kts, ktop(i)
            dp = (po_cup(i, k) - po_cup(i, k + 1))*100.
            !total KE dissiptaion estimate
            dts = dts - (dellu(i, k)*us(i, k) + dellv(i, k)*vs(i, k))*dp/c_grav
            !
            ! fpi needed for calcualtion of conversion to pot. energyintegrated
            fpi = fpi + sqrt(dellu(i, k)*dellu(i, k) + dellv(i, k)*dellv(i, k))*dp
         end do
         if (fpi .gt. 0.) then
            do k = kts, ktop(i)
               fp = sqrt((dellu(i, k)*dellu(i, k) + dellv(i, k)*dellv(i, k)))/fpi
               dellat(i, k) = dellat(i, k) + fp*dts*c_grav/real(c_cp)
            end do
         end if
      end do

   end subroutine keToHeating

   ! ------------------------------------------------------------------------------------
   subroutine getInCloudScChemUp(cumulus, fscav, mtp, se, se_cup, sc_up, pw_up, tot_pw_up_chem &
                               , z_cup, rho, po, po_cup, qrco, tempco, pwo, zuo, up_massentro, up_massdetro &
                               , vvel2d, vvel1d, start_level, k22, kbcon, ktop, klcl, ierr, xland, itf, ktf &
                               , its, ite, kts, kte)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'getInCloudScChemUp' ! subroutine name
   
      real, parameter :: p_scav_eff = 0.6  
      !! for smoke : Chuang et al. (1992) J. Atmos. Sci.
      real, parameter :: p_cte_w_upd = 10. 
      !! m/s
      !    real, parameter :: kc = 5.e-3  
      !! s-1
      real, parameter :: p_kc = 2.e-3
      !! autoconversion parameter in GF is lower than what is used in GOCART s-1

      !Variables (input, output, inout)
       integer, intent(in)  :: itf, ktf, its, ite, kts, kte, mtp

      integer, intent(in) :: ierr(:)
      integer, intent(in) :: kbcon(:)
      integer, intent(in) :: ktop(:)
      integer, intent(in) :: k22(:)
      integer, intent(in) :: klcl(:)
      integer, intent(in) :: start_level(:)

      real, intent(in) :: fscav(:)
      real, intent(in) :: se(:, :, :)
      real, intent(in) :: se_cup(:, :, :)
      real, intent(in) :: z_cup(:, :)
      real, intent(in) :: rho(:, :)
      real, intent(in) :: po_cup(:, :)
      real, intent(in) :: qrco(:, :)
      real, intent(in) :: tempco(:, :)
      real, intent(in) :: pwo(:, :)
      real, intent(in) :: zuo(:, :)
      real, intent(in) :: up_massentro(:, :)
      real, intent(in) :: up_massdetro(:, :)
      real, intent(in) :: po(:, :)
      real, intent(in) :: vvel2d(:, :)
      real, intent(in) :: vvel1d(:)
      real, intent(in) :: xland(:)

      character(len=*), intent(in) :: cumulus

      real, intent(out) :: sc_up(:, :, :)
      real, intent(out) :: pw_up(:, :, :)
      real, intent(out) :: tot_pw_up_chem(:, :)

      !Local variables:
      real, dimension(mtp, its:ite) ::  sc_b
      real, dimension(mtp) :: conc_mxr
      real :: x_add, dz, xzz, xzd, xze, denom, henry_coef, w_upd, fliq, dp
      integer :: i, k, ispc
      real, dimension(mtp, its:ite, kts:kte) ::  factor_temp

      !--initialization
      sc_up = se_cup
      pw_up = 0.0
      tot_pw_up_chem = 0.0

      if (USE_TRACER_SCAVEN == 2 .and. cumulus /= 'shallow') then
         factor_temp = 1.
         do i = its, itf
            if (ierr(i) /= 0) cycle
            do ispc = 1, mtp
               ! - if tracer is type "carbon" then set coefficient to 0 for hydrophobic
               if (trim(chem_name(ispc) (1:len_trim('OCphobic'))) == 'OCphobic') factor_temp(ispc, :, :) = 0.0

               ! - suppress scavenging most aerosols at cold T except BCn1 (hydrophobic), dust, and HNO3
               if (trim(chem_name(ispc) (1:len_trim('BCphobic'))) == 'BCphobic') then
                  where (tempco < 258.) factor_temp(ispc, :, :) = 0.0
               end if

               if (trim(chem_name(ispc)) == 'sulfur' .or. &
                   trim(chem_name(ispc) (1:len_trim('ss'))) == 'ss' .or. & ! 'seasalt'
                   trim(chem_name(ispc)) == 'SO2' .or. &
                   trim(chem_name(ispc)) == 'SO4' .or. &
                   trim(chem_name(ispc)) == 'nitrate' .or. &
                   trim(chem_name(ispc)) == 'bromine' .or. &
                   trim(chem_name(ispc)) == 'NH3' .or. &
                   trim(chem_name(ispc)) == 'NH4a') then

                  where (tempco < 258.) factor_temp(ispc, :, :) = 0.0
               end if

            end do
         end do
      end if

      do i = its, itf
         if (ierr(i) /= 0) cycle
         !start_level(i) = klcl(i)
         !start_level(i) = k22(i)

         do ispc = 1, mtp
            call getCloudBc(cumulus, kts, kte, ktf, xland(i), po(i, kts:kte), se_cup(ispc, i, kts:kte), sc_b(ispc, i), k22(i))
         end do
         do k = kts, start_level(i)
            sc_up(:, i, k) = sc_b(:, i)
            !sc_up   (:,i,k) = se_cup(:,i,k)
         end do
      end do

      do i = its, itf
         if (ierr(i) /= 0) cycle
         loopk: do k = start_level(i) + 1, ktop(i) + 1

            !-- entr,detr, mass flux ...
            xzz = zuo(i, k - 1)
            xzd = 0.5*up_massdetro(i, k - 1)
            xze = up_massentro(i, k - 1)
            denom = (xzz - xzd + xze)

            !-- transport + mixing
            if (denom > 0.) then
               sc_up(:, i, k) = (sc_up(:, i, k - 1)*xzz - sc_up(:, i, k - 1)*xzd + se(:, i, k - 1)*xze)/denom
            else
               sc_up(:, i, k) = sc_up(:, i, k - 1)
            end if

            !-- scavenging section
            if (USE_TRACER_SCAVEN == 0 .or. cumulus == 'shallow') cycle loopk
            dz = z_cup(i, k) - z_cup(i, k - 1)

            !-- in-cloud vert velocity for scavenging formulation 2
            !           w_upd = cte_w_upd
            !           w_upd = vvel1d(i)
            w_upd = vvel2d(i, k)

            do ispc = 1, mtp
               if (fscav(ispc) > 1.e-6) then ! aerosol scavenging

                  !--formulation 1 as in GOCART with RAS conv_par
                  if (USE_TRACER_SCAVEN == 1) pw_up(ispc, i, k) = max(0., sc_up(ispc, i, k)*(1.-exp(-fscav(ispc)*(dz/1000.))))

                  !--formulation 2 as in GOCART
                  if (USE_TRACER_SCAVEN == 2) pw_up(ispc, i, k) = max(0., sc_up(ispc, i, k) * (1.-exp(-chem_adj_autoc(ispc)*p_kc &
                                            * (dz/w_upd)))*factor_temp(ispc, i, k))

                  !--formulation 3 - orignal GF conv_par
                  if (USE_TRACER_SCAVEN == 3) then
                     !--- cloud liquid water tracer concentration
                     conc_mxr(ispc) = p_scav_eff*sc_up(ispc, i, k) !unit [kg(aq)/kg(air)]  for aerosol/smoke
                     !---   aqueous-phase concentration in rain water
                     pw_up(ispc, i, k) = conc_mxr(ispc)*pwo(i, k)/(1.e-8 + qrco(i, k))
                  end if

                  !---(in cloud) total mixing ratio in gas and aqueous phases
                  sc_up(ispc, i, k) = sc_up(ispc, i, k) - pw_up(ispc, i, k)

                  !
               elseif (hcts(ispc)%hstar > 1.e-6) then ! tracer gas phase scavenging

                  !--- equilibrium tracer concentration - Henry's law
                  henry_coef = Henry(ispc, tempco(i, k), rho(i, k))

                  if (USE_TRACER_SCAVEN == 3) then
                     !--- cloud liquid water tracer concentration
                     conc_mxr(ispc) = (henry_coef*qrco(i, k)/(1.+henry_coef*qrco(i, k)))*sc_up(ispc, i, k)
                     !
                     !---   aqueous-phase concentration in rain water
                     pw_up(ispc, i, k) = conc_mxr(ispc)*pwo(i, k)/(1.e-8 + qrco(i, k))

                  else

                     !-- this the 'alpha' parameter in Eq 8 of Mari et al (2000 JGR) = X_aq/X_total
                     fliq = henry_coef*qrco(i, k)/(1.+henry_coef*qrco(i, k))

                     !---   aqueous-phase concentration in rain water
                     pw_up(ispc, i, k) = max(0., sc_up(ispc, i, k) &
                                       *(1.-exp(-fliq*chem_adj_autoc(ispc)*p_kc*dz/w_upd)))!*factor_temp(ispc,i,k))

                  end if

                  !---(in cloud) total mixing ratio in gas and aqueous phases
                  sc_up(ispc, i, k) = sc_up(ispc, i, k) - pw_up(ispc, i, k)

                  !
                  !---(in cloud)  mixing ratio in aqueous phase
                  !sc_up_aq(ispc,i,k) = conc_mxr(ispc) !if using set to zero at the begin.
               end if
            end do
            !
            !-- total aerosol/gas in the rain water
            dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))

            tot_pw_up_chem(:, i) = tot_pw_up_chem(:, i) + pw_up(:, i, k)*dp/c_grav
         end do loopk
         !
         !----- get back the in-cloud updraft gas-phase mixing ratio : sc_up(ispc,k)
         !          do k=start_level(i)+1,ktop(i)+1
         !            do ispc = 1,mtp
         !             sc_up(ispc,i,k) = sc_up(ispc,i,k) - sc_up_aq(ispc,i,k)
         !            enddo
         !          enddo
      end do
   end subroutine getInCloudScChemUp

   !---------------------------------------------------------------------------------------------------
   function Henry(ispc, temp, rhoair) result(henry_coef)
      !! ## compute Henry's constant
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! —
      !! **Full description**:
      !!
      !! calculate Henry's constant for solubility of gases into cloud water
      !! inputs : ak0(ispc), dak(ispc),  hstar(ispc), dhr(ispc)
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !!  
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'henry' ! function name
   
      !Variables (input):
      integer, intent(in) :: ispc
      real, intent(in) :: temp
      real, intent(in) :: rhoair
   
      !Local variables:
      real :: henry_coef ! output

      real :: fct, tcorr, corrh

      ! aqueous-phase concentrations XXXa [mol/m3(air)]!
      ! gas-phase concentrations XXXg [mol/m3(air)]!
      ! Henry constants XXXh for scavenging [mol/(l*atm)]!
      ! converted to [(mol(aq)/m3(aq))/(mol(g)/m3(air))], i.e. dimensionless!
      ! in equilibrium XXXa = XXXh * LWC * XXXg!
      tcorr = 1./temp - c_temp0i

      !-P. Colarco corrected the expression below
      !fct   = conv7 * rgas_ * temp ! - for henry_coef in units 1/m3
      fct = c_rgas_atm*temp ! - for henry_coef dimensioless

      !-taking into account the acid dissociation constant
      ! ak=ak0*exp(dak*(1/t-1/298))
      corrh = 1.+hcts(ispc)%ak0*exp(hcts(ispc)%dak*tcorr)/c_hplus

      !-- for concentration in mol[specie]/mol[air] - Eq 5 in 'Compilation of Henry's law constants (version 4.0) for
      !-- water as solvent, R. Sander, ACP 2015'.
      henry_coef = hcts(ispc)%hstar*exp(hcts(ispc)%dhr*tcorr)*fct*corrh

   end function Henry

   ! ---------------------------------------------------------------------------------------------------
   subroutine getInCloudScChemDd(cumulus, fscav, mtp, se, se_cup, sc_dn, pw_dn, pw_up, sc_up, tot_pw_up_chem &
                               , tot_pw_dn_chem, z_cup, rho, po_cup, qrcdo, pwdo, pwevo, edto, zdo, dd_massentro &
                               , dd_massdetro, pwavo, pwo , jmin, ierr, itf, ktf, its, ite, kts, kte)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'getInCloudScChemDd' ! subroutine name
   
      !Variables (input, output, inout)
      integer, intent(in)  :: itf, ktf, its, ite, kts, kte, mtp

      integer, intent(in) :: ierr(:)
      integer, intent(in) :: jmin(:)

      real, intent(in) :: se(:, :, :)
      real, intent(in) :: se_cup(:, :, :)
      real, intent(in) :: pw_up(:, :, :)
      real, intent(in) :: sc_up(:, :, :)
      real, intent(in) :: fscav(:)
      real, intent(in) :: edto(:)
      real, intent(in) :: pwavo(:)
      real, intent(in) :: pwevo(:)
      real, intent(in) :: z_cup(:, :)
      real, intent(in) :: rho(:, :)
      real, intent(in) :: po_cup(:, :)
      real, intent(in) :: qrcdo(:, :)
      real, intent(in) :: pwdo(:, :)
      real, intent(in) :: zdo(:, :)
      real, intent(in) :: dd_massentro(:, :)
      real, intent(in) :: dd_massdetro(:, :)
      real, intent(in) :: pwo(:, :)
      real, intent(in) :: tot_pw_up_chem(:, :)

      character(len=*), intent(in)  :: cumulus

      real, intent(out) :: sc_dn(:, :, :)
      real, intent(out) :: pw_dn(:, :, :)
      real, intent(out) :: tot_pw_dn_chem(:, :)

      !Local variables:
      real, dimension(mtp) :: conc_mxr
      real :: x_add, dz, xzz, xzd, xze, denom, evaporate, pwdper, x1, frac_evap, dp, xkk
      integer :: i, k, ispc

      sc_dn = 0.0
      pw_dn = 0.0
      tot_pw_dn_chem = 0.0
      if (cumulus == 'shallow') return

      do i = its, itf
         if (ierr(i) /= 0) cycle

         !--- fration of the total rain that was evaporated
         frac_evap = -pwevo(i)/(1.e-16 + pwavo(i))

         !--- scalar concentration in-cloud - downdraft

         !--- at k=jmim
         k = jmin(i)
         pwdper = pwdo(i, k)/(1.e-16 + pwevo(i))*frac_evap  ! > 0
         if (USE_TRACER_EVAP == 0) pwdper = 0.0

         dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))

         do ispc = 1, mtp
            !--downdrafts will be initiate with a mixture of 50% environmental and in-cloud concentrations
            sc_dn(ispc, i, k) = se_cup(ispc, i, k)
            !sc_dn(ispc,i,k) = 0.9*se_cup(ispc,i,k)+0.1*sc_up(ispc,i,k)

            pw_dn(ispc, i, k) = -pwdper*tot_pw_up_chem(ispc, i)*c_grav/dp
            sc_dn(ispc, i, k) = sc_dn(ispc, i, k) - pw_dn(ispc, i, k)
            tot_pw_dn_chem(ispc, i) = tot_pw_dn_chem(ispc, i) + pw_dn(ispc, i, k)*dp/c_grav
         end do
         !
         !--- calculate downdraft mass terms
         do k = jmin(i) - 1, kts, -1
            xzz = zdo(i, k + 1)
            xzd = 0.5*dd_massdetro(i, k)
            xze = dd_massentro(i, k)

            denom = (xzz - xzd + xze)
            !-- transport + mixing
            if (denom > 0.) then
               sc_dn(:, i, k) = (sc_dn(:, i, k + 1)*xzz - sc_dn(:, i, k + 1)*xzd + se(:, i, k)*xze)/denom
            else
               sc_dn(:, i, k) = sc_dn(:, i, k + 1)
            end if
            !
            !-- evaporation term
            if (USE_TRACER_EVAP == 0) cycle

            dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))

            !-- fraction of evaporated precip per layer
            pwdper = pwdo(i, k)/(1.e-16 + pwevo(i))! > 0

            !-- fraction of the total precip that was actually evaporated at layer k
            pwdper = pwdper*frac_evap

            !-- sanity check
            pwdper = min(1., max(pwdper, 0.))

            do ispc = 1, mtp
               !-- amount evaporated by the downdraft from the precipitation
               pw_dn(ispc, i, k) = -pwdper*tot_pw_up_chem(ispc, i)*c_grav/dp ! < 0. => source term for the downdraft tracer concentration

               !if(ispc==1) print*,"pw=",pwdper,tot_pw_up_chem (ispc,i),pwevo(i)/pwavo(i),pwdo(i,k)/(1.e-16+pwo(i,k))

               !-- final tracer in the downdraft
               sc_dn(ispc, i, k) = sc_dn(ispc, i, k) - pw_dn(ispc, i, k) ! observe that -pw_dn is > 0.

               !-- total evaporated tracer
               tot_pw_dn_chem(ispc, i) = tot_pw_dn_chem(ispc, i) + pw_dn(ispc, i, k)*dp/c_grav

               !print*,"to=",k,tot_pw_dn_chem(ispc,i),pwdo(i,k)/(1.e-16+pwevo(i)),frac_evap,tot_pw_dn_chem(ispc,i)/tot_pw_up_chem (ispc,i)

            end do
         end do
         !
      end do
   end subroutine getInCloudScChemDd

   ! -----------------------------------------------------------------------------------------
   subroutine fct1d3(ktop, n, dt, z, tracr, massflx, trflx_in, del_out)
      !! ## modify a 1-D array of tracer fluxes for the purpose of maintaining monotonicity
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! modify a 1-D array of tracer fluxes for the purpose of maintaining
      !! monotonicity (including positive-definiteness) in the tracer field
      !! during tracer transport.
      !! the underlying transport equation is   (d tracr/dt) = - (d trflx/dz)
      !! where  dz = |z(k+1)-z(k)| (k=1,...,n) and  trflx = massflx * tracr
      !! note: tracr is carried in grid cells while z and fluxes are carried on
      !! interfaces. interface variables at index k are at grid location k-1/2.
      !! sign convention: mass fluxes are considered positive in +k direction.
      !! massflx and trflx_in  must be provided independently to allow the
      !! algorithm to generate an auxiliary low-order (diffusive) tracer flux
      !! as a stepping stone toward the final product trflx_out.
      !!
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'fct1d3' ! subroutine name
      
      logical, parameter :: p_hi_order = .false.
      real, parameter :: p_epsil = 1.e-22
      !! prevent division by zero
      real, parameter :: p_damp = 1.
      !! damper of antidff flux (1=no damping)

      !Variables (input, output, inout)
      integer, intent(in) :: n
      !! number of grid cells
      integer, intent(in) :: ktop

      real, intent(in) :: dt
      !! transport time step
      real, intent(in) :: z(n + 0)
      !! location of cell interfaces
      real, intent(in) :: tracr(n)
      !! the transported variable
      real, intent(in) :: massflx(n + 0)
      !! mass flux across interfaces
      real, intent(in) :: trflx_in(n + 0)
      !! original tracer flux
      
      real, intent(out):: del_out(n + 0)                
      !! modified tracr flux   
      
      !Local variables:
      real :: trflx_out(n + 0)                
      !! modified tracr flux
      integer :: k, km1, kp1
      logical :: nan, error = .false., vrbos = .false.
      real :: dtovdz(n), trmax(n), trmin(n), flx_lo(n + 0), antifx(n + 0), clipped(n + 0), soln_hi(n), totlin(n), totlout(n)
      real ::  soln_lo(n), clipin(n), clipout(n), arg

      nan(arg) = .not. (arg .ge. 0. .or. arg .le. 0.) ! NaN detector
      soln_lo(:) = 0.
      antifx(:) = 0.
      clipout(:) = 0.
      flx_lo(:) = 0.

      do k = 1, ktop
         dtovdz(k) = .01*dt/abs(z(k + 1) - z(k))                ! time step / grid spacing
         !     if (z(k).eq.z(k+1)) error=.true.
      end do
      if (vrbos .or. error) print '(a/(8es10.3))', '(fct1d) dtovdz =', dtovdz(1:ktop)

      do k = 2, ktop
         if (massflx(k) > 0.) then
            flx_lo(k) = massflx(k)*tracr(k - 1)              ! low-order flux, upstream
         else
            flx_lo(k) = massflx(k)*tracr(k)                ! low-order flux, upstream
         end if
         antifx(k) = trflx_in(k) - flx_lo(k)                ! antidiffusive flux
      end do
      flx_lo(1) = trflx_in(1)
      flx_lo(ktop + 1) = trflx_in(ktop + 1)
      antifx(1) = 0.
      antifx(ktop + 1) = 0.
      ! --- clip low-ord fluxes to make sure they don't violate positive-definiteness
      do k = 1, ktop
         totlout(k) = max(0., flx_lo(k + 1)) - min(0., flx_lo(k))         ! total flux out
         clipout(k) = min(1., tracr(k)/max(p_epsil, totlout(k))/(1.0001*dtovdz(k)))
      end do

      do k = 2, ktop
         if (massflx(k) .ge. 0.) then
            flx_lo(k) = flx_lo(k)*clipout(k - 1)
         else
            flx_lo(k) = flx_lo(k)*clipout(k)
         end if
      end do
      if (massflx(1) .lt. 0.) flx_lo(1) = flx_lo(1)*clipout(1)
      if (massflx(ktop + 1) .gt. 0.) flx_lo(ktop + 1) = flx_lo(ktop + 1)*clipout(ktop)

      ! --- a positive-definite low-order (diffusive) solution can now be  constructed
      do k = 1, ktop
         soln_lo(k) = tracr(k) - (flx_lo(k + 1) - flx_lo(k))*dtovdz(k)        ! low-ord solutn
         del_out(k) = -c_grav*(flx_lo(k + 1) - flx_lo(k))*dtovdz(k)/dt
      end do

      if (.not. p_hi_order) return

      soln_hi(:) = 0.
      clipin(:) = 0.
      trmin(:) = 0.
      trmax(:) = 0.
      clipped(:) = 0.
      trflx_out(:) = 0.

      do k = 1, ktop
         km1 = max(1, k - 1)
         kp1 = min(n, k + 1)
         trmax(k) = max(soln_lo(km1), soln_lo(k), soln_lo(kp1), tracr(km1), tracr(k), tracr(kp1)) ! upper bound
         trmin(k) = max(0., min(soln_lo(km1), soln_lo(k), soln_lo(kp1), tracr(km1), tracr(k), tracr(kp1)))  ! lower bound
      end do

      do k = 1, ktop
         totlin(k) = max(0., antifx(k)) - min(0., antifx(k + 1))                ! total flux in
         totlout(k) = max(0., antifx(k + 1)) - min(0., antifx(k))                ! total flux out

         clipin(k) = min(p_damp, (trmax(k) - soln_lo(k))/max(p_epsil, totlin(k))/(1.0001*dtovdz(k)))
         clipout(k) = min(p_damp, (soln_lo(k) - trmin(k))/max(p_epsil, totlout(k))/(1.0001*dtovdz(k)))

         if (nan(clipin(k))) print *, '(fct1d) error: clipin is NaN,  k=', k
         if (nan(clipout(k))) print *, '(fct1d) error: clipout is NaN,  k=', k

         if (clipin(k) .lt. 0.) then
            print 100, '(fct1d) error: clipin < 0 at k =', k, &
               'clipin', clipin(k), 'trmax', trmax(k), 'soln_lo', soln_lo(k), &
               'totlin', totlin(k), 'dt/dz', dtovdz(k)
            error = .true.
         end if
         if (clipout(k) .lt. 0.) then
            print 100, '(fct1d) error: clipout < 0 at k =', k, &
               'clipout', clipout(k), 'trmin', trmin(k), 'soln_lo', soln_lo(k), &
               'totlout', totlout(k), 'dt/dz', dtovdz(k)
            error = .true.
         end if
100      format(a, i3/(4(a10, "=", es9.2)))
      end do

      do k = 2, ktop
         if (antifx(k) .gt. 0.) then
            clipped(k) = antifx(k)*min(clipout(k - 1), clipin(k))
         else
            clipped(k) = antifx(k)*min(clipout(k), clipin(k - 1))
         end if
         trflx_out(k) = flx_lo(k) + clipped(k)
         if (nan(trflx_out(k))) then
            print *, '(fct1d) error: trflx_out is NaN,  k=', k
            error = .true.
         end if
      end do

      trflx_out(1) = trflx_in(1)
      trflx_out(ktop + 1) = trflx_in(ktop + 1)
      do k = 1, ktop
         soln_hi(k) = tracr(k) - (trflx_out(k + 1) - trflx_out(k))*dtovdz(k)
         del_out(k) = -c_grav*(trflx_out(k + 1) - trflx_out(k))*dtovdz(k)/dt
         !write(32,*)'3',k,soln_lo(k),soln_hi(k)
      end do

      if (vrbos .or. error) then
         do k = 2, ktop
            write (32, 99) k, &
               'tracr(k)', tracr(k), &
               'flx_in(k)', trflx_in(k), &
               'flx_in(k+1)', trflx_in(k + 1), &
               'flx_lo(k)', flx_lo(k), &
               'flx_lo(k+1)', flx_lo(k + 1), &
               'soln_lo(k)', soln_lo(k), &
               'trmin(k)', trmin(k), &
               'trmax(k)', trmax(k), &
               'totlin(k)', totlin(k), &
               'totlout(k)', totlout(k), &
               'clipin(k-1)', clipin(k - 1), &
               'clipin(k)', clipin(k), &
               'clipout(k-1)', clipout(k - 1), &
               'clipout(k)', clipout(k), &
               'antifx(k)', antifx(k), &
               'antifx(k+1)', antifx(k + 1), &
               'clipped(k)', clipped(k), &
               'clipped(k+1)', clipped(k + 1), &
               'flx_out(k)', trflx_out(k), &
               'flx_out(k+1)', trflx_out(k + 1), &
               'dt/dz(k)', dtovdz(k), &
               'final', tracr(k) - (trflx_out(k + 1) - trflx_out(k))*dtovdz(k)
99          format('(trc1d)   k =', i4/(3(a13, '=', es13.6)))
         end do
         if (error) stop '(fct1d error)'
      end if

   end subroutine fct1d3

   ! ---------------------------------------------------------------------------------------------------
   subroutine tridiag(m_size, aa, bb, cc, ff)
      !! ## solves the problem: aa*ff(k-1,t+1) + bb*ff(k,t+1) + cc*ff(k+1,t+1) = dd
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! This routine solves the problem: aa*ff(k-1,t+1) + bb*ff(k,t+1) + cc*ff(k+1,t+1) = dd
      !! an updated "ff" at time t+1 is the output
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'tridiag' ! subroutine name
   
      !Variables (input, output, inout)
      integer, intent(in) :: m_size

      real, intent(in) :: aa(:)
      real, intent(in) :: bb(:)

      real, intent(inout) :: cc(:)
      real, intent(inout) :: ff(:)

      !Local variables:
      real, dimension(m_size) :: qq
      integer :: k
      real :: pp

      cc(m_size) = 0.
      qq(1) = -cc(1)/bb(1)
      ff(1) = ff(1)/bb(1)
      do k = 2, m_size
         pp = 1./(bb(k) + aa(k)*qq(k - 1))
         qq(k) = -cc(k)*pp
         ff(k) = pp*(ff(k) - aa(k)*ff(k - 1))
      end do
      do k = m_size - 1, 1, -1
         ff(k) = ff(k) + qq(k)*ff(k + 1)
      end do

   end subroutine tridiag

   ! ---------------------------------------------------------------------------------------------------
   subroutine cupEnvClevChem(mtp, se_chem, se_cup_chem, ierr, itf, ktf, its, ite, kts, kte)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'cupEnvClevChem' ! subroutine name

      integer, parameter ::  p_clev_option = 2 
      !! use option 2
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte, mtp
      
      integer, intent(in) :: ierr(:)

      real, intent(in) :: se_chem(: , :, :)

      real, intent(out) :: se_cup_chem(: , :, :)

      !Local variables:
      integer :: i, k
      
      if (p_clev_option == 1) then
         !-- original version
         do i = its, itf
            if (ierr(i) /= 0) cycle
            do k = kts + 1, ktf
               se_cup_chem(1:mtp, i, k) = 0.5*(se_chem(1:mtp, i, k - 1) + se_chem(1:mtp, i, k))
            end do
            se_cup_chem(1:mtp, i, kts) = se_chem(1:mtp, i, kts)
            se_cup_chem(1:mtp, i, kte) = se_chem(1:mtp, i, ktf)
         end do
      else
         !-- version 2: se_cup (k+1/2) = se(k) => smoother profiles
         do i = its, itf
            if (ierr(i) /= 0) cycle
            do k = kts, ktf
               se_cup_chem(1:mtp, i, k) = se_chem(1:mtp, i, k)
            end do
         end do
      end if

   end subroutine cupEnvClevChem

   ! ---------------------------------------------------------------------------------------------------
   subroutine rainEvapBelowCloudBase(cumulus, itf, ktf, its, ite, kts, kte, ierr, kbcon, ktop, xmb, psur, xland &
                                 ,   qo_cup, t_cup, po_cup, qes_cup, pwavo, edto, pwevo, pwo, pwdo &
                                 ,   pre, prec_flx, evap_flx, outt, outq, outbuoy, evap_bcb)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'rainEvapBelowCloudBase' ! subroutine name

      real, parameter :: p_alpha1 = 5.44e-4 
      !! 1/sec
      real, parameter :: p_alpha2 = 5.09e-3 
      !! unitless
      real, parameter :: p_alpha3 = 0.5777 
      !! unitless
      real, parameter :: p_c_conv = 0.05
      !!conv fraction area, unitless
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte

      integer, intent(in) :: ierr(:)
      integer, intent(in) :: kbcon(:)
      integer, intent(in) :: ktop(:)

      real, intent(in) :: psur(:)
      real, intent(in) :: xland(:)
      real, intent(in) :: pwavo(:)
      real, intent(in) :: edto(:)
      real, intent(in) :: pwevo(:)
      real, intent(in) :: xmb(:)
      real, intent(in) :: po_cup(:, :)
      real, intent(in) :: qo_cup(:, :)
      real, intent(in) :: qes_cup(:, :)
      real, intent(in) :: pwo(:, :)
      real, intent(in) :: pwdo(:, :)
      real, intent(in) :: t_cup(:, :)

      character(len=*), intent(in) :: cumulus

      real, intent(inout) :: pre(:)
      real, intent(inout) :: outt(:, :)
      real, intent(inout) :: outq(:, :)
      real, intent(inout) :: outbuoy(:, :)
      real, intent(inout) :: prec_flx(:, :)
      real, intent(inout) :: evap_flx(:, :)

      real, intent(out)   :: evap_bcb(:, :)

      !Local variables:
      integer :: i, k
      real :: rh_cr, del_t, del_q, dp, q_deficit, pqsat, temp_pre
      real :: rh_cr_ocean, rh_cr_land
      real, dimension(its:ite) :: tot_evap_bcb, eff_c_conv

      if (trim(cumulus) == 'shallow') then
         rh_cr_ocean = 1.
         rh_cr_land = 1.
         eff_c_conv(:) = min(0.2, max(xmb(:), p_c_conv))
      else
         rh_cr_ocean = 0.95 !test 0.90
         rh_cr_land = 0.90
         eff_c_conv(:) = p_c_conv
      end if

      prec_flx = 0.0
      evap_flx = 0.0
      tot_evap_bcb = 0.0
      if (c0 < 1.e-6) return

      do i = its, itf

         if (ierr(i) /= 0) cycle

         !-- critical rel humidity  - check this, if the value is too small, not evapo will take place.
         rh_cr = rh_cr_ocean*xland(i) + rh_cr_land*(1.0 - xland(i))

         !if(xland(i)  < 0.90 ) then !- over land
         !  RH_cr = RH_cr_LAND
         !else
         !  RH_cr = RH_cr_OCEAN
         !endif

         do k = ktop(i), kts, -1

            dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))

            !p_liq_ice(i,k) = FractLiqF(tempco(i,k))

            !---rainfall evaporation below cloud base
            if (k <= kbcon(i)) then
               q_deficit = max(0., (rh_cr*qes_cup(i, k) - qo_cup(i, k)))
               !pqsat=SaturSpecHum(t_cup(i,k),po_cup(i,k))

               !--units here: kg[water]/kg[air}/sec
               evap_bcb(i, k) = eff_c_conv(i)*p_alpha1*q_deficit * (sqrt(po_cup(i, k)/psur(i))/p_alpha2*prec_flx(i, k + 1) &
                              / eff_c_conv(i))**p_alpha3

               !--units here: kg[water]/kg[air}/sec * kg[air]/m3 * m = kg[water]/m2/sec
               evap_bcb(i, k) = evap_bcb(i, k)*dp/c_grav

            else

               evap_bcb(i, k) = 0.0

            end if

            !-- before anything check if the evaporation already consumed all precipitation
            temp_pre = pre(i) - evap_bcb(i, k)
            if (temp_pre < 0.) evap_bcb(i, k) = pre(i)

            !-- get the net precitation flux after the local evaporation and downdraft
            prec_flx(i, k) = prec_flx(i, k + 1) - evap_bcb(i, k) + xmb(i)*(pwo(i, k) + edto(i)*pwdo(i, k))
            prec_flx(i, k) = max(0., prec_flx(i, k))

            evap_flx(i, k) = evap_flx(i, k + 1) + evap_bcb(i, k) - xmb(i)*edto(i)*pwdo(i, k)
            evap_flx(i, k) = max(0., evap_flx(i, k))

            tot_evap_bcb(i) = tot_evap_bcb(i) + evap_bcb(i, k)

            !-- feedback
            del_q = evap_bcb(i, k)*c_grav/dp          ! > 0., units: kg[water]/kg[air}/sec
            del_t = -evap_bcb(i, k)*c_grav/dp*(real(c_alvl)/real(c_cp)) ! < 0., units: K/sec

            outq(i, k) = outq(i, k) + del_q
            outt(i, k) = outt(i, k) + del_t
            !--- comment out 17nov
            !outbuoy(i,k) = outbuoy(i,k) + cp*del_t+xlv*del_q

            pre(i) = pre(i) - evap_bcb(i, k)

            !--for future use (rain and snow precipitation fluxes)
            !prec_flx_rain(k) = prec_flx(i,k)*(1.-p_liq_ice(k))
            !prec_flx_snow(k) = prec_flx(i,k)*    p_liq_ice(k)

         end do

         if (pre(i) < 0.) then
            print *, "prec evap neg for cumulus=", pre(i), trim(cumulus)
            call flush (6)
            !stop '@subroutine rain_evap_below_cloudbase'
         end if

      end do

   end subroutine rainEvapBelowCloudBase

   ! ---------------------------------------------------------------------------------------------------
   subroutine getPrecipFluxes(cumulus, klcl, kbcon, ktop, k22, ierr, xland, pre, xmb, pwo, pwavo, edto, pwevo, pwdo &
                            , t_cup, tempco, prec_flx, evap_flx, itf, ktf, its, ite, kts, kte)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'getPrecipFluxes' ! subroutine name
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte

      integer, intent(in) :: kbcon(:)
      integer, intent(in) :: ktop(:)
      integer, intent(in) :: k22(:)
      integer, intent(in) :: klcl(:)
      integer, intent(in) :: ierr(:)

      real, intent(in) :: xland(:)
      real, intent(in) :: pwavo(:)
      real, intent(in) :: pwevo(:)
      real, intent(in) :: edto(:)
      real, intent(in) :: pre(:)
      real, intent(in) :: xmb(:)
      real, intent(in) :: pwo(:, :)
      real, intent(in) :: pwdo(:, :)
      real, intent(in) :: t_cup(:, :)
      real, intent(in) :: tempco(:, :)

      character(len=*), intent(in) :: cumulus

      real, intent(out) :: prec_flx(:, :)
      !! units kg[water]/m2/s
      real, intent(out) :: evap_flx (:, :)
      !! units kg[water]/m2/s

      !Local variables:
      integer :: i, k

      prec_flx = 0.0
      evap_flx = 0.0
      if (c0 < 1.e-6) return

      do i = its, itf
         if (ierr(i) /= 0) cycle

         do k = ktop(i), kts, -1

            !--- precipitation flux (at 'cup' levels), units: kg[water]/m2/s
            prec_flx(i, k) = prec_flx(i, k + 1) + xmb(i)*(pwo(i, k) + edto(i)*pwdo(i, k))
            prec_flx(i, k) = max(0., prec_flx(i, k))

            !--- evaporation flux (at 'cup' levels), units: kg[water]/m2/s
            evap_flx(i, k) = evap_flx(i, k + 1) - xmb(i)*edto(i)*pwdo(i, k)
            evap_flx(i, k) = max(0., evap_flx(i, k))

            !
            !--for future use (rain and snow precipitation fluxes)
            !p_liq_ice(i,k) = FractLiqF(tempco(i,k))
            !prec_flx_rain(k) = prec_flx(i,k)*(1.-p_liq_ice(k))
            !prec_flx_snow(k) = prec_flx(i,k)*    p_liq_ice(k)
         end do

         !if(prec_flx   (i,kts) .ne. pre(i)) then
         !print*,"error=",100.*(prec_flx   (i,kts) - pre(i))/(1.e-16+pre(i)),pre(i),prec_flx   (i,kts)
         !STOP 'problem with water balance'
         !endif
      end do

   end subroutine getPrecipFluxes

   !------------------------------------------------------------------------------------
   function SaturSpecHum(pt, press) result(pqsat)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! —
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !!
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'SaturSpecHum' ! function name
   
      !Variables (input):
      real, intent(in) :: pt
      !! kelvin
      real, intent(in) :: press 
      !! hPa
      
      !Local variables:
      real :: pqsat !output

      real :: zew, zqs, zcor, foealfcu, foeewmcu
      
      foealfcu = min(1.0, ((max(c_rticecu, min(c_t00, pt)) - c_rticecu)*c_rtwat_rticecu_r)**2)
      foeewmcu = c_r2es*(foealfcu*exp(c_r3les*(pt - c_t00)/(pt - c_r4les)) &
               + (1.0 - foealfcu)*exp(c_r3ies*(pt - c_t00)/(pt - c_r4ies)))

      zew = foeewmcu
      zqs = zew/(100.*press)
      if (1.0 - c_retv*zqs > 0.) then
         zcor = 1.0/(1.0 - c_retv*zqs)  ! divide by zero
         pqsat = zqs*zcor
      else
         pqsat = c_max_qsat
      end if

   end function SaturSpecHum

   ! ---------------------------------------------------------------------------------------------------
   subroutine getJmin(cumulus, itf, ktf, its, ite, kts, kte, ierr, kdet, ktop, kbcon, jmin, ierrc, beta, depth_min, heso_cup &
                    , zo_cup, melting_layer)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'getJmin' ! subroutine name
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte

      integer, intent(in) :: ktop(:)
      integer, intent(in) :: kbcon(:)

      real, intent(in) :: heso_cup(:, :)
      real, intent(in) :: zo_cup(:, :)
      real, intent(in) :: melting_layer(:, :)
      real, intent(in) :: depth_min
      
      character(len=*), intent(in)    :: cumulus

      integer, intent(inout) :: ierr(:)
      integer, intent(inout) :: jmin(:)
      integer, intent(inout) :: kdet(:)
      
      real, intent(out) :: beta
      character(len=128), intent(out) :: ierrc(:)

      !Local variables:
      integer :: i, k, jmini, ki
      real :: dh, dz
      real, dimension(its:ite, kts:kte)  ::  hcdo
      logical :: keep_going

      if (trim(cumulus) == 'deep') beta = 0.05
      if (trim(cumulus) == 'mid') beta = 0.02

      if (trim(cumulus) == 'shallow') then
         beta = 0.02
         jmin(:) = 0
         return
      end if

      do i = its, itf
         if (ierr(i) /= 0) cycle

         if (trim(cumulus) == 'deep' .and. p_melt_glac) jmin(i) = max(jmin(i), maxloc(melting_layer(i, :), 1))
         !--- check whether it would have buoyancy, if there where
         !--- no entrainment/detrainment
         jmini = jmin(i)
         keep_going = .true.
         do while (keep_going)
            keep_going = .false.
            if (jmini - 1 .lt. kdet(i)) kdet(i) = jmini - 1
            if (jmini .ge. ktop(i) - 1) jmini = ktop(i) - 2
            ki = jmini
            hcdo(i, ki) = heso_cup(i, ki)
            dz = zo_cup(i, ki + 1) - zo_cup(i, ki)
            dh = 0.
            do k = ki - 1, 1, -1
               hcdo(i, k) = heso_cup(i, jmini)
               dz = zo_cup(i, k + 1) - zo_cup(i, k)
               dh = dh + dz*(hcdo(i, k) - heso_cup(i, k))
               if (dh .gt. 0.) then
                  jmini = jmini - 1
                  if (jmini .gt. 5) then
                     keep_going = .true.
                  else
                     ierr(i) = 9
                     ierrc(i) = "could not find jmini9"
                     exit
                  end if
               end if
            end do
         end do
         jmin(i) = jmini
         if (jmini .le. 5) then
            ierr(i) = 4
            ierrc(i) = "could not find jmini4"
         end if
      end do

      ! - must have at least depth_min m between cloud convective base and cloud top.
      do i = its, itf
         if (ierr(i) /= 0) cycle
         if (jmin(i) - 1 .lt. kdet(i)) kdet(i) = jmin(i) - 1
         if (-zo_cup(i, kbcon(i)) + zo_cup(i, ktop(i)) .lt. depth_min) then
            ierr(i) = 6
            ierrc(i) = "cloud depth very shallow"
         end if
      end do

   end subroutine getJmin

   ! ------------------------------------------------------------------------------------
   subroutine precipCwvFactor(itf, ktf, its, ite, kts, kte, ierr, t, po, qo, po_cup, cumulus, p_cwv_ave)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'precipCwvFactor' ! subroutine name

      real, parameter :: p_fpkup = 0.8  
      !! 90% of precip occurs above 80% of critical w
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte

      integer, intent(in) :: ierr(:)

      real, intent(in) :: t(:, :)
      real, intent(in) :: po(:, :)
      real, intent(in) :: qo(:, :)
      real, intent(in) :: po_cup(:, :)

      character(len=*), intent(in) :: cumulus

      real, intent(out)  :: p_cwv_ave(:)
      
      !Local variables:
      integer :: i, k
      real :: dp, trash
      real, dimension(its:ite) :: w_col, w_ccrit, t_troposph

      p_cwv_ave = 0.0
      if (trim(cumulus) /= 'deep') return

      !-- get the pickup of ensemble ave prec, following Neelin et al 2009.
      do i = its, itf
         w_col(i) = 0.
         w_ccrit(i) = 0.
         t_troposph(i) = 0.
         if (ierr(i) /= 0) cycle
         trash = 0.
         do k = kts, ktf
            if (po(i, k) .lt. 200.) exit

            dp = 100.*(po_cup(i, k) - po_cup(i, k + 1))
            trash = trash + dp/c_grav

            w_col(i) = w_col(i) + qo(i, k)*dp/c_grav ! unit mm
            t_troposph(i) = t_troposph(i) + t(i, k)*dp/c_grav
         end do
         !--average temperature
         t_troposph(i) = t_troposph(i)/(1.e-8 + trash)! unit K
         !
         !--- wcrit given by Neelin et al 2009.
         w_ccrit(i) = max(0., 56.2 + 2.1*(t_troposph(i) - 268.)) ! unit mm
         !
         !--- pickup (normalized by the factor 'a')
         !-- <p>=a[(w-w_c)/w_c]**beta, a=0.15, beta=0.23
         !
         p_cwv_ave(i) = (max(0., w_col(i) - p_fpkup*w_ccrit(i))/(1.e-8 + p_fpkup*w_ccrit(i)))**0.23
         p_cwv_ave(i) = max(0., min(1., p_cwv_ave(i)))

         !print*,"NEE=",i,w_col(i),t_troposph(i),w_ccrit(i),p_cwv_ave    (i)
         !print*,"=================================================="
      end do
   end subroutine precipCwvFactor

   !------------------------------------------------------------------------------------
   subroutine getWetbulb(jmin, qo_cup, t_cup, po_cup, q_wetbulb, t_wetbulb)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'getWetbulb' ! subroutine name
   
      !Variables (input, output, inout)
      integer, intent(in) :: jmin

      real, intent(in) :: qo_cup
      real, intent(in) :: t_cup
      real, intent(in) :: po_cup

      real, intent(inout) :: q_wetbulb
      real, intent(inout) :: t_wetbulb
   
      !Local variables:
      real ::  zqp, zcond, zcond1, zcor, zqsat
      real :: psp, pt, pq
      real :: z3es, z4es, z5alcp, zaldcp
      real :: ptare, evap
      real :: foedelta, f0eewmcu, f0ealfcu, f0edemcu, f0eldcpmcu

      !-- for testing
      !              PSP                   TEMP        Q                     ZCOND1
      ! input   85090.0000000000        289.140030372766     1.105078557441815E-002
      ! output  85090.0000000000        287.230570412846     1.181792062536557E-002 -2.761256206705639E-005
      ! PT  = 289.140030372766
      ! PQ  = 1.105078557441815E-002
      ! PSP = 85090.
      !----------------------

      !-- environmental values
      pt = t_cup       ! K
      pq = qo_cup      ! kg/kg
      psp = po_cup*100. ! hPa

      if (pt > c_t00) then
         z3es = c_r3les
         z4es = c_r4les
         z5alcp = c_r5alvcp
         zaldcp = c_ralvdcp
      else
         z3es = c_r3ies
         z4es = c_r4ies
         z5alcp = c_r5alscp
         zaldcp = c_ralsdcp
      end if

      !--- get wet bulb thermo properties --------------------------
      ptare = pt
      zqp = 1.0/psp

      f0ealfcu = min(1.0, ((max(c_rticecu, min(c_t00, ptare)) - c_rticecu)*c_rtwat_rticecu_r)**2)
      f0eewmcu = c_r2es*(f0ealfcu*exp(c_r3les*(ptare - c_t00)/(ptare - c_r4les)) + &
                       (1.0 - f0ealfcu)*exp(c_r3ies*(ptare - c_t00)/(ptare - c_r4ies)))
      zqsat = f0eewmcu*zqp

      zqsat = min(c_max_qsat, zqsat)
      zcor = 1.0/(1.0 - c_retv*zqsat)
      zqsat = zqsat*zcor

      f0edemcu = f0ealfcu*c_r5alvcp*(1.0/(ptare - c_r4les)**2) + &
                 (1.0 - f0ealfcu)*c_r5alscp*(1.0/(ptare - c_r4ies)**2)

      zcond = (pq - zqsat)/(1.0 + zqsat*zcor*f0edemcu)

      zcond = min(zcond, 0.0)

      f0eldcpmcu = f0ealfcu*c_ralvdcp + (1.0 - f0ealfcu)*c_ralsdcp
      pt = pt + f0eldcpmcu*zcond

      pq = pq - zcond

      !--update PTARE
      ptare = pt

      f0ealfcu = min(1.0, ((max(c_rticecu, min(c_t00, ptare)) - c_rticecu)*c_rtwat_rticecu_r)**2)
      f0eewmcu = c_r2es*(f0ealfcu*exp(c_r3les*(ptare - c_t00)/(ptare - c_r4les)) + &
                       (1.0 - f0ealfcu)*exp(c_r3ies*(ptare - c_t00)/(ptare - c_r4ies)))
      zqsat = f0eewmcu*zqp

      zqsat = min(0.5, zqsat)
      zcor = 1.0/(1.0 - c_retv*zqsat)
      zqsat = zqsat*zcor

      f0edemcu = f0ealfcu*c_r5alvcp*(1.0/(ptare - c_r4les)**2) + &
                 (1.0 - f0ealfcu)*c_r5alscp*(1.0/(ptare - c_r4ies)**2)
      zcond1 = (pq - zqsat)/(1.0 + zqsat*zcor*f0edemcu)

      if (zcond == 0.0) zcond1 = min(zcond1, 0.0)
      f0eldcpmcu = f0ealfcu*c_ralvdcp + (1.0 - f0ealfcu)*c_ralsdcp
      pt = pt + f0eldcpmcu*zcond1
      pq = pq - zcond1

      !-- set output --------------------------
      q_wetbulb = pq
      t_wetbulb = pt
      evap = -zcond1 != q_wetbulb-qo_cup, source for water vapor
   end subroutine getWetbulb

   ! ------------------------------------------------------------------------------------
   subroutine cupForcingEns3dShal(itf, ktf, its, ite, kts, kte, dtime, ichoice, ierrc, ierr, klcl, kpbl, kbcon, k22, ktop &
                              ,   xmb, tsur, cape, h_sfc_flux, le_sfc_flux, zws, po, hco, heo_cup, po_cup, t_cup, dhdt &
                              ,   rho, xff_shal2d, xf_dicycle, tke_pbl, wlpool, xf_coldpool)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'cupForcingEns3dShal' ! subroutine name

      real, parameter :: p_k1 = 1.2
      !! tuning numbers for the TKE-based closure for shallow convection
      real, parameter :: p_cloud_area = 0.15
      !! tuning numbers for the TKE-based closure for shallow convection
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte, ichoice

      integer, intent(in) :: klcl(:)
      integer, intent(in) :: kpbl(:)
      integer, intent(in) :: kbcon(:)
      integer, intent(in) :: k22(:)
      integer, intent(in) :: ktop(:)

      real, intent(in) :: dtime
      real, intent(in) :: tsur(:)
      real, intent(in) :: cape(:)
      real, intent(in) :: h_sfc_flux(:)
      real, intent(in) :: le_sfc_flux(:)
      real, intent(in) :: zws(:)
      real, intent(in) :: tke_pbl(:)
      real, intent(in) :: wlpool(:)
      real, intent(in) :: po(:, :)
      real, intent(in) :: hco(:, :)
      real, intent(in) :: heo_cup(:, :)
      real, intent(in) :: po_cup(:, :)
      real, intent(in) :: t_cup(:, :)
      real, intent(in) :: dhdt(:, :)
      real, intent(in) :: rho(:, :)

      integer, intent(inout):: ierr(:)

      real, intent(inout) :: xmb(:)
      real, intent(inout) :: xf_dicycle(:)
      real, intent(inout) :: xf_coldpool(:)

      character(len=128), intent(inout):: ierrc(:)

      real, intent(out)  :: xff_shal2d(:, :)

      !Local variables:
      real, dimension(its:ite)    :: xmbmax
      integer :: i, k, kbase
      real :: blqe, trash, tcold, fin, fsum, efic, thot, dp
      real, dimension(p_shall_closures)  :: xff_shal

      do i = its, itf
         xmb(i) = 0.
         xf_dicycle(i) = 0.
         if (ierr(i) /= 0) cycle

         xmbmax(i) = 100.*(po(i, kbcon(i)) - po(i, kbcon(i) + 1))/(c_grav*dtime)

         !- limiting the mass flux at cloud base
         xmbmax(i) = min(p_xmbmaxshal, xmbmax(i))

         !- cloud base
         kbase = kbcon(i)
         !kbase=klcl(i)

         !--- closure from Grant (2001): ichoice = 1
         xff_shal(1) = .030*zws(i)*rho(i, kpbl(i))
         xff_shal(2) = xff_shal(1)
         xff_shal(3) = xff_shal(1)

         !--- closure from the heat-engine principle : ichoice = 4
         !- Renno and Ingersoll(1996), Souza et al (1999)
         !- get the averaged environment temperature between cloud base
         !- and cloud top
         tcold = 0.
         do k = kbase, ktop(i)
            dp = po_cup(i, k) - po_cup(i, k + 1)
            tcold = tcold + t_cup(i, k)*dp
         end do
         tcold = tcold/(po_cup(i, kbase) - po_cup(i, ktop(i) + 1))

         !-surface temperature
         thot = tsur(i)  ! + ztexec(i)
         !- thermodynamic eficiency
         !efic = max(0.05, (thot-tcold)/thot )
         efic = max(0.0, (thot - tcold)/thot)

         !- total heat flux from surface
         fin = max(0.0, h_sfc_flux(i) + le_sfc_flux(i))

         !--- mass flux at cloud base
         !if(cape(i) > 0.0 .and. h_sfc_flux(i) >0.0 ) then
         if (cape(i) > 0.0) then
            xff_shal(4) = efic*fin/cape(i)
         else
            xff_shal(4) = 0.0
         end if
         xff_shal(5) = xff_shal(4)
         xff_shal(6) = xff_shal(4)

         !--- closure from boundary layer QE (Raymond 1995): ichoice = 7
         blqe = 0.
         trash = 0.
         if (k22(i) .lt. kpbl(i) + 1) then
            do k = kts, kbase
               blqe = blqe + 100.*dhdt(i, k)*(po_cup(i, k) - po_cup(i, k + 1))/c_grav
            end do
            trash = max((hco(i, kbase) - heo_cup(i, kbase)), 1.e1)
            xff_shal(7) = max(0., blqe/trash)
         else
            xff_shal(7) = 0.0
         end if
         xff_shal(8) = xff_shal(7)
         xff_shal(9) = xff_shal(7)

         !--- new closure based on the PBL TKE mean (Zheng et al, 2020 GRL): ichoice = 10
         !-- shallow cumulus active area is for now keept by 0.15 (Zheng 2021 p. commun.)
         !-- k1 is 'slope' of the curve between Wb x (TKE_PBL)**0.5
         !--        and varies between 1.2 (from lidar) to 1.6 (from WRF and SAM models)
         xff_shal(10) = p_cloud_area*rho(i, kbase)*p_k1*sqrt(tke_pbl(i))
         xff_shal(11) = xff_shal(10)
         xff_shal(12) = xff_shal(10)

         !--- store all closures for later.
         xff_shal2d(i, :) = xff_shal(:)

      end do
   end subroutine cupForcingEns3dShal

   ! ------------------------------------------------------------------------------------
   subroutine cupUpLightning(itf, ktf, its, ite, kts, kte, ierr, kbcon, ktop, xland, cape, zo, zo_cup, t_cup, t, tempco &
                           , qrco, po_cup, rho, prec_flx, lightn_dens)
      !! ## Lightning parameterization
      !!
      !! Author: Saulo Freitas [SRF]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>
      !!
      !! Date: 10-Aug-2019
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Lightning parameterization based on:
      !! "A Lightning Parameterization for the ECMWF Integrated Forecasting System"
      !! P. Lopez, 2016 MWR
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'cupUpLightning' ! subroutine name
      real, parameter :: p_v_graup = 3.0  
      !! m/s
      real, parameter :: p_v_snow = 0.5  
      !! m/s
      real, parameter :: p_beta_land = 0.70 
      !! 1
      real, parameter :: p_beta_ocean = 0.45 
      !! 1
      real, parameter :: p_alpha = 37.5 
      !! 1
      real, parameter :: p_t_initial = 0.0 + 273.15 
      !! K
      real, parameter :: p_t_final = -25.+273.15 
      !! K
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte

      integer, intent(in) :: ierr(:)
      integer, intent(in) :: kbcon(:)
      integer, intent(in) :: ktop(:)

      real, intent(in) :: cape(:)
      real, intent(in) :: xland(:)
      real, intent(in) :: po_cup(:, :)
      real, intent(in) :: zo_cup(:, :)
      real, intent(in) :: t_cup(:, :)
      real, intent(in) :: t(:, :)
      real, intent(in) :: tempco(:, :)
      real, intent(in) :: zo(:, :)
      real, intent(in) :: qrco(:, :)
      real, intent(in) :: rho(:, :)
      real, intent(in) :: prec_flx(:, :)

      real, intent(out) :: lightn_dens(:) 
      !! lightning flash density - rate (units: 1/km2/day)

      !Local variables:
      integer :: i, k, k_initial, k_final
      real :: q_r, z_base, beta, prec_flx_fr, dz
      real, dimension(kts:kte) :: p_liq_ice, q_graup, q_snow

      do i = its, itf
         lightn_dens(i) = 0.0
         if (ierr(i) /= 0) cycle

         beta = xland(i)*p_beta_ocean + (1.-xland(i))*p_beta_land

         q_graup(:) = 0.
         q_snow(:) = 0.

         do k = kts, ktop(i)

            p_liq_ice(k) = FractLiqF(tempco(i, k))

            prec_flx_fr = p_liq_ice(k)*prec_flx(i, k)/rho(i, k)

            q_graup(k) = beta*prec_flx_fr/p_v_graup ! - graupel mixing ratio (kg/kg)
            q_snow(k) = (1.-beta)*prec_flx_fr/p_v_snow  ! - snow    mixing ratio (kg/kg)

         end do

         k_initial = minloc(abs(tempco(i, kbcon(i):ktop(i)) - p_t_initial), 1) + kbcon(i) - 1
         k_final = minloc(abs(tempco(i, kbcon(i):ktop(i)) - p_t_final), 1) + kbcon(i) - 1

         q_r = 0.0
         do k = k_initial, k_final
            dz = zo(i, k) - zo(i, k - 1)
            q_r = q_r + dz*rho(i, k)*(q_graup(k)*(qrco(i, k) + q_snow(k)))
            !print*,"qr=",q_r,tempco(i,k)-273.15,k,tempco(i,k)-t_initial
         end do

         z_base = zo_cup(i, kbcon(i))/1000. ! km

         !---
         !--- lightning flash density (units: number of flashes/km2/day) - equation 5
         !--- (to compare with Lopez 2016's results, convert to per year: lightn_dens*365)
         !
         lightn_dens(i) = p_alpha*q_r*sqrt(max(0., cape(i)))*min(z_base, 1.8)**2
         !
      end do
   end subroutine cupUpLightning

   ! ------------------------------------------------------------------------------------
   subroutine getInterp(q_old, t_old, po_cup, q_new, t_new)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'getInterp' ! subroutine name
      !Variables (input, output, inout)
      real, intent(in) :: po_cup 
      !! original

      real, intent(inout) :: q_old
      real, intent(inout) :: t_old
      real, intent(inout) :: q_new
      real, intent(inout) :: t_new 
      !! extrapolated
   
      !Local variables:
      real ::  zqp, zcond1, zcor, zqsat
      real ::  psp, pt, pq, ptare
      real ::  foealfcu, foeewmcu, foedemcu, foeldcpmcu
      integer :: i

      pt = t_old       ! K
      pq = q_old       ! kg/kg
      psp = po_cup*100. ! hPa

      !-- for testing
      !              PSP                   TEMP        Q                     ZCOND1
      ! input    27940.0000000000        236.604976804749       3.220181796223121E-004
      ! output   27940.0000000000        236.361132108860       4.084506812610067E-004
      !  PT  = 236.604976804749      ! K
      !  PQ  = 3.220181796223121E-004       ! kg/kg
      !  PSP = 27940. ! hPa
      !----------------------
      !print*,"1",PSP,PT,PQ

      zqp = 1.0/psp
      do i = 1, 2
         ptare = pt

         foealfcu = min(1.0, ((max(c_rtice, min(c_t00, ptare)) - c_rtice)*c_rtwat_rtice_r)**2)
         foeewmcu = c_r2es*(foealfcu*exp(c_r3les*(ptare - c_t00)/(ptare - c_r4les)) + (1.0 - foealfcu)*exp(c_r3ies*(ptare - c_t00) &
                  / (ptare - c_r4ies)))
         zqsat = foeewmcu*zqp

         !    if(1.0-RETV  *ZQSAT == 0.) then
         !
         !      print*,"ZQSAT=",ZQP,FOEEWMCU,q_old,t_old,po_cup,q_new,t_new
         !3.5491847E-02   46.36052      0.5000000       249.8219
         !  0.2817549      0.5000000       249.8219
         !      call flush(6)
         !      stop 3333
         !    endif

         zcor = 1.0/(1.0 - c_retv*zqsat)
         zqsat = zqsat*zcor

         foedemcu = foealfcu*c_r5alvcp*(1.0/(ptare - c_r4les)**2) + (1.0 - foealfcu)*c_r5alscp*(1.0/(ptare - c_r4ies)**2)

         zcond1 = (pq - zqsat)/(1.0 + zqsat*zcor*foedemcu)

         foeldcpmcu = foealfcu*c_ralvdcp + (1.0 - foealfcu)*c_ralsdcp
         pt = pt + foeldcpmcu*zcond1
         pq = pq - zcond1
      end do
      !-- FINAL --------------------------
      q_new = pq
      t_new = pt
      !print*,"2",PSP,PT,PQ
      !print*,"E",100*(PT-236.361132108860)/236.361132108860,100*(PQ-4.084506812610067E-004)/4.084506812610067E-004
   end subroutine getInterp

   ! -----------------------------------------------------------------------------------------------------------------------
   subroutine sound(part, cumulus, int_time, dtime, ens4, itf, ktf, its, ite, kts, kte, xlats, xlons, jcol, whoami_all &
                    , z, qes, he, hes, t, q, po, z1, psur, zo, qeso, heo, heso, tn, qo, us, vs, omeg, xz, h_sfc_flux, le_sfc_flux &
                    , tsur, dx, stochastic_sig, zws, ztexec, zqexec, xland, kpbl, k22, klcl, kbcon, ktop, aa0, aa1, sig, xaa0, hkb &
                    , xmb, pre, edto, zo_cup, dhdt, rho, zuo, zdo, up_massentro, up_massdetro, outt, outq, outqc, outu, outv)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'sound' ! subroutine name
      real, parameter :: p_latsnd = -10., p_lonsnd = 301., p_deltx = 0.2
      !      real, parameter :: LATSND= -8.72, LONSND= 186.6, DELTX=0.2
   
      !Variables (input, output, inout)
      integer, intent(in) ::ens4, itf, ktf, its, ite, kts, kte, jcol, whoami_all, part

      integer, intent(in) :: k22(:)
      integer, intent(in) :: klcl(:)
      integer, intent(in) :: kbcon(:)
      integer, intent(in) :: ktop(:)

      real, intent(in) :: zo_cup(:, :)
      real, intent(in) :: zuo(:, :)
      real, intent(in) :: zdo(:, :)
      real, intent(in) :: up_massentro(:, :)
      real, intent(in) :: up_massdetro(:, :)
      real, intent(in) :: outt(:, :)
      real, intent(in) :: outq(:, :)
      real, intent(in) :: outqc(:, :)
      real, intent(in) :: outu(:, :)
      real, intent(in) :: outv(:, :)
      real, intent(in) :: stochastic_sig(:)
      real, intent(in) :: xland(:)
      real, intent(in) :: aa0(:)
      real, intent(in) :: aa1(:)
      real, intent(in) :: xaa0(:)
      real, intent(in) :: hkb(:)
      real, intent(in) :: xmb(:)
      real, intent(in) :: pre(:)
      real, intent(in) :: edto(:)
      real, intent(in) :: sig(:)

      real, intent(in) :: int_time, dtime

      character(len=*), intent(in)    :: cumulus

      integer, intent(inout) :: kpbl(:)

      real, intent(inout) :: h_sfc_flux(:)
      real, intent(inout) :: le_sfc_flux(:)
      real, intent(inout) :: tsur(:)
      real, intent(inout) :: dx(:)
      real, intent(inout) :: zws(:)
      real, intent(inout) :: ztexec(:)
      real, intent(inout) :: zqexec(:)
      real, intent(inout) :: xlats(:)      
      real, intent(inout) :: xlons(:)      
      real, intent(inout) :: z1(:)
      real, intent(inout) :: psur(:)
      real, intent(inout) :: qes(:, :)
      real, intent(inout) :: he(:, :)
      real, intent(inout) :: hes(:, :)
      real, intent(inout) :: t(:, :)
      real, intent(inout) :: q(:, :)
      real, intent(inout) :: po(:, :)
      real, intent(inout) :: zo(:, :)
      real, intent(inout) :: heo(:, :)
      real, intent(inout) :: heso(:, :)
      real, intent(inout) :: tn(:, :)
      real, intent(inout) :: qo(:, :)
      real, intent(inout) :: us(:, :)
      real, intent(inout) :: vs(:, :)
      real, intent(inout) :: dhdt(:, :)
      real, intent(inout) :: omeg(:, :, :)

      real, intent(out) :: z(:, :)
      real, intent(out) :: xz(:, :)
      real, intent(out) :: qeso(:, :)
      real, intent(out) :: rho(:, :)

      !---locals
      integer :: i, k, x_kte, x_i, x_jcol, x_k, x_whoami_all
      real :: x_time
      real, dimension(its:ite) :: x_stochastic_sig, x_xland
      
      character(len=200) :: lixo


      if (trim(rundata) == "NONE") then
         if (mod(int_time, 3600.) < dtime) then
            open (15, file="dataLXXX.dat_"//trim(cumulus), status='unknown', position="APPEND")
            if (part == 1) then
               do i = its, itf
                  if (xlats(i) > p_latsnd - p_deltx .and. xlats(i) < p_latsnd + p_deltx) then
                     if (xlons(i) > p_lonsnd - p_deltx .and. xlons(i) < p_lonsnd + p_deltx) then

                        print *, "==============================================="
                        print *, "00>", i, jcol, xlats(i), xlons(i), whoami_all, int_time/3600.
                        call flush (6)

                        write (15, *) "====begin====="
                        write (15, *) "i,jcol,xlats(i),xlons(i),int_time/3600."
                        write (15, *) i, jcol, xlats(i), xlons(i), int_time/3600.

                        write (15, *) "kte,z1(i),psur(i),tsur(i),xland(i)"
                        write (15, *) kte, z1(i), psur(i), tsur(i), xland(i)

                        write (15, *) "h_sfc_flux(i),le_sfc_flux(i),ztexec(i),zqexec(i)"
                        write (15, *) h_sfc_flux(i), le_sfc_flux(i), ztexec(i), zqexec(i)

                        write (15, *) "stochastic_sig(i), dx(i),zws(i),kpbl(i)"
                        write (15, *) stochastic_sig(i), dx(i), zws(i), kpbl(i)

                        write (15, *) "=>k zo po t tn-t q qo-q us vs qes he hes qeso-qes heo-he heso-hes dhdt omeg"
                        do k = kts, kte
                           write (15, 100) k, zo(i, k), po(i, k), t(i, k), tn(i, k) - t(i, k), q(i, k), qo(i, k) - q(i, k) &
                              , us(i, k), vs(i, k), qes(i, k), he(i, k), hes(i, k), qeso(i, k) - qes(i, k), heo(i, k) - he(i, k) &
                              , heso(i, k) - hes(i, k), dhdt(i, k), omeg(i, k, 1:ens4)
                        end do

                     end if
                  end if
               end do
            else
               do i = its, itf
                  if (xlats(i) > p_latsnd - p_deltx .and. xlats(i) < p_latsnd + p_deltx) then
                     if (xlons(i) > p_lonsnd - p_deltx .and. xlons(i) < p_lonsnd + p_deltx) then

                        write (15, *) "====outputs======="
                        write (15, *) "L=", i, jcol, xlats(i), xlons(i), whoami_all
                        write (15, *) "A=", aa0(i), aa1(i), xaa0(i), sig(i)
                        write (15, *) "K=", k22(i), klcl(i), kpbl(i), kbcon(i), ktop(i)
                        write (15, *) "Z=", zo_cup(i, k22(i)) - z1(i), zo_cup(i, klcl(i)) - z1(i), zo_cup(i, kpbl(i)) - z1(i) &
                           , zo_cup(i, kbcon(i)) - z1(i), zo_cup(i, ktop(i)) - z1(i)
                        write (15, *) "H=", hkb(i)/real(c_cp), edto(i)
                        write (15, *) "T=", maxval(outt(i, 1:ktop(i)))*86400., maxval(outq(i, 1:ktop(i)))*86400.*1000., &
                           minval(outt(i, 1:ktop(i)))*86400., minval(outq(i, 1:ktop(i)))*86400.*1000.
                        write (15, *) "P=", xmb(i)*1000., 'g/m2/s', 3600*pre(i), 'mm/h'
                        if (xmb(i) > 0.0) then
                           write (15, *) "=> k zo po zuo,zdo,up_massentro,up_massdetro,outt, outq,outqc,outu,outv"
                           do k = kts, kte
                              write (15, 101) k, zo(i, k), po(i, k), zuo(i, k), zdo(i, k), up_massentro(i, k), up_massdetro(i, k) &
                                            , outt(i, k)*86400., outq(i, k)*86400.*1000., outqc(i, k)*86400.*1000., outu(i, k) &
                                            *86400., outv(i, k)*86400.

                           end do
                        end if
                        write (15, *) "=====end=========="
                     end if
                  end if
               end do
            end if
            close (15)
         end if
      else
         if (part == 1) then
            open (15, file=trim(rundata), status='old')
            i = 1
            read (15, *) lixo
            read (15, *) lixo
            read (15, *) x_i, x_jcol, xlats(i), xlons(i), x_time
            read (15, *) lixo
            read (15, *) x_kte, z1(i), psur(i), tsur(i), x_xland(i)
            !-- check
            if (x_kte .ne. kte) stop " X_kte .ne. kte "
            read (15, *) lixo
            read (15, *) h_sfc_flux(i), le_sfc_flux(i), ztexec(i), zqexec(i)
            read (15, *) lixo
            read (15, *) x_stochastic_sig(i), dx(i), zws(i), kpbl(i)
            read (15, *) lixo
            do k = kts, kte
               read (15, 100) x_k, zo(i, k), po(i, k), t(i, k), tn(i, k), q(i, k), qo(i, k), us(i, k), vs(i, k) , qes(i, k) &
                            , he(i, k), hes(i, k), qeso(i, k), heo(i, k), heso(i, k), dhdt(i, k), omeg(i, k, 1:ens4)
            end do
            close (15)
            !---settings
            tn(i, :) = t(i, :) + tn(i, :) ! input is delta(T)
            qo(i, :) = q(i, :) + qo(i, :) ! input is delta(Q)
            qeso(i, :) = qes(i, :) + qeso(i, :) ! input is delta(Q)
            heo(i, :) = he(i, :) + heo(i, :) ! input is delta(H)
            heso(i, :) = hes(i, :) + heso(i, :) ! input is delta(HO)
            xz(i, :) = zo(i, :)
            z(i, :) = zo(i, :)
            rho(i, :) = 1.e2*po(i, :)/(c_rgas*t(i, :))
         else
            do i = its, itf
               if (xlats(i) > p_latsnd - p_deltx .and. xlats(i) < p_latsnd + p_deltx) then
                  if (xlons(i) > p_lonsnd - p_deltx .and. xlons(i) < p_lonsnd + p_deltx) then

                     print *, "====outputs======="
                     print *, "A=", aa0(i), aa1(i), xaa0(i), sig(i)
                     print *, "K=", k22(i), klcl(i), kpbl(i), kbcon(i), ktop(i)
                     print *, "Z=", zo_cup(i, k22(i)) - z1(i), zo_cup(i, klcl(i)) - z1(i), zo_cup(i, kpbl(i)) - z1(i) &
                        , zo_cup(i, kbcon(i)) - z1(i), zo_cup(i, ktop(i)) - z1(i)
                     print *, "H=", hkb(i)/real(c_cp), edto(i)
                     print *, "T=", maxval(outt(i, 1:ktop(i)))*86400., maxval(outq(i, 1:ktop(i)))*86400.*1000. &
                                  , minval(outt(i, 1:ktop(i)))*86400., minval(outq(i, 1:ktop(i)))*86400.*1000.
                     print *, "P=", xmb(i)*1000., 'g/m2/s', 3600*pre(i), 'mm/h'
                     if (xmb(i) > 0.0) then
                        print *, "=> k zo po zuo,zdo,up_massentro,up_massdetro,outt, outq,outqc,outu,outv"
                        do k = kts, kte
                           write (*, 101) k, zo(i, k), po(i, k) &
                              , zuo(i, k), zdo(i, k), up_massentro(i, k), up_massdetro(i, k), outt(i, k)*86400. &
                              , outq(i, k)*86400.*1000., outqc(i, k)*86400.*1000., outu(i, k)*86400., outv(i, k)*86400.
                        end do
                     end if
                  end if
               end if
            end do
         end if
      end if
100   format(1x, i4, 16e16.8)
101   format(1x, i4, 11e16.8)

   end subroutine sound

   ! ------------------------------------------------------------------------------------
   subroutine cloudDissipation(cumulus, itf, ktf, its, ite, kts, kte, ierr, kbcon, ktop, dtime, xmb, xland &
                           ,   qo_cup, qeso_cup, po_cup, outt, outq, outqc, zuo, vvel2d, rho_hydr &
                           ,   qrco, sig, tempco, qco, tn_cup, heso_cup, zo)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'cloudDissipation' ! subroutine name

      real, parameter :: p_cloud_lifetime = 1800.
      integer, parameter :: p_versionx = 2
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte

      integer, intent(in) :: ierr(:)
      integer, intent(in) :: kbcon(:)
      integer, intent(in) :: ktop(:)

      real, intent(in) :: xmb(:)
      real, intent(in) :: xland(:)
      real, intent(in) :: sig(:)
      real, intent(in) :: po_cup(:,:)
      real, intent(in) :: qo_cup(:,:)
      real, intent(in) :: qeso_cup(:,:)
      real, intent(in) :: zuo(:,:)
      real, intent(in) :: vvel2d(:,:)
      real, intent(in) :: rho_hydr(:,:)
      real, intent(in) :: tempco(:,:)
      real, intent(in) :: qco(:,:)
      real, intent(in) :: tn_cup(:,:)
      real, intent(in) :: heso_cup(:,:)
      real, intent(in) :: zo(:,:)
      real, intent(in) :: dtime

      character(len=*), intent(in) :: cumulus

      real, intent(inout) :: outt(:,:)
      real, intent(inout) :: outq(:,:)
      real, intent(inout) :: outqc(:,:)
      real, intent(inout) :: qrco(:,:)

      !Local variables:
      integer :: i, k
      real :: del_t, del_q, dp, frh
      real :: qrc_diss, fractional_area, outqc_diss, outq_mix, outt_diss, outt_mix, tempx, qvx

      do i = its, itf

         if (ierr(i) /= 0) cycle

         do k = ktop(i), kbcon(i), -1

            !--- cloud liq/ice remained in the convection plume
            qrc_diss = max(0., qrco(i, k) - outqc(i, k)*dtime)

            !dp  = 100.*(po_cup(i,k)-po_cup(i,k+1))

            !--- get relative humidity
            frh = 0. !min(qo_cup(i,k)/qeso_cup(i,k),1.)

            !--- estimation of the fractional area
            fractional_area = (xmb(i)/sig(i))*zuo(i, k)/(rho_hydr(i, k)*vvel2d(i, k))

            !--- source of enviroment moistening/cooling due to the 'remained' cloud dissipation into it.
            outqc_diss = (qrc_diss*(1.-frh))/p_cloud_lifetime

            if (p_versionx == 1 .or. p_coupl_mphysics .eqv. .false.) then

               outt_diss = -outqc_diss*(real(c_alvl)/real(c_cp)) !--- cooling

               !--- source of enviroment moistening/warming due to the 'remained' in-cloud water vapor mixing into it.
               !  qvx   = qco   (i,k)
               !  tempx = tempco(i,k)
               qvx = qeso_cup(i, k)
               tempx = (heso_cup(i, k) - c_grav*zo(i, k) - real(c_alvl)*qeso_cup(i, k))/real(c_cp)

               outq_mix = (qvx - qo_cup(i, k))/p_cloud_lifetime

               outt_mix = (tempx - tn_cup(i, k))/p_cloud_lifetime

               !-- feedback
               del_q = (outqc_diss + outq_mix)*USE_CLOUD_DISSIPATION*fractional_area ! units: kg[water]/kg[air}/sec
               del_t = (outt_diss + outt_mix)*USE_CLOUD_DISSIPATION*fractional_area ! units: K/sec

               outq(i, k) = outq(i, k) + del_q
               outt(i, k) = outt(i, k) + del_t

            else

               outqc(i, k) = outqc(i, k) + outqc_diss*fractional_area*USE_CLOUD_DISSIPATION

            end if

            !print*,"diss2=",k,real(outqc_diss*86400.*1000),real(sqrt(1.-sig(i)),4),real( fractional_area*100.,4)

            qrco(i, k) = max(0., qrco(i, k) - outqc_diss*USE_CLOUD_DISSIPATION*fractional_area*dtime)
            !if(qrco (i,k) <0.) print*,"qrc<0",trim(cumulus),qrco(i,k)

         end do
      end do

   end subroutine cloudDissipation

   ! -------------------------------------------------------------------------------------------------------
   subroutine gfConvParInit(mynum)
      !! ## read the namelist
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! read the namelist
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'gfConparInit' ! subroutine name
   
      !Variables (input, output, inout)
      integer, intent(in) :: mynum
      !!Number of current processor
      
      !Local variables:
      character(len=64) :: fn_nml = 'GF_ConvPar_nml'
      logical :: exists
      integer :: nlunit = 4
      
      !Code:
      
      namelist /GF_NML/ ICUMULUS_GF, CLOSURE_CHOICE, USE_SCALE_DEP, DICYCLE &
         , USE_TRACER_TRANSP, USE_TRACER_SCAVEN, USE_FLUX_FORM, USE_TRACER_EVAP, DOWNDRAFT, USE_FCT &
         , USE_REBCB, VERT_DISCR, SATUR_CALC, CLEV_GRID, APPLY_SUB_MP, ALP1 &
         , SGS_W_TIMESCALE, LIGHTNING_DIAG, AUTOCONV, BC_METH, OVERSHOOT, USE_WETBULB &
         , C1, C0_DEEP, QRC_CRIT, LAMBAU_DEEP, LAMBAU_SHDN, C0_MID &
         , CUM_MAX_EDT_LAND, CUM_MAX_EDT_OCEAN, CUM_HEI_DOWN_LAND &
         , CUM_HEI_DOWN_OCEAN, CUM_HEI_UPDF_LAND, CUM_HEI_UPDF_OCEAN &
         , CUM_ENTR_RATE, TAU_DEEP, TAU_MID &
         , USE_MOMENTUM_TRANSP, MOIST_TRIGGER, FRAC_MODIS &
         , CUM_USE_EXCESS, CUM_AVE_LAYER, ADV_TRIGGER, USE_SMOOTH_PROF &
         , USE_CLOUD_DISSIPATION, USE_SMOOTH_TEND, USE_GUSTINESS, USE_RANDOM_NUM &
         , DCAPE_THRESHOLD, BETA_SH, C0_SHAL, USE_LINEAR_SUBCL_MF, LIQ_ICE_NUMBER_CONC &
         , ALPHA_ADV_TUNING, CAP_MAXS, SIG_FACTOR, CUM_FADJ_MASSFLX, LCL_TRIGGER &
         , RH_DICYCLE, CUM_T_STAR, CONVECTION_TRACER, TAU_OCEA_CP, TAU_LAND_CP &
         , USE_MEMORY, ADD_COLDPOOL_PROP, MX_BUOY1, MX_BUOY2, MAX_TQ_TEND, CUM_ZUFORM &
         , ADD_COLDPOOL_CLOS, ADD_COLDPOOL_DIFF

      inquire (file=trim(fn_nml), exist=exists)
      if (.not. exists) then
         write (6, *) 'GF_convpar_nml :: namelist file: ', trim(fn_nml), ' does not exist'
         stop 31415
      else
         open (nlunit, file=fn_nml, status='old', form='formatted')
         read (nlunit, nml=GF_NML)
         close (nlunit)
      end if
      if (mynum == 1) then
         !- print the namelist
         print *, "           "
         print *, "------------- GF ConvPar namelist -------------"
         print *, "!---- the main controls"
         print *, 'icumulus_gf        ', ICUMULUS_GF
         print *, 'cum_entr           ', real(CUM_ENTR_RATE, 4)
         print *, 'closure_choice     ', CLOSURE_CHOICE
         print *, 'use_scale_dep      ', USE_SCALE_DEP
         print *, 'sig_factor         ', real(SIG_FACTOR, 4)
         print *, 'dicycle            ', DICYCLE
         print *, 't_star             ', real(CUM_T_STAR, 4)
         print *, 'rh_dicycle         ', RH_DICYCLE
         print *, 'alpha_adv_tuning   ', real(ALPHA_ADV_TUNING, 4)
         print *, 'cap_maxs           ', real(CAP_MAXS, 4)
         print *, 'moist_trigger      ', MOIST_TRIGGER
         print *, 'adv_trigger        ', ADV_TRIGGER
         print *, 'lcl_trigger        ', LCL_TRIGGER
         print *, 'dcape_threshold    ', real(DCAPE_THRESHOLD, 4)
         print *, 'tau_deep,tau_mid   ', real(TAU_DEEP, 4), real(TAU_MID, 4)
         print *, 'SGS_W_TIMESCALE    ', SGS_W_TIMESCALE
         print *, 'CONVECTION_TRACER  ', CONVECTION_TRACER
         print *, 'ADD_COLDPOOL_PROP  ', ADD_COLDPOOL_PROP
         print *, 'ADD_COLDPOOL_CLOS  ', ADD_COLDPOOL_CLOS
         print *, 'ADD_COLDPOOL_DIFF  ', ADD_COLDPOOL_DIFF
         print *, 'tau_ocea_cp        ', TAU_OCEA_CP
         print *, 'tau_land_cp        ', TAU_LAND_CP
         print *, 'mx_buoy1 - kJ/kg   ', MX_BUOY1*1.e-3
         print *, 'mx_buoy2 - kJ/kg   ', MX_BUOY2*1.e-3
         print *, 'USE_MEMORY         ', USE_MEMORY

         print *, '!--- controls rainfall evaporation'
         print *, 'USE_REBCB          ', USE_REBCB
         print *, 'downdraft          ', DOWNDRAFT
         print *, 'max_edt_land       ', real(CUM_MAX_EDT_LAND, 4)
         print *, 'max_edt_ocean      ', real(CUM_MAX_EDT_OCEAN, 4)

         print *, '!---- boundary condition specification'
         print *, 'BC_METH            ', BC_METH
         print *, 'cum_use_excess     ', CUM_USE_EXCESS
         print *, 'cum_ave_layer      ', real(CUM_AVE_LAYER, 4)

         print *, '!---- for mass flux profiles - (deep ,shallow ,congestus)'
         print *, 'CUM_ZUFORM         ', CUM_ZUFORM
         print *, 'hei_down_land      ', real(CUM_HEI_DOWN_LAND, 4)
         print *, 'hei_down_ocean     ', real(CUM_HEI_DOWN_OCEAN, 4)
         print *, 'hei_updf_land      ', real(CUM_HEI_UPDF_LAND, 4)
         print *, 'hei_updf_ocean     ', real(CUM_HEI_UPDF_OCEAN, 4)
         print *, 'beta_sh            ', real(BETA_SH, 4)
         print *, 'use_linear_subcl_mf', USE_LINEAR_SUBCL_MF
         print *, 'use_smooth_prof    ', USE_SMOOTH_PROF
         print *, 'use_smooth_tend    ', USE_SMOOTH_TEND
         print *, 'use_random_num     ', USE_RANDOM_NUM

         print *, '!---- the cloud microphysics'
         print *, 'AUTOCONV           ', AUTOCONV
         print *, 'C0_DEEP            ', real(C0_DEEP, 4)
         print *, 'C0_MID             ', real(C0_MID, 4)
         print *, 'C0_SHAL            ', real(C0_SHAL, 4)
         print *, 'c1                 ', real(C1, 4)
         print *, 'QRC_CRIT           ', real(QRC_CRIT, 4)

         print *, '!--- for momentum transport'
         PRINT *, 'USE_MOMENTUM_TRANS ', USE_MOMENTUM_TRANSP
         print *, 'lambau_deep        ', real(LAMBAU_DEEP, 4)
         print *, 'lambau_shdn        ', real(LAMBAU_SHDN, 4)

         print *, '!--- for tracer transport'
         print *, 'USE_TRACER_TRANSP  ', USE_TRACER_TRANSP
         print *, 'USE_TRACER_SCAVEN  ', USE_TRACER_SCAVEN
         print *, 'USE_FLUX_FORM      ', USE_FLUX_FORM
         print *, 'use_fct            ', USE_FCT
         print *, 'use_tracer_evap    ', USE_TRACER_EVAP
         print *, 'apply_sub_mp       ', APPLY_SUB_MP
         print *, 'ALP1               ', real(ALP1, 4)

         print *, '!---- couplings w/ other parameterizations'
         print *, 'LIGHTNING_DIAG     ', LIGHTNING_DIAG
         print *, 'OVERSHOOT          ', real(OVERSHOOT, 4)
         print *, 'liq_ice_number_conc', LIQ_ICE_NUMBER_CONC
         print *, 'use_gustiness      ', USE_GUSTINESS

         print *, '!----misc controls'
         print *, 'frac_modis         ', FRAC_MODIS
         print *, 'use_cloud_dissipat ', real(USE_CLOUD_DISSIPATION, 4)
         print *, 'USE_WETBULB        ', USE_WETBULB
         print *, 'CLEV_GRID          ', CLEV_GRID
         print *, 'vert_discr         ', VERT_DISCR
         print *, 'satur_calc         ', SATUR_CALC
         print *, 'max_tq_tend        ', real(MAX_TQ_TEND, 4)
         print *, 'cum_fadj_massflx   ', real(CUM_FADJ_MASSFLX, 4)
         print *, "========================================================================"
         call flush (6)
      end if
   
   end subroutine gfConvParInit

   ! ------------------------------------------------------------------------------------
   subroutine getLiqIceNumberConc(itf, ktf, its, ite, kts, kte, ierr, ktop, dtime, rho, outqc, tempco, outnliq, outnice)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'getLiqIceNumberConc' ! subroutine name
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte

      integer, intent(in) :: ierr(:) 
      integer, intent(in) :: ktop(:)

      real, intent(in) :: outqc(:,:)
      real, intent(in) :: tempco(:,:)
      real, intent(in) :: rho(:,:)
      real, intent(in) :: dtime

      real, intent(out) :: outnliq(:,:)
      real, intent(out) :: outnice(:,:)

      !Local variables:
      integer :: i, k
      real :: fr, tqliq, tqice, dtinv
      real, dimension(its:ite, kts:kte) :: nwfa   
      !! in the future set this as NCPL
      real, dimension(its:ite, kts:kte) :: nifa   
      !! in the future set this as NCPI

      nwfa(:, :) = 99.e7  ! in the future set this as NCPL
      nifa(:, :) = 0.     ! in the future set this as NCPI
      dtinv = 1./dtime
      do i = its, itf
         if (ierr(i) /= 0) cycle

         do k = kts, ktop(i) + 1

            fr = FractLiqF(tempco(i, k))
            tqliq = dtime*outqc(i, k)*rho(i, k)*fr
            tqice = dtime*outqc(i, k)*rho(i, k)*(1.-fr)

            outnice(i, k) = max(0.0, MakeIceNumber(tqice, tempco(i, k))/rho(i, k))
            outnliq(i, k) = max(0.0, MakeDropletNumber(tqliq, nwfa(i, k))/rho(i, k))

         end do
         !-- convert in tendencies
         outnice = outnice*dtinv ! unit [1/s]
         outnliq = outnliq*dtinv ! unit [1/s]
         !--- for update
         ! nwfa =nwfa + outnliq*dtime
         ! nifa =nifa + outnice*dtime

      end do

   end subroutine getLiqIceNumberConc

   ! -----------------------------------------------------------------------
   elemental function MakeIceNumber(q_ice, temp) result(ice_number)
      !! ## number_concentrations
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! —
      !! **Full description**:
      !!
      !! -! ---- module_mp_thompson_make_number_concentrations
      !! -!  Developed by H. Barnes @ NOAA/OAR/ESRL/GSL Earth Prediction Advancement Division
      !!
      !!     Q_ice              is cloud ice mixing ratio, units of kg/m3
      !!     Q_cloud            is cloud water mixing ratio, units of kg/m3
      !!     Q_rain             is rain mixing ratio, units of kg/m3
      !!     temp               is air temperature in Kelvin
      !!     make_IceNumber     is cloud droplet number mixing ratio, units of number per m3
      !!     MakeDropletNumber is rain number mixing ratio, units of number per kg of m3
      !!     make_RainNumber    is rain number mixing ratio, units of number per kg of m3
      !!     qnwfa              is number of water-friendly aerosols in number per kg
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !!  
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'MakeIceNumber' ! function name

      real, parameter:: c_ice_density = 890.0
   
      !Variables (input):
      real, intent(in) :: q_ice
      real, intent(in) :: temp
   
      !Local variables:
      real :: ice_number ! output

      integer :: idx_rei
      real :: corr, reice, deice
      double precision :: lambda

      !+---+-----------------------------------------------------------------+
      !..Table of lookup values of radiative effective radius of ice crystals
      !.. as a function of Temperature from -94C to 0C.  Taken from WRF RRTMG
      !.. radiation code where it is attributed to Jon Egill Kristjansson
      !.. and coauthors.
      !+---+-----------------------------------------------------------------+

      real, dimension(95), parameter:: p_retab = (/ &
                                       5.92779, 6.26422, 6.61973, 6.99539, 7.39234, &
                                       7.81177, 8.25496, 8.72323, 9.21800, 9.74075, 10.2930, &
                                       10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319, &
                                       15.0351, 15.8799, 16.7674, 17.6986, 18.6744, 19.6955, &
                                       20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125, &
                                       27.7895, 28.6450, 29.4167, 30.1088, 30.7306, 31.2943, &
                                       31.8151, 32.3077, 32.7870, 33.2657, 33.7540, 34.2601, &
                                       34.7892, 35.3442, 35.9255, 36.5316, 37.1602, 37.8078, &
                                       38.4720, 39.1508, 39.8442, 40.5552, 41.2912, 42.0635, &
                                       42.8876, 43.7863, 44.7853, 45.9170, 47.2165, 48.7221, &
                                       50.4710, 52.4980, 54.8315, 57.4898, 60.4785, 63.7898, &
                                       65.5604, 71.2885, 75.4113, 79.7368, 84.2351, 88.8833, &
                                       93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424, &
                                       124.954, 130.630, 136.457, 142.446, 148.608, 154.956, &
                                       161.503, 168.262, 175.248, 182.473, 189.952, 197.699, &
                                       205.728, 214.055, 222.694, 231.661, 240.971, 250.639/)

      if (q_ice == 0) then
         ice_number = 0
         return
      end if

      !+---+-----------------------------------------------------------------+
      !..From the model 3D temperature field, subtract 179K for which
      !.. index value of retab as a start.  Value of corr is for
      !.. interpolating between neighboring values in the table.
      !+---+-----------------------------------------------------------------+

      idx_rei = int(temp - 179.)
      idx_rei = min(max(idx_rei, 1), 94)
      corr = temp - int(temp)
      reice = p_retab(idx_rei)*(1.-corr) + p_retab(idx_rei + 1)*corr
      deice = 2.*reice*1.e-6

      !+---+-----------------------------------------------------------------+
      !..Now we have the final radiative effective size of ice (as function
      !.. of temperature only).  This size represents 3rd moment divided by
      !.. second moment of the ice size distribution, so we can compute a
      !.. number concentration from the mean size and mass mixing ratio.
      !.. The mean (radiative effective) diameter is 3./Slope for an inverse
      !.. exponential size distribution.  So, starting with slope, work
      !.. backwords to get number concentration.
      !+---+-----------------------------------------------------------------+

      lambda = 3.0/deice
      ice_number = q_ice*lambda*lambda*lambda/(c_pi*c_ice_density)

      !+---+-----------------------------------------------------------------+
      !..Example1: Common ice size coming from Thompson scheme is about 30 microns.
      !.. An example ice mixing ratio could be 0.001 g/kg for a temperature of -50C.
      !.. Remember to convert both into MKS units.  This gives N_ice=357652 per kg.
      !..Example2: Lower in atmosphere at T=-10C matching ~162 microns in retab,
      !.. and assuming we have 0.1 g/kg mixing ratio, then N_ice=28122 per kg,
      !.. which is 28 crystals per liter of air if the air density is 1.0.
      !+---+-----------------------------------------------------------------+
   end function MakeIceNumber

   ! ------------------------------------------------------------------------------
   elemental function MakeDropletNumber(q_cloud, qnwfa) result(droplet_number)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! —
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !!
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'MakeDropletNumber' ! function name
   
      real, parameter:: c_am_r = c_pi*1000./6.
      real, dimension(15), parameter:: c_g_ratio = (/24, 60, 120, 210, 336, 504, 720, 990, 1320, 1716, 2184, 2730, 3360, 4080 &
                                                   , 4896/)
      !Variables (input):
      real, intent(in):: q_cloud
      real, intent(in):: qnwfa
   
      !Local variables:
      real :: droplet_number ! out
      double precision:: lambda, qnc
      real:: q_nwfa, x1, xDc
      integer:: nu_c

      if (Q_cloud == 0) then
         droplet_number = 0
         return
      end if
      
      q_nwfa = max(99.e6, min(qnwfa, 5.e10))
      nu_c = max(2, min(nint(2.5e10/q_nwfa), 15))

      x1 = max(1., min(q_nwfa*1.e-9, 10.)) - 1.
      xDc = (30.-x1*20./9.)*1.e-6

      lambda = (4.0d0 + nu_c)/xDc
      qnc = Q_cloud/c_g_ratio(nu_c)*lambda*lambda*lambda/c_am_r
      droplet_number = SNGL(qnc)

      return
   end function MakeDropletNumber

   ! ----------------------------------------------------------------------------------------
   pure function IntFuncGamma(x, y) result(z)
      !! ## ???
      !!
      !! Author: Demerval Moreira [DSM]
      !!
      !! E-mail: <mailto:demerval.moreira@unesp.br>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! —
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !!
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'IntFuncGamma' ! function name
   
      !Variables (input):
      real, intent(in) :: x
      real, intent(in) :: y
   
      !Local variables:
      real :: z ! output
   
      !Code:
      z = x**(y - 1.0)*exp(-x)

   end function IntFuncGamma

   ! ---------------------------------------------------------------------------------------------
   function GammaBrams(a) result(g)
      !! ## ???
      !!
      !! Author: Demerval Moreira [DSM]
      !!
      !! E-mail: <mailto:demerval.moreira@unesp.br>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! —
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !!
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'GammaBrams' ! funtion name

      real, parameter :: p_small = 1.0e-4
      integer, parameter :: p_points = 100000
   
      !Variables (input):
      real, intent(in) :: a
   
      !Local variables:
      real :: g  !Output

      real :: infty, dx, p, sp(2, p_points), x
      integer :: i
      logical :: correction

      x = a

      correction = .false.
      ! value with x<1 gives \infty, so we use
      ! \Gamma(x+1) = x\Gamma(x)
      ! to avoid the problem
      if (x < 1.0) then
         correction = .true.
         x = x + 1
      end if

      ! find a "reasonable" infinity...
      ! we compute this integral indeed
      ! \int_0^M dt t^{x-1} e^{-t}
      ! where M is such that M^{x-1} e^{-M} ≤ \epsilon
      infty = 1.0e4
      do while (IntFuncGamma(infty, x) > p_small)
         infty = infty*10.0
      end do

      ! using simpson
      dx = infty/real(p_points)
      sp = 0.0
      forall (i=1:p_points/2 - 1) sp(1, 2*i) = IntFuncGamma(2.0*(i)*dx, x)
      forall (i=1:p_points/2) sp(2, 2*i - 1) = IntFuncGamma((2.0*(i) - 1.0)*dx, x)
      g = (IntFuncGamma(0.0, x) + 2.0*sum(sp(1, :)) + 4.0*sum(sp(2, :)) + &
           IntFuncGamma(infty, x))*dx/3.0

      if (correction) g = g/a

   end function GammaBrams

   ! ------------------------------------------------------------------------------------------
   function Ran1(idum) result(random_number)
      !! ## Random number generator
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! —
      !! **Full description**:
      !!
      !! This is contributed code standardized by Yong Wang
      !! Random number generator taken from Press et al.
      !!
      !! Returns numbers in the range 0-->1
      !!
      !! Their description...
      !! "Minimal" random number generator of Park and Miller with Bays-Durham
      !! shuffle and added safeguards. Returns a uniform deviate between 0.0 and 1.0
      !! (exclusive of the endpoint values). Call with idum a negative integer to
      !! initialize; thereafter, do not alter idum between successive calls in a
      !! sequence. RNMX should approximate the largest floating value that is less
      !! than 1.
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !!  
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'Ran1' ! function name

      integer(kind = kind_ib), parameter:: p_ntab = 32
      integer(kind = kind_ib), parameter:: p_iq = 127773
      integer(kind = kind_ib), parameter:: p_ia = 16807
      integer(kind = kind_ib), parameter:: p_ir = 2836
      integer(kind = kind_ib), parameter:: p_im = 2147483647
      integer(kind = kind_ib), parameter:: p_ndiv = 1 + (p_im - 1)/p_ntab
      real(kind = kind_rb), parameter:: p_am = 1.0/p_im
      real(kind = kind_rb), parameter:: p_eps = 1.2e-7
      real(kind = kind_rb), parameter:: p_rnmx = 1.0 - p_eps
   
      !Variables (input):
      integer(kind = kind_ib), intent(inout):: idum
   
      !Local variables:
      real(kind = kind_rb) :: random_number !output

      integer(kind = kind_ib):: iy
      integer(kind = kind_ib), dimension(p_ntab):: iv
      !save iv,iy
      data iv/p_ntab*0/, iy/0/
      integer(kind = kind_ib):: j, k

      if (idum .le. 0 .or. iy .eq. 0) then
         ! initialize
         idum = max(-idum, 1)
         do j = p_ntab + 8, 1, -1
            k = idum/p_iq
            idum = p_ia*(idum - k*p_iq) - p_ir*k
            if (idum .lt. 0) idum = idum + p_im
            if (j .le. p_ntab) iv(j) = idum
         end do
         iy = iv(1)
      end if
      !
      k = idum/p_iq
      ! compute idum = mod(ia*idum,im) without overflows by schrage's method
      idum = p_ia*(idum - k*p_iq) - p_ir*k
      if (idum .lt. 0) idum = idum + p_im
      ! j will be in the range 1-->ntab
      j = 1 + iy/p_ndiv
      ! output previously stored value and refill the shuffle table
      iy = iv(j)
      iv(j) = idum
      random_number = min(p_am*iy, p_rnmx)

   end function Ran1

  ! ----------------------------------------------------------------------
   subroutine getDelmix(cumulus, kts, kte, ktf, xland, subcl_level, po, ain, aout)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'getDelmix' ! subroutine namea
   
      !Variables (input, output, inout)
      integer, intent(in) :: kts, kte, ktf, subcl_level
      
      real, intent(in) :: ain(:)
      real, intent(in) :: po(:)
      real, intent(in) :: xland
      
      character(len=*), intent(in) :: cumulus

      real, intent(inout) :: aout(:)

      !Local variables:
      integer :: k
      real :: x1, x2, dp, del, qc

      qc = aout(kts)
 
      x2 = 0.
      x1 = 0.
      do k = kts, subcl_level
         dp = po(k + 1) - po(k)
         x2 = x2 + dp
         x1 = x1 + dp*ain(k)
      end do
      del = abs(qc - x1/(x2 + 1.e-12))
      aout(kts:subcl_level) = ain(kts:subcl_level) + del
 
   end subroutine getDelmix

   ! ----------------------------------------------------------------------------------------------
   subroutine getQadv(cumulus, itf, ktf, its, ite, kts, kte, ierr, dt, q, qo, qo_adv, po, po_cup, qeso, q_adv, col_sat_adv &
                    , alpha_adv, tau_bl, zo_cup, kbcon, ktop)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'getQadv' ! subroutine name

      real, parameter :: p_ptop = 60.
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte

      integer, intent(in) :: ierr(:)
      integer, intent(in) :: kbcon(:)
      integer, intent(in) :: ktop(:)

      real, intent(in) :: q(:,:)
      real, intent(in) :: qo(:,:)
      real, intent(in) :: qo_adv(:,:)
      real, intent(in) :: po_cup(:,:)
      real, intent(in) :: qeso(:,:)
      real, intent(in) :: po(:,:)
      real, intent(in) :: zo_cup(:,:)
      real, intent(in) :: tau_bl(:)
      real, intent(in) :: dt

      character(len=*), intent(in) :: cumulus

      real, intent(inout) :: q_adv(:)
      real, intent(inout) :: col_sat_adv(:)
      real, intent(inout) :: alpha_adv(:)
   
      !Local variables:
      integer :: i, k
      real :: dp, layer, H_cloud, dz

      !-- get the advective moisture tendency scaled with the relative humidity
      !--  Q_adv = integral( q/q*  DQv/Dt_adv dp), see Eq 1 Becker et al(2021 QJRMS)
      !-- units here are "J m^-3" _or_  "J kg^-1"

      do i = its, itf
         col_sat_adv(i) = 0.   !check if it needs be inout, perhavps only local var

         if (ierr(i) /= 0) cycle

         alpha_adv(i) = ALPHA_ADV_TUNING
         layer = 0.

         loopN: do k = kts, ktf
            if (po(i, k) < p_ptop) exit loopN

            !dp=100.*(po_cup(i,k+1)-po_cup(i,k)) ! dp < 0.
            dz = zo_cup(i, k + 1) - zo_cup(i, k)  ! dz > 0

            !-- integral over dp
            !Q_adv(i) = Q_adv(i) + dp*(qo(i,k)/qeso(i,k))*(qo_adv(i,k)-q(i,k))/dt

            !-- integral over dz
            q_adv(i) = q_adv(i) + dz*(qo(i, k)/qeso(i, k))*(qo_adv(i, k) - q(i, k))/dt

            col_sat_adv(i) = col_sat_adv(i) + dz*qo(i, k)/qeso(i, k)

            layer = layer + dz

         end do loopN
         !-- get the column-average saturation fraction
         col_sat_adv(i) = col_sat_adv(i)/(1.e-8 + layer)

         !--check if the col-ave saturation fraction is over the threshold
         if (col_sat_adv(i) > col_sat_adv_threshold) then

            alpha_adv(i) = 0.0
            cycle

         end if

         !-- check if cloud top _OR_cloud layer   !<<<< check this
         H_cloud = zo_cup(i, ktop(i)) - zo_cup(i, kbcon(i))

         !-- convert Q_adv to units as in Eq (1) => J m^-3
         !Q_adv(i) = - Q_adv(i) * tau_bl(i) * xlv / (g * H_cloud)

         !-- convert Q_adv to units as in cloud work function => J kg-1
         q_adv(i) = q_adv(i)*tau_bl(i)*real(c_alvl)/(H_cloud)

         !if(abs(q_adv(i))>1.) print*,"Qadv=",i,q_adv(i),Q_adv_dz(i)call flush(6)
      end do

   end subroutine getQadv

   ! ----------------------------------------------------------------------
   subroutine rhControls(itf, ktf, its, ite, kts, kte, ierr, temp, po, qo, qeso, po_cup, cumulus, rh_entr_factor, rh_dicycle_fct &
                       , entr_rate_input, entr_rate, xlons, dt)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'rhControls' ! subroutine name

      real, parameter :: p_ref_local_time = 8., p_ftun3 = 0.25
      logical, parameter :: p_free_troposphere = .true.

      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, ite, kts, kte

      integer, intent(in) :: ierr(:)

      real, intent(in) :: temp(:,:)
      real, intent(in) :: po(:,:)
      real, intent(in) :: qo(:,:)
      real, intent(in) :: po_cup(:,:)
      real, intent(in) :: qeso(:,:)
      real, intent(in) :: xlons(:)
      real, intent(in) :: entr_rate_input
      real, intent(in) :: dt

      character(len=*), intent(in) :: cumulus

      real, intent(inout) :: entr_rate(:)
      real, intent(inout) :: rh_entr_factor(:)
      real, intent(inout) :: rh_dicycle_fct(:)

      !Local variables:
      integer :: i, k
      real*8  :: y, x
      real :: dpg, trash, dayhr, p_start = 1000.
      real, dimension(its:ite) :: frh, dayhrr

      if (MOIST_TRIGGER /= 2 .and. RH_DICYCLE == 0) return

      !-- ave rh from 1000 -> 450 hPa, following Tian et al 2022 GRL.
      ! OR
      !-- ave rh from 800 -> 450 hPa accounts only for the ´free troposphere'
      if (p_free_troposphere) p_start = 800.
      do i = its, itf
         frh(i) = 0.
         trash = 0.

         do k = kts, ktf
            if (po(i, k) .gt. p_start .and. po(i, k) .lt. 450.) cycle
            dpg = 100.*(po_cup(i, k) - po_cup(i, k + 1))/c_grav
            trash = trash + dpg
            frh(i) = frh(i) + (qo(i, k)/qeso(i, k))*dpg
         end do

         !--average relative humidity
         frh(i) = 100.*frh(i)/(1.e-8 + trash) ! no unit
         frh(i) = max(1., min(100., frh(i)))
         !
         !--- this is for the moist_trigger = 2
         x = dble(frh(i))
         y = 9.192833d0 - 0.2529055d0*x + 0.002344832d0*x**2 - 0.000007230408d0*x**3
         rh_entr_factor(i) = real(y, 4)

         !--- local time
         dayhr = (time_in/3600.+float(itime1_in/100) + float(mod(itime1_in, 100))/60.)
         dayhrr(i) = mod(dayhr + xlons(i)/15.+24., 24.)

         !print*,"FRH=",i,frh(i),rh_dicycle_fct(i),dayhrr,time_in/3600.,xlons(i)
         !print*,"LONS=",i,dayhrr,time_in/3600.,xlons(i)
         !print*,"=================================================="
      end do
      if (MOIST_TRIGGER == 2) then
         entr_rate(:) = entr_rate_input*rh_entr_factor(:)
         !print*,"rh-entr-fac",minval(rh_entr_factor),maxval(rh_entr_factor)
      end if

      if (RH_DICYCLE == 1) then
         do i = its, itf
                        !--- ftun3 controls the domain of the dicycle closure
            !    ftun3 => 1 the closure is insensitive to the mean tropospheric RH
            !    ftun3 => 0 it depends on the RH, with mininum = ftun3
            !rh_dicycle_fct(i) = ftun3 +(1. - ftun3)*&
            !                 (1.-(atan((frh(i)-60.)/10.)+atan(50.))/3.1016)/0.9523154
            if (abs(dayhrr(i) - p_ref_local_time) < 1. .or. time_in < dt + 1.) rh_dicycle_fct(i) = p_ftun3 + (1.-p_ftun3) &
               * (1.-(atan((frh(i) - 55.)/10.) + atan(55.))/3.1016)
            !print*,"fac=",xlons(i),frh(i), rh_dicycle_fct(i);call flush(6)
         end do
      end if
      !-- code to test the atan function
      !  do i = 1 ,100 !relative humidity
      !      y = 0.25 +0.75*(1.-(atan((float(i)-60.)/10.)+atan(50.))/3.1016)/0.9523154
      !      print*,i,y
      !  enddo

      !print*,"FRH",maxval(frh),minval(frh),maxval(rh_dicycle_fct),minval(rh_dicycle_fct)
      !call flush(6)
   end subroutine rhControls

   ! ------------------------------------------------------------------------------------
   function ColdPoolStart(cnv_tr) result(cp_start_out)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! —
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !!
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'ColdPoolStart' ! function name

      real, parameter   :: p_width = 100. 
      !! orig 100
   
      !Variables (input):
      real, intent(in) :: cnv_tr
   
      !Local variables:
      real :: cp_start_out !output

      real :: f1

      f1 = min(MX_BUOY2, cnv_tr)
      !--- f1 > mx_buoy1 => cp_start_out ---> 1
      cp_start_out = (1.35 + atan((f1 - mx_buoy1)/p_width))/2.8
      cp_start_out = max(0.00, min(cp_start_out, 1.00))
      !cp_start_out =  max(0.05,min(cp_start_out,0.95))
   end function ColdPoolStart

   ! ------------------------------------------------------------------------------------
   function FractLiqF(temp2) result(r_fract_liq_f)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! —
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !!
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'FractLiqF' ! function name

      real, parameter :: p_max_temp = 46. 
      !! Celsius
   
      !Variables (input):
      real, intent(in) :: temp2 
      !! K
   
      !Local variables:
      real :: r_fract_liq_f !Output

      real :: temp, ptc

      select case (FRAC_MODIS)
      case (1)
         temp = temp2 - 273.16 !Celsius
         temp = min(p_max_temp, max(-p_max_temp, temp))
         ptc = 7.6725 + 1.0118*temp + 0.1422*temp**2 + 0.0106*temp**3 + 3.39e-4*temp**4 + 3.95e-6*temp**5
         r_fract_liq_f = 1./(1.+exp(-ptc))
         !WMP skew ice fraction for deep convective clouds
         !       fract_liq_f = fract_liq_f**4
         !WMP
      case default
         r_fract_liq_f = min(1., (max(0., (temp2 - c_t_ice))/(c_t00 - c_t_ice))**2)

      end select

   end function FractLiqF


   function initModConvParGF() result(is_init)
      !! # Initialize the module with values
      !!
      !! Author: Rodrigues, L.F. [LFR]
      !!
      !! E-mail: <mailto:luiz.rodrigues@inpe.br>
      !!
      !! Date: 26Janeiro2023 16:41
      !!
      !! #####Version: 0.1.0
      !!
      !! —
      !! **Full description**:
      !!
      !! Initialize all variables from module
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !!
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'initModConvParGF' ! function name
   
      !Local variables:
      integer :: is_init
   
      !Code:
      ! Se o módulo já foi inicializado retorna -1
      if(modConvParGF_initialized) then
         is_init = -1
         return
      endif

      !Inicializa as variáveis do módulo
      ICUMULUS_GF = (/1, 1, 1/) !ok
      CLOSURE_CHOICE = (/0, 7, 3/) 
      CUM_ENTR_RATE = (/ &
                                  9.00e-4 & !deep
                                  , 1.00e-3 & !shallow
                                  , 5.00e-4 & !mid
                                  /)
      CUM_ZUFORM = (/20, 20, 20/) 
      USE_TRACER_TRANSP = 1 
      USE_TRACER_SCAVEN = 2 
      USE_FLUX_FORM = 1 
      USE_FCT = 1 
      USE_TRACER_EVAP = 1 
      CONVECTION_TRACER = 0 
      USE_MEMORY = 2 
      ADD_COLDPOOL_PROP = 3 
      ADD_COLDPOOL_CLOS = 2 
      ADD_COLDPOOL_DIFF = 3 
      USE_SCALE_DEP = 1 
      DICYCLE = 1 
      RH_DICYCLE = 0 
      CLEV_GRID = 1 
      USE_REBCB = 1 
      VERT_DISCR = 1 
      SATUR_CALC = 1 
      SGS_W_TIMESCALE = 1 
      LIGHTNING_DIAG = 0 
      APPLY_SUB_MP = 0 
      ALP1 = 1 
      USE_WETBULB = 0 
      BC_METH = 1 
      OVERSHOOT = 0.
      AUTOCONV = 1     
      C0_DEEP = 1.0e-3
      C0_MID = 1.5e-3
      C0_SHAL = 0.    
      QRC_CRIT = 6.e-4 
      C1 = 0.0   
      USE_MOMENTUM_TRANSP = 1   
      LAMBAU_DEEP = 0.0 
      LAMBAU_SHDN = 2.0 
      DOWNDRAFT = 1   
      TAU_DEEP = 3600.  
      TAU_MID = 1200.  
      MAX_TQ_TEND = 300.   
      USE_SMOOTH_PROF = 1      
      USE_SMOOTH_TEND = 1      
      CUM_HEI_DOWN_LAND = (/0.40, 0.00, 0.35/)
      CUM_HEI_DOWN_OCEAN = (/0.35, 0.00, 0.35/)
      CUM_HEI_UPDF_LAND = (/0.55, 0.10, 0.55/)
      CUM_HEI_UPDF_OCEAN = (/0.55, 0.10, 0.55/)
      CUM_MAX_EDT_LAND = (/0.60, 0.00, 0.20/)
      CUM_MAX_EDT_OCEAN = (/0.60, 0.00, 0.20/)
      CUM_FADJ_MASSFLX = (/1.00, 1.00, 1.00/)
      CUM_AVE_LAYER = (/50., 75., 25./)
      CUM_T_STAR = (/15., -99., -99./)
      CUM_USE_EXCESS = (/1, 1, 1/)
      MOIST_TRIGGER = 0   
      FRAC_MODIS = 1   
      ADV_TRIGGER = 0   
      DCAPE_THRESHOLD = 70. 
      LCL_TRIGGER = 0   
      TAU_OCEA_CP = 7200. 
      TAU_LAND_CP = 7200. 
      MX_BUOY1 = (real(c_cp)*5.0 + real(c_alvl)*2.e-3)*0.025  
      MX_BUOY2 = (real(c_cp)*10.+real(c_alvl)*4.e-3) 
      USE_CLOUD_DISSIPATION = 0.   
      USE_GUSTINESS = 0    
      USE_RANDOM_NUM = 0.   
      BETA_SH = 2.2  
      USE_LINEAR_SUBCL_MF = 1    
      CAP_MAXS = 50.  
      LIQ_ICE_NUMBER_CONC = 1    
      ALPHA_ADV_TUNING = 0.8  
      SIG_FACTOR = 0.22 != exponential factor for the sigma determination (orig = 0.1)

      hei_down_land  = 0.     
      hei_down_ocean = 0.    
      hei_updf_land  = 0.    
      hei_updf_ocean = 0.    
      max_edt_land   = 0.    
      max_edt_ocean  = 0.    
      fadj_massflx  = 0.    
      t_star      = 0.    
      ave_layer      = 0.    
      c0             = 0.    
      col_sat_adv_threshold = 0.94 
      chem_adj_autoc = 0.
      time_in = 0.
      int_time = 0.
      ispc_co = 0
      whoami_all = 0
      jcol = 0
      itime1_in = 0
      nrec = 0
      ntimes = 0
       use_excess   = 0      
      output_sound = 0   
      ind_chem = 0
      chem_name_mask = 0
      chem_name_mask_evap = 0
      use_c1d = .false. 
      first_guess_w = .false. 
      wrtgrads = .false.
      chem_name = ""
      if(allocated(hcts)) then
         hcts%hstar = 0.
         hcts%dhr = 0.
         hcts%ak0 = 0.
         hcts%dak = 0.
      endif
      ! Informa que já inicializado
      modConvParGF_initialized = .true.
      ! Retorna 0, foi inicializado dessa vez
      is_init = 0

   end function initModConvParGF

   subroutine initDriverJ(its,itf, aot500, xland, sfc_press, temp2m, topt, kpbl, lons, lats, rtgt, ztexec &
                     ,   zqexec, last_ierr, fixout_qv, conprr, lightn_dens, var2d, revsu_gf_2d, prfil_gf_2d &
                     ,   var3d_agf_2d, var3d_bgf_2d, t_pert_2d, temp_tendqv, outt, outu, outv, outq, outqc, outnice &
                     ,   outnliq, outbuoy, outmpqi, outmpql, outmpcf, omeg, out_chem, ccn, xlandi, psur, tsur &
                     ,   ter11, kpbli, xlons, xlats)
      !! ## Init some variables and arrays for driver
      !!
      !! Author: Rodrigues, L.F [LFR]
      !!
      !! E-mail: <mailto:luizfrodrigues@protonmail.com>
      !!
      !! Date: 10Fevereiro2023 08:50
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Init some variables and arrays for driver
      !!
      !! ** History**:
      !!
      !! --- 
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!

      implicit none
      ! Parameters:
      character(len=*), parameter :: p_procedure_name = 'initDriver' 
      !! subroutine name

      ! Variables (input, output, inout)
      integer, intent(in) :: its, itf

      integer, intent(in) :: kpbl(:)

      real, intent(in) :: aot500(:)
      real, intent(in) :: xland(:)
      real, intent(in) :: sfc_press(:)
      real, intent(in) :: temp2m(:)
      real, intent(in) :: topt(:)
      real, intent(in) :: lons(:)
      real, intent(in) :: lats(:)

      real, intent(inout) :: outmpqi(:,:,:,:)
      real, intent(inout) :: outmpql(:,:,:,:)
      real, intent(inout) :: outmpcf(:,:,:,:)
      real, intent(inout) :: out_chem(:,:,:,:)
      real, intent(inout) :: omeg(:,:,:)
      real, intent(inout) :: ccn(:)

      integer, intent(out) :: kpbli (:)
      integer, intent(out) :: last_ierr(:)

      real, intent(out) :: rtgt(:)
      real, intent(out) :: ztexec(:)
      real, intent(out) :: zqexec(:)
      real, intent(out) :: fixout_qv(:)
      real, intent(out) :: conprr(:)
      real, intent(out) :: lightn_dens(:)
      real, intent(out) :: var2d(:)
      real, intent(out) :: revsu_gf_2d(:,:)
      real, intent(out) :: prfil_gf_2d(:,:)
      real, intent(out) :: var3d_agf_2d(:,:)
      real, intent(out) :: var3d_bgf_2d(:,:)
      real, intent(out) :: t_pert_2d(:,:)
      real, intent(out) :: temp_tendqv(:, :)
      real, intent(out) :: outt   (:,:,:)
      real, intent(out) :: outu(:,:,:)
      real, intent(out) :: outv(:,:,:)
      real, intent(out) :: outq(:,:,:)
      real, intent(out) :: outqc(:,:,:)
      real, intent(out) :: outnice(:,:,:)
      real, intent(out) :: outnliq(:,:,:)
      real, intent(out) :: outbuoy(:,:,:)
      real, intent(out) :: xlandi(:)
      real, intent(out) :: psur  (:)
      real, intent(out) :: tsur  (:)
      real, intent(out) :: ter11 (:)
      real, intent(out) :: xlons (:)
      real, intent(out) :: xlats (:)
      
      ! Local variables:
      integer :: i_cnt

      !-- initialization
      do i_cnt = its, itf
         rtgt     (i_cnt) = 1.0
         ztexec   (i_cnt) = 0.0
         zqexec   (i_cnt) = 0.0
         last_ierr(i_cnt) = -999
         fixout_qv(i_cnt) = 1.0
         !
         conprr     (i_cnt) = 0.0
         lightn_dens(i_cnt) = 0.0
         var2d      (i_cnt) = 0.0
         !--- (i,k)
         revsu_gf_2d (i_cnt, :) = 0.0
         prfil_gf_2d (i_cnt, :) = 0.0
         var3d_agf_2d(i_cnt, :) = 0.0
         var3d_bgf_2d(i_cnt, :) = 0.0
         t_pert_2d    (i_cnt, :) = 0.0
         !
         temp_tendqv(i_cnt, :) = 0.0
         !- tendencies (w/ maxiens)
         outt   (i_cnt, :, :) = 0.0
         outu   (i_cnt, :, :) = 0.0
         outv   (i_cnt, :, :) = 0.0
         outq   (i_cnt, :, :) = 0.0
         outqc  (i_cnt, :, :) = 0.0
         outnice(i_cnt, :, :) = 0.0
         outnliq(i_cnt, :, :) = 0.0
         outbuoy(i_cnt, :, :) = 0.0
      end do
      if (APPLY_SUB_MP == 1) then
         do i_cnt = its, itf
            !- tendencies (w/ nmp and maxiens)
            outmpqi(:, i_cnt, :, :) = 0.0
            outmpql(:, i_cnt, :, :) = 0.0
            outmpcf(:, i_cnt, :, :) = 0.0
         end do
      end if
      do i_cnt = its, itf
         omeg(i_cnt, :, :) = 0.0
      end do
      if (USE_TRACER_TRANSP == 1) then
         out_chem = 0.0
      end if
      !
      if (AUTOCONV == 2) then
         do i_cnt = its, itf
            ccn(i_cnt) = max(100., (370.37*(0.01 + max(0., aot500(i_cnt))))**1.555)
         end do
      else
         do i_cnt = its, itf
            ccn(i_cnt) = 100.
         end do
      end if
      do i_cnt = its, itf
         xlandi(i_cnt) = xland(i_cnt)!flag < 1 para land
         !flag  =1 para water
         psur (i_cnt) = sfc_press(i_cnt)*1.e-2 ! mbar
         tsur (i_cnt) = temp2m(i_cnt)
         ter11(i_cnt) = max(0., topt(i_cnt))
         kpbli(i_cnt) = kpbl(i_cnt)
         xlons(i_cnt) = lons(i_cnt)*180./3.14159
         xlats(i_cnt) = lats(i_cnt)*180./3.14159
      end do

   end subroutine initDriverJ

   subroutine initDriverK(its, itf, kts, ktf, kte, mtp, dt, rtgt, topt, zt, press, temp, rvap, curr_rvap, u_wind, v_wind &   
                        , w_wind, dm, buoy_exc, rth_advten, rqvften, flip, mp_ice, mp_liq, mp_cf, tracer, zo, po, temp_old  &  
                        , qv_old, qv_curr, rhoi, tkeg, rcpg, us, vs, dm2d, buoy_exc2d, temp_new_adv, qv_new_adv &
                        , mpqi, mpql, mpcf, se_chem, omeg)
      !! ## Init some other variables from driver
      !!
      !! Author: Rodrigues, L.F. [LFR]
      !!
      !! E-mail: <mailto:luizfrodrigues@protonmail.com>
      !!
      !! Date: 10Fevereiro2023 10:15
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Init some other variables from driver
      !!
      !! ** History**:
      !!
      !! --- 
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!
   
      implicit none
      !Parameters:
      character(len=*), parameter :: p_procedure_name = 'initDriverK' 
      !! subroutine name
   
      !Variables (input, output, inout)
      integer, intent(in) :: its, itf, kts, ktf, kte, mtp

      integer, intent(in) :: flip(:)

      real, intent(in) :: dt
      real, intent(in) :: rtgt      (:)
      real, intent(in) :: topt      (:)
      real, intent(in) :: zt        (:,:)
      real, intent(in) :: press     (:,:)
      real, intent(in) :: temp      (:,:)
      real, intent(in) :: rvap      (:,:)
      real, intent(in) :: curr_rvap (:,:)
      real, intent(in) :: u_wind    (:,:)
      real, intent(in) :: v_wind    (:,:)
      real, intent(in) :: w_wind    (:,:)
      real, intent(in) :: dm        (:,:)
      real, intent(in) :: buoy_exc  (:,:)
      real, intent(in) :: rth_advten(:,:)
      real, intent(in) :: rqvften   (:,:)
      real, intent(in) :: mp_ice    (:,:,:)
      real, intent(in) :: mp_liq    (:,:,:)
      real, intent(in) :: mp_cf     (:,:,:)
      real, intent(in) :: tracer    (:,:,:)

      real, intent(inout) :: mpqi   (:,:,:)
      real, intent(inout) :: mpql   (:,:,:)
      real, intent(inout) :: mpcf   (:,:,:)
      real, intent(inout) :: se_chem(:,:,:)

      real, intent(out) :: zo          (:,:)
      real, intent(out) :: po          (:,:)
      real, intent(out) :: temp_old    (:,:)
      real, intent(out) :: qv_old      (:,:)
      real, intent(out) :: qv_curr     (:,:)
      real, intent(out) :: rhoi        (:,:)
      real, intent(out) :: tkeg        (:,:)
      real, intent(out) :: rcpg        (:,:)
      real, intent(out) :: us          (:,:)
      real, intent(out) :: vs          (:,:)
      real, intent(out) :: dm2d        (:,:)
      real, intent(out) :: buoy_exc2d  (:,:)
      real, intent(out) :: temp_new_adv(:,:)
      real, intent(out) :: qv_new_adv  (:,:)
      real, intent(out) :: omeg(:,:,:)

      !Local variables:
      integer :: k_cnt, i_cnt, kr, ispc
   
      do k_cnt = kts, ktf
         do i_cnt = its, itf
            kr = k_cnt   !+1   !<<<< only kr=k (the input was already converted to the BRAMS vertical grid,
            !                see cup_grell3.f90 routine)
            !- heigths, current pressure, temp and water vapor mix ratio
            zo      (i_cnt, k_cnt) = zt(kr, i_cnt)*rtgt(i_cnt) + topt(i_cnt)
            po      (i_cnt, k_cnt) = press(kr, i_cnt)*1.e-2 !mbar
            temp_old(i_cnt, k_cnt) = temp(kr, i_cnt)
            qv_old  (i_cnt, k_cnt) = rvap(kr, i_cnt) ! @ begin of the timestep
            qv_curr (i_cnt, k_cnt) = curr_rvap(kr, i_cnt) ! current (after dynamics + physical processes called before GF)
            !- air density, TKE and cloud liq water mixing ratio
            rhoi    (i_cnt, k_cnt) = 1.e2*po(i_cnt, k_cnt)/(287.04*temp_old(i_cnt, k_cnt)*(1.+0.608*qv_old(i_cnt, k_cnt)))
            tkeg    (i_cnt, k_cnt) = c_tkmin_ms
            rcpg    (i_cnt, k_cnt) = 0.
            !- wind velocities
            us      (i_cnt, k_cnt) = u_wind(kr, i_cnt)
            vs      (i_cnt, k_cnt) = v_wind(kr, i_cnt)
            omeg(i_cnt, k_cnt, :) = w_wind(kr, i_cnt)
            dm2d(i_cnt, k_cnt) = dm(kr, i_cnt)
            !omeg   (i,k,:)= -g*rho(kr,i,j)*w(kr,i,j)
            !-buoyancy excess
            buoy_exc2d(i_cnt, k_cnt) = buoy_exc(kr, i_cnt)
            !- temp/water vapor modified only by advection
            temp_new_adv(i_cnt, k_cnt) = temp_old(i_cnt, k_cnt) + (rth_advten(kr, i_cnt))*dt
            qv_new_adv(i_cnt, k_cnt) = qv_old(i_cnt, k_cnt) + (rqvften(kr, i_cnt))*dt
         end do
      end do
      if (APPLY_SUB_MP == 1) then
         do k_cnt = kts, ktf
            do i_cnt = its, itf
               kr = k_cnt   !+1   !<<<< only kr=k
               !- microphysics ice and liq mixing ratio, and cloud fraction of the host model
               !- (only subsidence is applied)
               mpqi(:, i_cnt, k_cnt) = mp_ice(:, kr, i_cnt) ! kg/kg
               mpql(:, i_cnt, k_cnt) = mp_liq(:, kr, i_cnt) ! kg/kg
               mpcf(:, i_cnt, k_cnt) = mp_cf(:, kr, i_cnt) ! 1
            end do
         end do
      end if
      if (USE_TRACER_TRANSP == 1) then
         do k_cnt = kts, kte
            do i_cnt = its, itf
               kr = k_cnt !+1
               !- atmos composition
               do ispc = 1, mtp
                  se_chem(ispc, i_cnt, k_cnt) = max(p_mintracer, tracer(i_cnt, flip(kr), ispc))
               end do
            end do
         end do
      end if

   
   end subroutine initDriverK

end module modConvParGF

