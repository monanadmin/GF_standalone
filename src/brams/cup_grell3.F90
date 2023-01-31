!--------------------------------------------------------------------------------!
! Cumulus Parameterizations G3d and GF                                          !
! Implemented in BRAMS by Saulo Freitas @ Feb/2012 /  Feb/2021                               !
! Rafael Mello: included parallelization for spread/g3d_smoothh arrays               !
!--------------------------------------------------------------------------------!
MODULE CUPARM_GRELL3

  use ModNamelistFile, only: namelistFile

  use ModMessageSet, only: &
      PostRecvSendMsgs,    &
      WaitRecvMsgs

  use mem_basic         , only: basic_g
  use mem_tend          , only: tend
  use mem_cuparm        , only: confrq,cuparm_g,cuparm_g_sh&
                               ,nnqparm

  use node_mod          , only: mynum,   &   ! INTENT(IN)
             		        mxp,     &   ! INTENT(IN)
             		        myp,     &   ! INTENT(IN)
             		        mzp,     &   ! INTENT(IN)
             		        ia,      &   ! INTENT(IN)
             		        iz,      &   ! INTENT(IN)
             		        ja,      &   ! INTENT(IN)
             		        jz,      &   ! INTENT(IN)
             		        i0,      &   ! INTENT(IN)
             		        j0,      &   ! INTENT(IN)
			                 ibcon        ! INTENT(IN)

  use mem_grid          , only: time,    &   ! INTENT(IN)
            		   	initial, &   ! INTENT(IN)
            		   	dtlt,	 &   ! INTENT(IN)
            		   	itime1,  &   ! INTENT(IN)
            		   	ngrid,   &   ! INTENT(IN)
            		   	grid_g,  &   ! INTENT(IN)
            		   	dtlongn, &   ! INTENT(IN)
           		         deltaxn, &   ! INTENT(IN)
           		   	   deltayn, &   ! INTENT(IN)
				            npatch,  &   ! INTENT(IN)
				            ztn,     &   ! INTENT(IN)
				            zmn,     &   ! INTENT(IN)
				            akminvar,&   ! INTENT(IN)
                        nxtnest, &
                        naddsc

  use mem_varinit, only: nudlat

  use rconstants        , only: tkmin
  use mem_turb          , only: turb_g,akmin
  use mem_micro         , only: micro_g
  use io_params         , only: frqanl
  use mem_leaf          , only: leaf_g, isfcl
  use micphys           , only: level,mcphys_type

  use mem_grell_param, only: mgmxp,mgmyp,mgmzp,maxiens,ngrids_cp,&
       maxens_g3d ,                        & !INTENT(IN)
       maxens2_g3d,                        & !INTENT(IN)
       maxens3_g3d,                        & !INTENT(IN)
       ensdim_g3d ,                        & !INTENT(IN)
       maxens ,                            & !INTENT(IN)
       maxens2,                            & !INTENT(IN)
       maxens3,                            & !INTENT(IN)
       ensdim ,                            & !INTENT(IN)
       icoic

  use mem_scratch1_grell, only: ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d,kpbl4d,   &
                                kstabi4d,kstabm4d,xmb4d,edt4d,pwav4d,               &
                                zup5d, zdn5d,iruncon, pcup5d, prup5d,prdn5d,        &
				                    clwup5d,tup5d,enup5d,endn5d,deup5d,dedn5d,zcup5d,   &
                                up_massdetr5d, up_massentr5d,                       &
                                dd_massdetr5d, dd_massentr5d,                       &
				                    conv_cld_fr5d,sigma4d,klcl4d,cprr4d


  use mem_grell         , only: cuforc_g,cuforc_sh_g

  use mem_carma, only: carma

  use mem_radiate, only: ISWRTYP, ILWRTYP,radiate_g ! Intent(in)

  use mem_turb, only:  idiffk !INTENT(IN)

!-----------Grell G3d - GD-FIM - GF
  use module_cu_g3    , only:  G3DRV
  use module_cu_gf    , only:  GFDRV
  use module_cu_gf2   , only:  GFDRV2

  USE Phys_const, only: cp, p00, tcrit, g, cpor , XL, rm,rgas

!----------- GF - GEOS-5
  USE ConvPar_GF_GEOS5, only: GF_GEOS5_DRV, deep, shal, mid , nmp, lsmp , cnmp &
                            , GF_convpar_init,apply_sub_mp,icumulus_gf, liq_ice_number_conc &
                            , convection_tracer, tau_ocea_cp, tau_land_cp,use_gustiness

! use mem_scalar, only: scalar_g

!-----------

  use ccatt_start, only: ccatt
  use mem_chem1  , only: chemistry

  use mem_jules, only: jules_g


  use ModPostUtils, only : rs

  implicit none

  TYPE g3d_ens_vars
     REAL, POINTER, DIMENSION(:,:)  ::apr
     REAL, POINTER, DIMENSION(:,:)  ::accapr
     REAL, POINTER, DIMENSION(:,:)  ::weight
     !-----------
  END TYPE g3d_ens_vars
  TYPE (g3d_ens_vars)    , allocatable :: g3d_ens_g(:,:), g3d_ensm_g(:,:)

  TYPE g3d_vars
     !-- 2 dimensions
     REAL, POINTER, DIMENSION(:,:  )  ::xmb_deep
     REAL, POINTER, DIMENSION(:,:  )  ::xmb_deep_dd
     REAL, POINTER, DIMENSION(:,:  )  ::err_deep
     REAL, POINTER, DIMENSION(:,:  )  ::xmb_shallow
     REAL, POINTER, DIMENSION(:,:  )  ::rh_dicy_fct
     REAL, POINTER, DIMENSION(:,:  )  ::lightn_dens
     REAL, POINTER, DIMENSION(:,:  )  ::aa0
     REAL, POINTER, DIMENSION(:,:  )  ::aa1   
     REAL, POINTER, DIMENSION(:,:  )  ::aa1_bl
     REAL, POINTER, DIMENSION(:,:  )  ::aa1_adv
     REAL, POINTER, DIMENSION(:,:  )  ::aa1_radpbl
     REAL, POINTER, DIMENSION(:,:  )  ::vcpool
     REAL, POINTER, DIMENSION(:,:  )  ::umcl
     REAL, POINTER, DIMENSION(:,:  )  ::vmcl
     REAL, POINTER, DIMENSION(:,:  )  ::wlpool
     REAL, POINTER, DIMENSION(:,:  )  ::sgs_shf
     REAL, POINTER, DIMENSION(:,:  )  ::sgs_lhf

     !-- 3 dimensions     
     REAL, POINTER, DIMENSION(:,:,:)  ::cugd_ttens
     REAL, POINTER, DIMENSION(:,:,:)  ::cugd_qvtens
     REAL, POINTER, DIMENSION(:,:,:)  ::thsrc
     REAL, POINTER, DIMENSION(:,:,:)  ::rtsrc
     REAL, POINTER, DIMENSION(:,:,:)  ::clsrc
     REAL, POINTER, DIMENSION(:,:,:)  ::nlsrc
     REAL, POINTER, DIMENSION(:,:,:)  ::nisrc
     REAL, POINTER, DIMENSION(:,:,:)  ::usrc
     REAL, POINTER, DIMENSION(:,:,:)  ::vsrc
     REAL, POINTER, DIMENSION(:,:,:)  ::buoysrc
     REAL, POINTER, DIMENSION(:,:,:)  ::cnv_tr
     REAL, POINTER, DIMENSION(:,:,:)  ::qexcp
     REAL, POINTER, DIMENSION(:,:,:)  ::hexcp


     REAL, POINTER, DIMENSION(:,:,:)  ::mupdp
     REAL, POINTER, DIMENSION(:,:,:)  ::mdddp
     REAL, POINTER, DIMENSION(:,:,:)  ::mupsh
     REAL, POINTER, DIMENSION(:,:,:)  ::mupmd
     REAL, POINTER, DIMENSION(:,:,:)  ::eupdp
     REAL, POINTER, DIMENSION(:,:,:)  ::dupdp
  END TYPE g3d_vars

  TYPE (g3d_vars)       , allocatable :: g3d_g(:),g3dm_g(:)

  integer ::    ids,ide, jds,jde, kds,kde            &
               ,ims,ime, jms,jme, kms,kme            &
               ,ips,ipe, jps,jpe, kps,kpe            &
               ,its,ite, jts,jte, kts,kte

  integer           ::  imomentum=1
  integer           ::  ishallow_g3

  !- define if the training will be used or not
  INTEGER,PARAMETER:: training=0
  character(len=255) :: g3d_training_file

  !- define if the lateral subsidence spread will be done or not
  INTEGER :: g3d_spread ! 1=ON, 0=OFF
  INTEGER :: cugd_avedx

  !- define if the horizontal smoothing is to be done or not
  INTEGER :: g3d_smoothh! 1=ON, 0=OFF

  !- define if the vertical smoothing is to be done or not
  INTEGER :: g3d_smoothv! 1=ON, 0=OFF

  !- number of members of prec ensemble
  INTEGER,PARAMETER :: train_dim= 5

  CHARACTER(LEN=6),PARAMETER,DIMENSION(train_dim) :: pre_name=(/ &
      'apr_gr' & !
     ,'apr_w ' & !
     ,'apr_mc' & !
     ,'apr_st' & !
     ,'apr_as' & !
   /)

  INTEGER,PARAMETER :: apr_gr=001
  INTEGER,PARAMETER :: apr_w =002
  INTEGER,PARAMETER :: apr_mc=003
  INTEGER,PARAMETER :: apr_st=004
  INTEGER,PARAMETER :: apr_as=005


  integer,parameter :: CPTIME = 0. !orig: CPTIME = 7200.

  integer,parameter :: i_forcing = 1

  integer,parameter :: autoconv = 1!2 ! =1, Kessler
                                      ! =2, Berry
  integer,parameter :: aerovap = 1!3  ! =1, orig
                                    ! =2, mix orig+new
				    ! =3, new
  !- direct link cupar-microphysics
  LOGICAL,parameter :: do_cupar_mcphys_coupling = .true.
  !- read GF namelist
  LOGICAL :: read_GF_ConvPar_nml =.true. 

  character(len=*), parameter :: header="**(cup_grell3)**"

Contains
!-----------------------------------------
  subroutine nullify_grell3(g3d_ens,g3d,ndim_train)

    implicit none
    integer, intent(in) ::ndim_train
    type (g3d_ens_vars),dimension(ndim_train) :: g3d_ens
    type (g3d_vars) :: g3d
    integer i

    do i=1,ndim_train
      if (associated(g3d_ens(i)%apr))    nullify (g3d_ens(i)%apr)
      if (associated(g3d_ens(i)%accapr)) nullify (g3d_ens(i)%accapr)
      if (associated(g3d_ens(i)%weight)) nullify (g3d_ens(i)%weight)
    enddo

    if (associated(g3d%xmb_deep))    nullify (g3d%xmb_deep)
    if (associated(g3d%xmb_deep_dd)) nullify (g3d%xmb_deep_dd)
    if (associated(g3d%err_deep))    nullify (g3d%err_deep)
    if (associated(g3d%xmb_shallow)) nullify (g3d%xmb_shallow)
    if (associated(g3d%rh_dicy_fct)) nullify (g3d%rh_dicy_fct)
    if (associated(g3d%lightn_dens)) nullify (g3d%lightn_dens)
    if (associated(g3d%aa0))    nullify (g3d%aa0)
    if (associated(g3d%aa1))    nullify (g3d%aa1)
    if (associated(g3d%aa1_bl)) nullify (g3d%aa1_bl)
    if (associated(g3d%aa1_adv)) nullify (g3d%aa1_adv)
    if (associated(g3d%aa1_radpbl)) nullify (g3d%aa1_radpbl)
    if (associated(g3d%vcpool))    nullify (g3d%vcpool)
    if (associated(g3d%umcl  ))    nullify (g3d%umcl)
    if (associated(g3d%vmcl  ))    nullify (g3d%vmcl)
    if (associated(g3d%wlpool))    nullify (g3d%wlpool)
    if (associated(g3d%sgs_shf))   nullify (g3d%sgs_shf)
    if (associated(g3d%sgs_lhf))   nullify (g3d%sgs_lhf)
    if (associated(g3d%sgs_lhf))   nullify (g3d%sgs_lhf)
    if (associated(g3d%sgs_lhf))   nullify (g3d%sgs_lhf)
    if (associated(g3d%qexcp  ))   nullify (g3d%qexcp  )
    if (associated(g3d%hexcp  ))   nullify (g3d%hexcp  )


    if (associated(g3d%cugd_ttens))  nullify (g3d%cugd_ttens)
    if (associated(g3d%cugd_qvtens)) nullify (g3d%cugd_qvtens)

    if (associated(g3d%thsrc)) nullify (g3d%thsrc)
    if (associated(g3d%rtsrc)) nullify (g3d%rtsrc)
    if (associated(g3d%clsrc)) nullify (g3d%clsrc)
    if (associated(g3d%nlsrc)) nullify (g3d%nlsrc)
    if (associated(g3d%nisrc)) nullify (g3d%nisrc)
    if (associated(g3d%usrc )) nullify (g3d%usrc)
    if (associated(g3d%vsrc )) nullify (g3d%vsrc)
    if (associated(g3d%buoysrc )) nullify (g3d%buoysrc)
    if (associated(g3d%cnv_tr )) nullify (g3d%cnv_tr)

    if (associated(g3d%mupdp)) nullify (g3d%mupdp)
    if (associated(g3d%mupsh)) nullify (g3d%mupsh)
    if (associated(g3d%mdddp)) nullify (g3d%mdddp)
    if (associated(g3d%mupmd)) nullify (g3d%mupmd)
    
    if (associated(g3d%eupdp)) nullify (g3d%eupdp)
    if (associated(g3d%dupdp)) nullify (g3d%dupdp)

  end subroutine nullify_grell3
!-----------------------------------------
   subroutine alloc_grell3(g3d_ens,g3d, m1, m2, m3, ng,ndim_train)
    implicit none
    type (g3d_ens_vars),dimension(ndim_train) :: g3d_ens
    type (g3d_vars) :: g3d
    integer, intent(in) :: m1, m2, m3, ng,ndim_train
    integer :: i

    if(nnqparm(ng) == 3 .or. nnqparm(ng) == 6 .or. nnqparm(ng) == 5  &
       .or. nnqparm(ng) == 8) then
      do i=1,ndim_train
       allocate(g3d_ens(i)%apr   (m2,m3))   ; g3d_ens(i)%apr    =0.0
       allocate(g3d_ens(i)%accapr(m2,m3))   ; g3d_ens(i)%accapr =0.0
       allocate(g3d_ens(i)%weight(m2,m3))   ; g3d_ens(i)%weight =0.0
      enddo
      allocate (g3d%cugd_ttens (m1, m2, m3))  ;g3d%cugd_ttens =0.0
      allocate (g3d%cugd_qvtens(m1, m2, m3))  ;g3d%cugd_qvtens=0.0
    endif

    allocate (g3d%xmb_deep   (m2,m3))       ;g3d%xmb_deep   =0.0
    allocate (g3d%xmb_deep_dd(m2,m3))       ;g3d%xmb_deep_dd=0.0
    allocate (g3d%err_deep   (m2,m3))       ;g3d%err_deep   =0.0
    allocate (g3d%xmb_shallow(m2,m3))       ;g3d%xmb_shallow=0.0
    allocate (g3d%rh_dicy_fct(m2,m3))       ;g3d%rh_dicy_fct=0.0
    allocate (g3d%lightn_dens(m2,m3))       ;g3d%lightn_dens=0.0
    allocate (g3d%aa0        (m2,m3))       ;g3d%aa0        =0.0
    allocate (g3d%aa1        (m2,m3))       ;g3d%aa1        =0.0
    allocate (g3d%aa1_bl     (m2,m3))       ;g3d%aa1_bl     =0.0
    allocate (g3d%aa1_adv    (m2,m3))       ;g3d%aa1_adv    =0.0
    allocate (g3d%aa1_radpbl (m2,m3))       ;g3d%aa1_radpbl =0.0
    allocate (g3d%vcpool     (m2,m3))       ;g3d%vcpool     =0.0
    allocate (g3d%vmcl       (m2,m3))       ;g3d%vmcl       =0.0
    allocate (g3d%umcl       (m2,m3))       ;g3d%umcl       =0.0
    allocate (g3d%wlpool     (m2,m3))       ;g3d%wlpool     =0.0
    allocate (g3d%sgs_shf    (m2,m3))       ;g3d%sgs_shf    =0.0
    allocate (g3d%sgs_lhf    (m2,m3))       ;g3d%sgs_lhf    =0.0
    allocate (g3d%qexcp      (m1,m2,m3))    ;g3d%qexcp      =0.0
    allocate (g3d%hexcp      (m1,m2,m3))    ;g3d%hexcp      =0.0

    allocate (g3d%thsrc(m1, m2, m3))  ;g3d%thsrc=0.0
    allocate (g3d%rtsrc(m1, m2, m3))  ;g3d%rtsrc=0.0
    allocate (g3d%clsrc(m1, m2, m3))  ;g3d%clsrc=0.0
    allocate (g3d%nlsrc(m1, m2, m3))  ;g3d%nlsrc=0.0
    allocate (g3d%nisrc(m1, m2, m3))  ;g3d%nisrc=0.0

    if( (imomentum==0 .or. imomentum==1 ) .and. nnqparm(ng) >= 4) then
     allocate (g3d%usrc(m1, m2, m3))  ;g3d%usrc=0.0
     allocate (g3d%vsrc(m1, m2, m3))  ;g3d%vsrc=0.0
    endif
    allocate (g3d%mupsh(m1, m2, m3))  ;g3d%mupsh=0.0
    allocate (g3d%mupdp(m1, m2, m3))  ;g3d%mupdp=0.0
    allocate (g3d%mdddp(m1, m2, m3))  ;g3d%mdddp=0.0
    allocate (g3d%mupmd(m1, m2, m3))  ;g3d%mupmd=0.0

    allocate (g3d%eupdp(m1, m2, m3))  ;g3d%eupdp=0.0
    allocate (g3d%dupdp(m1, m2, m3))  ;g3d%dupdp=0.0
  
    if(nnqparm(ng) == 8 ) then ! and convection_tracer = 1
      allocate (g3d%buoysrc(m1, m2, m3))  ;g3d%buoysrc=0.0
      allocate (g3d%cnv_tr (m1, m2, m3))  ;g3d%cnv_tr =0.0
    endif 
  end subroutine alloc_grell3

!-----------------------------------------
  subroutine filltab_grell3(g3d_ens,g3d,g3d_ensm,g3dm,imean, m1, m2, m3, ng,ndim_train)
    use var_tables
    implicit none
    include "i8.h"

    type (g3d_ens_vars),dimension(ndim_train) :: g3d_ens,g3d_ensm
    type (g3d_vars) :: g3d,g3dm
    integer, intent(in) :: imean, m1, m2, m3, ng,ndim_train
    integer(kind=i8) :: npts
    integer :: i
    character (len=4) :: arrprop
    ! Fill pointers to arrays into variable tables

    npts=m2*m3
    do i=1,ndim_train
     if (associated(g3d_ens(i)%apr))  &
         call InsertVTab (g3d_ens(i)%apr   ,g3d_ensm(i)%apr    &
         ,ng, npts, imean,  &
         trim(pre_name(i))//' :2:hist:mpti:mpt3')

     if (associated(g3d_ens(i)%accapr))  &
         call InsertVTab (g3d_ens(i)%accapr   ,g3d_ensm(i)%accapr    &
         ,ng, npts, imean,  &
         'acc'//trim(pre_name(i)(2:len_trim(pre_name(i))))//' :2:hist:anal:mpti:mpt3')
    enddo

    do i=1,ndim_train
     if (associated(g3d_ens(i)%weight))  &
         call InsertVTab (g3d_ens(i)%weight   ,g3d_ensm(i)%weight    &
         ,ng, npts, imean,  &
         'weight'//trim(pre_name(i)(4:len_trim(pre_name(i))))//' :2:hist:anal:mpti:mpt3')

    enddo

    if (associated(g3d%xmb_deep))  &
         call InsertVTab (g3d%xmb_deep   ,g3dm%xmb_deep    &
         ,ng, npts, imean,  &
         'MFUP :2:hist:anal:mpti:mpt3')
    if (associated(g3d%xmb_deep_dd))  &
         call InsertVTab (g3d%xmb_deep_dd   ,g3dm%xmb_deep_dd    &
         ,ng, npts, imean,  &
         'MFDD :2:hist:anal:mpti:mpt3')
    if (associated(g3d%err_deep))  &
         call InsertVTab (g3d%err_deep   ,g3dm%err_deep    &
         ,ng, npts, imean,  &
         'XIERR :2:hist:anal:mpti:mpt3')
    if (associated(g3d%xmb_shallow))  &
         call InsertVTab (g3d%xmb_shallow   ,g3dm%xmb_shallow    &
         ,ng, npts, imean,  &
         'MFSH :2:hist:anal:mpti:mpt3')
    if (associated(g3d%rh_dicy_fct))  &
         call InsertVTab (g3d%rh_dicy_fct   ,g3dm%rh_dicy_fct    &
         ,ng, npts, imean,  &
         'RHDICY :2:hist:anal:mpti:mpt3')

    if (associated(g3d%lightn_dens))  &
         call InsertVTab (g3d%lightn_dens   ,g3dm%lightn_dens    &
         ,ng, npts, imean,  &
         'LIGHTN :2:hist:anal:mpti:mpt3')
 
    if (associated(g3d%aa0)) call InsertVTab (g3d%aa0   ,g3dm%aa0    &
         ,ng, npts, imean, 'AA0 :2:hist:anal:mpti:mpt3')
    if (associated(g3d%aa1)) call InsertVTab (g3d%aa1   ,g3dm%aa1    &
         ,ng, npts, imean, 'AA1 :2:hist:anal:mpti:mpt3')
    if (associated(g3d%aa1_bl)) call InsertVTab (g3d%aa1_bl   ,g3dm%aa1_bl    &
         ,ng, npts, imean, 'AA1_BL :2:hist:anal:mpti:mpt3')
    if (associated(g3d%aa1_radpbl)) call InsertVTab (g3d%aa1_radpbl   ,g3dm%aa1_radpbl    &
         ,ng, npts, imean, 'AA1_RADPBL :2:hist:anal:mpti:mpt3')
    if (associated(g3d%aa1_adv)) call InsertVTab (g3d%aa1_adv   ,g3dm%aa1_adv    &
         ,ng, npts, imean, 'AA1_ADV :2:hist:anal:mpti:mpt3')
    
    if (associated(g3d%vcpool)) call InsertVTab (g3d%vcpool   ,g3dm%vcpool    &
        ,ng, npts, imean, 'VCPOOL :2:hist:anal:mpti:mpt3')
    if (associated(g3d%vmcl)) call InsertVTab (g3d%vmcl   ,g3dm%vmcl    &
        ,ng, npts, imean, 'VMCL :2:hist:anal:mpti:mpt3')
    if (associated(g3d%umcl)) call InsertVTab (g3d%umcl   ,g3dm%umcl    &
        ,ng, npts, imean, 'UMCL :2:hist:anal:mpti:mpt3')
    if (associated(g3d%wlpool)) call InsertVTab (g3d%wlpool   ,g3dm%wlpool    &
        ,ng, npts, imean, 'WLPOOL :2:hist:anal:mpti:mpt3')

    if (associated(g3d%sgs_shf)) call InsertVTab (g3d%sgs_shf   ,g3dm%sgs_shf   &
        ,ng, npts, imean, 'SGS_SHF :2:hist:anal:mpti:mpt3')
    if (associated(g3d%sgs_lhf)) call InsertVTab (g3d%sgs_lhf   ,g3dm%sgs_lhf   &
        ,ng, npts, imean, 'SGS_LHF :2:hist:anal:mpti:mpt3')

 

    !- 3D Arrays
    npts=m1*m2*m3

    !- define if the arrays will exchange 1 row x 1 line (not in use anymore)
    arrprop=''

    if (associated(g3d%cugd_ttens))  &
         call InsertVTab (g3d%cugd_ttens     ,g3dm%cugd_ttens      &
         ,ng, npts, imean,  &
         'TTENS :3:hist:anal:mpti:mpt3'//trim(arrprop))

    if (associated(g3d%cugd_qvtens))  &
         call InsertVTab (g3d%cugd_qvtens    ,g3dm%cugd_qvtens     &
         ,ng, npts, imean,  &
         'QVTTENS :3:hist:anal:mpti:mpt3'//trim(arrprop))

    if (associated(g3d%thsrc))  &
         call InsertVTab (g3d%thsrc     ,g3dm%thsrc      &
         ,ng, npts, imean,  &
         'THSRC :3:hist:anal:mpti:mpt3'//trim(arrprop))

    if (associated(g3d%rtsrc))  &
         call InsertVTab (g3d%rtsrc     ,g3dm%rtsrc     &
         ,ng, npts, imean,  &
         'RTSRC :3:hist:anal:mpti:mpt3'//trim(arrprop))

    !- this array does not need to be parallelized (only column)
    if (associated(g3d%clsrc))  &
         call InsertVTab (g3d%clsrc     ,g3dm%clsrc     &
         ,ng, npts, imean,  &
         'CLSRC :3:hist:anal:mpti:mpt3')
 
    if (associated(g3d%nlsrc))  &
         call InsertVTab (g3d%nlsrc     ,g3dm%nlsrc     &
         ,ng, npts, imean,  &
         'NLSRC :3:hist:anal:mpti:mpt3')
    if (associated(g3d%nisrc))  &
         call InsertVTab (g3d%nisrc     ,g3dm%nisrc     &
         ,ng, npts, imean,  &
         'NISRC :3:hist:anal:mpti:mpt3')


    if(imomentum==1 .and. nnqparm(ng) >= 4) then
    !- these arrays does not need to be parallelized (only column)
      if (associated(g3d%usrc))  &
         call InsertVTab (g3d%usrc     ,g3dm%usrc     &
         ,ng, npts, imean,  &
         'USRC :3:hist:anal:mpti:mpt3')
      if (associated(g3d%vsrc))  &
         call InsertVTab (g3d%vsrc     ,g3dm%vsrc     &
         ,ng, npts, imean,  &
         'VSRC :3:hist:anal:mpti:mpt3')
    endif

 
    if (associated(g3d%buoysrc))  &
 
         call InsertVTab (g3d%buoysrc     ,g3dm%buoysrc     &
         ,ng, npts, imean,  &
         'BUOYSRC :3:hist:anal:mpti:mpt3'//trim(arrprop))

 !-----------
   if (associated(g3d%cnv_tr))  &
         call InsertVTab (g3d%cnv_tr     ,g3dm%cnv_tr     &
         ,ng, npts, imean,  &
         'CNV_TR :3:hist:anal:mpti:mpt3:mpt1'//trim(arrprop))
 !--------------------
   if (associated(g3d%qexcp)) &
         call InsertVTab (g3d%qexcp   ,g3dm%qexcp   &
        ,ng, npts, imean, 'QEXCP :3:hist:anal:mpti:mpt3:mpt1'//trim(arrprop))
   if (associated(g3d%hexcp)) &
         call InsertVTab (g3d%hexcp   ,g3dm%hexcp   &
        ,ng, npts, imean, 'HEXCP :3:hist:anal:mpti:mpt3:mpt1'//trim(arrprop))



    if (associated(g3d%mupsh))  &
         call InsertVTab (g3d%mupsh    ,g3dm%mupsh     &
         ,ng, npts, imean,  &
         'ZMFSH :3:hist:anal:mpti:mpt3')
    if (associated(g3d%mupdp))  &
         call InsertVTab (g3d%mupdp    ,g3dm%mupdp     &
         ,ng, npts, imean,  &
         'ZMFUP :3:hist:anal:mpti:mpt3')
    if (associated(g3d%mdddp))  &
         call InsertVTab (g3d%mdddp    ,g3dm%mdddp    &
         ,ng, npts, imean,  &
         'ZMFDD :3:hist:anal:mpti:mpt3')
    if (associated(g3d%mupmd))  &
         call InsertVTab (g3d%mupmd    ,g3dm%mupmd     &
         ,ng, npts, imean,  &
         'ZMFMD :3:hist:anal:mpti:mpt3')
    
    if (associated(g3d%eupdp))  &
         call InsertVTab (g3d%eupdp    ,g3dm%eupdp     &
         ,ng, npts, imean,  &
         'ENTUP :3:hist:anal:mpti:mpt3')
    if (associated(g3d%dupdp))  &
         call InsertVTab (g3d%dupdp    ,g3dm%dupdp     &
         ,ng, npts, imean,  &
         'DETUP :3:hist:anal:mpti:mpt3')
  end subroutine filltab_grell3
!-------------------------------------------------------------

subroutine CUPARM_GRELL3_CATT(OneGrid, iens,iinqparm,iinshcu)

    USE mem_radiate, ONLY: &
         ilwrtyp, iswrtyp        ! INTENT(IN)
    use ModGrid, only: &
         Grid
!ML -- In case you want to output massflux
    use mem_stilt         , only: imassflx

    use mem_jules

  implicit none

  include "i8.h"
  integer, intent(IN) :: iens,iinqparm,iinshcu
  type(Grid), pointer :: OneGrid ! intent(in)
  integer :: i,j,k
  real :: grid_length,theta2temp,tau_cp, src_cnvtr, snk_cnvtr

  REAL, DIMENSION( mxp , myp ) :: aot500,temp2m

!----GF-GEOS-5 -------------------------------->
  INTEGER, PARAMETER :: mtp=1 !tmp
  INTEGER, DIMENSION(mzp)   :: flip
  INTEGER, DIMENSION(mxp,myp) :: kpbl,do_this_column
  REAL   , DIMENSION(mtp)   :: FSCAV_INT
  REAL   , DIMENSION(mxp,myp) :: CNV_FRC,AA0,AA1,AA1_ADV,AA1_RADPBL,AA1_BL,AA1_CIN &
                              ,AA2,AA3,TAU_BL,TAU_EC
  REAL   , DIMENSION(mzp,mxp,myp ) ::  zm3d	&
 				      ,zt3d	&
        	          ,dm3d	&
 				      ,up	&
 				      ,vp	&
 				      ,wp	&
 				      ,rvap	&
 				      ,temp	&
 				      ,press	&
                      !,rhenv   &
                      , gsf_t	& ! grid-scale forcing for temp
 				      , gsf_q	& ! grid-scale forcing fo rv
 				      ,sgsf_t	& ! sub-grid scale forcing for temp
 				      ,sgsf_q	  ! sub-grid scale forcing for rv

  REAL,  DIMENSION(mzp , mxp, myp ) ::          &
		      advf_t	  &
		     ,REVSU_GF    &
		     ,PRFIL_GF    & 
		     ,VAR3d_aGF,VAR3d_bGF,VAR3d_cGF,VAR3d_dGF

  REAL,  DIMENSION(nmp, mzp , mxp, myp ) ::       &
					 mp_ice   &
					,mp_liq   &
					,mp_cf    

  REAL,  DIMENSION(nmp, mzp , mxp, myp ) ::       &
                                         SUB_MPQI & ! subsidence transport applied to ice mix ratio
                                        ,SUB_MPQL & ! subsidence transport applied to cloud mix ratio
                                        ,SUB_MPCF   ! subsidence transport applied to cloud fraction
  
  REAL   , DIMENSION(mxp,myp)  :: var2d,col_sat,stochastic_sig

  REAL   ,DIMENSION(mtp,mzp,mxp,myp)  :: SRC_CHEM
  !---
  REAL   ,DIMENSION(mxp,myp,mzp,mtp)  :: TRACER  !geos-5 data structure
  !---
  REAL   , DIMENSION(mxp,myp)  :: dx2d ,lons,lats,sfc_press,xland, tke_pbl

  INTEGER :: kr,n,i1,i2,j1,j2
  INTEGER, DIMENSION(mxp,myp,maxiens) ::     &
  	     ierr4d_tmp 		     &
  	    ,jmin4d_tmp 		     &
  	    ,klcl4d_tmp 		     &
  	    ,k224d_tmp  		     &
  	    ,kbcon4d_tmp		     &
  	    ,ktop4d_tmp 		     &
  	    ,kstabi4d_tmp		     &
  	    ,kstabm4d_tmp

  REAL,DIMENSION(mxp,myp,maxiens)     ::     &
  	     cprr4d_tmp 		     &
  	    ,xmb4d_tmp  		     &
  	    ,edt4d_tmp  		     &
  	    ,pwav4d_tmp 		     &
  	    ,sigma4d_tmp
  REAL,DIMENSION(mxp,myp,mzp,maxiens) ::     &
  	     pcup5d_tmp 		     &
  	    ,up_massentr5d_tmp    &
  	    ,up_massdetr5d_tmp    &
  	    ,dd_massentr5d_tmp    &
  	    ,dd_massdetr5d_tmp    &
  	    ,zup5d_tmp  		     &
  	    ,zdn5d_tmp  		     &
  	    ,prup5d_tmp 		     &
  	    ,prdn5d_tmp 		     &
  	    ,clwup5d_tmp		     &
  	    ,tup5d_tmp  		     &
  	    ,conv_cld_fr5d_tmp

!if(initial.eq.2.and.time.lt.cptime) return
!if(initial.eq.2.and.time.lt.dtlt) return

!if(time < 10800.) return
 
 !--- convection tracer 
 if( nnqparm(ngrid) == 8  .and. convection_tracer == 1) then  

      if(use_gustiness == 1) call sgs_coldpool_sfx_and_diff()

      call adv_convection_tracer()
 endif


 if(mod(time,confrq) < dtlt  .or. time < 0.01 .or. abs(time-cptime) < 0.01) then
 
    !-start convective transport of tracers
    iruncon=1
    g3d_g(ngrid)%thsrc       = 0.0
    g3d_g(ngrid)%rtsrc       = 0.0
    g3d_g(ngrid)%clsrc       = 0.0
    g3d_g(ngrid)%cugd_ttens  = 0.0
    g3d_g(ngrid)%cugd_qvtens = 0.0
    g3d_g(ngrid)%mupsh       = 0.0
    g3d_g(ngrid)%mupdp       = 0.0
    g3d_g(ngrid)%mdddp       = 0.0
    g3d_g(ngrid)%mupmd       = 0.0
    g3d_g(ngrid)%eupdp       = 0.0
    g3d_g(ngrid)%dupdp       = 0.0
    cuparm_g(ngrid)%conprr   = 0.0
  
    if(liq_ice_number_conc > 0) then 
      g3d_g(ngrid)%nlsrc     = 0.0
      g3d_g(ngrid)%nisrc     = 0.0
    endif
    if(imomentum == 1 .and. nnqparm(ngrid) >= 4) then
      g3d_g(ngrid)%usrc      = 0.0
      g3d_g(ngrid)%vsrc      = 0.0
    endif

    if(nnqparm(ngrid) == 8) then
      g3d_g(ngrid)%buoysrc   = 0.0
    endif

    ishallow_g3=0
    
    if(i_forcing /= 1) then
      !call check(mzp * myp * mzp,tend%THT ,cuforc_g(ngrid)%lsfth ,tend%RTT  ,cuforc_g(ngrid)%lsfrt)
      call atob(mxp * myp * mzp,tend%THT  ,cuforc_g(ngrid)%lsfth     )
      call atob(mxp * myp * mzp,tend%RTT  ,cuforc_g(ngrid)%lsfrt     )
    endif

    !- converting WRF setting to BRAMS
    ids=1   ;ide=mxp ;jds=1   ;jde=myp ;kds=1; kde=mzp
    ims=1   ;ime=mxp ;jms=1   ;jme=myp ;kms=1; kme=mzp
    ips=ia+1;ipe=iz-2;jps=ja+1;jpe=jz-2;kps=1; kpe=mzp
    its=ia  ;ite=iz  ;jts=ja  ;jte=jz  ;kts=1; kte=mzp-1

    grid_length=sqrt(deltaxn(ngrid)*deltayn(ngrid))

!-------------------------------------------------------------

    if(iinqparm==3) then  ! G3d scheme
      !
      !- lateral spreading
      if(g3d_spread == 0 )cugd_avedx=1
      if(g3d_spread == 1 )cugd_avedx=3

      if(ilwrtyp==4 .or. iswrtyp==4) THEN
        aot500(:,:)=carma(ngrid)%aot(:,:,11)
      else
        aot500(:,:)=0.0
      end if

      CALL G3DRV( mynum,i0,j0,time          &
              ,dtlt                         & !
              ,grid_length                  & !
              ,autoconv                     & !
              ,aerovap                      & !
              ,basic_g(ngrid)%dn0           & !
              ,cuparm_g(ngrid)%CONPRR       & !
              ,basic_g(ngrid)%up            & !
              ,basic_g(ngrid)%vp            & !
              ,basic_g(ngrid)%theta         & !
              ,basic_g(ngrid)%thp           & !
              ,basic_g(ngrid)%pp            & !
              ,basic_g(ngrid)%pi0           & !
              ,basic_g(ngrid)%wp            & !
              ,basic_g(ngrid)%rv            & !
              ,grid_g(ngrid)%RTGT           & !
              ,tend%PT                      & !
              ,XL			                    & !
              ,CP			    & !
              ,G			    & !
              ,rm                           &
              ,p00                          &
              ,cpor                         & !
              ,g3d_ens_g(apr_gr,ngrid)%apr  &
              ,g3d_ens_g(apr_w ,ngrid)%apr  &
              ,g3d_ens_g(apr_mc,ngrid)%apr  &
              ,g3d_ens_g(apr_st,ngrid)%apr  &
              ,g3d_ens_g(apr_as,ngrid)%apr  &
              ,g3d_g(ngrid)%xmb_deep        &
              ,g3d_g(ngrid)%xmb_shallow     &
!
              ,g3d_ens_g(apr_gr,ngrid)%weight &
              ,g3d_ens_g(apr_w ,ngrid)%weight &
              ,g3d_ens_g(apr_mc,ngrid)%weight &
              ,g3d_ens_g(apr_st,ngrid)%weight &
              ,g3d_ens_g(apr_as,ngrid)%weight &
!
              ,training &
!
              ,grid_g(ngrid)%topt              &
              ,leaf_g(ngrid)%patch_area        &
              ,npatch                          &
              ,radiate_g(ngrid)%rshort         &

              ,cugd_avedx					&
              ,imomentum          				&
              ,ensdim,maxiens,maxens,maxens2,maxens3,icoic      &
              ,ishallow_g3                                      &
              ,ids,ide, jds,jde, kds,kde                        &
              ,ims,ime, jms,jme, kms,kme                        &
              ,ips,ipe, jps,jpe, kps,kpe                        &
              ,its,ite, jts,jte, kts,kte                        &
              ,g3d_g(ngrid)%THSRC      & ! temp tendency
              ,g3d_g(ngrid)%RTSRC      & ! rv tendency
              ,g3d_g(ngrid)%CLSRC      & ! cloud/ice tendency
              ,g3d_g(ngrid)%cugd_ttens    &
              ,g3d_g(ngrid)%cugd_qvtens   &
              ! forcings -  for deep/shallow
              ,cuforc_g(ngrid)%	lsfth    & ! forcing for theta deep
              ,cuforc_g(ngrid)%	lsfrt    & ! forcing for rv deep
              ,cuforc_sh_g(ngrid)%lsfth   & ! forcing for theta shallow
              ,cuforc_sh_g(ngrid)%lsfrt   & ! forcing for rv shallow
              ,level                      &
              ,micro_g(ngrid)%rcp         & ! liquid water
              ,micro_g(ngrid)%rrp         & ! pristine
              ,micro_g(ngrid)%rpp         &
              ,micro_g(ngrid)%rsp         &
              ,micro_g(ngrid)%rap         &
              ,micro_g(ngrid)%rgp         &
              ,micro_g(ngrid)%rhp         &
              ,aot500                     & ! aot at 500nm
              )
!
!- exchange border information for parallel run
   if( g3d_spread == 1 .or. g3d_smoothh == 1) then
      call PostRecvSendMsgs(OneGrid%SendG3D, OneGrid%RecvG3D)
      call WaitRecvMsgs    (OneGrid%SendG3D, OneGrid%RecvG3D)
   endif
!
!
!- call routine to do the lateral spread, smooths and limiters/fixers
   CALL conv_grell_spread3d_brams(mzp,mxp,myp,ia,iz,ja,jz,dtlt,level,cugd_avedx&
	      ,XL			    &
	      ,CP			    &
	      ,G			    &
	      ,rm                           &
	      ,p00                          &
	      ,cpor                         &
	      ,cuparm_g(ngrid)%CONPRR       &!preci rate
              ,basic_g(ngrid)%theta         &
              ,basic_g(ngrid)%thp           &
              ,basic_g(ngrid)%pp            &
              ,basic_g(ngrid)%pi0           &
	           ,basic_g(ngrid)%rv            &
              ,tend%PT                      &
	           ,micro_g(ngrid)%rcp           & ! liquid water
    	        ,micro_g(ngrid)%rrp           & ! pristine
    	        ,micro_g(ngrid)%rpp           &
	           ,micro_g(ngrid)%rsp           &
    	        ,micro_g(ngrid)%rap           &
	           ,micro_g(ngrid)%rgp           &
    	        ,micro_g(ngrid)%rhp           &
!	      
              ,g3d_g(ngrid)%THSRC           & ! temp tendency
              ,g3d_g(ngrid)%RTSRC           & ! rv tendency
              ,g3d_g(ngrid)%CLSRC           & ! cloud/ice tendency
              ,g3d_g(ngrid)%cugd_ttens      &
              ,g3d_g(ngrid)%cugd_qvtens     &
              ,g3d_ens_g(apr_gr,ngrid)%apr  &
              ,g3d_ens_g(apr_w ,ngrid)%apr  &
              ,g3d_ens_g(apr_mc,ngrid)%apr  &
              ,g3d_ens_g(apr_st,ngrid)%apr  &
              ,g3d_ens_g(apr_as,ngrid)%apr  )

!-------------------------------------------------------------
 elseif(iinqparm==5) then  ! GF 2014 scheme
   !
   !- no lateral spreading
   cugd_avedx=1

   if(ilwrtyp==4 .or. iswrtyp==4) THEN
   	aot500(:,:)=carma(ngrid)%aot(:,:,11)
   else
   	aot500(:,:)=0.0
   endif

   if(iinshcu == 3) ishallow_g3=1

   CALL GFDRV( CCATT                        &
              ,mgmxp,mgmyp,mgmzp,ngrid,ngrids_cp,iens &
              ,mynum,i0,j0,time,mzp,mxp,myp &
              ,dtlt                           & !
              ,grid_length                    & !
              ,autoconv                       & ! Const
              ,aerovap                        & ! Const
              ,basic_g(ngrid)%dn0             & !3d ok
              ,cuparm_g(ngrid)%CONPRR         & !2d ok
              ,basic_g(ngrid)%up              & !3d ok
              ,basic_g(ngrid)%vp              & !3d ok
              ,basic_g(ngrid)%theta           & !3d ok
              ,basic_g(ngrid)%thp             & !3d ok
              ,basic_g(ngrid)%pp              & !3d ok
              ,basic_g(ngrid)%pi0             & !3d ok
              ,basic_g(ngrid)%wp              & !3d ok
              ,basic_g(ngrid)%rv              & !3d ok
              ,basic_g(ngrid)%rtp             & !3d ok
              ,grid_g(ngrid)%rtgt             & !2d ok
              ,tend%pt                        & !3d !*** borda
              ,xl                             & ! Const
              ,cp                             & ! Const
              ,g                              & ! Const
              ,rm                             & ! Const
              ,p00                            & ! Const
              ,cpor                           & ! Const
              ,rgas                           & ! Const
              ,zmn(:,ngrid)                   & !
              ,ztn(:,ngrid)                   & !
              ,g3d_ens_g(apr_gr,ngrid)%apr    & !2d ok
              ,g3d_ens_g(apr_w ,ngrid)%apr    & !2d ok
              ,g3d_ens_g(apr_mc,ngrid)%apr    & !2d ok
              ,g3d_ens_g(apr_st,ngrid)%apr    & !2d ok
              ,g3d_ens_g(apr_as,ngrid)%apr    & !2d ok
              ,g3d_g(ngrid)%xmb_deep          & !2d ok
              ,g3d_g(ngrid)%err_deep          & !2d ok
              ,g3d_g(ngrid)%xmb_shallow       & !2d ok
              ,g3d_ens_g(apr_gr,ngrid)%weight & !2d ok
              ,g3d_ens_g(apr_w ,ngrid)%weight & !2d ok
              ,g3d_ens_g(apr_mc,ngrid)%weight & !2d ok
              ,g3d_ens_g(apr_st,ngrid)%weight & !2d ok
              ,g3d_ens_g(apr_as,ngrid)%weight & !2d ok
              ,training                       & !
              ,grid_g(ngrid)%topt             & !2d ok
              ,leaf_g(ngrid)%patch_area       & !3d *** Borda
              ,npatch                         & !
              ,radiate_g(ngrid)%rshort        & !2d ok
              ,cugd_avedx                     & !
              ,imomentum                      & !
              ,ensdim_g3d                     & !
              ,maxiens                        & !
              ,maxens_g3d                     & !
              ,maxens2_g3d                    & !
              ,maxens3_g3d                    & !
              ,icoic                          & !
              ,ishallow_g3                    & !
              ,ids,ide, jds,jde, kds,kde      & !
              ,ims,ime, jms,jme, kms,kme      & !
              ,ips,ipe, jps,jpe, kps,kpe      & !
              ,its,ite, jts,jte, kts,kte      & !
              ,g3d_g(ngrid)%THSRC             & !3d ok ! temp tendency
              ,g3d_g(ngrid)%RTSRC             & !3d ok ! rv tendency
              ,g3d_g(ngrid)%CLSRC             & !3d ok ! cloud/ice tendency
              ,g3d_g(ngrid)%cugd_ttens        & !3d ok
              ,g3d_g(ngrid)%cugd_qvtens       & !3d ok
              ,cuforc_g(ngrid)%	lsfth        & !3d *** borda forcing for theta deep
              ,cuforc_g(ngrid)%	lsfrt        & !3d *** borda forcing for rv deep
              ,cuforc_sh_g(ngrid)%lsfth       & !3d *** borda forcing for theta shallow
              ,cuforc_sh_g(ngrid)%lsfrt       & !3d *** borda forcing for rv shallow
              ,level                          &
              ,micro_g(ngrid)%rcp             & !3d ok ! liquid water
              ,micro_g(ngrid)%rrp             & !3d ok ! pristine
              ,micro_g(ngrid)%rpp             & !3d ok
              ,micro_g(ngrid)%rsp             & !3d ok
              ,micro_g(ngrid)%rgp             & !3d ok
              ,aot500                         &! aot at 500nm
              ,turb_g(ngrid)%sflux_r          & !2d *** borda
              ,turb_g(ngrid)%sflux_t          & !2d *** borda
              ,turb_g(ngrid)%tkep             & !3d ok
              ,TKMIN                          &
              ,akmin(ngrid)                   &
!- for convective transport-start
              ,ierr4d                         & !4d *** borda
              ,jmin4d                         & !4d ok
              ,kdet4d                         & !4d ok
              ,k224d                          & !4d ok
              ,kbcon4d                        & !4d ok
              ,ktop4d                         & !4d ok
              ,kpbl4d                         & !4d ok
              ,kstabi4d                       & !4d ok
              ,kstabm4d                       & !4d ok
              ,xmb4d                          & !4d ok
              ,edt4d                          & !4d ok
              ,pwav4d                         & !4d ok
              ,pcup5d                         & !5d ok
              ,up_massentr5d                  & !5d ok
              ,up_massdetr5d                  & !5d ok
              ,dd_massentr5d                  & !5d ok
              ,dd_massdetr5d                  & !5d ok
              ,zup5d                          & !5d *** ERRO !!!!!
              ,zdn5d                          & !5d ok
              ,prup5d                         & !5d ok
              ,prdn5d                         & !5d ok
              ,clwup5d                        & !5d ok
              ,tup5d                          & !5d ok
!- for convective transport- end
          			     )
!-------------------------------------------------------------
 elseif(iinqparm==6) then  ! GF 2015 scheme
   !
   !- no lateral spreading
   cugd_avedx=1

   if(ilwrtyp==4 .or. iswrtyp==4) THEN
   	aot500(:,:)=carma(ngrid)%aot(:,:,11)
   else
   	aot500(:,:)=0.0
   endif

   if(isfcl == 5) then
     temp2m(:,:) = jules_g(ngrid)%t2mj(:,:)
   else
     temp2m(:,:) =0.5*( basic_g(ngrid)%theta(1,:,:)* &
                       (basic_g(ngrid)%pp(1,:,:)+basic_g(ngrid)%pi0(1,:,:))/cp + &
                        basic_g(ngrid)%theta(2,:,:)*&
		                 (basic_g(ngrid)%pp(2,:,:)+basic_g(ngrid)%pi0(2,:,:))/cp )
   endif

   if(iinshcu == 3) ishallow_g3=1

   CALL GFDRV2(mgmxp,mgmyp,mgmzp,ngrid,ngrids_cp,iens &
              ,mynum,i0,j0,time,mzp,mxp,myp &
              ,dtlt         		    & !
              ,grid_length                  & !
              ,autoconv                     & !
              ,aerovap                      & !
              ,basic_g(ngrid)%dn0           & !
	           ,cuparm_g(ngrid)%CONPRR       & !
              ,basic_g(ngrid)%up            & !
              ,basic_g(ngrid)%vp            & !
              ,basic_g(ngrid)%theta         & !
              ,basic_g(ngrid)%thp           & !
              ,basic_g(ngrid)%pp            & !
              ,basic_g(ngrid)%pi0           & !
	           ,basic_g(ngrid)%wp            & !
	           ,basic_g(ngrid)%rv            & !
	           ,basic_g(ngrid)%rtp           & !
              ,grid_g(ngrid)%RTGT           & !
              ,tend%PT                      & !
	           ,XL			                    & !
	           ,CP			                    & !
	           ,G			                    & !
	           ,rm                           &
	           ,p00                          &
	           ,cpor                         & !
	           ,rgas                         & !
	           ,zmn(:,ngrid)                 & !
	           ,ztn(:,ngrid)                 & !
              ,g3d_ens_g(apr_gr,ngrid)%apr  &
              ,g3d_ens_g(apr_w ,ngrid)%apr  &
              ,g3d_ens_g(apr_mc,ngrid)%apr  &
              ,g3d_ens_g(apr_st,ngrid)%apr  &
              ,g3d_ens_g(apr_as,ngrid)%apr  &
              ,g3d_g(ngrid)%xmb_deep        &
              ,g3d_g(ngrid)%err_deep        &
              ,g3d_g(ngrid)%xmb_shallow     &
	           ,g3d_ens_g(apr_gr,ngrid)%weight &
              ,g3d_ens_g(apr_w ,ngrid)%weight &
              ,g3d_ens_g(apr_mc,ngrid)%weight &
              ,g3d_ens_g(apr_st,ngrid)%weight &
              ,g3d_ens_g(apr_as,ngrid)%weight &
              ,training &
	           ,grid_g(ngrid)%topt              &
              ,leaf_g(ngrid)%patch_area        &
	           ,npatch                          &
              ,radiate_g(ngrid)%rshort         &
	           ,cugd_avedx					        &
	           ,imomentum          				  &
              ,ensdim_g3d,maxiens,maxens_g3d,maxens2_g3d,maxens3_g3d,icoic      &
              ,ishallow_g3                                      &
	           ,ids,ide, jds,jde, kds,kde                        &
              ,ims,ime, jms,jme, kms,kme                        &
              ,ips,ipe, jps,jpe, kps,kpe                        &
              ,its,ite, jts,jte, kts,kte                        &
              ,g3d_g(ngrid)%THSRC      & ! temp tendency
              ,g3d_g(ngrid)%RTSRC      & ! rv tendency
              ,g3d_g(ngrid)%CLSRC      & ! cloud/ice tendency
              ,g3d_g(ngrid)%USRC       & ! U tendency
              ,g3d_g(ngrid)%VSRC       & ! V tendency
              ,g3d_g(ngrid)%MUPDP      & ! updraft mass flux
	           ,cuforc_g(ngrid)%lsfth   & ! forcing for theta deep
	           ,cuforc_g(ngrid)%lsfrt   & ! forcing for rv deep
	           ,cuforc_sh_g(ngrid)%lsfth  & ! forcing for theta shallow
	           ,cuforc_sh_g(ngrid)%lsfrt  & ! forcing for rv shallow
              ,level                     &
	           ,micro_g(ngrid)%rcp        & ! liquid water
	           ,aot500                    &! aot at 500nm
	           ,temp2m                    &! aot at 500nm
  	           ,turb_g(ngrid)%sflux_r     &
              ,turb_g(ngrid)%sflux_t     &
              ,turb_g(ngrid)%tkep        &
              ,TKMIN                     &
	           ,akmin(ngrid)              &
	           ,do_cupar_mcphys_coupling  &
!- for convective transport-start
              ,ierr4d  		  &
	      ,jmin4d  		     &
	      ,kdet4d  		     &
	      ,k224d	           &
	      ,kbcon4d 		     &
	      ,ktop4d  		     &
	      ,kpbl4d  		     &
	      ,kstabi4d		     &
	      ,kstabm4d		     &
	      ,xmb4d		        &
	      ,edt4d		        &
	      ,pwav4d		        &
	      ,pcup5d  		     &
         ,up_massentr5d	     &
	      ,up_massdetr5d	     &
	      ,dd_massentr5d	     &
	      ,dd_massdetr5d	     &
	      ,zup5d		        &
	      ,zdn5d   		     &
	      ,prup5d  		     &
	      ,prdn5d  		     &
	      ,clwup5d 		     &
	      ,tup5d   		     &
!- for convective transport- end
          			     )

!--------------------------------------------------------------------------------------------------------------------------

 elseif(iinqparm==8) then  ! GF 2021 scheme

    IF(read_GF_ConvPar_nml) THEN
      !-- read the GF namelist
      call GF_convpar_init(mynum)
      read_GF_ConvPar_nml = .false.
      !scalar_g(1,1)%sclp(:,:,:) = 0.  ! Check this later
      !if(mynum==7)  g3d_g(ngrid)%cnv_tr=10. ! only for debugging
    ENDIF

    

    !--- these arrays must be set to zero every timestep.
    ierr4d_tmp        = 0.0
    jmin4d_tmp        = 0.0
    klcl4d_tmp        = 0.0
    k224d_tmp	       = 0.0
    kbcon4d_tmp       = 0.0
    ktop4d_tmp        = 0.0
    kstabi4d_tmp      = 0.0
    kstabm4d_tmp      = 0.0
    cprr4d_tmp        = 0.0
    xmb4d_tmp	       = 0.0
    edt4d_tmp	       = 0.0
    pwav4d_tmp        = 0.0
    sigma4d_tmp       = 0.0
    pcup5d_tmp        = 0.0
    up_massentr5d_tmp = 0.0
    up_massdetr5d_tmp = 0.0
    dd_massentr5d_tmp = 0.0
    dd_massdetr5d_tmp = 0.0
    zup5d_tmp	       = 0.0
    zdn5d_tmp	       = 0.0
    prup5d_tmp        = 0.0
    prdn5d_tmp        = 0.0
    clwup5d_tmp       = 0.0
    tup5d_tmp	       = 0.0
    conv_cld_fr5d_tmp = 0.0
    CNV_FRC (:,:)     = 0.0
    TRACER  (:,:,:,:) = 0.0 
    SRC_CHEM(:,:,:,:) = 0.0 

    call set_index_loops( ims,ime, jms,jme, kms,kme,    &
                          its,ite, jts,jte, kts,kte,    &
                          mxp,myp,mzp                   )

    FSCAV_INT(:)   = 0.1 

    do k=1,mzp
      flip   (k)   = k
    enddo
    if( idiffk(ngrid) /= 2 .and. idiffk(ngrid) /= 3) then 
        if(idiffk(ngrid) == 7 ) then          
          kpbl (:,:) = max(1,nint(turb_g(ngrid)%kpbl(:,:)))
        else
          do j=1,myp
            do i=1,mxp
               call get_zi_gf2018(mzp,tkmin,turb_g(ngrid)%tkep(:,i,j),zmn(:,ngrid) &
                                 ,grid_g(ngrid)%rtgt(i,j)                          &
	                              ,grid_g(ngrid)%topt(i,j),kpbl(i,j) )
	            kpbl (i,j) = max(1,min(kpbl (i,j),mzp-1))
            enddo
         enddo
        end if
        !--- mean pbl TKE for new shallow convection mass flux closure
        do j=1,myp
            do i=1,mxp
             call get_mean_tke(mzp,tkmin,turb_g(ngrid)%tkep(:,i,j) ,grid_g(ngrid)%rtgt(i,j) &
                              ,kpbl(i,j),basic_g(ngrid)%dn0(:,i,j), tke_pbl(i,j) )
            enddo
         enddo
         !print*,"tkemean=",maxval(tke_pbl),minval(tke_pbl),maxval(kpbl),minval(kpbl)
    else
      kpbl       = 5  ! check later (introduce better formulation for Zi )
      tke_pbl(:,:) = tkmin
    endif
    !
    do j=1,myp
     do i=1,mxp
         do_this_column(i,j)=0

         do k=1,mzp-1
	        kr=k+1
           zm3d   (k,i,j) = zmn(kr,ngrid)*grid_g(ngrid)%rtgt(i,j) !m - height above local terrain
           zt3d   (k,i,j) = ztn(kr,ngrid)*grid_g(ngrid)%rtgt(i,j) !m
           dm3d   (k,i,j) = basic_g(ngrid)%dn0  (kr,i,j) !kg/m3
           rvap   (k,i,j) = basic_g(ngrid)%rv   (kr,i,j) !kg/kg
           
	        theta2temp     = (basic_g(ngrid)%pp(kr,i,j)+basic_g(ngrid)%pi0(kr,i,j))/cp   !K
	        temp   (k,i,j) = basic_g(ngrid)%theta(kr,i,j)* theta2temp
           press  (k,i,j) = ((basic_g(ngrid)%pp(kr,i,j)+basic_g(ngrid)%pi0(kr,i,j))/cp)**cpor*p00 !Pa
          !rhenv  (k,i,j) = 100.*min(1., max (0., rvap(k,i,j)/ rs(press(k,i,j), temp(k,i,j))))

           up     (k,i,j) = basic_g(ngrid)%up(kr,i,j) !m/s
           vp     (k,i,j) = basic_g(ngrid)%vp(kr,i,j) !m/s
           wp     (k,i,j) = basic_g(ngrid)%wp(kr,i,j)*(-g*basic_g(ngrid)%dn0(kr,i,j)) ! omega Pa/s
            gsf_t (k,i,j) = (cuforc_g   (ngrid)%lsfth(kr,i,j) + radiate_g(ngrid)%fthrd(kr,i,j))* theta2temp ! Adv+Rad, K/s
            gsf_q (k,i,j) =  cuforc_g   (ngrid)%lsfrt(kr,i,j)              !kg/kg/s  Adv only
           sgsf_t (k,i,j) =  cuforc_sh_g(ngrid)%lsfth(kr,i,j) * theta2temp !K/s     PBL only 
           sgsf_q (k,i,j) =  cuforc_sh_g(ngrid)%lsfrt(kr,i,j)              !kg/kg/s PBL only 
           advf_t (k,i,j) =  cuforc_g   (ngrid)%lsfth(kr,i,j) * theta2temp ! advection only, see 'prepare_lsl' routine
    enddo;enddo;enddo

    IF(APPLY_SUB_MP == 1) THEN
      do j=1,myp
       do i=1,mxp
        do k=1,mzp-1
	      kr=k+1
	      mp_ice   (:,k,i,j) = 0. ! microg%ice(kr,i,j) in the future includes ice mix ratio 
	      mp_liq   (:,k,i,j) = 0. ! microg%liq(kr,i,j) in the future includes liq ratio 
	      mp_cf    (:,k,i,j) = 0. ! cloud fraction  
      enddo; enddo; enddo
    ENDIF

    do j=1,myp
     do i=1,mxp
       sfc_press(i,j) = 0.5*( ((basic_g(ngrid)%pp(1,i,j)+basic_g(ngrid)%pi0(1,i,j))/cp)**cpor*p00 +  &
                              ((basic_g(ngrid)%pp(2,i,j)+basic_g(ngrid)%pi0(2,i,j))/cp)**cpor*p00 ) !Pa

       xland(i,j) = leaf_g(ngrid)%patch_area(i,j,1); if(xland(i,j) < 0.98) xland(i,j)=0. ! water = 1, land < 1
       lons (i,j) = grid_g(ngrid)%glon(i,j)*3.14159/180. !- convert to rad
       lats (i,j) = grid_g(ngrid)%glat(i,j)*3.14159/180. !- convert to rad
       dx2d (i,j) = sqrt(1./(grid_g(ngrid)%dxt(i,j)*grid_g(ngrid)%dyt(i,j)))
       col_sat       (i,j) = 0.
       stochastic_sig(i,j) = 1.
       
       !-- stores precip of the previous timestep 
       AA2(i,j)   = g3d_ens_g(1,ngrid)%apr(i,j) !cprr4d_tmp(:,:,deep) -> previous time step
    enddo;enddo

    if(ilwrtyp==4 .or. iswrtyp==4) THEN
       aot500(:,:)=carma(ngrid)%aot(:,:,11)
     else
       aot500(:,:)=0.0
    endif

    if(isfcl == 5 .and.  time > dtlt ) then
          temp2m(:,:) = jules_g(ngrid)%t2mj(:,:) !K
    else
          temp2m(:,:) =0.5*(basic_g(ngrid)%theta(1,:,:)* &
                           (basic_g(ngrid)%pp(1,:,:)+basic_g(ngrid)%pi0(1,:,:))/cp + &
                            basic_g(ngrid)%theta(2,:,:)*&
		                     (basic_g(ngrid)%pp(2,:,:)+basic_g(ngrid)%pi0(2,:,:))/cp ) !Kelvin
    endif

   
    if(use_gustiness==2 .and. convection_tracer==1)  then
       call sgs_coldpool_sfx_and_diff()
    else
       g3d_g(ngrid)%sgs_shf(:,:) = turb_g(ngrid)%sflux_t(:,:)
       g3d_g(ngrid)%sgs_lhf(:,:) = turb_g(ngrid)%sflux_r(:,:)
    endif

    !- call the driver routine to apply the parameterization
    CALL GF_GEOS5_DRV(mxp,myp,mzp,mtp ,nmp, time, itime1 &
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
                    !,turb_g(ngrid)%sflux_r &
                    !,turb_g(ngrid)%sflux_t &
                     ,g3d_g(ngrid)%sgs_lhf  &
                     ,g3d_g(ngrid)%sgs_shf  &
                     ,g3d_g(ngrid)%qexcp    &
                     ,g3d_g(ngrid)%hexcp    &
                     ,g3d_g(ngrid)%wlpool   &
                     ,grid_g(ngrid)%topt    &
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
                     ,g3d_g(ngrid)%cnv_tr & 
                     , gsf_t   & ! forcing for theta adv+rad
		             , gsf_q   & ! forcing for rv    adv
                     ,advf_t   &
		             ,sgsf_t   & ! forcing for theta pbl
 		             ,sgsf_q   & ! forcing for rv    pbl
                     !---- output ----
                     ,cuparm_g(ngrid)%conprr  &
                     ,g3d_g(ngrid)%lightn_dens& 
                     ,g3d_g(ngrid)%rh_dicy_fct&
                     ,g3d_g(ngrid)%thsrc      & ! temp tendency
                     ,g3d_g(ngrid)%rtsrc      & ! rv tendency
                     ,g3d_g(ngrid)%clsrc      & ! cloud/ice  mass   mix ratio tendency
                     ,g3d_g(ngrid)%nlsrc      & ! cloud drop number mix ratio tendency
                     ,g3d_g(ngrid)%nisrc      & ! ice        number mix ratio tendency
                     ,g3d_g(ngrid)%usrc       & ! u tendency
                     ,g3d_g(ngrid)%vsrc       & ! v tendency
                     ,sub_mpqi    & 
                     ,sub_mpql    & 
                     ,sub_mpcf    & 
                     ,g3d_g(ngrid)%buoysrc    &
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
                     ,AA0,AA1,AA1_ADV,AA1_RADPBL,AA1_BL,AA2,AA3,AA1_CIN,TAU_BL,TAU_EC  &
                     ,VAR2d,VAR3d_aGF,VAR3d_bGF,VAR3d_cGF,VAR3d_dGF &
                     )
  !
  !-- outputs ....

   if( icumulus_gf(deep) == 1) then 
      g3d_g(ngrid)%aa0        (:,:) = AA0        (:,:)
      g3d_g(ngrid)%aa1        (:,:) = AA1        (:,:)
      g3d_g(ngrid)%aa1_bl     (:,:) = AA1_BL     (:,:)
      !g3d_g(ngrid)%aa1_radpbl (:,:) = AA1_RADPBL (:,:)
      g3d_g(ngrid)%aa1_adv    (:,:) = AA1_ADV    (:,:)
      !tmp
      g3d_g(ngrid)%aa1_radpbl (:,:) = AA3 (:,:)
   endif

   !-- saving the precip of each mode in the array g3d_ens_g(1,ngrid)%accapr(
   if( icumulus_gf(deep) == 1)  g3d_ens_g(1,ngrid)%apr(:,:)= cprr4d_tmp(:,:,deep)  
   if( icumulus_gf(mid)  == 1)  g3d_ens_g(3,ngrid)%apr(:,:)= cprr4d_tmp(:,:,mid)  
   if( icumulus_gf(shal) == 1)  g3d_ens_g(2,ngrid)%apr(:,:)= cprr4d_tmp(:,:,shal)  
    

   if( icumulus_gf(deep) == 1) then 
      ! updraft mass flux at cloud base
       g3d_g(ngrid)%xmb_deep   (:,:) = xmb4d_tmp       (:,:,deep)  
      ! downdraft mass flux at cloud base
       g3d_g(ngrid)%xmb_deep_dd(:,:) = xmb4d_tmp       (:,:,deep) * edt4d_tmp(:,:,deep)
       g3d_g(ngrid)%err_deep   (:,:) = float(ierr4d_tmp(:,:,deep))
   endif
   if( icumulus_gf(shal) == 1) then 
       g3d_g(ngrid)%xmb_shallow(:,:) = xmb4d_tmp       (:,:,shal)
   endif

   g3d_g(ngrid)%THSRC=EOSHIFT(g3d_g(ngrid)%THSRC, SHIFT=-1, BOUNDARY=g3d_g(ngrid)%THSRC(1,:,:), DIM=1)
   g3d_g(ngrid)%RTSRC=EOSHIFT(g3d_g(ngrid)%RTSRC, SHIFT=-1, BOUNDARY=g3d_g(ngrid)%RTSRC(1,:,:), DIM=1)
   g3d_g(ngrid)%CLSRC=EOSHIFT(g3d_g(ngrid)%CLSRC, SHIFT=-1, BOUNDARY=g3d_g(ngrid)%CLSRC(1,:,:), DIM=1)
   g3d_g(ngrid)%USRC =EOSHIFT(g3d_g(ngrid)%USRC,  SHIFT=-1, BOUNDARY=g3d_g(ngrid)%USRC (1,:,:), DIM=1)
   g3d_g(ngrid)%VSRC =EOSHIFT(g3d_g(ngrid)%VSRC,  SHIFT=-1, BOUNDARY=g3d_g(ngrid)%VSRC (1,:,:), DIM=1)


   if(liq_ice_number_conc > 0) then
      g3d_g(ngrid)%NLSRC=EOSHIFT(g3d_g(ngrid)%NLSRC, SHIFT=-1, BOUNDARY=g3d_g(ngrid)%NLSRC(1,:,:), DIM=1)
      g3d_g(ngrid)%NISRC=EOSHIFT(g3d_g(ngrid)%NISRC, SHIFT=-1, BOUNDARY=g3d_g(ngrid)%NISRC(1,:,:), DIM=1)
   endif

   !-- converting Dtemp/Dt to Dtheta/ Dt (temp = cp * theta/exner function) 
   g3d_g(ngrid)%THSRC = g3d_g(ngrid)%THSRC * cp / (basic_g(ngrid)%pp + basic_g(ngrid)%pi0)
 
   if( icumulus_gf(deep) == 1) then 
     
     do j=1,myp
       do i=1,mxp
           temp(:,i,j)= zup5d_tmp(i,j,:,deep) ! zup already includes the XMB
     enddo;enddo
     g3d_g(ngrid)%MUPDP=EOSHIFT(temp,  SHIFT=-1, BOUNDARY=temp (1,:,:), DIM=1)
     
     do j=1,myp
       do i=1,mxp
           temp(:,i,j)= zdn5d_tmp(i,j,:,deep)* edt4d_tmp(i,j,deep) ! zdo already includes the XMB
     enddo;enddo
     g3d_g(ngrid)%MDDDP=EOSHIFT(temp,  SHIFT=-1, BOUNDARY=temp (1,:,:), DIM=1)
   
     do j=1,myp
       do i=1,mxp
           temp(:,i,j)= up_massentr5d_tmp(i,j,:,deep) !  already includes the XMB
     enddo;enddo
     g3d_g(ngrid)%EUPDP=EOSHIFT(temp,  SHIFT=-1, BOUNDARY=temp (1,:,:), DIM=1)
     
     do j=1,myp
       do i=1,mxp
           temp(:,i,j)= up_massdetr5d_tmp(i,j,:,deep) !  already includes the XMB
     enddo;enddo
     g3d_g(ngrid)%DUPDP=EOSHIFT(temp,  SHIFT=-1, BOUNDARY=temp (1,:,:), DIM=1)

   endif

   if( icumulus_gf(shal) == 1) then 
     do j=1,myp
       do i=1,mxp
           temp(:,i,j)= zup5d_tmp(i,j,:,shal) ! zup already includes the XMB
     enddo;enddo
     g3d_g(ngrid)%MUPSH=EOSHIFT(temp,  SHIFT=-1, BOUNDARY=temp (1,:,:), DIM=1)
   endif
   
   if( icumulus_gf(mid) == 1) then 
     do j=1,myp
       do i=1,mxp
           temp(:,i,j)= zup5d_tmp(i,j,:,mid) ! zup already includes the XMB
     enddo;enddo
     g3d_g(ngrid)%MUPMD=EOSHIFT(temp,  SHIFT=-1, BOUNDARY=temp (1,:,:), DIM=1)
   endif

   !--- output for RRTM/CARMA and convective transport
   if((ilwrtyp>0 .or. iswrtyp>0) .or. chemistry >= 0) THEN
     do j=1,myp
      do i=1,mxp
        if(do_this_column(i,j)==0) cycle
	     do n=1,maxiens
             if( icumulus_gf(n) /= 1) cycle 
             xmb4d   (i,j,n,ngrid) = xmb4d_tmp   (i,j,n)
             ierr4d  (i,j,n,ngrid) = ierr4d_tmp  (i,j,n)
        enddo

	     do k=1,mzp
	         do n=1,maxiens
             if( icumulus_gf(n) /= 1) cycle 
             zup5d        (k,i,j,n,ngrid) =  zup5d_tmp        (i,j,k,n)
             clwup5d      (k,i,j,n,ngrid) =  clwup5d_tmp      (i,j,k,n)
             up_massdetr5d(k,i,j,n,ngrid) =  up_massdetr5d_tmp(i,j,k,n)
         enddo
        enddo
	     if( chemistry >= 0) THEN ! - for convective transport only
	       do n=1,maxiens
            if( icumulus_gf(n) /= 1) cycle 
	          jmin4d  (i,j,n,ngrid) = jmin4d_tmp  (i,j,n)
	          k224d   (i,j,n,ngrid) = k224d_tmp   (i,j,n)
	          kbcon4d (i,j,n,ngrid) = kbcon4d_tmp (i,j,n)
	          ktop4d  (i,j,n,ngrid) = ktop4d_tmp  (i,j,n)
	          kstabi4d(i,j,n,ngrid) = kstabi4d_tmp(i,j,n)
	          kstabm4d(i,j,n,ngrid) = kstabm4d_tmp(i,j,n)
	          edt4d   (i,j,n,ngrid) = edt4d_tmp   (i,j,n)
	          pwav4d  (i,j,n,ngrid) = pwav4d_tmp  (i,j,n)
          enddo

	       do k=1,mzp
	         do n=1,maxiens
             if( icumulus_gf(n) /= 1) cycle 
             zdn5d        (k,i,j,n,ngrid) =  zdn5d_tmp        (i,j,k,n)
             up_massentr5d(k,i,j,n,ngrid)	=  up_massentr5d_tmp(i,j,k,n)
             dd_massentr5d(k,i,j,n,ngrid)	=  dd_massentr5d_tmp(i,j,k,n)
             dd_massdetr5d(k,i,j,n,ngrid)	=  dd_massdetr5d_tmp(i,j,k,n)
             pcup5d       (k,i,j,n,ngrid)	=  pcup5d_tmp       (i,j,k,n)
	          prup5d       (k,i,j,n,ngrid)	=  prup5d_tmp       (i,j,k,n)
	          prdn5d       (k,i,j,n,ngrid)	=  prdn5d_tmp       (i,j,k,n)
	          tup5d        (k,i,j,n,ngrid)	=  tup5d_tmp        (i,j,k,n)
           enddo
          enddo
	     endif
     enddo;enddo
   endif

!--------------------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------------
 endif

!--- filling the output tendencies for the level k =1 and k=mzp
  g3d_g(ngrid)%THSRC(1  ,1:mxp,1:myp)= g3d_g(ngrid)%THSRC(2    ,1:mxp,1:myp)
  g3d_g(ngrid)%RTSRC(1  ,1:mxp,1:myp)= g3d_g(ngrid)%RTSRC(2    ,1:mxp,1:myp)
  g3d_g(ngrid)%CLSRC(1  ,1:mxp,1:myp)= g3d_g(ngrid)%CLSRC(2    ,1:mxp,1:myp)
  g3d_g(ngrid)%THSRC(mzp,1:mxp,1:myp)= g3d_g(ngrid)%THSRC(mzp-1,1:mxp,1:myp)
  g3d_g(ngrid)%RTSRC(mzp,1:mxp,1:myp)= g3d_g(ngrid)%RTSRC(mzp-1,1:mxp,1:myp)
  g3d_g(ngrid)%CLSRC(mzp,1:mxp,1:myp)= g3d_g(ngrid)%CLSRC(mzp-1,1:mxp,1:myp)
  if(imomentum==1 .and. nnqparm(ngrid) >= 4) then
   g3d_g(ngrid)%USRC(1  ,1:mxp,1:myp)= g3d_g(ngrid)%USRC (2    ,1:mxp,1:myp)
   g3d_g(ngrid)%VSRC(1  ,1:mxp,1:myp)= g3d_g(ngrid)%VSRC (2    ,1:mxp,1:myp)
   g3d_g(ngrid)%USRC(mzp,1:mxp,1:myp)= g3d_g(ngrid)%USRC (mzp-1,1:mxp,1:myp)
   g3d_g(ngrid)%VSRC(mzp,1:mxp,1:myp)= g3d_g(ngrid)%VSRC (mzp-1,1:mxp,1:myp)
  endif
  if(liq_ice_number_conc>0) then
   g3d_g(ngrid)%NLSRC(mzp,1:mxp,1:myp)= g3d_g(ngrid)%NLSRC(mzp-1,1:mxp,1:myp)
   g3d_g(ngrid)%NISRC(mzp,1:mxp,1:myp)= g3d_g(ngrid)%NISRC(mzp-1,1:mxp,1:myp)
  endif
!-------------------------------------------------------------
 endif! 002
!-------------------------------------------------------------
! stores precipitation rate for each closure, only for output/training
 if(nnqparm(ngrid) == 3 .and. training > 0) then
   do i=1,train_dim
     call update(mxp*myp, g3d_ens_g(i,ngrid)%accapr,g3d_ens_g(i,ngrid)%apr,dtlt)
   enddo
 endif
 !--- for output only 
 if(nnqparm(ngrid) == 8) then
   do i=1,train_dim
     call update(mxp*myp, g3d_ens_g(i,ngrid)%accapr,g3d_ens_g(i,ngrid)%apr,dtlt)
   enddo
 endif
!----------------------------------------------------------

 call update(mxp*myp, cuparm_g(ngrid)%aconpr   ,cuparm_g(ngrid)%conprr   ,dtlt)

 call accum(int(mxp*myp*mzp,i8), tend%tht, g3d_g(ngrid)%thsrc)
 call accum(int(mxp*myp*mzp,i8), tend%rtt, g3d_g(ngrid)%rtsrc)

 if(imomentum == 1 .and. nnqparm(ngrid) >= 4) then
  call accum(int(mxp*myp*mzp,i8), tend%ut, g3d_g(ngrid)%usrc)
  call accum(int(mxp*myp*mzp,i8), tend%vt, g3d_g(ngrid)%vsrc)
 endif

 if(do_cupar_mcphys_coupling) then
   call cupar2mcphysics(mzp,mxp,myp,ia,iz,ja,jz,ngrid,dtlt,&
                        g3d_g  (ngrid)%clsrc   ,&
			               basic_g(ngrid)%theta   ,&
			               basic_g(ngrid)%pp      ,&
			               basic_g(ngrid)%pi0     ,&
			               basic_g(ngrid)%dn0      )
 else
   !if there is not direct coupling, send cloud/ice source to rtotal tendency
   call accum(int(mxp*myp*mzp,i8), tend%rtt, g3d_g(ngrid)%clsrc)

 endif
!
!--------- Convective Transport based on mass flux scheme -
 if (CCATT == 1 .and. iruncon == 1 .and. (iinqparm==5 .or. iinqparm==6) ) then
     
     if(iinqparm==5 .and. iens .ne. 1 ) &
       stop 'conv transp with GF scheme version 2014 only for deep convection'

     !- this call convective transport for deep convection
     call trans_conv_mflx_GF(1,iinqparm)

     !- if shallow convection was solved by GF version 2015, call again
     !- the convective transport routine to include the transport
     !- by shallow convection scheme.
     if(iinqparm==6 .and. iinshcu == 3)  call trans_conv_mflx_GF(2,iinqparm)

 endif

 !- this calls convective transport for deep/mid/shallow convection
 if (CCATT == 1 .and. iruncon == 1 .and. (iinqparm==8) ) then
     do n = 1, maxiens
        if(n == shal .and. iinshcu == 0) cycle
        call trans_conv_mflx_GF(n,iinqparm)
     enddo
 endif


! [ML------------- Stilt - BRAMS coupling  ------------------
  if (imassflx == 1 ) then

       !-srf -  mass fluxes from deep convection
       if( iinqparm==5 .or. iinqparm==6 )                                 &
          call prep_convflx_to_stilt(mzp,mxp,myp,ia,iz,ja,jz              &
              ,mgmxp,mgmyp,mgmzp,maxiens,ngrid,ngrids_cp		  &
              ,ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d,kpbl4d 	  &
              ,kstabi4d,kstabm4d,xmb4d,edt4d				  &
              ,zcup5d,pcup5d,enup5d,endn5d,deup5d,dedn5d,zup5d,zdn5d	  &
              ,1)! = iens
       !-srf if shallow convection was solved by GF version 2015, call again
       !-    the convective transport routine to include the mass fluxes
       !-    from the shallow convection scheme.
       if(  iinqparm==6 .and. iinshcu == 3 )                              &
          call prep_convflx_to_stilt(mzp,mxp,myp,ia,iz,ja,jz              &
              ,mgmxp,mgmyp,mgmzp,maxiens,ngrid,ngrids_cp		  &
              ,ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d,kpbl4d 	  &
              ,kstabi4d,kstabm4d,xmb4d,edt4d				  &
              ,zcup5d,pcup5d,enup5d,endn5d,deup5d,dedn5d,zup5d,zdn5d	  &
              ,2)! = iens

   endif
! ------------- Stilt - BRAMS coupling  ------------------ ML]

 
!-- convection tracer update 
  if( nnqparm(ngrid) == 8 .and. convection_tracer == 1) then 
!      if (time>7200.) return 
!      if(naddsc > 0) then

         do j=1,myp
           do i=1,mxp
              tau_cp = 1.* tau_land_cp - (1.-1.)*tau_ocea_cp
              do k=1,mzp-1
                 kr=k+1

                 !- sink term (exp decay tau_cp hours)
                 snk_cnvtr = dtlt * abs( g3d_g(ngrid)%cnv_tr(k,i,j) )/tau_cp
   
                 !- source term
                 !- downdraft detrainment of buoyancy [ J/kg s^{-1}]
                 !- negative sign => source for updraft lifting
                 src_cnvtr = - dtlt * min(0.,g3d_g(ngrid)%buoysrc(k,i,j))
                 
                 !- 'continuity' equation = ADV + SRC - SINK
                 g3d_g(ngrid)%cnv_tr(k,i,j) = g3d_g(ngrid)%cnv_tr(k,i,j) + src_cnvtr - snk_cnvtr


 !               !-- send the convection_tracer out for ADV + DIFF
 !               snk_cnvtr = dtlt * abs( scalar_g(1,1)%sclp(kr,i,j) )/tau_cp
 !               scalar_g(1,1)%sclp(kr,i,j) = scalar_g(1,1)%sclp(kr,i,j) + src_cnvtr &
 !                                        - snk_cnvtr
 
         enddo;enddo;enddo
 !     else
 !        print*,"convection_tracer is on, but NADDSC is zero"
 !        stop 'wrong NADDSC for convection_tracer = 1'
 !     endif
   endif
!
end subroutine CUPARM_GRELL3_CATT
!
!-------------------------------------------------------------------------------------------------
!
subroutine init_weights(ng,n2,n3,nnqparm)
implicit none
integer, intent(in)::ng,n2,n3,nnqparm
integer :: it,i,j
real sumx,hweight
!- ordem dos pesos
!apr_gr=001
!apr_w =002
!apr_mc=003
!apr_st=004
!apr_as=005

if(training == 0) return

!-- training on closures
if(training == 1) then
   if(nnqparm==3) hweight = 0.2
   if(nnqparm==5) hweight = 0.25
   do j=1,n3
     do i=1,n2
      do it=1,train_dim

       g3d_ens_g(it,ng)%weight(i,j)=hweight
       !print*,'weights=', it,i,j, g3d_ens_g(it,ng)%weight(i,j)

       !if(it==apr_st) g3d_ens_g(it,ng)%weight(i,j)=0.175
       !if(it==apr_as) g3d_ens_g(it,ng)%weight(i,j)=0.25
       !if(it==apr_w ) g3d_ens_g(it,ng)%weight(i,j)=0.25
       !if(it==apr_mc) g3d_ens_g(it,ng)%weight(i,j)=0.25
       if(nnqparm==5) then

        !-special treatment over the ocean
        if(leaf_g(ng)%patch_area(i,j,1) .gt. 0.9 ) then ! water
           if(it==apr_as) g3d_ens_g(it,ng)%weight(i,j)=0.0
	        if(it==apr_st) g3d_ens_g(it,ng)%weight(i,j)=0.425
           if(it==apr_gr) g3d_ens_g(it,ng)%weight(i,j)=0.1667
           if(it==apr_w ) g3d_ens_g(it,ng)%weight(i,j)=0.1667
           if(it==apr_mc) g3d_ens_g(it,ng)%weight(i,j)=0.1667

        else ! land

           if(it==apr_as) g3d_ens_g(it,ng)%weight(i,j)=0.0
	        if(it==apr_st) g3d_ens_g(it,ng)%weight(i,j)=0.175
           if(it==apr_gr) g3d_ens_g(it,ng)%weight(i,j)=0.25
           if(it==apr_w ) g3d_ens_g(it,ng)%weight(i,j)=0.25
           if(it==apr_mc) g3d_ens_g(it,ng)%weight(i,j)=0.25
       endif
      endif
      !g3d_ens_g(it,ng)%weight(i,j)=float(i+j)*exp(-(float(it-2))**2)*float(i*j)

    enddo;enddo;enddo

!-- training on CAPS
elseif(training == 2) then

    do j=1,n3; do i=1,n2

       g3d_ens_g(apr_gr,ng)%weight(i,j)=0.3333
       g3d_ens_g(apr_w ,ng)%weight(i,j)=0.3333
       g3d_ens_g(apr_mc,ng)%weight(i,j)=0.3333
       g3d_ens_g(apr_st,ng)%weight(i,j)=0.0
       g3d_ens_g(apr_as,ng)%weight(i,j)=0.0

     enddo;enddo

endif



return! <<<<
if(training == 1) then
 !- normalize a 1
 do j=1,n3
    do i=1,n2
     sumx=0.
     do it=1,train_dim
       sumx=sumx+g3d_ens_g(it,ng)%weight(i,j)
     enddo
      do it=1,train_dim
      g3d_ens_g(it,ng)%weight(i,j) = g3d_ens_g(it,ng)%weight(i,j)/sumx
     enddo
 enddo;enddo
endif

end subroutine init_weights
!-------------------------------------------------------------
!-------------------------------------------------------------
SUBROUTINE conv_grell_spread3d_brams(m1,m2,m3,ia,iz,ja,jz,dt, &
               level,cugd_avedx,                              &
	       XLV,CP,G,r_v,p00,cpor,                         &
               conprr, theta,thetail,pp,pi0,                  &
	       rv,pt,rcp,rrp,rpp,rsp,rap,rgp,rhp,             &
	       RTHcuten,				      &
               RQVcuten,				      &
	       RQCcuten,				      &
               cugd_ttens,				      &
	       cugd_qvtens,				      &
               apr_gr,					      &
	       apr_w,					      &
	       apr_mc,					      &
	       apr_st,					      &
	       apr_as					      )

IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: m1,m2,m3,ia,iz,ja,jz,level,cugd_avedx
   REAL,         INTENT(IN   ) :: dt
   REAL,         INTENT(IN   ) :: XLV, R_v
   REAL,         INTENT(IN   ) :: CP,G, cpor, p00

   REAL, DIMENSION(m1,m2,m3),INTENT(IN   ) ::     &
 	  	 theta   ,&
          	 thetail ,&
          	 pp	 ,&
          	 pi0	 ,&
        	 pt	 ,&
        	 rv      ,rcp,rrp,rpp,rsp,rap,rgp,rhp


   REAL, DIMENSION(m2,m3),INTENT(INOUT) ::   &
               conprr,                       &
               apr_gr,			     &
	       apr_w ,			     &
	       apr_mc,			     &
	       apr_st,			     &
	       apr_as



   REAL, DIMENSION(m1,m2,m3),INTENT(INOUT) ::    &
                        RTHcuten,    &
                        RQVcuten,    &
                        RQCcuten,    &
                        cugd_ttens,  &
	                cugd_qvtens


  ! local var
   REAL  ::   exner,r_sol,r_liq,fxc,tempk,dfxcdt,outt
   INTEGER :: j,i,k,kk,jfs,jfe,ifs,ife,kts,kte,ii,jj
   INTEGER :: cugd_spread

   REAL, DIMENSION (m1,m2,m3) :: & ! orig (its-2:ite+2,kts:kte,jts-2:jte+2) ::     &
          RTHcuten_tmp, &  ! tmp RTHcuten
	  RQVcuten_tmp     ! tmp RQVcuten

   REAL, DIMENSION (m2,m3) :: & ! orig (its-2:ite+2,jts-2:jte+2) ::
          Qmem

   REAL   :: & ! orig (its-1:ite+1,jts-1:jte+1) ::
          smTT,smTQ

   REAL, DIMENSION (m1) :: & ! orig (kts:kte) ::
          conv_TRASHT,conv_TRASHQ

   REAL :: Qmem1,Qmem2,Qmemf,Thresh

  !-initial settings
  ! g3d_smoothh=0  ! 0 or 1: do horizontal smoothing
  ! g3d_smoothv=0  ! 0 or 1: do vertical smoothing
   cugd_spread=cugd_avedx/2 ! = 0/1 => no/do spreading

   RTHcuten_tmp  = 0.0
   RQVcuten_tmp  = 0.0
   Qmem       = 1.0
   smTT       = 0.0
   smTQ       = 0.0
   conv_TRASHT= 0.0
   conv_TRASHQ= 0.0
   jfs=ja
   jfe=jz
   ifs=ia
   ife=iz
   kts=2
   kte=m1-1 !check if this correct or should be kte=m1

   !if(g3d_smoothh ==1 .or. cugd_spread > 0) then
   !  jfs=1
   !  jfe=m3
   !  ifs=1
   !  ife=m2
   !endif

   !- store input tendencies
   ! *** jm note -- for smoothing this goes out one row/column beyond tile in i and j
   do j=1,m3
     do i=1,m2
         RTHcuten_tmp(:,i,j)=RTHcuten (:,i,j)
         RQVcuten_tmp(:,i,j)=RQVcuten (:,i,j)
     enddo
   enddo



! ---------------- spreading   section --------------
   do j=ja,jz
     do i=ia,iz
!
! for high res run, spread the subsidence
! this is tricky......only consider grid points where there was no rain,
! so cugd_tten and such are zero!
!
!      if do spreading
       if(cugd_spread > 0)then
         do k=kts,kte
	    do jj=j-1,j+1 ! only 3 neighboors
	      do ii=i-1,i+1 ! only 3 neighboors

               RTHcuten_tmp(k,i,j)=RTHcuten_tmp(k,i,j)     &
                                            +Qmem(ii,jj)*cugd_ttens(k,ii,jj)

               RQVcuten_tmp(k,i,j)=RQVcuten_tmp(k,i,j)     &
                                            +Qmem(ii,jj)*cugd_qvtens(k,ii,jj)
             enddo
           enddo
         enddo
!      end spreading
!
!      if not spreading
       elseif(cugd_spread == 0)then
         do k=kts,kte
           RTHcuten_tmp(k,i,j)=RTHcuten_tmp(k,i,j)+cugd_ttens (k,i,j)
           RQVcuten_tmp(k,i,j)=RQVcuten_tmp(k,i,j)+cugd_qvtens(k,i,j)
         enddo
       endif
!
     enddo  ! end i
   enddo  ! end j

! ----------------horizontal smoothing  section --------------

!- if not doing horizontal smoothing, get the final tendencies
  if(g3d_smoothh == 0)then
      do j=ja,jz
         do i=ia,iz
            do k=kts,kte
               RTHcuten(k,i,j)=RTHcuten_tmp(k,i,j)
               RQVcuten(k,i,j)=RQVcuten_tmp(k,i,j)
            enddo ! enf k
          enddo  ! end j
      enddo  ! end j

!- if doing horizontal smoothing ...
   else if(g3d_smoothh == 1)then
      do k=kts,kte
        do j=ja,jz
           do i=ia,iz

	    smTT = 0.0
            smTQ = 0.0
	    do jj=j-1,j+1 ! only 3 neighboors
	      do ii=i-1,i+1 ! only 3 neighboors
               smTT = smTT +RTHcuten_tmp(k,ii,jj)
	       smTQ = smTQ +RQVcuten_tmp(k,ii,jj)

              enddo  ! end ii
            enddo  ! end jj

            RTHcuten(k,i,j)=(3.*RTHcuten_tmp(k,i,j) + smTT)/12.
            RQVcuten(k,i,j)=(3.*RQVcuten_tmp(k,i,j) + smTQ)/12.

            enddo  ! end i
          enddo  ! end j
       enddo  ! end k

   endif  ! g3d_smoothh
  !
  ! - checking and limiting moistening/heating rates
  !
   do j=ja,jz
      do i=ia,iz
        !--- moistening section ------
	Qmemf  = 1.0
        Thresh = 1.e-20
        do k=kts,kte

	 if(RQVcuten(k,i,j) < 0.0) then
	    Qmem1 = rv(k,i,j)+RQVcuten(k,i,j)*dt
	    if(Qmem1 < Thresh)then
              Qmem1 = RQVcuten(k,i,j)
              Qmem2 = (Thresh-rv(k,i,j))/dt
              Qmemf = min(Qmemf,Qmem2/Qmem1)
              Qmemf = max(0.,Qmemf)
              Qmemf = min(1.,Qmemf)
            endif
          endif

         enddo  ! end k
         ! - limiting moistening
         do k=kts,kte
          RQVcuten   (k,i,j) = RQVcuten   (k,i,j)*Qmemf
          RQCcuten   (k,i,j) = RQCcuten   (k,i,j)*Qmemf
          RTHcuten   (k,i,j) = RTHcuten   (k,i,j)*Qmemf
	  cugd_ttens (k,i,j) = cugd_ttens (k,i,j)*Qmemf
          cugd_qvtens(k,i,j) = cugd_qvtens(k,i,j)*Qmemf

         enddo ! end k

	 !- limiting precip for consistency
	 conprr (i,j) = conprr (i,j)*Qmemf
         apr_gr (i,j) = apr_gr (i,j)*Qmemf
	 apr_w  (i,j) = apr_w  (i,j)*Qmemf
	 apr_mc (i,j) = apr_mc (i,j)*Qmemf
	 apr_st (i,j) = apr_st (i,j)*Qmemf
	 apr_as (i,j) = apr_as (i,j)*Qmemf
         ! no futuro inclua tambem o limting para o fluxo de massa
	 ! xmb(i,j)=xmb(i,j)*qmemf

	 !--- heating section ------
         Thresh=200. ! max heating/cooling rate allowed  K/day
!srf         Thresh=100. ! max heating/cooling rate allowed  K/day
         Qmemf=1.
         Qmem1=0.

	 do k=kts,kte
            Qmem1=abs(RTHcuten(k,i,j))*86400.

            if(Qmem1 > Thresh)then
              Qmem2 = Thresh/Qmem1
              Qmemf = min(Qmemf,Qmem2)
              Qmemf = max(0.,Qmemf)
            endif
         enddo

         ! - limiting heating/cooling
         do k=kts,kte
          RTHcuten   (k,i,j) = RTHcuten   (k,i,j)*Qmemf
          RQVcuten   (k,i,j) = RQVcuten   (k,i,j)*Qmemf
          RQCcuten   (k,i,j) = RQCcuten   (k,i,j)*Qmemf
	  cugd_ttens (k,i,j) = cugd_ttens (k,i,j)*Qmemf
          cugd_qvtens(k,i,j) = cugd_qvtens(k,i,j)*Qmemf
         enddo ! end k

	 !- limiting precip for consistency
	 conprr (i,j) = conprr (i,j)*Qmemf
         apr_gr (i,j) = apr_gr (i,j)*Qmemf
	 apr_w  (i,j) = apr_w  (i,j)*Qmemf
	 apr_mc (i,j) = apr_mc (i,j)*Qmemf
	 apr_st (i,j) = apr_st (i,j)*Qmemf
	 apr_as (i,j) = apr_as (i,j)*Qmemf
         ! no futuro inclua tambem o limting para o fluxo de massa
	 ! xmb(i,j)=xmb(i,j)*qmemf

     enddo  ! end i
   enddo  ! end j
  !
  ! ---  vertical smooth ------------
  !
   if (g3d_smoothv == 1)then

    do j=ja,jz
      do i=ia,iz

          do k=kts+2,kte-2
            conv_TRASHT(k)= .25*(RTHcuten(k-1,i,j)+2.*RTHcuten(k,i,j)+RTHcuten(k+1,i,j))
            conv_TRASHQ(k)= .25*(RQVcuten(k-1,i,j)+2.*RQVcuten(k,i,j)+RQVcuten(k+1,i,j))
          enddo
          do k=kts+2,kte-2
            RTHcuten(k,i,j)=conv_TRASHT(k)
            RQVcuten(k,i,j)=conv_TRASHQ(k)
          enddo
     enddo  ! end i
    enddo  ! end j

   endif

  ! Converte tend da temperatura (OUTT) em tend de theta (OUTTEM)
  ! cp*T=Pi*Theta => cp dT/dt = Theta*dPi/dt + Pi*dTheta/dt,

   !if(level <=2) then

    do j=ja,jz; do i=ia,iz; do k=kts,kte

	!if(RTHCUTEN (k,i,j) /= 0.0) then

           ! Converte tend da temperatura (OUTT) em tend de theta (OUTTEM)
           ! cp*T=Pi*Theta => cp dT/dt = Theta*dPi/dt + Pi*dTheta/dt,
           ! Exner's function = pp(k,i,j)+pi0(k,i,j)
            exner	   = pp(k,i,j) + pi0(k,i,j)
           ! tendencia do theta devida a conv profunda
           RTHcuten(k,i,J) = cp/exner * RTHCUTEN(k,i,J) !- theta(k,i,j)*pt(k,i,j)/exner
        !endif
	     !RQVcuten(k,i,J) = RQVCUTEN(k,i,J)+ RQCCUTEN(k,i,J)
        !RQCcuten(k,i,J) = 0.

    enddo; enddo; enddo


   !elseif(level > 2) then
   !
   ! do j=ja,jz; do i=ia,iz; do k=kts,kte
   !	    !
   !	 ! - tend na temperatura (para uso na convers�o do thetail
   !	 outt=RTHCUTEN (k,i,j)
   !	 ! Exner's function = pp(k,i,j)+pi0(k,i,j)
   !	 exner= pp(k,i,j) + pi0(k,i,j)
   !	 if(outt /= 0.0 ) then
   !	   !
   !	   ! converte tend da temperatura (outt) em tend de theta (outtem)
   !	   ! cp*T=Pi*Theta => cp dT/dt = Theta*dPi/dt + Pi*dTheta/dt,
   !
   !	   ! tendencia do theta  devida a conv profunda
   !	   RTHCUTEN (k,i,j) = cp/exner * RTHCUTEN(k,i,j) - theta(k,i,j)*pt(k,i,j)/exner
   !
   !	 endif
   !
   !	 ! tendencia do theta_il devida a conv profunda
   !	 r_liq= max(0.,rcp(k,i,j) + rrp(k,i,j))
   !
   !	 r_sol= max(0.,rsp(k,i,j)+rpp(k,i,j)+ &
   !		       rap(k,i,j)+rgp(k,i,j)+  &
   !		       rhp(k,i,j))
   !
   !	 tempk = theta(k,i,j)*(exner)/cp ! air temp (Kelvin)
   !
   !	 if(tempk.le.253) then
   !	   fxc =   (2.5e6*r_liq+2.83e6*r_sol)/(cp*amax1(tempk,253.))
   !
   !	   dfxcdt = 2.83e6*RQCCUTEN(k,i,J)/(cp*amax1(tempk,253.))
   !
   !	  RTHCUTEN (k,i,j) = (1./(1.+fxc))*( RTHCUTEN (k,i,j) - thetail(k,i,j)*dfxcdt )
   !
   !     else
   !
   !       fxc =   (2.5e6*r_liq+2.83e6*r_sol)/(cp*amax1(tempk,253.))
   !
   !!orig     dfxcdt = 2.5e6*OUTQC(I,K)*cuten(i)/(cp*amax1(tempk,253.)) - &
   !!orig         fxc/(cp*amax1(tempk,253.)) * cp * OUTT(I,K)
!  !
   !	  dfxcdt = 2.5e6*RQCCUTEN(k,i,J)/(cp*amax1(tempk,253.)) - &
   !          	   fxc/(cp*amax1(tempk,253.)) * cp * OUTT
   !
   !       RTHCUTEN (k,i,j) = (1./(1.+fxc))*( RTHCUTEN (k,i,j) - thetail(k,i,j)*dfxcdt )
   !
   !     endif
   !
   ! enddo; enddo; enddo
   !endif
  !- tendencies at boundaries
   RTHcuten(1,ia:iz,ja:jz)=RTHcuten(2,ia:iz,ja:jz)
   RQVcuten(1,ia:iz,ja:jz)=RQVcuten(2,ia:iz,ja:jz)
   RQCcuten(1,ia:iz,ja:jz)=RQCcuten(2,ia:iz,ja:jz)

   RTHcuten(m1,ia:iz,ja:jz)=RTHcuten(m1-1,ia:iz,ja:jz)
   RQVcuten(m1,ia:iz,ja:jz)=RQVcuten(m1-1,ia:iz,ja:jz)
   RQCcuten(m1,ia:iz,ja:jz)=RQCcuten(m1-1,ia:iz,ja:jz)

END SUBROUTINE conv_grell_spread3d_brams

!-------------------------------------------------------------

  subroutine StoreNamelistFileAtCup_grell3(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile

    g3d_spread = oneNamelistFile%g3d_spread
    g3d_smoothh = oneNamelistFile%g3d_smoothh
    g3d_smoothv = oneNamelistFile%g3d_smoothv

  end subroutine StoreNamelistFileAtCup_grell3

END MODULE CUPARM_GRELL3
!-----------------------------------------------------------------------------------

subroutine moveup(m1,A)
      implicit none
      integer, intent(in) :: m1
      real, dimension(m1) :: A,B

      real :: dummy
      integer :: k,kr

      B=A
      do k=1,m1-1
        kr=k+1
        A(kr) = B(k)
      enddo
      A(1) = A(2)
 end subroutine moveup
!-----------------------------------------------------------------------------------
 
 SUBROUTINE set_index_loops( ims,ime, jms,jme, kms,kme,    &
                             its,ite, jts,jte, kts,kte,    &
                             mxp,myp,mzp                   )

    IMPLICIT NONE
    INTEGER, INTENT(IN)         :: mxp,myp,mzp
    INTEGER, INTENT(INOUT)      :: ims,ime, jms,jme, kms,kme,&
                                   its,ite, jts,jte, kts,kte


    ims=1  ;ime=mxp ;jms=1  ;jme=myp ;kms=1 ;kme=mzp
    its=1  ;ite=mxp ;jts=1  ;jte=myp ;kts=1 ;kte=mzp

 END SUBROUTINE set_index_loops
!*************************************************************************************
subroutine check (m1,tht,ath,rtt,artt)
implicit none
integer, intent(in) :: m1
real, dimension(m1), intent(in) :: tht,ath,rtt,artt

integer k
do k=1,m1
print*,"check",k, tht(k),ath(k),rtt(k),artt(k)
enddo
print*,"mx1",maxval(tht),maxval(ath),minval(tht),minval(ath)
print*,"mx2",maxval(rtt),maxval(artt),minval(rtt),minval(artt)

end subroutine check
!------------------------------------------------------------
subroutine cupar2mcphysics(m1,m2,m3,ia,iz,ja,jz,ngrid,dtlt &
                          ,clsrc  ,theta,pp,pi0,dn0)
			  
  use micphys     ,only: level,mcphys_type
  use mem_micro   ,only: micro_g
  use mem_tend    ,only: tend
  use rconstants  ,only: cpi
  implicit none
  integer m1,m2,m3,ia,iz,ja,jz,k,i,j,ngrid
  real dtlt
  real, dimension(m1,m2,m3),intent(in) :: theta, pp, pi0,dn0
  real, dimension(m1,m2,m3),intent(in) :: clsrc! liquid/ice tendency from
                                    ! cumulus parameterization
   if(level < 2  .and. mcphys_type < 2 ) return

   if(level == 2 .and. mcphys_type < 2) then
       call mcphysics0(m1,m2,m3,ia,iz,ja,jz,dtlt &
	             ,clsrc                  &
                ,tend%rct(1)            &
		          ,tend%rtt(1)            )

   elseif(level == 3 .and. (mcphys_type >= 0))  then

       call mcphysics1(mcphys_type,m1,m2,m3,ia,iz,ja,jz,dtlt,cpi  &
            ,theta, pp, pi0,dn0 &
            ,clsrc          &! cumulus tendency
	         ,tend%rct(1)    &! cloud water mass mix ratio tendency 
            ,tend%rpt(1)    &! pristine mass mix ratio tendency 
            ,tend%rtt(1)    &! total water mass mix ratio tendency
            )
   
      if(mcphys_type == 2) &
       call mcphysics2(mcphys_type,m1,m2,m3,ia,iz,ja,jz,dtlt,cpi  &
            ,theta, pp, pi0,dn0 &
            ,clsrc          &! cumulus tendency
            ,tend%rct(1)    &! cloud water mass mix ratio tendency 
            ,tend%rpt(1)    &! pristine mass mix ratio tendency 
            ,tend%rtt(1)    &! total water mass mix ratio tendency
            ,tend%cpt(1)    &! pristine number conc tendency 
            )

      if(mcphys_type == 3) &
       call mcphysics3(mcphys_type,m1,m2,m3,ia,iz,ja,jz,dtlt,cpi  &
            ,theta, pp, pi0,dn0,micro_g(ngrid)%ccp &
            ,clsrc          &! cumulus tendency
            ,tend%rct(1)    &! cloud water mass mix ratio tendency 
            ,tend%rpt(1)    &! pristine mass mix ratio tendency 
            ,tend%rtt(1)    &! total water mass mix ratio tendency
            ,tend%cpt(1)    &! pristine number conc tendency 
            ,tend%cct(1)    &! cloud water  number conc tendency 
            )

    endif
return
end  subroutine cupar2mcphysics

!------------------------------------------------------------------------
subroutine mcphysics0(m1,m2,m3,ia,iz,ja,jz,dtlt,clsrc,rct,rtt)
implicit none
integer m1,m2,m3,ia,iz,ja,jz,k,i,j
real dtlt
real, dimension(m1,m2,m3),intent(in   ) :: clsrc
real, dimension(m1,m2,m3),intent(inout) :: rct,rtt

do j = ja,jz
   do i = ia,iz
     do k = 1,m1
       rct(k,i,j)=rct(k,i,j)+clsrc(k,i,j)
       rtt(k,i,j)=rtt(k,i,j)+clsrc(k,i,j)
     enddo
   enddo
enddo
end subroutine mcphysics0

!------------------------------------------------------------------------
subroutine mcphysics1(mcphys_type,m1,m2,m3,ia,iz,ja,jz,dtlt,cpi,theta,pp,pi0,dn0 &
                    ,clsrc,rct,rpt,rtt) 

use ConvPar_GF_GEOS5, only : fract_liq_f
implicit none
integer :: m1,m2,m3,ia,iz,ja,jz,k,i,j,mcphys_type
real dtlt,cpi
real, dimension(m1,m2,m3),intent(in)    :: clsrc
real, dimension(m1,m2,m3),intent(in)    :: theta,pp,pi0,dn0
real, dimension(m1,m2,m3),intent(inout) :: rct,rpt,rtt 
real ::tempk,tem1,tqice,tqliq,add_ncp,add_npp
real, parameter :: tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf)

do j = ja,jz
   do i = ia,iz
     do k = 1,m1
            tempk = theta(k,i,j)*(pp(k,i,j)+pi0(k,i,j))*cpi ! air temp (Kelvin)

            tem1 = fract_liq_f(tempk)
            
	         !- splitting cumulus tendency into water and ice tendencies
            rct(k,i,j) = rct(k,i,j)+clsrc(k,i,j)* tem1 ! cloud water

            rpt(k,i,j) = rpt(k,i,j)+clsrc(k,i,j)*(1.-tem1) ! pristine ice

            !- it must include also the ice/liq tendencies at rtt for
	         !- consistency, since rtt includes ice and liq mixing ratios
	         rtt(k,i,j) = rtt(k,i,j)+clsrc(k,i,j)

     enddo
   enddo
enddo

end subroutine mcphysics1
!------------------------------------------------------------------------
subroutine mcphysics2(mcphys_type,m1,m2,m3,ia,iz,ja,jz,dtlt,cpi,theta,pp,pi0,dn0 &
                    ,clsrc,rct,rpt,rtt,cpt)

use ConvPar_GF_GEOS5, only : make_IceNumber,fract_liq_f
implicit none
integer :: m1,m2,m3,ia,iz,ja,jz,k,i,j,mcphys_type
real dtlt,cpi
real, dimension(m1,m2,m3),intent(in)    :: clsrc
real, dimension(m1,m2,m3),intent(in)    :: theta,pp,pi0,dn0
real, dimension(m1,m2,m3),intent(inout) :: rct,rpt,rtt,cpt
real ::tempk,tem1,tqice,tqliq,add_ncp,add_npp
real, parameter :: tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf)

  do j = ja,jz
    do i = ia,iz
      do k = 1,m1
            tempk = theta(k,i,j)*(pp(k,i,j)+pi0(k,i,j))*cpi ! air temp (Kelvin)

            tem1 = fract_liq_f(tempk)

            !-- detrained pristine mass mixing ratio
	         tqice = (1.-tem1) * clsrc(k,i,j) * dn0(k,i,j)* dtlt

            !-- detrained ICN ice number concenration in the time "dtlt"
	         add_npp = max(0.0, make_IceNumber(tqice, tempk)/dn0(k,i,j))
            
	         !- update tendency 
            cpt(k,i,j) = cpt(k,i,j)+ add_npp/dtlt

      enddo
    enddo
  enddo
end subroutine mcphysics2

!------------------------------------------------------------------------
subroutine mcphysics3(mcphys_type,m1,m2,m3,ia,iz,ja,jz,dtlt,cpi,theta,pp,pi0,dn0 &
                    ,ccp,clsrc,rct,rpt,rtt,cpt,cct)

use ConvPar_GF_GEOS5, only : make_DropletNumber &
                            ,make_IceNumber     &
			    ,fract_liq_f
implicit none
integer :: m1,m2,m3,ia,iz,ja,jz,k,i,j,mcphys_type
real dtlt,cpi
real, dimension(m1,m2,m3),intent(in)    :: clsrc
real, dimension(m1,m2,m3),intent(in)    :: theta,pp,pi0,dn0,ccp
real, dimension(m1,m2,m3),intent(inout) :: rct,rpt,rtt,cpt,cct
real ::tempk,tem1,tqice,tqliq,add_ncp,add_npp
real, parameter :: tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf)

  do j = ja,jz
    do i = ia,iz
      do k = 1,m1
            tempk = theta(k,i,j)*(pp(k,i,j)+pi0(k,i,j))*cpi ! air temp (Kelvin)

            tem1 = fract_liq_f(tempk)

            !-- detrained pristine mass mixing ratio
            tqice = (1.-tem1) * clsrc(k,i,j) * dn0(k,i,j)* dtlt

            !-- detrained ICN ice number concenration in the time "dtlt"
            add_npp = max(0.0, make_IceNumber(tqice, tempk)/dn0(k,i,j))
            
            !- update tendency 
            cpt(k,i,j) = cpt(k,i,j)+ add_npp/dtlt

            !--cloud number concentration
	         tqliq =  tem1 * clsrc(k,i,j) * dn0(k,i,j) * dtlt

!--check if
            add_ncp = make_DropletNumber(tqliq, ccp(k,i,j))/dn0(k,i,j)
!or
!	         add_ncp = make_DropletNumber(tqliq, nwfa(k,i,j))/dn0(k,i,j)
           
	         !- update tendency 
	         cct(k,i,j) = cct(k,i,j)+ max(0.0, add_ncp/dtlt)

      enddo
    enddo
  enddo

end subroutine mcphysics3

!------------------------------------------------------------------------
subroutine prepare_lsf(nnqparm,nnshcu,iwork)

  use mem_grell   ,only: cuforc_g,cuforc_sh_g
  use mem_tend    ,only: tend
  !use mem_scratch ,only: scratch
  use mem_grid    ,only: time,ngrid,dtlt, dyncore_flag
  use mem_cuparm  ,only: confrq ,cuparm_g_sh
  use node_mod    ,only: mxp,myp,mzp ,ia,iz,ja,jz,mynum
  use mem_radiate, only: ilwrtyp, iswrtyp, radiate_g
  use mem_grid, only:ngrid, nzpmax, grid_g, dtlt, if_adap, jdim, time, &
       zt, zm, dzm, dzt, hw4,itopo
  use mem_basic, only: basic_g
  use mem_scratch, only : vctr1,vctr2

  implicit none
  include "i8.h"
  character(len=3) :: forcing
  integer,intent(IN) :: nnqparm,nnshcu,iwork
  !- scratchs (local arrays)
  real :: vt3da(mzp,mxp,myp)
  real :: vt3db(mzp,mxp,myp)
  real :: vt3dc(mzp,mxp,myp)
  real :: vt3dh(mzp,mxp,myp)
  real :: vt3dj(mzp,mxp,myp)
  real :: vt3dk(mzp,mxp,myp)
  real :: vt3di(mzp,mxp,myp)
  real :: vt3df(mzp,mxp,myp)
  real :: vt3dg(mzp,mxp,myp)
  real :: vt3de(mzp,mxp,myp)
  real :: vt3dd(mzp,mxp,myp)
 ! real :: vctr1(mzp)
 ! real :: vctr2(mzp)
  real :: scr1(mzp,mxp,myp)
  integer :: i,j,k
  !- parameter to define if include or not diffusion tendencies at forcing for deep convection
  logical,parameter :: forc_deep_pbl = .false.


 IF(mod(time,confrq).lt.dtlt .or. time .lt. dtlt+.01) then

       !-
       !  the forcing for shallow is only due to diffusion in PBL only (which is calculated in turb routines)
       !  the forcing for deep is due to radiation + 3dim advection
    if(iwork.eq.1) then

	 !----------- include radiation for theta
	      if(ilwrtyp + iswrtyp > 0  .and. nnqparm /= 8 ) then
	         cuforc_g(ngrid)%lsfth(1:mzp,1:mxp,1:myp)= radiate_g(ngrid)%fthrd(1:mzp,1:mxp,1:myp)
	      else
	         cuforc_g(ngrid)%lsfth(1:mzp,1:mxp,1:myp)= 0.
         endif
    
         !-reset lsf for water vapor
         cuforc_g(ngrid)%lsfrt(1:mzp,1:mxp,1:myp)= 0.

	      !----------- include advection for theta and rv (or should be rtp?)
         vt3dd=0.0
         vt3de=0.0
         vt3df=0.0
         vt3dg=0.0
         vt3dh=0.0
         vt3di=0.0
         vt3dj=0.0
         vt3dk=0.0
         vctr1=0.0
         vctr2=0.0
	      if(dyncore_flag == 0) then
          do j = 1,myp
            do i = 1,mxp
              do k = 1,mzp
                vt3da(k,i,j) = (basic_g(ngrid)%up(k,i,j)+ basic_g(ngrid)%uc(k,i,j))*dtlt*0.5
                vt3db(k,i,j) = (basic_g(ngrid)%vp(k,i,j)+ basic_g(ngrid)%vc(k,i,j))*dtlt*0.5
                vt3dc(k,i,j) = (basic_g(ngrid)%wp(k,i,j)+ basic_g(ngrid)%wc(k,i,j))*dtlt*0.5
              end do
            end do
          end do
         else
	        do j = 1,myp
            do i = 1,mxp
              do k = 1,mzp
                vt3da(k,i,j) = basic_g(ngrid)%uc(k,i,j)*dtlt
                vt3db(k,i,j) = basic_g(ngrid)%vc(k,i,j)*dtlt
                vt3dc(k,i,j) = basic_g(ngrid)%wc(k,i,j)*dtlt
              end do
            end do
          end do
         endif
         call fa_preptc(mzp,mxp,myp            &
	     ,vt3da	     ,vt3db	       &
	     ,vt3dc	     ,vt3dd	       &
	     ,vt3de	     ,vt3df	       &
	     ,vt3dh	     ,vt3di	       &
	     ,vt3dj	     ,vt3dk	       &
	     ,mynum			       )

	     if(dyncore_flag == 0) then
         !---- thp
	      scr1(1:mzp,1:mxp,1:myp) = basic_g(ngrid)%thp(1:mzp,1:mxp,1:myp)

         ! output: scr1,vt3dg
         call fa_xc(mzp,mxp,myp,ia,iz,1,myp,basic_g(ngrid)%thp,scr1,vt3da,vt3dd,vt3dg,vt3dh,vt3di,mynum)

         ! input: scalarp, scr1,vt3db,vt3de,vt3dj,vt3di
         ! output: scr1,vt3dg
         if (jdim .eq. 1)  &
              call fa_yc(mzp,mxp,myp,ia,iz,ja,jz,basic_g(ngrid)%thp,scr1,vt3db,vt3de,vt3dg,vt3dj,vt3di,jdim,mynum)

         ! input: scalarp, scr1,vt3dc,vt3df,vt3dk, vctr1,vctr2
         ! output: scr1,vt3dg
         call fa_zc(mzp,mxp,myp,ia,iz,ja,jz,basic_g(ngrid)%thp,scr1,vt3dc,vt3df,vt3dg,vt3dk,vctr1,vctr2,mynum)

         ! input:  thetap , lsfth,scr1, dtlt
         ! output: lsfth
         call advtndc(mzp,mxp,myp,ia,iz,ja,jz,basic_g(ngrid)%thp,scr1,cuforc_g(ngrid)%lsfth,dtlt,mynum)
         !
	     else
         !---- thc
         scr1(1:mzp,1:mxp,1:myp) = basic_g(ngrid)%thc(1:mzp,1:mxp,1:myp)

         ! output: scr1,vt3dg
         call fa_xc(mzp,mxp,myp,ia,iz,1,myp,basic_g(ngrid)%thc,scr1,vt3da,vt3dd,vt3dg,vt3dh,vt3di,mynum)

         ! input: scalarp, scr1,vt3db,vt3de,vt3dj,vt3di
         ! output: scr1,vt3dg
         if (jdim .eq. 1)  &
              call fa_yc(mzp,mxp,myp,ia,iz,ja,jz,basic_g(ngrid)%thc,scr1,vt3db,vt3de,vt3dg,vt3dj,vt3di,jdim,mynum)

         ! input: scalarp, scr1,vt3dc,vt3df,vt3dk, vctr1,vctr2
         ! output: scr1,vt3dg
         call fa_zc(mzp,mxp,myp,ia,iz,ja,jz,basic_g(ngrid)%thc,scr1,vt3dc,vt3df,vt3dg,vt3dk,vctr1,vctr2,mynum)

         ! input:  thetac , lsfth,scr1, dtlt
         ! output: lsfth
         call advtndc(mzp,mxp,myp,ia,iz,ja,jz,basic_g(ngrid)%thc,scr1,cuforc_g(ngrid)%lsfth,dtlt,mynum)
        endif

        !---- water vapor
	     scr1(1:mzp,1:mxp,1:myp) = basic_g(ngrid)%rv(1:mzp,1:mxp,1:myp)

        ! output: scr1,vt3dg
        call fa_xc(mzp,mxp,myp,ia,iz,1,myp,basic_g(ngrid)%rv,scr1,vt3da,vt3dd,vt3dg,vt3dh,vt3di,mynum)

        ! input: scalarp, scr1,vt3db,vt3de,vt3dj,vt3di
        ! output: scr1,vt3dg
        if (jdim .eq. 1)  &
              call fa_yc(mzp,mxp,myp,ia,iz,ja,jz,basic_g(ngrid)%rv,scr1,vt3db,vt3de,vt3dg,vt3dj,vt3di,jdim,mynum)

         ! input: scalarp, scr1,vt3dc,vt3df,vt3dk, vctr1,vctr2
         ! output: scr1,vt3dg
         call fa_zc(mzp,mxp,myp,ia,iz,ja,jz,basic_g(ngrid)%rv,scr1,vt3dc,vt3df,vt3dg,vt3dk,vctr1,vctr2,mynum)

         ! input: basic(ngrid)%rv, scalart,scr1, dtlt
         ! output: lsfrt = rad + adv
         call advtndc(mzp,mxp,myp,ia,iz,ja,jz,basic_g(ngrid)%rv,scr1,cuforc_g(ngrid)%lsfrt,dtlt,mynum)

         !- here the forcings contain rad+adv for temp and adv for water vapor
         !-end of inclusion of the advection forcings

         !- flag to include diffusion in pbl (only vertical) in deep convection forcing.
         !- the pbl forcing (cuforc_sh_g(ngrid)%lsfth and %lsfrt) has been calculated
         !- in the turbulence routines.
         if(forc_deep_pbl) then
          ! calcula o forcing para conveccao profunda = rad + pbl turb + adv
          ! for deep convection    LSF =  radiation + pbl_turb + advection
          cuforc_g(ngrid)%lsfth(:,:,:) = cuforc_g(ngrid)%lsfth(:,:,:)+cuforc_sh_g(ngrid)%lsfth(:,:,:)
          cuforc_g(ngrid)%lsfrt(:,:,:) = cuforc_g(ngrid)%lsfrt(:,:,:)+cuforc_sh_g(ngrid)%lsfrt(:,:,:)
         endif

     endif
     if(iwork.eq.2) then
       call accum(int(mxp*myp*mzp,i8), cuforc_g(ngrid)%lsfth(1,1,1), cuparm_g_sh(ngrid)%thsrc(1,1,1))
     	 call accum(int(mxp*myp*mzp,i8), cuforc_g(ngrid)%lsfrt(1,1,1), cuparm_g_sh(ngrid)%rtsrc(1,1,1))
     endif

ENDIF
end subroutine prepare_lsf

!-------------------------------------------------------------------
SUBROUTINE get_zi_gf2018(m1,tkmin,tkeg,z,rtgt,ztop,kzi)

  implicit none
  integer,intent(in):: m1
  integer :: kzimax,ktke_max,i,k
  real tkmin,tke_tmp
  real,intent(in), dimension(m1) :: tkeg,z
  real,intent(in)    :: ztop,rtgt
  integer,intent(out) :: kzi
  real, parameter :: rcpmin=1.e-5 , pblhmax=3000.

  kzi      = 1
  ktke_max = 1
  kzimax   = m1-1
  !---  max level for kzi
    do k=1,m1
    if(z(k).ge. pblhmax+ztop) then
       kzimax = min(k,m1-1)
       !if(j==8 .and. i==10) print*,"1",z(i,k), pblhmax,ztop(i),kzimax
       exit
    endif
  enddo
  !---
  do k=ktke_max,kzimax
    if(tkeg(k) .gt. 1.1*tkmin )  then
      kzi = k
      cycle
    else
       kzi = max(1,k-1)
       exit
    endif
  enddo
  kzi = max(1	  ,kzi)
  kzi = min(kzimax,kzi)
  !print*,"2",kzi(i),i;call flush(6)
  !pbl(i) = max( z(i,kzi(i))-ztop(i), z(i,1)-ztop(i) )
 END SUBROUTINE get_zi_gf2018
!-------------------------------------------------------------------
SUBROUTINE get_mean_tke(m1,tkmin,tke1d,rtgt,kzi,dn01d,tke_pbl)
  
  use mem_grid, only:   dzt ! intent(IN)
  implicit none
  integer,intent(in):: m1, kzi
  real, intent(in):: tkmin,rtgt
  real, intent(in), dimension(m1) :: tke1d,dn01d
  real, intent(out) :: tke_pbl
  integer :: k
  real :: dzpho,total_dz, shmf
  real :: k1 = 1.2 , cloud_area = 0.15

  tke_pbl  = 0. 
  total_dz = 0.

  do k = 2, kzi + 1
      dzpho    = rtgt/dzt(k) * dn01d(k)
      tke_pbl  = tke_pbl  + tke1d(k) * dzpho
      total_dz = total_dz + dzpho
  enddo
  tke_pbl = tke_pbl / (1.e-6 + total_dz) 
  tke_pbl = max(tkmin, tke_pbl)
 
  !-- just for checking 
    !-- potential closure for the mass flux shallow convection
    !shmf = cloud_area * dn01d(kzi) * k1 * sqrt(tke_pbl)
    !tke_pbl = shmf
  !--
   
 END SUBROUTINE get_mean_tke
!-------------------------------------------------------------------
!------------------------------------------------------------------------
subroutine adv_convection_tracer()
 
  use cuparm_grell3, only:  g3d_g
  use mem_tend    ,only: tend
  use mem_grid    ,only: time,ngrid,dtlt, dyncore_flag,dzt,ztn
  use mem_cuparm  ,only: confrq ,cuparm_g_sh
  use node_mod    ,only: mxp,myp,mzp ,ia,iz,ja,jz,mynum, nmachs
  use mem_grid    ,only: ngrid, nzpmax, grid_g, dtlt, if_adap, jdim, time, &
                         zt, zm, dzm, dzt, hw4,itopo
  use mem_basic   ,only: basic_g
  use mem_scratch ,only: vctr1,vctr2
  use mem_turb    ,only: turb_vars, turb_g , xkhkm
  use mem_scratch1_grell, only: jmin4d
  USE ConvPar_GF_GEOS5, only: add_coldpool_prop,mx_buoy2,use_gustiness, add_coldpool_diff
  use rconstants
  
  USE extras, ONLY:   &
         na_extra3d,  &    
         extra3d,     &         
         na_extra2d,  &
         extra2d

  use ModComm, only: commHaloAcou

  implicit none
  include "i8.h"
  real, parameter :: xlv = 2.5e6, Kfr = 0.9, epsx = 1.e2
  !- scratchs (local arrays)
  real :: vt3da(mzp,mxp,myp)
  real :: vt3db(mzp,mxp,myp)
  real :: vt3dc(mzp,mxp,myp)
  real :: vt3dh(mzp,mxp,myp)
  real :: vt3dj(mzp,mxp,myp)
  real :: vt3dk(mzp,mxp,myp)
  real :: vt3di(mzp,mxp,myp)
  real :: vt3df(mzp,mxp,myp)
  real :: vt3dg(mzp,mxp,myp)
  real :: vt3de(mzp,mxp,myp)
  real :: vt3dd(mzp,mxp,myp)
  real :: scr1(mzp,mxp,myp)
  real, dimension (mzp,mxp,myp) :: scalart
  real, dimension (mxp,myp)     :: sfcflx

  integer :: i,j,k, k_lfs,kr,n, k_amb_wind
  !- parameter to define if include or not vertical advection
  logical,parameter :: vert_adv = .true. ! must be true even for 2-D advection

 !real, parameter :: p2=300., p1=850., p3=900., p4=750.
  real, parameter :: p2=400., p1=750., p3=950., p4=750., p5=900.,p6=500.,p7=600.&
                    ,p1allj=980.,p2allj=800.
  real :: press, dzpho,total_dz, u, v, w, ullj,vllj,wllj,total_dz_llj
  real :: theta2temp,temp , rvap   ,zmsl , H_env ,dz,     MCL_speed, aux, dtlt_local 
  character(len=255) :: wind_type = 'lower_troposphere' !  'Corfidi'
  real, parameter  :: factor=0.6,  tuning=1.0
  real, parameter  :: alpha = 45 * 3.1416/180. ! 45 graus slope of cold pool head
  integer, parameter  :: nsubsteps = 1
  logical :: singleProcRun

  singleProcRun = nmachs == 1 !Avalia se está rodando com mais de 1 proc

  !do i = 1,mxp;  do j = 1,myp
  !if(jmin4d(i,j,1,1) > 0) print*,"jmin",i,j,jmin4d(i,j,1,1)
  !enddo;enddo
  !k_lfs=30
         
         dtlt_local = dtlt/float(nsubsteps)
         !-
         !----------- includes advection for the convection tracer 
         vt3dd=0.0
         vt3de=0.0
         vt3df=0.0
         vt3dg=0.0
         vt3dh=0.0
         vt3di=0.0
         vt3dj=0.0
         vt3dk=0.0
         vctr1=0.0
         vctr2=0.0

      if(trim(wind_type) == 'real_wind') then  !=============

         if(dyncore_flag == 0) then
          do j = 1,myp
            do i = 1,mxp
              do k = 1,mzp             
                k_lfs = k
                vt3da(k,i,j) = (basic_g(ngrid)%up(k_lfs,i,j)+ basic_g(ngrid)%uc(k_lfs,i,j))*dtlt_local*0.5
                vt3db(k,i,j) = (basic_g(ngrid)%vp(k_lfs,i,j)+ basic_g(ngrid)%vc(k_lfs,i,j))*dtlt_local*0.5
                vt3dc(k,i,j) = (basic_g(ngrid)%wp(k_lfs,i,j)+ basic_g(ngrid)%wc(k_lfs,i,j))*dtlt_local*0.5
              end do
            end do
          end do
        
         else

          do j = 1,myp
            do i = 1,mxp
              do k = 1,mzp
                 k_lfs = k
                vt3da(k,i,j) = basic_g(ngrid)%uc(k_lfs,i,j)*dtlt_local
                vt3db(k,i,j) = basic_g(ngrid)%vc(k_lfs,i,j)*dtlt_local
                vt3dc(k,i,j) = basic_g(ngrid)%wc(k_lfs,i,j)*dtlt_local
              end do
            end do
          end do
         
         endif

      elseif(trim(wind_type) == 'Corfidi') then              

          do j = 1,myp
            do i = 1,mxp
              
              !=============  Define the mean cloud motion wind
              total_dz = 0. ; u = 0. ; v = 0. ; w = 0.
              do k = 1,mzp
                press = 1.e-2*((basic_g(ngrid)%pp(k,i,j)+basic_g(ngrid)%pi0(k,i,j))/cp)**cpor*p00 !hPa
                if(press .le. p1 .and. press .ge. p2 ) then 

                  dzpho  = grid_g(ngrid)%rtgt(i,j)/dzt(k) * basic_g(ngrid)%dn0(k,i,j)
                  u = u + basic_g(ngrid)%uc(k,i,j)*dzpho
                  v = v + basic_g(ngrid)%vc(k,i,j)*dzpho
                  w = w + basic_g(ngrid)%wc(k,i,j)*dzpho
                  total_dz = total_dz + dzpho                
                endif 
              end do

              !=============  low level jet 
              total_dz_llj= 0.;  ullj= 0. ; vllj= 0. ;  wllj = 0.  
              do k = 1,mzp
                press = 1.e-2*((basic_g(ngrid)%pp(k,i,j)+basic_g(ngrid)%pi0(k,i,j))/cp)**cpor*p00 !hPa
                if(press .gt. p1 .and. press .le. p3 ) then 

                  dzpho  = grid_g(ngrid)%rtgt(i,j)/dzt(k) * basic_g(ngrid)%dn0(k,i,j)
                  ullj = ullj + basic_g(ngrid)%uc(k,i,j)*dzpho
                  vllj = vllj + basic_g(ngrid)%vc(k,i,j)*dzpho
                  wllj = wllj + basic_g(ngrid)%wc(k,i,j)*dzpho
                  total_dz_llj = total_dz_llj + dzpho                
                endif 
              end do
     
              !============= combines MCMW - LLJ (Corfidi 2003) => MCS propagation
              vt3da(:,i,j) = u/(total_dz+1.e-8)*dtlt_local +  ullj/(total_dz_llj+1.e-8)*dtlt_local
              vt3db(:,i,j) = v/(total_dz+1.e-8)*dtlt_local +  vllj/(total_dz_llj+1.e-8)*dtlt_local
              vt3dc(:,i,j) = w/(total_dz+1.e-8)*dtlt_local +  wllj/(total_dz_llj+1.e-8)*dtlt_local
              !print*,'U2=',i,j, maxval(vt3da(:,i,j))/dtlt_local,maxval(basic_g(ngrid)%uc(:,i,j))
            end do
          end do
      !
      else
      
         do j = 1,myp
            do i = 1,mxp
              
              !=============  Define the mean cloud motion wind by the lowest part of the troposphere
              total_dz = 0. ; u = 0. ; v = 0. !; w = 0.
              do k = 1,mzp
                press = 1.e-2*((basic_g(ngrid)%pp(k,i,j)+basic_g(ngrid)%pi0(k,i,j))/cp)**cpor*p00 !hPa
               
               !-artigo
               if(press .ge. p7 .and. press .le. p5 ) then !p5=900.,p7=600.
               !-Amz llj
               !if(press .ge. p2allj .and. press .le. p1allj ) then !p1allj=980.,p2allj=800.

                  dzpho    = grid_g(ngrid)%rtgt(i,j)/dzt(k) * basic_g(ngrid)%dn0(k,i,j)
                  u        = u + basic_g(ngrid)%uc(k,i,j)*dzpho
                  v        = v + basic_g(ngrid)%vc(k,i,j)*dzpho
                 !w        = w + basic_g(ngrid)%wc(k,i,j)*dzpho
                  total_dz = total_dz + dzpho
                  
                endif 
              end do
              !============= 
              vt3da(:,i,j) = u/(total_dz+1.e-8)*dtlt_local 
              vt3db(:,i,j) = v/(total_dz+1.e-8)*dtlt_local 
              vt3dc(:,i,j) = 0. !w/(total_dz+1.e-8)*dtlt_local 
              
            end do
          end do

      endif     !=============

      !-- save the mean cloud layer velocity
      g3d_g(ngrid)%umcl(:,:) = vt3da(2,:,:)/dtlt_local  !mean cloud layer U-wind
      g3d_g(ngrid)%vmcl(:,:) = vt3db(2,:,:)/dtlt_local  !mean cloud layer V-wind
      
      !--- add cold pool velocity
      if( add_coldpool_prop .ge. 0 ) then
         do j = 1,myp
             do i = 1,mxp

                g3d_g(ngrid)%wlpool(i,j) = 0.
                g3d_g(ngrid)%vcpool(i,j) = 0.
                total_dz = 0.
                
               loopK: do k = 2,mzp
                  press = 1.e-2*((basic_g(ngrid)%pp(k,i,j)+basic_g(ngrid)%pi0(k,i,j))/cp)**cpor*p00 !hPa
                     
                  !if( g3d_g(ngrid)%cnv_tr(k,i,j) < epsx .or. press    .lt. p7    ) cycle 
                  if( g3d_g(ngrid)%cnv_tr(k,i,j) < epsx .or. total_dz .gt. 2000. ) exit loopK 
                 
                  theta2temp = ( basic_g(ngrid)%pp   (k,i,j) + basic_g(ngrid)%pi0(k,i,j))/cp   !K
                  temp       =   basic_g(ngrid)%theta(k,i,j) * theta2temp  !K
                  rvap       =   basic_g(ngrid)%rv   (k,i,j)               !kg/kg
                  zmsl       =   ztn(k,ngrid)*grid_g(ngrid)%rtgt(i,j) + grid_g(ngrid)%topt(i,j)  !m
                  
                  H_env      = cp * temp + xlv * rvap + g* zmsl !J/kg
                  dz         = grid_g(ngrid)%rtgt(i,j)/dzt(k)  !m
                  total_dz   = total_dz + dz
                
                  !-- See eq 8.5 page 346 Storm and Cloud Dynamics, Cotton et al. 2nd edition.
                  g3d_g(ngrid)%vcpool(i,j) = g3d_g(ngrid)%vcpool(i,j) +  &
                                             (g*dz*min(mx_buoy2,g3d_g(ngrid)%cnv_tr(k,i,j))/H_env)
                  
                  !-- Max vertical velocity at the leading of the cold pool
                  !-- see "Estimating the Maximum Vertical Velocity at the Leading Edge 
                  !--      of a Density Current, Reif et al 2020"
                  g3d_g(ngrid)%wlpool(i,j) = g3d_g(ngrid)%wlpool(i,j) +                      &
                                             (g*dz*min(mx_buoy2,g3d_g(ngrid)%cnv_tr(k,i,j)) * & 
                                                   sin(alpha)**2 / H_env )

                  kr = k

                end do loopK
                !-- see page 346 of the same book.
                g3d_g(ngrid)%vcpool(i,j) = min(20., Kfr *sqrt(g3d_g(ngrid)%vcpool(i,j)))
                
                g3d_g(ngrid)%wlpool(i,j) = min(10., Kfr *sqrt(g3d_g(ngrid)%wlpool(i,j)))

                !if(g3d_g(ngrid)%vcpool(i,j) > 1.) &
                !print*,"cp=", i,j,g3d_g(ngrid)%vcpool(i,j) ,maxval( g3d_g(ngrid)%cnv_tr(2:k,i,j) )/1000.,kr 

             end do
         end do
         
         !-- replace the MCL wind by the cold pool velocity, keeping the MCL wind direction
         if( add_coldpool_prop == 1 ) then
          do j = 1,myp
             do i = 1,mxp
             
             MCL_speed= sqrt( (g3d_g(ngrid)%umcl(i,j))**2 + (g3d_g(ngrid)%vmcl(i,j))**2 )
             aux =  1./(MCL_speed+1.e-6)
             
             g3d_g(ngrid)%umcl(i,j) = aux * g3d_g(ngrid)%umcl(i,j) *  g3d_g(ngrid)%vcpool(i,j)
             g3d_g(ngrid)%vmcl(i,j) = aux * g3d_g(ngrid)%vmcl(i,j) *  g3d_g(ngrid)%vcpool(i,j)

             vt3da(:,i,j)=  g3d_g(ngrid)%umcl(i,j) * dtlt_local  !mean cloud layer U-wind is the cold pool propagation
             vt3db(:,i,j)=  g3d_g(ngrid)%vmcl(i,j) * dtlt_local  !mean cloud layer V-wind is the cold pool propagation
             end do
          end do       
         endif 

         !-- add cold pool to the MCL wind, keeping the original MCL wind direction
         !-- 'factor' = 0.62 is suggested by Simpson and Britter (1980) to add the cold pool velocity 
         !-- the to ambient wind (see page 346 Storm and Cloud Dynamics, Cotton et al. 2nd edition)
         if( add_coldpool_prop == 2 ) then
          do j = 1,myp
            do i = 1,mxp
           
             MCL_speed= sqrt( (g3d_g(ngrid)%umcl(i,j))**2 + (g3d_g(ngrid)%vmcl(i,j))**2 )
             aux =  (factor * MCL_speed + g3d_g(ngrid)%vcpool(i,j))/(MCL_speed+1.e-6)
             g3d_g(ngrid)%umcl(i,j) = aux * g3d_g(ngrid)%umcl(i,j)
             g3d_g(ngrid)%vmcl(i,j) = aux * g3d_g(ngrid)%vmcl(i,j)

             vt3da(:,i,j)=  g3d_g(ngrid)%umcl(i,j) * dtlt_local  !mean cloud layer U-wind with cold pool propagation
             vt3db(:,i,j)=  g3d_g(ngrid)%vmcl(i,j) * dtlt_local  !mean cloud layer V-wind with cold pool propagation
           
            end do
          end do       
         endif 


         !-- add the influence of ambient flow to the MCL wind
         !-- 'factor' = 0.62 is suggested by Simpson and Britter (1980) to add the cold pool velocity 
         !-- the to ambient wind (see page 346 Storm and Cloud Dynamics, Cotton et al. 2nd edition)
         if( add_coldpool_prop == 3 ) then
          k_amb_wind = 2
          do j = 1,myp
            do i = 1,mxp
           
             MCL_speed= sqrt( (g3d_g(ngrid)%umcl(i,j))**2 + (g3d_g(ngrid)%vmcl(i,j))**2 )
             aux =  (MCL_speed + g3d_g(ngrid)%vcpool(i,j))/(MCL_speed+1.e-6)
             g3d_g(ngrid)%umcl(i,j) = aux * g3d_g(ngrid)%umcl(i,j)
             g3d_g(ngrid)%vmcl(i,j) = aux * g3d_g(ngrid)%vmcl(i,j)

             !-- including the influence of ambient flow to the MCL wind

             g3d_g(ngrid)%umcl(i,j) = g3d_g(ngrid)%umcl(i,j) + factor*basic_g(ngrid)%uc(k_amb_wind,i,j)
             g3d_g(ngrid)%vmcl(i,j) = g3d_g(ngrid)%vmcl(i,j) + factor*basic_g(ngrid)%vc(k_amb_wind,i,j)

             !-- 
             vt3da(:,i,j)=  g3d_g(ngrid)%umcl(i,j) * dtlt_local  !mean cloud layer U-wind with cold pool propagation
             vt3db(:,i,j)=  g3d_g(ngrid)%vmcl(i,j) * dtlt_local  !mean cloud layer V-wind with cold pool propagation
           
            end do
          end do       
         endif 


      endif      

    call fa_preptc(mzp,mxp,myp  &
        ,vt3da      ,vt3db        &
        ,vt3dc      ,vt3dd        &
        ,vt3de      ,vt3df        &
        ,vt3dh      ,vt3di        &
        ,vt3dj      ,vt3dk        &
        ,mynum                    )

    !---- 2-d advection plus 3 or 2-d diffusion
    do n=1,nsubsteps

         !---- tracer to be advected 
         scr1(1:mzp,1:mxp,1:myp) = g3d_g(ngrid)%cnv_tr(1:mzp,1:mxp,1:myp)

         ! output: scr1,vt3dg
         call fa_xc(mzp,mxp,myp,ia,iz,1,myp,g3d_g(ngrid)%cnv_tr,scr1,vt3da,vt3dd,vt3dg,vt3dh,vt3di,mynum)

         ! input: scalarp, scr1,vt3db,vt3de,vt3dj,vt3di
         ! output: scr1,vt3dg
         if (jdim .eq. 1)  &
         call fa_yc(mzp,mxp,myp,ia,iz,ja,jz,g3d_g(ngrid)%cnv_tr,scr1,vt3db,vt3de,vt3dg,vt3dj,vt3di,jdim,mynum)
       
         ! input: scalarp, scr1,vt3dc,vt3df,vt3dk, vctr1,vctr2
         ! output: scr1,vt3dg
         if(vert_adv) &
         call fa_zc(mzp,mxp,myp,ia,iz,ja,jz,g3d_g(ngrid)%cnv_tr,scr1,vt3dc,vt3df,vt3dg,vt3dk,vctr1,vctr2,mynum)
         
         !--- get the adection tendency, if desired
         !call advtndc(mzp,mxp,myp,ia,iz,ja,jz,g3d_g(ngrid)%cnv_tr,scr1,dummy,dtlt,mynum)

         !----------------------
         !Aproveitando as rotinas de comunicação criadas para passar o Halo quando em uso da RK
         !para comunicar o scr1 aos demais processadores
         !if (mynum==1) then 
         !  scr1(5,25,25)=8.8
         !    write(88,*) scr1(5,25,25);call flush(88)
         !endif
         !if (mynum==2) then 
         !  scr1(5,1,25)=5.5
         !  write(89,*) scr1(5,1,25);call flush(89)
         !endif
         !----------------------
        
         !--- apply 3, 2 or 1-d diffusion
         if(add_coldpool_diff >= 1) then 
       
           !scr1 = scalar to be applied 3-D, 2-D or 1-D diffusion
           scalart = 0. 
           sfcflx  = 0.
          
           call sclr_diffusion(add_coldpool_diff, mzp, mxp, myp, ia, iz, ja, jz, dtlt,dzm      &
                              ,dzt,grid_g(ngrid)%rtgt,scr1, scalart, sfcflx, basic_g(ngrid)%dn0 &  
                              ,turb_g(ngrid)%vkh, turb_g(ngrid)%hkm,xkhkm(ngrid)                )

           !--- update with the diffusion process
           scr1 (:,:,:)= scr1(:,:,:)+scalart(:,:,:)*dtlt
         endif
         
         !--- avoid neg scalar
         where(scr1 < 0.)  scr1 = 0.0

         !--- parallel communication 
         if (.not. singleProcRun .and. n < nsubsteps) call commHaloAcou(scr1,mxp,myp,mzp,myNum,'scr1')
        
         !--- updated tracer CNV_TR = CNV_TR + ADV * dtlt + DIFF * dtlt
         g3d_g(ngrid)%cnv_tr = scr1

     enddo

     if(use_gustiness==1) then

      !--- water vapor mixing ratio
      scr1(1:mzp,1:mxp,1:myp) = g3d_g(ngrid)%qexcp(1:mzp,1:mxp,1:myp)

      do n=1,nsubsteps

         ! output: scr1,vt3dg
         call fa_xc(mzp,mxp,myp,ia,iz,1,myp,g3d_g(ngrid)%qexcp,scr1,vt3da,vt3dd,vt3dg,vt3dh,vt3di,mynum)

         ! input: scalarp, scr1,vt3db,vt3de,vt3dj,vt3di
         ! output: scr1,vt3dg
         if (jdim .eq. 1)  &
         call fa_yc(mzp,mxp,myp,ia,iz,ja,jz,g3d_g(ngrid)%qexcp,scr1,vt3db,vt3de,vt3dg,vt3dj,vt3di,jdim,mynum)
       
         ! input: scalarp, scr1,vt3dc,vt3df,vt3dk, vctr1,vctr2
         ! output: scr1,vt3dg
         if(vert_adv) &
         call fa_zc(mzp,mxp,myp,ia,iz,ja,jz,g3d_g(ngrid)%qexcp,scr1,vt3dc,vt3df,vt3dg,vt3dk,vctr1,vctr2,mynum)

         where(scr1 < 0.)  scr1 = 0.0
         
         if (.not. singleProcRun .and. n < nsubsteps) then
           call commHaloAcou(scr1,mxp,myp,mzp,myNum,'scr1')
         endif
         
         !--- updated tracer CNV_TR = CNV_TR + ADV * dtlt
         g3d_g(ngrid)%qexcp = scr1
       
      enddo
      
      !--- moist static energy
      scr1(1:mzp,1:mxp,1:myp) = g3d_g(ngrid)%hexcp(1:mzp,1:mxp,1:myp)

      do n=1,nsubsteps

         ! output: scr1,vt3dg
         call fa_xc(mzp,mxp,myp,ia,iz,1,myp,g3d_g(ngrid)%hexcp,scr1,vt3da,vt3dd,vt3dg,vt3dh,vt3di,mynum)

         ! input: scalarp, scr1,vt3db,vt3de,vt3dj,vt3di
         ! output: scr1,vt3dg
         if (jdim .eq. 1)  &
         call fa_yc(mzp,mxp,myp,ia,iz,ja,jz,g3d_g(ngrid)%hexcp,scr1,vt3db,vt3de,vt3dg,vt3dj,vt3di,jdim,mynum)
       
         ! input: scalarp, scr1,vt3dc,vt3df,vt3dk, vctr1,vctr2
         ! output: scr1,vt3dg
         if(vert_adv) &
         call fa_zc(mzp,mxp,myp,ia,iz,ja,jz,g3d_g(ngrid)%hexcp,scr1,vt3dc,vt3df,vt3dg,vt3dk,vctr1,vctr2,mynum)

         where(scr1 < 0.)  scr1 = 0.0
         
         if (.not. singleProcRun .and. n < nsubsteps) then
           call commHaloAcou(scr1,mxp,myp,mzp,myNum,'scr1')
         endif
         
         !--- updated tracer CNV_TR = CNV_TR + ADV * dtlt
         g3d_g(ngrid)%hexcp = scr1
       
      enddo
  endif ! - if gustiness = 1


end subroutine adv_convection_tracer

!---------------------------------------------------------------
subroutine sclr_diffusion(add_coldpool_diff, m1, m2, m3, ia, iz, ja, jz, dtlt,dzm,dzt &
                          ,rtgt, scp, sct, sfcflx, dn0, vkkh, hkkh, xkhkm )

  implicit none
  integer  , intent(in) :: m1, m2, m3, ia, iz, ja, jz,add_coldpool_diff
  real     , intent(in) :: dtlt,xkhkm
  real     , intent(in) :: scp(m1,m2,m3), vkkh(m1,m2,m3), hkkh(m1,m2,m3),dn0(m1,m2,m3)

  real, intent(in)    :: dzm(m1), dzt(m1)
  real, intent(in)    :: sfcflx(m2,m3),rtgt(m2,m3)
 
  real, intent(inout) :: sct(m1,m2,m3)

  real                :: vt3da(m1,m2,m3), vt3db(m1,m2,m3),  &
       vt3dc(m1,m2,m3), vt3df(m1,m2,m3), vt3dg(m1,m2,m3), vt3dj(m1,m2,m3), &
       vt3dk(m1,m2,m3), vt3do(m1,m2,m3), dn03i(m1,m2,m3), vt3dl(m1,m2,m3), &
       vt3dm(m1,m2,m3)
 
  real :: vt2db(m2,m3)

  real, dimension(m1) :: vctr1  
       
  integer :: i, j, k, ihorgrad,jd
  real :: c1
  real :: dtlti

  !          compute vertical diffusion matrix coefficients for scalars

  JD=1
  c1 = 0.5*dtlt
  !ihorgrad=0  !no horizontal difusion
  !nstbot=1
  !nsttop=1
  
  !vt3df=0.! check necessidade
  !vt3dg=0.! check necessidade

   
     do j=ja,jz
        do i=ia,iz
           do k=1,m1-1
              vctr1(k) = dzm(k)*(vkkh(k,i,j) + vkkh(k+1,i,j))
           enddo
 
           do k=2,m1-1
              
              vt3dj(k,i,j) = -c1 * dzt(k) * vctr1(k-1)
              vt3dk(k,i,j) = -c1 * dzt(k) * vctr1(k)
             
              vt3do(k,i,j) = dn0(k,i,j) - vt3dj(k,i,j) - vt3dk(k,i,j)
              dn03i(k,i,j) = 1./dn0(k,i,j)

           enddo
        enddo
     enddo

!     if (nstbot==1 .and. nsttop==1) then
        do j=ja,jz
           do i=ia,iz
              vt3dj(2,i,j) = 0.
              vt3do(2,i,j) = dn0(2,i,j) - vt3dk(2,i,j)
              vt3dk(m1-1,i,j) = 0.
              vt3do(m1-1,i,j) = dn0(m1-1,i,j) - vt3dj(m1-1,i,j)

              vt2db(i,j) = 1./vt3do(2,i,j)
              vt3dl(2,i,j) = vt3dk(2,i,j)*vt2db(i,j)
           enddo
        enddo
!     elseif (nstbot==1) then
!        do j=ja,jz
!           do i=ia,iz
!              vt3dj(2,i,j)    = 0.
!              vt3do(2,i,j)    = dn0(2,i,j) - vt3dk(2,i,j)
!              vt3dk(m1-1,i,j) = 0.
!
!         ! Bob (10/27/00):  Evaluate ci2i, cjp1, and cji terms here
!         ! (vt2db, vt3dl, vt3dm) for new tridiff2.
!         ! With this, no longer need to pass ci or cip1 (vt3do, vt3dk)
!         ! to tridiff2.
!
!         vt2db(i,j)      = 1./vt3do(2,i,j)
!         vt3dl(2,i,j)    = vt3dk(2,i,j)*vt2db(i,j)
!      enddo
!   enddo
!     elseif (nsttop==1) then
!   do j=ja,jz
!      do i=ia,iz
!         vt3dj(2,i,j)    = 0.
!         vt3dk(m1-1,i,j) = 0.
!         vt3do(m1-1,i,j) = dn0(m1-1,i,j) - vt3dj(m1-1,i,j)
!
!         ! Bob (10/27/00):  Evaluate ci2i, cjp1, and cji terms here
!         ! (vt2db, vt3dl, vt3dm) for new tridiff2.
!         ! With this, no longer need to pass ci or cip1 (vt3do, vt3dk)
!         ! to tridiff2.
!
!         vt2db(i,j)      = 1./vt3do(2,i,j)
!         vt3dl(2,i,j)    = vt3dk(2,i,j)*vt2db(i,j)
!      enddo
!   enddo
!     else
!   do j=ja,jz
!      do i=ia,iz
!         vt3dj(2,i,j)    = 0.
!         vt3dk(m1-1,i,j) = 0.
!
!         ! Bob (10/27/00):  Evaluate ci2i, cjp1, and cji terms here
!         ! (vt2db, vt3dl, vt3dm) for new tridiff2.
!         ! With this, no longer need to pass ci or cip1 (vt3do, vt3dk)
!         ! to tridiff2.
!
!         vt2db(i,j)      = 1./vt3do(2,i,j)
!         vt3dl(2,i,j)    = vt3dk(2,i,j)*vt2db(i,j)
!      enddo
!   enddo
!     endif
!
     do j=ja,jz
        do i=ia,iz
           do k=3,m1-1
              vt3dm(k,i,j) = 1./(vt3do(k,i,j) - vt3dj(k,i,j)*vt3dl(k-1,i,j))
              vt3dl(k,i,j) = vt3dk(k,i,j)*vt3dm(k,i,j)
           enddo
        enddo
     enddo

!  endif

if (add_coldpool_diff >= 2) then 
!------------------------- horizontal difusion------------------
!     compute 2 horizontal scalar gradients needed for dscp/dt
    call grad(m1, m2, m3,  1, iz, ja, jz, scp, vt3df, 'XDIR', 'TPNT')
    call grad(m1, m2, m3, ia, iz,  1, jz, scp, vt3dg, 'YDIR', 'TPNT')
!
    do j=ja,jz
      do i=1,iz
       do k=1,m1-1
          vt3df(k,i,j) = -vt3df(k,i,j)*0.5*(hkkh(k,i,j) + hkkh(k,i+1,j))* xkhkm
       enddo
      enddo
    enddo

    do j=1,jz
      do i=ia,iz
       do k=1,m1-1
          vt3dg(k,i,j) = -vt3dg(k,i,j)*0.5*(hkkh(k,i,j) + hkkh(k,i,j+jd))* xkhkm
       enddo
      enddo
    enddo
!
!       horizontal flux divergence for scalars
    call divcart(m1, m2, m3, ia, iz, ja, jz, vt3df, vt3da, 'XDIR', 'UPNT')
    call divcart(m1, m2, m3, ia, iz, ja, jz, vt3dg, vt3db, 'YDIR', 'VPNT')

    do j=ja,jz
      do i=ia,iz
        do k=2,m1-1
          sct(k,i,j) = sct(k,i,j)  - (vt3da(k,i,j) + vt3db(k,i,j))*dn03i(k,i,j)
        enddo
      enddo
    enddo
!---------------------------------------------------------------
 endif 
 
 if (add_coldpool_diff >= 1) then 
!---------------------------------------------------------------
  !  finish matrix coefficients for vertical diffusion

   do j=ja,jz
     do i=ia,iz
        do k=2,m1-1
           vt3dc(k,i,j) = scp(k,i,j)*dn0(k,i,j)  ! scalar being diffused
        enddo
   
           vt3dc(2,i,j) = scp(2,i,j)*dn0(2,i,j) + & !surface fluxes
                          sfcflx(i,j)*dtlt*dzt(2)/rtgt(i,j)

          ! vt3dc(2,i,j) = max(0.,vt3dc(2,i,j))!positividade

     enddo
   enddo

   do j=ja,jz
     do i=ia,iz

        vt3df(2,i,j) = vt3dc(2,i,j)*vt2db(i,j)

        do k=3,m1-1
           vt3df(k,i,j) = (vt3dc(k,i,j) - vt3dj(k,i,j)*vt3df(k-1,i,j))* &
                vt3dm(k,i,j)
        enddo

        do k=m1-2,2,-1
           vt3df(k,i,j) = vt3df(k,i,j) - vt3dl(k,i,j)*vt3df(k+1,i,j)
        enddo

     enddo
   enddo

   dtlti = 1.0/dtlt
   do j=ja,jz
     do i=ia,iz
        do k=2,m1-1
           sct(k,i,j) = sct(k,i,j) + vt3df(k,i,j)*dtlti - scp(k,i,j)*dtlti
        enddo
     enddo
   enddo
!--------------------------------------------------------------
 endif
end subroutine sclr_diffusion

!------------------------------------------------------------------------------------
subroutine sgs_coldpool_sfx_and_diff()

   use cuparm_grell3, only: g3d_g
   use mem_grid     , only: time,ng=>ngrid,dtlt, dyncore_flag,dzt,ztn
   use node_mod     , only: mxp,myp,mzp ,ia,iz,ja,jz,mynum, nmachs
   use mem_grid     , only: ngrid, nzpmax, grid_g, dtlt, if_adap, jdim, time, &
                            zt, zm, dzm, dzt, hw4,itopo,npatch
   use mem_basic    , only: basic_g
   use mem_leaf     , only: leaf_vars, leaf_g 
   use mem_turb     , only: turb_vars, turb_g , xkhkm

   use convpar_gf_geos5, only: coldpool_start,tau_land_cp,tau_ocea_cp, use_gustiness &
                              ,add_coldpool_diff

   use mem_jules       , only: jules_g
   use rconstants   
   !use leaf_coms, only: timefac_sst, niter_leaf, niter_can, dtll, dtll_factor, &
   !    dtlc_factor, dtlc, hcapcan, wcapcan, hcapveg, dtllohcc, dtllowcc, &
   !    dtlcohcc, dtlcowcc, dtlcohcv, z0fac_water, snowrough, ups, vps, ths, rvs, &
   !!    zts, pis, dens, prss, vels, gzotheta, tempk, snowfac, thetacan, &
    !   vels_pat, rdi, pcpgl,fracliq

   implicit none
   include "i8.h"
   integer ::  i,j,ip,niter_leaf, niter_can,k
   real :: factor,snk,src,tau_cp
   real                      :: ths,rvs,pis,dens,ups,vps,zts, gzotheta,prss,temps, dummy=0.&
                               ,vels,patch_area,  dtll, dtll_factor,dtllohcc,vels_pat &
                               ,dtlcohcc, dtlc_factor, dtlc, hcapcan , thetacan, can_rvap &
                               ,ustar0,tstar0, rstar0,sflux_u,sflux_v, sflux_w, sflux_t, sflux_r &
                               ,z0fac_water,vels_ratio
   real, dimension (mxp,myp) :: ustar,tstar, rstar,patch_rough
   real, dimension (mxp,myp) :: ths2,rvs2,pis2,dens2,ups2,vps2,zts2,vels1,vels_cp
   real, dimension (mxp,myp,npatch) :: l_patch_rough
   
   real, dimension (mzp,mxp,myp) :: scalart,scalarp,scr1
   real, dimension (mxp,myp)     :: sfcflx
   
   vels1=0.;vels_cp=0.
   do j=ja,jz
      do i=ia,iz
       dens= (basic_g(ng)%dn0(1,i,j)   + basic_g(ng)%dn0(2,i,j)) * .5
       ups = (basic_g(ng)%up (2,i-1,j) + basic_g(ng)%up(2,i,j) ) * .5
       vps = (basic_g(ng)%vp (2,i,j-1) + basic_g(ng)%vp(2,i,j) ) * .5

       vels1(i,j) = sqrt(ups ** 2 + vps ** 2) + 1.e-6

       !--- include the cold-pool gust front effects on the surface fluxes
       factor=  coldpool_start(g3d_g(ng)%CNV_TR(2,i,j))
       
       ups = (1.-factor)*ups+factor*g3d_g(ng)%umcl(i,j)
       vps = (1.-factor)*vps+factor*g3d_g(ng)%vmcl(i,j)

       vels_cp(i,j) = sqrt(ups ** 2 + vps ** 2)
       
       !print*,"LE",alvl* turb_g(ng)%sflux_r(i,j), alvl * ustar0*rstar0 * dens
       !print*," H",   cp*turb_g(ng)%sflux_t(i,j),cp * sflux_t
       
       vels_ratio = vels_cp(i,j)/vels1(i,j)

       !-- sensible heat flux K m/s
       g3d_g(ng)%sgs_shf(i,j) = turb_g(ng)%sflux_t(i,j)*vels_ratio !*dens*cp   = H W/m2
       
       !-- latent heat flux
       g3d_g(ng)%sgs_lhf(i,j) = turb_g(ng)%sflux_r(i,j)*vels_ratio !*dens*alvl = LE W/m2

   enddo; enddo

   if(use_gustiness  == 1) then
        
      tau_cp = 1.* tau_land_cp - (1.-1.)*tau_ocea_cp ! check latter 
        
      do j=ja,jz
         do i=ia,iz
             !-- apply the source term at surface

             ! a) water vapor mixing ratio 
             src= min(2000./alvl, (vels_cp(i,j)/vels1(i,j) - 1.) * turb_g(ng)%sflux_r(i,j)) &
                    * dtlt *dzt(2)/grid_g(ng)%rtgt(i,j)
             
             g3d_g(ng)%qexcp(2,i,j) = g3d_g(ng)%qexcp(2,i,j) + src 
                       

             ! b) moist static energy (using hexcp array)
             src = cp *  turb_g(ng)%sflux_t(i,j) +  alvl *  turb_g(ng)%sflux_r(i,j)
             src = min (3000., src) * (vels_cp(i,j)/vels1(i,j) - 1.) 
             src = src * dtlt *dzt(2)/grid_g(ng)%rtgt(i,j)

             g3d_g(ng)%hexcp(2,i,j) = g3d_g(ng)%hexcp(2,i,j) + src 
            

            !-- apply the sink term at all column
            !
            do k = 1, mzp
                snk= dtlt * abs(g3d_g(ng)%qexcp(k,i,j))/tau_cp
                g3d_g(ng)%qexcp(k,i,j) = g3d_g(ng)%qexcp(k,i,j) - snk

                snk= dtlt * abs(g3d_g(ng)%hexcp(k,i,j))/tau_cp
                g3d_g(ng)%hexcp(k,i,j) = g3d_g(ng)%hexcp(k,i,j) - snk
               
            enddo
         enddo
      enddo

      !--- apply the vertical diffusion term
      !
       if(add_coldpool_diff >= 1) then 
       
       do k =1,2
          if(k==1) scalarp = g3d_g(ng)%qexcp ! scalar to be applied vert diffusion
          if(k==2) scalarp = g3d_g(ng)%hexcp ! scalar to be applied vert diffusion
          scalart = 0. 
          sfcflx  = 0.

          call sclr_diffusion(add_coldpool_diff, mzp, mxp, myp, ia, iz, ja, jz, dtlt,dzm,dzt   &
                             ,grid_g(ng)%rtgt,scalarp, scalart, sfcflx                         &
                            , basic_g(ng)%dn0,  turb_g(ng)%vkh , turb_g(ng)%hkm, xkhkm(ng) )

          !--- update
          if(k==1) g3d_g(ng)%qexcp (:,:,:)= g3d_g(ng)%qexcp (:,:,:)+scalart(:,:,:)*dtlt
          if(k==2) g3d_g(ng)%hexcp (:,:,:)= g3d_g(ng)%hexcp (:,:,:)+scalart(:,:,:)*dtlt
        enddo
      endif
   
      !--- sanity check for q and MSE
      where (g3d_g(ng)%qexcp < 0.          ) g3d_g(ng)%qexcp = 0.
      where (g3d_g(ng)%qexcp > 20.e-3      ) g3d_g(ng)%qexcp = 20.e-3
      where (g3d_g(ng)%hexcp < 0.          ) g3d_g(ng)%hexcp = 0.
      where (g3d_g(ng)%hexcp > 25.e-3*alvl ) g3d_g(ng)%hexcp = 25.e-3*alvl

   endif ! if gustiness == 1




   !print*,"vels",maxval(abs(vels_cp)),maxval(abs(vels1)),minval(abs(vels_cp)), minval(abs(vels1))
   !print*,"HLE ",maxval(abs(g3d_g(ng)%sgs_shf))*1.1*cp,maxval(abs(g3d_g(ng)%sgs_lhf))*1.1*alvl
   !call flush(6)
   return

  !----- for future use-----------------------------------------------------------

  !-------------- check this later 
  ! Define leaf3 and canopy time-split timesteps here.  This ensures that leaf3
  ! will not use a timestep longer than about 40 seconds, and canopy will not
  ! use a timestep longer than about 15 seconds.  This allows values of
  ! hcapcan = 2.e4, wcapcan = 2.e1, and hcapveg = 3.e4 as are now defined below.

  niter_leaf  = max(1,nint(dtlt/40.+.4))
  niter_can   = max(1,nint(dtll/15.+.4))
  dtll_factor = 1. / float(niter_leaf)
  dtll        = dtlt * dtll_factor
  dtlc_factor = 1. / float(niter_can)
  dtlc        = dtll * dtlc_factor
  hcapcan = 2.0e4
  dtlcohcc = dtlc / hcapcan
  dtllohcc = dtll / hcapcan
!----------------------------
  dtll        = dtlt
 

   !z0fac_water = .016 / g
   ! Copy surface atmospheric variables into 2d arrays for input to leaf
   call sfc_fields(mzp,mxp,myp,ia,iz,ja,jz,jdim                                     &
          ,basic_g(ng)%theta(1,1,1) ,basic_g(ng)%rv (1,1,1) ,basic_g(ng)%up(1,1,1)  &
          ,basic_g(ng)%vp   (1,1,1) ,basic_g(ng)%dn0(1,1,1) ,basic_g(ng)%pp(1,1,1)  &
          ,basic_g(ng)%pi0  (1,1,1) ,grid_g (ng)%rtgt(1,1)  ,zt                     &
          ,ths2,rvs2,ups2,vps2,pis2,dens2,zts2)    

   do j=ja,jz
      do i=ia,iz
       ustar   (i,j) = 0.
       tstar   (i,j) = 0.          
       rstar   (i,j) = 0.
       patch_rough(i,j) = 0.

      ! For water surface (patch 1), compute surface saturation mixing ratio
      ! and roughness length based on previous ustar.
      ! For soil patches, compute roughness length based on vegetation and snow.

      !leaf%ground_rsat(i,j,ip) = rslif(prss,tempk(mzg))
       l_patch_rough(i,j,1) = max(z0fac_water * leaf_g(ng)%ustar(i,j,1) ** 2,.0001)

         
       l_patch_rough(i,j,2:npatch) = max(grid_g(ng)%topzo(i,j),leaf_g(ng)%soil_rough(i,j,2:npatch)  &
                                                              ,leaf_g(ng)%veg_rough (i,j,2:npatch))  
                        
       ustar0=leaf_g(ng)%ustar (i,j,1)
       tstar0=leaf_g(ng)%tstar (i,j,1)
       rstar0=leaf_g(ng)%rstar (i,j,1)



       do ip=1,npatch
        ustar (i,j) = ustar (i,j) + leaf_g(ng)%ustar (i,j,ip)*leaf_g(ng)%patch_area(i,j,ip)
        tstar (i,j) = tstar (i,j) + leaf_g(ng)%tstar (i,j,ip)*leaf_g(ng)%patch_area(i,j,ip)            
        rstar (i,j) = rstar (i,j) + leaf_g(ng)%rstar (i,j,ip)*leaf_g(ng)%patch_area(i,j,ip)     
        patch_rough(i,j) =patch_rough(i,j) +  l_patch_rough(i,j,ip)*leaf_g(ng)%patch_area(i,j,ip) 
        !print*," rough=",ip,patch_rough(i,j) ,l_patch_rough(i,j,ip), leaf_g(ng)%patch_area(i,j,ip) 
       enddo
      
       ups = ups2(i,j)
       vps = vps2(i,j)
       ths = ths2(i,j)
       rvs = rvs2(i,j)
       zts = zts2(i,j)
       pis = pis2(i,j)
       dens = dens2(i,j)
       gzotheta = g * zts / ths
       prss = pis ** cpor * p00
       temps = ths * pis 
   

      !--- include the cold-pool gust front effects on the surface fluxes
       factor=  coldpool_start(g3d_g(ng)%CNV_TR(2,i,j))
       
       ups = (1.-factor)*ups2(i,j)+factor*g3d_g(ng)%umcl(i,j)
       vps = (1.-factor)*vps2(i,j)+factor*g3d_g(ng)%vmcl(i,j)

       vels = sqrt(ups ** 2 + vps ** 2)
       
       patch_area = 1.

!jules_g(ng)%u10mj(ia:iz,ja:jz)=sf_diag%u10m(:,:)
!jules_g(ng)%v10mj(ia:iz,ja:jz)=sf_diag%v10m(:,:)
!jules_g(ng)%t2mj(ia:iz,ja:jz)=sf_diag%t1p5m(:,:)
!jules_g(ng)%rv2mj(ia:iz,ja:jz)=sf_diag%q1p5m(:,:)/(1-sf_diag%q1p5m(:,:))

       thetacan = jules_g(ng)%t2mj(i,j)/pis
       can_rvap = jules_g(ng)%rv2mj(i,j)

     !  call stars(ustar(i,j),tstar(i,j) ,rstar(i,j)     &
     !            ,ths,rvs,thetacan,can_rvap             &
     !            ,zts,patch_rough(i,j),patch_area       &
     !            ,vels,vels_pat,vonk,dummy,dens,dtll,dummy )

!print*,"us",ustar(i,j),ustar0,tstar(i,j),tstar0,ths-thetacan,can_rvap,rvs

       sflux_u = 0.; sflux_v = 0.;  sflux_w = 0.;  sflux_t = 0.; sflux_r = 0.
       call sfclmcv(i,j,ustar(i,j), tstar(i,j),                  &
           rstar(i,j), vels, vels_pat, ups, vps, gzotheta,       &
           patch_area, sflux_u,sflux_v, sflux_w, sflux_t, sflux_r)

       sflux_u = sflux_u * dens
       sflux_v = sflux_v * dens 
       sflux_w = sflux_w * dens
       sflux_t = sflux_t * dens
       sflux_r = sflux_r * dens 
      
       print*,"LE",alvl* turb_g(ng)%sflux_r(i,j), alvl * ustar0*rstar0 * dens
       !print*," H",   cp*turb_g(ng)%sflux_t(i,j),cp * sflux_t

       g3d_g(ng)%sgs_shf(i,j) = cp   * ustar0*tstar0 * dens !sflux_t
       g3d_g(ng)%sgs_lhf(i,j) = alvl * ustar0*rstar0 * dens !sflux_r


   enddo; enddo
!stop 3

end subroutine sgs_coldpool_sfx_and_diff
!-----------------------------------------------------------------------------------------
