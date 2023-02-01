   subroutine get_zi_gf(j, its, ite, kts, kte, istart, iend, ktf, ierr, kzi, pbl, tkeg, &
                        rcpg, z, ztop, tkmin)

      implicit none
      integer, intent(in):: j, its, ite, kts, kte, ktf, istart, iend
      integer :: kzimax, ktke_max, i, k
      real tkmin, tke_tmp
      real, dimension(its:ite, kts:kte) :: tkeg, rcpg, z
      real, dimension(its:ite)          :: ztop, pbl
      integer, dimension(its:ite)          :: kzi, ierr

      real, parameter :: rcpmin = 1.e-5, pblhmax = 3000.

      do i = istart, iend
         kzi(i) = 1
         !    if(ierr(i).eq.0)then
         !         tke_tmp = 0.
         ktke_max = 1
         kzimax = ktf - 1
         !---  max level for kzi
         do K = kts, ktf
            if (z(i, k) .ge. pblhmax + ztop(i)) then
               kzimax = min(k, ktf - 1)
               !if(j==8 .and. i==10) print*,"1",z(i,k), pblhmax,ztop(i),kzimax
               exit
            end if
         end do
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

         do k = ktke_max, kzimax + 1

            !           if(tkeg(i,k) .gt. 1.1*tkmin .and. rcpg(i,k) .lt. rcpmin )  then
            if (tkeg(i, k) .gt. 1.1*tkmin) then
               kzi(i) = k
               !if(j==8 .and. i==10) print*,"I",k,rcpg(i,k),tkeg(i,k),kzi(i),z(i,k)-ztop(i)
               cycle

            else
               kzi(i) = max(1, k - 1)
               !if(j==8 .and. i==10) print*,"F",k,rcpg(i,k),tkeg(i,k),kzi(i),z(i,k)-ztop(i)
               exit
            end if

         end do
         kzi(i) = max(1, kzi(i))
         !print*,"1",kzi(i),i
         kzi(i) = min(kzimax, kzi(i))
         !print*,"2",kzi(i),icall flush(6)
         pbl(i) = max(z(i, kzi(i)) - ztop(i), z(i, 1) - ztop(i))
      end do

   end subroutine get_zi_gf

      !------------------------------------------------------------------------------------

   subroutine get_zu_zd_pdf_orig(draft, ierr, kb, kt, zs, zuf, ztop, zu, kts, kte, ktf)

      implicit none
      integer, intent(in) ::kb, kt, kts, kte, ktf
      real, intent(in) :: Zs, Zuf, Ztop
      real, intent(inout) :: zu(kts:kte)
      integer, intent(inout) :: ierr
      character*(*), intent(in) ::draft

      !- local var
      integer :: add, i, nrec = 0, k, kb_adj
      real ::zumax, ztop_adj
      real ::beta, alpha, kratio, tunning

      !- kb cannot be at 1st level
      kb_adj = max(kb, 2)

      !-- fill zu with zeros
      zu = 0.0

      if (draft == "UP" .or. draft == "up") then
         if (kt <= kb_adj) then
            !stop "ktop must be larger than kbcon"
            ierr = 99
            return
         end if
         !beta=4.  !=> must larger than 1
         !=> higher makes the profile sharper
         !=> around the maximum zu
         add = 0     !=> additional levels above kbcon, where
         !=> the maximum zu will resides
         kb_adj = kb_adj + add
         kb_adj = max(10, kb_adj)

         !- this alpha constrains the location of the maximun ZU to be at
         !- "kb_adj" vertical level
         !alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))

         !- 2nd approach for beta and alpha parameters
         !- the tunning parameter must be between 0.5 (low  level max zu)
         !-                                   and 1.5 (high level max zu)
         tunning = 0.6
         beta = 2.0/tunning
         alpha = tunning*beta

         !- Beta PDF
         do k = kts, min(kte, kt + 1)
            kratio = float(k)/float(kt + 1)

            zu(k) = kratio**(alpha - 1.0)*(1.0 - kratio)**(beta - 1.0)
         end do

      elseif (draft == "DOWN" .or. draft == "DOWNM") then
         add = 0    !=> additional levels above kbcon, where
         !=> the maximum zu will resides
         beta = 4.  !=> must larger than 1
         !=> higher makes the profile sharper
         !=> around the maximum zu
         alpha = 0.25*beta

         !- for downdrafts kb = jmin(i)-levadj
         kb_adj = kb_adj + add

         !- 2nd approach for beta and alpha parameters
         !- the tunning parameter must be between 0.5 (low  level max zu)
         !-                                   and 1.5 (high level max zu)
         tunning = 1.
         beta = 2.0/tunning
         alpha = tunning*beta

         !- Beta PDF
         do k = kts, min(kte, kt)
            kratio = float(k)/float(kt)

            zu(k + 1) = kratio**(alpha - 1.0)*(1.0 - kratio)**(beta - 1.0)
         end do

      elseif (draft == "shallow" .or. draft == "SHALLOW") then

         alpha = 3.
         beta = 2.*alpha
         kb_adj = 1 ! level where mass flux starts

         !- Beta PDF
         do k = kts + kb_adj - 1, min(kte, kt + 1)
            kratio = float(k + 1 - kb_adj)/float(kt + 1)  !-kb_adj+1)

            zu(k) = kratio**(alpha - 1.0)*(1.0 - kratio)**(beta - 1.0)
         end do

      else
         print *, "unknown type of flow", draft
         stop "routine get_zu_zd"

      end if

      !- normalize ZU
      zu(kts:min(kte, kt + 1)) = zu(kts:min(kte, kt + 1))/maxval(zu(kts:min(kte, kt + 1)))

      !--- Sanity checks
      if (beta <= 1) stop "beta must be larger than 1"

      if (minval(zu(:)) < 0.0) then
         print *, " zu has negative values for ", draft
         stop " zu < zero"
      end if
      if (maxval(zu(:)) > 1.0) then
         print *, " zu has values greater than 1 for ", draft
         stop " zu  >  one"
      end if

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
   subroutine FLIPZ(flip, mzp)
      implicit none
      integer, intent(in) :: mzp
      integer, dimension(mzp), intent(inout) :: flip
      integer :: m, k
      m = mzp
      do k = 1, mzp
         flip(k) = m
         m = m - 1
      end do
   end subroutine FLIPZ

      !------------------------------------------------------------------------------------

   subroutine set_index_loops(ims, ime, jms, jme, kms, kme, &
                              its, ite, jts, jte, kts, kte, &
                              mxp, myp, mzp)

      implicit none
      integer, intent(in)         :: mxp, myp, mzp
      integer, intent(inout)      :: ims, ime, jms, jme, kms, kme, &
                                     its, ite, jts, jte, kts, kte

      ims = 1
      ime = mxp
      jms = 1
      jme = myp
      kms = 1
      kme = mzp
      its = 1
      ite = mxp
      jts = 1
      jte = myp
      kts = 1
      kte = mzp

   end subroutine set_index_loops

      !------------------------------------------------------------------------------------

   subroutine get_vars(LM, mxp, myp, Q, T, PLE, ZLE, ZLO, PLO, PK, TH)
      implicit none
      integer, intent(in) :: LM, mxp, myp
      real, intent(in), dimension(mxp, myp, 0:LM) :: PLE
      real, intent(in), dimension(mxp, myp, LM)   :: T, Q

      real, intent(out), dimension(mxp, myp, 0:LM) :: ZLE
      real, intent(out), dimension(mxp, myp, LM)   :: ZLO, PLO, PK, TH

      !-local var
      integer :: L
      real, dimension(mxp, myp, 0:LM) :: CNV_PLE, PKE

      CNV_PLE = PLE*.01
      PLO = 0.5*(CNV_PLE(:, :, 0:LM - 1) + CNV_PLE(:, :, 1:LM))
      PKE = (CNV_PLE/1000.)**(c_rgas2/c_cp2)
      PK = (PLO/1000.)**(c_rgas2/c_cp2)
      TH = T/PK

      ZLE(:, :, LM) = 0.
      do L = LM, 1, -1
         ZLE(:, :, L - 1) = TH(:, :, L)*(1.+c_vireps*Q(:, :, L))
         ZLO(:, :, L) = ZLE(:, :, L) + (c_cp2/c_grav)*(PKE(:, :, L) - PK(:, :, L))*ZLE(:, :, L - 1)
         ZLE(:, :, L - 1) = ZLO(:, :, L) + (c_cp2/c_grav)*(PK(:, :, L) - PKE(:, :, L - 1))*ZLE(:, :, L - 1)
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

   !----------------------------------------------------------------------
      subroutine alloc_grads_arr(n, mzp, task, jl)
      implicit none
      integer, intent(in)    :: n, mzp, task
      integer, intent(inout) :: jl
      integer :: nvar

      if (task == 1) then
         jl = n
         allocate (cupout(nvar_grads))
         do nvar = 1, nvar_grads
            allocate (cupout(nvar)%varp(n, mzp))
            allocate (cupout(nvar)%varn(3))
            cupout(nvar)%varp(:, :) = 0.0
            cupout(nvar)%varn(:) = "xxxx"
         end do
      else
         do nvar = 1, nvar_grads
            deallocate (cupout(nvar)%varp)
            deallocate (cupout(nvar)%varn)
         end do
         deallocate (cupout)
      end if

   end subroutine alloc_grads_arr

      !------------------------------------------------------------------------------------
   subroutine writetxt(mzp, t, ple, th1, pk, q1, u1, zle, zlo &
                       , DYNF_Q, DYNF_T, DYNF_PLE, DYNF_UA &
                       !
                       , theta, pp, rv, up, zm3d, zt3d, vp, omega &
                       , sflux_r, sflux_t, topt, xland, sfc_press, dx, kpbl, temp2m, dt_moist &
                       )

      implicit none
      integer, intent(in)    :: mzp
      real, dimension(mzp)   :: t, th1, pk, q1, u1, zlo &
                                , DYNF_Q, DYNF_T, DYNF_UA

      real, dimension(mzp)   :: theta, pp, rv, up, zm3d, zt3d, vp, omega
      real, dimension(0:mzp) :: ple, zle, DYNF_PLE
      real :: sflux_r, sflux_t, topt, xland, sfc_press, dx, temp2m, dt_moist

      integer :: k, kpbl

      write (8, *) "================================================"
      write (7, *) "================================================"
      write (7, *) kpbl, sflux_r, sflux_t, topt, xland, sfc_press, dx, temp2m, dt_moist
      do k = 1, mzp
         write (8, 10) k, ple(k), t(k), th1(k), pk(k), 1000.*q1(k), u1(k), zle(k), zlo(k), 86400.*DYNF_Q(k)
         write (7, 11) k, theta(k), pp(k), 1000.*rv(k), up(k), zm3d(k), zt3d(k), vp(k), omega(k)
      end do
      call flush (7)
      call flush (8)
10    format(1x, i4, 9f11.3)
11    format(1x, i4, 8f11.3)
   end subroutine writetxt

      !------------------------------------------------------------------------------------
   subroutine get_cloud_fraction(mzp, kts, ktf, &
                                 PPABS, PZZ, PT, PRV, QCO, QRCO, PMFLX, PCLDFR, PRC, PRI)
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
      real, parameter :: XP00 = 1.e5        ! reference pressure
      real, parameter :: XPI = 3.141592654 ! Pi
      real, parameter ::  XG = 9.80665     ! gravity constant
      real, parameter :: XMD = 28.9644e-3  ! molecular weight of dry air
      real, parameter :: XMV = 18.0153e-3  ! molecular weight of water vapor
      real, parameter :: XRD = 287.05967   ! gaz constant for dry air
      real, parameter :: XRV = 461.524993  ! gaz constant for water vapor
      real, parameter :: XCPD = 1004.708845 ! specific heat of dry air
      real, parameter :: XCPV = 1846.1      ! specific heat of water vapor
      real, parameter :: XRHOLW = 1000.       ! density of liquid water
      real, parameter :: XCL = 4218.       ! specific heat of liquid water
      real, parameter :: XCI = 2106.       ! specific heat of ice
      real, parameter :: XTT = 273.16      ! triple point temperature
      real, parameter :: XLVTT = 2.5008e6    ! latent heat of vaporisation at XTT
      real, parameter :: XLSTT = 2.8345e6    ! latent heat of sublimation at XTT
      real, parameter :: XLMTT = 0.3337e6    ! latent heat of melting at XTT
      real, parameter :: XESTT = 611.14      ! saturation pressure at XTT
      real, parameter :: XALPW = 60.22416    ! constants in saturation pressure over liquid water
      real, parameter :: XBETAW = 6822.459384
      real, parameter :: XGAMW = 5.13948
      real, parameter :: XALPI = 32.62116    ! constants in saturation pressure over ice
      real, parameter :: XBETAI = 6295.421
      real, parameter :: XGAMI = 0.56313
      logical, parameter :: LUSERI = .true. ! logical switch to compute both
      ! liquid and solid condensate (LUSERI=.TRUE.)
      ! or only liquid condensate (LUSERI=.FALSE.)
      !
      !*       0.1   Declarations of dummy arguments :
      !
      !
      integer, intent(in)   :: mzp     ! vertical dimension
      integer, intent(in)   :: kts     ! vertical  computations start at
      !                                             ! KTS that is at least 1
      integer, intent(in)   :: ktf     ! vertical computations can be
      ! limited to MZP + 1 - KTF
      ! default=1
      real, dimension(mzp), intent(in)    :: PPABS  ! pressure (Pa)
      real, dimension(mzp), intent(in)    :: PZZ    ! height of model levels (m)
      real, dimension(mzp), intent(in)    :: PT     ! grid scale T  (K)
      real, dimension(mzp), intent(in)    :: PRV    ! grid scale water vapor mixing ratio (kg/kg)
      real, dimension(mzp), intent(in)    :: PMFLX  ! convective mass flux (kg/(s m^2))
      real, dimension(mzp), intent(in)    :: QRCO   ! sub-grid scale liq water mixing ratio (kg/kg)
      real, dimension(mzp), intent(in)    :: QCO    ! in-cloud water mixing ratio (kg/kg)
      real, dimension(mzp), intent(inout), optional :: PRC    ! grid scale r_c mixing ratio (kg/kg)
      real, dimension(mzp), intent(inout), optional :: PRI    ! grid scale r_i (kg/kg)
      real, dimension(mzp), intent(out)   :: PCLDFR ! fractional cloudiness (between 0 and 1)
      !
      !
      !*       0.2   Declarations of local variables :
      !
      integer  ::  JKT, JKP, JKM, K     ! loop index
      real, dimension(mzp) :: ZTLK, ZRT       ! work arrays for T_l, r_t
      real, dimension(mzp) :: ZL              ! length-scale
      integer   :: ITPL    ! top levels of tropopause/highest inversion
      real   :: ZTMIN   ! min Temp. related to ITPL
      real, dimension(mzp) :: LOC_PRC, LOC_PRI
      !
      real :: ZTEMP, ZLV, ZLS, ZTL, ZPV, ZQSL, ZPIV, ZQSI, ZFRAC, ZCOND, ZCPD ! thermodynamics
      real :: ZLL, DZZ, ZZZ ! length scales
      real :: ZAH, ZA, ZB, ZSBAR, ZQ1, ZSIGMA, ZDRW, ZDTL ! related to computation of Sig_s
      real :: ZSIG_CONV, ZSIGMA_NOCONV, ZQ1_NOCONV      ! convective part of Sig_s
      !
      !*       0.3  Definition of constants :
      !
      !-------------------------------------------------------------------------------
      !
      real :: ZL0 = 600.        ! tropospheric length scale
      ! changed to 600 m instead of 900 m to give a consistent
      ! value (linear increase) in general 500 m deep oceanic
      ! mixed layer - but could be put back to 900 m if wished
      real :: ZCSIGMA = 0.2         ! constant in sigma_s parameterization
      real :: ZCSIG_CONV = 0.30e-2  ! scaling factor for ZSIG_CONV as function of mass flux
      !
      !
      logical :: ONLY_CONVECTIVE_CLOUD_FRACTION = .true. ! set .false. for the total cloud fraction
      !-------------------------------------------------------------------------------
      !RETURN
      !
      if (present(PRC)) then
         LOC_PRC(:) = PRC(:)
      else
         LOC_PRC(:) = 0.0
      end if
      if (present(PRI)) then
         LOC_PRI(:) = PRI(:)
      else
         LOC_PRI(:) = 0.0
      end if

      PCLDFR(:) = 0. ! Initialize values
      !

      JKT = MZP + 1 - KTS
      !-will limit the model vertical column to 60 hPa
      do K = KTF, KTS, -1
         if (PPABS(k) > 60.*100.) then
            JKT = k
            !PRINT*,"JKT=",K,MZP+1-KTS ;CALL FLUSH(6)
            exit
         end if
      end do

      do K = KTS, JKT
         ZTEMP = PT(k)
         !latent heat of vaporisation/sublimation
         ZLV = XLVTT + (XCPV - XCL)*(ZTEMP - XTT)
         ZLS = XLSTT + (XCPV - XCI)*(ZTEMP - XTT)

         !store temperature at saturation and total water mixing ratio
         ZRT(k) = PRV(k) + LOC_PRC(k) + LOC_PRI(k)
         ZCPD = XCPD + XCPV*PRV(k) + XCL*LOC_PRC(k) + XCI*LOC_PRI(k)
         ZTLK(k) = ZTEMP - ZLV*LOC_PRC(k)/ZCPD - ZLS*LOC_PRI(k)/ZCPD
      end do

      !-------------------------------------------------------------------------------
      ! Determine tropopause/inversion  height from minimum temperature

      ITPL = KTS + 1
      ZTMIN = 400.
      do k = KTS + 1, JKT - 1
         if (PT(k) < ZTMIN) then
            ZTMIN = PT(k)
            ITPL = K
         end if
      end do

      ! Set the mixing length scale - used for computing the "turbulent part" of Sigma_s

      ZL(:) = 20.
      do k = KTS + 1, JKT

         ! free troposphere
         ZL(k) = ZL0
         JKP = ITPL
         ZZZ = PZZ(k) - PZZ(KTS)
         ! approximate length for boundary-layer : linear increase
         if (ZL0 > ZZZ) ZL(k) = ZZZ
         ! gradual decrease of length-scale near and above tropopause/top inversion
         if (ZZZ > 0.9*(PZZ(JKP) - PZZ(KTS))) &
            ZL(k) = .6*ZL(K - 1)
      end do
      !-------------------------------------------------------------------------------

      do k = KTS + 1, JKT - 1
         JKP = k + 1
         JKM = k - 1
         ZTEMP = PT(k)

         !latent heat of vaporisation/sublimation
         ZLV = XLVTT + (XCPV - XCL)*(ZTEMP - XTT)
         ZLS = XLSTT + (XCPV - XCI)*(ZTEMP - XTT)

         ZCPD = XCPD + XCPV*PRV(k) + XCL*LOC_PRC(k) + XCI*LOC_PRI(k)
         !temperature at saturation
         ZTL = ZTEMP - ZLV*LOC_PRC(k)/ZCPD - ZLS*LOC_PRI(k)/ZCPD

         !saturated water vapor mixing ratio over liquid water
         ZPV = min(exp(XALPW - XBETAW/ZTL - XGAMW*log(ZTL)), 0.99*PPABS(k))
         ZQSL = XRD/XRV*ZPV/(PPABS(k) - ZPV)

         !saturated water vapor mixing ratio over ice
         ZPIV = min(exp(XALPI - XBETAI/ZTL - XGAMI*log(ZTL)), 0.99*PPABS(k))
         ZQSI = XRD/XRV*ZPIV/(PPABS(k) - ZPIV)

         !interpolate between liquid and solid as function of temperature
         ! glaciation interval is specified here to 20 K
         ZFRAC = (ZTL - 250.16)/(XTT - 250.16)  ! liquid/solid fraction
         ZFRAC = max(0., min(1., ZFRAC))

         if (.not. LUSERI) ZFRAC = 1.
         ZQSL = (1.-ZFRAC)*ZQSI + ZFRAC*ZQSL
         ZLV = (1.-ZFRAC)*ZLS + ZFRAC*ZLV

         !coefficients a and b
         ZAH = ZLV*ZQSL/(XRV*ZTL**2)*(XRV*ZQSL/XRD + 1.)
         !orig  ZAH  = ZLV * ZQSL / ( XRV * ZTL**2 )

         ZA = 1./(1.+ZLV/ZCPD*ZAH)
         ZB = ZAH*ZA

         !- parameterize Sigma_s with first_order closure
         DZZ = PZZ(JKP) - PZZ(JKM)
         ZDRW = ZRT(JKP) - ZRT(JKM)
         ZDTL = ZTLK(JKP) - ZTLK(JKM) + XG/ZCPD*DZZ
         ZLL = ZL(k)

         !- standard deviation due to convection
         ZSIG_CONV = ZCSIG_CONV*PMFLX(k)/ZA

         !- turb + conv
         ZSIGMA = sqrt(max(1.e-25, ZCSIGMA*ZCSIGMA*ZLL*ZLL/(DZZ*DZZ)*( &
                           ZA*ZA*ZDRW*ZDRW - 2.*ZA*ZB*ZDRW*ZDTL &
                           + ZB*ZB*ZDTL*ZDTL) &
                           + ZSIG_CONV*ZSIG_CONV))

         !- zsigma should be of order 4.e-4 in lowest 5 km of atmosphere
         ZSIGMA = max(ZSIGMA, 1.e-10)

         !- normalized saturation deficit
         ZSBAR = ZA*(ZRT(k) - ZQSL)
         !- "Q1" parameter
         ZQ1 = ZSBAR/ZSIGMA

         !- total cloud fraction
         PCLDFR(k) = max(0., min(1., 0.5 + 0.36*atan(1.55*ZQ1)))

         if (ONLY_CONVECTIVE_CLOUD_FRACTION) then
            !- get cloud fraction associated with ONLY the sub-grid scale convective part
            !- this sigma does not include the sub-grid scale convective part
            ZSIGMA_NOCONV = sqrt(max(1.e-25, ZCSIGMA*ZCSIGMA*ZLL*ZLL/(DZZ*DZZ)*( &
                                     ZA*ZA*ZDRW*ZDRW - 2.*ZA*ZB*ZDRW*ZDTL &
                                     + ZB*ZB*ZDTL*ZDTL)))
            !- zsigma should be of order 4.e-4 in lowest 5 km of atmosphere
            ZSIGMA_NOCONV = max(ZSIGMA_NOCONV, 1.e-10)
            ZQ1_NOCONV = ZSBAR/ZSIGMA_NOCONV

            !- cloud fraction associated with ONLY convective part ("total-turb")
            PCLDFR(k) = 0.36*(atan(1.55*ZQ1) - atan(1.55*ZQ1_NOCONV))

            PCLDFR(k) = max(0., min(1., PCLDFR(k)))

         end if
         !- newer formulation, see GMD 2015
         !PCLDFR(k) = MAX( 0., MIN(1.,0.5+0.34*ATAN(1.85*ZQ1+2.33)) )
         !- this is area fraction of cloud cores
         !PCLDFR(k) = MAX( 0., MIN(1.,0.292/ZQ1**2) )

         cycle
         !total condensate diagnostic (not being used)
         if (ZQ1 > 0. .and. ZQ1 <= 2.) then
            !orig   ZCOND =     EXP(-1.)+.66*ZQ1+.086*ZQ1*ZQ1
            ZCOND = min(exp(-1.) + .66*ZQ1 + .086*ZQ1**2, 2.) ! We use the MIN function for continuity
         else if (ZQ1 > 2.) then
            ZCOND = ZQ1
         else
            ZCOND = exp(1.2*ZQ1 - 1.)
         end if
         ZCOND = ZCOND*ZSIGMA

         if (zcond < 1.e-12) then
            zcond = 0.
            pcldfr(k) = 0.
         end if
         if (pcldfr(k) == 0.) then
            zcond = 0.
         end if

         LOC_PRC(k) = ZFRAC*ZCOND ! liquid condensate
         if (LUSERI) then
            LOC_PRI(k) = (1.-ZFRAC)*ZCOND   ! solid condensate
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
   function auto_rk(n, step, aux, xexp, qrc1) result(PW)
      integer, intent(in) :: n
      real, intent(in) :: step, aux, qrc1, xexp
      real             :: PW

      PW = step*qrc1*(1.0 - exp(-aux**xexp))/float(n)

   end function auto_rk

   !---------------------------------------------------------------------------------------------------
   subroutine interface_aerchem(mtp, itrcr, aer_chem_mech, cnames, qnames, fscav, fscav_int)
      implicit none
      integer, intent(in) :: mtp, ITRCR
      character(len=*), intent(in)   :: AER_CHEM_MECH
      character(len=*), dimension(mtp), intent(in)   :: CNAMES, QNAMES
      real, dimension(ITRCR), intent(in) :: FSCAV

      real, dimension(mtp), intent(out)   :: FSCAV_INT
      !-local vars
      integer :: ispc, len_ACM, len_spc, irun = 0
      character(len=100) :: TMP_AER_NAME
      ispc_co = 1
      chem_name_mask = 1
      chem_name_mask_evap = 1
      chem_adj_autoc = 1.0

      !- GOCART + PCHEM section
      if (AER_CHEM_MECH == 'GOCART') then
         len_ACM = len(trim(AER_CHEM_MECH))
         do ispc = 1, mtp
            FSCAV_INT(ispc) = fscav(ispc)
            len_spc = len(trim(cnames(ispc)))
            chem_name(ispc) = trim(qnames(ispc))!(len_ACM+1:len_spc))
            if (trim(chem_name(ispc)) == 'CO') ispc_co = ispc

            !-the tracers below are not being transported :
            if (trim(chem_name(ispc)) == "NCPI") chem_name_mask(ispc) = 0
            if (trim(chem_name(ispc)) == "NCPL") chem_name_mask(ispc) = 0
            if (trim(chem_name(ispc)) == "QGRAUPEL") chem_name_mask(ispc) = 0
            if (trim(chem_name(ispc)) == "QRAIN") chem_name_mask(ispc) = 0
            if (trim(chem_name(ispc)) == "QSNOW") chem_name_mask(ispc) = 0
            !if(trim(chem_name(ispc)) == "QW"      ) CHEM_NAME_MASK     (ispc) = 0
            !if(trim(chem_name(ispc)) == "QW"      ) CHEM_NAME_MASK_EVAP(ispc) = 0

            if (trim(chem_name(ispc) (1:3)) == "du0") chem_adj_autoc(ispc) = 1.0
            if (trim(chem_name(ispc) (1:3)) == "ss0") chem_adj_autoc(ispc) = 1.0
            if (trim(chem_name(ispc) (1:3)) == "OCp") chem_adj_autoc(ispc) = 1.0
            if (trim(chem_name(ispc) (1:3)) == "BCp") chem_adj_autoc(ispc) = 1.0
            if (trim(chem_name(ispc) (1:3)) == "SO4") chem_adj_autoc(ispc) = 1.0
         end do
         !-----------------------------------------------------------------------------------------
         !---temporary section to fill Henrys cts for N2O and CH4 of PCHEM chemical mechanism
         do ispc = 1, mtp
            if (trim(chem_name(ispc)) == 'N2O' .or. &
                trim(chem_name(ispc)) == 'CH4' &
                ) then

               call getHenryLawCts(trim(chem_name(ispc)), &
                                     hcts(ispc)%hstar, hcts(ispc)%dhr, hcts(ispc)%ak0, hcts(ispc)%dak)
            end if
         end do
         !***   all tracers not having wet deposition _must_ have all Hcts reset to zero here.
         !.. WHERE( Hcts(:)%hstar == -1.) Hcts(:)%hstar = 0.
         !.. WHERE( Hcts(:)%dhr   == -1.) Hcts(:)%dhr   = 0.
         !.. WHERE( Hcts(:)%ak0   == -1.) Hcts(:)%ak0   = 0.
         !.. WHERE( Hcts(:)%dak   == -1.) Hcts(:)%dak   = 0.
         !-----------------------------------------------------------------------------------------
      end if
      return
      !IF( MAPL_AM_I_ROOT() .and. irun == 0)THEN
      irun = 1
      print *, "========================================================================="; call flush (6)
      print *, " the following tracers will be transport by GF scheme                    "; call flush (6)

      write (10, *) "================= table of tracers in the GF conv transport ==================="
      write (10, *) " SPC,  CHEM_NAME,  FSCAV          - the four Henrys law cts  -   Transport Flag - kc adjust"
      do ispc = 1, mtp
         write (10, 121) ispc, trim(chem_name(ispc)), FSCAV_INT(ispc), hcts(ispc)%hstar &
            , hcts(ispc)%dhr, hcts(ispc)%ak0, hcts(ispc)%dak, chem_name_mask(ispc), chem_adj_autoc(ispc)
         if (chem_name_mask(ispc) == 1) then
            print *, "GF is doing transp and wet removal of: ", trim(chem_name(ispc))
            call flush (6)
         end if
      end do
      print *, "========================================================================="; call flush (6)
121   format(1x, I4, A10, 5f14.5, I4, F14.5)
      !ENDIF
   end subroutine interface_aerchem

      !---------------------------------------------------------------------------------------------------
   subroutine bidiag(m, b, c, f)
      !-- this routine solves the problem:  bb*f(k,t+1) + cc*f(k+1,t+1) = dd
      !-- an updated "f" at time t+1 is the output
      implicit none
      integer, intent(in) :: m
      real, dimension(m), intent(inout) :: b, c
      real, dimension(m), intent(inout) :: f
      !--locals
      real, dimension(m) :: q
      integer :: k
      real :: p

      c(m) = 0.
      q(1) = -c(1)/b(1)
      f(1) = f(1)/b(1)
      do k = 2, m
         p = 1./b(k)
         q(k) = -c(k)*p
         f(k) = f(k)*p
      end do
      do k = m - 1, 1, -1
         f(k) = f(k) + q(k)*f(k + 1)
      end do

   end subroutine bidiag


   !------------------------------------------------------------------------------------
   subroutine cup_up_rain(cumulus, klcl, kbcon, ktop, k22, ierr, xland &
                          , zo_cup, qco, qrco, pwo, pwavo, po, p_cup, t_cup, tempco &
                          , zuo, up_massentr, up_massdetr, vvel2d, rho &
                          , qrr &
                          , itf, ktf, its, ite, kts, kte)

      implicit none
      character*(*), intent(in) :: cumulus
      integer, intent(in) :: itf, ktf, its, ite, kts, kte
      integer, dimension(its:ite), intent(in) :: kbcon, ktop, k22, klcl, ierr
      real, dimension(its:ite), intent(in) :: xland, pwavo
      real, dimension(its:ite, kts:kte), intent(in)  :: &
         zo_cup, qco, qrco, pwo, po, p_cup, t_cup, zuo &
         , up_massentr, up_massdetr, vvel2d, tempco, rho

      !--for future use (rain water mixing ratio)
      real, dimension(its:ite, kts:kte), intent(out) :: qrr      !-- units kg[water]/kg[air]

      !-- locals
      integer :: i, k
      real :: tmp
      integer, dimension(its:ite) :: start_level
      real :: dz, z1, zrnew, zc, zd, zint, z2, zrold, denom, fall_fact, wup, exp1, R_vr
      real, dimension(kts:kte) :: prec_flx_rain, prec_flx_snow
      real, dimension(kts:kte) :: pw, p_liq_ice ! - rainfall source
      real, parameter :: rho1000mb = 1.2, rhow = 1000., N_r = 0.1 ! cm^-3, rainfall drops number concen
      real, parameter :: exp_KR = 1./5. &! Kuo & Raymond 1983
                         , exp_SB = 2./3.  ! Seifert & Beheng 2006 eq 27

      qrr = 0.
      if (c0 < 1.e-6) return

      !--- rain water mixing ratio
      do i = its, itf
         if (ierr(i) /= 0) cycle

         do k = ktop(i), kts, -1

            p_liq_ice(k) = fract_liq_f(tempco(i, k))

            !--- transport + mixing
            denom = zuo(i, k + 1) - .5*up_massdetr(i, k) + up_massentr(i, k)
            if (denom > 0.) then
               qrr(i, k) = (qrr(i, k + 1)*zuo(i, k + 1) - .5*up_massdetr(i, k)*qrr(i, k + 1))/denom
            else
               qrr(i, k) = qrr(i, k + 1)
            end if

            !--- rain source
            pw(k) = pwo(i, k)/(1.e-16 + zuo(i, k))

            !-- rainfall sedimentation
            !-- Kuo & Raymond 1983
            !-- fallout of rain (21.18 * qrr^0.2 have m/s as units with qrr in kg/kg)
            !---                                half velocity for ice phase
            fall_fact = 21.18*(p_liq_ice(k) + 0.5*(1.-p_liq_ice(k)))
            exp1 = exp_KR
            !
            !-- Seifert & Beheng 2006 eq 27, units= m/s
            ! fall_fact = 159.*sqrt(rho1000mb/rho(i,k))*( p_liq_ice(k) + 0.5*(1.-p_liq_ice(k) ))
            ! exp1      = exp_SB

            !-- Kogan 2013
            R_vr = (4.*c_pi*rhow/(3.*rho(i, k)))**(-1./3.)*(qrr(i, k) + pw(k))**(1./3.)* &
                   (N_r)**(-1./3.)! cm, N_r is the rainfall drops number concentration, and not CCN or N_c

            R_vr = max(40., R_vr*1.e-2*1.e+6)   ! micrometer
            fall_fact = 1.e-2*(2.4*R_vr - 62.0)        ! m/s
            exp1 = 0.

            wup = min(15., max(2., vvel2d(i, k)))
            z2 = fall_fact/wup

            !--exact solution
            if (qrr(i, k) + pw(k) > 0.) then
               !-- this is the sedimentation speed divided by W_up
               zd = z2*(qrr(i, k) + pw(k))**exp1
               zint = exp(-zd)
               zrold = qrr(i, k)
               zc = pw(k)
               zrnew = zrold*zint + zc/zd*(1.-zint)
               zrnew = max(0.0, min(qrr(i, k) + pw(k), zrnew))
            else
               zrnew = 0.
            end if
            qrr(i, k) = zrnew

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
         end do
      end do

   end subroutine cup_up_rain


   !------------------------------------------------------------------------------------
   subroutine get_condensation(q_old, t_old, po_cup, q_new, t_new)

      !-- calculate condensation and adjust t and q accordingly
      implicit none
      real, intent(in) :: po_cup, q_old, t_old ! before condensation
      real, intent(inout) ::        q_new, t_new ! after  condensation

      !---locals
      real ::  zqp, zcond, zcond1, zcor, zqsat, zi, zl, zf
      real :: psp, pt, pq
      real :: z3es, z4es, z5alcp, zaldcp
      real :: ptare, cond
      real :: foeewmcu, foealfcu, foedemcu, foeldcpmcu

      real, parameter :: &
         RD = 287.06 &
         , RV = 461.52 &
         , RCPD = 1004.71 &
         , RTT = 273.16 &
         , RHOH2O = 1000. &
         , RLVTT = 2.5008e+6 &
         , RLSTT = 2.8345e+6 &
         , RETV = RV/RD - 1.0 &
         , RLMLT = RLSTT - RLVTT &
         , RCPV = 4.*RV &
         , R2ES = 611.21*RD/RV &
         , R3LES = 17.502 &
         , R3IES = 22.587 &
         , R4LES = 32.19 &
         , R4IES = -0.7 &
         , R5LES = R3LES*(RTT - R4LES) &
         , R5IES = R3IES*(RTT - R4IES) &
         , R5ALVCP = R5LES*RLVTT/RCPD &
         , R5ALSCP = R5IES*RLSTT/RCPD &
         , RALVDCP = RLVTT/RCPD &
         , RALSDCP = RLSTT/RCPD &
         , RALFDCP = RLMLT/RCPD &
         , RTWAT = RTT &
         , RTBER = RTT - 5. &
         , RTBERCU = RTT - 5.0 &
         , RTICE = RTT - 23. &
         , RTICECU = RTT - 23. &
         , RTWAT_RTICE_R = 1./(RTWAT - RTICE) &
         , RTWAT_RTICECU_R = 1./(RTWAT - RTICECU) &
         , RVTMP2 = RCPV/RCPD - 1. &
         , ZQMAX = 0.5

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
      PT = t_old       ! K
      PQ = q_old       ! kg/kg
      PSP = po_cup*100. ! hPa

      !--- get condensation in moist ascent --------------------------
      PTARE = PT
      ZQP = 1.0/PSP

      ZL = 1.0/(PT - R4LES)
      ZI = 1.0/(PT - R4IES)

      FOEALFCU = min(1.0, ((max(RTICECU, min(RTWAT, PTARE)) - RTICECU)*RTWAT_RTICECU_R)**2)
      ZQSAT = R2ES*(FOEALFCU*exp(R3LES*(PTARE - RTT)*ZL) + &
                    (1.0 - FOEALFCU)*exp(R3IES*(PTARE - RTT)*ZI))

      ZQSAT = ZQSAT*ZQP
      ZQSAT = min(0.5, ZQSAT)
      ZCOR = 1.0 - RETV*ZQSAT

      ZF = FOEALFCU*R5ALVCP*ZL**2 + (1.0 - FOEALFCU)*R5ALSCP*ZI**2
      ZCOND = (PQ*ZCOR**2 - ZQSAT*ZCOR)/(ZCOR**2 + ZQSAT*ZF)

      if (ZCOND > 0.0) then
         FOELDCPMCU = FOEALFCU*RALVDCP + (1.0 - FOEALFCU)*RALSDCP
         PT = PT + FOELDCPMCU*ZCOND
         PTARE = PT
         PQ = PQ - ZCOND

         ZL = 1.0/(PT - R4LES)
         ZI = 1.0/(PT - R4IES)

         FOEALFCU = min(1.0, ((max(RTICECU, min(RTWAT, PTARE)) - RTICECU)*RTWAT_RTICECU_R)**2)
         ZQSAT = R2ES*(FOEALFCU*exp(R3LES*(PT - RTT)*ZL) + &
                       (1.0 - FOEALFCU)*exp(R3IES*(PT - RTT)*ZI))

         ZQSAT = ZQSAT*ZQP
         ZQSAT = ZQSAT - 0.5*(abs(0.5 - ZQSAT) - (0.5 - ZQSAT))

         ZCOR = 1.0 - RETV*ZQSAT
         ZF = FOEALFCU*R5ALVCP*ZL**2 + (1.0 - FOEALFCU)*R5ALSCP*ZI**2

         ZCOND1 = (PQ*ZCOR**2 - ZQSAT*ZCOR)/(ZCOR**2 + ZQSAT*ZF)
         if (ZCOND == 0.0) ZCOND1 = 0.0
         FOELDCPMCU = FOEALFCU*RALVDCP + (1.0 - FOEALFCU)*RALSDCP
         PT = PT + FOELDCPMCU*ZCOND1
         PQ = PQ - ZCOND1
      end if

      !-- FINAL --------------------------
      q_new = PQ
      t_new = PT
      cond = -ZCOND1 != q_old-qnew, source for the liquid water
   end subroutine get_condensation

   !---------------------------------------------------------------------------------
   elemental real function make_RainNumber(Q_rain, temp)

      implicit none

      real, intent(in):: q_rain, temp
      double precision:: lambda, n0, qnr
      real, parameter:: am_r = c_pi*1000./6.

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
      elseif (temp .gt. 271.15 .and. temp .lt. 273.15) then
         N0 = 8.*10**(279.15 - temp)
      end if

      lambda = sqrt(sqrt(N0*am_r*6.0/Q_rain))
      qnr = Q_rain/6.0*lambda*lambda*lambda/am_r
      make_RainNumber = SNGL(qnr)

      return
   end function make_RainNumber
