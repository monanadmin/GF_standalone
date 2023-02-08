module modConstants
   !! ## Module with constants to be used by MONAN Sources
   !!
   !! ![](https://i.ibb.co/LNqGy3S/logo-Monan-Color-75x75.png)
   !! ## MONAN
   !!
   !! Author: Rodrigues, L. F. [LFR]
   !!
   !! E-mail: <mailto:luiz.rodrigues@inpe.br>
   !!
   !! Date: 28 September 2018 09:00
   !!
   !! #####Version: 0.1.0
   !!
   !! ---
   !! **Full description**:
   !!
   !! Module with constants to be used by MONAN Sources
   !!
   !! ** History**:
   !!
   !! --- 
   !! ** Licence **:
   !!
   !!  <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
   !!
   !!  This program is free software: you can redistribute it and/or modify
   !!  it under the terms of the GNU General Public License as published by
   !!  the  Free  Software  Foundation, either version 3 of the License, or
   !!  (at your option) any later version.
   !!
   !!  This program is distributed in the hope that it  will be useful, but
   !!  ** WITHOUT  ANY  WARRANTY **;  without  even  the   implied   warranty  of
   !!  **MERCHANTABILITY** or **FITNESS FOR A  PARTICULAR PURPOSE**.  See  the, GNU
   !!  GNU General Public License for more details.
   !!
   !!  You should have received a copy  of the GNU General  Public  License
   !!  along with this program.  If not, see [GNU Public License](https://www.gnu.org/licenses/gpl-3.0.html).
   !!
   
   implicit none
   character(len=*), parameter :: sourceName = 'modConstants.F90' 
   !! Source code name 
   character(len=*), parameter :: moduleName = 'modConstants' 
   !! module name 

   !integer kinds
   integer, parameter :: kind_ib = selected_int_kind(13)
   !! 8 byte integer
   integer, parameter :: kind_im = selected_int_kind(6)
   !! 4 byte integer
   integer, parameter :: kind_in = kind(1)
   !! native integer

   !real kinds
   integer, parameter :: kind_rb = selected_real_kind(12)
   !! 8 byte real
   integer, parameter :: kind_rm = selected_real_kind(6)
   !! 4 byte real
   integer, parameter :: kind_rn = kind(1.0)
   !! native real

   character(len=*), parameter :: c_modelVersion='Rev. 0.1.0'
   character(len=*), parameter :: modelVersion='Rev. 0.1.0'
   character(len=*), parameter :: c_license='CC Attribution-ShareAlike 4.0 International'
   !! Last version of model

   !! Phys & others constants
   real(kind=kind_rb), parameter :: c_pi       = 3.1415926535897932384626433
   real(kind=kind_rb), parameter :: c_rgas     = 287.
   !! J K-1 kg-1
   real(kind=kind_rb), parameter :: c_rgas_ = 8.205e-2 
   !! atm M^-1 K^-1 ! 8.314 gas constant [J/(mol*K)]
   real(kind=kind_rb), parameter :: c_avogad = 6.022e23
   !! Avogadro constant [1/mol]
   real(kind=kind_rb), parameter:: c_rhoH2O = 999.9668
   !! density of water [kg/m3]
   real(kind=kind_rb), parameter:: c_hplus = 1.175e-4     
   !!  for cloud water. pH is asuumed to be 3.93: pH=3.93 =>hplus=10**(-pH)
   real, parameter:: c_temp0 =    298.15
   !! standard temperature [K]
   real, parameter:: c_temp0i= 1./c_temp0
   !! inverse of standard temperature [K]
   real(kind=kind_rb), parameter :: c_cp       = 1004.
   !! J K-1 kg-1
   real(kind=kind_rb), parameter :: c_rv = 461.
   !! J K-1 kg-1
   real(kind=kind_rb), parameter :: c_cv       = 717.
   real(kind=kind_rb), parameter :: c_rm       = 461.
   real(kind=kind_rb), parameter :: c_p00      = 1.e5
   !! hPa
   real(kind=kind_rb), parameter :: c_t00      = 273.16
   !Absolute temperature
   real(kind=kind_rb), parameter :: c_tcrit = 258.
   !! K
   real(kind=kind_rb), parameter :: c_xlv = 2.5e6
   !! J kg-1
   real(kind=kind_rb), parameter :: c_akmin = 1.0
   !! #
   real(kind=kind_rb), parameter :: c_pi180    = c_pi / 180.
   real(kind=kind_rb), parameter :: c_i_pi180  = 1./c_pi180
   real(kind=kind_rb), parameter :: c_pi4      = c_pi * 4.
   real(kind=kind_rb), parameter :: c_spcon    = 111120.
   real(kind=kind_rb), parameter :: c_erad     = 6367000.
   real(kind=kind_rb), parameter :: c_arc      = c_erad*c_pi180
   real(kind=kind_rb), parameter :: c_vonk     = 0.40
   real(kind=kind_rb), parameter :: c_tkmin    = 5.e-4
   !! Minimum TKE [J/kg]
   real(kind=kind_rb), parameter :: c_tkmin2 = 1.e-5
   !! m+2 s-2
   real(kind=kind_rb), parameter :: c_grav      = 9.80665
   !! Gravity acceleration [m/s]
   real(kind=kind_rb), parameter :: c_ccnclean = 250.
   !! # cm-3
   real(kind=kind_rb), parameter :: c_t_0 = 273.16
   !! K
   real(kind=kind_rb), parameter :: c_t_ice = 235.16
   !! K
   real(kind=kind_rb), parameter :: c_xlf = 0.333e6
   !! latent heat of freezing (J kg-1)
   real(kind=kind_rb), parameter :: c_max_qsat = 0.5
   !! kg/kg
   real(kind=kind_rb), parameter :: c_alvl     = 2.50e6
   real(kind=kind_rb), parameter :: c_alvi     = 2.834e6
   real(kind=kind_rb), parameter :: c_alli     = 0.334e6
   real(kind=kind_rb), parameter :: c_alvl2    = 6.25e12
   real(kind=kind_rb), parameter :: c_alvi2    = 8.032e12
   real(kind=kind_rb), parameter :: c_solar    = 1.3533e3
   real(kind=kind_rb), parameter :: c_stefan   = 5.6696e-8
   real(kind=kind_rb), parameter :: c_cww      = 4218.
   real(kind=kind_rb), parameter :: c_c0       = 752.55 * 4.18684e4
   real(kind=kind_rb), parameter :: c_viscos   = .15e-4
   real(kind=kind_rb), parameter :: c_rowt     = 1.e3
   real(kind=kind_rb), parameter :: c_dlat     = 111120.
   real(kind=kind_rb), parameter :: c_omega    = 7.292e-5
   real(kind=kind_rb), parameter :: c_rocp     = c_rgas / c_cp
   real(kind=kind_rb), parameter :: c_p00i     = 1. / c_p00
   real(kind=kind_rb), parameter :: c_cpor     = c_cp / c_rgas
   real(kind=kind_rb), parameter :: c_rocv     = c_rgas / c_cv
   real(kind=kind_rb), parameter :: c_cpi      = 1. / c_cp
   real(kind=kind_rb), parameter :: c_cpi4     = 4. * c_cpi
   real(kind=kind_rb), parameter :: c_cp253i   = c_cpi / 253.
   real(kind=kind_rb), parameter :: c_allii    = 1. / c_alli
   real(kind=kind_rb), parameter :: c_aklv     = c_alvl / c_cp
   real(kind=kind_rb), parameter :: c_akiv     = c_alvi / c_cp
   real(kind=kind_rb), parameter :: c_gama     = c_cp / c_cv
   real(kind=kind_rb), parameter :: c_gg       = .5 * c_grav
   real(kind=kind_rb), parameter :: c_ep       = c_rgas / c_rm
   real(kind=kind_rb), parameter :: c_airmw    = 28.965
   !! effective molecular mass of air [kg/kmol]
   real(kind=kind_rb), parameter :: c_h2omw = 18.015 
   !! molecular mass of water [kg/kmol]
   real(kind=kind_rb), parameter :: c_runiv = 8314.47
   !! J/(Kmole K)
   real(kind=kind_rb), parameter :: c_rdry = c_runiv/c_airmw  
   !! j/(kg k) 
   real(kind=kind_rb), parameter :: c_cpdry = 3.5*c_rdry          !
   !! j/(kg k) 
   real(kind=kind_rb), parameter :: c_kappa = c_rdry/c_cpdry   !
   !! (2.0/7.0) 
   real(kind=kind_rb), parameter :: c_epsilon = c_h2omw/c_airmw  
   !! -- 
   real(kind=kind_rb), parameter :: c_rgas2 = c_rdry              
   !! j/(kg k) (deprecated) 
   real(kind=kind_rb), parameter :: c_cp2 = c_rgas/c_kappa  
   !! j/(kg k) (deprecated) 
   real(kind=kind_rb), parameter :: c_vireps = 1.0/c_epsilon - 1.0   
   !!          (deprecated)
   real(kind=kind_rb), parameter :: c_p00k     = 26.870941
   !!  = p00 ** rocp
   real(kind=kind_rb), parameter :: c_p00ki    = 1. / c_p00k

   real(kind=kind_rb), parameter :: c_onethird  = 1./3.
   !! 1/3
   real(kind=kind_rb), parameter :: c_half  = 1./2.
   !! 1/2

   real(kind=kind_rb), parameter :: c_retv = c_rm/c_rgas - 1.0 
   real(kind=kind_rb), parameter :: c_rlmlt = c_alvi - c_alvl
   real(kind=kind_rb), parameter :: c_rcpv = 4.*c_rm 
   real(kind=kind_rb), parameter :: c_r2es = 611.21*c_ep
   real(kind=kind_rb), parameter :: c_r3les = 17.502
   real(kind=kind_rb), parameter :: c_r3ies = 22.587
   real(kind=kind_rb), parameter :: c_r4les = 32.19
   real(kind=kind_rb), parameter :: c_r4ies = -0.7
   real(kind=kind_rb), parameter :: c_r5les = c_r3les*(c_t00 - c_r4les)
   real(kind=kind_rb), parameter :: c_r5ies = c_r3ies*(c_t00 - c_r4ies)
   real(kind=kind_rb), parameter :: c_r5alvcp = c_r5les*c_alvl/c_cp
   real(kind=kind_rb), parameter :: c_r5alscp = c_r5ies*c_alvi/c_cp
   real(kind=kind_rb), parameter :: c_ralvdcp = c_alvl/c_cp
   real(kind=kind_rb), parameter :: c_ralsdcp = c_alvi/c_cp
   real(kind=kind_rb), parameter :: c_ralfdcp = c_rlmlt/c_cp
   real(kind=kind_rb), parameter :: c_rtber = c_t00 - 5.
   real(kind=kind_rb), parameter :: c_rtIce = c_t00 - 23.
   real(kind=kind_rb), parameter :: c_rtwat_rtIce_r = 1./(c_t00 - c_rtIce)
   real(kind=kind_rb), parameter :: c_rvtmp2 = c_rcpv/c_cp - 1.
   real(kind=kind_rb), parameter :: p_zqmax = 0.5


   ! Lower bounds for turbulence-related variables                                         !
   real(kind=kind_rb), parameter :: c_sigwmin     = 1.e-4
   !! Minimum sigma-w                     [m/s]
   real(kind=kind_rb), parameter :: c_abslmomin   = 1.e-4
   !! Minimum abs value of Obukhov length [m]
   real(kind=kind_rb), parameter :: c_ltscalemax  = 1.e5
   !! Maximum Lagrangian timescale        [s]
   real(kind=kind_rb), parameter :: c_abswltlmin  = 1.e-4
   !! Minimum abs value of Theta*         [K m/s]
   real(kind=kind_rb), parameter :: c_lturbmin    = 1.e-3
   !! Minimum abs value of turb. lenght   [m]

   !Conversion factors
   real, parameter:: p_m_2_cm = 100.   
   !!  [m]   to [cm]
   real, parameter :: p_m3_2_l = 1000.     
   !!  [m^3]      to [l]
   real, parameter:: p_l_2_m3 = 1/p_m3_2_l   
   !!  [l]    to [m^3]
   real, parameter:: p_pa_2_atm = 1./101325.  
   !!  [Pa]      to [atm]

   real,parameter :: c_scday=86400.0
   !! Seconds by day

   !Tiny numbers
   real(kind=kind_rb), parameter :: c_tinyReal = 1.17549435E-38
   !! Tiny real number
   !real(kind=kind_rb), parameter :: tinyDouble = 2.2250738585072014E-308
   !! Tiny double precision number

   !Errors type for dump functions
   integer, parameter :: p_noError=0
   !! No Error, just write
   integer, parameter :: p_notice=1
   !! Notice Error
   integer, parameter :: p_warning=2
   !! warning Error
   integer, parameter :: p_fatal=3
   !! Fatal error
   integer, parameter :: p_kill=4
   !! Kill the model Signal
   integer, parameter :: p_continue=5
   !! Continue run Signal

   logical, parameter :: p_yes=.true.
   !! Yes is the .true. value
   logical, parameter :: p_no=.false.
   !! No is the false value

   character(len=*), parameter :: p_empty=''
   !! empty string

   !Colors strings for print and write commands
   character(len=*), parameter :: p_darkGrey   =achar(27)//'[90m'
   character(len=*), parameter :: p_peach      =achar(27)//'[91m'
   character(len=*), parameter :: p_lightGreen =achar(27)//'[92m'
   character(len=*), parameter :: p_lightYellow=achar(27)//'[93m'
   character(len=*), parameter :: p_lightBlue  =achar(27)//'[94m'
   character(len=*), parameter :: p_pink       =achar(27)//'[95m'
   character(len=*), parameter :: p_lightAqua  =achar(27)//'[96m'
   character(len=*), parameter :: p_pearlWhite =achar(27)//'[97m'
   character(len=*), parameter :: p_black      =achar(27)//'[30m'
   character(len=*), parameter :: p_red        =achar(27)//'[31m'
   character(len=*), parameter :: p_green      =achar(27)//'[32m'
   character(len=*), parameter :: p_yellow     =achar(27)//'[33m'
   character(len=*), parameter :: p_blue       =achar(27)//'[34m'
   character(len=*), parameter :: p_purple     =achar(27)//'[35m'
   character(len=*), parameter :: p_aqua       =achar(27)//'[36m'
   character(len=*), parameter :: p_blink      =achar(27)//'[31;5;95;38;5;214m'
   character(len=*), parameter :: p_noColor    =achar(27)//'[0m'
   character(len=*), parameter :: p_Iblack     =achar(27)//'[40m'
   character(len=*), parameter :: p_Ired       =achar(27)//'[41m'
   character(len=*), parameter :: p_Igreen     =achar(27)//'[42m'
   character(len=*), parameter :: p_Iyellow    =achar(27)//'[43m'
   character(len=*), parameter :: p_Iblue      =achar(27)//'[44m'
   character(len=*), parameter :: p_IMagenta   =achar(27)//'[45m'
   character(len=*), parameter :: p_Icyan      =achar(27)//'[46m'
   character(len=*), parameter :: p_underline  =achar(27)//'[4m'
   character(len=*), parameter :: p_strike     =achar(27)//'[9m'
   character(len=*), parameter :: p_bold       =achar(27)//'[1m'
   character(len=*), parameter :: p_normal     =achar(27)//'[0m'
   character(len=*), parameter :: p_blinking   =achar(27)//'[5m'
   character(len=*), parameter :: p_reverse    =achar(27)//'[7m'
   character(len=*), parameter :: p_inverted   =achar(27)//'[30m'//achar(27)//'[47m'
   character(len=*), parameter :: p_bg_black   =achar(27)//'[40m'
   character(len=*), parameter :: p_bg_red     =achar(27)//'[41m'
   character(len=*), parameter :: p_bg_green   =achar(27)//'[42m'
   character(len=*), parameter :: p_bg_brown   =achar(27)//'[43m'
   character(len=*), parameter :: p_bg_blue    =achar(27)//'[44m'
   character(len=*), parameter :: p_bg_purple  =achar(27)//'[45m'
   character(len=*), parameter :: p_bg_cyan    =achar(27)//'[46m'
   character(len=*), parameter :: p_bg_lgray   =achar(27)//'[47m' 


   integer, parameter :: p_tty=6
   !! Default TTY (terminal) output
!!ifdef splitlog
!   integer, parameter :: logUnit=69
!!else
   integer, parameter :: p_logUnit=6
!!endif
   !! Number of file to log
   integer :: iErrNumber
   !! integer to use with dumps

   real(kind=kind_rb), parameter :: p_adjust=0.01
   !! to convert from Pascal to mbar

   character(len=3), parameter, dimension(12) :: p_monthName=(/'jan','feb','mar' &
                                                            , 'apr','may','jun' &
                                                            , 'jul','aug','sep' &
                                                            , 'oct','nov','dec'/)

   character(len=3), parameter, dimension(7)  ::  p_weekName=(/'sun','mon','tue' &
                                                            , 'wed','thu','fri' &
                                                            , 'sat'/)

end module modConstants

