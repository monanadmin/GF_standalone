module modGate
  IMPLICIT NONE

  !-for BRAMS runs, set use_gate=.false.
  LOGICAL, PARAMETER :: p_use_gate = .true. 
  
  !- Here are the place for data related with the GATE soundings
  INTEGER, PARAMETER :: p_gate = 1   ! flag to turn on/off : 1/0
  INTEGER, PARAMETER :: p_klon = 161 ! number of soundings for gate
  INTEGER, PARAMETER :: p_klev = 41  ! number of vertical levels
  INTEGER, PARAMETER :: p_ktrac= 2   ! number of chemical tracers
  INTEGER, PARAMETER :: p_levs=p_klev  
  INTEGER, PARAMETER :: p_nvar_grads=300
  
  TYPE t_cupout_vars
    REAL             ,POINTER     :: varp(:,:)
    CHARACTER(LEN=80),ALLOCATABLE :: varn(:)
  END  TYPE t_cupout_vars
  
  TYPE (t_cupout_vars), ALLOCATABLE :: cupout(:)
  
  REAL, DIMENSION(p_klon,p_klev):: pgeo,ppres,ptemp,pq,pu,pv,pvervel, &
                               zrvten,ztten,zq1,zq2,zqr,zadvt,zadvq

  INTEGER :: JL,KLEV_SOUND
  character (len=128) :: runname, runlabel,rundata="NONE" 

end module modGate
