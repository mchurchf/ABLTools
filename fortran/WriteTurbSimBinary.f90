program wrtTurbSim


!SUBROUTINE WrBinTURBSIM( V, NumGrid_Y, NumGrid_Z, GridRes_Y, GridRes_Z, TimeStep, UHub, HubHt, Z1, DescStr, RootName )
   ! -------------------------------------------------------------------------------------------------------
   ! This is for dimension 2 of V(:,:,:)....
   !
   ! The grid of points on the Cartesian plane is numbered in the following way (notice that the first
   ! height starts at the bottom of the grid):
   !
   !               Y(1)                        Y(2)             Y(3)             ... Y(NumGrid_Y)
   !              -------------------------------------------------------------------------------------------
   ! Z(NumGrid_Z):|V(NumGrid_Y*(NumGrid_Z-1)+1) ...                                  V(NTot) |
   ! ...          |...                                                                                      |
   ! Z(2)        :|V(NumGrid_Y + 1)            V(NumGrid_Y + 2) V(NumGrid_Y + 3) ... V(2*NumGrid_Y)         |
   ! Z(1)        :|V(1)                        V(2)             V(3)             ... V(NumGrid_Y)           |
   !              -------------------------------------------------------------------------------------------
   !
   ! Z(i) < Z(i+1) for all integers i, 1 <= i < NumGrid_Z
   ! Y(j) < Y(j+1) for all integers j, 1 <= j < NumGrid_Y
   !
   ! If an extra hub point is necessary because the point does not fall on the grid,
   ! then it is added immediately following the regular grid points, i.e.
   ! Hub point = NumGrid_Y * NumGrid_Z + 1.
   !
   ! If the tower wind file output is selected, those extra points (in a single vertical
   ! line) are added at the end, after the hub point.
   ! --------------------------------------------------------------------------------------------------------

USE NWTC_Library

IMPLICIT                  NONE


! !..............................................................................................
! ! these are defined in the NWTC Subroutine library:
! !
! INTEGER, PARAMETER              :: B1Ki     = SELECTED_INT_KIND(  2 )           ! Kind for one-byte whole numbers
! INTEGER, PARAMETER              :: B2Ki     = SELECTED_INT_KIND(  4 )           ! Kind for two-byte whole numbers
! INTEGER, PARAMETER              :: B4Ki     = SELECTED_INT_KIND(  9 )           ! Kind for four-byte whole numbers
! INTEGER, PARAMETER              :: B8Ki     = SELECTED_INT_KIND( 18 )           ! Kind for eight-byte whole numbers
!
! INTEGER, PARAMETER              :: QuKi     = SELECTED_REAL_KIND( 20, 500 )     ! Kind for 16-byte, floating-point numbers
! INTEGER, PARAMETER              :: R8Ki     = SELECTED_REAL_KIND( 14, 300 )     ! Kind for eight-byte floating-point numbers
! INTEGER, PARAMETER              :: SiKi     = SELECTED_REAL_KIND(  6,  30 )     ! Kind for four-byte, floating-point numbers
!
!
!       ! The default kinds for reals and integers:
!
! INTEGER, PARAMETER              :: IntKi    = B4Ki                              ! Default kind for integers
! INTEGER, PARAMETER              :: ReKi     = SiKi                              ! Default kind for floating-point numbers
! INTEGER, PARAMETER              :: DbKi     = R8Ki                              ! Default kind for double floating-point numbers
! !..............................................................................................


INTEGER(4)                :: Sttus                      ! Status returned by an attempted allocation.
REAL(Reki), ALLOCATABLE   :: V(:,:,:)                   ! An array containing the velocity matrix: (NumSteps,NTot,3).
REAL(ReKi)                :: Z1                         ! The first element of the Z(:) array: vertical locations of the points., i.e., the lowest height on the grid
CHARACTER                 :: DescStr*7                  ! String used to describe the file

   ! local data

REAL(SiKi), PARAMETER     :: IntMax   =  32767.0
REAL(SiKi), PARAMETER     :: IntMin   = -32768.0
REAL(SiKi), PARAMETER     :: IntRng   = IntMax - IntMin ! Max Range of 2-byte integer

REAL(SiKi)                :: UOff                       ! Offset for the U component
REAL(SiKi)                :: UScl                       ! Slope  for the U component
REAL(ReKi)                :: VMax(3)                    ! Maximum value of the 3 wind components
REAL(ReKi)                :: VMin(3)                    ! Minimum value of the 3 wind components
REAL(SiKi)                :: VOff                       ! Offset for the V component
REAL(SiKi)                :: VScl                       ! Slope  for the V component
REAL(SiKi)                :: WOff                       ! Offset for the W component
REAL(SiKi)                :: WScl                       ! Slope  for the W component

INTEGER, PARAMETER        :: DecRound  = 3              ! Number of decimal places to round to
INTEGER                   :: IC                         ! counter for the velocity component of V
INTEGER                   :: II                         ! counter for the point on the grid/tower
INTEGER                   :: IT                         ! counter for the timestep
INTEGER(B4Ki)             :: LenDesc                    ! Length of the description string
INTEGER(B4Ki)             :: NumGrid                    ! Number of points on the grid
INTEGER(B4Ki)             :: NumTower                   ! Number of points on the tower
INTEGER(B4Ki)             :: TwrStart                   ! First index of a tower point
INTEGER(B4Ki)             :: TwrTop                     ! The index of top of the tower (it could be on the grid instead of at the end)

INTEGER(B4Ki)             :: IP
INTEGER(B2Ki),ALLOCATABLE :: TmpVarray(:)               ! This array holds the normalized velocities before being written to the binary file

INTEGER                   :: AllocStat

!...................................................
! variables from TSSubs, that I've set here:
!...................................................

LOGICAL, PARAMETER        :: WrADTWR    = .FALSE.       ! Flag to determine if tower points should be added
LOGICAL, PARAMETER        :: ExtraTwrPt = .FALSE.       ! Flag to indicate if the tower is on the regular grid or if an extra point must be added
LOGICAL, PARAMETER        :: ExtraHubPT = .FALSE.       ! Flag to indicate if the hub is on the regular grid or if an extra point must be added

INTEGER, PARAMETER        :: UAFFW    = 9               ! I/O unit for AeroDyn FF data (*.bts file).
INTEGER                   :: NTot                       ! Number of points in grid, plus the hub center.
INTEGER                   :: NumOutSteps                ! Number of output time steps.
REAL(ReKi)                :: Z(1)                       ! Vertical locations of the points


integer                   :: i
integer                   :: j 
integer                   :: k
character                 :: fileID*10

REAL(ReKi)                :: tmp_u                      
REAL(ReKi)                :: tmp_v                      
REAL(ReKi)                :: tmp_w

REAL(ReKi)                :: tmp_x                   
REAL(ReKi)                :: tmp_y                      
REAL(ReKi)                :: tmp_z


REAL(ReKi)                :: pi_deg

real(ReKi), ALLOCATABLE   :: xpos(:)
real(ReKi), ALLOCATABLE   :: ypos(:)
real(ReKi), ALLOCATABLE   :: zpos(:)


!!!  input data to read in
REAL(ReKi)                :: theta                      ! inflow angle
INTEGER                   :: NumGrid_Y                  ! Grid dimension. (in horizontal direction)
INTEGER                   :: NumGrid_Z                  ! Grid dimension. (in vertical direction)
REAL(ReKi)                :: GridRes_Y                  ! Distance between two consecutive horizontal points on the grid (Horizontal resolution)
REAL(ReKi)                :: GridRes_Z                  ! Distance between two consecutive vertical points on the grid (Vertical resolution)
REAL(ReKi)                :: TimeStep                   ! Time step.

integer                   :: nt                         ! number of time-series data planes
REAL(ReKi)                :: UHub                       ! Hub-height (total) wind speed (m/s)
REAL(ReKi)                :: HubHt                      ! Hub height.
CHARACTER                 :: timestr*10
REAL(ReKi)                :: curr_time
CHARACTER                 :: RootName*16                ! Root name of the I/O files
CHARACTER                 :: fileName*24                ! Root name of the I/O files


!!! read file
open(35, file='tsConv.inp')
  read(35,*) theta
  read(35,*) NumGrid_Y
  read(35,*) NumGrid_Z
  read(35,*) GridRes_Z
  read(35,*) GridRes_Y
  read(35,*) TimeStep
  read(35,*) nt
  read(35,*) UHub
  read(35,*) HubHt
  read(35,*) RootName
  read(35,*) fileName
  read(35,*) curr_time
close(35)


NTot = NumGrid_Y*NumGrid_Z
pi_deg = 2.0*asin(1.0)
theta = theta*pi_deg/180.0

NumOutSteps = nt
Z(1)        = 0.0   !Z1
NumTower = 0


!!! allocate dynamic memory !!!

ALLOCATE ( V(nt,NTot,3) , STAT=Sttus )
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error allocating memory for the V array.' )
ENDIF

ALLOCATE ( xpos(NTot) , STAT=Sttus )
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error allocating memory for the xpos array.' )
ENDIF

ALLOCATE ( ypos(NTot) , STAT=Sttus )
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error allocating memory for the ypos array.' )
ENDIF

ALLOCATE ( zpos(NTot) , STAT=Sttus )
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error allocating memory for the zpos array.' )
ENDIF


CALL OpenBin ( UAFFW, TRIM(RootName)//'.bts', 2 )


do i = 1, nt
  
  write(fileID,"(f9.3)") curr_time
  write(*,*) fileID
  
  open(25, file=TRIM(fileID)//'/'//TRIM(fileName))

  do k = 1, NumGrid_Z 
    do j = 1, NumGrid_Y 

      read(25,*) xpos(j+NumGrid_Y*(k-1)), ypos(j+NumGrid_Y*(k-1)), zpos(j+NumGrid_Y*(k-1)), tmp_u, tmp_v, tmp_w     
     
      V(i,j+NumGrid_Y*(k-1),1) = cos(theta)*tmp_u - sin(theta)*tmp_v
      V(i,j+NumGrid_Y*(k-1),2) = sin(theta)*tmp_u + cos(theta)*tmp_v
      V(i,j+NumGrid_Y*(k-1),3) = tmp_w
      
    enddo   
  enddo
  close(25)
  
  curr_time = curr_time + TimeStep
  
enddo



! subroutine from TurbSim:

      ! Find the range of our velocity

   DO IC=1,3

         ! Initialize the Min/Max values

      VMin(IC) = V(1,1,IC)
      VMax(IC) = V(1,1,IC)

      DO II=1,NTot   ! Let's check all of the points
         DO IT=1,NumOutSteps  ! Use only the number of timesteps requested originally

            IF ( V(IT,II,IC) > VMax(IC) ) THEN

               VMax(IC) = V(IT,II,IC)

            ELSEIF ( V(IT,II,IC) < VMin(IC) ) THEN

               VMin(IC) = V(IT,II,IC)

            ENDIF

         ENDDO !IT
      ENDDO !II

   ENDDO !IC


      ! Calculate the scaling parameters for each component


   IF ( VMax(1) == VMin(1) ) THEN
      UScl = 1
   ELSE
      UScl = IntRng/REAL( VMax(1) - VMin(1) , SiKi )
   ENDIF

   IF ( VMax(2) == VMin(2) ) THEN
      VScl = 1
   ELSE
      VScl = IntRng/REAL( VMax(2) - VMin(2) , SiKi )
   ENDIF

   IF ( VMax(3) == VMin(3) ) THEN
      WScl = 1
   ELSE
      WScl = IntRng/REAL( VMax(3) - VMin(3) , SiKi )
   ENDIF


   UOff = IntMin - UScl*REAL( VMin(1)    , SiKi )
   VOff = IntMin - VScl*REAL( VMin(2)    , SiKi )
   WOff = IntMin - WScl*REAL( VMin(3)    , SiKi )


      ! Find the first tower point

   NumGrid  = NTot

   IF ( WrADTWR ) THEN

      TwrStart = NumGrid + 1

      IF ( ExtraHubPT ) THEN
         TwrStart = TwrStart + 1
      ENDIF

      IF ( ExtraTwrPt ) THEN
         TwrTop   = TwrStart
         TwrStart = TwrStart + 1
      ELSE
         TwrTop = INT(NumGrid_Y / 2) + 1      ! The top tower point is on the grid where Z = 1
      ENDIF

      NumTower = Ntot - TwrStart + 2

   ELSE

      NumTower = 0

   ENDIF


   DescStr = 'testing1'


   LenDesc = LEN_TRIM( DescStr )             ! Length of the string that contains program name, version, date, and time


   !CALL WrScr ( ' Generating AeroDyn binary time-series file "'//TRIM( RootName )//'.bts"' )


      ! Write the header



   WRITE (UAFFW)   INT(   7                , B2Ki )          ! TurbSim format

   WRITE (UAFFW)   INT( NumGrid_Z          , B4Ki )          ! the number of grid points vertically
   WRITE (UAFFW)   INT( NumGrid_Y          , B4Ki )          ! the number of grid points laterally
   WRITE (UAFFW)   INT( NumTower           , B4Ki )          ! the number of tower points
   WRITE (UAFFW)   INT( NumOutSteps        , B4Ki )          ! the number of time steps

   WRITE (UAFFW)  REAL( GridRes_Z          , SiKi )          ! grid spacing in vertical direction, in m
   WRITE (UAFFW)  REAL( GridRes_Y          , SiKi )          ! grid spacing in lateral direction, in m
   WRITE (UAFFW)  REAL( TimeStep           , SiKi )          ! grid spacing in delta time, in m/s
   WRITE (UAFFW)  REAL( UHub               , SiKi )          ! the mean wind speed in m/s at hub height
   WRITE (UAFFW)  REAL( HubHt              , SiKi )          ! the hub height, in m
   WRITE (UAFFW)  REAL( Z(1)               , SiKi )          ! the height of the grid bottom, in m

   WRITE (UAFFW)  REAL( UScl               , SiKi )          ! the U-component slope for scaling
   WRITE (UAFFW)  REAL( UOff               , SiKi )          ! the U-component offset for scaling
   WRITE (UAFFW)  REAL( VScl               , SiKi )          ! the V-component slope for scaling
   WRITE (UAFFW)  REAL( VOff               , SiKi )          ! the V-component offset for scaling
   WRITE (UAFFW)  REAL( WScl               , SiKi )          ! the W-component slope for scaling
   WRITE (UAFFW)  REAL( WOff               , SiKi )          ! the W-component offset for scaling

   WRITE (UAFFW)   INT( LenDesc            , B4Ki )          ! the number of characters in the string, max 200

   DO II=1,LenDesc

      WRITE (UAFFW)  INT( IACHAR( DescStr(II:II) ), B1Ki )   ! converted ASCII characters

   ENDDO

      ALLOCATE ( TmpVarray( 3*(NumGrid_Z*NumGrid_Y + NumTower) ) , STAT=AllocStat )

      !IF ( AllocStat /= 0 )  THEN
      !   CALL TS_Abort ( 'Error allocating memory for temporary wind speed array.' )
      !ENDIF

      ! Loop through time.

   DO IT=1,NumOutSteps  !Use only the number of timesteps requested originally

         ! Write out grid data in binary form. II = (IZ - 1)*NumGrid_Y + IY, IY varies most rapidly

      IP = 1

      DO II=1,NumGrid

         TmpVarray(IP)   =  NINT( Max( Min( REAL(UScl*V(IT,II,1) + UOff, SiKi), IntMax ),IntMin) , B2Ki )
         TmpVarray(IP+1) =  NINT( Max( Min( REAL(VScl*V(IT,II,2) + VOff, SiKi), IntMax ),IntMin) , B2Ki )
         TmpVarray(IP+2) =  NINT( Max( Min( REAL(WScl*V(IT,II,3) + Woff, SiKi), IntMax ),IntMin) , B2Ki )

         IP = IP + 3
      ENDDO ! II


      IF ( WrADTWR ) THEN

            ! Write out the tower data in binary form

            ! Value at the top of the tower (bottom of grid)
         TmpVarray(IP)   =  NINT( Max( Min( REAL(UScl*V(IT,TwrTop,1) + UOff, SiKi), IntMax ),IntMin) , B2Ki )
         TmpVarray(IP+1) =  NINT( Max( Min( REAL(VScl*V(IT,TwrTop,2) + VOff, SiKi), IntMax ),IntMin) , B2Ki )
         TmpVarray(IP+2) =  NINT( Max( Min( REAL(WScl*V(IT,TwrTop,3) + Woff, SiKi), IntMax ),IntMin) , B2Ki )

         IP = IP + 3
         DO II=TwrStart,NTot
                ! Values of tower data
            TmpVarray(IP)   =  NINT( Max( Min( REAL(UScl*V(IT,II,1) + UOff, SiKi), IntMax ),IntMin) , B2Ki )
            TmpVarray(IP+1) =  NINT( Max( Min( REAL(VScl*V(IT,II,2) + VOff, SiKi), IntMax ),IntMin) , B2Ki )
            TmpVarray(IP+2) =  NINT( Max( Min( REAL(WScl*V(IT,II,3) + Woff, SiKi), IntMax ),IntMin) , B2Ki )

            IP = IP + 3
         ENDDO ! II

      ENDIF

      WRITE ( UAFFW ) TmpVarray(:)
   ENDDO ! IT

   CLOSE ( UAFFW )

   IF ( ALLOCATED( TmpVarray ) ) DEALLOCATE( TmpVarray )

   
   IF ( ALLOCATED( V ) ) DEALLOCATE( V ) 

END Program wrtTurbSim 
