MODULE SharedDataTypes
!
!----------------------------------------------------------------------------------------------------

   USE                              Precision

   IMPLICIT NONE

   !-------------------------------------------------------------------------------------------------
   ! Public types and subroutines
   !-------------------------------------------------------------------------------------------------
   TYPE, PUBLIC :: ProgDesc !bjj: do we want to add date as well?
      CHARACTER(20)              :: Name
      CHARACTER(99)              :: Ver
   END TYPE ProgDesc


   TYPE, PUBLIC :: Marker
      REAL(ReKi)                 :: Position(3)
      REAL(ReKi)                 :: Orientation(3,3)     ! Direction Cosine Matrix (DCM)
      REAL(ReKi)                 :: TranslationVel(3)
      REAL(ReKi)                 :: RotationVel(3)
   END TYPE Marker

   TYPE, PUBLIC :: Load
      REAL(ReKi)                 :: Force(3)
      REAL(ReKi)                 :: Moment(3)
      REAL(ReKi)                 :: AddedMass(6,6)
   END TYPE Load


      ! Data types developed for AeroDyn

   TYPE, PUBLIC :: AllAeroMarkers
      TYPE(Marker), ALLOCATABLE  :: Blade(:,:)
      TYPE(Marker), ALLOCATABLE  :: Hub(:)
      TYPE(Marker), ALLOCATABLE  :: RotorFurl(:)
      TYPE(Marker), ALLOCATABLE  :: Nacelle(:)
      TYPE(Marker), ALLOCATABLE  :: Tower(:)
      TYPE(Marker), ALLOCATABLE  :: Tail(:)
   END TYPE AllAeroMarkers

   TYPE, PUBLIC :: AllAeroLoads
      TYPE(Load),  ALLOCATABLE   :: Blade(:,:)
      TYPE(Load),  ALLOCATABLE   :: Hub(:)
      TYPE(Load),  ALLOCATABLE   :: RotorFurl(:)
      TYPE(Load),  ALLOCATABLE   :: Nacelle(:)
      TYPE(Load),  ALLOCATABLE   :: Tower(:)
      TYPE(Load),  ALLOCATABLE   :: Tail(:)
   END TYPE AllAeroLoads


   TYPE, PUBLIC :: AeroConfig
      TYPE(Marker), ALLOCATABLE  :: Blade(:)
      REAL(ReKi)                 :: BladeLength
      TYPE(Marker)               :: Hub
      TYPE(Marker)               :: RotorFurl
      TYPE(Marker)               :: Nacelle
      TYPE(Marker)               :: TailFin
      TYPE(Marker)               :: Tower
!      TYPE(Marker)               :: Substructure
!      TYPE(Marker)               :: Foundation
   END TYPE AeroConfig


!      ! Data types developed for HydroDyn
!
!   TYPE, PUBLIC :: AllHydroMarkers
!      TYPE(Marker), ALLOCATABLE  :: Substructure (:)  
!   END TYPE
!   
!   TYPE, PUBLIC :: AllHydroLoads
!      TYPE(Load),   ALLOCATABLE  :: Substructure (:)  
!   END TYPE
!   
!   
!   TYPE, PUBLIC :: HydroConfig
!      TYPE(Marker)               :: Substructure
!   END TYPE HydroConfig
   
   
   
   
!====================================================================================================
!CONTAINS
!====================================================================================================
!====================================================================================================

END MODULE SharedDataTypes
