! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to declare a structure for output fields.
!
! Description:
!   This module contains the declaration of the structure
!   used to store output fields for the radiation code.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiance Core
!
!------------------------------------------------------------------------------
MODULE def_out

USE realtype_rd

IMPLICIT NONE


TYPE StrOut

! Fluxes and radiances calculated
  REAL (RealK), ALLOCATABLE :: flux_direct(:, :, :)
!   Direct flux
  REAL (RealK), ALLOCATABLE :: flux_down(:, :, :)
!   Total downward flux
  REAL (RealK), ALLOCATABLE :: flux_up(:, :, :)
!   Upward flux
  REAL (RealK), ALLOCATABLE :: flux_direct_clear(:, :, :)
!   Clear-sky direct flux
  REAL (RealK), ALLOCATABLE :: flux_down_clear(:, :, :)
!   Clear-sky downward flux
  REAL (RealK), ALLOCATABLE :: flux_up_clear(:, :, :)
!   Clear-sky upward flux
  REAL (RealK), ALLOCATABLE :: radiance(:, :, :, :)
!   Radiances
  REAL (RealK), ALLOCATABLE :: photolysis(:, :, :)
!   Rates of photolysis
  REAL (RealK), ALLOCATABLE :: solar_tail_flux(:)
!   Solar tail flux considered in LW region

! Diagnostic fluxes
  REAL (RealK), ALLOCATABLE :: flux_up_tile(:, :, :)
!   Upward fluxes at tiled surface points
  REAL (RealK), ALLOCATABLE :: flux_up_blue_tile(:, :, :)
!   Upward blue fluxes at tiled surface points
  REAL (RealK), ALLOCATABLE :: flux_direct_blue_surf(:)
!   Direct blue flux at the surface
  REAL (RealK), ALLOCATABLE :: flux_down_blue_surf(:)
!   Total downward blue flux at the surface
  REAL (RealK), ALLOCATABLE :: flux_up_blue_surf(:)
!   Upward blue flux at the surface
  REAL (RealK), ALLOCATABLE :: flux_direct_diag(:, :, :)
!   Diagnostic direct flux weighted by band
  REAL (RealK), ALLOCATABLE :: flux_down_diag(:, :, :)
!   Diagnostic total downward flux weighted by band
  REAL (RealK), ALLOCATABLE :: flux_up_diag(:, :, :)
!   Diagnostic upward flux weighted by band
  REAL (RealK), ALLOCATABLE :: flux_down_diag_surf(:)
!   Diagnostic total downward flux at the surface weighted by band
  REAL (RealK), ALLOCATABLE :: flux_down_clear_diag_surf(:)
!   Diagnostic clear-sky downward flux at the surface weighted by band
  REAL (RealK), ALLOCATABLE :: flux_direct_band(:, :, :)
!   Direct flux per band
  REAL (RealK), ALLOCATABLE :: flux_down_band(:, :, :)
!   Total downward flux per band
  REAL (RealK), ALLOCATABLE :: flux_up_band(:, :, :)
!   Upward flux per band
  REAL (RealK), ALLOCATABLE :: flux_direct_clear_band(:, :, :)
!   Clear-sky direct flux per band
  REAL (RealK), ALLOCATABLE :: flux_down_clear_band(:, :, :)
!   Clear-sky downward flux per band
  REAL (RealK), ALLOCATABLE :: flux_up_clear_band(:, :, :)
!   Clear-sky upward flux per band

! Cloud diagnostics
  REAL (RealK), ALLOCATABLE :: tot_cloud_cover(:)
!   Total cloud cover
  REAL (RealK), ALLOCATABLE :: cloud_absorptivity(:, :)
!   Absorptivity of cloud weighted by cloud fraction
!   and upward clear-sky infra-red flux.
  REAL (RealK), ALLOCATABLE :: cloud_weight_absorptivity(:, :)
!   Weights to be applied to absorptivies.
  REAL (RealK), ALLOCATABLE :: ls_cloud_absorptivity(:, :)
!   Absorptivity of layer cloud weighted by cloud fraction
!   and upward clear-sky infra-red flux.
  REAL (RealK), ALLOCATABLE :: ls_cloud_weight_absorptivity(:, :)
!   Weights to be applied to layer cld. absorptivies.
  REAL (RealK), ALLOCATABLE :: cnv_cloud_absorptivity(:, :)
!   Absorptivity of conv.cloud weighted by cloud fraction
!   and upward clear-sky infra-red flux.
  REAL (RealK), ALLOCATABLE :: cnv_cloud_weight_absorptivity(:, :)
!   Weights to be applied to conv.cld absorptivies.
  REAL (RealK), ALLOCATABLE :: cloud_extinction(:, :)
!   Absorptivity of cloud weighted by cloud fraction
!   and downward clear-sky solar flux.
  REAL (RealK), ALLOCATABLE :: cloud_weight_extinction(:, :)
!   Weights to be applied to extinctions.
  REAL (RealK), ALLOCATABLE :: ls_cloud_extinction(:, :)
!   Absorptivity of layer cloud weighted by cloud fraction
!   and downward clear-sky solar flux.
  REAL (RealK), ALLOCATABLE :: ls_cloud_weight_extinction(:, :)
!   Weights to be applied to layer cld. extinctions.
  REAL (RealK), ALLOCATABLE :: cnv_cloud_extinction(:, :)
!   Absorptivity of conv.cloud weighted by cloud fraction
!   and downward clear-sky solar flux.
  REAL (RealK), ALLOCATABLE :: cnv_cloud_weight_extinction(:, :)
!   Weights to be applied to conv. cld. extinctions.

! Aerosol diagnostics
  REAL (RealK), ALLOCATABLE :: aerosol_absorption_band(:, :, :)
!   Total aerosol absorption per band
  REAL (RealK), ALLOCATABLE :: aerosol_scattering_band(:, :, :)
!   Total aerosol scattering per band
  REAL (RealK), ALLOCATABLE :: aerosol_asymmetry_band(:, :, :)
!   Total aerosol asymmetry (weighted by scattering) per band

END TYPE StrOut


CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE allocate_out(radout, control, dimen, sp)

USE def_control,  ONLY: StrCtrl
USE def_dimen,    ONLY: StrDim
USE def_spectrum, ONLY: StrSpecData

IMPLICIT NONE

TYPE (StrOut),      INTENT(INOUT) :: radout
TYPE (StrCtrl),     INTENT(IN)    :: control
TYPE (StrDim),      INTENT(IN)    :: dimen
TYPE (StrSpecData), INTENT(IN)    :: sp

IF (.NOT. ALLOCATED(radout%flux_direct))                                       &
  ALLOCATE(radout%flux_direct                  ( dimen%nd_flux_profile,        &
                                                 0: dimen%nd_layer,            &
                                                 dimen%nd_channel            ))

IF (.NOT. ALLOCATED(radout%flux_down))                                         &
  ALLOCATE(radout%flux_down                    ( dimen%nd_flux_profile,        &
                                                 0: dimen%nd_layer,            &
                                                 dimen%nd_channel            ))

IF (.NOT. ALLOCATED(radout%flux_up))                                           &
  ALLOCATE(radout%flux_up                      ( dimen%nd_flux_profile,        &
                                                 0: dimen%nd_layer,            &
                                                 dimen%nd_channel            ))

IF (.NOT. ALLOCATED(radout%flux_direct_clear))                                 &
  ALLOCATE(radout%flux_direct_clear            ( dimen%nd_2sg_profile,         &
                                                 0: dimen%nd_layer,            &
                                                 dimen%nd_channel            ))

IF (.NOT. ALLOCATED(radout%flux_down_clear))                                   &
  ALLOCATE(radout%flux_down_clear              ( dimen%nd_2sg_profile,         &
                                                 0: dimen%nd_layer,            &
                                                 dimen%nd_channel            ))

IF (.NOT. ALLOCATED(radout%flux_up_clear))                                     &
  ALLOCATE(radout%flux_up_clear                ( dimen%nd_2sg_profile,         &
                                                 0: dimen%nd_layer,            &
                                                 dimen%nd_channel            ))

IF (.NOT. ALLOCATED(radout%radiance))                                          &
  ALLOCATE(radout%radiance                     ( dimen%nd_radiance_profile,    &
                                                 dimen%nd_viewing_level,       &
                                                 dimen%nd_direction,           &
                                                 dimen%nd_channel            ))

IF (.NOT. ALLOCATED(radout%photolysis))                                        &
  ALLOCATE(radout%photolysis                   ( dimen%nd_j_profile,           &
                                                 dimen%nd_viewing_level,       &
                                                 dimen%nd_channel            ))

IF (.NOT. ALLOCATED(radout%solar_tail_flux))                                   &
  ALLOCATE(radout%solar_tail_flux              ( dimen%nd_profile            ))


IF (.NOT. ALLOCATED(radout%flux_up_tile))                                      &
  ALLOCATE(radout%flux_up_tile                 ( dimen%nd_point_tile,          &
                                                 dimen%nd_tile,                &
                                                 dimen%nd_channel            ))

IF (.NOT. ALLOCATED(radout%flux_up_blue_tile))                                 &
  ALLOCATE(radout%flux_up_blue_tile            ( dimen%nd_point_tile,          &
                                                 dimen%nd_tile,                &
                                                 dimen%nd_channel            ))

IF (.NOT. ALLOCATED(radout%flux_direct_blue_surf))                             &
  ALLOCATE(radout%flux_direct_blue_surf        ( dimen%nd_flux_profile       ))

IF (.NOT. ALLOCATED(radout%flux_down_blue_surf))                               &
  ALLOCATE(radout%flux_down_blue_surf          ( dimen%nd_flux_profile       ))

IF (.NOT. ALLOCATED(radout%flux_up_blue_surf))                                 &
  ALLOCATE(radout%flux_up_blue_surf            ( dimen%nd_flux_profile       ))

IF (control%l_flux_direct_diag) THEN
  IF (.NOT. ALLOCATED(radout%flux_direct_diag))                                &
    ALLOCATE(radout%flux_direct_diag           ( dimen%nd_flux_profile,        &
                                                 0: dimen%nd_layer,            &
                                                 dimen%nd_channel            ))
END IF

IF (control%l_flux_down_diag) THEN
  IF (.NOT. ALLOCATED(radout%flux_down_diag))                                  &
    ALLOCATE(radout%flux_down_diag             ( dimen%nd_flux_profile,        &
                                                 0: dimen%nd_layer,            &
                                                 dimen%nd_channel            ))
END IF

IF (control%l_flux_up_diag) THEN
  IF (.NOT. ALLOCATED(radout%flux_up_diag))                                    &
    ALLOCATE(radout%flux_up_diag               ( dimen%nd_flux_profile,        &
                                                 0: dimen%nd_layer,            &
                                                 dimen%nd_channel            ))
END IF

IF (control%l_flux_down_diag_surf) THEN
  IF (.NOT. ALLOCATED(radout%flux_down_diag_surf))                             &
    ALLOCATE(radout%flux_down_diag_surf        ( dimen%nd_flux_profile       ))
END IF

IF (control%l_flux_down_clear_diag_surf) THEN
  IF (.NOT. ALLOCATED(radout%flux_down_clear_diag_surf))                       &
    ALLOCATE(radout%flux_down_clear_diag_surf  ( dimen%nd_flux_profile       ))
END IF

IF (control%l_flux_direct_band) THEN
  IF (.NOT. ALLOCATED(radout%flux_direct_band))                                &
    ALLOCATE(radout%flux_direct_band           ( dimen%nd_flux_profile,        &
                                                 0: dimen%nd_layer,            &
                                                 sp%dim%nd_band              ))
END IF

IF (control%l_flux_down_band) THEN
  IF (.NOT. ALLOCATED(radout%flux_down_band))                                  &
    ALLOCATE(radout%flux_down_band             ( dimen%nd_flux_profile,        &
                                                 0: dimen%nd_layer,            &
                                                 sp%dim%nd_band              ))
END IF

IF (control%l_flux_up_band) THEN
  IF (.NOT. ALLOCATED(radout%flux_up_band))                                    &
    ALLOCATE(radout%flux_up_band               ( dimen%nd_flux_profile,        &
                                                 0: dimen%nd_layer,            &
                                                 sp%dim%nd_band              ))
END IF

IF (control%l_flux_direct_clear_band) THEN
  IF (.NOT. ALLOCATED(radout%flux_direct_clear_band))                          &
    ALLOCATE(radout%flux_direct_clear_band     ( dimen%nd_2sg_profile,         &
                                                 0: dimen%nd_layer,            &
                                                 sp%dim%nd_band              ))
END IF

IF (control%l_flux_down_clear_band) THEN
  IF (.NOT. ALLOCATED(radout%flux_down_clear_band))                            &
    ALLOCATE(radout%flux_down_clear_band       ( dimen%nd_2sg_profile,         &
                                                 0: dimen%nd_layer,            &
                                                 sp%dim%nd_band              ))
END IF

IF (control%l_flux_up_clear_band) THEN
  IF (.NOT. ALLOCATED(radout%flux_up_clear_band))                              &
    ALLOCATE(radout%flux_up_clear_band         ( dimen%nd_2sg_profile,         &
                                                 0: dimen%nd_layer,            &
                                                 sp%dim%nd_band              ))
END IF

IF (.NOT. ALLOCATED(radout%tot_cloud_cover))                                   &
  ALLOCATE(radout%tot_cloud_cover              ( dimen%nd_profile            ))

IF (control%l_cloud_absorptivity) THEN
  IF (.NOT. ALLOCATED(radout%cloud_absorptivity))                              &
    ALLOCATE(radout%cloud_absorptivity         ( dimen%nd_profile,             &
                                                 dimen%nd_layer              ))
  IF (.NOT. ALLOCATED(radout%cloud_weight_absorptivity))                       &
    ALLOCATE(radout%cloud_weight_absorptivity  ( dimen%nd_profile,             &
                                                 dimen%nd_layer              ))
END IF

IF (control%l_ls_cloud_absorptivity) THEN
  IF (.NOT. ALLOCATED(radout%ls_cloud_absorptivity))                           &
    ALLOCATE(radout%ls_cloud_absorptivity      ( dimen%nd_profile,             &
                                                 dimen%nd_layer              ))
  IF (.NOT. ALLOCATED(radout%ls_cloud_weight_absorptivity))                    &
    ALLOCATE(radout%ls_cloud_weight_absorptivity ( dimen%nd_profile,           &
                                                   dimen%nd_layer            ))
END IF

IF (control%l_cnv_cloud_absorptivity) THEN
  IF (.NOT. ALLOCATED(radout%cnv_cloud_absorptivity))                          &
    ALLOCATE(radout%cnv_cloud_absorptivity     ( dimen%nd_profile,             &
                                                 dimen%nd_layer              ))
  IF (.NOT. ALLOCATED(radout%cnv_cloud_weight_absorptivity))                   &
    ALLOCATE(radout%cnv_cloud_weight_absorptivity ( dimen%nd_profile,          &
                                                    dimen%nd_layer           ))
END IF

IF (control%l_cloud_extinction) THEN
  IF (.NOT. ALLOCATED(radout%cloud_extinction))                                &
    ALLOCATE(radout%cloud_extinction           ( dimen%nd_profile,             &
                                                 dimen%nd_layer              ))
  IF (.NOT. ALLOCATED(radout%cloud_weight_extinction))                         &
    ALLOCATE(radout%cloud_weight_extinction    ( dimen%nd_profile,             &
                                                 dimen%nd_layer              ))
END IF

IF (control%l_ls_cloud_extinction) THEN
  IF (.NOT. ALLOCATED(radout%ls_cloud_extinction))                             &
    ALLOCATE(radout%ls_cloud_extinction        ( dimen%nd_profile,             &
                                                 dimen%nd_layer              ))
  IF (.NOT. ALLOCATED(radout%ls_cloud_weight_extinction))                      &
    ALLOCATE(radout%ls_cloud_weight_extinction ( dimen%nd_profile,             &
                                                 dimen%nd_layer              ))
END IF

IF (control%l_cnv_cloud_extinction) THEN
  IF (.NOT. ALLOCATED(radout%cnv_cloud_extinction))                            &
    ALLOCATE(radout%cnv_cloud_extinction       ( dimen%nd_profile,             &
                                                 dimen%nd_layer              ))
  IF (.NOT. ALLOCATED(radout%cnv_cloud_weight_extinction))                     &
    ALLOCATE(radout%cnv_cloud_weight_extinction( dimen%nd_profile,             &
                                                 dimen%nd_layer              ))
END IF

IF (control%l_aerosol_absorption_band) THEN
  IF (.NOT. ALLOCATED(radout%aerosol_absorption_band))                         &
    ALLOCATE(radout%aerosol_absorption_band    ( dimen%nd_profile,             &
                                                 dimen%nd_layer,               &
                                                 sp%dim%nd_band              ))
END IF

IF (control%l_aerosol_scattering_band) THEN
  IF (.NOT. ALLOCATED(radout%aerosol_scattering_band))                         &
    ALLOCATE(radout%aerosol_scattering_band    ( dimen%nd_profile,             &
                                                 dimen%nd_layer,               &
                                                 sp%dim%nd_band              ))
END IF

IF (control%l_aerosol_asymmetry_band) THEN
  IF (.NOT. ALLOCATED(radout%aerosol_asymmetry_band))                          &
    ALLOCATE(radout%aerosol_asymmetry_band     ( dimen%nd_profile,             &
                                                 dimen%nd_layer,               &
                                                 sp%dim%nd_band              ))
END IF

END SUBROUTINE allocate_out
!------------------------------------------------------------------------------
SUBROUTINE deallocate_out(radout)

IMPLICIT NONE

TYPE (StrOut), INTENT(INOUT) :: radout

IF (ALLOCATED(radout%aerosol_asymmetry_band)) &
    DEALLOCATE(radout%aerosol_asymmetry_band)
IF (ALLOCATED(radout%aerosol_scattering_band)) &
    DEALLOCATE(radout%aerosol_scattering_band)
IF (ALLOCATED(radout%aerosol_absorption_band)) &
    DEALLOCATE(radout%aerosol_absorption_band)
IF (ALLOCATED(radout%cnv_cloud_weight_extinction)) &
    DEALLOCATE(radout%cnv_cloud_weight_extinction)
IF (ALLOCATED(radout%cnv_cloud_extinction)) &
    DEALLOCATE(radout%cnv_cloud_extinction)
IF (ALLOCATED(radout%ls_cloud_weight_extinction)) &
    DEALLOCATE(radout%ls_cloud_weight_extinction)
IF (ALLOCATED(radout%ls_cloud_extinction)) &
    DEALLOCATE(radout%ls_cloud_extinction)
IF (ALLOCATED(radout%cloud_weight_extinction)) &
    DEALLOCATE(radout%cloud_weight_extinction)
IF (ALLOCATED(radout%cloud_extinction)) &
    DEALLOCATE(radout%cloud_extinction)
IF (ALLOCATED(radout%cnv_cloud_weight_absorptivity)) &
    DEALLOCATE(radout%cnv_cloud_weight_absorptivity)
IF (ALLOCATED(radout%cnv_cloud_absorptivity)) &
    DEALLOCATE(radout%cnv_cloud_absorptivity)
IF (ALLOCATED(radout%ls_cloud_weight_absorptivity)) &
    DEALLOCATE(radout%ls_cloud_weight_absorptivity)
IF (ALLOCATED(radout%ls_cloud_absorptivity)) &
    DEALLOCATE(radout%ls_cloud_absorptivity)
IF (ALLOCATED(radout%cloud_weight_absorptivity)) &
    DEALLOCATE(radout%cloud_weight_absorptivity)
IF (ALLOCATED(radout%cloud_absorptivity)) &
    DEALLOCATE(radout%cloud_absorptivity)
IF (ALLOCATED(radout%tot_cloud_cover)) &
    DEALLOCATE(radout%tot_cloud_cover)
IF (ALLOCATED(radout%flux_up_clear_band)) &
    DEALLOCATE(radout%flux_up_clear_band)
IF (ALLOCATED(radout%flux_down_clear_band)) &
    DEALLOCATE(radout%flux_down_clear_band)
IF (ALLOCATED(radout%flux_direct_clear_band)) &
    DEALLOCATE(radout%flux_direct_clear_band)
IF (ALLOCATED(radout%flux_up_band)) &
    DEALLOCATE(radout%flux_up_band)
IF (ALLOCATED(radout%flux_down_band)) &
    DEALLOCATE(radout%flux_down_band)
IF (ALLOCATED(radout%flux_direct_band)) &
    DEALLOCATE(radout%flux_direct_band)
IF (ALLOCATED(radout%flux_down_clear_diag_surf)) &
    DEALLOCATE(radout%flux_down_clear_diag_surf)
IF (ALLOCATED(radout%flux_down_diag_surf)) &
    DEALLOCATE(radout%flux_down_diag_surf)
IF (ALLOCATED(radout%flux_up_diag)) &
    DEALLOCATE(radout%flux_up_diag)
IF (ALLOCATED(radout%flux_down_diag)) &
    DEALLOCATE(radout%flux_down_diag)
IF (ALLOCATED(radout%flux_direct_diag)) &
    DEALLOCATE(radout%flux_direct_diag)
IF (ALLOCATED(radout%flux_up_blue_surf)) &
    DEALLOCATE(radout%flux_up_blue_surf)
IF (ALLOCATED(radout%flux_down_blue_surf)) &
    DEALLOCATE(radout%flux_down_blue_surf)
IF (ALLOCATED(radout%flux_direct_blue_surf)) &
    DEALLOCATE(radout%flux_direct_blue_surf)
IF (ALLOCATED(radout%flux_up_blue_tile)) &
    DEALLOCATE(radout%flux_up_blue_tile)
IF (ALLOCATED(radout%flux_up_tile)) &
    DEALLOCATE(radout%flux_up_tile)
IF (ALLOCATED(radout%solar_tail_flux)) &
    DEALLOCATE(radout%solar_tail_flux)
IF (ALLOCATED(radout%photolysis)) &
    DEALLOCATE(radout%photolysis)
IF (ALLOCATED(radout%radiance)) &
    DEALLOCATE(radout%radiance)
IF (ALLOCATED(radout%flux_up_clear)) &
    DEALLOCATE(radout%flux_up_clear)
IF (ALLOCATED(radout%flux_down_clear)) &
    DEALLOCATE(radout%flux_down_clear)
IF (ALLOCATED(radout%flux_direct_clear)) &
    DEALLOCATE(radout%flux_direct_clear)
IF (ALLOCATED(radout%flux_up)) &
    DEALLOCATE(radout%flux_up)
IF (ALLOCATED(radout%flux_down)) &
    DEALLOCATE(radout%flux_down)
IF (ALLOCATED(radout%flux_direct)) &
    DEALLOCATE(radout%flux_direct)

END SUBROUTINE deallocate_out
!------------------------------------------------------------------------------

END MODULE def_out
