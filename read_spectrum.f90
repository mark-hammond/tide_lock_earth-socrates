! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutines to read a Spectral File.
!
! Description:
!   The file is opened with checks. Logical arrays are set.
!   The file is now read until a line beginning a block is found.
!   A routine is called to read the block which in turn calls
!   an appropriate subroutine depending on the type, sub-type
!   and version of the block.
!
!------------------------------------------------------------------------------
SUBROUTINE read_spectrum(file_spectral, Sp)

USE realtype_rd, ONLY: RealK
USE def_spectrum, ONLY: StrSpecData, allocate_spectrum
USE rad_pcf
USE dimensions_spec_ucf
USE gas_list_pcf
USE file_manager, ONLY: assign_file_unit, release_file_unit
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


! Dummy variables.
CHARACTER (LEN=*), INTENT(IN) :: file_spectral
!   Name of spectral file
TYPE (StrSpecData) :: Sp
!   Spectral data


! Local variables.
CHARACTER (LEN=80) :: line
!   Line read from file
CHARACTER (LEN=80) :: char_dum
!   Dummy charcater variable
CHARACTER (LEN=256) :: spectral_k
!   Name of extended spectral file
CHARACTER (LEN=256) :: spectral_var
!   Name of spectral variability file
CHARACTER (LEN=errormessagelength) :: iomessage
!   I/O error message
INTEGER :: ios
!   I/O error status
INTEGER :: ierr = i_normal
!   Error flag
INTEGER :: iu_spc, iu_spc1, iu_spc2
!   Unit numbers for I/O of the spectral file
INTEGER :: i_type
!   Type of block read in
INTEGER :: i_subtype
!   Subtype of block
INTEGER :: i_version
!   Version for type and subtype
INTEGER :: i
!   Loop variable
LOGICAL :: l_exist_k
!   Existence flag for file

! Local values for dimensions: used to shorten declarations later
INTEGER :: nd_band
!   Size allocated for spectral bands
INTEGER :: nd_exclude
!   Size allocated for excluded bands
INTEGER :: nd_k_term
!   Size allocated for k-terms
INTEGER :: nd_species
!   Size allocated for gaseous species
INTEGER :: nd_scale_variable
!   Size allocated for scaling variables
INTEGER :: nd_continuum
!   Size allocated for continua
INTEGER :: nd_drop_type
!   Size allocated for drop types
INTEGER :: nd_ice_type
!   Size allocated for ice crystal types
INTEGER :: nd_aerosol_species
!   Size allocated for aerosol species
INTEGER :: nd_thermal_coeff
!   Size allocated for thermal coefficients
INTEGER :: nd_cloud_parameter
!   Size allocated for cloud parameters
INTEGER :: nd_humidity
!   Size allocated for humidities
INTEGER :: nd_aod_wavel
!   Number of wavelengths for aerosol optical depths
INTEGER :: nd_phase_term
!   Size allocated for terms in the phase function
INTEGER :: nd_tmp
!   Number of reference temperature for k-terms
INTEGER :: nd_pre
!   Number of reference pressures for k-terms
INTEGER :: nd_mix
!   Number of eta for mixture absorbing species
INTEGER :: nd_band_mix_gas
!   Number of bands where mixed species exist
INTEGER :: nd_sub_band
!   Size allocated for spectral sub-bands (for spectral variability)
INTEGER :: nd_times
!   Size allocated for times (for spectral variability)

CHARACTER (LEN=errormessagelength) :: cmessage
CHARACTER (LEN=*), PARAMETER       :: RoutineName = 'READ_SPECTRUM'


! Initialise array sizes. Some cannot at present conveniently be inferred
! from the spectral file: eventually these should all be determined
! from a block of dimensions or by other dynamic means.
nd_band = 1
nd_exclude = npd_exclude
nd_k_term = npd_k_term
nd_species = 1
nd_scale_variable = npd_scale_variable
nd_continuum = npd_continuum
nd_drop_type = npd_drop_type
nd_ice_type  = npd_ice_type
nd_aerosol_species = 1
nd_thermal_coeff = 1
nd_cloud_parameter  = npd_cloud_parameter
nd_humidity  = npd_humidities
nd_aod_wavel = 0
nd_phase_term = npd_phase_term
nd_tmp = 0
nd_pre = 0
nd_mix = 0
nd_band_mix_gas = 0
nd_sub_band = 0
nd_times = 0

Sp%Dim%nd_type = npd_type

! Initialise other integer values to zero
Sp%Basic%n_band = 0
Sp%Gas%n_absorb = 0
Sp%Planck%n_deg_fit = 0
Sp%Aerosol%n_aerosol = 0
Sp%Aerosol%n_aerosol_mr = 0
Sp%Aerosol%n_aod_wavel = 0

! It is important to know which gas is water vapour for some
! applications: here we initialize INDEX_WATER to 0 to
! produce an error in the code if it is not reset to a legal value.
! This should guard against omitting water vapour when it is needed.
Sp%Cont%index_water=0

! Initialization of logical variables.
ALLOCATE(Sp%Basic%l_present(0:Sp%Dim%nd_type))
Sp%Basic%l_present = .FALSE.
Sp%Planck%l_planck_tbl = .FALSE.

! Get a unit to read the file.
CALL assign_file_unit(file_spectral, iu_spc, handler="fortran")

! Open the file for reading
OPEN(UNIT=iu_spc, FILE=file_spectral, IOSTAT=ios, STATUS='OLD', &
     ACTION='READ', IOMSG=iomessage)
IF (ios /= 0) THEN
  cmessage = 'Spectral file could not be opened: ' // TRIM(iomessage)
  ierr=i_err_fatal
  CALL ereport(RoutineName, ierr, cmessage)
END IF

! Check existence of extended file of k-term data and open.
i=INDEX(file_spectral, ' ') - 1
spectral_k = file_spectral(1:i) // '_k'
INQUIRE(FILE=spectral_k, EXIST=l_exist_k)

IF (l_exist_k) THEN
  CALL assign_file_unit(spectral_k, iu_spc1, handler="fortran")
  OPEN(UNIT=iu_spc1,FILE=spectral_k,IOSTAT=ios,STATUS='OLD',               &
       ACTION='READ',IOMSG=iomessage)
  IF (ios /= 0) THEN
    cmessage = 'Extended spectral file exists but could not be opened: '// &
                TRIM(iomessage)
    ierr=i_warning
    CALL ereport(RoutineName, ierr, cmessage)
    l_exist_k=.FALSE.
  END IF
  ! Skip over head of the file
  DO
    READ(iu_spc1, '(a1)') char_dum
    IF ( char_dum(1:1) == '*' ) EXIT
  END DO
END IF


! Read through the file processing the blocks of data as they are encountered.
! Each line is read into an internal file and then processed.
DO
  READ(iu_spc, '(a80)', IOSTAT=ios) line

  IF (ios /= 0) EXIT

  ! Locate a block header.
  IF (line(1:6) == '*BLOCK') THEN
    READ(line, FMT='(a15, i4, 12x, i4, 12x, i4)', IOSTAT=ios, IOMSG=iomessage) &
      char_dum, i_type, i_subtype, i_version
    IF (ios /= 0) THEN
      cmessage = 'Block header in spectrum is incorrect: ' // TRIM(iomessage)
      ierr=i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    END IF

    ! Read in the rest of the block.
    CALL read_block_int
    IF (ierr /= i_normal) THEN
      CALL ereport(RoutineName, ierr, cmessage)
    END IF

    ! Read in the termination statement
    READ(iu_spc, '(a4)') char_dum
    IF (char_dum(1:4) /= '*END') THEN
      WRITE(cmessage, '(a48, i3, a12, i3, a11, i3)') &
        'Spectral block is incorrectly terminated. Type =', i_type, &
        ': Sub-type =', i_subtype, &
        ': Version =', i_version
      ierr=i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    END IF

    ! The block has been properly read: record the data as present.
    Sp%Basic%l_present(i_type)= .TRUE.

  END IF

END DO

CLOSE(iu_spc)
CALL release_file_unit(iu_spc, handler="fortran")

IF (l_exist_k) THEN
  CLOSE(iu_spc1)
  CALL release_file_unit(iu_spc1, handler="fortran")
END IF

! Initialise number of aerosol mixing ratios to aerosols present
Sp%Dim%nd_aerosol_mr = nd_aerosol_species
Sp%Aerosol%n_aerosol_mr = Sp%Aerosol%n_aerosol

! Set the index of water.
DO i=1, Sp%Gas%n_absorb
  IF (Sp%Gas%type_absorb(i) == ip_h2o) Sp%Cont%index_water=i
END DO

! Set dimensions in spectrum type from local values.
Sp%Dim%nd_band = nd_band
Sp%Dim%nd_exclude = nd_exclude
Sp%Dim%nd_k_term = nd_k_term
Sp%Dim%nd_species = nd_species
Sp%Dim%nd_scale_variable = nd_scale_variable
Sp%Dim%nd_continuum = nd_continuum
Sp%Dim%nd_drop_type = nd_drop_type
Sp%Dim%nd_ice_type = nd_ice_type
Sp%Dim%nd_aerosol_species = nd_aerosol_species
Sp%Dim%nd_thermal_coeff = nd_thermal_coeff
Sp%Dim%nd_cloud_parameter = nd_cloud_parameter
Sp%Dim%nd_humidity = nd_humidity
Sp%Dim%nd_aod_wavel = nd_aod_wavel
Sp%Dim%nd_phase_term = nd_phase_term
Sp%Dim%nd_tmp = nd_tmp
Sp%Dim%nd_pre = nd_pre
Sp%Dim%nd_mix = nd_mix
Sp%Dim%nd_band_mix_gas = nd_band_mix_gas
Sp%Dim%nd_sub_band = nd_sub_band
Sp%Dim%nd_times = nd_times

! Allocate spectrum arrays that remain unallocated
CALL allocate_spectrum(Sp)


CONTAINS


SUBROUTINE read_block_int

IMPLICIT NONE

! Local variables
LOGICAL :: l_block_read
!   Flag for correctly read block
CHARACTER (LEN=80) :: line_temp
!   Temporary string to hold the current line
!   (this is only temporary: old spectral files should be processed
!   thoroughly eventually)


! Start from the position that the block has not been read, then
! use the flag to check for errors.
l_block_read= .FALSE.

! Depending on the value of I_TYPE, the appropriate subroutine is called.
IF (i_type == 0) THEN
  IF (i_subtype == 0) THEN
    IF (i_version == 1) THEN
      CALL read_block_0_0_1_int
      l_block_read= .TRUE.
    ELSE IF (i_version == 2) THEN
      CALL read_block_0_0_2_int
      l_block_read= .TRUE.
    END IF
  END IF
ELSE IF (i_type == 1) THEN
  IF (i_subtype == 0) THEN
    IF (i_version == 0) THEN
      CALL read_block_1_0_0_int
      l_block_read= .TRUE.
    END IF
  END IF
ELSE IF (i_type == 2) THEN
  IF (i_subtype == 0) THEN
    IF (i_version == 0) THEN
      CALL read_block_2_0_0_int
      l_block_read= .TRUE.
    ELSE IF (i_version == 1) THEN
      CALL read_block_2_0_1_int
      l_block_read= .TRUE.
    END IF
  ELSE IF (i_subtype == 1) THEN
    IF (i_version == 0) THEN
      CALL read_block_2_1_0_int
      l_block_read= .TRUE.
    END IF
  END IF
ELSE IF (i_type == 3) THEN
  IF (i_subtype == 0) THEN
    IF (i_version == 0) THEN
      CALL read_block_3_0_0_int
      l_block_read= .TRUE.
    END IF
  ELSE IF (i_subtype == 1) THEN
    IF (i_version == 0) THEN
      CALL read_block_3_1_0_int
      l_block_read= .TRUE.
    END IF
  END IF
ELSE IF (i_type == 4) THEN
  IF (i_subtype == 0) THEN
    IF (i_version == 0) THEN
      CALL read_block_4_0_0_int
      l_block_read= .TRUE.
    ELSE IF (i_version == 1) THEN
      CALL read_block_4_0_1_int
      l_block_read= .TRUE.
    END IF
  END IF
ELSE IF (i_type == 5) THEN
  IF (i_subtype == 0) THEN
    IF (i_version == 0) THEN
      CALL read_block_5_0_0_int
      l_block_read= .TRUE.
    ELSE IF (i_version == 1) THEN
      CALL read_block_5_0_1_int
      l_block_read= .TRUE.
    END IF
  ELSE IF (i_subtype == 1) THEN
    IF (i_version == 0) THEN
      CALL read_block_5_1_0_int
      l_block_read= .TRUE.
    END IF
  END IF
ELSE IF (i_type == 6) THEN
  IF (i_subtype == 0) THEN
    IF (i_version == 0) THEN
      CALL read_block_6_0_0_int
      l_block_read= .TRUE.
    ELSE IF (i_version == 1) THEN
      CALL read_block_6_0_1_int
      l_block_read= .TRUE.
    END IF
  END IF
ELSE IF (i_type == 8) THEN
  IF (i_subtype == 0) THEN
    IF (i_version == 0) THEN
      CALL read_block_8_0_0_int
      l_block_read= .TRUE.
    END IF
  END IF
ELSE IF (i_type == 9) THEN
  IF (i_subtype == 0) THEN
    IF (i_version == 0) THEN
      CALL read_block_9_0_0_int
      l_block_read= .TRUE.
    END IF
  ELSE IF (i_subtype == 1) THEN
    IF (i_version == 0) THEN
      CALL read_block_9_1_0_int
      l_block_read= .TRUE.
    END IF
  END IF
ELSE IF (i_type == 10) THEN
  IF (i_subtype == 0) THEN
    IF (i_version == 0) THEN
      CALL read_block_10_0_0_int
      l_block_read= .TRUE.
    ELSE IF (i_version == 1) THEN
      CALL read_block_10_0_1_int
      l_block_read= .TRUE.
    ELSE IF (i_version == 2) THEN
      CALL read_block_10_0_2_int
      l_block_read= .TRUE.
    END IF
  END IF
ELSE IF (i_type == 11) THEN
  IF (i_subtype == 0) THEN
    IF (i_version == 1) THEN
      CALL read_block_11_0_1_int
      l_block_read= .TRUE.
    ELSE IF (i_version == 2) THEN
      CALL read_block_11_0_2_int
      l_block_read= .TRUE.
    END IF
  ELSE IF (i_subtype == 1) THEN
    IF (i_version == 0) THEN
      CALL read_block_11_1_0_int
      l_block_read= .TRUE.
    ELSE IF (i_version == 1) THEN
      CALL read_block_11_1_1_int
      l_block_read= .TRUE.
    ELSE IF (i_version == 2) THEN
      CALL read_block_11_1_2_int
      l_block_read= .TRUE.
    END IF
  END IF
ELSE IF (i_type == 12) THEN
  IF (i_subtype == 0) THEN
    IF (i_version == 0) THEN
      CALL read_block_12_0_0_int
      l_block_read= .TRUE.
    ELSE IF (i_version == 1) THEN
      CALL read_block_12_0_1_int
      l_block_read= .TRUE.
    ELSE IF (i_version == 2) THEN
      CALL read_block_12_0_2_int
      l_block_read= .TRUE.
    END IF
  END IF
ELSE IF (i_type == 13) THEN
  IF (i_subtype == 0) THEN
    IF (i_version == 0) THEN
      CALL read_block_13_0_0_int
      l_block_read= .TRUE.
    END IF
  END IF
ELSE IF (i_type == 14) THEN
  IF (i_subtype == 0) THEN
    IF (i_version == 0) THEN
      CALL read_block_14_0_0_int
      l_block_read= .TRUE.
    END IF
  END IF
ELSE IF (i_type == 15) THEN
  IF (i_subtype == 0) THEN
    CALL read_block_15_0_0
    l_block_read= .TRUE.
  ELSE IF (i_subtype == 1 ) THEN
    CALL read_block_15_1_0
    l_block_read= .TRUE.
  END IF
ELSE IF (i_type == 16) THEN
  IF (i_subtype == 0) THEN
    IF (i_version == 0) THEN
      CALL read_block_16_0_0_int
      l_block_read= .TRUE.
    END IF
  END IF
ELSE IF (i_type == 17) THEN
  IF (i_subtype == 0) THEN
    IF (i_version == 0) THEN
      CALL read_block_17_0_0
      l_block_read= .TRUE.
    END IF
  END IF
END IF

IF (ierr /= i_normal) THEN
  l_block_read = .FALSE.
  RETURN
END IF

IF (.NOT. l_block_read) THEN
  ! This block is not of a supported type.
  WRITE(cmessage, '(a43, i3, a12, i3, a11, i3)') &
    'This sort of block is not permitted. Type =', i_type, &
    ': Sub-type =', i_subtype, &
    ': Version =', i_version
  ierr=i_warning
  CALL ereport(RoutineName, ierr, cmessage)
  ! Read through to the end of the block.
  DO
    READ(iu_spc, '(a)') line_temp
    IF (line_temp(1:4) == '*END') THEN
      BACKSPACE(iu_spc)
      EXIT
    END IF
  END DO
  RETURN
END IF

END SUBROUTINE read_block_int



SUBROUTINE read_block_0_0_1_int

IMPLICIT NONE

! Local variables.
CHARACTER :: chdum
!   Dummy character
INTEGER :: idum
!   Dummy integer


! Skip over the header.
READ(iu_spc, *)

! Read in the number of spectral intervals, the number of
! gaseous absorbers and the number of aerosols.
READ(iu_spc, '(27x, i5)', IOSTAT=ios, IOMSG=iomessage) Sp%Basic%n_band
IF (ios /= 0) THEN
  cmessage = '*** Error in subroutine read_block_0_0_1.\n' // &
    'Number of bands could not be read: ' // TRIM(iomessage)
  ierr=i_err_fatal
  RETURN
END IF
nd_band=Sp%Basic%n_band
ALLOCATE(Sp%Basic%n_band_exclude(nd_band))
! This must be zeroed lest block 14 should not be present and
! the array be filled with random values.
Sp%Basic%n_band_exclude(1:nd_band) = 0

READ(iu_spc, '(36x, i5)', IOSTAT=ios, IOMSG=iomessage) Sp%Gas%n_absorb
IF (ios /= 0) THEN
  cmessage = '*** Error in subroutine read_block_0_0_1.\n' // &
    'Number of absorbers could not be read: '// TRIM(iomessage)
  ierr=i_err_fatal
  RETURN
END IF
nd_species=Sp%Gas%n_absorb

READ(iu_spc, '(27x, i5)', IOSTAT=ios, IOMSG=iomessage) Sp%Aerosol%n_aerosol
IF (ios /= 0) THEN
  cmessage = '*** Error in subroutine read_block_0_0_1. ' // &
    'Number of aerosols could not be read: '// TRIM(iomessage)
  ierr=i_err_fatal
  RETURN
END IF
nd_aerosol_species=MAX(Sp%Aerosol%n_aerosol, 1)


! Read over the headers and the list of absorbers.
ALLOCATE(Sp%Gas%type_absorb(nd_species))
READ(iu_spc, '(/)')
DO i=1, Sp%Gas%n_absorb
  READ(iu_spc, '(i5, 7x, i5, 7x, a)') &
    idum, Sp%Gas%type_absorb(i), chdum
END DO

! Read over the headers and the list of aerosols.
ALLOCATE(Sp%Aerosol%type_aerosol(nd_aerosol_species))
READ(iu_spc, '(/)')
DO i=1, Sp%Aerosol%n_aerosol
  READ(iu_spc, '(i5, 7x, i5, 7x, a)') &
    idum, Sp%Aerosol%type_aerosol(i), chdum
END DO

END SUBROUTINE read_block_0_0_1_int



SUBROUTINE read_block_0_0_2_int

IMPLICIT NONE

! Local variables.
INTEGER :: desc_end
!   Position of equals sign to delimit end of item description
CHARACTER :: chdum
!   Dummy character
INTEGER :: idum
!   Dummy integer

DO
  READ(iu_spc, '(a80)', IOSTAT=ios) line
  IF (line(1:4) == '*END') THEN
    BACKSPACE(iu_spc)
    EXIT
  END IF

  desc_end=SCAN(line,'=',.TRUE.)
  IF (desc_end==0) desc_end=LEN_TRIM(line)+1
  SELECT CASE (line(1:desc_end-1))

  CASE ('Number of spectral bands','nd_band')
    READ(line(desc_end+1:),*,IOSTAT=ios, IOMSG=iomessage) Sp%Basic%n_band
    nd_band=Sp%Basic%n_band
  CASE ('Total number of gaseous absorbers','nd_species')
    READ(line(desc_end+1:),*,IOSTAT=ios, IOMSG=iomessage) Sp%Gas%n_absorb
    nd_species=Sp%Gas%n_absorb
  CASE ('Total number of aerosols','nd_aerosol_species')
    READ(line(desc_end+1:),*,IOSTAT=ios, IOMSG=iomessage) &
                                                     Sp%Aerosol%n_aerosol
    nd_aerosol_species=MAX(Sp%Aerosol%n_aerosol, 1)
  CASE ('List of indexing numbers and absorbers.')
    ALLOCATE(Sp%Gas%type_absorb(nd_species))
    READ(iu_spc, *)
    DO i=1, Sp%Gas%n_absorb
      READ(iu_spc, '(i5, 7x, i5, 7x, a)') &
        idum, Sp%Gas%type_absorb(i), chdum
    END DO
  CASE ('List of indexing numbers of aerosols.')
    ALLOCATE(Sp%Aerosol%type_aerosol(nd_aerosol_species))
    READ(iu_spc, *)
    DO i=1, Sp%Aerosol%n_aerosol
      READ(iu_spc, '(i5, 7x, i5, 7x, a)') &
        idum, Sp%Aerosol%type_aerosol(i), chdum
    END DO
  CASE ('Maximum number of k-terms in a band','nd_k_term')
    READ(line(desc_end+1:),*,IOSTAT=ios, IOMSG=iomessage) nd_k_term

  END SELECT
  IF (ios /= 0) THEN
    cmessage = 'Error in subroutine read_block_0_0_2_int: ' //           &
               TRIM(iomessage)
    ierr=i_err_fatal
    RETURN
  END IF

END DO

ALLOCATE(Sp%Basic%n_band_exclude(nd_band))
! This must be zeroed lest block 14 should not be present and
! the array be filled with random values.
Sp%Basic%n_band_exclude(1:nd_band) = 0

END SUBROUTINE read_block_0_0_2_int



SUBROUTINE read_block_1_0_0_int

IMPLICIT NONE

! Local variables.
INTEGER :: idum
!   Dummy integer


! Skip over the headers.
READ(iu_spc, '(//)')

! Read in the limits on the intervals in the spectrum
ALLOCATE(Sp%Basic%wavelength_short(nd_band))
ALLOCATE(Sp%Basic%wavelength_long(nd_band))
DO i=1, Sp%Basic%n_band
  READ(iu_spc, FMT='(i5, 7x, 1pe16.9, 4x, 1pe16.9)', IOSTAT=ios,         &
       IOMSG=iomessage) &
    idum, Sp%Basic%wavelength_short(i), &
          Sp%Basic%wavelength_long(i)
  IF (ios /= 0) THEN
    cmessage = '*** Error in subroutine read_block_1_0_0.\n' // &
      'Wavelength limits of bands could not be read: ' // TRIM(iomessage)
    ierr=i_err_fatal
    RETURN
  END IF
END DO

END SUBROUTINE read_block_1_0_0_int



SUBROUTINE read_block_2_0_0_int

IMPLICIT NONE

! Local variables.
INTEGER :: idum
!   Dummy integer


! Skip over the headers.
READ(iu_spc, '(/)')

! Read in the limits on the intervals in the spectrum
ALLOCATE(Sp%Solar%solar_flux_band(nd_band))
DO i=1, Sp%Basic%n_band
  READ(iu_spc, FMT='(i5, 7x, 1pe16.9)', IOSTAT=ios, IOMSG=iomessage) &
    idum, Sp%Solar%solar_flux_band(i)
  IF (ios /= 0) THEN
    cmessage = '*** Error: solar spectral data are not correct: ' // &
               TRIM(iomessage)
    ierr=i_err_fatal
    RETURN
  END IF
END DO

END SUBROUTINE read_block_2_0_0_int



SUBROUTINE read_block_2_0_1_int

IMPLICIT NONE

! Local variables.
INTEGER :: idum
!   Dummy integer


! Skip over the headers.
READ(iu_spc, '(/)')

! Read in the limits on the intervals in the spectrum
ALLOCATE(Sp%Solar%solar_flux_band(nd_band))
ALLOCATE(Sp%Solar%weight_blue(nd_band))
DO i=1, Sp%Basic%n_band
  READ(iu_spc, FMT='(i5, 7x, 1pe16.9, 4x, 1pe16.9)', IOSTAT=ios,         &
       IOMSG=iomessage) &
    idum, Sp%Solar%solar_flux_band(i), &
    Sp%Solar%weight_blue(i)
  IF (ios /= 0) THEN
    cmessage = '*** Error: solar spectral data are not correct: '//  &
               TRIM(iomessage)
    ierr=i_err_fatal
    RETURN
  END IF
END DO

END SUBROUTINE read_block_2_0_1_int



SUBROUTINE read_block_2_1_0_int

IMPLICIT NONE

! Local variables.
INTEGER :: iss, k
!   Loop variables
INTEGER :: i_band
!   Number of band
INTEGER :: n_term
!   Number of terms
INTEGER :: i_solar_src = 3
!   Solar spectrum to use (hardwired here for now)

ALLOCATE( Sp%Solar%solar_flux_band_ses(nd_k_term, nd_band) )

! Read up to the solar spectrum to be used
DO iss=1, i_solar_src
  ! Skip over the headers.
  READ(iu_spc, '(/)')

  ! Read in the limits on the intervals in the spectrum
  DO i=1, Sp%Basic%n_band
    READ(iu_spc, FMT='(I5, 7X, I5)', IOSTAT=ios) i_band, n_term
    READ(iu_spc, FMT='(4(1X,1PE16.9))',IOSTAT=ios, IOMSG=iomessage) &
      (Sp%Solar%solar_flux_band_ses(k, i),k=1, n_term)
    IF (ios /= 0) THEN
      cmessage = '*** Error: solar spectral data are not correct: ' &
                 // TRIM(iomessage)
      ierr=i_err_fatal
      RETURN
    END IF
  END DO
END DO
! Read through to the end of the block.
DO
  READ(iu_spc, '(a)') char_dum
  IF (char_dum(1:4) == '*END') THEN
    BACKSPACE(iu_spc)
    EXIT
  END IF
END DO

END SUBROUTINE read_block_2_1_0_int



SUBROUTINE read_block_3_0_0_int

IMPLICIT NONE

! Local variables.
INTEGER :: idum
!   Dummy integer


! Skip over the headers.
READ(iu_spc, '(//)')

! Rayleigh scattering coefficients are for the total gas
Sp%Rayleigh%i_rayleigh_scheme = ip_rayleigh_total

! Read in the limits on the intervals in the spectrum
ALLOCATE(Sp%Rayleigh%rayleigh_coeff(nd_band))
DO i=1, Sp%Basic%n_band
  READ(iu_spc, FMT='(i5, 7x, 1pe16.9)', IOSTAT=ios, IOMSG=iomessage) &
    idum, Sp%Rayleigh%rayleigh_coeff(i)
  IF (ios /= 0) THEN
    cmessage = '*** Error: rayleigh scattering data are not correct: ' // &
               TRIM(iomessage)
    ierr=i_err_fatal
    RETURN
  END IF
END DO

END SUBROUTINE read_block_3_0_0_int


SUBROUTINE read_block_3_1_0_int

IMPLICIT NONE


! Skip over the headers.
READ(iu_spc, *)

! Rayleigh scattering coefficients are tabulated for each gas
Sp%Rayleigh%i_rayleigh_scheme = ip_rayleigh_custom

! Get number of gases for which Rayleigh scattering coefficients are tabulated
READ(iu_spc, '(38X, I5)') Sp%Rayleigh%n_gas_rayleigh
IF (ios /= 0) THEN
  cmessage = '*** Error: rayleigh scattering data are not correct: ' // &
             TRIM(iomessage)
  ierr=i_err_fatal
  RETURN
END IF

! Read indexing numbers of gases
ALLOCATE(Sp%Rayleigh%index_rayleigh(nd_species))
READ(iu_spc, FMT='(/, 15(2X, I3))', IOSTAT=ios, IOMSG=iomessage) &
  Sp%Rayleigh%index_rayleigh(1:Sp%Rayleigh%n_gas_rayleigh)
IF (ios /= 0) THEN
  cmessage = '*** Error: rayleigh scattering data are not correct: ' // &
             TRIM(iomessage)
  ierr=i_err_fatal
  RETURN
END IF

! Read in the limits on the intervals in the spectrum
ALLOCATE(Sp%Rayleigh%rayleigh_coeff_gas(nd_species,nd_band))
DO i=1, Sp%Basic%n_band
  READ(iu_spc, FMT='(/, (4(3X, 1PE16.9)))', IOSTAT=ios, IOMSG=iomessage) &
    Sp%Rayleigh%rayleigh_coeff_gas(1:Sp%Rayleigh%n_gas_rayleigh,i)
  IF (ios /= 0) THEN
    cmessage = '*** Error: rayleigh scattering data are not correct: ' // &
               TRIM(iomessage)
    ierr=i_err_fatal
    RETURN
  END IF
END DO

END SUBROUTINE read_block_3_1_0_int



SUBROUTINE read_block_4_0_0_int

IMPLICIT NONE

! Local variables.
INTEGER :: idum
!   Dummy integer
INTEGER :: j
!   Loop variable


! Skip over the headers.
READ(iu_spc, '(////)')

! Read in the list of absorbers in each band.
ALLOCATE(Sp%Gas%n_band_absorb(nd_band))
ALLOCATE(Sp%Gas%index_absorb(nd_species, nd_band))
DO i=1, Sp%Basic%n_band
  READ(iu_spc, FMT='(i5, 7x, i5)', IOSTAT=ios, IOMSG=iomessage) &
    idum, Sp%Gas%n_band_absorb(i)
  IF (ios /= 0) THEN
    cmessage = '*** Error in subroutine read_block_4_0_0.\n' // &
      'The list of absorbers is not correct: '// TRIM(iomessage)
    ierr=i_err_fatal
    RETURN
  END IF
  IF (Sp%Gas%n_band_absorb(i) > 0) THEN
    READ(iu_spc, '(5x, 4(2x, i3))', IOSTAT=ios, IOMSG=iomessage) &
      ( Sp%Gas%index_absorb(j, i), &
        j=1, Sp%Gas%n_band_absorb(i) )
  END IF
  IF (ios /= 0) THEN
    cmessage = '*** Error in subroutine read_block_4_0_0. ' // &
      'The index of absorbers is not correct: ' // TRIM(iomessage)
    ierr=i_err_fatal
    RETURN
  END IF
END DO

END SUBROUTINE read_block_4_0_0_int



SUBROUTINE read_block_4_0_1_int

IMPLICIT NONE

! Local variables.
INTEGER :: idum
!   Dummy integer
INTEGER :: j
!   Loop variable


! Skip over the headers.
READ(iu_spc, '(////)')

! Read in the list of absorbers in each band.
ALLOCATE(Sp%Gas%n_band_absorb(nd_band))
ALLOCATE(Sp%Gas%index_absorb(nd_species, nd_band))
DO i=1, Sp%Basic%n_band
  READ(iu_spc, FMT='(i5, 7x, i5)', IOSTAT=ios, IOMSG=iomessage) &
    idum, Sp%Gas%n_band_absorb(i)
  IF (ios /= 0) THEN
    cmessage = '*** Error in subroutine read_block_4_0_1. ' // &
      'The list of absorbers is not correct: '// TRIM(iomessage)
    ierr=i_err_fatal
    RETURN
  END IF
  IF (Sp%Gas%n_band_absorb(i) > 0) THEN
    READ(iu_spc, '(5x, 8(2x, i3))', IOSTAT=ios, IOMSG=iomessage) &
      ( Sp%Gas%index_absorb(j, i), &
        j=1, Sp%Gas%n_band_absorb(i) )
  END IF
  IF (ios /= 0) THEN
    cmessage = '*** Error in subroutine read_block_4_0_1. ' // &
      'The index of absorbers is not correct: ' // TRIM(iomessage)
    ierr=i_err_fatal
    RETURN
  END IF
END DO

END SUBROUTINE read_block_4_0_1_int



SUBROUTINE read_block_5_0_0_int

IMPLICIT NONE

! Local variables.
INTEGER :: idum_band
!   Dummy integer
INTEGER :: idum_species, jspecies(nd_species, nd_band)
!   Dummy integer
INTEGER :: idum_scale
!   Dummy integer
INTEGER :: idum_fnc
!   Dummy integer
INTEGER :: number_term
!   Number of ESFT terms
INTEGER :: j, i_term, l, it, ip
!   Loop variables
LOGICAL :: l_lookup


! Allocate space for the arrays of k-terms.
ALLOCATE(Sp%Gas%i_band_k(nd_band, nd_species))
ALLOCATE(Sp%Gas%i_scale_k(nd_band, nd_species))
ALLOCATE(Sp%Gas%i_scale_fnc(nd_band, nd_species))
Sp%Gas%i_scale_fnc=0.0_RealK
ALLOCATE(Sp%Gas%k(nd_k_term, nd_band, nd_species))
ALLOCATE(Sp%Gas%w(nd_k_term, nd_band, nd_species))
ALLOCATE(Sp%Gas%p_ref(nd_species, nd_band))
ALLOCATE(Sp%Gas%t_ref(nd_species, nd_band))
ALLOCATE(Sp%Gas%scale(nd_scale_variable, nd_k_term, &
  nd_band, nd_species))
ALLOCATE(Sp%Gas%i_scat(nd_k_term, nd_band, nd_species))
ALLOCATE(Sp%Gas%num_ref_p(nd_species, nd_band))
ALLOCATE(Sp%Gas%num_ref_t(nd_species, nd_band))
Sp%Gas%num_ref_p=0
Sp%Gas%num_ref_t=0
l_lookup=.FALSE.

! Skip over the headers.
READ(iu_spc, '(///)')

! Read in the number of k-terms in each band.
DO i=1, Sp%Basic%n_band
  DO j=1, Sp%Gas%n_band_absorb(i)
    READ(iu_spc, FMT='(i5, 5(7x, i5))', IOSTAT=ios, IOMSG=iomessage) &
      idum_band, idum_species, number_term, idum_scale, idum_fnc
    IF (ios /= 0) THEN
      cmessage = '*** Error in subroutine read_block_5_0_0.\n' //       &
        'k-distribution data are not consistent with the summary: ' // &
        TRIM(iomessage)
      ierr=i_err_fatal
      RETURN
    END IF
    jspecies(j,i) = idum_species
    IF ( (idum_fnc /= IP_scale_power_law)    .AND. &
         (idum_fnc /= IP_scale_power_quad)   .AND. &
         (idum_fnc /= IP_scale_doppler_quad) .AND. &
         (idum_fnc /= IP_scale_wenyi)        .AND. &
         (idum_fnc /= IP_scale_dbl_pow_law)  .AND. &
         (idum_fnc /= IP_scale_dbl_pow_quad) .AND. &
         (idum_fnc /= IP_scale_dbl_dop_quad) .AND. &
         (idum_fnc /= IP_scale_lookup)       .AND. &
         (idum_fnc /= IP_scale_fnc_null) ) THEN
      cmessage = '*** Error in subroutine read_block_5_0_0. ' // &
        'an illegal scaling function has been specified.'
      ierr=i_err_fatal
      RETURN
    END IF
    IF (number_term > nd_k_term) THEN
      cmessage = '*** Error in subroutine read_block_5_0_0. ' // &
        'Too many esft terms have been given. ' // &
        'Increase npd_k_term and recompile.'
      ierr=i_err_fatal
      RETURN
    END IF
    Sp%Gas%i_band_k(idum_band, idum_species)=number_term
    Sp%Gas%i_scale_k(idum_band, idum_species)=idum_scale
    Sp%Gas%i_scale_fnc(idum_band, idum_species)=idum_fnc

    ! Read the reference temperature and pressure.
    IF (idum_fnc == ip_scale_lookup) THEN
      IF (.NOT. l_exist_k) THEN
        cmessage = 'Extended spectral file cannot be opened, ' // &
          'k-term lookup table for BLOCK 5 cannot be read.'
        ierr=i_err_fatal
        RETURN
      END IF
      READ(iu_spc, '(2(7x, i5))') &
        Sp%Gas%num_ref_p(idum_species, idum_band), &
        Sp%Gas%num_ref_t(idum_species, idum_band)
      l_lookup=.TRUE.
    ELSE
      READ(iu_spc, '(2(6x, 1pe16.9))') &
        Sp%Gas%p_ref(idum_species, idum_band), &
        Sp%Gas%t_ref(idum_species, idum_band)
    END IF
    ! For each band read in the k-terms and weights.
    DO i_term=1, number_term
      READ(iu_spc, '(2(3x, 1pe16.9),i3,:,(t42, 1pe16.9,3x,1pe16.9))' &
        , IOSTAT=ios, IOMSG=iomessage) &
          Sp%Gas%k(i_term, idum_band, idum_species), &
          Sp%Gas%w(i_term, idum_band, idum_species), &
          Sp%Gas%i_scat(i_term, idum_band, idum_species), &
          (Sp%Gas%scale(l, i_term, idum_band, idum_species), &
          l=1, n_scale_variable(idum_fnc))
      IF (ios /= 0) THEN
        cmessage = '*** Error in subroutine read_block_5_0_0.\n' // &
          'k-distribution data are not consistent with the summary: ' // &
          TRIM(iomessage)
        ierr=i_err_fatal
        RETURN
      END IF
    END DO

  END DO
END DO

IF (l_lookup) THEN
  nd_pre = MAXVAL(Sp%Gas%num_ref_p)
  nd_tmp = MAXVAL(Sp%Gas%num_ref_t)

  ALLOCATE(Sp%Gas%p_lookup( nd_pre ))
  ALLOCATE(Sp%Gas%t_lookup( nd_tmp, nd_pre ))
  ALLOCATE(Sp%Gas%k_lookup( nd_tmp, nd_pre, &
                            nd_k_term, nd_species, nd_band ))

  ! Skip over the headers.
  READ(iu_spc1, '(/)')
  DO ip=1, nd_pre
    READ(iu_spc1, '(6(1PE13.6))', IOSTAT=ios, IOMSG=iomessage) &
      Sp%Gas%p_lookup(ip), (Sp%Gas%t_lookup(it, ip), it=1, nd_tmp)
    IF (ios /= 0) THEN
      WRITE(cmessage,'(a, i4, A)') &
        '*** Error in subroutine read_block_5_0_0.\n' // &
        'P/T look-up table entry:', ip, 'Error: ' // TRIM(iomessage)
      ierr=i_err_fatal
      RETURN
    END IF
    Sp%Gas%p_lookup(ip)=LOG(Sp%Gas%p_lookup(ip))
  END DO
  DO i=1, Sp%Basic%n_band
    DO j=1, Sp%Gas%n_band_absorb(i)
      idum_species = jspecies(j,i)
      IF ( Sp%Gas%i_scale_fnc(i, idum_species) == &
           ip_scale_lookup ) THEN
        IF (Sp%Gas%num_ref_p(idum_species, i) /= nd_pre .OR. &
            Sp%Gas%num_ref_t(idum_species, i) /= nd_tmp) THEN
          cmessage = '*** Error in subroutine read_block_5_0_0: '// &
            'P/T lookup table size is not consistent.'
          ierr=i_err_fatal
          RETURN
        END IF
        ! Skip over the headers.
        READ(iu_spc1, '(/)')
        DO i_term=1, Sp%Gas%i_band_k(i, idum_species)
          DO ip=1, nd_pre
            READ(iu_spc1, '(6(1PE13.6))', IOSTAT=ios, IOMSG=iomessage) &
              (Sp%Gas%k_lookup(it,ip,i_term,idum_species,i), &
               it=1, nd_tmp)
            IF (ios /= 0) THEN
              WRITE(cmessage,'(a, 4i4, A)') &
                '*** Error in subroutine read_block_5_0_0:\n' // &
                'Look-up table entry:', i, i_term, idum_species, ip, &
                '\nError: ' // TRIM(iomessage)
              ierr=i_err_fatal
              RETURN
            END IF
          END DO
        END DO
      END IF
    END DO
  END DO
END IF

END SUBROUTINE read_block_5_0_0_int



SUBROUTINE read_block_5_0_1_int

IMPLICIT NONE

! Local variables.
INTEGER :: idum_band
!   Dummy integer
INTEGER :: idum_species, jspecies(nd_species, nd_band)
!   Dummy integer
INTEGER :: idum_scale
!   Dummy integer
INTEGER :: idum_fnc
!   Dummy integer
INTEGER :: number_term
!   Number of ESFT terms
INTEGER :: j, i_term, l, it, ip
!   Loop variables
LOGICAL :: l_lookup
!   True if a k-table is used
LOGICAL :: l_k_table_exists
!   True if the k-table is present in the extended spectral file


! Allocate space for the arrays of k-terms.
ALLOCATE(Sp%Gas%i_band_k(nd_band, nd_species))
ALLOCATE(Sp%Gas%i_scale_k(nd_band, nd_species))
ALLOCATE(Sp%Gas%i_scale_fnc(nd_band, nd_species))
Sp%Gas%i_scale_fnc=0.0_RealK
ALLOCATE(Sp%Gas%k(nd_k_term, nd_band, nd_species))
ALLOCATE(Sp%Gas%w(nd_k_term, nd_band, nd_species))
ALLOCATE(Sp%Gas%p_ref(nd_species, nd_band))
ALLOCATE(Sp%Gas%t_ref(nd_species, nd_band))
ALLOCATE(Sp%Gas%scale(nd_scale_variable, nd_k_term, &
  nd_band, nd_species))
ALLOCATE(Sp%Gas%i_scat(nd_k_term, nd_band, nd_species))
ALLOCATE(Sp%Gas%num_ref_p(nd_species, nd_band))
ALLOCATE(Sp%Gas%num_ref_t(nd_species, nd_band))
Sp%Gas%num_ref_p=0
Sp%Gas%num_ref_t=0
l_lookup=.FALSE.

! Skip over the headers.
READ(iu_spc, '(///)')

! Read in the number of k-terms in each band.
DO i=1, Sp%Basic%n_band
  DO j=1, Sp%Gas%n_band_absorb(i)
    READ(iu_spc, FMT='(i5, 5(7x, i5))', IOSTAT=ios, IOMSG=iomessage) &
      idum_band, idum_species, number_term, idum_scale, idum_fnc
    IF (ios /= 0) THEN
      cmessage = '*** Error in subroutine read_block_5_0_1.\n' // &
        'k-distribution data are not consistent with the summary: ' // &
        '\nError: ' // TRIM(iomessage)
      ierr=i_err_fatal
      RETURN
    END IF
    jspecies(j,i) = idum_species
    IF ( (idum_fnc /= IP_scale_power_law)    .AND. &
         (idum_fnc /= IP_scale_power_quad)   .AND. &
         (idum_fnc /= IP_scale_doppler_quad) .AND. &
         (idum_fnc /= IP_scale_wenyi)        .AND. &
         (idum_fnc /= IP_scale_dbl_pow_law)  .AND. &
         (idum_fnc /= IP_scale_dbl_pow_quad) .AND. &
         (idum_fnc /= IP_scale_dbl_dop_quad) .AND. &
         (idum_fnc /= IP_scale_lookup)       .AND. &
         (idum_fnc /= IP_scale_fnc_null) ) THEN
      cmessage = '*** Error in subroutine read_block_5_0_1. ' // &
        'an illegal scaling function has been specified.'
      ierr=i_err_fatal
      RETURN
    END IF
    IF (number_term > nd_k_term) THEN
      cmessage = '*** Error in subroutine read_block_5_0_1. ' // &
        'Too many esft terms have been given. ' // &
        'Increase npd_k_term and recompile.'
      ierr=i_err_fatal
      RETURN
    END IF
    Sp%Gas%i_band_k(idum_band, idum_species)=number_term
    Sp%Gas%i_scale_k(idum_band, idum_species)=idum_scale
    Sp%Gas%i_scale_fnc(idum_band, idum_species)=idum_fnc

    ! Read the reference temperature and pressure.
    IF (idum_fnc == ip_scale_lookup) THEN
      IF (.NOT. l_exist_k) THEN
        cmessage = 'Extended spectral file cannot be opened, ' // &
          'k-term lookup table for BLOCK 5 cannot be read.'
        ierr=i_err_fatal
        RETURN
      END IF
      READ(iu_spc, '(2(7x, i5))') &
        Sp%Gas%num_ref_p(idum_species, idum_band), &
        Sp%Gas%num_ref_t(idum_species, idum_band)
      l_lookup=.TRUE.
    ELSE
      READ(iu_spc, '(2(6x, 1pe16.9))') &
        Sp%Gas%p_ref(idum_species, idum_band), &
        Sp%Gas%t_ref(idum_species, idum_band)
    END IF
    ! For each band read in the k-terms and weights.
    DO i_term=1, number_term
      READ(iu_spc, '(2(3x, 1pe16.9),i3,:,(t42, 1pe16.9,3x,1pe16.9))' &
        , IOSTAT=ios, IOMSG=iomessage) &
          Sp%Gas%k(i_term, idum_band, idum_species), &
          Sp%Gas%w(i_term, idum_band, idum_species), &
          Sp%Gas%i_scat(i_term, idum_band, idum_species), &
          (Sp%Gas%scale(l, i_term, idum_band, idum_species), &
          l=1, n_scale_variable(idum_fnc))
      IF (ios /= 0) THEN
        cmessage = '*** Error in subroutine read_block_5_0_1.\n' // &
          'k-distribution data are not consistent with the summary: ' // &
          '\nError: ' // TRIM(iomessage)
        ierr=i_err_fatal
        RETURN
      END IF
    END DO

  END DO
END DO

IF (l_lookup) THEN
  nd_pre = MAXVAL(Sp%Gas%num_ref_p)
  nd_tmp = MAXVAL(Sp%Gas%num_ref_t)

  ALLOCATE(Sp%Gas%p_lookup( nd_pre ))
  ALLOCATE(Sp%Gas%t_lookup( nd_tmp, nd_pre ))
  ALLOCATE(Sp%Gas%k_lookup( nd_tmp, nd_pre, &
                            nd_k_term, nd_species, nd_band ))

  ! Locate correct block in extended spectral file
  l_k_table_exists=.FALSE.
  REWIND(iu_spc1)
  DO
    READ(iu_spc1, '(a80)') char_dum
    IF ( char_dum(1:6) == '*BLOCK' ) THEN
      IF ( char_dum(9:16) == 'k-table' ) THEN
        l_k_table_exists=.TRUE.
        EXIT
      END IF
    END IF
  END DO

  ! Return with error if the k-table does not exist
  IF (.NOT. l_k_table_exists) THEN
    cmessage = 'No k-table in extended spectral file.'
    ierr = i_err_fatal
    RETURN
  END IF

  ! Skip over the headers.
  READ(iu_spc1, '(/)')
  DO ip=1, nd_pre
    READ(iu_spc1, '(6(1PE13.6))', IOSTAT=ios, IOMSG=iomessage) &
      Sp%Gas%p_lookup(ip), (Sp%Gas%t_lookup(it, ip), it=1, nd_tmp)
    IF (ios /= 0) THEN
      WRITE(cmessage,'(a, i4, A)') &
        '*** Error in subroutine read_block_5_0_1.\n' // &
        'P/T look-up table entry:', ip, '\nError: ' // TRIM(iomessage)
      ierr=i_err_fatal
      RETURN
    END IF
    Sp%Gas%p_lookup(ip)=LOG(Sp%Gas%p_lookup(ip))
  END DO
  DO i=1, Sp%Basic%n_band
    DO j=1, Sp%Gas%n_band_absorb(i)
      idum_species = jspecies(j,i)
      IF ( Sp%Gas%i_scale_fnc(i, idum_species) == &
           ip_scale_lookup ) THEN
        IF (Sp%Gas%num_ref_p(idum_species, i) /= nd_pre .OR. &
            Sp%Gas%num_ref_t(idum_species, i) /= nd_tmp) THEN
          cmessage = '*** Error in subroutine read_block_5_0_1: '// &
            'P/T lookup table size is not consistent.'
          ierr=i_err_fatal
          RETURN
        END IF
        ! Skip over the headers.
        READ(iu_spc1, '(/)')
        DO i_term=1, Sp%Gas%i_band_k(i, idum_species)
          DO ip=1, nd_pre
            READ(iu_spc1, '(6(1PE13.6))', IOSTAT=ios, IOMSG=iomessage) &
              (Sp%Gas%k_lookup(it,ip,i_term,idum_species,i), &
               it=1, nd_tmp)
            IF (ios /= 0) THEN
              WRITE(cmessage,'(a, 4i4, A)') &
                '*** Error in subroutine read_block_5_0_1:\n' // &
                'Look-up table entry:', i, i_term, idum_species, ip, &
                '\nError: ' // TRIM(iomessage)
              ierr=i_err_fatal
              RETURN
            END IF
          END DO
        END DO
      END IF
    END DO
  END DO
END IF

END SUBROUTINE read_block_5_0_1_int



SUBROUTINE read_block_5_1_0_int

IMPLICIT NONE

! Local variables.
INTEGER :: idum_band
!   Dummy integer
INTEGER :: i_band_k_ses
!   Number of k-terms in band
INTEGER :: num_ref_p
!   Number of refrence pressures
INTEGER :: num_mix
!   Number of binary parameter for interpolation of absorption
!   coefficient for mixture of two species
REAL (RealK) :: f_mix
!   Mixing ratio of mixed absorber amount
INTEGER :: k, it, ip, i_gas, i_term
!   Loop variables


IF (.NOT. l_exist_k) THEN
  cmessage = 'Extended spectral file cannot be opened, ' // &
    'k-terms for BLOCK 5 cannot be read.'
  ierr=i_err_fatal
  RETURN
END IF

! Allocate space for the arrays of k-terms.
ALLOCATE(Sp%Gas%i_band_k_ses(nd_band))
ALLOCATE(Sp%Gas%num_ref_p(1, nd_band))
ALLOCATE(Sp%Gas%num_mix(nd_band))
ALLOCATE(Sp%Gas%f_mix(nd_band))
ALLOCATE(Sp%Gas%w_ses(nd_k_term, nd_band))
ALLOCATE(Sp%Gas%i_scale_fnc(nd_band, nd_species))
ALLOCATE(Sp%Gas%i_scale_k(nd_band, nd_species))

Sp%Gas%i_scale_fnc = ip_scale_ses2
Sp%Gas%i_scale_k   = ip_scale_null

READ(iu_spc, '(/)')

! Read in the number of k-terms in each band.
DO i=1, Sp%Basic%n_band
  READ(iu_spc, FMT='(4I5, 1x, 1PE16.9 )', IOSTAT=ios, IOMSG=iomessage) &
    idum_band, i_band_k_ses, num_ref_p, num_mix, f_mix
  IF (ios /= 0) THEN
    cmessage = '*** 1st error in subroutine read_block_5_1_0.\n'// &
      ' ESFT data is not consistent with summary: ' // TRIM(iomessage)
    ierr=i_err_fatal
    RETURN
  END IF
  Sp%Gas%i_band_k_ses(idum_band) = i_band_k_ses
  Sp%Gas%num_ref_p(1, idum_band) = num_ref_p
  Sp%Gas%num_mix(idum_band) = num_mix
  Sp%Gas%f_mix(idum_band) = f_mix

  ! For each band read the values of the coefficients.
  READ(iu_spc, '(4(1X,1PE16.9))', IOSTAT=ios, IOMSG=iomessage) &
    (Sp%Gas%w_ses(k, idum_band),k=1, i_band_k_ses)
  IF (ios /= 0) THEN
    cmessage = '*** 2nd error in subroutine read_block_5_1_0.\n'// &
      ' ESFT data is not consistent with summary: ' // TRIM(iomessage)
    ierr=i_err_fatal
    RETURN
  END IF
END DO
nd_pre = MAXVAL(Sp%Gas%num_ref_p(1, 1:Sp%Basic%n_band))
nd_tmp = 5  ! Currently hardwired
nd_mix = MAXVAL(Sp%Gas%num_mix(1:Sp%Basic%n_band))
nd_band_mix_gas = COUNT( Sp%Gas%num_mix(1:Sp%Basic%n_band) > 1 )

ALLOCATE(Sp%Gas%k_lookup( nd_tmp, nd_pre, &
                          nd_k_term, nd_species, nd_band ))

DO i=1, Sp%Basic%n_band

  ! Skip over the headers.
  READ(iu_spc1, '(/)')

  DO i_gas=1, Sp%Gas%n_band_absorb(i)
    DO i_term=1, Sp%Gas%i_band_k_ses(i)
      DO ip=1, Sp%Gas%num_ref_p(1, i)
        READ(iu_spc1, '(5(1PE12.6,1x))', IOSTAT=ios, IOMSG=iomessage) &
          (Sp%Gas%k_lookup(it,ip,i_term,i_gas,i), it=1, nd_tmp)
        IF (ios /= 0) THEN
          WRITE(cmessage, '(a, 4i4, A)') &
            '*** 3rd error in subroutine read_block_5_1_0: ', &
            i,i_term,i_gas,ip, TRIM(iomessage)
          ierr=i_err_fatal
          RETURN
        END IF
      END DO
    END DO
  END DO

END DO

END SUBROUTINE read_block_5_1_0_int



SUBROUTINE read_block_6_0_0_int

IMPLICIT NONE

! Local variables.
INTEGER :: k
!   Loop variable
INTEGER :: i_band
!   Number of band


READ(iu_spc, '(/, 23x, i5, 26x, 1pe16.9)') &
    Sp%Planck%n_deg_fit, Sp%Planck%t_ref_planck
nd_thermal_coeff = Sp%Planck%n_deg_fit + 1
ALLOCATE(Sp%Planck%thermal_coeff(0:nd_thermal_coeff-1, nd_band))

READ(iu_spc, '(/)')
DO i=1, Sp%Basic%n_band
  READ(iu_spc, '(i5, 7x, (t13, 3(1pe16.9, 4x)))', IOSTAT=ios,            &
       IOMSG=iomessage) &
    i_band, (Sp%Planck%thermal_coeff(k, i), &
    k=0, Sp%Planck%n_deg_fit)
  IF (ios /= 0) THEN
    cmessage = '*** Error in subroutine read_block_6_0_0. ' // &
      'The data for the thermal source function could not be read: ' // &
      TRIM(iomessage)
    ierr=i_err_fatal
    RETURN
  END IF
END DO

END SUBROUTINE read_block_6_0_0_int



SUBROUTINE read_block_6_0_1_int

IMPLICIT NONE

! Local variables.
INTEGER :: k, it
!   Loop variables
INTEGER :: i_band
!   Number of band
CHARACTER(LEN=10) :: planck_data_type
!   String with type of Plack data
LOGICAL :: l_planck_table_exists


READ(iu_spc, '(/,14x, a10)') planck_data_type
IF (planck_data_type(1:5) == 'table') THEN
  Sp%Planck%l_planck_tbl=.TRUE.
ELSE
  Sp%Planck%l_planck_tbl=.FALSE.
END IF

READ(iu_spc, '(24x, i5, 26x, 1pe16.9)') &
  Sp%Planck%n_deg_fit, Sp%Planck%t_ref_planck
IF (Sp%Planck%l_planck_tbl) THEN
  Sp%Planck%n_deg_fit=Sp%Planck%n_deg_fit-1
END IF
nd_thermal_coeff = Sp%Planck%n_deg_fit + 1
ALLOCATE(Sp%Planck%thermal_coeff( &
  0:nd_thermal_coeff-1, nd_band))


IF (Sp%Planck%l_planck_tbl) THEN

  ALLOCATE(Sp%Planck%theta_planck_tbl(0:nd_thermal_coeff-1))

  ! Locate correct block in extended spectral file
  l_planck_table_exists=.FALSE.
  REWIND(iu_spc1)
  DO
    READ(iu_spc1, '(a80)') char_dum
    IF ( char_dum(1:6) == '*BLOCK' ) THEN
      IF ( char_dum(9:20) == 'Planck table' ) THEN
        l_planck_table_exists=.TRUE.
        EXIT
      END IF
    END IF
  END DO

  ! Return with error if the Planck table does not exist
  IF (.NOT. l_planck_table_exists) THEN
    cmessage = 'No table of Planck function in extended spectral file.'
    ierr = i_err_fatal
    RETURN
  END IF

  ! Read temperatures
  READ(iu_spc1,'(/,5(1PE16.8E3))', IOSTAT=ios, IOMSG=iomessage) &
    (Sp%Planck%theta_planck_tbl(it), &
    it=0, Sp%Planck%n_deg_fit)
  IF (ios /= 0) THEN
    cmessage = '*** Error in subroutine read_block_6_0_1.\n' // &
      'Planck table temperatures could not be read: ' // TRIM(iomessage)
    ierr=i_err_fatal
    RETURN
  END IF

  ! Read Planck table
  DO i=1, Sp%Basic%n_band
    READ(iu_spc1, '(/)')
    READ(iu_spc1, '(5(1PE16.8E3))', IOSTAT=ios, IOMSG=iomessage) &
      (Sp%Planck%thermal_coeff(it,i), &
       it=0, Sp%Planck%n_deg_fit)
    IF (ios /= 0) THEN
      cmessage = '*** Error in subroutine read_block_6_0_1.\n' // &
        'Planck table could not be read: ' // TRIM(iomessage)
      ierr=i_err_fatal
      RETURN
    END IF
  END DO
ELSE
  READ(iu_spc, '(/)')
  DO i=1, Sp%Basic%n_band
    READ(iu_spc, '(i5, 7x, (t13, 3(1pe16.9, 4x)))', IOSTAT=ios,          &
         IOMSG=iomessage) &
      i_band, (Sp%Planck%thermal_coeff(k, i), &
      k=0, Sp%Planck%n_deg_fit)
    IF (ios /= 0) THEN
      cmessage = '*** Error in subroutine read_block_6_0_1. ' // &
        'The data for the thermal source function could not be read: ' &
        // TRIM(iomessage)
      ierr=i_err_fatal
      RETURN
    END IF
  END DO
END IF

END SUBROUTINE read_block_6_0_1_int



SUBROUTINE read_block_8_0_0_int

IMPLICIT NONE

! Local variables.
INTEGER :: idum
!   Dummy integer
INTEGER :: j
!   Loop variable


! Allocate the continuum arrays:
ALLOCATE(Sp%Cont%n_band_continuum(nd_band))
ALLOCATE(Sp%Cont%index_continuum(nd_band, nd_continuum))
! Skip over the headers.
READ(iu_spc, '(////)')

! Read in the limits on the intervals in the spectrum
DO i=1, Sp%Basic%n_band
  READ(iu_spc, FMT='(i5, 7x, i5)', IOSTAT=ios, IOMSG=iomessage) &
    idum, Sp%Cont%n_band_continuum(i)
  IF (ios /= 0) THEN
    cmessage = '*** Error in subroutine read_block_8_0_0:\n' // &
      'the list of continua is not correct: '// TRIM(iomessage)
    ierr=i_err_fatal
    RETURN
  END IF
  IF (Sp%Cont%n_band_continuum(i) > nd_continuum) THEN
    cmessage = '*** Error in subroutine read_block_8_0_0.\n' // &
      'There are too many continua: ' // &
      'increase npd_continuum and recompile: '// TRIM(iomessage)
    ierr=i_err_fatal
    RETURN
  END IF
  IF (Sp%Cont%n_band_continuum(i) > 0) THEN
    READ(iu_spc, '(5x, 4(2x, i3))') &
      ( Sp%Cont%index_continuum(i, j), &
        j=1, Sp%Cont%n_band_continuum(i) )
  END IF
  IF (ios /= 0) THEN
    cmessage = '*** Error in subroutine read_block_8_0_0. ' // &
      'The list of continua is not correct.'
    ierr=i_err_fatal
    RETURN
  END IF
END DO

! Read in the indices of gases forming the continuum species.
READ(iu_spc, '(/, 22x, i5)') Sp%Cont%index_water

END SUBROUTINE read_block_8_0_0_int



SUBROUTINE read_block_9_0_0_int

IMPLICIT NONE

! Local variables.
INTEGER :: idum_band
!   Dummy integer
INTEGER :: idum_continuum
!   Dummy integer
INTEGER :: idum_fnc
!   Dummy integer
INTEGER :: j
!   Loop variable
INTEGER :: l
!   Loop variable


! Allocate space for the variables.
ALLOCATE(Sp%Cont%k_cont(nd_band, nd_continuum))
ALLOCATE(Sp%Cont%i_scale_fnc_cont(nd_band, nd_continuum))
ALLOCATE(Sp%Cont%scale_cont(nd_scale_variable, &
  nd_band, nd_continuum))
ALLOCATE(Sp%Cont%t_ref_cont(nd_continuum, nd_band))
ALLOCATE(Sp%Cont%p_ref_cont(nd_continuum, nd_band))

! Skip over the headers.
READ(iu_spc, '(//)')

DO i=1, Sp%Basic%n_band
  DO j=1, Sp%Cont%n_band_continuum(i)
    READ(iu_spc, FMT='(i5, 2(7x, i5))', IOSTAT=ios, IOMSG=iomessage) &
      idum_band, idum_continuum, idum_fnc
    IF (ios /= 0) THEN
      cmessage = '*** Error in subroutine read_block_9_0_0:\n' // &
        'continua in band could not be read: '// TRIM(iomessage)
      ierr=i_err_fatal
      RETURN
    END IF
    Sp%Cont%i_scale_fnc_cont(idum_band, idum_continuum)=idum_fnc
    ! Read the reference temperature and pressure.
    READ(iu_spc, '(2(6x, 1pe16.9))') &
      Sp%Cont%p_ref_cont(idum_continuum, idum_band), &
      Sp%Cont%t_ref_cont(idum_continuum, idum_band)
    ! For each band read the values of the coefficients.
    READ(iu_spc, '(6x, 1pe16.9, (t23, 2(6x, 1pe16.9)))', &
        IOSTAT=ios, IOMSG=iomessage) &
      Sp%Cont%k_cont(idum_band, idum_continuum), &
      (Sp%Cont%scale_cont(l, idum_band, idum_continuum), &
       l=1, n_scale_variable(idum_fnc))
    IF (ios /= 0) THEN
      cmessage = '*** Error in subroutine read_block_9_0_0: ' // &
        'continuum data could not be read: ' // TRIM(iomessage)
      ierr=i_err_fatal
      RETURN
    END IF
  END DO
END DO

END SUBROUTINE read_block_9_0_0_int



SUBROUTINE read_block_9_1_0_int

IMPLICIT NONE

! Local variables.
INTEGER :: idum_band
!   Dummy integer
INTEGER :: idum_continuum
!   Dummy integer
INTEGER :: n_term, num_t
!   Dummy integer
INTEGER :: j, i_term, ip, it
!   Loop variables


! Allocate space for the variables.
ALLOCATE(Sp%Cont%k_h2oc(nd_pre, nd_tmp, nd_k_term, nd_band))
ALLOCATE(Sp%Cont%k_cont_ses(nd_k_term, nd_tmp, nd_band, nd_continuum))
ALLOCATE(Sp%Cont%i_scale_fnc_cont(nd_band, nd_continuum))

Sp%Cont%i_scale_fnc_cont = ip_scale_ses2

! Skip over the headers.
READ(iu_spc, '(/)')

DO i=1, Sp%Basic%n_band
  DO j=1, Sp%Cont%n_band_continuum(i)
    READ(iu_spc, *, IOSTAT=ios, IOMSG=iomessage) &
      idum_band, idum_continuum, n_term, num_t
    IF (ios /= 0) THEN
      cmessage = '*** Error in subroutine read_block_9_1_0:\n' // &
        'continua in band could not be read: ' // TRIM(iomessage)
      ierr=i_err_fatal
      RETURN
    END IF

    IF (idum_continuum == 1) THEN
      READ(iu_spc1, '(/)')
      DO i_term=1, n_term
        DO ip=1, 21
          READ(iu_spc1, '(5(1PE12.6,1x))', IOSTAT=ios, IOMSG=iomessage) &
            (Sp%Cont%k_h2oc(ip,it,i_term,i), it=1, nd_tmp)
          IF (ios /= 0) THEN
            WRITE(cmessage, '(a, 3i4,A)') &
              '*** 2nd Error in subroutine read_block_9_1_0:\n', &
              i,i_term,ip, TRIM(iomessage)
            ierr=i_err_fatal
            RETURN
          END IF
        END DO
      END DO

    ELSE
      DO i_term=1, n_term
        READ(iu_spc, '(6X, 5(1x, 1PE16.9))', IOSTAT=ios, IOMSG=iomessage) &
          (Sp%Cont%k_cont_ses(i_term, it, idum_band, j), it=1, num_t)
      END DO
    END IF
  END DO
END DO

END SUBROUTINE read_block_9_1_0_int



SUBROUTINE read_block_10_0_0_int

IMPLICIT NONE

! Local variables.
INTEGER :: i_drop
!   Type of droplet
INTEGER :: i_parametrization_drop
!   Dummy index of parameter scheme
INTEGER :: n_parameter
!   Number of parameters
INTEGER :: k
!   Loop variable
INTEGER :: i_dummy
!   Dummy reading variable



! Allocate storage for droplet data if entering this block for the
! first time.
IF ( .NOT. ALLOCATED(Sp%Drop%l_drop_type) ) THEN
  ALLOCATE(Sp%Drop%l_drop_type(nd_drop_type))
  Sp%Drop%l_drop_type(1:nd_drop_type)=.FALSE.
  ALLOCATE(Sp%Drop%i_drop_parm(nd_drop_type))
  ALLOCATE(Sp%Drop%parm_min_dim(nd_drop_type))
  ALLOCATE(Sp%Drop%parm_max_dim(nd_drop_type))
  ALLOCATE(Sp%Drop%n_phf(nd_drop_type))
  ALLOCATE(Sp%Drop%parm_list(nd_cloud_parameter, nd_band, &
    nd_drop_type))
END IF

! Read the headers.
READ(iu_spc, '(/, 27x, i5, /, 34x, i5, 27x, i5)') &
  i_drop, i_parametrization_drop, i_dummy
IF (i_drop > nd_drop_type) THEN
  cmessage = '*** Error in subroutine read_block_10_0_0. ' // &
    'Indexing number of droplet exceeds maximum permitted value: ' // &
    'Increase npd_drop_type and recompile.'
  ierr=i_err_fatal
  RETURN
END IF

IF ( (i_parametrization_drop == IP_slingo_schrecker)  .OR. &
     (i_parametrization_drop == IP_ackerman_stephens) .OR. &
     (i_parametrization_drop == IP_drop_pade_2) ) THEN
  ! Data are parametrized.


  ! Coding for backward compatibility: default settings of 0
  ! introduced.
  Sp%Drop%parm_min_dim(i_drop)=0.0_RealK
  Sp%Drop%parm_max_dim(i_drop)=0.0_RealK
  ! Only a parametrization of the asymmetry can be accommodated.
  Sp%Drop%n_phf(i_drop)=1


  n_parameter=i_dummy

  DO i=1, Sp%Basic%n_band
    READ(iu_spc, FMT='(/, (4(4x, 1pe12.5)))', IOSTAT=ios, IOMSG=iomessage) &
      (Sp%Drop%parm_list(k, i, i_drop), &
       k=1, n_parameter)
    ! For each band read the values of the parameters.
    IF (ios /= 0) THEN
      cmessage = '*** Error in subroutine read_block_10_0_0:\n' // &
        'data for droplets are not in the correct format: ' // &
        TRIM(iomessage)
      ierr=i_err_fatal
      RETURN
    END IF
  END DO
ELSE
  ! Illegal parametrization scheme encountered.
  cmessage = '*** Error in subroutine read_block_10_0_0:\n' // &
    'an unknown parametrization scheme has been specified.'
  ierr=i_err_fatal
  RETURN
END IF

! Record the presence of the drop type and the index
! of the parametrization
Sp%Drop%l_drop_type(i_drop)= .TRUE.
Sp%Drop%i_drop_parm(i_drop)=i_parametrization_drop

END SUBROUTINE read_block_10_0_0_int



SUBROUTINE read_block_10_0_1_int

IMPLICIT NONE

! Local variables.
INTEGER :: i_drop
!   Type of droplet
INTEGER :: i_parametrization_drop
!   Dummy index of parameter scheme
INTEGER :: n_parameter
!   Number of parameters
INTEGER :: i
!   Loop variable
INTEGER :: k
!   Loop variable
INTEGER :: ios
!   I/O error flag
INTEGER :: i_dummy
!   Dummy reading variable



! Allocate storage for droplet data if entering this block for
! the first time.
IF ( .NOT. ALLOCATED(Sp%Drop%l_drop_type) ) THEN
  ALLOCATE(Sp%Drop%l_drop_type(nd_drop_type))
  Sp%Drop%l_drop_type(1:nd_drop_type)=.FALSE.
  ALLOCATE(Sp%Drop%i_drop_parm(nd_drop_type))
  ALLOCATE(Sp%Drop%parm_min_dim(nd_drop_type))
  ALLOCATE(Sp%Drop%parm_max_dim(nd_drop_type))
  ALLOCATE(Sp%Drop%n_phf(nd_drop_type))
  ALLOCATE(Sp%Drop%parm_list(nd_cloud_parameter, nd_band, &
    nd_drop_type))
END IF

! Read the headers.
READ(iu_spc, '(/, 27x, i5, /, 34x, i5, 27x, i5)') &
  i_drop, i_parametrization_drop, i_dummy
IF (i_drop > nd_drop_type) THEN
  cmessage = '*** Error in subroutine read_block_10_0_1.\n' // &
    'The indexing number of a droplet exceeds the '// &
    'maximum permitted value: ' // &
    'increase npd_drop_type and recompile.'
  ierr=i_err_fatal
  RETURN
END IF

IF ( (i_parametrization_drop == IP_slingo_schrecker)  .OR. &
     (i_parametrization_drop == IP_ackerman_stephens) .OR. &
     (i_parametrization_drop == IP_drop_pade_2) ) THEN
  ! Data are parametrized.

  ! Settings for backward compatibility:
  ! Only a parametrization of the asymmetry can be accommodated.
  Sp%Drop%n_phf(i_drop)=1

  READ(iu_spc, '(39x, 1pe12.5, 4x, 1pe12.5)') &
    Sp%Drop%parm_min_dim(i_drop), &
    Sp%Drop%parm_max_dim(i_drop)

  n_parameter=i_dummy

  DO i=1, Sp%Basic%n_band
    READ(iu_spc, FMT='(/, (4(4x, 1pe12.5)))', IOSTAT=ios, IOMSG=iomessage) &
      (Sp%Drop%parm_list(k, i, i_drop), k=1, n_parameter)
    ! For each band read the values of the parameters.
    IF (ios /= 0) THEN
      cmessage = '*** Error in subroutine read_block_10_0_1.\n' // &
        'Data for droplets are not in the correct format: ' // &
        TRIM(iomessage)
      ierr=i_err_fatal
      RETURN
    END IF
  END DO
ELSE
  ! Illegal parametrization scheme encountered.
  cmessage = '*** Error in subroutine read_block_10_0_1.\n' // &
    'An unknown parametrization scheme has been specified.'
  ierr=i_err_fatal
  RETURN
END IF

! Record the presence of the drop type and the index
! of the parametrization
Sp%Drop%l_drop_type(i_drop)= .TRUE.
Sp%Drop%i_drop_parm(i_drop)=i_parametrization_drop

END SUBROUTINE read_block_10_0_1_int



SUBROUTINE read_block_10_0_2_int

IMPLICIT NONE

! Local variables.
INTEGER :: i_drop
!   Type of droplet
INTEGER :: i_parametrization_drop
!   Dummy index of parameter scheme
INTEGER :: n_parameter
!   Number of parameters
INTEGER :: i
!   Loop variable
INTEGER :: k
!   Loop variable
INTEGER :: ios
!   I/O error flag
INTEGER :: i_dummy
!   Dummy reading variable



! Allocate storage for droplet data if entering this block for
! the first time.
IF ( .NOT. ALLOCATED(Sp%Drop%l_drop_type) ) THEN
  ALLOCATE(Sp%Drop%l_drop_type(nd_drop_type))
  Sp%Drop%l_drop_type(1:nd_drop_type)=.FALSE.
  ALLOCATE(Sp%Drop%i_drop_parm(nd_drop_type))
  ALLOCATE(Sp%Drop%parm_min_dim(nd_drop_type))
  ALLOCATE(Sp%Drop%parm_max_dim(nd_drop_type))
  ALLOCATE(Sp%Drop%n_phf(nd_drop_type))
  ALLOCATE(Sp%Drop%parm_list(nd_cloud_parameter, nd_band, &
    nd_drop_type))
END IF

! Read the headers.
READ(iu_spc, '(/, 27x, i5, /, 34x, i5, 27x, i5)') &
  i_drop, i_parametrization_drop, i_dummy
IF (i_drop > nd_drop_type) THEN
  cmessage = '*** Error in subroutine read_block_10_0_2.\n' // &
    'The indexing number of a droplet exceeds the maximum ' // &
    'permitted value: increase npd_drop_type and recompile.'
  ierr=i_err_fatal
  RETURN
END IF

IF ( (i_parametrization_drop == IP_slingo_schrecker)  .OR. &
     (i_parametrization_drop == IP_Slingo_Schr_PHF )  .OR. &
     (i_parametrization_drop == IP_ps_size_PHF )      .OR. &
     (i_parametrization_drop == IP_ackerman_stephens) .OR. &
     (i_parametrization_drop == IP_drop_pade_2) ) THEN
  ! Data are parametrized.

  READ(iu_spc, '(42x, i5)') Sp%Drop%n_phf(i_drop)

  READ(iu_spc, '(39x, 1pe12.5, 4x, 1pe12.5)') &
    Sp%Drop%parm_min_dim(i_drop), &
    Sp%Drop%parm_max_dim(i_drop)

  n_parameter=i_dummy

  DO i=1, Sp%Basic%n_band
    READ(iu_spc, FMT='(/, (4(4x, 1pe12.5)))', IOSTAT=ios, IOMSG=iomessage) &
      (Sp%Drop%parm_list(k, i, i_drop), k=1, n_parameter)
    ! For each band read the values of the parameters.
    IF (ios /= 0) THEN
      cmessage = '*** Error in subroutine read_block_10_0_2.\n' // &
        'Data for droplets are not in the correct format: ' // &
        TRIM(iomessage)
      ierr=i_err_fatal
      RETURN
    END IF
  END DO
ELSE
  ! Illegal parametrization scheme encountered.
  cmessage = '*** Error in subroutine read_block_10_0_2.\n' // &
    'An unknown parametrization scheme has been specified.'
  ierr=i_err_fatal
  RETURN
END IF

! Record the presence of the drop type and the index
! of the parametrization
Sp%Drop%l_drop_type(i_drop)= .TRUE.
Sp%Drop%i_drop_parm(i_drop)=i_parametrization_drop

END SUBROUTINE read_block_10_0_2_int



SUBROUTINE read_block_11_0_1_int

IMPLICIT NONE

! Local variables.
INTEGER :: i_species
!   Index of species
INTEGER :: idum
!   Dummy reading variable


! Allocate space for aerosol data if entering this block for the
! first time.
IF ( .NOT. ALLOCATED(Sp%Aerosol%l_aero_spec) ) THEN
  ALLOCATE(Sp%Aerosol%l_aero_spec(nd_aerosol_species))
  Sp%Aerosol%l_aero_spec(1:nd_aerosol_species) = .FALSE.
  ALLOCATE(Sp%Aerosol%i_aerosol_parm(nd_aerosol_species))
  ALLOCATE(Sp%Aerosol%nhumidity(nd_aerosol_species))
  ALLOCATE(Sp%Aerosol%humidities(nd_humidity, nd_aerosol_species))
  ALLOCATE(Sp%Aerosol%n_aerosol_phf_term(nd_aerosol_species))
  ALLOCATE(Sp%Aerosol%abs(nd_humidity, nd_aerosol_species, nd_band))
  ALLOCATE(Sp%Aerosol%scat(nd_humidity, nd_aerosol_species, nd_band))
  ALLOCATE(Sp%Aerosol%phf_fnc(nd_humidity, nd_phase_term, &
    nd_aerosol_species, nd_band))
END IF
READ(iu_spc, '(/, 19x, i5, //)') i_species
IF (i_species > nd_aerosol_species) THEN
  WRITE(cmessage, '(a,i2,a)') &
    '*** Error in subroutine read_block_11_0_1. Aerosol index ', &
    i_species, ' is greater than the total number listed in BLOCK 0'
  ierr=i_err_fatal
  RETURN
END IF

! Read in the scattering parameters for each band
DO i=1, Sp%Basic%n_band
  READ(iu_spc, *, IOSTAT=ios, IOMSG=iomessage) idum, &
   Sp%Aerosol%abs(1, i_species, i), &
   Sp%Aerosol%scat(1, i_species, i), &
   Sp%Aerosol%phf_fnc(1, 1, i_species, i)
  IF (ios /= 0) THEN
    cmessage = '*** Error in subroutine read_block_11_0_1.\n' // &
      'Dry aerosol scattering data are not in the correct format: ' //&
      TRIM(iomessage)
    ierr=i_err_fatal
    RETURN
  END IF
END DO
Sp%Aerosol%nhumidity(i_species)=0
Sp%Aerosol%i_aerosol_parm(i_species)=IP_aerosol_param_dry
Sp%Aerosol%n_aerosol_phf_term(i_species)=1

! After successful reading the presence of this species is recorded.
Sp%Aerosol%l_aero_spec(i_species)= .TRUE.

END SUBROUTINE read_block_11_0_1_int



SUBROUTINE read_block_11_0_2_int

IMPLICIT NONE

! Local variables.
INTEGER :: l
!   Loop variable
INTEGER :: i_species
!   Index of species
INTEGER :: idum
!   Dummy reading variable


! Allocate space for aerosol data if entering
! this block for the first time.
IF ( .NOT. ALLOCATED(Sp%Aerosol%l_aero_spec) ) THEN
  ALLOCATE(Sp%Aerosol%l_aero_spec(nd_aerosol_species))
  Sp%Aerosol%l_aero_spec(1:nd_aerosol_species) = .FALSE.
  ALLOCATE(Sp%Aerosol%i_aerosol_parm(nd_aerosol_species))
  ALLOCATE(Sp%Aerosol%nhumidity(nd_aerosol_species))
  ALLOCATE(Sp%Aerosol%n_aerosol_phf_term(nd_aerosol_species))
  ALLOCATE(Sp%Aerosol%humidities(nd_humidity, nd_aerosol_species))
END IF

READ(iu_spc, '(/, 19x, i5)') i_species
IF (i_species > nd_aerosol_species) THEN
  WRITE(cmessage, '(a,i2,a)') &
    '*** Error in subroutine read_block_11_0_2. Aerosol index ', &
    i_species, ' is greater than the total number listed in BLOCK 0'
  ierr=i_err_fatal
  RETURN
END IF
READ(iu_spc, '(36x, i5, //)') &
  Sp%Aerosol%n_aerosol_phf_term(i_species)

IF ( .NOT. ALLOCATED(Sp%Aerosol%abs) ) THEN
  ALLOCATE(Sp%Aerosol%abs(nd_humidity, &
    nd_aerosol_species, nd_band))
  ALLOCATE(Sp%Aerosol%scat(nd_humidity, &
    nd_aerosol_species, nd_band))
  ALLOCATE(Sp%Aerosol%phf_fnc(nd_humidity, nd_phase_term, &
    nd_aerosol_species, nd_band))
END IF

! Read in the scattering parameters for each band
DO i=1, Sp%Basic%n_band
  READ(iu_spc, *, IOSTAT=ios, IOMSG=iomessage) idum, &
    Sp%Aerosol%abs(1, i_species, i), &
    Sp%Aerosol%scat(1, i_species, i), &
    (Sp%Aerosol%phf_fnc(1, l, i_species, i), &
     l=1, Sp%Aerosol%n_aerosol_phf_term(i_species))
  IF (ios /= 0) THEN
    cmessage = '*** Error in subroutine read_block_11_0_2.\n' // &
      'Dry aerosol scattering data are not in the correct format: ' //&
      TRIM(iomessage)
    ierr=i_err_fatal
    RETURN
  END IF
END DO
Sp%Aerosol%nhumidity(i_species)=0
IF (Sp%Aerosol%n_aerosol_phf_term(i_species) == 1) THEN
  Sp%Aerosol%i_aerosol_parm(i_species)=IP_aerosol_param_dry
ELSE IF (Sp%Aerosol%n_aerosol_phf_term(i_species) > 1) THEN
  Sp%Aerosol%i_aerosol_parm(i_species)=IP_aerosol_param_phf_dry
END IF

! After successful reading the presence of this species is recorded.
Sp%Aerosol%l_aero_spec(i_species)= .TRUE.

END SUBROUTINE read_block_11_0_2_int



SUBROUTINE read_block_11_1_0_int

IMPLICIT NONE

! Local variables.
INTEGER :: k
!   Loop variable
INTEGER :: i_component
!   Index of component


! Allocate arrays for aerosols.
IF ( .NOT. ALLOCATED(Sp%Aerosol%l_aero_spec) ) THEN
  ALLOCATE(Sp%Aerosol%l_aero_spec(nd_aerosol_species))
  Sp%Aerosol%l_aero_spec(1:nd_aerosol_species) = .FALSE.
  ALLOCATE(Sp%Aerosol%i_aerosol_parm(nd_aerosol_species))
  ALLOCATE(Sp%Aerosol%n_aerosol_phf_term(nd_aerosol_species))
  ALLOCATE(Sp%Aerosol%nhumidity(nd_aerosol_species))
  ALLOCATE(Sp%Aerosol%humidities(nd_humidity,nd_aerosol_species))
END IF

READ(iu_spc, '(/, 19x, i5 )') i_component
IF (i_component > nd_aerosol_species) THEN
  WRITE(cmessage, '(a,i2,a)') &
    '*** Error in subroutine read_block_11_1_0. Aerosol index ', &
    i_component, ' is greater than the total number listed in BLOCK 0'
  ierr=i_err_fatal
  RETURN
END IF
READ(iu_spc, '(28x, i3)') Sp%Aerosol%nhumidity(i_component)

IF ( .NOT.  ALLOCATED(Sp%Aerosol%abs) ) THEN
  ALLOCATE(Sp%Aerosol%abs(nd_humidity, nd_aerosol_species, &
    nd_band))
  ALLOCATE(Sp%Aerosol%scat(nd_humidity, nd_aerosol_species, &
    nd_band))
  ALLOCATE(Sp%Aerosol%phf_fnc(nd_humidity, nd_phase_term, &
    nd_aerosol_species, nd_band))
END IF

! Read in the scattering parameters for each band
DO i=1, Sp%Basic%n_band
  READ(iu_spc, '(//)')
  DO k=1, Sp%Aerosol%nhumidity(i_component)
    READ(iu_spc, *, IOSTAT=ios, IOMSG=iomessage) &
      Sp%Aerosol%humidities(k, i_component), &
      Sp%Aerosol%abs(k, i_component, i), &
      Sp%Aerosol%scat(k, i_component, i), &
      Sp%Aerosol%phf_fnc(k, 1, i_component, i)
    IF (ios /= 0) THEN
      cmessage = '*** Error in subroutine read_block_11_1_0.\n' // &
        'Moist aerosol scattering data are not in the correct format: ' &
        // TRIM(iomessage)
      ierr=i_err_fatal
      RETURN
    END IF
  END DO
END DO

Sp%Aerosol%i_aerosol_parm(i_component)=IP_aerosol_param_moist
Sp%Aerosol%n_aerosol_phf_term(i_component)=1

! After successful reading the presence of this component is recorded.
Sp%Aerosol%l_aero_spec(i_component)= .TRUE.

END SUBROUTINE read_block_11_1_0_int



SUBROUTINE read_block_11_1_1_int

IMPLICIT NONE

! Local variables.
INTEGER :: k
!   Loop variable
INTEGER :: i_species
!   Index of component


IF ( .NOT. ALLOCATED(Sp%Aerosol%l_aero_spec) ) THEN
  ALLOCATE(Sp%Aerosol%l_aero_spec(nd_aerosol_species))
  Sp%Aerosol%l_aero_spec(1:nd_aerosol_species) = .FALSE.
  ALLOCATE(Sp%Aerosol%i_aerosol_parm(nd_aerosol_species))
  ALLOCATE(Sp%Aerosol%n_aerosol_phf_term(nd_aerosol_species))
  ALLOCATE(Sp%Aerosol%nhumidity(nd_aerosol_species))
  ALLOCATE(Sp%Aerosol%humidities(nd_humidity,nd_aerosol_species))
END IF

READ(iu_spc, '(/, 19x, i5 )') i_species
IF (i_species > nd_aerosol_species) THEN
  WRITE(cmessage, '(a,i2,a)') &
    '*** Error in subroutine read_block_11_1_1. Aerosol index ', &
    i_species, ' is greater than the total number listed in BLOCK 0'
  ierr=i_err_fatal
  RETURN
END IF
READ(iu_spc, '(28x, i3)') Sp%Aerosol%nhumidity(i_species)

IF ( .NOT. ALLOCATED(Sp%Aerosol%abs) ) THEN
  ALLOCATE(Sp%Aerosol%abs(nd_humidity, nd_aerosol_species, &
    nd_band))
  ALLOCATE(Sp%Aerosol%scat(nd_humidity, nd_aerosol_species, &
    nd_band))
  ALLOCATE(Sp%Aerosol%phf_fnc(nd_humidity, nd_phase_term, &
    nd_aerosol_species, nd_band))
END IF

! Read in the scattering parameters for each band: only the
! asymmetry is used here.
DO i=1, Sp%Basic%n_band
  READ(iu_spc, '(//)')
  DO k=1, Sp%Aerosol%nhumidity(i_species)
    READ(iu_spc, *, IOSTAT=ios, IOMSG=iomessage) &
      Sp%Aerosol%humidities(k, i_species), &
      Sp%Aerosol%abs(k, i_species, i), &
      Sp%Aerosol%scat(k, i_species, i), &
      Sp%Aerosol%phf_fnc(k, 1, i_species, i)
    IF (ios /= 0) THEN
      cmessage = '*** Error in subroutine read_block_11_1_1.\n' // &
        'Moist aerosol scattering data are not in the correct format: ' &
        // TRIM(iomessage)
      ierr=i_err_fatal
      RETURN
    END IF
  END DO
END DO

Sp%Aerosol%i_aerosol_parm(i_species)=IP_aerosol_param_moist
Sp%Aerosol%n_aerosol_phf_term(i_species)=1

! After sucessful reading the presence of this species is specified.
Sp%Aerosol%l_aero_spec(i_species)= .TRUE.

END SUBROUTINE read_block_11_1_1_int



SUBROUTINE read_block_11_1_2_int

IMPLICIT NONE

! Local variables.
INTEGER :: k
!   Loop variable
INTEGER :: l
!   Loop variable
INTEGER :: i_species
!   Index of component


! Allocate space for aerosol arrays.
IF ( .NOT. ALLOCATED(Sp%Aerosol%l_aero_spec) ) THEN
  ALLOCATE(Sp%Aerosol%l_aero_spec(nd_aerosol_species))
  Sp%Aerosol%l_aero_spec(1:nd_aerosol_species) = .FALSE.
  ALLOCATE(Sp%Aerosol%i_aerosol_parm(nd_aerosol_species))
  ALLOCATE(Sp%Aerosol%n_aerosol_phf_term(nd_aerosol_species))
  ALLOCATE(Sp%Aerosol%nhumidity(nd_aerosol_species))
  ALLOCATE(Sp%Aerosol%humidities(nd_humidity,nd_aerosol_species))
END IF

READ(iu_spc, '(/, 19x, i5 )') i_species
IF (i_species > nd_aerosol_species) THEN
  WRITE(cmessage, '(a,i2,a)') &
    '*** Error in subroutine read_block_11_1_2. Aerosol index ', &
    i_species, ' is greater than the total number listed in BLOCK 0'
  ierr=i_err_fatal
  RETURN
END IF
READ(iu_spc, '(28x, i3)') Sp%Aerosol%nhumidity(i_species)
READ(iu_spc, '(36x, i5)') Sp%Aerosol%n_aerosol_phf_term(i_species)

IF ( .NOT. ALLOCATED(Sp%Aerosol%abs) ) THEN
  ALLOCATE(Sp%Aerosol%abs(nd_humidity, nd_aerosol_species, &
    nd_band))
  ALLOCATE(Sp%Aerosol%scat(nd_humidity, nd_aerosol_species, &
    nd_band))
  ALLOCATE(Sp%Aerosol%phf_fnc(nd_humidity, nd_phase_term, &
    nd_aerosol_species, nd_band))
END IF

! Read in the scattering parameters for each band.
DO i=1, Sp%Basic%n_band
  READ(iu_spc, '(//)')
  DO k=1, Sp%Aerosol%nhumidity(i_species)
    READ(iu_spc, *, IOSTAT=ios, IOMSG=iomessage) &
      Sp%Aerosol%humidities(k, i_species), &
      Sp%Aerosol%abs(k, i_species, i), &
      Sp%Aerosol%scat(k, i_species, i), &
      (Sp%Aerosol%phf_fnc(k, l, i_species, i), &
       l=1, Sp%Aerosol%n_aerosol_phf_term(i_species))
    IF (ios /= 0) THEN
      cmessage = '*** Error in subroutine read_block_11_1_2.\n' // &
        'Moist aerosol scattering data are not in the correct format:' //&
        TRIM(iomessage)
      ierr=i_err_fatal
      RETURN
    END IF
  END DO
END DO

IF (Sp%Aerosol%n_aerosol_phf_term(i_species) == 1) THEN
  Sp%Aerosol%i_aerosol_parm(i_species)=IP_aerosol_param_moist
ELSE IF (Sp%Aerosol%n_aerosol_phf_term(i_species) > 1) THEN
  Sp%Aerosol%i_aerosol_parm(i_species)=IP_aerosol_param_phf_moist
END IF

! After sucessful reading the presence of this species is recorded.
Sp%Aerosol%l_aero_spec(i_species)= .TRUE.

END SUBROUTINE read_block_11_1_2_int



SUBROUTINE read_block_12_0_0_int

IMPLICIT NONE

! Local variables.
INTEGER :: i_ice
!   Type of ice crystal
INTEGER :: i_parametrization_ice
!   Dummy index of parameter scheme
INTEGER :: n_parameter
!   Number of parameters
INTEGER :: k
!   Loop variable
INTEGER :: i_dummy
!   Dummy reading variable


! Allocate storage for ice data.
IF ( .NOT. ALLOCATED(Sp%Ice%l_ice_type) ) THEN
  ALLOCATE(Sp%Ice%l_ice_type(nd_ice_type))
  Sp%Ice%l_ice_type(1:nd_ice_type)=.FALSE.
  ALLOCATE(Sp%Ice%i_ice_parm(nd_ice_type))
  ALLOCATE(Sp%Ice%parm_min_dim(nd_ice_type))
  ALLOCATE(Sp%Ice%parm_max_dim(nd_ice_type))
  ALLOCATE(Sp%Ice%n_phf(nd_ice_type))
  ALLOCATE(Sp%Ice%parm_list(nd_cloud_parameter, nd_band, &
    nd_ice_type))
END IF

! Read the headers.
READ(iu_spc, '(/, 31x, i5, /, 34x, i5, 27x, i5)') &
  i_ice, i_parametrization_ice, i_dummy
IF (i_ice > nd_ice_type) THEN
  cmessage = '*** Error in subroutine read_block_12_0_0.\n' // &
    'Type of ice crystal exceeds maximum permitted value: ' // &
    'Increase npd_ice_type and recompile.'
  ierr=i_err_fatal
  RETURN
END IF

IF ( (i_parametrization_ice == IP_slingo_schrecker_ice) .OR. &
     (i_parametrization_ice == IP_ice_adt)              .OR. &
     (i_parametrization_ice == IP_ice_adt_10)           .OR. &
     (i_parametrization_ice == IP_ice_fu_solar)         .OR. &
     (i_parametrization_ice == IP_ice_fu_ir)            .OR. &
     (i_parametrization_ice == IP_sun_shine_vn2_vis)    .OR. &
     (i_parametrization_ice == IP_sun_shine_vn2_ir) ) THEN
  ! Data are parametrized.
  n_parameter=i_dummy


  ! Code for backward compatibility: set the range of validity
  ! of the parametrization to 0 to flag an unset range.
  Sp%Ice%parm_min_dim(i_ice)=0.0_RealK
  Sp%Ice%parm_max_dim(i_ice)=0.0_RealK
  ! Only one moment of the asymmetry can be accommodated here.
  Sp%Ice%n_phf(i_ice)=1


  DO i=1, Sp%Basic%n_band
    READ(iu_spc, '()')
    READ(iu_spc, *, IOSTAT=ios, IOMSG=iomessage) &
      (Sp%Ice%parm_list(k, i, i_ice), k=1, n_parameter)
    ! For each band read the values of the parameters.
    IF (ios /= 0) THEN
      cmessage = '*** Error in subroutine read_block_12_0_0:\n' // &
        'Data for ice crystals are not in the correct format: ' // &
        TRIM(iomessage)
      ierr=i_err_fatal
      RETURN
    END IF
  END DO
ELSE
  ! Illegal parametrization scheme encountered.
  cmessage = '*** Error in subroutine read_block_12_0_0:\n' // &
    'An unknown parametrization scheme has been specified.'
  ierr=i_err_fatal
  RETURN
END IF

! Record the presence of the ice crystal type and the index
! of the parametrization
Sp%Ice%l_ice_type(i_ice)= .TRUE.
Sp%Ice%i_ice_parm(i_ice)=i_parametrization_ice

END SUBROUTINE read_block_12_0_0_int



SUBROUTINE read_block_12_0_1_int

IMPLICIT NONE

! Local variables.
INTEGER :: i_ice
!   Type of ice crystal
INTEGER :: i_parametrization_ice
!   Dummy index of parameter scheme
INTEGER :: n_parameter
!   Number of parameters
INTEGER :: k
!   Loop variable
INTEGER :: i_dummy
!   Dummy reading variable


! Allocate storage for ice data.
IF ( .NOT. ALLOCATED(Sp%Ice%l_ice_type) ) THEN
  ALLOCATE(Sp%Ice%l_ice_type(nd_ice_type))
  Sp%Ice%l_ice_type(1:nd_ice_type)=.FALSE.
  ALLOCATE(Sp%Ice%i_ice_parm(nd_ice_type))
  ALLOCATE(Sp%Ice%parm_min_dim(nd_ice_type))
  ALLOCATE(Sp%Ice%parm_max_dim(nd_ice_type))
  ALLOCATE(Sp%Ice%n_phf(nd_ice_type))
  ALLOCATE(Sp%Ice%parm_list(nd_cloud_parameter, nd_band, &
    nd_ice_type))
END IF

! Read the headers.
READ(iu_spc, '(/, 31x, i5, /, 34x, i5, 27x, i5)') &
  i_ice, i_parametrization_ice, i_dummy
IF (i_ice > nd_ice_type) THEN
  cmessage = '*** Error in subroutine read_block_12_0_1.\n' // &
    'Type of ice crystal exceeds maximum permitted value: ' // &
    'Increase npd_ice_type and recompile.'
  ierr=i_err_fatal
  RETURN
END IF

IF ( (i_parametrization_ice == IP_slingo_schrecker_ice) .OR. &
     (i_parametrization_ice == IP_ice_adt)              .OR. &
     (i_parametrization_ice == IP_ice_adt_10)           .OR. &
     (i_parametrization_ice == IP_ice_fu_solar)         .OR. &
     (i_parametrization_ice == IP_ice_fu_ir)            .OR. &
     (i_parametrization_ice == IP_sun_shine_vn2_vis)    .OR. &
     (i_parametrization_ice == IP_sun_shine_vn2_ir) ) THEN
  ! Data are parametrized.
  n_parameter=i_dummy

  ! Options for backward compatibility:
  ! Only one moment of the asymmetry can be accommodated here.
  Sp%Ice%n_phf(i_ice)=1

  READ(iu_spc, '(39x, 1pe12.5, 4x, 1pe12.5)') &
    Sp%Ice%parm_min_dim(i_ice), &
    Sp%Ice%parm_max_dim(i_ice)

  DO i=1, Sp%Basic%n_band
    READ(iu_spc, '()')
    READ(iu_spc, *, IOSTAT=ios, IOMSG=iomessage) &
      Sp%Ice%parm_list(1:n_parameter, i, i_ice)
    ! For each band read the values of the parameters.
    IF (ios /= 0) THEN
      cmessage = '*** Error in subroutine read_block_12_0_1:\n' // &
        'Data for ice crystals are not in the correct format: ' // &
        TRIM(iomessage)
      ierr=i_err_fatal
      RETURN
    END IF
  END DO
ELSE
  ! Illegal parametrization scheme encountered.
  cmessage = '*** Error in subroutine read_block_12_0_1:\n' // &
    'An unknown parametrization scheme has been specified.'
  ierr=i_err_fatal
  RETURN
END IF

! Record the presence of the ice crystal type and the index
! of the parametrization
Sp%Ice%l_ice_type(i_ice)= .TRUE.
Sp%Ice%i_ice_parm(i_ice)=i_parametrization_ice

END SUBROUTINE read_block_12_0_1_int



SUBROUTINE read_block_12_0_2_int

IMPLICIT NONE

! Local variables.
INTEGER :: i_ice
!   Type of ice crystal
INTEGER :: i_parametrization_ice
!   Dummy index of parameter scheme
INTEGER :: n_parameter
!   Number of parameters
INTEGER :: k
!   Loop variable
INTEGER :: i_dummy
!   Dummy reading variable


! Allocate storage for ice data.
IF ( .NOT. ALLOCATED(Sp%Ice%l_ice_type) ) THEN
  ALLOCATE(Sp%Ice%l_ice_type(nd_ice_type))
  Sp%Ice%l_ice_type(1:nd_ice_type)=.FALSE.
  ALLOCATE(Sp%Ice%i_ice_parm(nd_ice_type))
  ALLOCATE(Sp%Ice%parm_min_dim(nd_ice_type))
  ALLOCATE(Sp%Ice%parm_max_dim(nd_ice_type))
  ALLOCATE(Sp%Ice%n_phf(nd_ice_type))
  ALLOCATE(Sp%Ice%parm_list(nd_cloud_parameter, nd_band, &
    nd_ice_type))
END IF

! Read the headers.
READ(iu_spc, '(/, 31x, i5, /, 34x, i5, 27x, i5)') &
  i_ice, i_parametrization_ice, i_dummy
IF (i_ice > nd_ice_type) THEN
  cmessage = '*** Error in subroutine read_block_12_0_2.\n' // &
    'Type of ice crystal exceeds maximum permitted value: ' // &
    'Increase npd_ice_type and recompile.'
  ierr=i_err_fatal
  RETURN
END IF

IF ( (i_parametrization_ice == IP_slingo_schrecker_ice) .OR. &
     (i_parametrization_ice == IP_ice_adt)              .OR. &
     (i_parametrization_ice == IP_ice_adt_10)           .OR. &
     (i_parametrization_ice == IP_ice_fu_solar)         .OR. &
     (i_parametrization_ice == IP_ice_fu_ir)            .OR. &
     (i_parametrization_ice == IP_slingo_schr_ice_phf)  .OR. &
     (i_parametrization_ice == IP_ice_fu_phf)           .OR. &
     (i_parametrization_ice == IP_sun_shine_vn2_vis)    .OR. &
     (i_parametrization_ice == IP_sun_shine_vn2_ir)     .OR. &
     (i_parametrization_ice == ip_ice_t_iwc)            .OR. &
     (i_parametrization_ice == ip_ice_iwc_only)         .OR. &
     (i_parametrization_ice == ip_ice_baran) ) THEN
  ! Data are parametrized.
  n_parameter=i_dummy

  READ(iu_spc, '(42x, i5)') Sp%Ice%n_phf(i_ice)

  READ(iu_spc, '(39x, 1pe12.5, 4x, 1pe12.5)') &
    Sp%Ice%parm_min_dim(i_ice), &
    Sp%Ice%parm_max_dim(i_ice)

  DO i=1, Sp%Basic%n_band
    READ(iu_spc, '()')
    READ(iu_spc, *, IOSTAT=ios, IOMSG=iomessage) &
      Sp%Ice%parm_list(1:n_parameter, i, i_ice)
    ! For each band read the values of the parameters.
    IF (ios /= 0) THEN
      cmessage = '*** Error in subroutine read_block_12_0_2:\n' // &
        'Data for ice crystals are not in the correct format: ' // &
        TRIM(iomessage)
      ierr=i_err_fatal
      RETURN
    END IF
  END DO
ELSE
  ! Illegal parametrization scheme encountered.
  cmessage = '*** Error in subroutine read_block_12_0_2:\n' // &
    'An unknown parametrization scheme has been specified.'
  ierr=i_err_fatal
  RETURN
END IF

! Record the presence of the ice crystal type and the index
! of the parametrization
Sp%Ice%l_ice_type(i_ice)= .TRUE.
Sp%Ice%i_ice_parm(i_ice)=i_parametrization_ice

END SUBROUTINE read_block_12_0_2_int



SUBROUTINE read_block_13_0_0_int

IMPLICIT NONE

! Local variables.
INTEGER :: idum
!   Dummy reading variable
LOGICAL :: l_dum
!   Dummy logical
REAL  (RealK) :: r_dum
!   Dummy real

ALLOCATE(Sp%Gas%l_doppler  (nd_species))
ALLOCATE(Sp%Gas%doppler_cor(nd_species))

READ(iu_spc, '(/)')

! Read in the Doppler information for each gas.
DO i=1, Sp%Gas%n_absorb
  READ(iu_spc, FMT='(i5, 7x, l6, 9x, 1pe12.5)', IOSTAT=ios,        &
       IOMSG=iomessage) &
    idum, l_dum, r_dum
  IF (ios /= 0) THEN
    cmessage = '*** Error: doppler broadening data are not in ' // &
      'the correct format: ' // TRIM(iomessage)
    ierr=i_err_fatal
    RETURN
  END IF
  Sp%Gas%l_doppler(idum)=l_dum
  Sp%Gas%doppler_cor(idum)=r_dum
END DO

END SUBROUTINE read_block_13_0_0_int



SUBROUTINE read_block_14_0_0_int

IMPLICIT NONE

! Local variables.
INTEGER :: idum
!   Dummy integer
INTEGER :: j
!   Loop variable


! Allocate space for arrys dealing with exclusions: n_band_exclude
! has been allocated earlier.
ALLOCATE(Sp%Basic%index_exclude(nd_exclude, nd_band))
! Skip over the headers.
READ(iu_spc, '(//)')

! Read in the list of excluded bands for each band in turn.
DO i=1, Sp%Basic%n_band
  READ(iu_spc, FMT='(i5, 7x, i5)', IOSTAT=ios, IOMSG=iomessage) &
    idum, Sp%Basic%n_band_exclude(i)
  IF (ios /= 0) THEN
    cmessage = '*** Error in subroutine read_block_14_0_0.\n' // &
      'The list of excluded bands is not correct: ' // TRIM(iomessage)
    ierr=i_err_fatal
    RETURN
  END IF
  IF (Sp%Basic%n_band_exclude(i) > 0) THEN
    READ(iu_spc, '(14x, 8(3x, i5))') &
      (Sp%Basic%index_exclude(j, i), &
       j=1, Sp%Basic%n_band_exclude(i) )
  END IF
  IF (ios /= 0) THEN
    cmessage = '*** Error in subroutine read_block_14_0_0. ' // &
      'The index of excluded bands is not correct.'
    ierr=i_err_fatal
    RETURN
  END IF
END DO

END SUBROUTINE read_block_14_0_0_int


! Dry aerosol optical depth
SUBROUTINE read_block_15_0_0

IMPLICIT NONE

! Local variables.
INTEGER :: idum_band
!   Dummy integer
INTEGER :: j
!   AOD type
INTEGER :: k
!   Loop variable

READ(iu_spc, '(/14X, I5)') Sp%Aerosol%n_aod_wavel

IF ( .NOT. ALLOCATED(Sp%Aerosol%i_aod_type) ) THEN
  IF (Sp%Aerosol%n_aod_wavel /= 6) THEN
    cmessage = '*** Error in subroutine read_block_15_0_0. ' // &
      'This version of block 15 requires 6 AOD wavelengths'
    ierr=i_err_fatal
    RETURN
  END IF

  nd_aod_wavel = Sp%Aerosol%n_aod_wavel
  ALLOCATE(Sp%Aerosol%i_aod_type(nd_aerosol_species))
  ALLOCATE(Sp%Aerosol%aod_wavel(nd_aod_wavel))
  ALLOCATE(Sp%Aerosol%aod_abs(nd_humidity, &
    nd_aerosol_species, nd_aod_wavel))
  ALLOCATE(Sp%Aerosol%aod_scat(nd_humidity, &
    nd_aerosol_species, nd_aod_wavel))

  Sp%Aerosol%aod_wavel(1:6) = 1.0e-9 * &
    (/ 380.0, 440.0, 550.0, 670.0, 865.0, 1020.0 /)
END IF

READ(iu_spc, '(19x,I5, 15x, I5)') i, j
Sp%Aerosol%i_aod_type(i) = j

READ(iu_spc, *)
DO k = 1 , Sp%Aerosol%n_aod_wavel
  READ(iu_spc, FMT='(I5,2(4X, 1PE12.5))', IOSTAT=ios, IOMSG=iomessage) &
     idum_band, Sp%Aerosol%aod_abs(1, i, k), &
     Sp%Aerosol%aod_scat(1, i, k)
  IF (ios /= 0) THEN
    cmessage = '*** Error in subroutine read_block_15_0_0.\n' // &
      'AOD data is not consistent with summary: ' // TRIM(iomessage)
    ierr=i_err_fatal
    RETURN
  END IF
END DO

END SUBROUTINE read_block_15_0_0


! Moist aerosol optical depth
SUBROUTINE read_block_15_1_0

IMPLICIT NONE

! Local variables.
INTEGER :: idum_band
!   Dummy integer
INTEGER :: j
!   AOD type
INTEGER :: k
!   Loop variable

READ(iu_spc, '(/14X, I5)') Sp%Aerosol%n_aod_wavel

IF ( .NOT. ALLOCATED(Sp%Aerosol%i_aod_type) ) THEN
  IF (Sp%Aerosol%n_aod_wavel /= 6) THEN
    cmessage = '*** Error in subroutine read_block_15_1_0.\n' // &
      'This version of block 15 requires 6 AOD wavelengths'
    ierr=i_err_fatal
    RETURN
  END IF

  nd_aod_wavel = Sp%Aerosol%n_aod_wavel
  ALLOCATE(Sp%Aerosol%i_aod_type(nd_aerosol_species))
  ALLOCATE(Sp%Aerosol%aod_wavel(nd_aod_wavel))
  ALLOCATE(Sp%Aerosol%aod_abs(nd_humidity, &
    nd_aerosol_species, nd_aod_wavel))
  ALLOCATE(Sp%Aerosol%aod_scat(nd_humidity, &
    nd_aerosol_species, nd_aod_wavel))

  Sp%Aerosol%aod_wavel(1:6) = 1.0e-9 * &
    (/ 380.0, 440.0, 550.0, 670.0, 865.0, 1020.0 /)
END IF

READ(iu_spc, '(19x,I5, 15x, I5)') i, j
Sp%Aerosol%i_aod_type(i) = j

DO k = 1 , Sp%Aerosol%n_aod_wavel
  READ(iu_spc, '(/)')
  DO j=1, Sp%Aerosol%nhumidity(i)
    READ(iu_spc, '(2(4X, 1PE12.5))', IOSTAT=ios, IOMSG=iomessage) &
      Sp%Aerosol%aod_abs(j, i,k), &
      Sp%Aerosol%aod_scat(j, i, k)
    IF (ios /= 0) THEN
      cmessage = '*** Error in subroutine ses_block_15_1_0.\n' // &
        'AOD data is not consistent with summary: ' // TRIM(iomessage)
      ierr=i_err_fatal
      RETURN
    END IF
  END DO

END DO

END SUBROUTINE read_block_15_1_0



SUBROUTINE read_block_16_0_0_int

IMPLICIT NONE

! Local variables.
INTEGER :: idum_band
!   Dummy integer
INTEGER :: i_band
!   Number of band
INTEGER :: j, i_term, i_mix, ip, it
!   Loop variables


! Allocate space.
ALLOCATE(Sp%Gas%mix_gas_band(nd_band))
ALLOCATE(Sp%Gas%n_mix_gas(nd_band))
ALLOCATE(Sp%Gas%index_mix_gas(2, nd_band_mix_gas))
! Note maximum number of mixed gases currently hardwired to 2.
ALLOCATE(Sp%Gas%k_mix_gas(nd_pre, nd_tmp, nd_mix, nd_k_term, &
                          nd_band_mix_gas))

READ(iu_spc, '(/)')

! Read in the number of ESFT terms in each band.
i = 0
Sp%Gas%mix_gas_band(1:nd_band)=0
DO i_band = 1 , Sp%Basic%n_band
  READ(iu_spc, FMT='(I5, 7X, I5)', IOSTAT=ios, IOMSG=iomessage) &
    idum_band, Sp%Gas%n_mix_gas(i_band)
  IF (ios /= 0) THEN
    cmessage = '*** Error in subroutine read_block_16_0_0.\n' // &
      'ESFT data is not consistent with summary: ' // TRIM(iomessage)
    ierr=i_err_fatal
    RETURN
  END IF

  IF ( Sp%Gas%n_mix_gas(i_band)  /=  0 ) THEN
    i=i+1

    Sp%Gas%mix_gas_band(i_band) = i

    READ(iu_spc,*)(Sp%Gas%index_mix_gas(j,i), &
      j = 1, Sp%Gas%n_mix_gas(i_band))

    READ(iu_spc1, '(/)')

    DO i_term = 1, Sp%Gas%i_band_k_ses(i_band)
      DO i_mix=1, Sp%Gas%num_mix(i_band)
        DO ip = 1, Sp%Gas%num_ref_p(1, i_band)
          READ(iu_spc1, *, err=1000, IOSTAT=ios, IOMSG=iomessage) &
          (Sp%Gas%k_mix_gas(ip,it,i_mix,i_term,i), it=1, nd_tmp)
        END DO
      END DO
    END DO

    1000    IF (ios /= 0) THEN
      WRITE(cmessage, '(a53, i3, a10, i3, a11, i3, a11, i3, A)') &
       'Error reading extended spectral file, block 16: ip = ', ip, &
       ': i_mix = ', i_mix, &
       ': i_term = ', i_term, &
       ': i_band = ', i_band, TRIM(iomessage)
      ierr=i_err_fatal
      RETURN
    END IF
  END IF
END DO

END SUBROUTINE read_block_16_0_0_int


! Spectral variability
SUBROUTINE read_block_17_0_0

IMPLICIT NONE

! Local variables.
INTEGER :: desc_end
!   Position of equals sign to delimit end of item description
INTEGER :: idum
!   Dummy integer
LOGICAL :: l_exist_var
!   True if spectral variability file exists and is readable
REAL (RealK), ALLOCATABLE :: rayleigh_coeff(:)

DO
  READ(iu_spc, '(a80)', IOSTAT=ios) line
  IF (line(1:4) == '*END') THEN
    BACKSPACE(iu_spc)
    EXIT
  END IF

  desc_end=SCAN(line,'=',.TRUE.)
  IF (desc_end==0) desc_end=LEN_TRIM(line)+1

  SELECT CASE (line(1:desc_end-1))
  CASE ('Number of spectral sub-bands','nd_sub_band')

    ! Read number of sub-bands
    READ(line(desc_end+1:),*,IOSTAT=ios, IOMSG=iomessage) Sp%Var%n_sub_band
    nd_sub_band=Sp%Var%n_sub_band

    ALLOCATE(Sp%Var%index_sub_band(2, nd_sub_band))
    ALLOCATE(Sp%Var%wavelength_sub_band(2, nd_sub_band))
    ALLOCATE(rayleigh_coeff(nd_sub_band))

    IF (Sp%Var%n_sub_band > Sp%Basic%n_band) THEN
      ! Skip header line
      READ(iu_spc, '(a80)', IOSTAT=ios, IOMSG=iomessage) line

      ! Read sub-band information
      DO i=1, Sp%Var%n_sub_band
        READ(iu_spc, '(3(i5, 2x), 2x, 1pe16.9, 2(4x, 1pe16.9))', &
          IOSTAT=ios, IOMSG=iomessage) &
          idum, Sp%Var%index_sub_band(:,i), Sp%Var%wavelength_sub_band(:,i), &
          rayleigh_coeff(i)
      END DO
    ELSE IF (Sp%Var%n_sub_band == Sp%Basic%n_band) THEN
      ! Sub-bands are equal to full-bands
      DO i=1, Sp%Var%n_sub_band
        Sp%Var%index_sub_band(1,i) = i
        Sp%Var%index_sub_band(2,i) = 0
        Sp%Var%wavelength_sub_band(1,i) = Sp%Basic%wavelength_short(i)
        Sp%Var%wavelength_sub_band(2,i) = Sp%Basic%wavelength_long(i)
        rayleigh_coeff(i) = Sp%Rayleigh%rayleigh_coeff(i)
      END DO
    ELSE
      cmessage = 'Not enough sub-bands in block 17 (should be >= n_band).'
      ierr=i_err_fatal
      RETURN
    END IF

    i=INDEX(file_spectral, ' ') - 1
    spectral_var = file_spectral(1:i) // '_var'

    ! Get a unit to read the spectral variability file.
    CALL assign_file_unit(spectral_var, iu_spc2, handler="fortran")

    ! Open the file for reading
    OPEN(UNIT=iu_spc2, FILE=spectral_var, IOSTAT=ios, STATUS='OLD', &
         ACTION='READ', IOMSG=iomessage)
    IF (ios == 0) THEN
      l_exist_var = .TRUE.
    ELSE
      l_exist_var = .FALSE.
      ios = 0
    END IF

    Sp%Var%n_times          = 0
    nd_times                = 0
    Sp%Var%n_repeat_times   = 0
    Sp%Var%n_rayleigh_coeff = 0
    IF (l_exist_var) THEN
      ! Read number of times / dates and skip over header
      DO
        READ(iu_spc2, '(a80)', IOSTAT=ios, IOMSG=iomessage) line
        IF (ios /= 0) THEN
           cmessage = 'Error reading spectral variability file: ' // &
             TRIM(iomessage)
           ierr=i_err_fatal
           RETURN
        END IF
        desc_end=SCAN(line,'=',.TRUE.)
        IF ( line(1:desc_end-1) == 'Number of times in look-up table' .OR. &
             line(1:desc_end-1) == 'nd_times' ) THEN
          READ(line(desc_end+1:),*,IOSTAT=ios, IOMSG=iomessage) Sp%Var%n_times
          nd_times=Sp%Var%n_times
        END IF
        IF ( line(1:desc_end-1) == 'Number of times for periodic repetition' ) &
          READ(line(desc_end+1:),*,IOSTAT=ios, IOMSG=iomessage) &
            Sp%Var%n_repeat_times
        IF ( line(1:desc_end-1) == 'Number of Rayleigh coefficients given' ) &
          READ(line(desc_end+1:),*,IOSTAT=ios, IOMSG=iomessage) &
            Sp%Var%n_rayleigh_coeff
        IF ( line(1:6) == '*BEGIN' ) EXIT
      END DO
      
      ! Read look-up table of spectral variability data
      ALLOCATE(Sp%Var%time(4, nd_times))
      ALLOCATE(Sp%Var%total_solar_flux(nd_times))
      ALLOCATE(Sp%Var%solar_flux_sub_band(nd_sub_band, nd_times))
      ALLOCATE(Sp%Var%rayleigh_coeff(nd_sub_band, 0:nd_times))
      DO i=1, Sp%Var%n_times
        READ(iu_spc2, '(4(i6),4x,1pe16.9)') &
          Sp%Var%time(:, i), Sp%Var%total_solar_flux(i)
        READ(iu_spc2, '(5(1pe16.9))') Sp%Var%solar_flux_sub_band(:, i)
        IF (Sp%Var%n_rayleigh_coeff > 0) READ(iu_spc2, '(5(1pe16.9))') &
          Sp%Var%rayleigh_coeff(1:Sp%Var%n_rayleigh_coeff, i)
      END DO
      
      CLOSE(iu_spc2)
    ELSE
      ALLOCATE(Sp%Var%rayleigh_coeff(nd_sub_band, 0:nd_times))
    END IF

    Sp%Var%rayleigh_coeff(:,0) = rayleigh_coeff
    DO i=1, Sp%Var%n_times
      Sp%Var%rayleigh_coeff(Sp%Var%n_rayleigh_coeff+1:,i) = &
             rayleigh_coeff(Sp%Var%n_rayleigh_coeff+1:)
    END DO
    DEALLOCATE(rayleigh_coeff)

    CALL release_file_unit(iu_spc2, handler="fortran")
  END SELECT

  IF (ios /= 0) THEN
    cmessage = 'Error in subroutine read_block_17_0_0: ' // &
               TRIM(iomessage)
    ierr=i_err_fatal
    RETURN
  END IF
END DO

END SUBROUTINE read_block_17_0_0


END SUBROUTINE read_spectrum
