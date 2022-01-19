PROGRAM xtometa
!>-------------------------------------------------------------------
!> vtk to raw converter
!
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!> Date:    12.09.2021
!> LastMod: 20.01.2022
!>-------------------------------------------------------------------

USE ISO_FORTRAN_ENV
USE global_std
USE user_interaction
USE meta
USE vtk_meta_data
USE raw_binary  

IMPLICIT NONE
! Parameter
INTEGER(KIND=ik), PARAMETER :: debug = 1

! Std variables
CHARACTER(LEN=scl) :: type_in, type_out, binary, restart, restart_cmd_arg
CHARACTER(LEN=mcl) :: filename=''

INTEGER(KIND=ik) :: hdr
INTEGER(KIND=ik), DIMENSION(3) :: dims
REAL   (KIND=rk), DIMENSION(3) :: spcng, origin
REAL   (KIND=rk) :: start, end

! Binary blob variables
REAL   (KIND=REAL32), DIMENSION(:,:,:), ALLOCATABLE :: rry_rk4
REAL   (KIND=REAL64), DIMENSION(:,:,:), ALLOCATABLE :: rry_rk8
INTEGER(KIND=INT16) , DIMENSION(:,:,:), ALLOCATABLE :: rry_ik2
INTEGER(KIND=INT32) , DIMENSION(:,:,:), ALLOCATABLE :: rry_ik4

LOGICAL :: stp = .FALSE.

CALL CPU_TIME(start)

!------------------------------------------------------------------------------
! Parse the command arguments
!------------------------------------------------------------------------------
CALL get_cmd_args(binary, filename, stp, restart, restart_cmd_arg)
IF(stp) GOTO 1001

IF (filename=='') THEN
    CALL usage(binary)    

    !------------------------------------------------------------------------------
    ! On std_out since file of std_out is not spawned
    !------------------------------------------------------------------------------
    CALL print_err_stop(6, "No input file given", 1)
END IF

CALL meta_create_new(TRIM(ADJUSTL(filename)))

!------------------------------------------------------------------------------
! Redirect std_out into a file in case std_out is not useful by environment.
! Place these lines before handle_lock_file :-)
!------------------------------------------------------------------------------
std_out = determine_stout()

!------------------------------------------------------------------------------
! Spawn standard out after(!) the basename is known
!------------------------------------------------------------------------------
IF(std_out/=6) CALL meta_start_ascii(std_out, '.std_out')

CALL show_title()

IF(debug >=0) WRITE(std_out, FMT_MSG) "Post mortem info probably in ./datasets/.temporary.std_out"

!------------------------------------------------------------------------------
! Read VTK file header
!------------------------------------------------------------------------------
WRITE(std_out, FMT_TXT) 'Reading the vtk header.'
CALL read_vtk_meta (filename, hdr, dims, origin, spcng, type_in) 

!------------------------------------------------------------------------------
! Determine data types
! Data types are converted to integer, since the provided precision on the
! data is sufficient.
!------------------------------------------------------------------------------
WRITE(std_out, FMT_TXT) 'Determining the output data types.'
SELECT CASE(type_in)
    CASE('rk4') ; type_out = 'ik4'
    CASE('rk8') ; type_out = 'ik4'
    CASE('ik4') ; type_out = 'ik4'
    CASE('ik2') ; type_out = 'ik2'
    CASE('uik2'); type_out = 'ik4'
END SELECT

!------------------------------------------------------------------------------
! Write meta data 
!------------------------------------------------------------------------------
WRITE(std_out, FMT_TXT) 'Writing meta information to *.meta file.'

WRITE(fhmeo, '(A)') "# HLRS|NUM Dataset Meta Information"
WRITE(fhmeo, '(A)') ""
WRITE(fhmeo, '(A)') "* GENERAL_INFORMATION"
CALL meta_write (fhmeo, 'CT_SCAN', in%dataset)
CALL meta_write (fhmeo, 'OWNER',         "TBD by user")
CALL meta_write (fhmeo, 'OWNER_CONTACT', "TBD by user")
CALL meta_write (fhmeo, 'DATE_CREATED',  "TBD by user")
CALL meta_write (fhmeo, 'INTERNAL_ID',   "TBD by user")
CALL meta_write (fhmeo, 'HISTORY', '(-)' , 1)

! Write original data type to file for documentaroy purposes
CALL meta_write (fhmeo, 'TYPE_IMPORT', TRIM(ADJUSTL(type_in)))
CALL meta_write (fhmeo, 'TYPE_RAW', TRIM(ADJUSTL(type_out)))

!------------------------------------------------------------------------------
! Data always are written as little endian; only vtk requires big endian.
!------------------------------------------------------------------------------
CALL meta_write (fhmeo, 'DATA_BYTE_ORDER'  , "LittleEndian")
CALL meta_write (fhmeo, 'DIMENSIONALITY'   , '(-)'  , 3)
CALL meta_write (fhmeo, 'DIMENSIONS'       , '(-)'  , dims)
CALL meta_write (fhmeo, 'NO_SCALAR_CMPNNTS', '(-)'  , 1)
CALL meta_write (fhmeo, 'SPACING'          , '(mm)' , spcng)
CALL meta_write (fhmeo, 'ORIGIN_SHIFT_GLBL', '(mm)' , [0._rk, 0._rk, 0._rk])
CALL meta_write (fhmeo, 'ORIGIN'           , '(-)'  , [0, 0, 0])
CALL meta_write (fhmeo, 'FIELD_OF_VIEW'    , '(mm)' , dims*spcng)
CALL meta_write (fhmeo, 'ENTRIES'          , '(-)'  , PRODUCT(dims))

FLUSH(fhmeo)

!------------------------------------------------------------------------------
! Read binary part of the vtk file - basically a *.raw file
!------------------------------------------------------------------------------
WRITE(std_out, FMT_TXT) 'Reading binary information of *.vtk file.'


SELECT CASE(type_in)
    CASE('rk4'); CALL ser_read_raw(fh_vtk, TRIM(in%p_n_bsnm)//vtk_suf, rry_rk4, 'SWAP')
    CASE('rk8'); CALL ser_read_raw(fh_vtk, TRIM(in%p_n_bsnm)//vtk_suf, rry_rk8, 'SWAP')
    CASE('ik4'); CALL ser_read_raw(fh_vtk, TRIM(in%p_n_bsnm)//vtk_suf, rry_ik4, 'SWAP')
    CASE('ik2', 'uik2') 
        CALL ser_read_raw(fh_vtk, TRIM(in%p_n_bsnm)//vtk_suf, rry_ik2, 'SWAP')
        IF(type_in=='uik2') THEN
            CALL uik2_to_ik4(rry_ik2, rry_ik4)
            DEALLOCATE(rry_ik2)
        END IF 
END SELECT

!------------------------------------------------------------------------------
! Convert arrays and write raw data
!------------------------------------------------------------------------------
WRITE(std_out, FMT_TXT) 'Converting and writing binary information to *.raw file.'

SELECT CASE(type_out)
    CASE('ik2') 
        CALL ser_write_raw(fh_vtk, TRIM(out%p_n_bsnm)//raw_suf, rry_ik2)
    CASE('ik4') 
        SELECT CASE(type_in)
            CASE('rk4') 
                rry_ik4 = INT(rry_rk4, KIND=INT32)
            CASE('rk8') 
                rry_ik4 = INT(rry_rk8, KIND=INT32)
        END SELECT

        CALL ser_write_raw(fh_vtk, TRIM(out%p_n_bsnm)//raw_suf, rry_ik4)
END SELECT

!------------------------------------------------------------------------------
! Only used in specific cases to finish more gracefully. (grep -i "GOTO")
!------------------------------------------------------------------------------
1001 Continue

!------------------------------------------------------------------------------
! Finish program
!------------------------------------------------------------------------------
CALL CPU_TIME(end)

WRITE(std_out, FMT_TXT_xAF0) 'Finishing the program took', end-start,'seconds.'
WRITE(std_out, FMT_TXT_SEP)

CALL meta_signing(binary)
CALL meta_close(INT(1, KIND=INT32)) ! MPI legacy

IF (std_out/=6) CALL meta_stop_ascii(fh=std_out, suf='.std_out')

END PROGRAM xtometa
