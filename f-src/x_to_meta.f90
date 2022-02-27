PROGRAM xtometa
!>-------------------------------------------------------------------
!> vtk to raw converter
!
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!> Date:    12.09.2021
!> LastMod: 28.12.2021
!>-------------------------------------------------------------------

USE ISO_FORTRAN_ENV
USE global_std
USE user_interaction
USE meta
USE MPI
USE vtk_meta_data
USE raw_binary  

IMPLICIT NONE
! Parameter
INTEGER(KIND=ik), PARAMETER :: debug = 1

! MPI variables
INTEGER(KIND=mik) :: ierr, my_rank, size_mpi

! Std variables
CHARACTER(LEN=mcl), DIMENSION(:), ALLOCATABLE :: m_rry      
CHARACTER(LEN=scl) :: type_in, type_out, binary,  restart_cmd_arg, dtrep=''
CHARACTER(LEN=mcl) :: filename='', cmd_arg_history='', file_type_in = ''
CHARACTER(LEN=4  ) :: suffix=''

INTEGER(KIND=ik) :: hdr
INTEGER(KIND=mik), DIMENSION(3) :: sections
INTEGER(KIND=ik), DIMENSION(3) :: dims, rry_dims, rank_section, sections_ik
INTEGER(KIND=ik), DIMENSION(3) :: remainder_per_dir, dims_reduced, subarray_origin
REAL   (KIND=rk), DIMENSION(3) :: spcng, origin
REAL   (KIND=rk) :: start, end

! Binary blob variables1
REAL   (KIND=REAL32), DIMENSION(:,:,:), ALLOCATABLE :: rry_rk4
REAL   (KIND=REAL64), DIMENSION(:,:,:), ALLOCATABLE :: rry_rk8
INTEGER(KIND=INT16) , DIMENSION(:,:,:), ALLOCATABLE :: rry_ik2
INTEGER(KIND=INT32) , DIMENSION(:,:,:), ALLOCATABLE :: rry_ik4

LOGICAL :: stp = .FALSE.

!------------------------------------------------------------------------------
! Invoke MPI 
!------------------------------------------------------------------------------
CALL mpi_init(ierr)
CALL print_err_stop(std_out, "MPI_INIT didn't succeed", INT(ierr, KIND=ik))

CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
CALL print_err_stop(std_out, "MPI_COMM_RANK couldn't be retrieved", INT(ierr, KIND=ik))

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size_mpi, ierr)
CALL print_err_stop(std_out, "MPI_COMM_SIZE couldn't be retrieved", INT(ierr, KIND=ik))

IF (size_mpi < 2) CALL print_err_stop(std_out, "We need at least 2 MPI processes to execute this program.", 1)

!------------------------------------------------------------------------------
! Rank 0 -- Init (Master) Process and broadcast init parameters 
!------------------------------------------------------------------------------
IF (my_rank==0) THEN

    CALL CPU_TIME(start)

    !------------------------------------------------------------------------------
    ! Parse the command arguments
    ! restart_cmd_arg not used in fmf, since lock file handling is not relevant.
    ! Implemented here to sustain api compatibility (maybe a bad API :-)
    !------------------------------------------------------------------------------
    CALL get_cmd_args(binary, filename, stp, restart_cmd_arg, cmd_arg_history)
    IF(stp) GOTO 1001
    
    IF (filename=='') THEN
        CALL usage(binary)    

        !------------------------------------------------------------------------------
        ! On std_out since file of std_out is not spawned
        !------------------------------------------------------------------------------
        CALL print_err_stop(6, "No input file given", 1)
    END IF

    !------------------------------------------------------------------------------
    ! Redirect std_out into a file in case std_out is not useful by environment.
    !------------------------------------------------------------------------------
    std_out = determine_stout()

    !------------------------------------------------------------------------------
    ! Spawn standard out after(!) the basename is known
    !------------------------------------------------------------------------------
    IF(std_out/=6) CALL meta_start_ascii(std_out, '.std_out')

    CALL show_title()
 
    IF(debug >=0) WRITE(std_out, FMT_MSG) "Post mortem info probably in ./datasets/temporary.std_out"
    WRITE(std_out, FMT_TXT) "Program invocation:"//TRIM(cmd_arg_history)          

    !------------------------------------------------------------------------------
    ! Check whether to convert from vtk to meta or the opposite way.
    !------------------------------------------------------------------------------
    IF (filename(LEN_TRIM(filename, KIND=ik)-2:LEN_TRIM(filename, KIND=ik)) == "vtk") THEN
        CALL meta_create_new(TRIM(ADJUSTL(filename)))
 
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
            CASE('rk4', 'rk8', 'ik4', 'uik2') ; type_out = 'ik4'
            CASE('ik2') ; type_out = 'ik2'
        END SELECT

        !------------------------------------------------------------------------------
        ! Write meta data 
        !------------------------------------------------------------------------------
        WRITE(std_out, FMT_TXT) 'Writing meta information to *.meta file.'

        WRITE(fhmeo, '(A)') "# HLRS|NUM Dataset Meta Information"
        WRITE(fhmeo, '(A)') ""
        WRITE(fhmeo, '(A)') "* GENERAL_INFORMATION"
        CALL meta_write('CT_SCAN', in%dataset)
        CALL meta_write('OWNER',         "TBD by user")
        CALL meta_write('OWNER_CONTACT', "TBD by user")
        CALL meta_write('DATE_CREATED',  "TBD by user")
        CALL meta_write('INTERNAL_ID',   "TBD by user")
        CALL meta_write('HISTORY', '(-)' , 1)
        
        ! Write original data type to file for documentaroy purposes
        CALL meta_write('TYPE_IMPORT', TRIM(ADJUSTL(type_in)))
        CALL meta_write('TYPE_RAW', TRIM(ADJUSTL(type_out)))
        
        !------------------------------------------------------------------------------
        ! Data always are written as little endian; only vtk requires big endian.
        !------------------------------------------------------------------------------
        CALL meta_write('DATA_BYTE_ORDER'  , "LittleEndian")
        CALL meta_write('DIMENSIONALITY'   , '(-)'  , 3)
        CALL meta_write('DIMENSIONS'       , '(-)'  , dims)
        CALL meta_write('NO_SCALAR_CMPNNTS', '(-)'  , 1)
        CALL meta_write('SPACING'          , '(mm)' , spcng)
        CALL meta_write('ORIGIN_SHIFT_GLBL', '(mm)' , [0._rk, 0._rk, 0._rk])
        CALL meta_write('ORIGIN'           , '(-)'  , [0, 0, 0])
        CALL meta_write('FIELD_OF_VIEW'    , '(mm)' , dims*spcng)
        CALL meta_write('ENTRIES'          , '(-)'  , PRODUCT(dims))
    
        FLUSH(fhmeo)

        file_type_in = "vtk"
    ELSE IF (filename(LEN_TRIM(filename, KIND=ik)-3:LEN_TRIM(filename, KIND=ik)) == "meta") THEN
        
        !------------------------------------------------------------------------------
        ! Special treatment of meta format, not a regular i/o procedure.
        !------------------------------------------------------------------------------
        in%full = TRIM(ADJUSTL(filename))
        in%p_n_bsnm = filename(1:LEN_TRIM(ADJUSTL(filename))-5)

        out = in

        CALL meta_invoke(m_rry)

        CALL meta_read('TYPE_RAW'  , m_rry, type_out)
        CALL meta_read('SPACING'   , m_rry, spcng)
        CALL meta_read('ORIGIN'    , m_rry, origin)
        CALL meta_read('DIMENSIONS', m_rry, dims)

        CALL write_vtk_struct_points_header &
            (fh_vtk, TRIM(out%p_n_bsnm)//vtk_suf, type_out, spcng, origin, dims)

        CLOSE (fhmei)

        file_type_in = "meta"
        hdr = 0_ik
        type_in = type_out

    ELSE
        mssg = "No valid input file (*.meta or *.vtk) given."
        CALL print_err_stop(std_out, mssg, 1_ik)
    END IF 

END IF ! my_rank==0

!------------------------------------------------------------------------------
! Send required variables
!------------------------------------------------------------------------------
CALL MPI_BCAST (type_in     , INT(scl, KIND=mik)     , MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST (type_out    , INT(scl, KIND=mik)     , MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST (in%p_n_bsnm , INT(meta_mcl, KIND=mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST (out%p_n_bsnm, INT(meta_mcl, KIND=mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST (file_type_in, INT(meta_mcl, KIND=mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST (hdr , 1_mik, MPI_INTEGER8, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST (dims, 3_mik, MPI_INTEGER8, 0_mik, MPI_COMM_WORLD, ierr)

!------------------------------------------------------------------------------
! Get dimensions for each domain. Every processor reveives its own domain.
! Therefore, each my_rank calculates its own address/dimensions/parameters.
!
! The decomposition of an image depends on the task to be done. Therefore, 
! these calculations are not encapsuled in a routine. However, this may
! follow in an update during refactoring the 3D scalar image filter.
!
! Allocation of subarray memory is done in the read_raw routines.
!------------------------------------------------------------------------------
sections=0
CALL MPI_DIMS_CREATE (size_mpi, 3_mik, sections, ierr)
CALL get_rank_section(INT(my_rank, KIND=ik), INT(sections, KIND=ik), rank_section)

sections_ik = INT(sections, KIND=ik)
remainder_per_dir = MODULO(dims, sections_ik)

dims_reduced   = dims - remainder_per_dir

rry_dims  = (dims_reduced / sections_ik)

subarray_origin = (rank_section-1_ik) * (rry_dims)

! Add the remainder to the last domains of each dimension
IF(rank_section(1) == sections_ik(1)) rry_dims(1) = rry_dims(1) + remainder_per_dir(1)
IF(rank_section(2) == sections_ik(2)) rry_dims(2) = rry_dims(2) + remainder_per_dir(2)
IF(rank_section(3) == sections_ik(3)) rry_dims(3) = rry_dims(3) + remainder_per_dir(3)

!------------------------------------------------------------------------------
! Read binary part of the vtk file - basically a *.raw file
!------------------------------------------------------------------------------
IF(my_rank==0) WRITE(std_out, FMT_TXT) 'Reading binary information of *.'//TRIM(file_type_in)//' file.'

IF (file_type_in == "vtk") THEN
    suffix = vtk_suf
    dtrep = 'EXTERNAL32'
ELSE IF (file_type_in == "meta") THEN
    suffix = raw_suf
    dtrep = 'NATIVE'
END IF 
 
SELECT CASE(type_in)
    CASE('rk4') 
        ! MPI_OFFSET_KIND needs ik=8 in this case.
        CALL mpi_read_raw(TRIM(in%p_n_bsnm)//suffix, INT(hdr, KIND=8), dims, rry_dims, subarray_origin, rry_rk4, TRIM(dtrep))
    CASE('rk8') 
        CALL mpi_read_raw(TRIM(in%p_n_bsnm)//suffix, INT(hdr, KIND=8), dims, rry_dims, subarray_origin, rry_rk8, TRIM(dtrep))
    CASE('ik4') 
        CALL mpi_read_raw(TRIM(in%p_n_bsnm)//suffix, INT(hdr, KIND=8), dims, rry_dims, subarray_origin, rry_ik4, TRIM(dtrep))
    CASE('ik2', 'uik2') 
        CALL mpi_read_raw(TRIM(in%p_n_bsnm)//suffix, INT(hdr, KIND=8), dims, rry_dims, subarray_origin, rry_ik2, TRIM(dtrep))
        IF(type_in=='uik2') THEN
            CALL uik2_to_ik4(rry_ik2, rry_ik4)
            DEALLOCATE(rry_ik2)
        END IF 
END SELECT

!------------------------------------------------------------------------------
! Convert arrays and write raw data
!------------------------------------------------------------------------------
IF(my_rank==0) WRITE(std_out, FMT_TXT) 'Converting and writing binary information to *.raw file.'

!------------------------------------------------------------------------------
! interchange suffixes to convert vtk/meta or meta/vtk
!------------------------------------------------------------------------------
IF (file_type_in == "vtk") THEN
    suffix = raw_suf
    hdr = 0_8
    dtrep = 'NATIVE'
ELSE IF (file_type_in == "meta") THEN
    suffix = vtk_suf
    dtrep = 'EXTERNAL32'

    !------------------------------------------------------------------------------
    ! Get the size of the meta header
    !------------------------------------------------------------------------------
    INQUIRE(FILE=TRIM(out%p_n_bsnm)//suffix, SIZE=hdr)
END IF 

SELECT CASE(type_out)
    CASE('ik2') 
        CALL mpi_write_raw(TRIM(out%p_n_bsnm)//suffix, INT(hdr, KIND=8), dims, &
            rry_dims, subarray_origin, rry_ik2, TRIM(dtrep))
    CASE('ik4') 
        SELECT CASE(type_in)
            CASE('rk4') 
                rry_ik4 = INT(rry_rk4, KIND=INT32)
            CASE('rk8') 
                rry_ik4 = INT(rry_rk8, KIND=INT32)
        END SELECT

        CALL mpi_write_raw(TRIM(out%p_n_bsnm)//suffix, INT(hdr, KIND=8), dims, &
            rry_dims, subarray_origin, rry_ik4, TRIM(dtrep))
END SELECT

!------------------------------------------------------------------------------
! Only used in specific cases to finish more gracefully. (grep -i "GOTO")
!------------------------------------------------------------------------------
1001 Continue

!------------------------------------------------------------------------------
! Finish program
!------------------------------------------------------------------------------
IF(my_rank == 0) THEN
    CALL CPU_TIME(end)

    WRITE(std_out, FMT_TXT_xAF0) 'Finishing the program took', end-start,'seconds.'
    WRITE(std_out, FMT_TXT_SEP)

    IF (file_type_in == "vtk") THEN
        CALL meta_signing(binary)
        CALL meta_close(size_mpi)

    ELSE IF (file_type_in == "meta") THEN
        CALL write_vtk_struct_points_footer (fh_vtk, TRIM(out%p_n_bsnm)//suffix)
    
    END IF

    IF (std_out/=6) CALL meta_stop_ascii(fh=std_out, suf='.std_out')

END IF ! (my_rank == 0)

Call MPI_FINALIZE(ierr)
CALL print_err_stop(std_out, "MPI_FINALIZE didn't succeed", INT(ierr, KIND=ik))

END PROGRAM xtometa
