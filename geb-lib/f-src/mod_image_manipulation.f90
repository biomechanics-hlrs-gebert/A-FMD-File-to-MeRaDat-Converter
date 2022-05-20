!------------------------------------------------------------------------------
! MODULE: image_manipulation
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @description: 
!> Collection of image manipulating routines.
!------------------------------------------------------------------------------
MODULE image_manipulation

USE ISO_FORTRAN_ENV
USE global_std

IMPLICIT NONE

INTERFACE downscale
    MODULE PROCEDURE downscale_ik2
    MODULE PROCEDURE downscale_ik4
END INTERFACE downscale

CONTAINS

!------------------------------------------------------------------------------
! SUBROUTINE: downscale_ik2
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Reduce the resolution of a given image by a given integer factor.
!
!> @param[in] in_array Input array
!> @param[in] scale_factor Factor of high resolution voxels to get a mean ave.
!> @param[out] out_array Output array
!------------------------------------------------------------------------------  
SUBROUTINE downscale_ik2(in_array, scale_factor, out_array)

    INTEGER(KIND=INT16), DIMENSION(:,:,:), INTENT(OUT) :: out_array
    INTEGER(KIND=INT16), DIMENSION(:,:,:), INTENT(IN) :: in_array
    INTEGER(KIND=ik), DIMENSION(3), INTENT(IN) :: scale_factor

    INTEGER(KIND=ik), DIMENSION(3) :: shp_in
    INTEGER(KIND=ik) :: ii, jj, kk, ll, mm, nn
    REAL(KIND=rk) :: scale_volume

    scale_volume = REAL(PRODUCT(scale_factor), KIND=rk)
    shp_in = SHAPE(in_array)

    nn = 1_ik
    DO kk=1, SIZE(in_array, DIM=3)-scale_factor(3), scale_factor(3)
        IF (shp_in(3) < kk+scale_factor(3)) CYCLE

        mm = 1_ik
        DO jj=1, SIZE(in_array, DIM=2)-scale_factor(2), scale_factor(2)
            IF (shp_in(2) < jj+scale_factor(2)) CYCLE

            ll = 1_ik
            DO ii=1, SIZE(in_array, DIM=1)-scale_factor(1), scale_factor(1)
                IF (shp_in(1) < ii+scale_factor(1)) CYCLE

                !------------------------------------------------------------------------------  
                ! DO NOT interchange SUM and REAL as the sum will overflow without conversion.
                !------------------------------------------------------------------------------  
                out_array(ll, mm, nn) = INT(SUM(REAL(in_array(&
                    ii:ii+scale_factor(1)-1_ik, &
                    jj:jj+scale_factor(2)-1_ik, &
                    kk:kk+scale_factor(3)-1_ik), KIND=rk)) / scale_volume, KIND=INT16)
                ll = ll + 1_ik
            END DO

            mm = mm + 1_ik
        END DO

        nn = nn + 1_ik
    END DO
END SUBROUTINE downscale_ik2

!------------------------------------------------------------------------------
! SUBROUTINE: downscale_ik4
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Reduce the resolution of a given image by a given integer factor.
!
!> @param[in] in_array Input array
!> @param[in] scale_factor Factor of high resolution voxels to get a mean ave.
!> @param[out] out_array Output array
!------------------------------------------------------------------------------  
SUBROUTINE downscale_ik4(in_array, scale_factor, out_array)

    INTEGER(KIND=INT32), DIMENSION(:,:,:), INTENT(OUT) :: out_array
    INTEGER(KIND=INT32), DIMENSION(:,:,:), INTENT(IN) :: in_array
    INTEGER(KIND=ik), DIMENSION(3), INTENT(IN) :: scale_factor

    INTEGER(KIND=ik), DIMENSION(3) :: shp_in
    INTEGER(KIND=ik) :: ii, jj, kk, ll, mm, nn
    REAL(KIND=rk) :: scale_volume

    scale_volume = REAL(PRODUCT(scale_factor), KIND=rk)
    shp_in = SHAPE(in_array)

    nn = 1_ik
    DO kk=1, SIZE(in_array, DIM=3)-scale_factor(3), scale_factor(3)
        IF (shp_in(3) < kk+scale_factor(3)) CYCLE

        mm = 1_ik
        DO jj=1, SIZE(in_array, DIM=2)-scale_factor(2), scale_factor(2)
            IF (shp_in(2) < jj+scale_factor(2)) CYCLE

            ll = 1_ik
            DO ii=1, SIZE(in_array, DIM=1)-scale_factor(1), scale_factor(1)
                IF (shp_in(1) < ii+scale_factor(1)) CYCLE

                !------------------------------------------------------------------------------  
                ! DO NOT interchange SUM and REAL as the sum will overflow without conversion.
                !------------------------------------------------------------------------------  
                out_array(ll, mm, nn) = INT(SUM(REAL(in_array(&
                    ii:ii+scale_factor(1)-1_ik, &
                    jj:jj+scale_factor(2)-1_ik, &
                    kk:kk+scale_factor(3)-1_ik), KIND=rk)) / scale_volume, KIND=INT32)
                ll = ll + 1_ik
            END DO

            mm = mm + 1_ik
        END DO

        nn = nn + 1_ik
    END DO
END SUBROUTINE downscale_ik4

END MODULE image_manipulation