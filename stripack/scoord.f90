SUBROUTINE SCOORD(PX, PY, PZ, PLAT, PLON, PNRM)
    REAL PX, PY, PZ, PLAT, PLON, PNRM
!
!***********************************************************
!   This subroutine converts a point P from Cartesian coor-
! dinates to spherical coordinates.
!
!
! On input:
!
!       PX,PY,PZ = Cartesian coordinates of P.
!
! Input parameters are not altered by this routine.
!
! On output:
!
!       PLAT = Latitude of P in the range -PI/2 to PI/2, or
!              0 if PNRM = 0.  PLAT should be scaled by
!              180/PI to obtain the value in degrees.
!
!       PLON = Longitude of P in the range -PI to PI, or 0
!              if P lies on the Z-axis.  PLON should be
!              scaled by 180/PI to obtain the value in
!              degrees.
!
!       PNRM = Magnitude (Euclidean norm) of P.
!
!***********************************************************
!
    PNRM = SQRT(PX * PX + PY * PY + PZ * PZ)
    IF (PX .NE. 0. .OR. PY .NE. 0.) THEN
        PLON = ATAN2(PY, PX)
    ELSE
        PLON = 0.
    END IF
    IF (PNRM .NE. 0.) THEN
        PLAT = ASIN(PZ / PNRM)
    ELSE
        PLAT = 0.
    END IF
    RETURN
END
