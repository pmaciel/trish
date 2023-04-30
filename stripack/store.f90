REAL FUNCTION STORE(X)
    REAL X
!
!***********************************************************
!   This function forces its argument X to be stored in a
! memory location, thus providing a means of determining
! floating point number characteristics (such as the machine
! precision) when it is necessary to avoid computation in
! high precision registers.
!
!
! On input:
!
!       X = Value to be stored.
!
! X is not altered by this function.
!
! On output:
!
!       STORE = Value of X after it has been stored and
!               possibly truncated or rounded to the single
!               precision word length.
!
! Modules required by STORE:  None
!
!***********************************************************
!
    REAL Y
    COMMON / STCOM / Y
    Y = X
    STORE = Y
    RETURN
END
