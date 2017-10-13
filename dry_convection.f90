MODULE dry_convection_mod

  IMPLICIT NONE


  PUBLIC :: dry_convection

CONTAINS


  SUBROUTINE dry_convection(tin, p_full, &
       Tref)
    REAL, INTENT(in), DIMENSION(:) :: tin, p_full
    REAL, INTENT(out), DIMENSION(:) :: Tref
    REAL :: T1, P1, P2, T2, Tbar, DTref, deltap1, deltap2, pfact,Rcp
    INTEGER :: n, k

    Rcp = 0.285
    n = SIZE(tin)
    Tref = tin *1.
       !Downward pass
    DO k=1,n-1
       T1 = Tref(k)
       T2 = Tref(k+1)
       P1 = p_full(k)
       P2 = p_full(k+1)

       pfact = (P1/P2)**Rcp


       IF (T1 .LT. T2*pfact) THEN

          Tbar = 0.5 * (T1+T2)
          T2 = 2.*Tbar/(1.+pfact)
          T1 = T2*pfact
          Tref(k) = T1
          Tref(k+1) = T2

       END IF
    END DO

       !Upward pass
    DO k=n-1,1,-1
       T1 = Tref(k)
       T2 = Tref(k+1)
       P1 = p_full(k)
       P2 = p_full(k+1)

       pfact = (P1/P2)**Rcp


       IF (T1 .LT. T2*pfact) THEN

          Tbar = 0.5 * (T1+T2)
          T2 = 2.*Tbar/(1.+pfact)
          T1 = T2*pfact
          Tref(k) = T1
          Tref(k+1) = T2

       END IF
    END DO
  END SUBROUTINE dry_convection


END MODULE dry_convection_mod
