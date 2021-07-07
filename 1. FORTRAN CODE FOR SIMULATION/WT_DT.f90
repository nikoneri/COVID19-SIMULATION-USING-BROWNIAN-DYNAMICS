SUBROUTINE WDT(CXR,CYR,KI,KI2,DT,NFREC2,NN1,NN2,N)
!*********************************************************************************************************************
!   1. AUTHORS AND E-MAIL CONTACT:
!   -Isaac Neri Gomez Sarmiento (isaacneri.gs@gmail.com)
!   -César Omar Ramírez Álvarez (cesaromarramirezalvarez@gmail.com)
!   -Jonas Valenzuela Terán     (24jonass@gmail.com)
!   -We thank the support and advice of our molecular simulations teacher Dr. Laura Lorenia Yeomans Reyna
!*********************************************************************************************************************
!   2. DESCRIPTION:
!       THIS SUBROUTINE ALLOWS US TWO THINGS:
!           1.Compute the mean square displacement function W(t),which gives us an insight on how much area the particles spread out.
!           2.Diffussion Coefficient D(t), which gives us an insight on the mobility (speed) of particles.
!*********************************************************************************************************************

!********************
!VARIABLE DECLARATION
!********************
IMPLICIT NONE

INTEGER,INTENT(IN) :: N !Number of particles
INTEGER,INTENT(IN) :: NN1 !Maximum row size for particles' properties matrices: NN1=N+10
INTEGER,INTENT(IN) :: NN2 !Number of saved configurations or simulation steps
INTEGER,INTENT(IN) :: KI !Counter for the number of saved simulation steps for position matrices with periodic boundary conditions
INTEGER,INTENT(IN) :: KI2 !Counter for the number of saved simulation steps for the position matrices that will be used for computing g(r), W(t) and D(t)
INTEGER,INTENT(IN) :: NFREC2 !Frecuency at which we'll be saving configurations for computing W(t) and D(t) (i.e same as NFREC)
REAL,DIMENSION(NN1,NN2), INTENT(IN) :: CXR, CYR !Two dimensional arrays. Rows: particles, Columns:position at saved steps (no periodic boundary conditions)
REAL, INTENT(IN) :: DT !Timestep
INTEGER:: i, L, j, NTMAX !Counters
REAL(8) :: TIM !Time equivalent for the frecuency for saving configurations (TIM=NFrec2*DT)
REAL(8) :: TIME !Time equivalent to simulation steps TIME=TIM*i
REAL(8) :: WTX, WTY, WT !WTX, WTY: Mean square displacement values for X and Y coordinates. WT: "Magnitude" of mean square displacement
REAL(8) :: DIF !Diffusion coefficient values

!Creating files where we'll store W(t) and D(t)
OPEN(28,FILE='WT.dat',STATUS='UNKNOWN')
OPEN(12,FILE='DT.dat',STATUS='UNKNOWN')


TIM = REAL (NFrec2)*DT

!Iteraring over all simulation steps
DO i=1,KI-1

    NTMAX = KI-i

    WTX = 0.d0
    WTY = 0.d0
    WT = 0.d0

    !Iterating over all particles
    DO L=1,N

        !Iterating over simulation steps from j=1 to j=NTMAX=KI-i
        DO j=1, NTMAX
            WTX = WTX + (CXR(L,i+j)-CXR(L,j))**2
            WTY = WTY + (CYR(L,i+j)-CYR(L,j))**2
        END DO

    END DO

    TIME = TIM*Real(i)

    WT = (WTX + WTY)/Real(NTMAX)/Real(N)/6.d0 !Computing W(t)

    DIF = WT/TIME !Computing D(t)

    WRITE(28,*) TIME,WT !Saving t and W(t)
    WRITE(12,*) TIME,DIF !Saving t and D(t)

END DO

CLOSE(28)
CLOSE(12)

END SUBROUTINE
