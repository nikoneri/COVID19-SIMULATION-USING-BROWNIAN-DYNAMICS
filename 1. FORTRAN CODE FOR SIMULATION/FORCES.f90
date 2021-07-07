SUBROUTINE FUERZAS(L,X,Y,FX,FY,BoxL,Rcut,r0,v,ICOV,SUK,IC,N,NN1,DINF)
!*********************************************************************************************************************
!   1. AUTHORS AND E-MAIL CONTACT:
!   -Isaac Neri Gomez Sarmiento (isaacneri.gs@gmail.com)
!   -César Omar Ramírez Álvarez (cesaromarramirezalvarez@gmail.com)
!   -Jonas Valenzuela Terán     (24jonass@gmail.com)
!   -We thank the support and advice of our molecular simulations teacher Dr. Laura Lorenia Yeomans Reyna
!*********************************************************************************************************************
!   2. DESCRIPTION:
!       THIS SUBROUTINE ALLOWS US THREE THINGS:
!           1. COMPUTE FORCES BETWEEN PARTICLES THAT WILL ALLOW US TO MOVE PARTICLES
!           2. COMPUTE SYSTEM'S ENERGY
!           3. COMPUTE THE NUMBER OF INFECTED PEOPLE IN EACH SAVED SIMULATION STEP
!*********************************************************************************************************************

!********************
!VARIABLE DECLARATION
!********************
IMPLICIT NONE
INTEGER,INTENT(IN) :: N !Number of particles
INTEGER,INTENT(IN) :: NN1 !!Maximum  size for particles' properties arrays: NN1=N+10
INTEGER,INTENT(IN) :: L !Simulation step
INTEGER,INTENT(IN) :: SUK !If SUK=1, we'll count infected people, else we won't
INTEGER:: I, J !Iterators
INTEGER:: CONTA, IC !These are counters for the number of new infected people per simulation step
REAL,DIMENSION(NN1), INTENT(IN) :: X, Y !Arrays that hold position coordinates
INTEGER,DIMENSION(NN1) :: ICOV !Array that holds the state of infection of each person
REAL,DIMENSION(NN1) :: FX, FY !Arrays that hold the X and Y components of the net force applied to each particle
REAL, INTENT(IN) :: BOXL !Box size
REAL, INTENT(IN) :: RCUT !Cut radius (forces and energy between particles after that distance are not computed)
REAL, INTENT(IN) :: r0 !Recommended distance to keep from other people.
REAL, INTENT(IN) :: v  !Keep your Distance level of intensity. When v=1, people don't keep their distance. As v increases, people keep their distance more frequently.
REAL(8) :: ENPOT !Total potential energy of the system
REAL(8) :: FXI,FYI !X and Y components of the total force applied to the i-th particle
REAL(8) :: FXIJ,FYIJ !Forces between a pair of particles i and j
REAL(8) :: XIJ,YIJ !Difference between two particles' coordinates
REAL(8) :: RIJ !Distance between two particles
REAL(8) :: U,U2 !U = v/(RIJ**2), U2 = U*((r0/RIJ))**v
REAL :: DINF !Infection distance. If a particle A is infected and particle B is not and if the distance between their center is less than DINF, then B will get infected.

conta=0
ENPOT=0.d0

DO i=1,N
    fx(i)=0.d0
    fy(i)=0.d0
END DO

!Creating new file where we'll store the mean potential energy of each particle
OPEN(UNIT=41,FILE='EP.dat',STATUS='UNKNOWN')


DO i=1,N-1

    FXI=fx(i)
    FYI=fy(i)

    DO J=i+1,N

        XIJ=X(i)-X(j)
        YIJ=Y(i)-Y(j)


        !Using periodic boundary conditions (minimim-image convention)
        XIJ=XIJ-BOXL*DNINT(XIJ/BOXL)
        YIJ=YIJ-BOXL*DNINT(YIJ/BOXL)

        RIJ=SQRT(XIJ**2+YIJ**2)


        IF (SUK == 1) THEN
            IF (RIJ < DINF .AND. ICOV(i)==1 .AND. ICOV(j)==0) THEN
                ICOV(j)=1
                CONTA = CONTA + 1
            ELSE IF (RIJ < DINF .AND. ICOV(i)==0 .AND. ICOV(j)==1) THEN
                ICOV(i)=1
                CONTA = CONTA + 1
            END IF
            IC = CONTA
        END IF

        !Inicia la especificación del modelo de potencial del cual se
        ! deriva la fuerza de interacción entre partículas

        IF (RIJ < RCUT) THEN

            U = v/(RIJ**2)
            U2 = U*((r0/RIJ))**v
            ENPOT = ((r0/RIJ))**v + ENPOT

            FXIJ = (XIJ)*U2
            FYIJ = (YIJ)*U2

            FXI = FXI+FXIJ
            FYI = FYI+FYIJ

            FX(J) = FX(J)-FXIJ
            FY(J) = FY(J)-FYIJ

        ENDIF

    END DO

    FX(I) = FXI
    FY(I) = FYI

END DO

    WRITE(41,*)L,ENPOT/REAL(N)

END SUBROUTINE
