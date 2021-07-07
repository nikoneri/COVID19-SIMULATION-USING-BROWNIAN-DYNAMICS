PROGRAM COVID19_BROWNIAN_DYNAMICS
!*********************************************************************************************************************
!   1. AUTHORS AND E-MAIL CONTACT:
!   -Isaac Neri Gomez Sarmiento (isaacneri.gs@gmail.com)
!   -César Omar Ramírez Álvarez (cesaromarramirezalvarez@gmail.com)
!   -Jonas Valenzuela Terán     (24jonass@gmail.com)
!   -We thank the support and advice of our molecular simulations teacher Dr. Laura Lorenia Yeomans Reyna
!*********************************************************************************************************************
!   2. DESCRIPTION:
!   *In this program we obtained a simple contagion model, when making certain analogies between a human population and
!   a two-dimensional coloidal suspension system under brownian motion.
!   We consider Hermosillo's population density (a city located in Sonora, Mexico).
!   In this model we include sanitary measurements such as "Quedate en Casa" or Stay Home
!   and "Sana Distancia" or Keep your Distance.
!
!   *We can compare different scenarios such as :
!   1. Not implementing any sanitary measurements (the parameters that allow this scenario are set by default in this code).
!   2. Implementing just Stay Home (with a gradual intensity).
!   3. Implementing just Keep your Distance (with a gradual intensity).
!   4. Implementing both Stay Home and Keep your Distance (with a gradual intensity).
!
!****************************************************************************************************************************************
!   3. SOME REMARKS:
!   *We use the Ermak algorithm to move particles (people) and solve overdamped Langevin equations or in the diffusive regime.
!   *This work is originally in Spanish, so some of the variables and subroutines are in that language.
!   *This simple contagion model doesn't include recovered people or inmmune people. Including this in the future will improve the model.
!   *We consider people to be 2D colloids with a diameter of 1 m.
!   *We work with reduced units, also called dimensionless units.
!   *We work with a repulsive potential of the type U=(r0/r)^v.
!   *Besides computing the number of infected people, this program can also compute structural and dynamical properties:
!       1. Radial distribution function g(r), which gives us an insight on the most likely average distance between particles.
!       2. Mean Square Displacement W(t), which gives us an insight on how much area the particles spread out.
!       3. Diffussion Coefficient D(t), which gives us an insight on the mobility (speed) of particles.
!   *This work has been presented at two congress:
!       1. LXIII National Congress of Physics in Mexico.
!       2. XXXI National Week of Research and Teaching in Mathematics in Hermosillo, Sonora, Mexico.
!****************************************************************************************************************************************

!***********************
!0. VARIABLE DECLARATION
!***********************
IMPLICIT NONE !Only the variables that we declare will work.
INTEGER :: N !Total number of people (this number doesn't have to be the real number of people in a real population, what matters is that the population density matches)
REAL ::  fi !Initial fraction of infected people.
REAL :: NI !Initial number of infected people: NI=N*fi
REAL :: Dens ! Population density (reduced units): (Real Number of People)*(diameter of people in m)^2/(Real Area in m^2)
REAL :: r0 !Recommended distance to keep from other people.
REAL :: QC !Stay Home level of intensity. When QC=1, people don't stay home. As QC increases, people have less mobility.
REAL :: v !Keep your Distance level of intensity. When v=1, people don't keep their distance. As v increases, people keep their distance more frequently.
REAL :: DINF !Infection distance. If a particle A is infected and particle B is not and if the distance between their center is less than DINF, then B will get infected.
REAL :: DT !Time step
REAL :: Var! Varianza: Var=SQRT(2.D0*DT)
REAL :: BoxL !Size of the simulation square
REAL :: RCut !Cut radius for particle-particle interactions
REAL :: PI !Value of pi
INTEGER :: NN1 !Maximum row size for particles' properties matrices: NN1=N+10
INTEGER :: NN2 !Number of saved configurations or simulation steps
INTEGER :: NN3 !Number of histogram bins to compute g(r).
INTEGER :: Nener !Step or configuration where we consider the system has reached thermalization (system energy is almost constant)
INTEGER :: NStep !Number of simulation steps or configurations
INTEGER :: NFREC !Frecuency at which we'll be saving configurations for computing g(r) (i.e every 200 configs.)
INTEGER :: NFrec2 !Frecuency at which we'll be saving configurations for computing W(t) and D(t) (i.e same as NFREC)

INTEGER :: Opcion !Variable to choose between an initial regular configuration or a random initial configuration.
INTEGER :: i,j,L !Iterators for for loops (i is for particles, j is for saved simulation steps, L is for every simulation step, either saved or not saved).
INTEGER :: KI !Counter for the number of saved simulation steps for position matrices with periodic boundary conditions
INTEGER :: KI2 !Counter for the number of saved simulation steps for the position matrices that will be used for computing g(r), W(t) and D(t)
INTEGER :: KI3 !Counter for the number of saved simulation steps for the matrix of infected people.
REAL :: RN1, RN2 !Variables that will hold uniform random numbers
REAL :: AX, AY !Variables that will hold gaussian random numbers

REAL,ALLOCATABLE,DIMENSION(:) :: X,Y !One dimensional arrays that will hold X, Y coordinates of every particle at each simulation step (periodic boundary conditions)
REAL,ALLOCATABLE,DIMENSION(:,:) :: CX,CY !Two dimensional arrays. Rows=particles, Columns:position at saved steps (periodic boundary conditions)
REAL,ALLOCATABLE,DIMENSION(:) :: XR,YR !One dimensional arrays that will hold X, Y coordinates of every particle at each simulation step (no periodic boundary conditions)
REAL,ALLOCATABLE,DIMENSION(:,:) :: CXR,CYR !Two dimensional arrays. Rows: particles, Columns:position at saved steps (no periodic boundary conditions)
REAL,ALLOCATABLE,DIMENSION(:) :: FX,FY !One dimensional arrays that will hold X, Y force components applied to each particle at each simulation step


INTEGER,ALLOCATABLE,DIMENSION(:):: ICOV !One dimensional array that will hold the state of infection for every person (1 is infected and 0 is not infected).
INTEGER,ALLOCATABLE,DIMENSION(:,:):: ICOVi !Two dimensional array. Rows: particles, Columns: State of infection in every saved configuration.
INTEGER :: IT !Number of accumulated infected people
INTEGER :: IC !Number of infected people in each configuration

!NOTE: CX, CY and ICOVi matrices will be used for creating the animation GIF in Python.

!*********************
!1. INITIAL CONDITIONS
!*********************

N=529 !We choose this number of people because at first we were using a initial regular configuration in which SQRT(529)=23, so 23 people would be at each side of the simulation square.
dens=0.004829D0  !This population density is for Hermosillo. It was computed as dens=(812229 people)*(1 m^2)/(168.2*10**6 m^2)
fi=0.1!10% of N will get infected

!The following parameters simulate the case where people don't keep their distance, neither they stay home.
r0=1.D0 !People keep a 1 meter distance from each other's center, which means they are touching, since they disks of 1 meter of diameter
v=1.D0 !v=1 means people are not keeping their distance. Increase this value so people keep their distance more frequently.
QC=1.D0 !QC=1 means people are not staying home. Increase this value if you want people to reduce their mobility.
DINF=2.5 !Maximum distance of infection. If distance between an infected person and a non-infected person is less than DINF, then the non-infected person will get infected. Otherwise, it won't.

Nener=15000
NFREC=200
NFrec2 = 200
Nstep=215000
NN1 = N+10
NN2=int((NStep - Nener)/NFrec2)
NN3 =1000000

!*********************************
!2. DEFINING THE DIMENSION FOR ARRAYS
!*********************************

ALLOCATE (X(NN1),Y(NN1))
ALLOCATE (XR(NN1),YR(NN1))
ALLOCATE (FX(NN1),FY(NN1))
ALLOCATE (CX(NN1,NN2),CY(NN1,NN2))
ALLOCATE (CXR(NN1,NN2),CYR(NN1,NN2))
ALLOCATE (ICOV(NN1))
ALLOCATE (ICOVi(NN1,NN2))

!************************************************************
!3. INITIALIZING ARRAYS, VARIABLES AND CREATING FILES WE'LL NEED
!************************************************************
PI = 4.D0*DATAN(1.D0) !Value for pi=3.14159...
DT = 0.0004D0
BoxL = ((1.D0*N)/Dens)**(1./2.)
RCut = BoxL/2.D0
Var = SQRT(2.D0*DT)
KI = 0
KI2 = 0
KI3=0
NI=INT(REAL(N)*fi)
IC=0
IT=0
write(*,*) "The initial number of infected people will be: ", NI

DO i = 1, NN1
    ICOV(i) = 0 !Initializing with 0's the array of infected people
END DO

DO i = 1, NN1
    DO j=1,NN2
        ICOVi(i,j) = 0 !Initializing with 0's the matrix of infected people
    END DO
END DO

OPEN (12, FILE = "CFin.dat", STATUS = 'unknown')
OPEN (1, FILE = "X.dat", STATUS = 'unknown')
OPEN (2, FILE = "Y.dat", STATUS = 'unknown')
OPEN (3, FILE = "Infectadas.dat", STATUS = 'unknown')
OPEN (50, FILE = "MATRIZ_INFECCION.dat", STATUS = 'unknown')

!***************************************************************
!4. INITIALIZING PEOPLE'S POSITION IN THE SIMULATION BOX AND FORCES
!***************************************************************
Opcion=1 !1 for random initial configuration. 0 for regular initial configuration.

IF (Opcion == 0) THEN
    PRINT*, "Regular initial configuration was chosen"
    CALL CONFIGINIREG(X,Y,N,BoxL)
ELSE
    PRINT*, "Random initial configuration was chosen"
    CALL CONFIGINIALE(X,Y,N,BoxL)
END IF


!The 0 in the subroutine FUERZAS indicates that we are not counting infected people yet
!If it was 1, then we'll count the number of infected people.
!Since we compute distances between people to get forces and potential energy,
!this distance will allow us to know if two persons approach each other enough to get infected.
L=1
CALL FUERZAS(L,X,Y,FX,FY,BoxL,RCut,r0,v,ICOV,0,IC,N,NN1,DINF)



!*************************************************************************************
!5. MOVING PARTICLES (PEOPLE) UNTIL THERMALIZATION (EQUILIBRATION OF POTENTIAL ENERGY)
!THERE IS NOT INFECTED PEOPLE YET
!*************************************************************************************

!ITERATING OVER CONFIGURATIONS
 DO L = 1, INT(nener)-1
    IF (MOD (L, 1000) == 0 ) THEN
        print*, "Configuration", L
    END IF

    CALL initRandomSeed() !Here we're reinitializing the seed of the random number generator

!ITERATING OVER EACH PARTICLE
    DO i = 1, N

!HERE WE'RE GENERATING RANDOM GAUSSIAN NUMBERS AX AND AY THAT WILL BE USED TO MOVE PARTICLES
      10 CALL RANDOM_NUMBER(RN1)
       CALL RANDOM_NUMBER(RN2)
       IF (RN1 <= 0.0001) THEN
            GO TO 10
       END IF
       AX=SQRT(-2.D0*LOG(RN1))*COS(2.D0*PI*RN2)

      11 CALL RANDOM_NUMBER(RN1)
       CALL RANDOM_NUMBER(RN2)
       IF (RN1 <= 0.0001) THEN
            GO TO 11
       END IF
       AY=SQRT(-2.D0*LOG(RN1))*COS(2.D0*PI*RN2)


!ERMAK ALGORITHM TO MOVE PARTICLES

         X(i) = X(i) + FX(i)*DT + Var*AX
         Y(i) = Y(i) + FY(i)*DT + Var*AY

         XR(i) = XR(i) + FX(i)*DT + Var*AX
         YR(i) = YR(i) + FY(i)*DT + Var*AY

!APPLYING PERIODIC BOUNDARY CONDITIONS (SO PARTICLES DON'T GET OUT OF THE BOX: IF ONE DISSAPPEARS FROM THE RIGHT, IT WILL APPEAR ON THE LEFT)
         X(i) = X(i) - BoxL*ANINT(X(i)/BoxL)
         Y(i) = Y(i) - BoxL*ANINT(Y(i)/BoxL)

    END DO


!DECIDING IF SAVING POSITIONS AT CERTAIN CONFIGURATION (PERIODIC BOUNDARY CONDITIONS)
    IF (MOD (L, NFrec) == 0 .AND. L > Nener) THEN
        KI = KI + 1
        DO i = 1,N,1
            CX(i,KI) = X(i)
            CY(i,KI) = Y(i)
        END DO
    END IF

!DECIDING IF SAVING POSITIONS AT CERTAIN CONFIGURATION (NO PERIODIC BOUNDARY CONDITIONS)
    IF (MOD (L, NFrec2) == 0 .AND. L > Nener) THEN
        KI2 = KI2 + 1
        DO i = 1,N,1
            CXR(I,KI2)=XR(I)
            CYR(I,KI2)=YR(I)
        END DO
    END IF

!COMPUTING NEW FORCES
    CALL FUERZAS(L,X,Y,FX,FY,BoxL,RCut,r0,v,ICOV,0,IC,N,NN1,DINF)

END DO

!***************************************
!6.INFECTING PEOPLE AFTER THERMALIZATION
!***************************************

IC=NI
IT=IC
WRITE (*,*) "THE FIRST", NI, "PEOPLE ARE NOW INFECTED"

DO I=1, INT(NI)
    ICOV(i)=1
END DO

!************************************************************
!7. MOVING PARTICLES (PEOPLE), BOTH INFECTED AND NON-INFECTED.
!************************************************************

!ITERATING OVER EACH CONFIGURATION, AFTER THERMALIZATION
DO L=NENER,NSTEP
    IF (MOD (L, 1000) == 0 ) THEN
        print*, "Configuration", L
    END IF
     CALL initRandomSeed() !Here we're reinitializing the seed of the random number generator

!ITERATING OVER EACH PERSON
    DO I=1,N

!HERE WE'RE GENERATING RANDOM GAUSSIAN NUMBERS AX AND AY THAT WILL BE USED TO MOVE PARTICLES
      28 CALL RANDOM_NUMBER(RN1)
       CALL RANDOM_NUMBER(RN2)
       IF (RN1 <= 0.0001) THEN
            GO TO 28
       END IF
       AX=SQRT(-2.D0*LOG(RN1))*COS(2.D0*PI*RN2)


      29 CALL RANDOM_NUMBER(RN1)
       CALL RANDOM_NUMBER(RN2)
       IF (RN1 <= 0.0001) THEN
            GO TO 29
       END IF
       AY=SQRT(-2.D0*LOG(RN1))*COS(2.D0*PI*RN2)


!ERMAK ALGORITHM TO MOVE PARTICLES. NOW IT INCLUDES THE STAY HOME PAREMETER "QC".

         X(i) = X(i) + FX(i)*DT/QC + Var*AX/SQRT(QC)
         Y(i) = Y(i) + FY(i)*DT/QC + Var*AY/SQRT(QC)

         XR(i) = XR(i) + FX(i)*DT + Var*AX
         YR(i) = YR(i) + FY(i)*DT + Var*AY

!APPLYING PERIODIC BOUNDARY CONDITIONS (SO PARTICLES DON'T GET OUT OF THE BOX: IF ONE DISSAPPEARS FROM THE RIGHT, IT WILL APPEAR ON THE LEFT)
         X(i) = X(i) - BoxL*ANINT(X(i)/BoxL)
         Y(i) = Y(i) - BoxL*ANINT(Y(i)/BoxL)

    END DO

!DECIDING IF SAVING ARRAY OF INFECTED PEOPLE AT CERTAIN CONFIGURATION
    IF (MOD (L, NFrec) == 0) THEN
        KI3 = KI3 + 1
        DO i=1,N,1
            ICOVi(I,KI3) = ICOV(I)
        END DO
    END IF

!COMPUTING NEW FORCES AND TRACKING NEW INFECTIONS INSIDE THE SAME SUBROUTINE "FORCES"
    CALL FUERZAS(L,X,Y,FX,FY,BoxL,RCut,r0,v,ICOV,1,IC,N,NN1,DINF)


    WRITE(3,*)L,IT !WRITING IN SCREEN NUMBER OF ACUMULATED INFECTIONS PER CONFIGURATION
    IT=IT+IC !Adding IC newn infecttions to the acumulated infections IT

!DECIDING IF SAVING POSITIONS AT CERTAIN CONFIGURATION (PERIODIC BOUNDARY CONDITIONS)
     IF (MOD (L, NFrec) == 0 .AND. L > Nener) THEN
        KI = KI + 1
        DO i = 1,N,1
            CX(i,KI) = X(i)
            CY(i,KI) = Y(i)
        END DO
    END IF

!DECIDING IF SAVING POSITIONS AT CERTAIN CONFIGURATION (NO PERIODIC BOUNDARY CONDITIONS). THESE POSITIONS WILL BE USED TO COMPUTE g(r), W(t) and D(t).
    IF (MOD (L, NFrec2) == 0 .AND. L > Nener) THEN
        KI2 = KI2 + 1
        DO i = 1,N,1
            CXR(I,KI2)=XR(I)
            CYR(I,KI2)=YR(I)
        END DO
    END IF
END DO


!***************************************
!8.SAVING POSITIONS AND INFECTION MATRIX
!***************************************
DO i = 1, N, 1
    WRITE(12,*)X(I),Y(I) !THIS IS FOR THE FINAL CONFIGURATION
    WRITE(1,*) (SNGL(CX(I,L)),L=1,NN2)
    WRITE(2,*) (SNGL(CY(I,L)),L=1,NN2)
    WRITE(50,*) (ICOVi(I,L),L=1,NN2)
END DO

WRITE(*,*)"TOTAL INFECTED PEOPLE (INCLUDING THE INITIAL ONES)", IT
WRITE(*,*)"COMPUTING RADIAL DISTRIBUTION FUNCTION G(R)..."
CALL GDR (N,NN1,KI,NN3,Dens,BoxL,CX,CY)

PRINT*, "COMPUTING MEAN SQUARE DISPLACEMENT W(T) AND DIFUSSION COEFFICIENT FUNCTIONS D(T)"
CALL WDT(CXR,CYR,KI,KI2,DT,NFREC2,NN1,NN2,N)


WRITE(*,*)"ALL IS DONE!"


!THIS IS TO RELEASE MEMORY THAT ARRAYS USED
DEALLOCATE(X, Y)
DEALLOCATE(XR,YR)
DEALLOCATE(CX, CY)
DEALLOCATE(CXR,CYR)
DEALLOCATE(FX,FY)
DEALLOCATE(ICOVi)
DEALLOCATE(ICOV)

END PROGRAM COVID19_BROWNIAN_DYNAMICS
