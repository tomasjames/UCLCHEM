!Simple physics module. Models points along a 1d line from the centre to edge of a cloud. Assuming the cloud is spherical you can average
!over the points to get a 1d average and then assume the rest of sphere is the same.

MODULE physics
    IMPLICIT NONE
    !Use main loop counters in calculations so they're kept here
    integer :: dstep,points
    !Switches for processes are also here, 1 is on/0 is off.
    integer :: collapse,switch,first,phase
    integer :: h2desorb,crdesorb,crdesorb2,uvcr,desorb

    !evap changes evaporation mode (see chem_evaporate), ion sets c/cx ratio (see initializeChemistry)
    !Flags let physics module control when evap takes place.flag=0/1/2 corresponding to not yet/evaporate/done
    integer :: evap,ion,solidflag,monoflag,volcflag,coflag,tempindx

    !variables either controlled by physics or that user may wish to change    
    double precision :: initialDens,timeInYears,initialTime,targetTime,currentTime,currentTimeold,finalDens,finalTime,d,dMin,dMax
    double precision :: cloudSize,rout,rin,baseAv,bc,tstart,maxTemp,dCool,lambda
    double precision, allocatable :: av(:),coldens(:),temp(:),dens(:),pressure(:)

    !Everything should be in cgs units. Helpful constants and conversions below
    double precision,parameter ::pi=3.141592654,mh=1.67e-24,kbolt=1.38d-23,mu=1
    double precision, parameter :: year=3.16455d-08,pc=3.086d18,km=1.d5,seconds=3.154d7

    character(2) ::filename
    character(1)  ::densint

    !Cshock specific parameters
    !*******************************************************************
    double precision :: initialTemp, z2,vs,v0,zn,vn,at,z3,targetTime0,tsat
    double precision :: ucm,z1,dv,vi,tempi,vn0,zn0,vA,dlength
    double precision :: grainRadius5,grainRadius,dens6
    double precision, allocatable :: tn(:),ti(:),tgc(:),tgr(:),tg(:)
    !variables for the collisional and radiative heating of grains
    double precision :: mun,tgc0,Frs,tgr0,tgr1,tgr2,tau100,trs0,G0
    double precision :: coshinv1,coshinv2,zmax,a1

    integer :: inrad
    double precision, parameter::nu0=3.0d15,kb2=1.38d-16,bm0=1.e-6,bt=6.
    !*******************************************************************

CONTAINS
!THIS IS WHERE THE REQUIRED PHYSICS ELEMENTS BEGIN. YOU CAN CHANGE THEM TO REFLECT YOUR PHYSICS BUT THEY MUST BE NAMED ACCORDINGLY.
    
    SUBROUTINE initializePhysics
    !Any initialisation logic steps go here
        allocate(av(points),coldens(points),temp(points),dens(points))
        cloudSize=(rout-rin)*pc

        if (collapse .eq. 1) THEN
            write(*,*) "Initialising physics."
            if (phase .eq. 2) THEN
                write(*,*) "Cannot have collapse on during jshock (phase=2)"
                Write(*,*) "setting collapse=0 and continuing"
                collapse=0
                dens=initialDens
            ELSE
                dens=1.001*initialDens
            END IF
        ELSE
            dens=initialDens
        ENDIF 

        write(*,*) "Setting temperature and cooling time."
        temp = initialTemp
        ! tstepCalc = real(vs*1d5)*real(points)/real(size)

        !calculate initial column density as distance from core edge to current point * density
        DO dstep=1,points
            coldens(dstep)=real(points-dstep+1)*cloudSize/real(points)*initialDens
            write(*,*) "Determined column density as ", coldens(dstep)
        END DO

        IF (phase .eq. 2) THEN
            allocate(tn(points),ti(points),tgc(points),tgr(points),tg(points))
            mun=2*mh
            grainRadius5=grainRadius/4.e-5
            dens6=dens(dstep)/1.e6
            targetTime0=0
        END IF

        !maxxtemp set by vs and pre-shock density, polynomial fits to values taken from Draine et al. 1983
        !have been made and coefficients placed here. Tested with log(dens)>3 <6
        IF (initialDens .gt. 10**4.5) THEN
            maxTemp=(2.91731*vs*vs)-(23.78974*vs)+225.204167337
        ELSE
            maxTemp=(0.47258*vs*vs)+(40.44161*vs)-128.635455216
        END IF    
        temp=initialTemp

        !tsat proportional to 1/pre-shock density. Fit to tsats from Jimenez-Serra 2008.
        tsat=(-15.38729*vs*vs*vs)+(2069.56962*vs*vs)-(90272.826991*vs)+1686858.54278
        tsat=tsat/initialDens

        write(*,*)tsat,maxTemp

    END SUBROUTINE initializePhysics

    !This is the time step for outputs from UCL_CHEM NOT the timestep for the integrater. DLSODE sorts that out based on chosen error
!tolerances (RTOL/ATOL) and is simply called repeatedly until it outputs a time >= targetTime. targetTime in seconds for DLSODE, timeInYears in
!years for output.

    SUBROUTINE updateTargetTime
        IF (phase .eq. 1) THEN
            IF (timeInYears .gt. 1.0d6) THEN
                targetTime=(timeInYears+1.0d5)/year
            ELSE IF (timeInYears .gt. 10000) THEN
                targetTime=(timeInYears+1000.0)/year
            ELSE IF (timeInYears .gt. 1000) THEN
                targetTime=(timeInYears+100.0)/year
            ELSE IF (timeInYears .gt. 0.0) THEN
                targetTime=(timeInYears*10)/year
            ELSE
                targetTime=3.16d7*10.d-8
            ENDIF
        ELSE
            IF (timeInYears .gt. 1.0d5) THEN
                targetTime=(timeInYears+100)/year
            ELSE IF (timeInYears.gt. 1.0d4) THEN
                targetTime=(timeInYears+10.)/year
            ELSE IF (timeInYears .gt. 1000) THEN
                targetTime=(timeInYears+1.)/year
            ELSE IF (timeInYears .gt. 10) THEN
                targetTime=(timeInYears+0.1)/year
            ELSE IF (timeInYears .gt. 1) THEN
                targetTime=(timeInYears+0.01)/year
            ELSE IF (timeInYears .gt. 0.05) THEN
                targetTime=(timeInYears+0.001)/year
            ELSE IF  (timeInYears.gt.0.0) THEN
                targetTime=(timeInYears+0.0001)/year
            ELSE
                targetTime=3.16d1
            ENDIF
        END IF
    END SUBROUTINE updateTargetTime

    !This routine is formed for every parcel at every time step.
    !update any physics here. For example, set density
    SUBROUTINE updatePhysics
        !calculate column density. Remember dstep counts from edge of core in to centre
        IF (dstep .lt. points) THEN
            !column density of current point + column density of all points further out
            coldens(dstep)=(cloudSize/real(points))*dens(dstep)
            coldens(dstep)=coldens(dstep)+sum(coldens(dstep:points))
        ELSE
            coldens(dstep)=cloudSize/real(points)*dens(dstep)
        END IF
                !calculate the Av using an assumed extinction outside of core (baseAv), depth of point and density
        av(dstep)= baseAv +coldens(dstep)/1.6d21

        ! phase=2 for the J-shock specific calculations
        IF (phase .eq. 2) THEN
            ! write(*,*) "In phase 2."
            ! Determine the cooling length (of the order of the mean free path)
            dCool = 1/((2**0.5)*initialDens*(pi*(2*5.291d-09)**2))
            ! Determine initial distance
            dMin = initialTime*vs*1d5
            ! Determine the final distance
            dMax = finalTime*vs*1d5
            ! Determine the current distance
            d = currentTime*vs*1d5
            lambda = LOG(maxTemp/initialTemp)

            ! Determine whether shock is still increasing the temperature
            ! Or whether it is in the post-shock cooling phase
            ! Or whether the temperature is now constant
            ! write(*,*) "Beginning the temperature and density calculations."
            IF (d - dMin .lt. dCool) THEN
                ! write(*,*) "The condition d - initialTime*vs .lt. dCool has been met."
                ! If it is in the initial temperature rise to max temp, then
                tn(dstep) = ((d/dCool)**3)*(maxTemp - initialTemp) + initialTemp
                dens = ((d/dCool)**2)*(4*initialDens) + initialDens
                ! write(*,*) "d: ", d, ' cm'
                ! write(*,*) "d/dCool: ", d/dCool
                ! write(*,*) "(d/dCool)**2: ", (d/dCool)**2
                ! write(*,*) "((d/dCool)**2)*(4*initialDens): ", ((d/dCool)**2)*(4*initialDens)
            ELSE IF (d - dMin .gt. dCool .AND. d - dMin .lt. dCool*1d6) THEN
                ! Otherwise we're in the cooling phase
                ! Determine the decay constant for the cooling
                tn(dstep) = maxTemp*EXP(-lambda*(d/(dCool*1d6)))
                dens = 4*initialDens
            ELSE 
                tn(dstep) = initialTemp
                dens = 4*initialDens
            END IF
            ! write(*,*) "Temperature, tn: ", tn(dstep)
            temp(dstep)=tn(dstep)

            IF (timeInYears .gt. 0) THEN
                write(92,1234) tn(dstep),dens(dstep),timeInYears
                1234 format(3(1x,ES12.4e2))
            ENDIF

            !At tsat, all mantle species evaporated. These flags make chem module aware of it.
            IF (timeInYears .gt. tsat .and. coflag .eq. 0) THEN
                evap=2
                coflag=1
            ENDIF
        ENDIF

    END SUBROUTINE updatePhysics

    !This FUNCTION works out the time derivative of the density, allowing DLSODE to update density with the rest of our ODEs
    !It get's called by F, the SUBROUTINE in chem.f90 that sets up the ODEs for DLSODE
    !Currently set to Rawlings 1992 freefall.
    pure FUNCTION densdot()
        double precision :: densdot

    !Rawlings et al. 1992 freefall collapse. With factor bc for B-field etc
        IF (dens(dstep) .lt. finalDens) THEN
             densdot=bc*(dens(dstep)**4./initialDens)**0.33*&
             &(8.4d-30*initialDens*((dens(dstep)/initialDens)**0.33-1.))**0.5
        ELSE
            densdot=1.0d-30       
        ENDIF
    END FUNCTION densdot

END MODULE physics

!REQUIRED PHYSICS ENDS HERE, ANY ADDITIONAL PHYSICS CAN BE ADDED BELOW.
