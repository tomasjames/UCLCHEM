!WORK IN PROGRESS
!THIS IS A SERIES OF SUBROUTINES THAT CAN BE COMPILED WITH F2PY TO PRODUCE A PYTHON MODULE
!EACH SUBROUTINE BECOMES A PYTHON FUNCTION

SUBROUTINE uclchem(init_dens,final_dens,shock_vel,phase_flag,outFile,startFile)
    USE physics
    USE chemistry
    IMPLICIT NONE
    double precision, INTENT(IN) :: init_dens,final_dens,shock_vel
    integer, INTENT(IN) :: phase_flag
    character(LEN=*), INTENT(IN) :: outFile,startFile
    character(LEN=100):: abundFile,outputFile,columnFile

    !f2py intent(in) init_dens,final_dens,max_temp,shock_vel,phase_flag,outFile,startFile
    include 'defaultparameters.f90'
    close(10)
    close(11)
    close(7)

    open(10,file=outFile,status='unknown') 
    columnFlag=.False.
    open(7,file=startFile,status='unknown') 

    initialDens=init_dens
    finalDens=final_dens
    vs=shock_vel
    phase=phase_flag

    IF (phase_flag .eq. 1) THEN
        switch=1
        readAbunds=0
        evap=0
        fr=1.0
        desorb=0
        h2desorb=0
        crdesorb=0
        uvcr=0
        IF (init_dens .gt. 1e4) THEN
            baseAv=20
        ELSE
            baseAv=10
        END IF
        rout=(baseAv*(1.6e21)/(init_dens))/(3.086d18)

    ELSE
        switch=0
        readAbunds=1
        evap=1
        fr=0.0
        desorb=1
        h2desorb=1
        crdesorb=1
        uvcr=1
        rout=rout
        baseAv=baseAv
    END IF

    CALL initializePhysics
    CALL initializeChemistry

    dstep=1
    currentTime=0.0
    timeInYears=0.0
    !loop until the end condition of the model is reached 
    DO WHILE ((switch .eq. 1 .and. dens(1) < finalDens) .or. (switch .eq. 0 .and. timeInYears < finalTime))

        !store current time as starting point for each depth step
        IF (points .gt. 1) THEN
            currentTimeold=targetTime
            currentTime=currentTimeold
        END IF
        !Each physics module has a subroutine to set the target time from the current time
        CALL updateTargetTime

        !loop over parcels, counting from centre out to edge of cloud
        DO dstep=1,points

            dens=abund(nspec+1,dstep)
            !update physics
            CALL updatePhysics
            !update chemistry
            CALL updateChemistry
            
            !set time to the final time of integrator rather than target     
            targetTime=currentTime
            !reset target for next depth point
            if (points .gt. 1)currentTime=currentTimeold
            !get time in years for output
            timeInYears= currentTime*year
            !write this depth step
            CALL output
        END DO
    END DO 

END SUBROUTINE uclchem
