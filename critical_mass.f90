PROGRAM CriticalMass
    implicit none
    real :: massU, R, sigT, xBank(10000000, 2), yBank(10000000, 2), &
            zBank(10000000, 2), nU, ran1, L, theta, phi, x, y, z
    real, parameter :: rhoU = 19.E3, sigS = 4.2E-28, sigA = 0.1E-28, &
                       sigF = 1.3E-28, prob2 = 0.56, prob3 = 0.44
    integer :: n, i, j, numNeutrons, iBank, numBank, current, next

    ! Open the data file
    open(unit = 2, file = 'critical_mass.dat', status = 'replace')

    ! Establish constants

    ! Mass of Uranium-235 in kg
    massU = 1.8E1
    ! Number of starting neutrons
    numNeutrons = 1.0E9
    ! Radius of sphere
    R = (3. * massU / (4. * acos(-1.0) * rhoU)) ** (1. / 3.)
    ! Number density of Uranium-235
    nU = rhoU / (235. * 1.661E-27)
    ! Total cross section
    sigT = sigS + sigA + sigF

    print *, "MassU:", massU, "Radius:", R
    write(2, *) 0, numNeutrons

    ! For 200 generations of neutrons,
    do j = 1, 200
        ! Reset neutron bank index
        iBank = 0

        ! Alternate between the database column indicies
        if(mod(j, 2).eq.0) then
            current = 2
            next = 1
        else
            current = 1
            next = 2
        endif

        ! For every neutron,
        do n = 1, numNeutrons
            if(j.eq.1) then
                ! Emit neutrons on first generation
                call emitNeutron(x, y, z, theta, phi)
            else
                ! Load locations from bank otherwise
                x = xBank(n, current)
                y = yBank(n, current)
                z = zBank(n, current)

                ! Get new directions
                call scatter(theta, phi)
            endif

            ! While neutron is still inside of the sphere,
            do while(sqrt(x*x + y*y + z*z).le.R)
                ! Generate random
                call random_number(ran1)

                ! Find distance traveled
                L = -log(ran1) / (nU * sigT)
                ! And update the neutron positions
                x = x + L * sin(theta) * cos(phi)
                y = y + L * sin(theta) * sin(phi)
                z = z + L * cos(theta)

                ! If still inside of the sphere,
                if(sqrt(x*x + y*y + z*z).le.R) then
                    ! Generate random
                    call random_number(ran1)

                    ! and find interaction type from the probabilities
                    if(ran1.lt.(sigA / sigT)) then
                        ! Neutron absorbed, go to next neutron
                        goto 2
                    elseif(ran1.lt.((sigA + sigS) / sigT)) then
                        ! Neutron scattered
                        call scatter(theta, phi)
                    else
                        ! Fission event, find number of neutrons produced
                        call random_number(ran1)
                        if(ran1.lt.prob2) then
                            ! 2 produced
                            numBank = 2
                        else
                            ! 3 produced
                            numBank = 3
                        endif

                        ! Bank neutrons' positions
                        do i = 1, numBank
                            iBank = iBank + 1

                            ! Exit if the bank is full
                            if (iBank.gt.9999999) goto 3

                            xBank(iBank, next) = x
                            yBank(iBank, next) = y
                            zBank(iBank, next) = z
                        enddo
                    endif
                endif
            enddo

2           continue
        enddo

        print *, "Gen:", j, "Number:", iBank
        ! update the number of neutrons for next generation
        numNeutrons = iBank
        ! Save to file
        write(2, *) j, numNeutrons
    enddo
    ! Close the file
3   close(2)
    print *, "Done."
end PROGRAM CriticalMass

SUBROUTINE emitNeutron(x, y, z, theta, phi)
!   Creates a photon packet by initializing its coordinates and direction vector
!   x       R4 x-coordinate
!   y       R4 y-coordinate
!   z       R4 z-coordinate
!   theta   R4 angle
!   phi     R4 angle
    implicit none
    real, intent(out) :: x, y, z, theta, phi

    ! Start a center of sphere
    x = 0.0
    y = 0.0
    z = 0.0

    call scatter(theta, phi)
end SUBROUTINE emitNeutron

SUBROUTINE scatter(theta, phi)
!   Scatters the photon isotropically by giving it new direction vectors
!   theta   R4 angle
!   phi     R4 angle
    implicit none
    real, intent(out) :: theta, phi
    real :: ran

    call random_number(ran)
    theta = acos(2.0 * ran - 1.0)

    call random_number(ran)
    phi = 2.0 * ran * acos(-1.0)  ! arccos(-1) = pi
end SUBROUTINE scatter
