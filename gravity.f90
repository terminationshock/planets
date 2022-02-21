! This file is part of Planets.
!
! Planets is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Planets is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Planets.  If not, see <http://www.gnu.org/licenses/>.

module gravity
    implicit none
    real(kind=8), parameter :: pi = 3.141592653589793d0
    real(kind=8), parameter :: G = 8.890125158720449d-10
    real(kind=8), parameter :: drmax = 5.d-4

    contains

    subroutine acceleration(iplanet, nplanets, x, xplanet, m, r, a, collide_any)
        implicit none
        integer(kind=4),                     intent(in)  :: iplanet, nplanets
        real(kind=8), dimension(nplanets,2), intent(in)  :: x
        real(kind=8), dimension(2),          intent(in)  :: xplanet
        real(kind=8), dimension(nplanets),   intent(in)  :: m, r
        real(kind=8), dimension(2),          intent(out) :: a
        logical,                             intent(out) :: collide_any
        real(kind=8), dimension(nplanets,2)              :: acc
        real(kind=8)                                     :: dr
        integer(kind=4)                                  :: i, j

        a(:) = 0.d0
        collide_any = .false.

!$OMP PARALLEL DO PRIVATE(i, j, dr)
        do i = 1, nplanets
            acc(i,:) = 0.d0
            if (i .ne. iplanet .and. m(i) .gt. 0.d0) then
                dr = sqrt((xplanet(1) - x(i,1))**2 + (xplanet(2) - x(i,2))**2)
                if (dr .gt. r(iplanet) + r(i)) then
                    do j = 1, 2
                        acc(i,j) = -G * m(i) * (xplanet(j) - x(i,j)) / dr**3
                    enddo
                else
                    collide_any = .true.
                endif
            endif
        enddo
!$OMP END PARALLEL DO

        a = sum(acc, 1)
    end subroutine acceleration

    subroutine rk4_step(iplanet, nplanets, x, v, m, r, dt, xnew, vnew, collide_any)
        implicit none
        integer(kind=4),                     intent(in)  :: iplanet, nplanets
        real(kind=8), dimension(nplanets,2), intent(in)  :: x, v
        real(kind=8), dimension(nplanets),   intent(in)  :: m, r
        real(kind=8),                        intent(in)  :: dt
        real(kind=8), dimension(2),          intent(out) :: xnew, vnew
        logical,                             intent(out) :: collide_any
        real(kind=8), dimension(4,2)                     :: kx, kv
        real(kind=8), dimension(2)                       :: a, xx, vv
        integer(kind=4)                                  :: j
        logical                                          :: bump

        collide_any = .false.

        do j = 1, 2
            xx(j) = x(iplanet,j)
            vv(j) = v(iplanet,j)
        enddo
        call acceleration(iplanet, nplanets, x, xx, m, r, a, bump)
        if (bump) collide_any = .true.
        do j = 1, 2
            kx(1,j) = dt * vv(j)
            kv(1,j) = dt * a(j)
            xx(j) = x(iplanet,j) + 0.5d0 * kx(1,j)
            vv(j) = v(iplanet,j) + 0.5d0 * kv(1,j)
        enddo
        call acceleration(iplanet, nplanets, x, xx, m, r, a, bump)
        if (bump) collide_any = .true.
        do j = 1, 2
            kx(2,j) = dt * vv(j)
            kv(2,j) = dt * a(j)
            xx(j) = x(iplanet,j) + 0.5d0 * kx(2,j)
            vv(j) = v(iplanet,j) + 0.5d0 * kv(2,j)
        enddo
        call acceleration(iplanet, nplanets, x, xx, m, r, a, bump)
        if (bump) collide_any = .true.
        do j = 1, 2
            kx(3,j) = dt * vv(j)
            kv(3,j) = dt * a(j)
            xx(j) = x(iplanet,j) + kx(3,j)
            vv(j) = v(iplanet,j) + kv(3,j)
        enddo
        call acceleration(iplanet, nplanets, x, xx, m, r, a, bump)
        if (bump) collide_any = .true.
        do j = 1, 2
            kx(4,j) = dt * vv(j)
            kv(4,j) = dt * a(j)
            xnew(j) = x(iplanet,j) + (kx(1,j) + 2.d0*kx(2,j) + 2.d0*kx(3,j) + kx(4,j)) / 6.d0
            vnew(j) = v(iplanet,j) + (kv(1,j) + 2.d0*kv(2,j) + 2.d0*kv(3,j) + kv(4,j)) / 6.d0
        enddo
    end subroutine rk4_step

    subroutine collide(iplanet, nplanets, x, v, m, r)
        implicit none
        integer(kind=4),                     intent(in)    :: iplanet, nplanets
        real(kind=8), dimension(nplanets,2), intent(inout) :: x, v
        real(kind=8), dimension(nplanets),   intent(inout) :: m, r
        real(kind=8)                                       :: dr
        integer(kind=4)                                    :: i, j

        do i = 1, nplanets
            if (i .ne. iplanet .and. m(iplanet) .gt. 0.d0) then
                dr = sqrt((x(iplanet,1) - x(i,1))**2 + (x(iplanet,2) - x(i,2))**2)
                if (dr .le. r(iplanet) + r(i) .and. m(i) .gt. 0.d0 .and. m(iplanet) .ge. m(i)) then
                    do j = 1, 2
                        x(iplanet,j) = 0.5d0 * (x(iplanet,j) + x(i,j))
                        v(iplanet,j) = (m(iplanet) * v(iplanet,j) + m(i) * v(i,j)) / (m(iplanet) + m(i))
                    enddo
                    m(iplanet) = m(iplanet) + m(i)
                    r(iplanet) = (r(iplanet)**3 + r(i)**3)**(1.d0/3.d0)
                    m(i) = 0.d0
                endif
            endif
        enddo

    end subroutine collide

    subroutine integrate(nplanets, ndt, dt, x, v, m, r, dtnew)
        implicit none
        integer(kind=4),                     intent(in)    :: nplanets, ndt
        real(kind=8),                        intent(in)    :: dt
        real(kind=8), dimension(nplanets,2), intent(inout) :: x, v
        real(kind=8), dimension(nplanets),   intent(inout) :: m, r
        real(kind=8),                        intent(out)   :: dtnew
        real(kind=8), dimension(2)                         :: xx, vv
        real(kind=8)                                       :: dr, ddr
        integer(kind=4)                                    :: i, j, n
        logical                                            :: collide_any

        dr = 0.d0
        do n = 1, ndt
            do i = 1, nplanets
                if (m(i) .gt. 0.d0) then
                    call rk4_step(i, nplanets, x, v, m, r, dt, xx, vv, collide_any)
                    ddr = 0.d0
                    do j = 1, 2
                        ddr = ddr + (x(i,j) - xx(j))**2
                        x(i,j) = xx(j)
                        v(i,j) = vv(j)
                    enddo
                    dr = max(dr, sqrt(ddr))
                endif
            enddo

            if (collide_any) then
                do i = 1, nplanets
                    call collide(i, nplanets, x, v, m, r)
                enddo
            endif
        enddo
        if (dr .gt. 0.d0) then
            dtnew = drmax / dr * dt
        else
            dtnew = 1.d0
        endif
    end subroutine integrate

    subroutine orbit(nplanets, x, v, m, a, b, e, alpha, perihelion, aphelion, T, isun)
        implicit none
        integer(kind=4),                        intent(in)  :: nplanets
        real(kind=8),    dimension(nplanets,2), intent(in)  :: x, v
        real(kind=8),    dimension(nplanets),   intent(in)  :: m
        real(kind=8),    dimension(nplanets),   intent(out) :: a, b, e, alpha, perihelion, aphelion, T
        integer(kind=4), dimension(nplanets),   intent(out) :: isun
        real(kind=8)                                        :: r, vv, h, p, mu, cosalpha, energy, energy_min
        integer(kind=4)                                     :: i, j

        a(:) = 0.d0
        b(:) = 0.d0
        e(:) = -1.d0
        alpha(:) = 0.d0
        perihelion(:) = 0.d0
        aphelion(:) = 0.d0
        T(:) = 0.d0
        isun(:) = 0

!$OMP PARALLEL DO PRIVATE(i, j, vv, mu, r, energy, energy_min)
        do i = 1, nplanets
            if (m(i) .gt. 0.d0) then
                vv = sqrt(v(i,1)**2 + v(i,2)**2)
                energy_min = 0.d0
                do j = 1, nplanets
                    if (m(i) .lt. m(j)) then
                        mu = G * (m(i) + m(j))
                        r = sqrt((x(i,1)-x(j,1))**2 + (x(i,2)-x(j,2))**2)
                        energy = 0.5d0 * vv**2 - mu / r
                        if (energy .lt. energy_min) then
                            isun(i) = j
                            energy_min = energy
                        endif
                    endif
                enddo
            endif
        enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(i, r, vv, mu, h, p, cosalpha)
        do i = 1, nplanets
            if (isun(i) .gt. 0) then
                r = sqrt((x(i,1)-x(isun(i),1))**2 + (x(i,2)-x(isun(i),2))**2)
                vv = sqrt(v(i,1)**2 + v(i,2)**2)
                mu = G * (m(i) + m(isun(i)))
                h = (x(i,1)-x(isun(i),1)) * v(i,2) - (x(i,2)-x(isun(i),2)) * v(i,1)
                p = h**2 / mu

                if (mu .eq. 0.5d0 * r * vv**2) cycle
                a(i) = 1.d0 / (2.d0 / r - vv**2 / mu)

                if (p .gt. a(i)) cycle
                e(i) = sqrt(1.d0 - p / a(i))
                b(i) = a(i) * sqrt(1.d0 - e(i)**2)
                perihelion(i) = a(i) * (1.d0 - e(i))
                aphelion(i) = a(i) * (1.d0 + e(i))
                T(i) = 2.d0 * pi * sqrt(a(i)**3 / (G * (m(isun(i)) + m(i))))

                cosalpha = (r - a(i) * (1.d0 - e(i)**2)) / (e(i) * r)
                if (abs(cosalpha) .gt. 1.d0) cycle
                alpha(i) = sign(acos(cosalpha), (x(i,1)-x(isun(i),1))*v(i,1) + (x(i,2)-x(isun(i),2))*v(i,2))
            endif
        enddo
!$OMP END PARALLEL DO
    end subroutine orbit

end module gravity
