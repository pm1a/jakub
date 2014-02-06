!> \file geometry.f90
!! \brief typy a procedury pro geometrii
!<

!> module pro popis geometrie
module geometry
    use pmatypy
    implicit none
    !> bod diskretizace
    type, public :: Point
        !> x-ova souradnice
        real(kind=rkind) :: x
        !> y-ova souradnice
        real(kind=rkind) :: y
        !> jestli lezi na hranici
        logical :: onBoundary = .false.
        !> hodnota v bode (realna i imaginarni cast)
        real(kind=rkind),dimension(1:2) :: value = 0
        logical :: used = .false.
    end type Point

    !> trojuhelnik v triangulaci
    type, public :: Triangle
        !> cisla vrcholu
        integer(kind=ikind) :: p1,p2,p3
        !> label druh materialu - cislo
        integer(kind=ikind) :: lab
        !> plocha trojuhelnika
        real(kind=rkind) :: plocha
        !> status rozdeleni
        integer(kind=ikind) :: typ
    end type Triangle

    !> jedna usecka na hranici
    type, public :: Bound
        !> koncove body
        integer(kind=ikind) :: p1,p2
        !> label pro druh okrajove podminky
        integer(kind=ikind) :: lab
        !> status rozdeleni
        integer(kind=ikind) :: typ
    end type Bound

    !> Label
    type, public :: Tlabel
        real(kind=rkind) :: Re = 0, Im = 1
        !> pocet vyskytu
        integer(kind=ikind) :: cnt = 0
    end type Tlabel

    !> label okrajove podminky
    type, public :: Blabel
        !> hodnota okrajove podminky
        real(kind=rkind) :: val = 0
        integer(kind=ikind) :: cnt = 0
    end type Blabel


    !integer, dimension(0:22) :: data1 = (/1,1,2,3,3,3,0,0,3,3,3,&
    !    1,2,4,5,6,7,8,9,10,11,12,13 /)
    integer, dimension(-95:101) :: data2 = 0

    real(kind=rkind), dimension(0:2,1:2) :: bvdata

    ! pocet aktualne uzitych bodu
    integer(kind=ikind) :: LastUsedPoint = 0
    ! konstanty pro pocitani
    real(kind=rkind), parameter :: eps0 = 8.85E-12_rkind
    real(kind=rkind), dimension(0:13) :: epsr=(/ 1.0_rkind, 1.0_rkind, &
        4.49_rkind, 1.0_rkind, 1.0_rkind, 1.0_rkind, 1.0_rkind, 1.0_rkind, &
        1.0_rkind, 1.0_rkind, 1.0_rkind, 1.0_rkind, 1.0_rkind, 1.0_rkind/)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!  tohle je origos
    !real(kind=rkind), dimension(0:13) :: gama=(/ 0.0_rkind, 0.0_rkind, &
        !1.0e-4_rkind, 57.0e6_rkind, 0.0_rkind, 0.0_rkind, 0.0_rkind, 0.0_rkind, &
        !0.0_rkind, 0.0_rkind, 0.0_rkind, 0.0_rkind, 0.0_rkind, 0.0_rkind/)
    ! tohle je modifikace
    real(kind=rkind), dimension(0:13) :: gama=(/ 0.0_rkind, 0.0_rkind, &
        1.0e-4_rkind, 57.0_rkind, 0.0_rkind, 0.0_rkind, 0.0_rkind, 0.0_rkind, &
        0.0_rkind, 0.0_rkind, 0.0_rkind, 0.0_rkind, 0.0_rkind, 0.0_rkind/)

    real(kind=rkind), parameter :: omega = 2*3.14159265_rkind *6.0e6_rkind


contains
end module geometry
