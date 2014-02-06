!> \file typy.f90
!! \brief definice kind konstant
!<


!> module pro definice typovych konstant
module typy
   !> integery - pro citace iteraci a indexovani
   integer, parameter,public :: ikind = SELECTED_INT_KIND(10)
   !> pro realna cisla
   integer, parameter, public :: rkind = SELECTED_REAL_KIND(25,99)
   !integer, parameter, public :: rkind = SELECTED_REAL_KIND(15,99)
   !>tohle je pro plplot - radeji neuzivat
   integer, parameter, public :: r8 = 8
end module typy
