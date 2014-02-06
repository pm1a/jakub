!> \file kreslitko.f90
!! \brief nastroje pro kresleni
!!
!<



!> \brief nastroje pro kresleni
!!
module kreslitko

    public :: plotstart
    public :: plotstop
    public :: plotriang
    public :: plotbound
    public :: plotboundi
    public :: plotdomain
    public :: plotdomain1
    public :: meze
    logical, private :: IsInit = .false.
contains

    !> \brief inicializace kresleni
    !!
    !!
    subroutine plotstart()
        implicit none
        print *,"nastavuji"
    end subroutine plotstart

    !> \brief ukonci kresleni
    !!
    !!
    subroutine plotstop()
        print *, "koncim"
    end subroutine plotstop

    !> \brief vykresli triangulaci
    !!
    !! \param U data ulohy
    !!
    subroutine plotriang(U)
        use zobraz
        use geometry
        use pmatypy
        use global
        implicit none

        type(Uloha), intent(in) :: U
        type(Picture_Type) :: P
        real(kind=rpl) :: xmin, xmax,ymin,ymax,sx,sy
        real(kind=rpl) :: wx,wy
        real(kind=rpl), dimension(1:4) :: tx, ty
        type(Point2d), dimension(1:4) :: Path
        integer(kind=ikind) :: i,j

        ! napred najdeme oramovani
        call meze(U,xmin,xmax,ymin,ymax)
        sx = (xmax-xmin)/100
        sy = (ymax-ymin)/100
        xmin = xmin - sx
        xmax = xmax + sx
        ymin = ymin - sy
        ymax = ymax + sy
        sx = 100
        sy = sx*(ymax-ymin)/(xmax-xmin)
        print *, sx,sy,xmin,xmax,ymin,ymax
        call Create(P,sx,sy,xmin,xmax,ymin,ymax)
        do i=0, U%nt-1
            Path(1)%x = U%PointList(U%TriList(i)%p1)%x
            Path(1)%y = U%PointList(U%TriList(i)%p1)%y
            Path(2)%x = U%PointList(U%TriList(i)%p2)%x
            Path(2)%y = U%PointList(U%TriList(i)%p2)%y
            Path(3)%x = U%PointList(U%TriList(i)%p3)%x
            Path(3)%y = U%PointList(U%TriList(i)%p3)%y
            Path(4)%x = Path(1)%x
            Path(4)%y = Path(1)%y
            call Line2d(P,Path,Thickness = 0.00001_rpl)
        end do
        call Close(P)
    end subroutine plotriang


    !> \brief vykresli hranici
    !!
    !! \param U uloha
    !!
    subroutine plotbound(U)
        use pmatypy
        use global
        use geometry
        use zobraz
        use pmatools
        implicit none
        type(Uloha), intent(in) :: U
        real(kind=rpl) :: xmin, xmax,ymin,ymax
        integer :: wl
        integer(kind=ikind) :: i
        type(Point2d), dimension(1:2) :: Path
        type(Picture_Type) :: P

        call Pockej("Kreslim hranici.")
        call meze(U,xmin,xmax,ymin,ymax)
        call Create(P,100.0_rpl,100*(ymax-ymin)/(xmax-xmin),&
        xmin,xmax,ymin,ymax)
        !call Axis(P,20,20)
        do i=0,U%nb-1
            wl = U%data(U%BList(i)%lab)
            if (wl>0) then
                Path(1)%x = U%PointList(U%BList(i)%p1)%x
                Path(2)%x = U%PointList(U%BList(i)%p2)%x
                Path(1)%y = U%PointList(U%BList(i)%p1)%y
                Path(2)%y = U%PointList(U%BList(i)%p2)%y
                call Line2d(P,Path,Thickness=0.00001_rpl)
            end if
        end do
        call Close(P)
    end subroutine plotbound

    subroutine plotboundi(U,ind,kam)

        use pmatypy
        use global
        use geometry
        use zobraz
        implicit none
        type(Uloha), intent(in) :: U
        integer, intent(in) :: ind
        integer, intent(in) :: kam
        real(kind=rpl) :: xmin, xmax,ymin,ymax
        real(kind=rpl) :: x1, x2,y1,y2
        integer :: wl
        integer(kind=ikind) :: i, cnt
        call meze(U,xmin,xmax,ymin,ymax)
        cnt = 0
        do i=0,ubound(U%BList,1)
            wl = U%BList(i)%lab
            if (wl == ind) then
                cnt = cnt + 1
                x1 = U%PointList(U%BList(i)%p1)%x
                x2 = U%PointList(U%BList(i)%p2)%x
                y1 = U%PointList(U%BList(i)%p1)%y
                y2 = U%PointList(U%BList(i)%p2)%y
            end if
        end do
        print *, "pocet usecek=",cnt
    end subroutine plotboundi

    !> \brief nakresli oblast
    !!
    !! \param U uloha
    !!
    subroutine plotdomain(U)
        use pmatypy
        use global
        use geometry
        use zobraz
        use pmatools
        implicit none
        type(Uloha), intent(in) :: U
        type(Picture_Type) :: P
        real(kind=rpl) :: xmin, xmax,ymin,ymax
        type(Point2d), dimension(1:3) :: Path
        type(RGB_Color_Type), dimension(:), pointer :: CMap
        integer :: wl,nc, mm
        integer(kind=ikind) :: i, cnt
        call meze(U,xmin,xmax,ymin,ymax)
        call Create(P,100.0_rpl,100*(ymax-ymin)/(xmax-xmin),&
        xmin,xmax,ymin,ymax)
        mm = lbound(U%permList,1)
        nc = ubound(U%permList,1)-mm+1
        print *,nc,mm
        call Create_Colormap(CMap,nc,1)
        CMap(1) = RGB_Color_Type(0,1,0) !vzduch
        CMap(11) = RGB_Color_Type(0,1,0)!vzduch
        CMap(2) = RGB_Color_Type(0,0.5,0.5) ! kartit
        cnt = 0
        do i=0,U%nt-1
            wl = U%TriList(i)%lab
            if ( .true.) then
                cnt = cnt + 1
                Path(1)%x = U%PointList(U%TriList(i)%p1)%x
                Path(2)%x = U%PointList(U%TriList(i)%p2)%x
                Path(3)%x = U%PointList(U%TriList(i)%p3)%x
                Path(1)%y = U%PointList(U%TriList(i)%p1)%y
                Path(2)%y = U%PointList(U%TriList(i)%p2)%y
                Path(3)%y = U%PointList(U%TriList(i)%p3)%y
                call Fill2d(P,Path,CMap(wl+1-mm))
                call Line2d(P,Path,Thickness=0.00001_rpl)
            end if
        end do
        call Close(P)
        print *, "vykresleno ",cnt, " domen"
    end subroutine plotdomain

    subroutine plotdomaini(U,ind,kam, mapovat)

        use pmatypy
        use global
        use geometry
        use zobraz
        implicit none
        type(Uloha), intent(in) :: U
        integer, intent(in) :: ind
        integer, intent(in) :: kam
        logical, intent(in) :: mapovat
        real(kind=rpl) :: xmin, xmax,ymin,ymax
        real(kind=rpl) :: x(1:3),y(1:3)
        real(kind=rkind) :: pl
        integer :: wl
        integer(kind=ikind) :: i, cnt
        call meze(U,xmin,xmax,ymin,ymax)
        cnt = 0
        pl = 0
        do i=0,U%nt - 1
            wl = U%TriList(i)%lab
            if (mapovat) then
                wl = U%data1(wl)
            end if

            if (wl == ind) then
                cnt = cnt + 1
                x(1) = U%PointList(U%TriList(i)%p1)%x
                x(2) = U%PointList(U%TriList(i)%p2)%x
                x(3) = U%PointList(U%TriList(i)%p3)%x
                y(1) = U%PointList(U%TriList(i)%p1)%y
                y(2) = U%PointList(U%TriList(i)%p2)%y
                y(3) = U%PointList(U%TriList(i)%p3)%y
                pl = pl + U%TriList(i)%Plocha
            end if
        end do
        print *, "pocet zobrazenych trojuhelniku:", cnt, &
            " z ",ubound(U%TriList,1)+1, &
            " celkova plocha=",pl
    end subroutine plotdomaini


    subroutine meze(U,xmin,xmax,ymin,ymax)

        use pmatypy
        use geometry
        use zobraz
        use global
        implicit none
        type(Uloha), intent(in) :: U
        real(kind=rpl), intent(out) :: xmin,xmax,ymin,ymax
        real(kind=rpl) :: wx,wy
        integer(kind=ikind) :: i,j
        Xmin = U%PointList(1)%x
        Ymin = U%PointList(1)%y
        Xmax = Xmin
        Ymax = Ymin
        do i=0,U%np-1
            wx = U%PointList(i)%x
            wy = U%PointList(i)%y
            if (wx < Xmin) Xmin = wx
            if (wx > Xmax) Xmax = wx
            if (wy < Ymin) Ymin = wy
            if (wy > Ymax) Ymax = wy
        end do
    end subroutine meze


end module kreslitko

