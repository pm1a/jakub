module global
    use geometry
    use pmatypy

    !> lokalni matice a dodatecne udaje
    !!
    !! Resena uloha je
    !! \f[
    !! -\nabla((\alpha+i*\beta)\nabla(u+i*v)=0
    !! \f]
    !! to lze prepsat na soustavu dvou rovnic
    !! \f{eqnarray*}{
    !!-\alpha\triangle u+\beta\triangle v & = & f\\
    !!+\beta\triangle u+\alpha\triangle v & = & 0
    !! \f}
    !<
    type, public :: Loc_Data
        !> koeficient u realne casti matice
        real(kind=rkind) :: alfa
        !> koeficient u imaginarni casti matice
        real(kind=rkind) :: beta
        !> diskretizace \f$ -\triangle u \f$
        real(kind=rkind), dimension(3,3) :: A
        !> prava strana, \n 1...realna cast\n 2...imaginarni cast
        real(kind=rkind), dimension(3,2) :: b
        !> cisla vrcholu - ve smyslu promennych
        integer(kind=ikind), dimension(3,2) :: map
        !> cislo trojuhelniku
        integer(kind=ikind) :: tid
    end type Loc_Data





    !> \brief drzi data o uloze
    !!
    !! !!!! Vzhledem k tomu,
    !! ze data pochazi z C nastroju, tak indexovani zacina nulou !!!
    type, public :: Uloha
        !> pocet bodu v diskretizaci
        integer(kind=ikind) :: np
        !> Seznam bodu v diskretizaci
        type(Point), public, dimension(:), allocatable    :: PointList
        !> pocet trojuhelniku v diskretizaci
        integer(kind=ikind) :: nt
        !> Seznam trojuhelnicku
        type(Triangle), public, dimension(:), allocatable :: TriList
        !> Pocet hranicnich hran
        integer(kind=ikind) :: nb
        !> seznam hran na okraji
        type(Bound), public, dimension(:), allocatable    :: BList
        !> zobrazi label hranice na typ hranice
        !> 0 ... vnitrek
        !> 1 ... dirichlet = 0
        !> 2 ... dirichlet = 1
        integer, dimension(:),allocatable :: data
        !> zobrazi label trojuhelniku na typ materialu
        !> 0 ... vzduch
        !> 1 ... cartit
        !> 2+ .... material(cislo-1)
        integer, dimension(:),allocatable :: data1

        type(BLabel), public, dimension(:), allocatable    :: BVList

        !> mapa znacek oblasti na mapu materialu
        type(TLabel), public, dimension(:), allocatable    :: permList
        !> ?
        !> ?
        integer(kind=ikind), dimension(:), allocatable, public :: Bindex
        !> pole lokalnich matic
        type(Loc_Data),dimension(:), allocatable :: LocList
    end type Uloha


    public :: ReadAgros
contains
    !> nacte geometrii z agrosu

    subroutine ReadAgros(basename, U)
        use pmatypy
        use pmatools
        implicit none
        character(len=*), intent(in) :: basename
        type(Uloha), intent(in out)   :: U
        integer :: fd
        integer(kind=ikind) :: np
        integer(kind=ikind) :: nt
        integer(kind=ikind) :: nb
        integer(kind=ikind) :: nbcnt
        integer(kind=ikind) :: i, i1,i2,i3,i4
        integer(kind=ikind) :: p1
        integer(kind=ikind) :: p2
        integer(kind=ikind) :: p3
        integer(kind=ikind) :: lab
        integer(kind=ikind) :: lab1
        integer(kind=ikind) :: lmin
        integer(kind=ikind) :: lmax
        integer(kind=ikind) :: pmax


        real(kind=rkind) :: x,y
        character :: znak
        character(*), parameter  :: n1 = ".edge"
        character(*), parameter  :: n2 = ".ele"
        character(*), parameter  :: n3 = ".node"
        character(*), parameter  :: n4 = ".poly"

        fd = 20
        open(unit=fd,file=basename // n3,action='READ', err = 1001)
        goto 1002
    1001 print *,"chyba pri otevirani souboru"
         stop "zachycena chyba pri otevirani" 
    1002 continue


        read (unit=fd,fmt=*) np
        U%np = np
        print *, "pocet bodu=",np
        allocate(U%PointList(0:np-1))
        do i=0,np-1
            read(unit=fd,fmt=*) i1,x,y,i2
            U%PointList(i)%x = x
            U%PointList(i)%y = y
            U%PointList(i)%onBoundary = .false.
            U%PointList(i)%value = 0.0_rkind
            U%PointList(i)%used = .false.
        end do
        close(fd)
        pmax = 0
        open(unit=fd,file=basename//n2,action='READ')
        read (unit=fd,fmt=*) nt
        print *, "pocet trojuhelniku=",nt
        U%nt = nt
        allocate(U%TriList(0:nt-1))
        lmin = nt
        lmax = 0
        do i=0,nt-1
            read(unit=fd,fmt=*) i1, p1,p2,p3, i2,i3,i4,lab
            p1 = p1
            p2 = p2
            p3 = p3
            U%TriList(i)%p1 = p1
            U%TriList(i)%p2 = p2
            U%TriList(i)%p3 = p3
            U%TriList(i)%lab = lab
            if (lab < lmin) lmin=lab
            if (lab > lmax) lmax=lab
            if (p1 > pmax ) pmax = p1
            if (p2 > pmax ) pmax = p2
            if (p2 > pmax ) pmax = p3
            U%PointList(p1)%used = .true.
            U%PointList(p2)%used = .true.
            U%PointList(p3)%used = .true.
        end do
        close(fd)
        print *,"labely od ",lmin," do ",lmax
        allocate(U%permlist(lmin:lmax))
        allocate(U%data1(lmin:lmax))
        ! nastavit data1 - zobrazi labely trojukelniku  na typ materialu
        ! tj. relativni permitivitu a gamma
        ! pro ted vsechno na vzduch
        U%data1 = 0
        do i=0,nt-1
            U%permlist(U%TriList(i)%lab)%cnt = &
                U%permlist(U%TriList(i)%lab)%cnt + 1
        end do
        print "(i6,f25.17,f25.17,i6)", (i,U%permlist(i), i=lmin,lmax)
        call pockej("po trojuhelnicich")

        open(unit=fd,file=basename//n1,action='READ')

        read (unit=fd,fmt=*) nb
        print *, "pocet hranic=",nb
        allocate(U%BList(0:nb-1))
        lmin = 0
        lmax = 0
        U%nb = nb
        do i=0,nb-1
            read(unit=fd,fmt=*) i1,p1,p2,lab
            U%BList(i)%p1 = p1
            U%BList(i)%p2 = p2
            U%BList(i)%lab = lab
            if (lab < lmin) lmin=lab
            if (lab > lmax) lmax=lab
        end do
        print *,"labely od ",lmin," do ",lmax
        call pockej
        allocate(U%BVlist(lmin:lmax))
        allocate(U%data(lmin:lmax))
        U%data = 0
        U%data(11:16) = 1 ! vnejsi hranice
        U%data(17:20) = 2 ! prvni elektroda
        U%data((/21,28,29,30/)) = 3 ! druha elektroda
        U%data((/22,31,32,33/)) = 4 ! treti elektroda
        U%data((/23,34,35,36/)) = 5 ! ctvrta elektroda
        U%data((/24,37,38,39/)) = 6 ! pata elektroda
        U%data((/25,40,41,42/)) = 7 ! sesta elektroda
        U%data((/26,43,44,45/)) = 8 ! sedma elektroda
        U%data((/27,46,47,48/)) = 9 ! osma elektroda
        U%data(51:130) = 10 !deska

        do i=0,nb-1
            U%BVlist(U%BList(i)%lab)%cnt = U%BVlist(U%BList(i)%lab)%cnt + 1
            U%BVlist(U%BList(i)%lab)%val = 0
        end do
        print "(i6,f25.17,i6)", (i,U%BVlist(i), i=lmin,lmax)
        close(fd)
        print *,  U%data
        print *, "pmax=", pmax
        LastUsedPoint = pmax



    end subroutine ReadAgros

end module global

