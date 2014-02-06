!> \file mtx.f90
!! obecna matice a operace
!<


!> definuje typy a oprace potrebne pro praci s maticemi
!! Typ matice je vlastnosti radku.
!! importuje module defs ( soubor defs.f90) - definice kind konstant.
!<
module mtx_jakub
    use pmatypy
    implicit none

    !> podrobnost protokolu
    integer :: debuglevel = 1

    !> typ definujici rozsah ulozenych indexu
    !!
    !! from ... poloha prvniho objektu \n
    !! to   ... poloha posledniho objektu \n
    !<  ..
    type, public :: rozsah
        !> prvni index
        integer(kind=ikind) :: from
        !> posledni index
        integer(kind=ikind) :: to
    end type rozsah



    !> typ definujici radek v matici
    !!
    !!
    !! umoznuje reprezentovat prazdny, ridky, obrysovy a plny radek.
    !! Vzdy se alokuji jen nezbytne casti.
    !! Objekty tohoto typu se vzdy odkazuji na matici za ktere pochazeji.
    !! Z techto duvodu musi byt privatni a jejich uziti obezretne.
    !<
    type,private :: Row_Type
        !> typ reprezentovaneho radku \n 0 .... prazdny radek \n
        !! 1 .... ridky radek \n
        !! 2 .... plny interval
        !<
        integer                                    :: druh = 0
        !> delka radku - t.j. nejvyssi povoleny index.
        !! Normalne pocet sloupcu matice
        !<
        integer(kind=ikind)                 :: len  = 0
        !> Pocet ulozenych prvku.
        integer(kind=ikind)                 :: nz   = 0
        !> index prvniho nenuloveho prvku
        integer(kind=ikind)                 :: cmin = 0
        !> index posledniho nenuloveho prvku
        integer(kind=ikind)                 :: cmax = 0

        !> rozsah indexu pro hodnoty
        type(rozsah) :: vrange = rozsah(0,0)
        !> rozsah pro sloupcove indexy
        type(rozsah) :: cirange  = rozsah(0,0)
        !> zpetny odkaz na values vlastnika
        real(kind=rkind), dimension(:), pointer :: ownerv
        !> zpetny odkaz na sloupcove indexy vlastnika
        integer(kind=ikind), dimension(:), pointer :: ownerc
    end type


    !> definice obecne matice
    type, public :: Mtx_Type
        private
        !> pocet radku
        integer(kind=ikind) :: n = 0
        !> pocet sloupcu
        integer(kind=ikind) :: m = 0
        !> hodnoty prvku matice
        real(kind=rkind), dimension(:), pointer :: v => null()
        !> hodnoty sloupcovych indexu
        integer(kind=ikind), dimension(:), pointer :: ci => null()
        !> pole radku
        type(Row_Type), dimension(:),pointer :: Rows => null()
        !> mnoziny der v hodnotach
        type(rozsah), dimension(50) :: freev
        !> mnoziny der v indexech sloupcu
        type(rozsah), dimension(50) :: freec
        !> pocet aktivnich der hodnot
        integer(kind=ikind) :: nfv
        !> pocet aktivnich der indexu
        integer(kind=ikind) :: nfc

    end type

    public :: mtx_test
    public :: init
    public :: Print_Mtx
    public :: Nonzero
    public :: Set_Elem
    private :: Set_Matrix_test
    private :: Get_Row
    private :: Set_Row
    private :: Set_Elem_Row
    private :: Get_Elem_Row
contains

    !> tester pro matice
    subroutine mtx_test
        use pmatypy
        implicit none
        type(Mtx_Type) :: A
        logical :: skoncit
        integer :: choice
        integer(kind=ikind) :: n,m

        skoncit = .false.
        print *," zacina testovani matic"
        do while ( .not. skoncit)
            print *, "vyber akci"
            print *, " 0 ... konec"
            print *, " 1 ... inicializace matice"
            print *, " 2 ... tisk matice"
            print *, " 3 ... zadavani matice"
            read(*,*) choice
            select case(choice)
                case (0)
                    skoncit = .true.
                case (1)
                    print *, "zadej pocet radku"
                    read(*,*) n
                    print *, "zadej pocet sloupcu"
                    read(*,*) m
                    call init(A,n,m)
                case (2)
                    call Print_Mtx(A,3,"pokusna matice")
                case (3)
                    call Set_Matrix_test(A)
                case default
                    print *,"chybna volba"
            end select
        end do
    end subroutine mtx_test


    !> zakladni konstruktor
    !!
    !! vytvori prazdnou  matici typu (n x m) s prostorem
    !! pro aspon nz nenulovych prvku
    !!
    !<
    subroutine init(A,n,m,nz)
        implicit none
        !> matice - dojde k pripadne realokaci
        type(Mtx_Type), intent(in out)           :: A
        !> pocet radku
        integer(kind=ikind),intent(in)           :: n
        !> pocet sloupcu - neni-li zadano je to pocet radku
        integer(kind=ikind),intent(in), optional :: m
        !>  pripravi prostor pro nz nenulovych prvku (zatim neuzito)
        integer(kind=ikind),intent(in), optional :: nz


        integer(kind=ikind) :: lm,lnz,i
        !> na zacatku vyresit nepovinne parametry
        if (present(m)) then
            lm = m
        else
            lm = n
        end if
        if (present(nz)) then
            lnz = nz
        else
            lnz =max(n,m)
        end if

        if (debuglevel > 0 ) then
            print *,"mtx :: Init"
            print *, "n=",n," m=",lm," nz=",lnz
        end if;
        !> napred udelat pripadanou realokaci
        !print *,"jdu testovat"
        if (associated(A%Rows)) then
            deallocate(A%Rows)
            print *, "uvolnuji"
        end if
        !print *, "jdu alokovat"
        allocate(A%Rows(1:n))
        A%n = n
        A%m = lm
        !print *,"vracim se"
        !> tady pripadne pripravit pracovni prostory
        A%nfc = 1
        A%nfv = 1
        if (associated(A%v)) then
            deallocate(A%v)
        end if
        allocate(A%v(1:lnz))
        if (associated(A%ci)) then
            deallocate(A%ci)
        end if
        allocate(A%ci(1:lnz))
        do i = 1,n
            A%Rows(i)%ownerv => A%v
            A%Rows(i)%ownerc => A%ci
        end do
        A%freev(1)%from = 1
        A%freev(1)%to   = lnz
        A%freec(1)%from = 1
        A%freec(1)%to   = lnz
    end subroutine init


    !> ziska radek z matice
    subroutine Get_Row(A,I,Row)
        type(Mtx_Type), intent(in) :: A
        integer(kind=ikind), intent(in) :: I
        type(Row_Type), intent(in out) :: Row
    end subroutine Get_Row

    !> nahradi radek matice
    subroutine Set_Row(A,I,Row)
        type(Mtx_Type), intent(in) :: A
        integer(kind=ikind), intent(in) :: I
        type(Row_Type), intent(in out) :: Row
    end subroutine Set_Row

    !> tiskne matici
    subroutine Print_Mtx(A,ncol,Caption)
        !> tistena matice
        type(Mtx_Type), intent(in)    :: A
        !> pocet sloupcu ve vypisu
        integer, intent(in), optional :: ncol
        !> popisek k matici
        character(len=*), intent(in), optional :: Caption

        integer :: lncol

        if (present(ncol)) then
            lncol = ncol
        else
            lncol = 5
        end if

        if (present(Caption)) then
            print *, "matice ", Caption
        else
            print *, "matice"
        end if

        print *, "pocet radku=",A%n," pocet sloupcu=",A%m, " nenul=", Nonzero(A)
    end subroutine Print_Mtx

    !> textova simulace funkce spy
    subroutine Spy(A)
        !> tistena matice
        type(Mtx_Type), intent(in) :: A
    end subroutine Spy

    function Nonzero(A) result(Nz)
        use pmatypy
        implicit none
        type(Mtx_Type), intent(in) :: A
        integer(kind=ikind) :: Nz

        Nz = 0
    end function Nonzero

    subroutine Set_Elem_Row(R,I,Val)
        use pmatypy
        implicit none
        type(Row_Type), intent(in out) :: R
        Integer(kind=ikind), intent(in) :: I
        Real(kind=rkind), intent(in) :: Val
    end subroutine Set_Elem_Row

    function Get_Elem_Row(R,I) result(Val)
        use pmatypy
        implicit none
        type(Row_Type), intent(in) :: R
        Integer(kind=ikind), intent(in) :: I
        Real(kind=rkind) :: Val

        Val = 0
    end function Get_Elem_Row


    subroutine Set_Matrix_test(A)
        use pmatypy
        implicit none
        type(Mtx_Type), intent(in out) :: A

        logical :: skoncit
        integer :: choice
        integer(kind=ikind) :: row, col
        real(kind=rkind) :: value

        skoncit = .false.
        print *, "cteni matice"
        do while (.not. skoncit)
            print *, "vyber cinnost"
            print *, " 0 ... skoncit"
            print *, " 1 ... zadavat z klavesnice"
            print *, " 2 ... precist ze souboru"
            read(*,*) choice
            select case(choice)
                case (0)
                    skoncit = .true.
                case (1)
                    row = 1
                    col = 1
                    value = 0
                    do while (row > 0)
                        print *, "zadej radek"
                        read(*,*) row
                        print *, "zadej sloupec"
                        read(*,*) col
                        print *,"zadej hodnotu prvku"
                        read(*,*) value
                        call Set_Elem(A,value,row,col)
                    end do
                case (2)
                case default
                    print *, "chybna volba"
            end select
        end do

    end subroutine Set_Matrix_test

    subroutine Set_Elem(A,value,row,col)
        use pmatypy
        implicit none
        type(Mtx_Type), intent(in out) :: A
        real(kind=rkind), intent(in) :: value
        integer(kind=ikind), intent(in) :: row
        integer(kind=ikind), intent(in) :: col

        integer(kind=ikind) :: ic,iv
        ! 1. zkus ho najit
        call Locate(A,row,col,iv,ic)
        ! 2. pokud tam je prepis
        if ( iv > 0 ) then ! tohle znamena nasel jsem
            if ( value == 0) then
                ! nulova hodnota znamena vymazani prvku z radku
                call Delete(A,row,col,iv,ic)
            else
                ! pokud to neni nula, tak jen prepis hodnotu
                A%v(iv) = value
            end if
        else
        ! 3. pokud neni vloz
           call Insert(A,row,col,value)
        end if
    end subroutine Set_Elem

    subroutine Locate(A,row,col,iv,ic)
       use pmatypy
       implicit none
       type(Mtx_Type),intent(in) :: A
       integer(kind=ikind),intent(in) :: row
       integer(kind=ikind),intent(in) :: col
       integer(kind=ikind),intent(out) :: iv
       integer(kind=ikind),intent(out) :: ic

       type(Row_Type) :: R

       iv = 0
       ic = 0
       ! otesuju jestli je v matici, napred radky
       if ((row < 1) .or. (row > A%n)) then
          return
       end if
       ! otesuju jestli je v matici, potom sloupce
       if ((col < 1) .or. (col > A%m)) then
          return
       end if
       ! vezmu si radek a zacnu zkoumat
       R = A%Rows(row)
       select case(R%druh)
          case (0)
             return !neni co delat, neni to tam
          case (1) ! budu hledat, je to ridky
             if ( (col < R%cmin) .or. (col > R%cmax)) then
                return !prece jen to tam neni
             end if
             ! je to sice utrideny, ale zatim to nevyuzivame
             iv = R%vrange%From
             ic = R%cirange%From
             do while (col /= A%ci(ic))
                iv = iv + 1
                ic = ic + 1
             end do
          case (2) !je to plny
             if ( (col < R%cmin) .or. (col > R%cmax)) then
                return !prece jen to tam neni
             end if
             ! je tam, ted ho oznacit
             ! ic ponechat nula jako znacku plnosti
             iv = R%vrange%from + (col - R%cmin)
          case default
             print *, "chybny typ radku"
             stop "Velka chyba implementace"
       end select

    end subroutine Locate



    subroutine Insert(A,row,col,value)
       use pmatypy
       implicit none
       type(Mtx_Type), intent(in out) :: A
       integer(kind=ikind),intent(in) :: row
       integer(kind=ikind),intent(in) :: col
       real(kind=rkind), intent(in)   :: value
       ! uz vim, ze tam neni - budu urcite vkladat
       type(Row_Type) :: R
       type(Row_Type) :: Rnew
       R = A%Rows(row)

       select case (R%druh)
         case (0) ! je prazdnej - proste vytvorime

         case (1) ! je ridkej - musime ho prodlouzit
         case (2) ! je plnej - prodlouzit, pripadne preves na ridkej
       end select

    end subroutine Insert

    !> najde volne misto v dirach, vrati ho a opravi diry
    !! volny rozsah doda v R
    subroutine Najdi_Misto(Delka,Diry,Nd,R)
       use pmatypy
       implicit none
       !> pozadovana delka pole
       Integer(kind=ikind), intent(in) :: Delka
       !> seznam der
       type(rozsah), dimension(:), intent(in out) :: Diry
       !> pocet der
       Integer(kind=ikind), intent(in out)            :: Nd
       !> nalezeny interval, pokud neuspeje
       !! vrati v polozce from: \n
       !! >0 ... prvni index hledaneho prostoru
       !!  0 ... zadna dira neni, ale v souctu se vejde
       !! <0 ... kolik je treba pridat
       type(Rozsah), intent(out)                  :: R

       Integer(kind=ikind) :: i
       Integer(kind=ikind) :: suma
       type(Rozsah) :: Wrk
       suma = 0
       do i=1,Nd
          Wrk = Diry(i)
          if (Wrk%To - Wrk%From + 1>= Delka) then
             ! mam to misto
             R%From = Wrk%From              ! tak to nastav
             R%To   = Wrk%From + Delka - 1
             Diry(i)%From = R%To + 1        ! a uprav diru
             return                        ! je hotovo - zpatky
          else
             ! nemam, tak aspon opravim celkem dostupne
             suma = suma + Wrk%To - Wrk%From + 1
          end if
       end do
       ! pokud jsem tady tak to nemam
       if (suma >= Delka) then
          R%From = 0
       else
          R%From = suma - Delka -1
       end if
    end subroutine  Najdi_Misto

    subroutine Delete(A,row,col,iv,ic)
       use pmatypy
       implicit none
       type(Mtx_Type),intent(in) :: A
       integer(kind=ikind),intent(in) :: row
       integer(kind=ikind),intent(in) :: col
       integer(kind=ikind),intent(out) :: iv
       integer(kind=ikind),intent(out) :: ic

    end subroutine Delete

end module mtx_jakub
