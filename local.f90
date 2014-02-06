!> \file local.f90
!! \brief vytvari lokalni matice tuhosti


!> module pro lokalni matice
module local
    use pmatypy
    implicit none

    public  :: CreateLocal
    private :: locmat

contains
    subroutine locmat(U,Id)
        use pmatypy
        use geometry
        use global
        implicit none
        type(Uloha), intent(in out) :: U
        integer(kind=ikind), intent(in) :: Id
        integer(kind=ikind) :: lab, permid,i,p
        real(kind=rkind) :: det, plocha
        real(kind=rkind), dimension(3,3) :: A,B
        character :: zn

        lab = U%TriList(id)%lab
        permid = U%data1(lab)
        U%LocList(id)%alfa = gama(permid)
        U%LocList(id)%beta = omega*eps0*epsr(permid)
        U%LocList(id)%b = 0.0_rkind !vetsinou je prava strana nula
        U%Loclist(id)%map(:,1)= &
            (/ U%TriList(id)%p1, U%TriList(id)%p2,U%TriList(id)%p3 /)
        U%Loclist(id)%map(:,2)= &
            (/ U%TriList(id)%p1, U%TriList(id)%p2,U%TriList(id)%p3 /)
        ! ted jsou mapy stejne a vetsinou zustanou
        do i=1,3
            A(i,3)=1
            A(i,1)= U%PointList(U%LocList(id)%map(i,1))%x
            A(i,2)= U%PointList(U%LocList(id)%map(i,1))%y
        end do
        det = (A(2,1)-A(1,1))*(A(3,2)-A(1,2))-(A(3,1)-A(1,1))*(A(2,2)-A(1,2))
        ! jeste vydelit plochou
        plocha = abs((A(2,1)-A(1,1))*(A(3,2)-A(1,2))-(A(3,1)-A(1,1))*(A(2,2)-A(1,2)))/2
        B(1,1) =  (A(2,2)-A(3,2))/det
        B(1,2) = -(A(1,2)-A(3,2))/det
        B(1,3) =  (A(1,2)-A(2,2))/det
        B(2,1) = -(A(2,1)-A(3,1))/det
        B(2,2) =  (A(1,1)-A(3,1))/det
        B(2,3) = -(A(1,1)-A(2,1))/det
        B(3,1) =  (A(2,1)*A(3,2)-A(3,1)*A(2,2))/det
        B(3,2) = -(A(1,1)*A(3,2)-A(3,1)*A(1,2))/det
        B(3,3) =  (A(1,1)*A(2,2)-A(2,1)*A(1,2))/det
        ! ve sloupcich B jsou koeficienty jednotlivych lokalnich bazovych funkci
        ! poradi je radek 1. x  2. y  3. konstanta
        U%LocList(id)%A(1,1) = B(1,1)*B(1,1)+B(2,1)*B(2,1)
        U%LocList(id)%A(1,2) = B(1,1)*B(1,2)+B(2,1)*B(2,2)
        U%LocList(id)%A(1,3) = B(1,1)*B(1,3)+B(2,1)*B(2,3)
        U%LocList(id)%A(2,1) = U%LocList(id)%A(1,2)
        U%LocList(id)%A(2,2) = B(1,2)*B(1,2)+B(2,2)*B(2,2)
        U%LocList(id)%A(2,3) = B(1,2)*B(1,3)+B(2,2)*B(2,3)
        U%LocList(id)%A(3,1) = U%LocList(id)%A(1,3)
        U%LocList(id)%A(3,2) = U%LocList(id)%A(2,3)
        U%LocList(id)%A(3,3) = B(1,3)*B(1,3)+B(2,3)*B(2,3)
        ! jsou to integraly z konstant pres trojuhelnik s nejakou plochou
        U%LocList(id)%A=U%LocList(id)%A*plocha
        U%TriList(Id)%plocha = plocha
        if (det == 0 ) then
            print *, id, U%LocList(id)
            print *,"A"
            print *, A
            print *,"B"
            print *, B

            stop
        end if
        !ted doresit rhs
        do i=1,3
            p=U%Loclist(id)%map(i,1)
            if (U%PointList(p)%onBoundary) then
                U%LocList(Id)%b(:,1) = U%LocList(Id)%b(:,1) - &
                    U%LocList(Id)%A(:,i)*&
                    U%PointList(p)%value(1)*U%LocList(id)%beta -&
                    U%LocList(Id)%A(:,i)*&
                    U%PointList(p)%value(2)*U%LocList(id)%alfa
                    U%LocList(Id)%b(:,2) = U%LocList(Id)%b(:,2) - &
                    U%LocList(Id)%A(:,i)*U%PointList(p)%value(1) &
                    *U%LocList(id)%alfa +&
                    U%LocList(Id)%A(1:3,i)*U%PointList(p)%value(2)!                    LocList(id)%beta
                U%Loclist(id)%map(i,1) = 0
            end if
        end do
    end subroutine locmat


    subroutine CreateLocal(U)
        use pmatypy
        use geometry
        use global
        implicit none

        type(Uloha), intent(inout) :: U
        integer(kind=ikind) :: nt,i,j,k
        integer :: fd
        print *, "create local"
        nt = ubound(U%TriList,1)+1
        allocate(U%LocList(0:nt-1))
        do i=0,nt-1
            call locmat(U,i)
        end do
        print *,"jdu to nasypat do souboru"
        fd = 22
        open(unit=fd,file='locm.txt',action='WRITE')
        write(unit=fd,fmt=*) nt
        do i=0,nt-1
            write(unit=fd, fmt=*) U%Loclist(i)%map(:,1)
            write(unit=fd, fmt=*) U%Loclist(i)%map(:,2)
            do j=1,3
                do k=1,3
                    write(unit=fd,fmt=*) &
                        U%Loclist(i)%map(j,1), &
                        U%Loclist(i)%map(k,1), &
                        U%Loclist(i)%A(j,k)
                end do
                write(unit=fd, fmt=*) &
                    U%Loclist(i)%map(j,1), &
                    U%Loclist(i)%b(j,1), &
                    U%Loclist(i)%b(j,2)
            end do
            write(unit=fd, fmt=*) &
                U%Loclist(i)%alfa,U%Loclist(i)%beta, U%TriList(i)%lab
        end do
        ! ted jeste zapis okrajovou podminku
        write( unit = fd, fmt = *),ubound(U%Blist,1)
        do i=1, ubound(U%Blist,1)
            write( unit = fd, fmt = *),U%Blist(i)
        end do

        close(fd)

    end subroutine CreateLocal

end module local
