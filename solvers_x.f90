!> \file solvers.f90
!! jednoduche resice

module solvers 
    public :: SDS
    public :: CGS
    public :: Axpy
    public :: AAxpy
    public :: getRhs
    public :: printMTX
contains 
    !> metoda nejvetsiho spadu pro A*A, A je synetricka, indefinitni\n
    !! resi soustavu A*x=b\n
    !! tak,ze ji prepise na (A*A)*x=A*b
    subroutine SDS(x,b)
        use typy
        implicit none
        !> aproximace reseni
        real(kind=rkind),dimension(0:,:), intent(inout) :: x
        !> prava strana
        real(kind=rkind),dimension(0:,:), intent(in) :: b

        integer(kind=ikind) :: n,cnt,rcnt, pcnt ,stepcnt
        real(kind=rkind),dimension(0:ubound(x,1),1:2) :: y,w1,w2,xx
        real(kind=rkind) :: rn, rar, rnold, rnold1, rmin, rn1, rar1
      
        n = ubound(x,1)
        x(0,:) = 0
        cnt = 0
        rcnt = 0
        stepcnt = 1000
        pcnt = stepcnt
        do
            cnt = cnt + 1
            y=0
            call Axpy(x,y)
            y = b - y  ! y = b - A*x === residum obycejne
            w1 = 0
            call Axpy(y,w1)  ! w1 = A*(b-A*x) == residuum z definice problemu
            rn = DOT_PRODUCT(w1(:,1),w1(:,1))+DOT_PRODUCT(w1(:,2),w1(:,2))
            rn1 = DOT_PRODUCT(y(:,1),y(:,1))+DOT_PRODUCT(y(:,2),y(:,2))
            if ( rn1 == 0 ) return
            if (pcnt == stepcnt) then
                xx = x
                print "(a,i8,a,1pg22.15,a,i4,a,1pg22.15,1pg22.15)",&
                "SDS iterace=",cnt," norma residua=", sqrt(rn1),&
                " rcnt=",rcnt," integral = ",integral(xx)
                print *,  rar1, sqrt(rn1)
                pcnt = 1
            else
                pcnt = pcnt + 1
            end if
            rnold1 = rnold
            rnold = rn1
            if (cnt < 2) rmin = rnold
            if (rnold < rmin) then
                rcnt = 0
                rmin = rnold
            else
                rcnt = rcnt + 1
            end if
            if (rcnt == 2000) return
            y = 0
            w2 = 0
            call Axpy(w1,y)
            !call Axpy(y,w2)
            !rar = DOT_PRODUCT(w1(:,1),w2(:,1))+DOT_PRODUCT(w1(:,2),w2(:,2))
            rar1 = DOT_PRODUCT(y(:,1),y(:,1))+DOT_PRODUCT(y(:,2),y(:,2))
            rn = rn/rar1  ! tohle staci
            x = x + rn*w1
        end do
    end subroutine SDS
   

    !> realizuje operaci y = A*x+y
    !! nasobeni provadi pomoci lokalnich matic
    subroutine Axpy(x,y)
        use typy
        use geometry
        use local
        implicit none
        !> x
        real(kind=rkind), dimension(0:,:), intent(in) :: x
        !> y
        real(kind=rkind), dimension(0:,:), intent(inout) :: y
        integer(kind=ikind) :: nt,i,j
        real(kind=rkind), dimension(3,2) :: lx, ly
        real(kind=rkind), dimension(3,3) :: lA
        real(kind=rkind) :: alfa, beta
        integer(kind=ikind) , dimension(1:3) :: mapa
      
        character :: zn
      
        nt = ubound(LocList,1)
        if ( (x(0,1) /= 0) .or. (x(0,2) /= 0) ) then
            print *, " chybna data pro Axpy"
            stop
        end if
        !print *,"axpy"
        do i=1,nt
            alfa = LocList(i)%alfa
            beta = LocList(i)%beta
            la   = LocList(i)%A
            mapa = LocList(i)%map(:,1)
            ly = x(mapa,:)
            do j=1,3
                if(mapa(j)==0)  then
                    ly(j,1) = 0
                    ly(j,2) = 0
                    la(j,:) = 0
                    la(:,j) = 0
                else
                    if (Pointlist(mapa(j))%onboundary) then
                        print *, "axpy  - potiz, pro i=",i
                        ly(j,1) = 0
                        ly(j,2) = 0
                        la(j,:) = 0
                        la(:,j) = 0
                    end if
                end if
            end do
            lx(:,1) = matmul(la,ly(:,1)) !to mam jen souciny s matici, ted jeste koeficienty
            lx(:,2) = matmul(la,ly(:,2)) !to mam jen souciny s matici, ted jeste koeficienty
            ly(:,1) =   beta*lx(:,1)+alfa*lx(:,2)
            ly(:,2) =   alfa*lx(:,1)-beta*lx(:,2)
            !print *,ly
            do j=1,3
                y(mapa(j),1) = y(mapa(j),1) + ly(j,1)
                y(mapa(j),2) = y(mapa(j),2) + ly(j,2)
            end do
            y(0,:) = 0
        end do
        y(0,:)=0
        y(Bindex,:)=0
        do i=1,ubound(y,1)
            if ( .not. Pointlist(i)%used) y(i,:) = 0
            if ( Pointlist(i)%onboundary) y(i,:) = 0
        end do
    end subroutine Axpy


    subroutine AAxpy(x,y)
        use typy
        use geometry
        use local
        implicit none
        real(kind=rkind), dimension(0:,:), intent(in) :: x
        real(kind=rkind), dimension(0:,:), intent(inout) :: y
        real(kind=rkind), dimension(0:ubound(x,1),1:2) :: wrk

        wrk = 0
        call Axpy(x,wrk)
        call Axpy(wrk,y)
    end subroutine AAxpy

   
    subroutine getRhs(x)
        use typy
        use geometry
        use local
        implicit none
      
        real(kind=rkind), dimension(0:,:), intent(out) :: x
        integer(kind=ikind) :: nt,i
        nt = ubound(Trilist,1)
        x=0
        do i=1,nt
            x(LocList(i)%map(:,1),:) =&
            x(LocList(i)%map(:,1),:) + LocList(i)%b
        end do
        x(0,:) = 0
        do i = 1, ubound(x,1)
            if ( .not. PointList(i)%used ) then
                x(i,:)=0
            end if
            if (PointList(i)%onboundary) then
                x(i,:) = 0
                print *, " set rhs - chyba", i
            end if
        end do
        x(0,:)  = 0
        ! pro pokus to uplne vymazeme
        !x = 0
        !x(3,:) = 1
    end subroutine getRhs
   
   
    function integral(x) result (y)
        use typy
        use geometry
        use local
        implicit none
        real(kind=rkind), dimension(0:,:), intent(inout) :: x
        real(kind=rkind), dimension(1:2) :: y
        real(kind=rkind) :: wrk
        real(kind=rkind), dimension(1:3) :: xl,x2
        integer(kind=ikind) :: nt,nb,i
      
        y = 0
        nt = ubound(LocList,1)
        nb = ubound(BIndex,1)
        !napsat do x spravne hodnoty na hranici
        do i=1,nb
            X(BIndex(I),:) = PointList(Bindex(i))%value
        end do
        do i=1,nt
            xl = x(LocList(i)%map(:,2),1)
            x2 = x(LocList(i)%map(:,2),2)
            wrk = DOT_PRODUCT(xl,matmul(LocList(i)%A,xl))*LocList(i)%beta
            y(1) = y(1) + wrk
            wrk = DOT_PRODUCT(x2,matmul(LocList(i)%A,x2))*LocList(i)%alfa
            y(2) = y(2) + wrk
           !print *, i, wrk
        end do
        !zase je vynulovat
        x(Bindex,:)=0
    end function integral
   
   
   
    subroutine CGS(x,b)
        use typy
        implicit none
        real(kind=rkind),dimension(0:,:), intent(inout) :: x
        real(kind=rkind),dimension(0:,:), intent(in)    :: b

        integer(kind=ikind) :: n,cnt,rcnt, pcnt ,stepcnt, n0
        real(kind=rkind),dimension(:,:), allocatable :: y,w1,w2,sm,xx
        real(kind=rkind) :: rn, rnr, rar, rnold, rnold1, rmin, resnorm
      
        n = ubound(x,1)
        n0 = lbound(x,1)
        if (n0 /= 0 ) print *, "CGS podivne"
        allocate(y(0:n,1:2))
        allocate(xx(0:n,1:2))
        allocate(w1(0:n,1:2))
        allocate(w2(0:n,1:2))
        allocate(sm(0:n,1:2))
        cnt = 0
        rcnt = 0
        stepcnt = 100
        pcnt = stepcnt
        y = 0
        do
            cnt = cnt + 1
        
            w2=0
            call Axpy(x,w2)
            w2 = b - w2
            resnorm  = DOT_PRODUCT(w2(:,1),w2(:,1))+DOT_PRODUCT(w2(:,2),w2(:,2))
            if (resnorm == 0) return
            w1 = 0
            call Axpy(w2,w1)  ! w1 = A*(b-A*x) == residuum
            if (cnt == 1) then
                sm = w1
            else
                ! udelame ortogonalizaci misto ponechani sm = w1
                !y = 0
                w2 = 0
                call Axpy(w1,w2)  ! w2 = A*r
                rn  = DOT_PRODUCT(y(:,1),w2(:,1))+ DOT_PRODUCT(y(:,2),w2(:,2)) ! rn=sm'*A*A*r
                rnr = rar    !rar=sm'*A'*A*sm
                rn = rn/rnr
                if (rnr == 0) then
                    print *,"problem"
                    stop
                end if
                sm = w1 - rn*sm
            end if
            rn  = DOT_PRODUCT(sm(:,1),sm(:,1))+DOT_PRODUCT(sm(:,2),sm(:,2)) ! rn  = sm'*sm
            rnr = DOT_PRODUCT(w1(:,1),w1(:,1))+DOT_PRODUCT(w1(:,2),w1(:,2)) ! rnr = r'*r

            if ( (pcnt == stepcnt) .or. (rcnt > 0) ) then
                xx = x
                print "(a,i8,a,i8,a,1pg22.15,a,i6,a,1pg22.15,1pg22.15,a,1pg22.15,a,1pg22.15,1pg22.15)",&
                "CGS iterace=",cnt," z celkem ",2*n," norma residua=", sqrt(resnorm),&
                " rcnt=",rcnt," integral = ",integral(xx)," delka smeru=", sqrt(rn), "nejmensi:",rmin
                pcnt = 1
            else
                pcnt = pcnt + 1
            end if
            rnold1 = rnold
            rnold = sqrt(resnorm)
            if (cnt < 2) rmin = rnold
            if (rnold < rmin) then
                rcnt = 0
                rmin = rnold
            else
                rcnt = rcnt + 1
            end if
            if (rcnt == 20) return
            y = 0
            call Axpy(sm,y)   ! y=A*sm
            rn =  DOT_PRODUCT(w1(:,1),sm(:,1))+DOT_PRODUCT(w1(:,2),sm(:,2)) ! rn = r'*sm
            rar = DOT_PRODUCT(y(:,1),y(:,1))+DOT_PRODUCT(y(:,2),y(:,2)) !rar=sm'*A'*A*sm
            rn = rn/rar
            x = x + rn*sm
        end do
    end subroutine CGS
   
    subroutine printMTX
        use typy
        use geometry
        use local
        implicit none

        integer :: fd
        integer(kind=ikind) :: n,i,j,pocet
        real(kind=rkind),dimension(:,:), allocatable :: y,x

        n = LastUsedPoint
        fd = 23
        allocate(x(0:n,1:2))
        allocate(y(0:n,1:2))
        pocet = 0
        open(unit=fd,file='matice.txt',action='WRITE')
        write(unit=fd,fmt=*) n
        do i=1,n
            x = 0
            y = 0
            x(i,1) = 1
            call Axpy(x,y)
            do j = 1,n
                if ( y(j,1) /= 0 ) then
                    write(unit=fd,fmt=*) i,j, y(j,1)
                    pocet = pocet + 1
                end if
                if ( y(j,2) /= 0 ) then
                    write(unit=fd,fmt=*) i,j+n, y(j,2)
                    pocet = pocet + 1
                end if
            end do
            print *, i,pocet
        end do
        do i=1,n
            x = 0
            y = 0

            x(i,2) = 1
            call Axpy(x,y)
            do j = 1,n
                if ( y(j,1) /= 0 ) then
                    write(unit=fd,fmt=*) i+n,j, y(j,1)
                    pocet = pocet + 1
                end if
                if ( y(j,2) /= 0 ) then
                    write(unit=fd,fmt=*) i+n,j+n, y(j,2)
                    pocet = pocet + 1
                end if
            end do
            print *, i+n,pocet
        end do
        close(fd)
    end subroutine printMTX

end module solvers

