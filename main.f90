!> \file main.f90
!! \brief hlavni program
!<



!> hlavni programek
!!
!! \mainpage Nastroje pro pocitani cidla
!!  obsahuje konstrukci modelu
!!  a jeho iterativni upravu pro vypocet odezvy na elektrodach.
!! Cilem je odhadnout tvar pozorovaneho objektu
!! - \subpage model "Model"
!! - \subpage diskr "Diskretizace"
!! - \subpage implem "Implementace matic"
!!  
!<
program jakub
    use pmatypy
    use local
    use geometry
    use solvers_jakub
    use kreslitko
    use mtx_jakub
    use pmatools
    use global
    use mtx
    use sparsematrix
    implicit none

    !!! vsude predpoklad, ze vektor(0,:)=0 pro funkci axpy
    real(kind=rkind), dimension(:,:), allocatable :: x,y,z
    real(kind=rkind), dimension(:), allocatable :: x1,y1,z1
    integer(kind=ikind) :: np,i,cnt, n,m
    integer(kind=ikind) :: ii,jj,gi,gj
    real(kind=rkind) :: val
    character :: zn
    character(*), parameter :: jmeno1=&
!        "/home/pmayer/Ubuntu One/triangulace_site/model2/generator_geometrie"
        "C:\Users\pmayer\SkyDrive\test\fortran\triangulace_site\model2\generator_geometrie"
    type(Uloha) :: U
    type(smtx)  :: A

    ! napred otestuj matice
    call plotstart()



    print *,"cau v"
    call ReadAgros(jmeno1,U)   !precte data z Agrosu
    print *, "pocet bodu", U%np
    print *, "pocet trojuhelniku", U%nt
    print *, "pocet hran", U%nb

    ! ted vykreslit jednotlive oblasti
    call plotriang(U)
    print *,"triangulace nakreslena"
    call plotbound(U)
    print *,"hranice nakreslena"
    call plotdomain(U)
    print *,"domeny nakresleny"

    call CreateLocal(U) !vyrobi lokalni matice


    !A ted udelam assembly
    call A%init(U%np,U%np)
    do i=1,U%nt
        do ii=1,3
            do jj=1,3
                val = U%LocList(i-1)%A(ii,jj)
                gi = U%LocList(i-1)%map(ii,1) + 1
                gj = U%LocList(i-1)%map(jj,1) + 1
                call A%add(val,gi,gj)
            end do
        end do
    end do
    call A%spy
    n = A%getn()
    m = A%getm()
    allocate(x1(n),y1(m))
    x1 = 1
    y1 = A%mul(x1)
    call pockej("po nasobeni")
    do i = 1,n
        print *,i,x1(i),y1(i)
    end do

    !   do ii = lbound(U%data1,1), ubound(U%data1,1)
    !      print *, "label=",ii
    !      call plotdomaini(U,ii,7,.false.)
    !      !read(*,*) zn
    !   end do
    !
    !   do ii = lbound(data2,1), ubound(data2,1)
    !      print *, "label=",ii
    !      call plotboundi(U,ii,8)
    !      !read(*,*) zn
    !   end do

    !   !jde se pocitat
    !   np = ubound(U%PointList,1)+1
    !   print *,"pocet bodu=",np
    !   allocate(x(0:np,1:2))
    !   allocate(y(0:np,1:2))
    !   allocate(z(0:np,1:2))
    !   x=0
    !   y=0
    !   z=0
    !   print *,np
    !   call GetRhs(U,z)
    !   call Axpy(U,x,y)
    !   print *,"nula?", Dot_Product(y(:,1),y(:,1))+Dot_Product(y(:,2),y(:,2))
    !   x = 1
    !   x(0,:)=0
    !   call Axpy(U,x,y)
    !   !z = y  ! jen pro otestovani, ze to funguje
    !   cnt = 0
    !   do i=0,np
    !      cnt = cnt + 1
    !      print *, i,y(i,:)
    !      if ( cnt == 20) then
    !         cnt = 0
    !         !read(*,*) zn
    !      end if
    !   end do
    !
    !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   ! je to moc zdlouhave - jedine az s pomoci openmp
    !   call printMTX(U)
    !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !   !! jdeme to spocitat
    !   x = x*0.2
    !   !call SDS(x,z)
    !   call CGS(U,x,z)
    !   call SDS(U,x,z)
    !   print *, "neco zmackni a enter"
    !  read(*,*) zn


    call plotstop()

end program jakub
