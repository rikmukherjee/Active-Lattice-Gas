! Asymmetric Exclusion Process code 
program asymmetric_exclusion_process
    implicit none
  integer,parameter::N=512
  integer,parameter::pe=12
    integer, parameter :: Nx=N, Ny=N,   EnsMax =1 
    real,parameter:: dt = 0.1 , dtt =  N , tFinal = 100*N   ! tFinal = Max Time , dtt =  saving Interval
    integer, parameter ::  tMax = int(tFinal/dtt)+1  ! Loop Index for saving , these many points will be saved
    integer,parameter:: navg = int(N**2/dt)*dtt
    real,parameter::  lnd = 4 
    real(kind=8):: p = 2.0*dt , q = 2.0*pe*(lnd/N)*dt , r = 2.0*(lnd**2.0/N**2.0)*dt    
    real, dimension(Nx) :: x
    real, dimension(Ny) :: y
    real, dimension(Nx,Ny) :: pIni
    integer,dimension(Nx,Ny)::u,uIni,up,un
    integer:: i,j, ti,tn
    real(kind=8) :: ran  
    character*132 :: fname
    !=====================================================================================
    !Norm = p+q+r
    !p = p/Norm
    !q = q/Norm
    !r = r/Norm
    print*,'p=',p,'q=',q,'r=',r,'Pe=',q/(sqrt(p*r)),'L_nd=',N*sqrt(r/p)
    print*,'4*p+q+r=',4*p+r
    print*,'tMAx=',tMax
    !-------------------------------------------------------------------------------------
    !=====================================================================================
    ! Initialize arrays
    ! Choosing Initial Condition
    ! Gaussian
    x = [(i, i=1,Nx)]
    y = [(i, i=1,Ny)]
    !x2 = spread(x, 1, size(y))
    !y2 = spread(y, 2, size(x))
    pIni = 0
    do i=1,Nx
        do j=1,Ny
         pIni(i,j) = exp(-((i-Nx/2)**2.0+(j-Ny/2)**2.0)/(0.6*N)**2.0)
         !pIni(i,:) = 0.5+0.4*exp(-((i-Nx/2)**2)/((0.35*N)**2))
 
 
        enddo
     enddo
 

    ! Step Function
    !
    !pIni(int(N/2-3*N/8):int(N/2+4.0*N/8),int(N/2-3*N/8):int(N/2+3.5*N/8))=1
    ! Sinusoidal
    !pIni = 0.2 + 0.8*sin(pi*y/N)
    ! Generate initial configuration uIni
    uIni(:,:)=0
    do i = 1, Nx
        do j=1,Ny
        call random_number(ran)
        if (ran < pIni(i,j) ) then
        call random_number(ran)
            if(ran < 0.5 )then
            uIni(i,j) = 1.0
            else
                uIni(i,j)=-1
        end if
        end if
        end do
    end do
    ! do i = 1, Nx
    !     do j=1,Ny
    !     call random_number(ran)
    !     if (ran < pIni(i,j).and.i<int(Nx/2)) then
    !     uIni(i,j) = 1.0
    !     end if
    !     if (ran < pIni(i,j).and.i>=int(Nx/2)) then
    !     uIni(i,j) = -1.0
    !     end if
    !     end do
    ! end do
   
    ! Output initial configuration
    open(1, file="initial.dat")
    do i = 1, Nx
        write(1,*) (uIni(i,j),j=1,Ny)
    end do
    close(1)
    print*,'Fraction=',sum(abs(real(uIni)))/real(N**2)

    

    open(unit=55,file='up.dat',status='unknown')
    open(unit=56,file='un.dat',status='unknown')
    ! Run the simulation
    u=uIni
    ! Time Loop
    do ti =1,tMax
        !--------------------------------------------
        do i=1,Nx
            do j=1,Ny
            if (u(i,j)==1) then
              up(i,j)=1 
            else
              up(i,j)=0
            end if
            if (u(i,j)==-1) then
              un(i,j)=1
            else
              un(i,j)=0
            end if
          end do
        end do
        !-----------------------------------
        do i=1,Nx
        write(100+ti,*)(u(i,j),j=1,Ny)
        end do
        !-----------------------------------------
        print*,'Time=',ti
        close(ti+100)
        do tn=1,navg
            ! Randomly selcting a site
            call random_number(ran)
            i = 1 + int(ran*Nx)
            call random_number(ran)
            j = 1 + int(ran*Ny)
            ! Site selected
            !---------------------------------
            call random_number(ran)
            !----------------------------------
            ! Diffusion + Activity
            !---------------------------------

            if(0<=ran .and. ran < p-q/2.0) then
                if(u(i,j)==1 .and. u(modulo(i-2, N) + 1,j)==0)then
                    u(modulo(i-2, N) + 1,j)= 1
                    u(i,j) = 0
                end if
                if(u(i,j)==-1 .and. u(modulo(i, N) + 1,j)==0)then
                    u(modulo(i, N) + 1,j)= -1
                    u(i,j) = 0
                end if
            end if
            
            
            if(p-q/2.0<=ran .and. ran < 2*p) then
                if(u(i,j)==1 .and. u(modulo(i, N) + 1,j)==0)then
                    u(modulo(i, N) + 1,j)= 1
                    u(i,j) = 0
                end if
                if(u(i,j)==-1 .and. u(modulo(i-2, N) + 1,j)==0)then
                    u(modulo(i-2, N) + 1,j)= -1
                    u(i,j) = 0
                end if
            end if
            
            
            if(2*p<=ran .and. ran < 3*p) then
                if(u(i,j) .ne. 0 .and. u(i,modulo(j-2,Ny)+1)==0)then
                    u(i,modulo(j-2, Ny) + 1)= u(i,j)
                    u(i,j) = 0
                end if
            end if
            
            
            if(3*p<=ran .and. ran < 4*p) then
                if(u(i,j) .ne. 0 .and. u(i,modulo(j, Ny) + 1)==0)then
                    u(i,modulo(j, Ny) + 1)= u(i,j)
                    u(i,j) = 0
                end if
            end if








                !if(u(i,j)==1 .and. u(i,modulo(j, N) + 1)==0)then
                !    u(i,modulo(j, N) + 1)= 1
                !    u(i,j) = 0
                !end if
                !if(u(i,j)==-1 .and. u(i,modulo(j-2, N) + 1)==0)then
                !    u(i,modulo(j-2, N) + 1)= -1
                !    u(i,j) = 0
                !end if
!
            
            ! Tumbling
            if(4*p<=ran .and. ran<4*p+r)then
                u(i,j) = -u(i,j)
            end if


        end do


    end do
   

end program asymmetric_exclusion_process
