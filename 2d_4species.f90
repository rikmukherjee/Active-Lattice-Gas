! Asymmetric Exclusion Process code 
program asymmetric_exclusion_process
    implicit none
  integer,parameter::N=512
  integer,parameter::pe=14
  real,parameter::  lnd = 8 
    integer, parameter :: Nx=N, Ny=N,   EnsMax =1 
    real,parameter:: dt = 0.1 , dtt =  0.1*(N/lnd)**2.0   ! tFinal = Max Time , dtt =  saving Interval
    integer, parameter ::  tMax = 100  ! Loop Index for saving , these many points will be saved
    integer,parameter:: navg = int(dtt/dt)
   
    real(kind=8):: p = 1.0*dt , q = 1.0*pe*(lnd/N)*dt , r = 1.0*(lnd**2/N**2)*dt    
    real, dimension(Nx) :: x
    real, dimension(Ny) :: y
    real, dimension(Nx,Ny) :: pIni
    integer,dimension(Nx,Ny)::u,uIni,up,un
    integer:: i,j, ti,tn , count ,ii
    real(kind=8) :: ran  
    !character*132 :: fname
    !=====================================================================================
    !Norm = p+q+r
    !p = p/Norm
    !q = q/Norm
    !r = r/Norm
    print*,'p=',p,'q=',q,'r=',r,'Pe=',q/(sqrt(p*r)),'L_nd=',N*sqrt(r/p)
    print*,' 4*p - q/2.0 + 3*r=', 4*p - q/2.0 + 3*r
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
    do i=1,Nx
       do j=1,Ny
        pIni(i,j) = exp(-((i-Nx/2)**2.0+(j-Ny/2)**2.0)/(0.62*N)**2.0)
        !pIni(i,:) = 0.5+0.4*exp(-((i-Nx/2)**2)/((0.35*N)**2))


       enddo
    enddo

    ! Step Function
    !pIni = 0
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
            if(ran<0.25               ) uIni(i,j) = 1
            if(ran>0.25 .and. ran<0.50) uIni(i,j)= -1
            if(ran>0.50 .and. ran<0.75) uIni(i,j)=  2
            if(ran>0.75 .and. ran<1.00) uIni(i,j)= -2
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
    count =0
    do i=1,Nx
        do j=1,Nx
            if(uIni(i,j) .ne. 0)count=count + 1
        end do
    end do
    print*,'Fraction=',count/real(N**2)

    

    open(unit=55,file='up.dat',status='unknown')
    open(unit=56,file='un.dat',status='unknown')
    ! Run the simulation
    u=uIni
    ! Time Loop
    do ti =1,tMax
        !--------------------------------------------
        ! do i=1,Nx
        !     do j=1,Ny
        !     if (u(i,j)==1) then
        !       up(i,j)=1 
        !     else
        !       up(i,j)=0
        !     end if
        !     if (u(i,j)==-1) then
        !       un(i,j)=1
        !     else
        !       un(i,j)=0
        !     end if
        !   end do
        ! end do
        !-----------------------------------
        do i=1,Nx
        write(100+ti,*)(u(i,j),j=1,Ny)
        end do
        print*,'Time=',ti
        close(ti+100)
        do tn = 1,navg
            do ii=1,N**2
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
                

                !-----------------Case I ------------------------------
                ! 1 left        -1  right       2   left        -2  right 
                ! Probability p-q/2
                !------------------------------------------------------ 


                if(0<=ran .and. ran < p-q/2.0) then
                    if(u(i,j)== 1 .and. u(modulo(i-2, N) + 1,j)==0)then
                        u(modulo(i-2, N) + 1,j)= 1
                        u(i,j) = 0
                    end if
                    if(u(i,j)==-1 .and. u(modulo(i, N) + 1,j)==0)then
                        u(modulo(i, N) + 1,j) = -1
                        u(i,j) = 0
                    end if
                    if(u(i,j)== 2 .and. u(modulo(i-2, N) + 1,j)==0)then
                        u(modulo(i-2, N) + 1,j) = 2
                        u(i,j) = 0
                    end if

                    if(u(i,j)==-2 .and. u(modulo(i , N) + 1,j)==0)then
                        u(modulo(i , N) + 1,j) = -2
                        u(i,j) = 0
                    end if

                end if
                
                !------------------ Case II ----------------------------------
                ! 1 bottom      -1 top        2 right      -2 left 
                ! Probability p-q/2
                !--------------------------------------------------------------

                if(p-q/2.0<=ran .and. ran < 2*p-q) then
                    if(u(i,j)== 1 .and. u(i,modulo(j-2, N) + 1)==0)then
                        u(i,modulo(j-2, N) + 1)= 1
                        u(i,j) = 0
                    end if
                    if(u(i,j)==-1 .and. u(i,modulo(j, N) + 1)==0)then
                        u(i,modulo(j, N) + 1) = -1
                        u(i,j) = 0
                    end if
                    if(u(i,j)== 2 .and. u(modulo(i, N) + 1,j)==0)then
                        u(modulo(i,N) + 1,j)= 2
                        u(i,j) = 0
                    end if

                    if(u(i,j)==-2 .and. u(modulo(i -2, N) + 1,j)==0)then
                        u(modulo(i -2, N) + 1,j)= -2
                        u(i,j) = 0
                    end if

                end if


                !--------------------- Case III -----------------------------------------
                ! 1 top     -1 bottom       2 bottom        -2 up
                ! Probability p-q/2
                !------------------------------------------------------------------------
                if(2*p-q<=ran .and. ran < 3*p-3*q/2.0) then
                    if(u(i,j)== 1 .and. u(i , modulo(j , N) + 1)==0)then
                        u(i , modulo(j , N) + 1) = 1
                        u(i,j) = 0
                    end if
                    if(u(i,j)==-1 .and. u(i,modulo(j-2, N) + 1)==0)then
                        u(i,modulo(j-2, N) + 1) = -1
                        u(i,j) = 0
                    end if
                    if(u(i,j)== 2 .and. u(i,modulo(j-2, N) + 1)==0)then
                        u(i,modulo(j-2,N) + 1)= 2
                        u(i,j) = 0
                    end if

                    if(u(i,j)==-2 .and. u(i,modulo(j , N) + 1)==0)then
                        u(i , modulo(j , N) + 1)= -2
                        u(i,j) = 0
                    end if

                end if




                !------------------- ---Case IV -----------------------------------------
                ! 1  right      -1 left       2 top          -2 bottom
                ! Probability = p+q/2
                !-----------------------------------------------------------------------
                
                if(3*p-3*q/2.0<=ran .and. ran < 4*p-q) then
                    if(u(i,j)== 1 .and. u(modulo(i, N) + 1,j)==0)then
                        u(modulo(i, N) + 1,j)= 1
                        u(i,j) = 0
                    end if
                    if(u(i,j)==-1 .and. u(modulo(i-2, N) + 1,j)==0)then
                        u(modulo(i-2, N) + 1,j)= -1
                        u(i,j) = 0
                    end if
                    if(u(i,j)== 2 .and. u(i,modulo(j , N) + 1)==0)then
                        u(i,modulo(j, N) + 1) = 2
                        u(i,j) = 0
                    end if
                    if(u(i,j)==-2 .and. u(i,modulo(j-2, N) + 1)==0)then
                        u(i,modulo(j-2, N) + 1)= -2
                        u(i,j) = 0
                    end if
                end if
                
                
                
            
            
            
                !if( 4.0*p - q  + 2.0*r <= ran .and. ran < 4.0*p - q + 3.0*r  ) u(i,j) = -u(i,j)
                

                ! New Piece of code :
                if(4.0*p - q         <= ran .and. ran < 4.0*p - q + r    .and. u(i,j)== 1 ) u(i,j) = -1
                if(4.0*p - q + r     <= ran .and. ran < 4.0*p - q + 2*r  .and. u(i,j)== 1 ) u(i,j) =  2
                if(4.0*p - q + 2*r   <= ran .and. ran < 4.0*p - q + 3*r  .and. u(i,j)== 1 ) u(i,j) = -2

                if(4.0*p - q + 3*r   <= ran .and. ran < 4.0*p - q + 4*r  .and. u(i,j)==-1 ) u(i,j) =  1
                if(4.0*p - q + 4*r   <= ran .and. ran < 4.0*p - q + 5*r  .and. u(i,j)==-1 ) u(i,j) =  2
                if(4.0*p - q + 5*r   <= ran .and. ran < 4.0*p - q + 6*r  .and. u(i,j)==-1 ) u(i,j) = -2
                
                if(4.0*p - q + 6*r   <= ran .and. ran < 4.0*p - q + 7*r  .and. u(i,j)== 2 ) u(i,j) =  1
                if(4.0*p - q + 7*r   <= ran .and. ran < 4.0*p - q + 8*r  .and. u(i,j)== 2 ) u(i,j) = -1 
                if(4.0*p - q + 8*r   <= ran .and. ran < 4.0*p - q + 9*r  .and. u(i,j)== 2 ) u(i,j) = -2 
                
                if(4.0*p - q + 9*r   <= ran .and. ran < 4.0*p - q + 10*r .and. u(i,j)==-2 ) u(i,j) =  1 
                if(4.0*p - q + 10*r  <= ran .and. ran < 4.0*p - q + 11*r .and. u(i,j)==-2 ) u(i,j) = -1 
                if(4.0*p - q + 11*r  <= ran .and. ran < 4.0*p - q + 12*r .and. u(i,j)==-2 ) u(i,j) =  2 



                



            end do
        end do


    end do
   

end program asymmetric_exclusion_process
