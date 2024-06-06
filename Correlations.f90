! Asymmetric Exclusion Process code 
program asymmetric_exclusion_process
    implicit none
  integer,parameter::N=2**8
  real,parameter:: pe = 1.0 , lnd = 1 
    integer,parameter:: Nx=N, Ny=1,   EnsMax =1
    real,parameter:: dt = 0.05 , dtt =  0.0001*(N/lnd)**2.0  ! tFinal = Max Time , dtt =  saving Interval
    integer, parameter ::   tMax = 50000    ! Outer Loop Index for saving , these many points will be saved
    integer,parameter:: navg = int(2*dtt/dt)
    real:: p = 0.5*dt , q = 100*0.5*pe*(lnd/N)*dt , r = 10000*0.5*(lnd**2/N**2)*dt , epsa = 5.0*dt
    real, dimension(Nx) :: x
    real, dimension(Ny) :: y
    real, dimension(Nx,Ny) :: pIni
    integer,dimension(Nx,Ny)::u,uIni,up,un
    integer:: i,j, ti,tn, tj , ii ,  count = 0
    real :: ran  
    integer,parameter::rMax=N/4,tgap=500
    real,dimension(rMax)::cor,corRhoM,AvgCor=0,AvgCorRhoM=0 
    real,dimension(N) :: AvgRho_0 = 0
    character*132 :: fname
    !=====================================================================================
    !Norm = p+q+r
    !p = p/Norm
    !q = q/Norm
    !r = r/Norm
    print*,'p=',p,'q=',q,'r=',r,'Pe=',q/(sqrt(p*r)),'L_nd=',N*sqrt(r/p)
    print*,'2*p+q+r=',2*p+q+r+2*epsa
    print*,'navg=',navg 
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
    pIni(i,j) = 0.7*exp(-((i-Nx/2)**2+(j-Ny/2)**2)/((0.1*N)**2.0))

        enddo
    enddo

    ! Step Function
    !pIni = 0
    !pIni(int(N/8):int(N),int(N/8):int(N))=1
    ! Sinusoidal
    !pIni = 0.2 + 0.8*sin(pi*y/N)
    ! Generate initial configuration uIni
    uIni(:,:)=0
    do i = 1, Nx
        do j=1,Ny
        call random_number(ran)
        if (ran < pIni(i,j).and.i<int(Nx/2)) then
        uIni(i,j) = 1.0
        end if
        if (ran < pIni(i,j).and.i>=int(Nx/2)) then
        uIni(i,j) = -1.0
        end if
        end do
    end do
   
    ! Output initial configuration
    open(1, file="initial.dat")
    do i = 1, Nx
        write(1,*) (uIni(i,j),j=1,Ny)
    end do
    close(1)

    

    open(unit=55,file='up.dat',status='unknown')
    open(unit=56,file='un.dat',status='unknown')
    ! Run the simulation
    u=uIni
    ! Time Loop
    do ti =1,tMax

        if(mod(ti,tgap)==0)print*,'Time=',ti ,'  ' , 'NumParticle=' ,sum(abs(u))

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
        !!======================================================
              !=============================================
        if(ti>200)then 
            do ii=1,rMax
            cor(ii)    = sum((up(:,1)+un(:,1))*cshift((up(:,1)+un(:,1)),ii-rMax/2))/real(N)
            corRhoM(ii)= sum((up(:,1)+un(:,1))*cshift((up(:,1)-un(:,1)),ii-rMax/2))/real(N)
            
            end do
            !========================
            AvgCor     = AvgCor     + cor
            AvgCorRhoM = AvgCorRhoM + corRhoM
            AvgRho_0   = AvgRho_0   + up(:,1)+un(:,1)

            count  = count  + 1 
        end if
      !=============================================
       
       
        do tn=1,navg
            do tj=1,N
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
                if(2*p<=ran .and. ran < 2*p + epsa) then
                    if(u(i,j) .ne. 0 .and. u(i,modulo(j-2,Ny)+1)==0)then
                        u(i,modulo(j-2, Ny) + 1)= u(i,j)
                        u(i,j) = 0
                    end if
                end if
                if(2*p+epsa<=ran .and. ran < 2*p+2*epsa) then
                    if(u(i,j) .ne. 0 .and. u(i,modulo(j, Ny) + 1)==0)then
                        u(i,modulo(j, Ny) + 1)= u(i,j)
                        u(i,j) = 0
                    end if
                end if
            
            if(2*p+2*epsa <=ran .and. ran<2*p+2*epsa+r)then
                u(i,j) = -u(i,j)
            end if

        enddo


        end do


    end do
    open(unit=102,file='Corr.dat',status='unknown')
    open(unit=103,file='CorrRhoM.dat',status='unknown')
  
        print*,'Avg_Rho_0=',sum(AvgRho_0)/real(count*N)
        
        do ii=1,rMax
         if(ii .ne. rMax/2)then
         write(102,*)ii/real(N),AvgCor(ii)/count !- AvgCor(rMax)/count ! sum(abs(u))**2.0/real(N**2.0) !
         end if
        end do
        close(102)
  
        do ii=1,rMax
          if(ii .ne. rMax/2)then
          write(103,*)ii/real(N),AvgCorRhoM(ii)/count !- AvgCorRhoM(rMax)/count!-(sum(abs(u))/real(N))**2.0
        end if
        end do
         close(103)
   

end program asymmetric_exclusion_process
