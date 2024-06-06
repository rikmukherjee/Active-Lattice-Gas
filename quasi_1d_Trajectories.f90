! Asymmetric Exclusion Process code 
module mod_wasep
   implicit none
  integer,parameter::N=256
  real,parameter::pe=10
  real,parameter:: lnd = 6 
  integer,parameter::Nx = N , Ny = 2
    ! tFinal = Max Time , dtt =  saving Interval
    real,parameter:: dt = 0.1 , dtt =   1  , tFinal = 200*N  
    ! Loop Index for saving , these many points will be saved
    integer, parameter ::  tMax = int(tFinal/dtt)           
    integer,parameter:: navg = int(2*dtt/dt)
    !!************************************************************************************
    real:: p = 1.0*dt , q = 1.0*pe*(lnd/N)*dt , r = 1*(lnd**2/N**2)*dt , epsa = 3.0*dt/N**0
    real, dimension(Nx) :: x
    real, dimension(Ny) :: y
    real, dimension(Nx,Ny) :: pIni
    integer,dimension(Nx,Ny)::u,uIni,up,un,udummy
    integer,dimension(:,:),allocatable::loc
    integer,dimension(:,:,:),allocatable::loctime
    integer,dimension(Nx,Ny,tMax)::statestime
    integer,dimension(:),allocatable::one_d_loc
    integer:: i,j, ii,ti,tn ,tii, site ,NumParticle,ix
    integer,dimension(1):: locindex
    real :: ran  
    character*132 :: fname
end module 




program asymmetric_exclusion_process
    use mod_wasep
    !=====================================================================================!
    print*,'p=',p,'q=',q,'r=',r,'Pe=',q/(sqrt(p*r)),'L_nd=',N*sqrt(r/p)
    print*,'2*p+2*epsa+q+r=',2*p+2*epsa+r
    print*,'tMax=',tMax
    !-------------------------------------------------------------------------------------!
    !=====================================================================================!
    ! Initialize arrays
    x = [(i, i=1,Nx)]
    y = [(i, i=1,Ny)]
    !======================================================================================!
    ! Choosing Initial Condition 
    !======================================================================================!

    ! Step Function
    
    ! call stepIni

    ! Gaussian
    
    call GaussianIni
    
    ! Now Generate initial configuration From the probablity distribution

    NumParticle=0
    uIni(:,:)=0

    call ConfigIni
    !=======================================================================================!
    ! Initial Configuration Done
    !=======================================================================================!

    print*,'NumParticle = ',NumParticle,'Fraction=',0.5*NumParticle/real(N)
    allocate(loc(NumParticle,2))
    allocate(loctime(NumParticle,2,tMax))

    call IniLocation
  
    ! Output initial configuration
    call writeInitials
    !---------------------------------------------------------------------------------------!
    open(20,form='unformatted' ,file='Trajectories.bin', status='unknown')
    open(21,form='unformatted' ,file='States.bin'      , status='unknown')
    !!======================================================================================!

    
    !!======================================================================================!
    !                  Main Simulation Starts
    !=======================================================================================!
    u=uIni
    ! Time Loop
    do ti =1,tMax
        !-------------------- Write Results-------------------------------------------------!
        !call write_up_un
        !write(20,*)(loc(ii,:),ii=1,NumParticle)
        !write(21,*)(abs(u(ii,1)),ii=1,Nx)
        loctime(:,:,ti) = loc
        statestime(:,:,ti) = u
        !print*,loc(:,2)
        !-----------------------------------------------------------------------------------!
        !print*,'Time=',ti

        do tn=1,navg


            do tii=1,Nx
        
                ! Randomly selcting a site
                call random_number(ran)
                i = 1 + int(ran*Nx)
                call random_number(ran)
                j = 1 + int(ran*Ny)
                ! Site selected
                !-----------------------------------
                call random_number(ran)
                !-----------------------------------
                ! Diffusion + Activity
                !-----------------------------------
                if(0<=ran .and. ran < p-q/2.0) then
                    
                    if(u(i,j)==1 .and. u(modulo(i-2, N) + 1,j)==0)then

                        locindex        =   pack([(ix,ix=1,size(loc(:,1)))],loc(:,2)==j.and.loc(:,1)==i)
                        !locindex        =   findloc(loc,i)
                        site            =   modulo(i-2, N) + 1
                        u(site,j)       =   1
                        u(i,j)          =   0
                        loc(locindex,1) =   site    ! update location
                    end if


                    if(u(i,j)==-1 .and. u(modulo(i, N) + 1,j)==0)then
                        site            =   modulo(i, N) + 1
                        locindex        =   pack([(ix,ix=1,size(loc(:,1)))],loc(:,2)==j.and.loc(:,1)==i)
                        !locindex        =   findloc(loc,i)
                        u(site,j)       =   -1
                        u(i,j)          =   0
                        loc(locindex,1) =   site    
                    end if


                end if


                if(p-q/2.0<=ran .and. ran < 2*p) then


                    if(u(i,j)==1 .and. u(modulo(i, N) + 1,j)==0)then
                        locindex        =   pack([(ix,ix=1,size(loc(:,1)))],loc(:,2)==j.and.loc(:,1)==i)
                        !locindex        =   findloc(loc,i)
                        site            =   modulo(i, N) + 1
                        u(site,j)       =   1
                        u(i,j)          =   0
                        loc(locindex,1)   =  site    ! update location
                    end if


                    if(u(i,j)==-1 .and. u(modulo(i-2, N) + 1,j)==0)then
                        locindex        =   pack([(ix,ix=1,size(loc(:,1)))],loc(:,2)==j.and.loc(:,1)==i)
                        !locindex        =   findloc(loc,i)
                        site            =   modulo(i-2, N) + 1
                        u(site,j)       =   -1
                        u(i,j)          =    0
                        loc(locindex,1) =   site      ! update location
                    end if


                end if


                if(2*p<=ran .and. ran < 2*p + epsa) then
                    
                    
                    if(u(i,j) .ne. 0 .and. u(i,modulo(j-2,Ny)+1)==0)then
                        locindex        =   pack([(ix,ix=1,size(loc(:,1)))],loc(:,2)==j.and.loc(:,1)==i)
                        site            =   modulo(j-2, Ny) + 1
                        u(i,site)       =   u(i,j)
                        u(i,j)          =   0
                        loc(locindex,2) =   site
                    end if


                end if


                if(2*p+epsa<=ran .and. ran < 2*p+2*epsa) then

                    if(u(i,j) .ne. 0 .and. u(i,modulo(j, Ny) + 1)==0)then
                        locindex        =   pack([(ix,ix=1,size(loc(:,1)))],loc(:,2)==j.and.loc(:,1)==i)
                        site            =   modulo(j, Ny) + 1
                        u(i,site)       =   u(i,j)
                        u(i,j)          =   0
                        loc(locindex,2) =   site

                    end if

                end if
                
                ! Tumbling
                if(2*p+2*epsa <=ran .and. ran<2*p+2*epsa+r)then

                    u(i,j) = -u(i,j) 
                
                end if
            enddo
            

        end do
           


    end do
   
   write(20)loctime
   write(21)statestime
   close(20)
   close(21)

end program asymmetric_exclusion_process

!subroutine swap(a,b)
!    real::a,b,temp
!    temp =  a 
!    a    =  b
!    b    =  temp

 
subroutine stepIni()
    use mod_wasep
    pIni = 0
    pIni(int(N/2-N/8):int(N/2+N/8),1)=1
end subroutine stepIni

subroutine GaussianIni()
    use mod_wasep
    real::sigma=0.35*N
    pIni = 0
    do i=1,Nx
        do j=1,Ny
            pIni(i,j) = 0.5+0.4*exp(-((i-Nx/2)**2+(j-int(Ny/2.0))**2)/(sigma**2))
        enddo
    enddo
end subroutine GaussianIni


subroutine ConfigIni()
    use mod_wasep
    do i = 1, Nx
        do j=1,Ny
        call random_number(ran)
        if (ran < pIni(i,j).and.i<int(Nx/2)) then
        uIni(i,j) = 1.0
        NumParticle=NumParticle+1
        end if
        if (ran < pIni(i,j).and.i>=int(Nx/2)) then
        uIni(i,j) = -1.0
        NumParticle=NumParticle+1
        end if
        end do
    end do
end subroutine ConfigIni

    
subroutine IniLocation()
    use mod_wasep
    locindex=1  
    do i=1,Nx 
        do j=1,Ny
            if(uIni(i,j) .ne. 0)then
                loc(locindex,1)=i
                loc(locindex,2)=j 
                locindex=locindex+1! Initialise location
            !print*,loc(i),uIni(i,1)
            endif
        enddo
    enddo
end subroutine IniLocation


subroutine writeInitials()
    use mod_wasep
    open(1, file="initial.dat")
    do i = 1, Nx
        write(1,*) (uIni(i,j),j=1,Ny)
    end do
    close(1)
    open(1, file="initial-location.dat",status='unknown')
    do i=1,NumParticle
        write(1,*) loc(i,:)
    enddo
    close(1)
end subroutine writeInitials

    
subroutine write_up_un()
    use mod_wasep
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

    
    !open(100+ti,status='unknown')
    !    do i=1,Nx
    !        write(100+ti,*)(u(i,j),up(i,j),un(i,j),j=1,Ny)
    !    end do
    !close(100+ti)

end subroutine write_up_un


