%% Burgers
%---------------------------------------
% This code solves hydrodynamics for Quasi-1d Lattice Gas model
% Pseudospectral method. Time marching in RK4 with
% linear term explicitly integrated.
%---------------------------------------
clear ; clc ; close
Pe = 10 ;
%----------- Domain ----------------------
L		= 	8                ;
N		= 	2^12             ;
dx 		= 	L/N              ;
x 		= 	(0:dx:L-dx)'     ;
%---------- Parameters -------------------
nu   	=  1           		;  % Diffusion = 1 
rho_0  	= .40          		;  % Mean density
%-------------Time------------------------------
tMax 	= 10             	;
dt   	= 10.0^-2       	;
dT   	= .1             	;
%------------------------------------------------

t  = 0:dt:tMax-dt       ;
T  = 0:dT:tMax          ;  
interval = fix(dT/dt)   ;
TStep = length(T)       ;
tStep = length(t)       ;
%------------- Wave Vector--------------------
k = 2*pi/L*(-N/2:N/2-1) ;
k = fftshift(k');
%--------------------------------------


%--------------Initial condition------------------
%------------- Choose Initial Condition here---------------------

%---------------Half Gaussian------------------
u0 = zeros(N,1);
v0 = zeros(N,1);

%------------ Both gaussian--------------------
%u0 =0
sigma = L/2;
u0 = 0.2*exp(-(x-L/2).^2/sigma^2) + 0.0*rand(N,1); 
v0 = 0.02*sin(x/L*(2*pi)) ;

%========================
%nn = int32(N/8);
%u0(N/2-nn:N/2+nn-1) = importdata('uIni.txt');
%v0(N/2-nn:N/2+nn-1) = importdata("vIni.txt");
%u0(1:N/4-1) = 0.1230*ones(N/4-1,1);
%u0(N/2+N/4:end) = 0.1230*ones(N/4+1,1);



%=========================================================================
%------------ Proper Mean Density-----------------------------------------
mean_den = mean(u0)    ;    
u0 = rho_0/mean_den*u0 ;
%v0 = rho_0/mean_den*v0;
% --------------Plotting Initial Profile--------------------------------
%figure()
plot(x,u0,LineWidth=1.5); hold on      ;   
plot(x,v0,LineWidth=1.5); hold on      ;
%plot(x,u0_1+v0_1,LineWidth=1.5,Color='k')  ;
legend('$\rho$','$m$',interpreter='latex',fontsize=14);
%------------------------------------------------------------------------
% Creating array for time evolution
w0 = cat(1,fft(u0),fft(v0)) ;
% ---------------Initialization of Arrays ----------------------------
uu_hat = zeros(N,TStep); uu = zeros(N,TStep);
vv_hat = zeros(N,TStep); vv = zeros(N,TStep);
%--------------------------------------------------------------------
uu(:,1) = u0 ; 
vv(:,1) = v0 ;
%-----------------------------------
w = w0;
%% ==================== Main Loop =======================================
for i = 1:tStep
    % --------- Writing ---------------
    if mod(i,interval)==0 ; j = fix(i/interval)+1 ;
       uu_hat(:,j) = w(1:N) ;
       uu(:,j) = real(ifft((w(1:N,1)))) ;
       vv_hat(:,j) = w(N+1:2*N,1) ;
       vv(:,j) = real(ifft((w(N+1:2*N,1)))) ;
    end
    w = RK4(w,dt,k,nu,Pe,N);
    %f_k = cat(1,fft(randn(N,1)) , randn(N,1))   ; 
    %w = w + .1*sqrt(dt)*f_k ;
end 
%%
writematrix(uu,'rho_MetaStable.dat')
writematrix(vv,'m_MetaStable.dat')

%% Function
function dw_dt = rhs(w,k,nu,Pe,N)
u = w(1:N)     ;
v = w(N+1:2*N) ;
% Here since D=1 we don't need dealiasing if we use sufficient 
% grid points : If D != 1 use the 2/3 dealiasing function defined 
% at the end 
d_dx_u  = 1i*k.*fft(real(ifft(v)).*(1-real(ifft(u)))) ;  
d_dx_v  = 1i*k.*fft(real(ifft(u)).*(1-real(ifft(u)))) ;
du2_dx2 = real(ifft(-k.^2.*u)) ;
dv2_dx2 = real(ifft(-k.^2.*v)) ; 
u_real  = real(ifft(u))        ;
v_real  = real(ifft(v))        ;
du_dt =  -Pe*d_dx_u                                                       ;
dv_dt =  -Pe*d_dx_v  - fft(u_real.*dv2_dx2) + fft(v_real.*du2_dx2) - 2*v  ;   
dw_dt = cat(1,du_dt , dv_dt)                                              ;
size(dw_dt);
end 


 %--------------------------------------------------
  function rk4 = RK4(u_hat,dt,k,nu,Pe,N)
    nu1 = nu ; nu2 = nu;
    semi_g_half1 = exp(-nu1*k.^2*dt/2) ; 
    semi_g_half2 = exp(-nu2*k.^2*dt/2) ;
    semi_g_half = cat(1,semi_g_half1,semi_g_half2);
    semi_g1 = exp(-nu1*k.^2*dt) ;
    semi_g2 = exp(-nu2*k.^2*dt) ;
    semi_g = cat(1,semi_g1,semi_g2);
    u_hat_old = u_hat ;  
    k1 = rhs(u_hat_old,k,nu,Pe,N) ; 
    u_hat = semi_g_half.*(u_hat + .5*dt*k1) ;
    k2 = rhs(u_hat,k,nu,Pe,N);
    u_hat = semi_g_half.*(u_hat_old) + .5*dt*k2 ;
    k3 = rhs(u_hat,k,nu,Pe,N);
    u_hat = semi_g.*(u_hat_old) + semi_g_half.*k3*dt ;
    k4 = rhs(u_hat,k,nu,Pe,N);
    u_hat = semi_g.*u_hat_old + dt*(semi_g.*k1 + 2*semi_g_half.*k2 + 2*semi_g_half.*k3 + k4)/6 ;
    rk4 = u_hat ; % This is finally returned
  end
  
  
function prod = de_2_3(u,v)
N = length(u);
u(N/2-ceil(N/6):N/2+ceil(N/6))=0;
v(N/2-ceil(N/6):N/2+ceil(N/6))=0;
prod = fft(real(ifft(u)).*real(ifft(v)));
prod(N/2-ceil(N/6):N/2+ceil(N/6))=0;
end
