%EM cloaking FDTD 2D, diagonalized constitutive parameter tensor
%UPML, TFSF, no loss, planewave
%Oliver Csernyava BME Project Laboratory 1. \mail: oliver.csernyava@sch.bme.hu

% All rights reserved

x = (1:n_x)*dx; %spatial vectors for plotting
z = (1:n_z)*dz;
%_______________________________________________________________________

% TIME _ COURANT LIMIT__________________________________________________
lambda_min = min(min(min(la1(la1>0))),min(min(eps_yy(eps_yy>0))));
dt = 0.9*dx*lambda_min / (c*sqrt(2));
%_______________________________________________________________________

%scale factor for the expected max value of the E -> for the PLOTting range
Q = dx/dt/c/2; %% SET the max value of E to the timestep

%Iteration No.----------------------------- 
    n_t = N_pass*(n_x-2*D_pml)*dx/dt/c; %number of total iterations 
%------------------------------------------

%_______________________________________________________________________
%UPML___________________________________________________________________
D_pml = ceil(n_x*q_pml*0.01); %number of the spatial points in the PML (width)
if D_pml < 15 %experienced minimum value for the PML width
    msg = 'D_pml should be greater than 15';
    error(msg) 
end
eta = sqrt( mu_0*mu_r/(eps_0*eps_r));
sig_max = -(m+1)*log(Rerr)/(2*eta*D_pml*dx);
sig_x = zeros(1,n_x);
sig_pmlE = ((D_pml-(1:D_pml))/(D_pml)).^m * sig_max;
sig_x(1:D_pml) = sig_pmlE; %set the parameters for the left side
for i = n_x-D_pml+1:n_x
    sig_x(i) = sig_pmlE(n_x-i+1); %set the parameters for the right side
end
                           %same parameters for the E and H fields as well!
sig_z = sig_x;
kap_x = ones(1,n_x);

kapp = (1+(kapp_max-1)*((D_pml-(1:D_pml))/(D_pml)).^m);
kap_x(1:D_pml) = kapp;
for i = n_x-D_pml+1:n_x
    kap_x(i) = kapp(n_x-i+1);
end
kap_z = kap_x;
%_______________________________________________________________________

%Auxiliary field coefficients(for the planewave source)-----------------
Ey_aux = zeros(1,n_x); %the electric field auxiliary 1 for 1Dim
Hz_aux = zeros(1,n_x-1); %the magnetic field auxiliary
K_M1 = (c*dt/sqrt(eps_r)-dx)/(c*dt/sqrt(eps_r)+dx); % Mur
K_M2 = 2*dx/(c*dt/sqrt(eps_r)+dx);
Ey_n = zeros(1,4);
K_E1 = dt/(eps_0*eps_r*dx);
K_M = dt/(mu_0*mu_r*dx);
%-----------------------------------------------------------------------

%Field vector coefficients initialization ------------------------------
Ey = zeros(n_x,n_z);
Hx = zeros(n_x-1,n_z-1);
Hz = zeros(n_x-1,n_z-1);
Dy = zeros(n_x,n_z);
Bx = zeros(n_x-1,n_z-1);
Bz = zeros(n_x-1,n_z-1);
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
%Parameters for the PML + Homogenous media(no loss, no dispersion)------
K_Dy1 = zeros(1,n_x);
K_Dy2 = zeros(1,n_x);
K_Dy3 = zeros(1,n_x);

K_Bz1 = zeros(1,n_x);
K_Bz2 = zeros(1,n_x);

K_Hx2 = zeros(1,n_x);
K_Hx3 = zeros(1,n_x);

K_Bx1 = zeros(1,n_z);
K_Bx2 = zeros(1,n_z);

K_Ey1 = zeros(1,n_z);
K_Ey2 = zeros(1,n_z);
K_Ey3 = zeros(1,n_z);
K_Ey12 = zeros(1,n_z);
K_Ey22 = zeros(1,n_z);
K_Ey32 = zeros(1,n_z);

K_Hz1 = zeros(1,n_z);
K_Hz2 = zeros(1,n_z);
K_Hz3 = zeros(1,n_z);

C_Eyz = dt/(eps_0*eps_r*dx);
C_Eyx = dt/(eps_0*eps_r*dz);
C_Hz = dt/(mu_0*mu_r*dx);
C_Hx = dt/(mu_0*mu_r*dz);

for i = 1:n_x 
K_Dy1(i) = -(dt*sig_x(i) - 2*kap_x(i)*eps_0)/(dt*sig_x(i) + 2*kap_x(i)*eps_0);
K_Dy2(i) = (2*dt*eps_0)/(dt*sig_x(i) + kap_x(i)*2*eps_0)/dx;
K_Dy3(i) = (2*dt*eps_0)/(dt*sig_x(i) + kap_x(i)*2*eps_0)/dz;

K_Bz1(i) = -(dt*sig_x(i) - 2*kap_x(i)*eps_0)/(dt*sig_x(i) + 2*kap_x(i)*eps_0);
K_Bz2(i) = (2*dt*eps_0)/(dt*sig_x(i) + kap_x(i)*2*eps_0)/dz;

K_Hx2(i) = dt*((dt*sig_x(i) + kap_x(i)*2*eps_0)/(2*dt*eps_0))/(mu_r*mu_0);
K_Hx3(i) = dt*((dt*sig_x(i) - kap_x(i)*2*eps_0)/(2*dt*eps_0))/(mu_r*mu_0);

K_Bx1(i) = -(dt*sig_z(i) - 2*kap_z(i)*eps_0)/(dt*sig_z(i) + 2*kap_z(i)*eps_0);
K_Bx2(i) = (2*dt*eps_0)/(dt*sig_z(i) + kap_z(i)*2*eps_0)/dz;

K_Hz1(i) = 1;
K_Hz2(i) = dt*((dt*sig_z(i) + kap_z(i)*2*eps_0)/(2*dt*eps_0))/(mu);
K_Hz3(i) = -dt*((dt*sig_z(i) - kap_z(i)*2*eps_0)/(2*dt*eps_0))/(mu);

K_Ey1(i) = -(dt*sig_z(i) - 2*kap_z(i)*eps_0)/(dt*sig_z(i) + 2*kap_z(i)*eps_0);
K_Ey2(i) = (2*dt*eps_0)/(dt*sig_z(i) + kap_z(i)*2*eps_0)/(dt*eps);
K_Ey3(i) = -(2*dt*eps_0)/(dt*sig_z(i) + kap_z(i)*2*eps_0)/(dt*eps);
end
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
%ANISOTROPIC COEFFICIENTS ------------------------------------------------
        
        A1n = eps_inf/dt^2 - omp2.^2/4;
        A1p = eps_inf/dt^2 + omp2.^2/4;
        
        A2n = mu_inf/dt^2 - omp.^2/4;
        A2p = mu_inf/dt^2 + omp.^2/4;
        B2n = (mu_inf.*t2.^2+la2.*t1.^2)/dt^2 - omp.^2.*t2.^2/4;
        B2p = (mu_inf.*t2.^2+la2.*t1.^2)/dt^2 + omp.^2.*t2.^2/4;
        C2n = t1.*t2.*(mu_inf-la2)/dt^2 - omp.^2.*t2.*t1/4;
        C2p = t1.*t2.*(mu_inf-la2)/dt^2 + omp.^2.*t2.*t1/4;
        
        A3n = mu_inf/dt^2 - omp.^2/4;
        A3p = mu_inf/dt^2 + omp.^2/4;
        B3n = (mu_inf.*t1.^2+la2.*t2.^2)/(dt^2) - omp.^2.*t1.^2/4;
        B3p = (mu_inf.*t1.^2+la2.*t2.^2)/(dt^2) + omp.^2.*t1.^2/4;
        C3n = t1.*t2.*(mu_inf-la2)/dt^2 - omp.^2.*t1.*t2/4;
        C3p = t1.*t2.*(mu_inf-la2)/dt^2 + omp.^2.*t1.*t2/4;
        
    for i = 1:2*D_ax
    for j = 1:2*D_az
        A_Ey1(i,j) = 2*A1n(i,j)/A1p(i,j);
        A_Ey2(i,j) = 1/(eps_0*dt^2*A1p(i,j));

        A_Hz1(i,j) = 2*A2n(i,j)/A2p(i,j);
        A_Hz2(i,j) = 1;
        A_Hz3(i,j) = B2p(i,j) / (mu_0*la2(i,j)*A2p(i,j));
        A_Hz4(i,j) = -2*B2n(i,j) / (mu_0*la2(i,j)*A2p(i,j));
        A_Hz5(i,j) = B2p(i,j) / (mu_0*la2(i,j)*A2p(i,j));
        A_Hz6(i,j) = C2p(i,j) / (mu_0*la2(i,j)*A2p(i,j));
        A_Hz7(i,j) = -2*C2n(i,j) / (mu_0*la2(i,j)*A2p(i,j));
        A_Hz8(i,j) = C2p(i,j) / (mu_0*la2(i,j)*A2p(i,j));
        
        A_Hx1(i,j) = 2*A3n(i,j)/A3p(i,j);
        A_Hx2(i,j) = 1;
        A_Hx3(i,j) = B3p(i,j) / (mu_0*la2(i,j)*A3p(i,j));
        A_Hx4(i,j) = -2*B3n(i,j) / (mu_0*la2(i,j)*A3p(i,j));
        A_Hx5(i,j) = B3p(i,j) / (mu_0*la2(i,j)*A3p(i,j));
        A_Hx6(i,j) = C3p(i,j) / (mu_0*la2(i,j)*A3p(i,j));
        A_Hx7(i,j) = -2*C3n(i,j) / (mu_0*la2(i,j)*A3p(i,j));
        A_Hx8(i,j) = C3p(i,j) / (mu_0*la2(i,j)*A3p(i,j));
    end
    end
%----------------------------------------------------------------

