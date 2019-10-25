%EM cloaking FDTD 2D, diagonalized constitutive parameter tensor
%UPML, TFSF, no loss, planewave
%Oliver Csernyava BME Project Laboratory 1. \mail: oliver.csernyava@sch.bme.hu

%DEFINE PARAMETERS
eps_0 = 8.854*1e-12;
mu_0  = 4*pi*1e-7;
c   = 1/sqrt(eps_0*mu_0); %the speed of the light in the medium

eps_r = 1; %default vacuum
mu_r  = 1;

n_z = 400; %UNIFORM GRID (number of spatial data points)
n_x = n_z;
%sppw = app.SampleNoEditField.Value; %minimum sample points per wavelength

%SET REQUIRED SPATIAL STEP (dx) (iteration)--------------------------------
%Before defining the required spatial and temporal resolution, the
%constitutional parameters of the cloak should be calculated. However, for
%these parameters we need a starting resolution, which is an estimation.
%This estimation should be further approximated with an iteration. The
%iteration stops, if the sample points per wavelength (sppw) for the very extremes
%of the constitutional parameter values are above a defined value (sppw).
%Thus, no undersampling happens.
spat_flag = 1; %if 0 -> iteration stops
while (spat_flag) %initialize required spatial step
    if exist('dx','var') %if the dx is already defined, then recalculate it with the results from the previous iteration
        n_z = ceil(lambda_free/dx/q/0.01); %set the spatial sample points to fit with the displayed wavelength requirements (eg. wavelength is the q % of the range)
        n_x = n_z;
    end        

% q = app.WavelengthEditField.Value;%[%] wavelength compared to the computational domain (percentage)
% q_pml = app.PMLwidthEditField.Value;%[%] PML width compared to the computational domain (percentage)

z_a = floor(n_z/2); %center of the anisotropic medium
x_a = floor(n_x/2);

E0 =1; %Amplitude of the incident planewave
% freq = app.FrequencyHzEditField.Value; % Hz, frequency of the incident wave
lambda_free = c / (freq); %vacuum wavelength

r = ceil(n_x*q*0.01*PEC_radius*0.01);  % radius relative to the half width of the anisotropic domain
R = r + ceil(n_x*q*0.01*CLOAK_width*0.01);


D_ax = R+2; %correction needed for calculation: the calculation uses the radius, and not the diameter!
D_az = R+2; %correction needed for calculation

orrf = 1; %outer ring flag, 1->cloak outer ring is visible
k = 1; %ratio of the ellipsoid ---> fix 1 ! the algorithm is not tested for ellipse cloak

% N_pass = app.PassNoEditField.Value; %Wavefront pass No. (the numbert of the passage of the wavfront on the simu range until the simu time)

m = 3.7; % PML parameters
Rerr = 1.0e-16; % error tolerance
kapp_max = 10;

diag_err = 0; % for diagonalization error

%CALCULATE PARAMETERS_____________________________________________________
eps = eps_r*eps_0;
mu = mu_r*mu_0;


%PARAMETERS OF THE ANISOTROPIC MEDIUM---------------------------
    A_Hz1 = zeros(2*D_ax+1,2*D_az+1); %initialization of variables_________
    A_Hz2 = zeros(2*D_ax+1,2*D_az+1);
    A_Hz3 = zeros(2*D_ax+1,2*D_az+1);
    A_Hz4 = zeros(2*D_ax+1,2*D_az+1);
    A_Hz5 = zeros(2*D_ax+1,2*D_az+1);
    A_Hz6 = zeros(2*D_ax+1,2*D_az+1);
    A_Hz7 = zeros(2*D_ax+1,2*D_az+1);
    A_Hz8 = zeros(2*D_ax+1,2*D_az+1);
    
    A_Hx1 = zeros(2*D_ax+1,2*D_az+1);
    A_Hx2 = zeros(2*D_ax+1,2*D_az+1);
    A_Hx3 = zeros(2*D_ax+1,2*D_az+1);
    A_Hx4 = zeros(2*D_ax+1,2*D_az+1);
    A_Hx5 = zeros(2*D_ax+1,2*D_az+1);
    A_Hx6 = zeros(2*D_ax+1,2*D_az+1);
    A_Hx7 = zeros(2*D_ax+1,2*D_az+1);
    A_Hx8 = zeros(2*D_ax+1,2*D_az+1);
    
    A_Ey1 = zeros(2*D_ax+1,2*D_az+1);
    A_Ey2 = zeros(2*D_ax+1,2*D_az+1);
    
    mu_xx = 1*ones(2*D_ax+1,2*D_az+1);
    mu_zx = 0*ones(2*D_ax+1,2*D_az+1);
    mu_xz = 0*ones(2*D_ax+1,2*D_az+1);
    mu_zz = 1*ones(2*D_ax+1,2*D_az+1);
    eps_yy = 1*ones(2*D_ax+1,2*D_az+1);
    eps_inf = eps_yy; 
    sig_y = zeros(2*D_ax,2*D_az);
    
    la1 = ones(2*D_ax+1,2*D_az+1);
    la2 = la1;
    mu_inf = ones(2*D_ax+1,2*D_az+1);
    t1 = ones(2*D_ax+1,2*D_az+1);
    t2 = zeros(2*D_ax+1,2*D_az+1);
    ksi1 = ones(2*D_ax+1,2*D_az+1);
    ksi2 = zeros(2*D_ax+1,2*D_az+1);
    bet = zeros(2*D_ax+1,2*D_az+1);
    alpha = zeros(2*D_ax+1,2*D_az+1);  %end of initialization block__________
    
    %--------------------------------------------------------------------
    %CLOAKING MATERIAL---------------------------------------------------
    for j = z_a-D_az : z_a+D_az
    for i = x_a-D_ax : x_a+D_ax
        i_x = i-(x_a-D_ax)+1;
        j_z = j-(z_a-D_az)+1;
        if (sqrt((z_a-j)^2+(x_a-i)^2) > r) &&  (sqrt((z_a-j)^2+(x_a-i)^2) <= R)
        rad = sqrt((z_a-j)^2+(x_a-i)^2);
        RAD = sqrt((z_a-j)^2+(x_a-i)^2);

        rr = (rad-r)/rad;
        rn = 1/rr;
        
        tet = atan((i-x_a)/(j-z_a));
            
        %CLOAK parameters---------------------------------
        if CLOAK_cylinder    
            mu_zz(i_x,j_z) = rr*cos(tet)^2 + sin(tet)^2/rr;
            mu_xx(i_x,j_z) = rr*sin(tet)^2 + cos(tet)^2/rr;
            mu_xz(i_x,j_z) = rr*sin(tet)*cos(tet) - sin(tet)*cos(tet)/rr;
            mu_zx(i_x,j_z) = mu_xz(i_x,j_z);
            eps_yy(i_x,j_z) = (R/(R-r))^2*rr;
        else %COSTUME parameters-------------------------------
            mu_zz(i_x,j_z) = 1;
            mu_xx(i_x,j_z) = 1;
            mu_xz(i_x,j_z) = 0;
            mu_zx(i_x,j_z) = mu_xz(i_x,j_z);
            eps_yy(i_x,j_z) = 1;
        end

        %EIGENVALUES-----
         la1(i_x,j_z) = ((mu_zz(i_x,j_z) + mu_xx(i_x,j_z)) - sqrt((mu_zz(i_x,j_z) + mu_xx(i_x,j_z))^2-4*(mu_xx(i_x,j_z)*mu_zz(i_x,j_z)-mu_xz(i_x,j_z)^2)))/2;
         la2(i_x,j_z) = ((mu_zz(i_x,j_z) + mu_xx(i_x,j_z)) + sqrt((mu_zz(i_x,j_z) + mu_xx(i_x,j_z))^2-4*(mu_xx(i_x,j_z)*mu_zz(i_x,j_z)-mu_xz(i_x,j_z)^2)))/2;                

        mu_inf(i_x,j_z) = max(1,la1(i_x,j_z));
        eps_inf(i_x,j_z) = max(1,eps_yy(i_x,j_z));

        bet(i_x,j_z) = sqrt((la2(i_x,j_z) - mu_xx(i_x,j_z))^2 + mu_xz(i_x,j_z)^2);
        if bet(i_x,j_z)<1e-12   %manage calculation error
            bet(i_x,j_z)=0;     
        end

        t1(i_x,j_z) = mu_xz(i_x,j_z)/bet(i_x,j_z);
        t2(i_x,j_z) = (la2(i_x,j_z)-mu_xx(i_x,j_z))/bet(i_x,j_z);

        %SET the minimum value of lambda in the simulation
            %small values-->long time simulation run
%              la1_min = app.MinpermeabilityEditField.Value;
             if la1(i_x,j_z) < la1_min
                 la1(i_x,j_z) = 1;
                 la2(i_x,j_z) = 1;
                 mu_xx(i_x,j_z) = 1;
                 mu_zz(i_x,j_z) = 1;
                 mu_xz(i_x,j_z) = 0;
                 mu_zx(i_x,j_z) = 0;
                 eps_yy(i_x,j_z) = 1;
                 t1(i_x,j_z) = 1;
                 t2(i_x,j_z) = 0;
             end
             
        %MANAGE singular values of t1, t2 ---> continous values required
        if (isnan(t1(i_x,j_z)))
            t1(i_x,j_z) = 1;
        end
        if (isnan(t2(i_x,j_z))) || (isinf(t2(i_x,j_z)))
           t2(i_x,j_z) = 0;
        end
        end 
    end
    end
    %____________________________________________________________________
    
    %calculate plasmonic frequencies of the drude model
    omp = 2*pi*freq*sqrt(mu_inf - la1); %for the mu (lambda_1)
    omp2 = 2*pi*freq*sqrt(eps_inf - eps_yy); % for the epsilon_yy

%--------------------------------------------
%CHECK TRANSFORM MATRICES -> check if the diagonalization is OK
for j = z_a-D_az+1 : z_a+D_az
for i = x_a-D_ax+1 : x_a+D_ax
    i_x = i-(x_a-D_ax+1)+1; %relative indexes appropriate for the parameters only calculated for the inner anysotropic medium range
    j_z = j-(z_a-D_az+1)+1;
    if (sqrt((z_a-j)^2+k^2*(x_a-i)^2) >= r) &&  (sqrt((z_a-j)^2+k^4*(x_a-i)^2) <= R) %only inside the cloak
    AM = [la1(i_x,j_z) 0; 0 la2(i_x,j_z)]; %diagonal mx
    P = [t1(i_x,j_z) t2(i_x,j_z); -t2(i_x,j_z) t1(i_x,j_z)]; %the invertible projection mx
    EPS = P*AM*transpose(P); %calculate the original tensor from the diagonal mx
    Err = EPS - [mu_zz(i_x,j_z) mu_zx(i_x,j_z); mu_xz(i_x,j_z) mu_xx(i_x,j_z)]; %%calculate the diffeence
    if max(max(Err)) > 1e-10 || (abs(min(min(Err))) > 1e-10)
        diag_err = 1;
    end
    end
end
end
%--------------------------------------------------

% SPATIAL RESOLUTION-------------------------------
%SET spatial resolution
omega = 2*pi*freq;
lambda_min = c / (sqrt( max(max(la2))*max(max(eps_yy) ))*freq); % minimum wavelength in the simulation domain

if ~exist('dx','var') %first iteration
    lambda0_rel = sppw; % sample points per wavelength
    dx = 1*lambda_min / lambda0_rel; %spatial step
    %check the PML width requirement
    n_xi = ceil(lambda_free/dx/q/0.01);
    D_pml = ceil(n_xi*q_pml*0.01); %number of the spatial points in the PML (width)
    if D_pml < 15 %experienced minimum value for the PML width (cells)
        n_xi = 15/(q_pml*0.01); %required datapoints for the minimum PML width
        dx = lambda_free/n_xi/q/0.01; %required dx for the minimum PML width
    end
else %further iteration steps
    if (lambda_min / dx) >= sppw %minimum sample points in a wavelength
        spat_flag = 0; %---> stop initializing iteration
        if diag_err == 1
            msg = 'Inproper transformation matrix';
            error(msg)
        end
    else
        lambda0_rel = sppw; % sample points per wavelength
        dx = 1*lambda_min / lambda0_rel; %spatial step
    end
end

dz = dx;
end

%END of the spatial initialization---------------------------------------------------


