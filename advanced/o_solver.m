%EM cloaking FDTD 2D, diagonalized constitutive parameter tensor
%UPML, TFSF, no loss, planewave
%Oliver Csernyava BME Project Laboratory 1. \mail: oliver.csernyava@sch.bme.hu

% All rights reserved

%MAIN LOOP##############################################################
%initialization
%%
t = 0;
Bz_o = zeros(n_x,n_z);
Bx_o = zeros(n_x,n_z);
Dy_o = zeros(n_x,n_z);
Hz_1 = zeros(n_x,n_z);
Hz_2 = zeros(n_x,n_z);
Hz_3 = zeros(n_x,n_z);
Hx_1 = zeros(n_x,n_z);
Hx_2 = zeros(n_x,n_z);
Hx_3 = zeros(n_x,n_z);
Ey_1 = zeros(n_x,n_z);
Ey_2 = zeros(n_x,n_z);
Ey_3 = zeros(n_x,n_z);
%%

while (t < n_t*dt)
    if ~ishandle(fig) %if figure is closed, then brake the loop
        break
    end
    %%
    %old values stored for the FDTD algorithm
    Dy_oo = Dy_o;
    Dy_o = Dy; 
    Bz_oo = Bz_o;
    Bz_o = Bz;
    Bx_oo = Bx_o;
    Bx_o = Bx;
    Hx_o = Hx;
    Hz_o = Hz;
    Ey_o = Ey;
    %%
% AUX FIELD excitation and MUR ABC components-------------------------
g = E0*sin(2*pi*(freq)*t); %excitation function
Ey_aux(D_pml-1) = Ey_aux(D_pml-1) + g;

Ey_n_1 = Ey_n; %the electric field components required for the MUR ABC
Ey_n(1:2)=Ey_aux(1:2);
Ey_n(3:4)=Ey_aux(n_x-1:n_x);
%---------------------------------------------------------------------

%AUX H FIELD-------------------------------------------------------
Hz_aux = Hz_aux - K_M*diff(Ey_aux); %update aux field
%%
%ELECTRIC FIELD_______________________________________________________
    %UPML-------------------------------------------------------------
    for j = 2:n_z-1
    for i = [2:D_pml,n_x-D_pml+1:n_x-1]
    Dy(i,j) = K_Dy1(i)*Dy(i,j) - (K_Dy2(i)*(Hz(i,j)-Hz(i-1,j)) - K_Dy3(i)*(Hx(i,j)-Hx(i,j-1)));
    Ey(i,j) = K_Ey1(j)*Ey(i,j) + (K_Ey2(j)*Dy(i,j) + K_Ey3(j)*Dy_o(i,j));
    end
    end

    for j = [2:D_pml,n_z-D_pml+1:n_z-1]
    for i = D_pml+1:n_x-D_pml
    Dy(i,j) = K_Dy1(i)*Dy(i,j) - (K_Dy2(i)*(Hz(i,j)-Hz(i-1,j)) - K_Dy3(i)*(Hx(i,j)-Hx(i,j-1)));
    Ey(i,j) = K_Ey1(j)*Ey(i,j) + (K_Ey2(j)*Dy(i,j) + K_Ey3(j)*Dy_o(i,j));
    end
    end
    %--------------------------------------------------------------------------

    %HOMOGENOUS SPACE + TF/SF -------------------------------------------------
    for j = D_pml+1:n_z-D_pml
    for i = [D_pml+1:x_a-D_ax-1 , x_a+D_ax+1:n_x-D_pml]
        Ey(i,j) = Ey(i,j) - (C_Eyz*(Hz(i,j)-Hz(i-1,j)) - C_Eyx*(Hx(i,j)-Hx(i,j-1)));
    end
    end

    for j = [D_pml+1:z_a-D_az-1 , z_a+D_az+1:n_z-D_pml]
    for i = x_a-D_ax : x_a+D_ax
        Ey(i,j) = Ey(i,j) - (C_Eyz*(Hz(i,j)-Hz(i-1,j)) - C_Eyx*(Hx(i,j)-Hx(i,j-1)));
    end
    end        

    i = D_pml+2; %TFSF
    for j = D_pml+2:n_z-D_pml-1
        Ey(i,j) = Ey(i,j) + C_Eyz*Hz_aux(i-1);
    end

    i = n_x-D_pml-1; %TFSF 
    for j = D_pml+2:n_z-D_pml-1
        Ey(i,j) = Ey(i,j) - C_Eyz*Hz_aux(i-1);
    end
    %----------------------------------------------------------------    

    %--------------------------------------------------------------------
    %ANISOTROPIC FIELD --------------------------------------------------
    for j = z_a-D_az : z_a+D_az %first update the constitutive equ
    for i = x_a-D_ax : x_a+D_ax
        Dy(i,j) = K_Dy1(i)*Dy(i,j) - (K_Dy2(i)*(Hz(i,j)-Hz(i-1,j)) - K_Dy3(i)*(Hx(i,j)-Hx(i,j-1)));
    end %important to calculate the D values before the H values, for the interpolation
    end

    for j = z_a-D_az : z_a+D_az %then update the Maxwell equ
    for i = x_a-D_ax : x_a+D_ax
        xi = i - (x_a-D_ax) +1;
        zj = j - (z_a-D_az) +1;
        if (sqrt((z_a-j)^2+(x_a-i)^2) > r) &&  (sqrt((z_a-j)^2+(x_a-i)^2) <= R)
           %inside the cloaking material 
            Ey_1(i,j) = Ey(i,j);
            Ey(i,j) = A_Ey1(xi,zj)*Ey_2(i,j) - Ey_3(i,j) + A_Ey2(xi,zj)*( Dy(i,j) - 2*Dy_o(i,j) + Dy_oo(i,j) );
            Ey_2(i,j) = Ey(i,j);
            Ey_3(i,j) = Ey_1(i,j);

        end
        if (sqrt((z_a-j)^2+(x_a-i)^2) > R) %outside the cloaking material (vacuum)
            Ey(i,j) = K_Ey1(j)*Ey(i,j) + (K_Ey2(j)*Dy(i,j) + K_Ey3(j)*Dy_o(i,j));
        end

        if (sqrt((z_a-j)^2+(x_a-i)^2) <= r) %the inner object
            if PEC_cylinder
                Ey(i,j) = 0; %--->PEC
            else
                Ey(i,j) = K_Ey1(j)*Ey(i,j) + (K_Ey2(j)*Dy(i,j) + K_Ey3(j)*Dy_o(i,j));
            end

        end
    end
    end
    %---------------------------------------------------------------   
%___________________________________________________________________
%%
% AUX E FIELD for the plane wave TF/SF--------------------------------
 Ey_aux(2:n_x-1) = Ey_aux(2:n_x-1) - K_E1*diff(Hz_aux(1:n_x-1));  % aux

    %Mur ABC 
   %at x = 0;
    Ey_aux(1)   = -Ey_n_1(2) + K_M1*(Ey_aux(2) + Ey_n_1(1)) + K_M2*(Ey_n(1) + Ey_n(2));
   %at x = h
    Ey_aux(n_x) = -Ey_n_1(3) + K_M1*(Ey_aux(n_x-1) + Ey_n_1(4)) + K_M2*(Ey_n(4) + Ey_n(3));
%---------------------------------------------------------------------
%%
%MAGNETIC FIELD _________________________________________________
    % PML ------------------------------------------------------------
    for j = 1:n_z-1
    for i = [1:D_pml,n_x-D_pml+1:n_x-1]
    Bx(i,j) = K_Bx1(j)*Bx(i,j) + K_Bx2(j)*(Ey(i,j+1)-Ey(i,j));
    Hx(i,j) = Hx(i,j) + (K_Hx2(i)*Bx(i,j) + K_Hx3(i)*Bx_o(i,j));
    Bz(i,j) = K_Bz1(i)*Bz(i,j) - K_Bz2(i)*(Ey(i+1,j)-Ey(i,j));      
    Hz(i,j) = K_Hz1(j)*Hz(i,j) + (K_Hz2(j)*Bz(i,j) - K_Hz3(j)*Bz_o(i,j));
    end
    end

    for j = [1:D_pml,n_z-D_pml+1:n_z-1]
    for i = D_pml+1:n_x-D_pml
    Bx(i,j) = K_Bx1(j)*Bx(i,j) + K_Bx2(j)*(Ey(i,j+1)-Ey(i,j));
    Hx(i,j) = Hx(i,j) + (K_Hx2(i)*Bx(i,j) + K_Hx3(i)*Bx_o(i,j));
    Bz(i,j) = K_Bz1(i)*Bz(i,j) - K_Bz2(i)*(Ey(i+1,j)-Ey(i,j));      
    Hz(i,j) = K_Hz1(j)*Hz(i,j) + (K_Hz2(j)*Bz(i,j) - K_Hz3(j)*Bz_o(i,j));
    end
    end
    %------------------------------------------------------------------

    %HOMOGENOUS FIELD + TF/SF- ----------------------------------------
    for j = D_pml+1:n_z-D_pml
    for i = [D_pml+1:x_a-D_ax-1 , x_a+D_ax+1:n_x-D_pml]
    Hx(i,j) = Hx(i,j) + C_Hx*(Ey(i,j+1)-Ey(i,j));
    Hz(i,j) = Hz(i,j) - C_Hz*(Ey(i+1,j)-Ey(i,j));      
    end
    end

    for j = [D_pml+1:z_a-D_az-1 , z_a+D_az+1:n_z-D_pml]
    for i = x_a-D_ax : z_a+D_az
    Hx(i,j) = Hx(i,j) + C_Hx*(Ey(i,j+1)-Ey(i,j));
    Hz(i,j) = Hz(i,j) - C_Hz*(Ey(i+1,j)-Ey(i,j));      
    end
    end

    i = D_pml+1; %TFSF
    for j = D_pml+2:n_z-D_pml-1
        Hz(i,j) = Hz(i,j) + C_Hz*Ey_aux(i+1);
    end

    i = n_x-D_pml-2; %TFSF
    for j = D_pml+2:n_z-D_pml-1
        Hz(i,j) = Hz(i,j) - C_Hz*Ey_aux(i+1);
    end

    j = D_pml+1; %TFSF
    for i = D_pml+2:n_x-D_pml-2
        Hx(i,j) = Hx(i,j) - C_Hx*Ey_aux(i);
    end

    j = n_x-D_pml-1; %TFSF
    for i = D_pml+2:n_x-D_pml-2
        Hx(i,j) = Hx(i,j) + C_Hx*Ey_aux(i);
    end
    %----------------------------------------------------------------
    %ANISOTROPIC FIELD ----------------------------------------------
    for j = z_a-D_az : z_a+D_az
    for i = x_a-D_ax : x_a+D_ax
        Bx(i,j) = K_Bx1(j)*Bx(i,j) + K_Bx2(j)*(Ey(i,j+1)-Ey(i,j));
        Bz(i,j) = K_Bz1(i)*Bz(i,j) - K_Bz2(i)*(Ey(i+1,j)-Ey(i,j));
    end
    end

    for j = z_a-D_az : z_a+D_az
    for i = x_a-D_ax : x_a+D_ax
        ix = i-(x_a-D_ax)+1;
        jz = j-(z_a-D_az)+1;
        %field outside the cloaking annulus
        if (sqrt((z_a-j)^2+(x_a-i)^2) <= r) ||  (sqrt((z_a-j)^2+(x_a-i)^2) > R)
            Hx(i,j) = Hx(i,j) + (K_Hx2(i)*Bx(i,j) + K_Hx3(i)*Bx_o(i,j));    
            Hz(i,j) = K_Hz1(j)*Hz(i,j) + (K_Hz2(j)*Bz(i,j) - K_Hz3(j)*Bz_o(i,j));
        end
        %field inside the cloaking annulus
        if (sqrt((z_a-j)^2+(x_a-i)^2) > r) &&  (sqrt((z_a-j)^2+(x_a-i)^2) <= R)
            %INTERPOLATION from the nearby data points
            Bx1 = ( Bx(i,j)+Bx(i,j-1)+Bx(i+1,j)+Bx(i+1,j-1) )/4;
            Bx2 = ( Bx_o(i,j)+Bx_o(i,j-1)+Bx_o(i+1,j)+Bx_o(i+1,j-1) )/4;
            Bx3 = ( Bx_oo(i,j)+Bx_oo(i,j-1)+Bx_oo(i+1,j)+Bx_oo(i+1,j-1) )/4;
            Bz1 = ( Bz(i,j)+Bz(i,j+1)+Bz(i-1,j+1)+Bz(i-1,j) )/4;
            Bz2 = ( Bz_o(i,j)+Bz_o(i,j+1)+Bz_o(i-1,j+1)+Bz_o(i-1,j) )/4;
            Bz3 = ( Bz_oo(i,j)+Bz_oo(i,j+1)+Bz_oo(i-1,j+1)+Bz_oo(i-1,j) )/4;

            Hz_1(i,j) = Hz(i,j); %Save the old value
            Hz(i,j) = A_Hz1(ix,jz)*Hz_2(i,j) - Hz_3(i,j) + A_Hz3(ix,jz)*Bz(i,j) + A_Hz4(ix,jz)*Bz_o(i,j) + A_Hz5(ix,jz)*Bz_oo(i,j) + A_Hz6(ix,jz)*Bx1 + A_Hz7(ix,jz)*Bx2 + A_Hz8(ix,jz)*Bx3;
            Hz_2(i,j) = Hz(i,j); %Save the recent value
            Hz_3(i,j) = Hz_1(i,j); %Save the old old value
            Hx_1(i,j) = Hx(i,j); %Save the old value
            Hx(i,j) = A_Hx1(ix,jz)*Hx_2(i,j) - Hx_3(i,j) + A_Hx3(ix,jz)*Bx(i,j) + A_Hx4(ix,jz)*Bx_o(i,j) + A_Hx5(ix,jz)*Bx_oo(i,j) + A_Hx6(ix,jz)*Bz1 + A_Hx7(ix,jz)*Bz2 + A_Hx8(ix,jz)*Bz3;
            Hx_2(i,j) = Hx(i,j); %Save the recent value
            Hx_3(i,j) = Hx_1(i,j); %Save the old old value
        end    
    end
    end
%%
%UPDATING PLOTS-------------------------------------------- 
title(['t =' num2str(t,'%.2e') '[sec]' ], 'FontSize', 20);
sp.ZData = Ey;
%pause(eps);
drawnow
%-----------------------------------------------------------
%%
%CAPTURING THE ANIMATION ------------------------------------
if (capture_on == 1) && (capture_aux == capture_div)   
  frame(capture_iter) = getframe(fig); %saving the frames --> Make_GIF.m
  count_frame=count_frame+1;
    capture_iter = capture_iter + 1;
    capture_aux = 1;
end 
capture_aux = capture_aux+1;
%---------------------------------------------------------------------
%%
t = t + dt; %timestepping
end
