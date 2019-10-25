%EM cloaking FDTD 2D, diagonalized constitutive parameter tensor
%UPML, TFSF, no loss, planewave
%Oliver Csernyava BME Project Laboratory 1. \mail: oliver.csernyava@sch.bme.hu

%INITIALIZE PLOT-------------------------------------------------
[X,Z] = meshgrid(x,z); %create the mesh for the visualization
figure(1)
%figure('Name','Cloak FDTD','NumberTitle','off');
fig = figure(1);
set(fig,'Name','Cloak FDTD','NumberTitle','off')
sp = surf(X,Z,Ey);
hold;
%plot the structure
tt = linspace(0,2*pi,100);
stt = size(tt);
zt = ones(1, stt(2));
xt = n_x*dx/2 + r*dx*sin(tt);
yt = n_x*dx/2 + k*r*dx*cos(tt);
if PEC_cylinder
line(xt,yt,zt*E0*Q,'Color','white','LineWidth',1.5) %PEC radius line
end
xt = n_x*dx/2 + R*dx*sin(tt);
yt = n_x*dx/2 + k*R*dx*cos(tt);
if CLOAK_cylinder
line(xt,yt,zt*E0*Q,'Color','white','LineStyle',':','LineWidth',1.5) %Cloak outer radius line
end
[~,j_aux] = min(min(la1));%find the inner radius of the cloak
[~,i_aux] = min(la1(1:end,j_aux)); 
R_v = sqrt((D_az-j_aux)^2+k^4*(D_ax-i_aux)^2);
xt = n_x*dx/2 + R_v*dx*sin(tt);
yt = n_x*dx/2 + k*R_v*dx*cos(tt);
if CLOAK_cylinder
line(xt,yt,zt*E0*Q,'Color','white','LineStyle',':','LineWidth',1.5) %Cloak inner radius line
end
hold off;
sp.EdgeColor = 'none';
axis equal
colormap(jet(256))
axis([0 n_x*dx 0 n_x*dx -E0*Q*3 E0*Q*3])
caxis([-E0*Q E0*Q])
view(2)
pause on
%----------------------------------------------------------

%CAPTURE PROPERTIES ---------------------------------------
%     capture_on = app.SavedataCheckBox.Value;
    capture_iter = 1;
    capture_aux = 1;
    

%     q_div = app.FramesperwavelengthEditField.Value; %captured frames while the wavefront passes one wavelength
    Q_div = q_div / (q*0.01); % captured frames while the wavefront passes the simulation domain once
    N_i = n_t/N_pass; % iterations while the wavefront passes the simulation domain once
    capture_div = ceil(N_i/Q_div); % frequency of frame sampling for the data visualization
    if capture_div <2 && capture_on == 1
        capture_div = 2;
        q_div = floor(N_i / capture_div *q*0.01);
        msg = sprintf('app.capture_div should be greater than 2 -> Frames per wavelength is set to %d', sppw);
        warning(msg)
    end
    N_d = ceil(n_t/capture_div)/N_pass; %corrected captured frames while the wavefront passes the simulation domain once
    count_frame = 0;
    frame(ceil(n_t/capture_div)) = struct('cdata',[],'colormap',[]);
%---------------------------------------------------------