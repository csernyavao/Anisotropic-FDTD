%EM cloaking FDTD 2D, diagonalized constitutive parameter tensor
%UPML, TFSF, no loss, planewave
%Oliver Csernyava BME Project Laboratory 1. \mail: oliver.csernyava@sch.bme.hu

%run after you got the frame array from some simulation
filename = append(Filename,'.gif');
for k = 1:size(frame,2)
      if k > count_frame || Make_GIF == 0
          break
      end
      dT = Pass_Time;
      delay = dT/N_d;
      colorres = 256; % 0 - 256
      [imind,cm] = rgb2ind(frame2im(frame(k)),colorres); %colormap resolution (1-256), if the frame contains more colour levels, then the resultiing gif will be full of sparkling dots 
       %Write to the GIF File 
      if k == 1
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',delay); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',delay); 
      end
end