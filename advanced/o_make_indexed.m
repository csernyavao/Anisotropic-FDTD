%EM cloaking FDTD 2D, diagonalized constitutive parameter tensor
%UPML, TFSF, no loss, planewave
%Oliver Csernyava BME \mail: csernyava.oliver@sch.bme.hu

% All rights reserved

%Saves the frames to a 4D indexed graphics file
%fig.IM -> Indices (at the size of the figure)
%fig.M -> Color map of the figure (max 256x3 matrix)
%open the resulting .mat file:
    % example = matfile('sample.mat');
    % fig = example.fig;
%run after you got the frame array from some simulation

filename = append(Filename,'.mat');

[imind,cm] = rgb2ind(frame2im(frame(1)),256); %get the size of the image from the first frame

fig = struct;

fig.IM = zeros(size(imind,1), size(imind,2));
fig.IM = uint8(fig.IM);
fig.M = zeros(45,3);

for k = 1:size(frame,2)
      if k > count_frame || Make_indexed == 0
          break
      end
      colorres = 256; % 0 - 256
      [imind,cm] = rgb2ind(frame2im(frame(k)),colorres); %colormap resolution (1-256), if the frame contains more colour levels, then the resultiing gif will be full of sparkling dots 
       
      fig.IM(:,:,k) = imind;
      for i = 1:size(cm,1)
        fig.M(i,:,k) = cm(i,:);
      end
end

save(filename,'fig');

