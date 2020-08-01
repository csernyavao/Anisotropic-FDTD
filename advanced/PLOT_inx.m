%EM cloaking FDTD 2D, diagonalized constitutive parameter tensor
%UPML, TFSF, no loss, planewave
%Oliver Csernyava BME Project Laboratory 1. \mail: csernyava.oliver@sch.bme.hu

% All rights reserved

exa = matfile('Sample.mat');
fig = exa.fig;

domain = 98; % [%]

k = size(fig.M,3);
frame = ceil(domain*0.01*k);

f = figure(1);
imshow(I,fig.M(:,:,frame));
