%EM cloaking FDTD 2D, diagonalized constitutive parameter tensor
%UPML, TFSF, no loss, planewave
%Oliver Csernyava BME Project Laboratory 1. \mail: csernyava.oliver@sch.bme.hu
% All rights reserved
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SET PARAMETERS FOR THE SIMULATION
Frequency   = 5e+10;    %[Hz]
Wavelength  = 20;       %[%]
PEC_radius  = 50;       %[%]
CLOAK_width = 50;       %[%]
PML_width   = 10;       %[%]
Pass_No     = 4;        %[pcs.]
Sample_No   = 10;       %[pcs.]
Min. permeability = 0.3;% MINIMUM CALCULATED PERMEABILITY

PEC_cylinder = 0;   %ACTIVATE PEC CYLINDER
CLOAK_cylinder = 0; %ACTIVATE CLOAK ANNULUS

Save_data = 1;      % ENABLE SAVING DATA
Frames_per_wavelength = 14; %FRAMES CAPTURED PER WAVELENGTH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_________________________________________________________________________%
%                                                                         %
%_________________________________________________________________________%
%                                                                         %
%_________________________________________________________________________%
%                                                                         %
%_________________________________________________________________________%
%                                                                         %
%_________________________________________________________________________%
%                                                                         %
%_________________________________________________________________________%
%                                                                         %
%_________________________________________________________________________%
%                                                                         %
%ALGORITHM STARTS HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name in the algorithm <-> Name in the data setup

freq = Frequency; %[Hz]
q = Wavelength; %[%]
PEC_radius = PEC_radius; %[%]
CLOAK_width = CLOAK_width; %[%]
q_pml = PML_width; %[%]
N_pass = Pass_No; %[pcs.]
sppw = Sample_No; %[pcs.]
la1_min = Min. permeability; %

PEC_cylinder = PEC_cylinder;
CLOAK_cylinder = CLOAK_cylinder;

capture_on = Save_data;
q_div = Frames_per_wavelength;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RUN algorithm

o_initialize
o_set
o_plotting
figure(1)
fig.MenuBar = 'none';
o_solver




