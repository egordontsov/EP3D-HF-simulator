%EP3D HF simulator for the geometry with symmetric stress barriers
%https://www.sciencedirect.com/science/article/abs/pii/S0013794415002404
%Egor Dontsov

clear all;clc;   
    
%Input problem parameters
Input = struct('K1C',1,...%MPa*m^1/2, fracture toughness
               'H',20,...%m, thickness of the reservoir layer
               'dsig',0.75,...%MPa, magnitude of the stress barrier
               'E',9.5,...%GPa, Young's modulus (the same for all layers)
               'nu',0.2,...%Poisson's ratio (the same for all layers)       
               'Cl',0.521/2*1e-2,...%mm/s^1/2, Carter's leak-off coefficient
               'mu',0.1,...%Pa*s, fluid viscosity
               'Q0',10,...%liters/s, injection rate
               't1',20,...%s, initial condition time
               't2',3600,...%s, end simulation time
               'tsave',[100 400 900 1800 3600],...%s, time instants at which fracture data is saved and plotted
               'Nx',50,...%number of grid points in x
               'Nz',50,...%number of grid points in z (for flux calculation)
               'Nt',2e2,...%number of time steps
               'alpha',1/2,...%time step distribution dt~t^alpha
               'animation',0,...%0 - no animation, 1 - do animation
               'save',0,...%0 - no saving, 1 - save animation
               'Nframes',100);%number of frames in animation
   
%run HF simulator
Results = EP3D_Simulator(Input);

%plotting options
%1 - plot, 0 - do not plot
options = struct('Wxz',0,...%w(x,z) at t(end)
                 'footpr',1,...%footprints, width and pressure cross-sections at times tsave
                 'Wxzt',0,...%3D polts at times tsave
                 'phLwt',1,...%time histories of wellbore width, pressure, fracture heigth, length, and volume
                 'plotrad',0,...%plot radial solution for comparison
                 'col','b-');%line color

%plot results
EP3D_plot(Results,options);

%save results
%save('EP3D_results.mat','Results');
