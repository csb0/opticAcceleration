function [sphCurl_F, sphDiv_F] = del_in_spherical(dTheta_F,dRho_F,thetaRes,rhoRes,rhoScale)
%% spherical curl and divergence based on:
% https://en.wikipedia.org/wiki/Del_in_cylindrical_and_spherical_coordinates

%% inputs
% dTheta_F and dRho_F are scattered interpolants of dRho and dTheta values
% at cartesian coordinates (x,y)
% res are the number of steps sweeping from 0:2*pi for theta, and
% 0:rhoScale (1/2 of grid width for interpolant) for rho


%% construct polar grid
thetaSpacing = 1/(thetaRes-1);
rhoSpacing = 1/(rhoRes-1);
[tt rr] = meshgrid(linspace(thetaSpacing,2*pi,thetaRes),linspace(rhoSpacing,rhoScale,rhoRes));




%% get dTheta and dRho polar grids
[xx yy] = pol2cart(tt,rr);
% xx = xx+rhoScale;
% yy = yy+rhoScale;

% scattered interpolant functions in cartesian space
dTheta_grid = dTheta_F(xx,yy);
dRho_grid = dRho_F(xx,yy);

%% compute grid of sin(rho) vals for spherical scaling

sinGrid = sin(pi-rr/rhoScale*pi/2);

%% compute gradients, also scale by step size

% non sin scaled
[dvTheta_dTheta dvTheta_dRho ] = gradient(dTheta_grid);
[dvRho_dTheta dvRho_dRho ] = gradient(dRho_grid);

% scale by step size
dvTheta_dRho = dvTheta_dRho/rhoSpacing;
dvRho_dRho = dvRho_dRho/rhoSpacing;
dvRho_dTheta = dvRho_dTheta/thetaSpacing;
dvTheta_dTheta = dvTheta_dTheta/thetaSpacing;

% sin scaled
[dvThetaSin_dTheta dvThetaSin_dRho] = gradient(dTheta_grid.*sinGrid);
[dvRhoSin_dTheta dvRhoSin_dRho ] = gradient(dRho_grid.*sinGrid);

% scale by step size
dvThetaSin_dRho = dvThetaSin_dRho/rhoSpacing;
dvRhoSin_dRho = dvRhoSin_dRho/rhoSpacing;
dvRhoSin_dTheta = dvRhoSin_dTheta/thetaSpacing;
dvThetaSin_dTheta = dvThetaSin_dTheta/thetaSpacing;

%% compute divergence
% divergence is (1/sin(rho))*(dvRhoSin_dRho + dvTheta_dTheta) 
      
sphDiv = 1./sinGrid.*(dvRhoSin_dRho + dvTheta_dTheta);
sphDiv_F = scatteredInterpolant(xx(:),yy(:),sphDiv(:));

%% compute curl
% curl is (1/sin(rho))*(dvThetaSin_dRho + dvRho_dTheta) 

sphCurl = 1./sinGrid.*(dvThetaSin_dRho + dvRho_dTheta);
sphCurl_F = scatteredInterpolant(xx(:),yy(:),sphCurl(:));










end

