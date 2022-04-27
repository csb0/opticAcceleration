function [zz] = gaussianPlot(xx, yy, xc, yc, sigma)
%GAUSSIANPLOT takes in x and y meshes, x and y centers, a sigma, and spits out a a big ol' normalized bump of a gaussian
%equation from https://www.mathworks.com/matlabcentral/answers/13020-2d-gaussian-function#comment_250586

%%
exponent = ((xx-xc).^2 + (yy-yc).^2)./(2*sigma^2);
amplitude = 1 / (sigma * sqrt(2*pi));  

zz       = amplitude  * exp(-exponent);

zz = zz./max(zz(:)); %normalize that bad boy

debug = false;
if debug
figure(5734)
clf
surface(xx,yy,zz); 
view(3)

end

