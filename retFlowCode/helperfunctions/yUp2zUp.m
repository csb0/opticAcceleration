function [outZup] = yUp2zUp(inYup, vargin)
%YUP2ZUP take in XYZ data where Y is up (e.g. Motion Shadow data) and spit
%it back with Z as the vertical dimention

if nargin == 2
    debug = vargin;
else
    debug = false;
end

outZup(:,1) = inYup(:,1);
outZup(:,2) = -inYup(:,3);
outZup(:,3) = inYup(:,2);


if debug
    figure(753)
    plot3(inYup(:,1), inYup(:,2), inYup(:,3),'r')
    hold on
    plot3(outZup(:,1), outZup(:,2), outZup(:,3),'b')
    axis equal
    grid on
end
end

