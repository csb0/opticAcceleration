function circlePoints3d = circleToPoints3d(varargin)
%DRAWCIRCLE3D Draw a 3D circle
%
%   Possible calls for the function:
%   drawCircle3d([XC YC ZC R THETA PHI])
%   drawCircle3d([XC YC ZC R], [THETA PHI])
%
%   where XC, YC, ZY are coordinates of circle center, R is the circle
%   radius, PHI and THETA are 3D angles in degrees of the normal to the
%   plane containing the circle:
%   * THETA between 0 and 180 degrees, corresponding to the colatitude
%       (angle with Oz axis).
%   * PHI between 0 and 360 degrees corresponding to the longitude (angle
%       with Ox axis)
%   
%   drawCircle3d([XC YC ZC R THETA PHI PSI])
%   drawCircle3d([XC YC ZC R], [THETA PHI PSI])
%   drawCircle3d([XC YC ZC R], THETA, PHI)
%   drawCircle3d([XC YC ZC], R, THETA, PHI)
%   drawCircle3d([XC YC ZC R], THETA, PHI, PSI)
%   drawCircle3d([XC YC ZC], R, THETA, PHI, PSI)
%   drawCircle3d(XC, YC, ZC, R, THETA, PHI)
%   drawCircle3d(XC, YC, ZC, R, THETA, PHI, PSI)
%   Are other possible syntaxes for this function.
%   
%   H = drawCircle3d(...)
%   return handle on the created LINE object
%
%   Example
%     % display 3 mutually orthogonal 3D circles
%     figure; hold on; 
%     drawCircle3d([10 20 30 50  0  0], 'LineWidth', 2, 'Color', 'b');
%     drawCircle3d([10 20 30 50 90  0], 'LineWidth', 2, 'Color', 'r');
%     drawCircle3d([10 20 30 50 90 90], 'LineWidth', 2, 'Color', 'g');
%     axis equal;
%     axis([-50 100 -50 100 -50 100]);
%     view([-10 20])
% 
%     % Draw several circles at once
%     center = [10 20 30];
%     circ1 = [center 50  0  0];
%     circ2 = [center 50 90  0];
%     circ3 = [center 50 90 90];
%     figure; hold on;
%     drawCircle3d([circ1 ; circ2 ; circ3]);
%     axis equal;
%
%   See also:
%   circles3d, drawCircleArc3d, drawEllipse3d, drawSphere
%
%   ------
%   Author: David Legland
%   e-mail: david.legland@grignon.inra.fr
%   Created: 2005-02-17
%   Copyright 2005 INRA - CEPIA Nantes - MIAJ (Jouy-en-Josas).

%   HISTORY
%   14/12/2006 allows unspecified PHI and THETA
%   04/01/2007 update doc, add todo for angle convention
%   19/06/2009 use localToGlobal3d, add drawing options
%   08/03/2010 use drawPolyline3d
%   2011-06-20 use angles in degrees, support several circles, update doc


%   Possible calls for the function, with number of arguments :
%   drawCircle3d([XC YC ZC R THETA PHI])            1
%   drawCircle3d([XC YC ZC R THETA PHI PSI])        1
%   drawCircle3d([XC YC ZC R], [THETA PHI])         2
%   drawCircle3d([XC YC ZC R], [THETA PHI PSI])     2
%   drawCircle3d([XC YC ZC R], THETA, PHI)          3
%   drawCircle3d([XC YC ZC], R, THETA, PHI)         4
%   drawCircle3d([XC YC ZC R], THETA, PHI, PSI)     4
%   drawCircle3d([XC YC ZC], R, THETA, PHI, PSI)    5
%   drawCircle3d(XC, YC, ZC, R, THETA, PHI)         6
%   drawCircle3d(XC, YC, ZC, R, THETA, PHI, PSI)    7


% extract drawing options
if verLessThan('matlab', '7.8')
    ind = find(cellfun('isclass', varargin, 'char'), 1, 'first');
else
    ind = find(cellfun(@ischar, varargin), 1, 'first');
end
options = {};
if ~isempty(ind)
    options = varargin(ind:end);
    varargin(ind:end) = [];
end

% Extract circle data
if length(varargin) == 1
    % get center and radius
    circle = varargin{1};
    xc = circle(:,1);
    yc = circle(:,2);
    zc = circle(:,3);
    r  = circle(:,4);
    
    % get colatitude of normal
    if size(circle, 2) >= 5
        theta = circle(:,5);
    else
        theta = zeros(size(circle, 1), 1);
    end

    % get azimut of normal
    if size(circle, 2)>=6
        phi     = circle(:,6);
    else
        phi = zeros(size(circle, 1), 1);
    end
    
    % get roll
    if size(circle, 2)==7
        psi = circle(:,7);
    else
        psi = zeros(size(circle, 1), 1);
    end
    


else
    error('drawCircle3d: please specify center and radius');
end

% circle parametrisation (by using N=60, some vertices are located at
% special angles like 45°, 30°...)
Nt  = 60;
t   = linspace(pi/8, 2*pi-pi/8, Nt+1);

nCircles = length(xc);
h = zeros(nCircles, 1);

for i = 1:nCircles
    % compute position of circle points
    x       = r(i) * cos(t)';
    y       = r(i) * sin(t)';
    z       = zeros(length(t), 1);
    circle0 = [x y z];

    % compute transformation from local basis to world basis
    trans   = localToGlobal3d_1(xc(i), yc(i), zc(i), theta(i), phi(i), psi(i));

    % compute points of transformed circle
    circle  = transformPoint3d(circle0, trans);

    % draw the curve of circle points
%     h(i) = drawPolyline3d(circle, options{:});
end

circlePoints3d = circle;

if nargout > 0
    varargout = {h};
end
