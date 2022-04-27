function [backEyeProjPoints_xyz_id_rot, retinaProjPoints_fr_X_Y_id,headX,headY] = defineRetinalPolarCoordinates(thisEyeStruct,back_proj_heading,fr, debug, spotcheck)
% [backEyeProjPoints_xyz_id_rot, retinaProjPoints_fr_theta_rho_id,
% retinaProjPoints_fr_X_Y_id,fixX,fixY] =
% defineRetinalPolarCoordinates(thisEyeStruct,back_proj_heading,fr, debug,
% spotcheck) 5/23


varNames = fieldnames(thisEyeStruct);
for i=1:length(varNames)
    eval([varNames{i} '=thisEyeStruct.' varNames{i} ';']);
end



thisFixPoint = fixXYZ(fr,:);

backEyeProjPoints_xyz_id_z = nan(size(backEyeProjPoints_xyz_id));
%zero projected points and foveaXYZ to eyeball center
backEyeProjPoints_xyz_id_z(:,1) = backEyeProjPoints_xyz_id(:,1) - eyeXYZ(fr,1);
backEyeProjPoints_xyz_id_z(:,2) = backEyeProjPoints_xyz_id(:,2) - eyeXYZ(fr,2);
backEyeProjPoints_xyz_id_z(:,3) = backEyeProjPoints_xyz_id(:,3) - eyeXYZ(fr,3);
backEyeProjPoints_xyz_id_z(:,4) = backEyeProjPoints_xyz_id(:,4);

% do the same thing to back_proj_heading

back_proj_heading(1) = back_proj_heading(1) - eyeXYZ(fr,1);
back_proj_heading(2) = back_proj_heading(2) - eyeXYZ(fr,2);
back_proj_heading(3) = back_proj_heading(3) - eyeXYZ(fr,3);


foveaXYZ_z = foveaXYZ - eyeXYZ(fr,:);

xAxisVec = [eyeRadius 0 0];

%     rotFovToX= createRotationVector3d(foveaXYZ_z, xAxisVec); %returns a 4x4 transformation matrix that puts the foveaXYZ_z vector onto the X-axis
%     rotFovToX = rotFovToX(1:3,1:3); % drop everything but the 3x3 bit
%
%
rotFovToX = rotXToGz^-1;

foveaXYZ_rot = rotFovToX*foveaXYZ_z';
assert(sum(foveaXYZ_rot - xAxisVec')<1e-10) %% make sure the rotation matrix works!

backEyeProjPoints_xyz_id_rot = nan(size(backEyeProjPoints_xyz_id_z));

for vv = 1:length(backEyeProjPoints_xyz_id_z)
    backEyeProjPoints_xyz_id_rot(vv,1:3) = rotFovToX*backEyeProjPoints_xyz_id_z(vv,1:3)';
    backEyeProjPoints_xyz_id_rot(vv,4) = backEyeProjPoints_xyz_id_z(vv,4);
end
% Again, do the same thing to the heading
back_proj_heading_rot = rotFovToX*back_proj_heading';

%%%%% get 2D polar projections of backEyeProjPoints
[backEyeTheta, backEyeRho] = cart2pol(backEyeProjPoints_xyz_id_rot(:,2), backEyeProjPoints_xyz_id_rot(:,3));
backEyeID = backEyeProjPoints_xyz_id_rot(:,4);
%%%%% get 2D polar projections of back_proj_heading
[backEyeTheta_fix, backEyeRho_fix] = cart2pol(back_proj_heading_rot(2), back_proj_heading_rot(3));

% %%% calculate "great circle distance" of each point, to avoid
% underestimating eccentricity due to sphere-plane projection
% convert projected points to spherical cooridinate
%% 5.23 comment out those lines 
[backEyeProjPointsAz, backEyeProjPointsEl, backEyeProjPointsR] = cart2sph(backEyeProjPoints_xyz_id_z(:,1),backEyeProjPoints_xyz_id_z(:,2), backEyeProjPoints_xyz_id_z(:,3));    % find distance between each point and (0,0,r), i.e. foveaXYZ - formula
[fovAz, fovEl, fovR] = cart2sph(foveaXYZ_z(1),foveaXYZ_z(2),foveaXYZ_z(3));
%heading
[backEyeProjPointsAz_fix, backEyeProjPointsEl_fix, backEyeProjPointsR_fix] = cart2sph(back_proj_heading(1),back_proj_heading(2), back_proj_heading(3));  

% equation from https://en.wikipedia.org/wiki/Great-circle_distance,taking
% point1 to be (0,0), phi to be elevation, and lambda to be azimuth
% test1: acos(sin(0)*sin(0) + cos(0)*cos(pi/2).*cos(pi/2)).*eyeRadius = (2*pi*eyeRadius)/4
% test2: acos(sin(0)*sin(0) + cos(0)*cos(pi).*cos(0)).*eyeRadius = (2*pi*eyeRadius)/2

nums  = sin(fovEl).*sin(backEyeProjPointsEl) + cos(fovEl).*cos(backEyeProjPointsEl).*cos(fovAz-backEyeProjPointsAz);
%heading
nums_heading = sin(fovEl)*sin(backEyeProjPointsEl_fix) + cos(fovEl)*cos(backEyeProjPointsEl_fix)*cos(fovAz-backEyeProjPointsAz_fix);

nums(abs(nums)>=1) = nan;
greatCircleDist = acos(nums).*eyeRadius;
%heading
greatCircleDist_fix = acos(nums_heading).*eyeRadius;
retinaProjPoints_fr_theta_rho_id= [backEyeTheta greatCircleDist backEyeID];

[fixX,fixY] = pol2cart(backEyeTheta_fix, greatCircleDist_fix);
% [X, Y] = pol2cart(backEyeTheta, greatCircleDist); in pixels 5/26
X = rad2deg(greatCircleDist.*cos(backEyeTheta)/eyeRadius)/2;
Y = rad2deg(greatCircleDist.*sin(backEyeTheta)/eyeRadius)/2;
retinaProjPoints_fr_X_Y_id = [X Y backEyeID];
headX = rad2deg(greatCircleDist_fix*cos(backEyeTheta_fix)/eyeRadius)/2;
headY = rad2deg(greatCircleDist_fix*sin(backEyeTheta_fix)/eyeRadius)/2;
%keyboard
%% replace the projection method with a more obvious one
% fixX = back_proj_heading_rot(2);
% fixY = back_proj_heading_rot(3);
% heading_deg_x = atand(fixX/eyeRadius + fixX.*back_proj_heading_rot(1)./(eyeRadius*(eyeRadius - back_proj_heading_rot(1))))
% heading_deg_y = atand(fixY/eyeRadius + fixY.*back_proj_heading_rot(1)./(eyeRadius*(eyeRadius - back_proj_heading_rot(1))))
% %% if instead we want things in degrees instaed of in mm
% 
% theta_x = atand(backEyeProjPoints_xyz_id_rot(:,2)/eyeRadius + backEyeProjPoints_xyz_id_rot(:,2).*backEyeProjPoints_xyz_id_rot(:,1)./(eyeRadius*(eyeRadius - backEyeProjPoints_xyz_id_rot(:,1))));
% theta_y = atand(backEyeProjPoints_xyz_id_rot(:,3)/eyeRadius + backEyeProjPoints_xyz_id_rot(:,3).*backEyeProjPoints_xyz_id_rot(:,1)./(eyeRadius*(eyeRadius - backEyeProjPoints_xyz_id_rot(:,1))));
% retinaProjPoints_fr_X_Y_id = [theta_x theta_y backEyeID];
% keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%debug plots!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if debug
    
    figure(334);clf
    subplot(1,3,1)
    hold on;
    axis equal;
    xlim([-eyeRadius*1.2 eyeRadius*1.2])
    ylim([-eyeRadius*1.2 eyeRadius*1.2])
    zlim([-eyeRadius*1.2 eyeRadius*1.2])
    view(3)
    box on
    plot3(backEyeProjPoints_xyz_id_z(:,1), backEyeProjPoints_xyz_id_z(:,2), backEyeProjPoints_xyz_id_z(:,3),'mo','MarkerFaceColor','m')
    hold on
    plot3(foveaXYZ_z(1), foveaXYZ_z(2), foveaXYZ_z(3),'kp','MarkerFaceColor','y')
    
    plot3(backEyeProjPoints_xyz_id_rot(:,1), backEyeProjPoints_xyz_id_rot(:,2), backEyeProjPoints_xyz_id_rot(:,3),'co','MarkerFaceColor','c')
    plot3(foveaXYZ_rot(1), foveaXYZ_rot(2), foveaXYZ_rot(3),'yp','MarkerFaceColor','y','MarkerSize',16)
    
    plot3(-eyeRadius,0,0, 'ko','MarkerFaceColor','g')% fovea
    plot3(eyeRadius,0,0, 'ko','MarkerFaceColor','r') % pupil
    s1 = drawSphere([0 0 0 eyeRadius]);
    s1.FaceColor = [1 1 1]*.5;
    
    %%% plot XYZ (RGB) axes
    plot3([0 eyeRadius*1.5], [0 0], [ 0 0], 'r:','LineWidth',3)
    plot3([0 0], [0 eyeRadius*1.5], [ 0 0], 'g:','LineWidth',3)
    plot3([0 0], [0 0], [ 0 eyeRadius*1.5], 'b:','LineWidth',3)
    plot3(-[0 eyeRadius*1.5], [0 0], [ 0 0], 'r:','LineWidth',3)
    plot3([0 0], -[0 eyeRadius*1.5], [ 0 0], 'g:','LineWidth',3)
    plot3([0 0], [0 0], -[ 0 eyeRadius*1.5], 'b:','LineWidth',3)
    
    %%% draw eye axes
    drawVector3d([0 0 0] , eyeXvec*eyeRadius*1.5,'Color','r','LineWidth',3)
    drawVector3d([0 0 0] , eyeYvec*eyeRadius*1.5,'Color','g','LineWidth',3)
    drawVector3d([0 0 0] , eyeZvec*eyeRadius*1.5,'Color','b','LineWidth',3)
    drawVector3d([0 0 0] , -eyeXvec*eyeRadius*1.5,'Color','r','LineWidth',3)
    drawVector3d([0 0 0] , -eyeYvec*eyeRadius*1.5,'Color','g','LineWidth',3)
    drawVector3d([0 0 0] , -eyeZvec*eyeRadius*1.5,'Color','b','LineWidth',3)
    
    axis equal
    lighting gouraud
    view(50,21)
    
    subplot(1,3,2)
    polarplot(retinaProjPoints_fr_theta_rho_id(:,1), greatCircleDist,'ro','MarkerFaceColor','r', 'MarkerSize',4,'DisplayName','Rho as Great Circle Dist')
    hold on
    polarplot(backEyeTheta, backEyeRho,'ko','MarkerFaceColor','k', 'MarkerSize',4,'DisplayName','Sphere Plane Projection')
    
    legend
    
    subplot(1,3,3)
    plot(X, Y,'ko','MarkerFaceColor','k', 'MarkerSize',4,'DisplayName','Cartesian conversion')
    hold on
    m = max(max(greatCircleDist));
    plot([-m m],[0 0],'r','HandleVisibility','off')
    plot([0 0],[-m m],'r','HandleVisibility','off')
    viscircles([0 0],m);
    axis equal
    legend
    
    if spotcheck
        dbstack
        keyboard
    end
end