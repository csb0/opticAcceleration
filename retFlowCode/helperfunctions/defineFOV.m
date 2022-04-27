function [fovRadiusDeg,fovRadiusRad, fovGroundPoints] = defineFOV(fovRadiusDeg, distCutOff, thisEyeStruct, fr,debug, spotcheck)

varNames = fieldnames(thisEyeStruct);
for i=1:length(varNames)
    eval([varNames{i} '=thisEyeStruct.' varNames{i} ';']);
end



%draw a circle centered on the YZ plane centered on the xAxis
fovRadiusRad = deg2rad(fovRadiusDeg); % fov radius in radians



fovRes = 1000;

eastCardInd = 1;
northCardInd = round(fovRes/4);
westCardInd = round(fovRes/2);
southCardInd = round(3*fovRes/4);

xCir = (ones(fovRes,1)*eyeRadius); %draw this circle at one eyeRadius in front of the pupil
yCir = (tan(fovRadiusRad)*eyeRadius) * cos(linspace(0,2*pi,fovRes));
zCir = (tan(fovRadiusRad)*eyeRadius) * sin(linspace(0,2*pi,fovRes));

% This is to make circle behind the front wall
% xCir = xCir*12.5;
% yCir = yCir*12.5;
% zCir = zCir*12.5;

%rotate all the circle points by rotXToGz <<-- rotation matrix to put the xUnit vector onto the gaze line

fovXYZ = nan(fovRes,3);

for cc = 1:length(xCir)
    fovXYZ(cc,:) = rotXToGz * [xCir(cc); yCir(cc); zCir(cc);];
end

fovXYZ(:,1) =     fovXYZ(:,1)*1e4 + pupilXYZ(1);
fovXYZ(:,2) =     fovXYZ(:,2)*1e4 + pupilXYZ(2);
fovXYZ(:,3) =     fovXYZ(:,3)*1e4 + pupilXYZ(3);

% fovXYZ(fovXYZ(:,3)>0,3) = 0; %fovValues aren't allowed to be larger than 0 (i.e. above the groundplane

%use edge ([x1 y1 z1 x2 y2 z2]) instead of infinite lines, b/c intersectLinePlane will find intersections *behind* eyeball
fovProjEdges = [ones(size(fovXYZ(:,1)))*pupilXYZ(1) ones(size(fovXYZ(:,1)))*pupilXYZ(2) ones(size(fovXYZ(:,1)))*pupilXYZ(3) fovXYZ];

fovGroundPoints = intersectEdgePlane(fovProjEdges, groundPlane);
fovGroundPoints(isnan(fovGroundPoints(:,1)),:) = [];%delete non intersecting points
%%%remove FOV points that are more than DistCutOff away ( to avoid plotting
%%%a bagillion dottos when the skel looks towards the horizon)
% [fovTh, fovR] = cart2pol(fovGroundPoints(:,1) - pupilXYZ(1), fovGroundPoints(:,2) - pupilXYZ(2));
% 
% fovR(fovR > distCutOff )  = distCutOff;
% 
% [fovGroundPoints(:,1), fovGroundPoints(:,2)] = pol2cart(fovTh, fovR);
% 
% %  % Mengjian Hua changed here 2/24
% fovGroundPoints(:,1) = fovGroundPoints(:,1) + pupilXYZ(1);
% fovGroundPoints(:,2) = fovGroundPoints(:,2) + pupilXYZ(2);


%% find crosshairs

% [vertCrossHairsIdx,dV] = knnsearch(fovGroundPoints, thisEyeStruct.verticalEyeAxisGroundIntersect);
% v1 = sortrows([dV vertCrossHairsIdx]); %sort by minimum distance
% v2 = unique(v1(:,2),'stable');%remove duplicate indices, keeping 'stable' order from previous sort
% 
% while pdist(fovGroundPoints(v2(1:2,:),:))<1 %make sure you aren't just picking two points on the same side of the fov
%     v2(2,:) = [];
% end
% verticalEyeCrossHairGroundPoints = fovGroundPoints(v2(1:2),:);
% 
% [horizCrossHairsIdx,dH] = knnsearch(fovGroundPoints, thisEyeStruct.horizontalEyeAxisGroundIntersect);
% h1 = sortrows([dH horizCrossHairsIdx]); %sort by minimum distance
% h2 = unique(h1(:,2),'stable');%remove duplicate indices, keeping 'stable' order from previous sort
% 
% while pdist(fovGroundPoints(h2(1:2,:),:))<1 %make sure you aren't just picking two points on the same side of the fov
%     h2(2,:) = [];
% end
% horizontalEyeCrossHairGroundPoints = fovGroundPoints(h2(1:2,:),:);

% 
% %%%debug
% figure(35730);clf
% drawPoint( thisEyeStruct.verticalEyeAxisGroundIntersect(:,1:2));
% hold on 
% drawPoint( thisEyeStruct.verticalEyeAxisGroundIntersect([1 end],1:2),'rp');
% drawPoint( thisEyeStruct.horizontalEyeAxisGroundIntersect(:,1:2));
% drawPoint( thisEyeStruct.horizontalEyeAxisGroundIntersect([1 end],1:2),'gp');
% drawPoint(fovGroundPoints(:,1:2))
% plot(horizontalEyeCrossHairGroundPoints(:,1),horizontalEyeCrossHairGroundPoints(:,2),'rp-','LineWidth',3)
% plot(verticalEyeCrossHairGroundPoints(:,1),verticalEyeCrossHairGroundPoints(:,2),'bp-','LineWidth',3)

%%



if debug
    figure(3029);clf
    hold on
    
    drawSphere(eyeSphere,'FaceColor',[.5 .5 .5], 'FaceAlpha',.8); %eyeball, at origin
    drawSphere([pupilXYZ 30],'FaceColor','k'); %pupil
    plot3(eyeXYZ(fr,1), eyeXYZ(fr,2), eyeXYZ(fr,3),'kp','MarkerSize', 12,'MarkerFaceColor','y')
    
    plot3(fovXYZ(:,1), fovXYZ(:,2), fovXYZ(:,3),'k.-','LineWidth',3)
    
    thisFixPoint = fixXYZ(fr,:);
    fx1 = drawPoint3d(thisFixPoint);
    fx1.Marker = 'p';
    fx1.MarkerSize = 12;
    fx1.MarkerFaceColor = 'k';
    fx1.Color = 'k';
    
    %%% draw gaze line (between eye center and fixated point)
    gl1 = drawLine3d(gazeLine);
    gl1.Color = 'm';
    gl1.LineWidth = 3;
    
    %%%%draw ground points
    gr1 = scatter3(groundPoints_xyz_id(:,1),groundPoints_xyz_id(:,2),groundPoints_xyz_id(:,3),'filled');
    gr1.Marker = 'o';
    gr1.SizeData = 8;
    gr1.CData = dotColors(groundPoints_xyz_id(:,4),:);
    
    %%%plot ground plane
    %[gx, gy] = meshgrid(-1e5:1e3:1e5);
    
   % g1 = surface(gx, gy, zeros(size(gy)));
    %g1.EdgeColor = 'none';
    %g1.FaceColor = 'c';
    
    %%% draw projection lines of FOV and points where they intersect grounplane
    f1 = drawPoint3d(fovGroundPoints);
    f1.MarkerFaceColor = 'k';

    plot(horizontalEyeCrossHairGroundPoints(:,1),horizontalEyeCrossHairGroundPoints(:,2),'rp-','LineWidth',3)
    plot(verticalEyeCrossHairGroundPoints(:,1),verticalEyeCrossHairGroundPoints(:,2),'bp-','LineWidth',3)
    
    xlim([min(fovGroundPoints(:,1))-2e3 max(fovGroundPoints(:,1))+2e3])
    ylim([min(fovGroundPoints(:,2))-2e3 max(fovGroundPoints(:,2))+2e3])
    zlim([-1e3 eyeXYZ(fr,3)+1.5e4]) % Mengjian Hua changed here 2/24
    view(3)
    drawnow
    
    if spotcheck
        dbstack
        
        keyboard
    end
end
