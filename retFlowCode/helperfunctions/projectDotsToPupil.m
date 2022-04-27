function [backEyeProjPoints_xyz_id,back_proj_heading] = projectDotsToPupil(thisEyeStruct, fr, debug, spotcheck)
% updated Mar 5 by Mengjian Hua 
% added some code to compute the projection of the fixation point on the
% back of the eye
varNames = fieldnames(thisEyeStruct);
for i=1:length(varNames)
    eval([varNames{i} '=thisEyeStruct.' varNames{i} ';']);
end

%%
% "project" points in point cloud onto back of cone by finding intersection
% of line from pupil through point with the back of cone

groundPoints_xyz_id_proj_lines = createLine3d(groundPoints_xyz_id(:, 1:3), pupilXYZ);
groundPoints_xyz_id_proj = intersectLinePlane(groundPoints_xyz_id_proj_lines, groundPlane);

% dotsInFov = inpolygon(groundPoints_xyz_id(:,2), groundPoints_xyz_id(:,3), fovGroundPoints(:,2), fovGroundPoints(:,3)); %5/23 very important!!!!
dotsInFov = inpolygon(groundPoints_xyz_id_proj(:,2), groundPoints_xyz_id_proj(:,3), fovGroundPoints(:,2), fovGroundPoints(:,3));
projPoint_id = groundPoints_xyz_id(dotsInFov,4); %id# of the points that are being projected on to the eye

% if(fr == 50)
%     keyboard;
% end
%%% debug plots
% plot(groundPoints_xyz_id(:,2),groundPoints_xyz_id(:,3),".");
% hold on
% plot(groundPoints_xyz_id(dotsInFov,2),groundPoints_xyz_id(dotsInFov,3),".");

%%%find points relevant to the occlusion detection, i.e. dots inthe
%%% FOV + dots between the eye positoin and the start of the FOV
%%% (if you just look at the FOV dots, you'll get dots "peeking
%%% out" from underneath. If you look at the whole groundplane, it
%%% takes foreeeeeever to calculate)
% eyeBoundx = [fovGroundPoints(:,1); pupilXYZ(1)];
% eyeBoundy = [fovGroundPoints(:,2); pupilXYZ(2)];
% 
% [c, idx] = unique([ eyeBoundx eyeBoundy],'rows','stable'); %to get rid of annoying "duplicate datapoints detected" message
% 
% [eyeBoundTRI] = delaunayTriangulation(eyeBoundx(idx), eyeBoundy(idx));
% eyeBoundConvHull = convexHull(eyeBoundTRI);  %works better than  "freeBoundary(eyeBoundTRI)." Don't ask me why, or what either of those things means...

% occMeshDots = inpolygon(groundPoints_xyz_id(:,1), groundPoints_xyz_id(:,2), eyeBoundTRI.Points(eyeBoundConvHull,1), eyeBoundTRI.Points(eyeBoundConvHull,2));%the dots defining the mesh that we'll use to determine occlusions
% %         %%%debug
%         triplot(eyeBoundTRI);hold on
%         plot(eyeBoundTRI.Points(eyeBoundConvHull,1), eyeBoundTRI.Points(eyeBoundConvHull,2),'-or','LineWidth',2)
%         plot(groundPoints_xyz_id(occMeshDots,1), groundPoints_xyz_id(occMeshDots,2),'ro')

projectionLines = createLine3d(groundPoints_xyz_id(dotsInFov,1:3), pupilXYZ); %full infinite lines (for intersectLineSphere)

%     projectionEdges = createEdge3d(groundPoints_xyz_id(dotsInFov,1:3), pupilXYZ); % 'edges' (i.e. line segments) only between pupil and dots (to be used to clean up back projected points found by intersectLineSphere)


%create mesh for dots in FOV using delaunay triangulation
% groundMesh.Vertices =(groundPoints_xyz_id(occMeshDots,[1:3]));
% 
% groundMesh.Faces  = delaunay(groundMesh.Vertices (:,[3 2])); % generate the mesh based only on the 2d projection, to avoid spurious connections between non-adjacent dots
% 
% inters = cell(length(projectionLines),1);
% occluded = false(length(projectionLines),1);%preallocating logicals

%% find projection lines that intersect with ground mesh before hitting their dot (i.e. find occlusions)

% if sum(abs(zz(:))) ~= 0 %skip the occulsion detector for flat ground (saves a LOT of time)
%     
%     for gg = 1:length(projectionLines(:,1))
%         [inters{gg}] = intersectLineMesh3d(projectionLines(gg,:), groundMesh.Vertices , groundMesh.Faces);
%         
%         thisPoint = groundPoints_xyz_id(projPoint_id(gg),1:3); % the dotto in question
%         thisPointDist = pdist([pupilXYZ; thisPoint]); %distance of this dot from the pupil (occlusions must be *closer* to pupil than dotto)
%         
%         dropThese = round(inters{gg},10) == round(thisPoint,10); %remove intersections that are the same as the dot (i.e. non-occlusion intersections)
%         
%         inters{gg}(dropThese(:,1),:) = [];
%         
%         occDists = sqrt( (pupilXYZ(1)-inters{gg}(:,1)).^2 + (pupilXYZ(2)-inters{gg}(:,2)).^2 + (pupilXYZ(3)-inters{gg}(:,3)).^2);
%         
%         dropTheseToo = occDists>thisPointDist;
%         
%         inters{gg}(dropTheseToo,:) = []; %remove intersections farther from pupil than dot in question (i.e. intersections *behind* the dot aren't occlusions)
%         occluded(gg) = ~isempty(inters{gg});
%         
%     end
%     
%     occIntersections = cell2mat(inters(occluded));
%     allPointsInFOV = groundPoints_xyz_id(dotsInFov,1:3);
%     occPoints = allPointsInFOV(occluded,:);
%     
%     occLines = projectionLines(occluded,:) ;
%     projectionLines(occluded,:) = []; %remove projection lines with occluded dots
%     projPoint_id(occluded) = [];
% end
%% find projection of ground grid onto eyeball

eyeProjIntersections_xyz_id = nan(length(projectionLines(:,1))*2, 4);

for ll = 1:length(projectionLines)%doing this in a loop in order to make it easier to maintain the ID of the projected points intersections

    eyeProjIntersections_xyz_id(ll*2-1:ll*2,1:3) = intersectLineSphere(projectionLines(ll,:), eyeSphere); %find intersections between projected line and eyeball sphere
    eyeProjIntersections_xyz_id(ll*2-1:ll*2,4) = projPoint_id(ll); %save the ID of the ground points that caused these intersections
    
end
%% find the projection of the heading vector onto eyeball

project_line_heading = createLine3d(heading, pupilXYZ);
two_intersection_points = intersectLineSphere(project_line_heading, eyeSphere);

%% Identify which intersection points are on the *back* of the eye -
%%%% defined as "points that are farther from the fixation point than the average betwee the Pupil and the eye center, i.e. the back 75% of the eye (USED TO BE - mid eye point (eyeXYZ))

thisFixPoint = fixXYZ(fr,:);
% thisEyeXYZ = eyeStruct.(thisEyeD).eyeXYZ; %setting this to be the eye
% center limits the potential FOV substantially. Replacing it with
% 'pupilXYZ', since the retina extevnds beyond the eye midline
thisEyeXYZ = mean([pupilXYZ; eyeXYZ(fr,:)]); %setting this to be the midpoint between the eye center and the pupil

projPointDists = sqrt( (thisFixPoint(1) - eyeProjIntersections_xyz_id(:,1)).^2 +...
    (thisFixPoint(2) - eyeProjIntersections_xyz_id(:,2)).^2 +...
    (thisFixPoint(3) - eyeProjIntersections_xyz_id(:,3)).^2);

eyeFixDist = sqrt( (thisFixPoint(1) - thisEyeXYZ(1)).^2 +...
    (thisFixPoint(2) - thisEyeXYZ(2)).^2 +...
    (thisFixPoint(3) - thisEyeXYZ(3)).^2);

projPointDists_fix = sqrt( (thisFixPoint(1) - eyeProjIntersections_xyz_id(1))^2 +...
    (thisFixPoint(2) - eyeProjIntersections_xyz_id(2))^2 +...
    (thisFixPoint(3) - eyeProjIntersections_xyz_id(3))^2);

projPointDists_heading = sqrt( (thisFixPoint(1) - two_intersection_points(:,1)).^2 +...
    (thisFixPoint(2) - two_intersection_points(:,2)).^2 +...
    (thisFixPoint(3) - two_intersection_points(:,3)).^2);

backEyeProjPoints_xyz_id = eyeProjIntersections_xyz_id(projPointDists > eyeFixDist, :);

back_proj_heading = two_intersection_points(projPointDists_heading > eyeFixDist, :); % not sure about this
if length(back_proj_heading) == 0 % CSB added, if heading is backwards, use other intersection: back of eye ball
    back_proj_heading = two_intersection_points(projPointDists_heading < eyeFixDist, :); % not sure about this
end

%% %plot it up!
%debug = true; %uncomment this to show debug plot
if debug 
    figure(5672); clf
    hold on
    axis equal
    drawSphere(eyeSphere,'FaceColor',[.5 .5 .5], 'FaceAlpha',.8); %eyeball, at origin
    drawSphere([pupilXYZ 30],'FaceColor','k'); %pupil
    plot3(thisEyeXYZ(1), thisEyeXYZ(2), thisEyeXYZ(3),'kp','MarkerSize', 12,'MarkerFaceColor','y')
    
    
    %%%%%draw fixated point
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
    gr1 = scatter3(groundPoints_xyz_id(dotsInFov,1),groundPoints_xyz_id(dotsInFov,2),groundPoints_xyz_id(dotsInFov,3),'filled');
    gr1.Marker = 'o';
    gr1.MarkerFaceColor = 'c';
    hold on
    
    
    
    %%%plot ground plane
    g1 = surface(xx, yy, zz);
    g1.EdgeColor =  'k';
    g1.EdgeAlpha = .5;
    g1.FaceColor = [.4 .8 .6];
    g1.FaceAlpha = 1;
    
    if exist('occLines')
        %%draw projection lines
        %                 ol1 = drawLine3d(occLines);
        
        %%%draw intersections
        i1 = scatter3(occIntersections(:,1), occIntersections(:,2),occIntersections(:,3),'filled');
        i1.MarkerFaceColor = 'm';
        i1.SizeData = 12;
        
        %%%draw occluded p9oints
        o1 = scatter3(occPoints(:,1), occPoints(:,2),occPoints(:,3),'filled');
        o1.MarkerFaceColor = 'r';
    end
    
    plot3(fovGroundPoints(:,1), fovGroundPoints(:,2), fovGroundPoints(:,3),'k.-')
    
    xlim([min(fovGroundPoints(:,1)) max(fovGroundPoints(:,1))])
    ylim([min(fovGroundPoints(:,2)) max(fovGroundPoints(:,2))])
    zlim([-1e3 max(eyeXYZ(:,3))+1e3])
    
    view(3)
    drawnow
    
    if spotcheck
        dbstack
        keyboard
    end
    
end