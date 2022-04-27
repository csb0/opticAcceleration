function [eyeMesh,...
    eyeSphere,...
    thisFixPoint,...
    gazeLine,...
    pupilXYZ,...
    foveaXYZ,...
    rotXToGzEul,...
    rotXToGz,...
    eyeXvec,...
    eyeYvec,...
    eyeZvec,...
    coronalEyePlane,...
    saggitalEyePlane,...
    transverseEyePlane ,...
    coronalEyeCircle,...
    saggitalEyeCircle,...
    transverseEyeCircle,...
    eyeCardinalPointNorth,...
    eyeCardinalPointSouth,...
    eyeCardinalPointEast,...
    eyeCardinalPointWest,...
    saggitalEyeCirclePoints, ...
    coronalEyeCirclePoints,...
    transverseEyeCirclePoints,...
    verticalEyeAxisGroundIntersect,...
    horizontalEyeAxisGroundIntersect] =...
    defineEyeball(thisEyeStruct,  fr, debug, spotcheck)


varNames = fieldnames(thisEyeStruct);
for i=1:length(varNames)
    eval([varNames{i} '=thisEyeStruct.' varNames{i} ';']);
end

thisEyeXYZ = thisEyeStruct.eyeXYZ;

eyeX = thisEyeXYZ(fr,1);
eyeY = thisEyeXYZ(fr,2);
eyeZ = thisEyeXYZ(fr,3);

eyeSphere = [eyeX eyeY eyeZ eyeRadius];

[v, f] = sphereMesh(eyeSphere);
f = triangulateFaces(f);
eyeMesh.faces = f;
eyeMesh.vertices = v;



thisFixPoint = fixXYZ(fr,:);

%% define "pupil" and "fovea" as the point on the front and back of eye sphere (respectively) that intersect with line between fixated point and eye center

%%%%%define line connecting eyeball center and fixated point
gazeLine = createLine3d(thisEyeXYZ(fr,:), thisFixPoint);


gazeEyeIntersections = intersectLineSphere(gazeLine, eyeSphere);


gazeEyeIntersectionDists = sqrt( (thisFixPoint(1) - gazeEyeIntersections(:,1)).^2 +...
    (thisFixPoint(2) - gazeEyeIntersections(:,2)).^2 +...
    (thisFixPoint(3) - gazeEyeIntersections(:,3)).^2);

assert(length(gazeEyeIntersectionDists) == 2,'Gaze intersects eyeballs at more (or less) than 2 points! Your eyeballs are screwy!!') %% if there are not 2 points, somethings screwy with your eyeballs!!

pupilXYZ = gazeEyeIntersections(gazeEyeIntersectionDists < mean(gazeEyeIntersectionDists),:);%the pupil is closer to the fixation point than the fovea
foveaXYZ = gazeEyeIntersections(gazeEyeIntersectionDists > mean(gazeEyeIntersectionDists),:);%cringing as I name this point "foveaXYZ," but I mean... what else would I call it?!



%%%
%% define coronal(midEyePlane, splitting Front and Back), saggital (splitting left and right), and transverse planes of the eye (splitting top and bottom)

%begin with the eye at the origin, pointing towards the "fixation point"
gaze_z = pupilXYZ - thisEyeXYZ(fr,:);

%%%% define some world coordinates
origin = [ 0 0 0];
xUnit = [1 0 0]; %x unit vector
yUnit = [0 1 0]; %y unit vector
zUnit = [0 0 1]; %z unit vector

xyPlane = createPlane(origin, xUnit, yUnit); %the XY plane contains the origin, the xUnit vector and the yUnit vector (could also have defined it using zUnit as the normal vector to the origin)
xzPlane = createPlane(origin, xUnit, zUnit); %etc
yzPlane = createPlane(origin, yUnit, zUnit); %etc

%%%get rotation matrix that moves xUnit vector onto gaze line
%         rotXToGz= createRotationVector3d(xUnit, gaze_z); %returns a 4x4 transformation matrix that puts the foveaXYZ_z vector onto the X-axis
%         rotXToGz = rotXToGz(1:3,1:3); % drop everything but the 3x3 bit
%
[xAz, xEl, xRho] = cart2sph(xUnit(1), xUnit(2), xUnit(3));
[gAz, gEl, gRho] = cart2sph(gaze_z(1), gaze_z(2), gaze_z(3));

rotXToGzEul = [gAz, -gEl, 0];
rotXToGz = eul2rotm(rotXToGzEul); %if you want to add torsion, replace that 0 with your desired torsion value (radians)


%%%% define coronal plane by gaze vector that cuts eye in half between front and back (to better define the "retina" )
eyeXvec = [rotXToGz * xUnit']';
coronalEyePlane = createPlane(origin, eyeXvec);
coronalEyePlane(1:3) = thisEyeXYZ(fr,:);
coronalEyeCircle = intersectPlaneSphere(coronalEyePlane, eyeSphere); %a circle defining the coronal plane of the eye
coronalEyeCirclePoints = circleToPoints3d(coronalEyeCircle);

%%%% define saggital plane by gaze vector that cuts eye in half between left and right
eyeYvec = [rotXToGz* yUnit']';
saggitalEyePlane = createPlane(origin, eyeYvec);
saggitalEyePlane(1:3) = thisEyeXYZ(fr,:);
saggitalEyeCircle = intersectPlaneSphere(saggitalEyePlane, eyeSphere); %a circle defining the transverse plane of the eye
saggitalEyeCirclePoints = circleToPoints3d(saggitalEyeCircle);

%%%% define transverse plane by gaze vector that cuts eye in half between top and bottom
eyeZvec = [rotXToGz* zUnit']';
transverseEyePlane = createPlane(origin, eyeZvec);
transverseEyePlane(1:3) = thisEyeXYZ(fr,:);
transverseEyeCircle = intersectPlaneSphere(transverseEyePlane, eyeSphere); %a circle defining the transverse plane of the eye
transverseEyeCirclePoints = circleToPoints3d(transverseEyeCircle);




%%%% define eye cardinal points
eyeCardinalPointNorth = eyeZvec *eyeRadius;
eyeCardinalPointSouth = -eyeZvec *eyeRadius;
eyeCardinalPointEast = eyeYvec *eyeRadius;
eyeCardinalPointWest = -eyeYvec *eyeRadius;


        %                 Get ground proj points for eye crosshairs
        
        %use edge ([x1 y1 z1 x2 y2 z2]) instead of infinite lines, b/c intersectLinePlane will find intersections *behind* eyeball
        for ll = 1:size(saggitalEyeCirclePoints)
            verticalEyeAxisProjEdge(ll,:) = [saggitalEyeCirclePoints(ll,:) saggitalEyeCirclePoints(ll,:)-pupilXYZ];            
            horizontalEyeAxisProjEdge(ll,:) = [transverseEyeCirclePoints(ll,:) transverseEyeCirclePoints(ll,:)-pupilXYZ];
        end

        verticalEyeAxisGroundIntersect= intersectLinePlane(verticalEyeAxisProjEdge, groundPlane);
        horizontalEyeAxisGroundIntersect = intersectLinePlane(horizontalEyeAxisProjEdge, groundPlane);

        %% 

        for dd = 1:3
            verticalEyeAxisGroundIntersect1(:,dd) = interp1(verticalEyeAxisGroundIntersect(:,dd),...
                1:.01:length(verticalEyeAxisGroundIntersect));
            horizontalEyeAxisGroundIntersect1(:,dd) = interp1(horizontalEyeAxisGroundIntersect(:,dd),...
                1:.01:length(horizontalEyeAxisGroundIntersect));
        end
        verticalEyeAxisGroundIntersect = verticalEyeAxisGroundIntersect1 ;
        horizontalEyeAxisGroundIntersect = horizontalEyeAxisGroundIntersect1;
        
        verticalEyeAxisGroundIntersect([1 end],:) = [];
        horizontalEyeAxisGroundIntersect1([1 end],:) = [];
        %%
assert(sum(dot(eyeYvec, eyeZvec))+ sum(dot(eyeYvec, eyeZvec)) + sum(dot(eyeYvec, eyeZvec)) <1e-10)%make sure eye vectors are orthogonal (i.e. the dot product between them is 0ish)

if debug
    figure(743);clf
    
    subplot(1,2,1);
    hold on
    view(3)
    axis equal
    drawSphere([0 0 0 eyeRadius],'FaceColor',[.5 .5 .5], 'FaceAlpha',.8); %eyeball, at origin
    drawSphere([pupilXYZ - thisEyeXYZ(fr,:) eyeRadius*.1],'FaceColor','k'); %pupil
    plot3(0,0,0,'kp','MarkerSize', 12,'MarkerFaceColor','y')
    
    drawPlane3d(xyPlane,'FaceColor','w', 'FaceAlpha',.6)
    drawPlane3d(xzPlane,'FaceColor','w', 'FaceAlpha',.6)
    drawPlane3d(yzPlane,'FaceColor','w', 'FaceAlpha',.6)
    
    drawVector3d(origin,eyeXvec*eyeRadius*1.5,'Color','m','LineWidth',3)
    drawPlane3d([origin coronalEyePlane(4:end)],'FaceColor','m', 'FaceAlpha',.3)
    drawCircle3d([0 0 0 coronalEyeCircle(4:end)],'LineWidth', 5,'Color', 'm');%coronal eye circle, at origin
    
    drawVector3d(origin, eyeYvec*eyeRadius*1.5,'Color','y','LineWidth',3)
    drawPlane3d([origin  saggitalEyePlane(4:end)],'FaceColor','y', 'FaceAlpha',.3)
    drawCircle3d([0 0 0  saggitalEyeCircle(4:end)],'LineWidth', 5,'Color', 'y');%Cyan eye circle, at origin
    
    drawVector3d(origin, eyeZvec*eyeRadius*1.5,'Color','c','LineWidth',3)
    drawPlane3d([origin transverseEyePlane(4:end)],'FaceColor','c', 'FaceAlpha',.3)
    drawCircle3d([0 0 0 transverseEyeCircle(4:end)],'LineWidth', 5,'Color', 'c');%Cyan eye circle, at origin
    
    
    drawVector3d(origin,xUnit*eyeRadius*1.5,'Color','r','LineWidth',3)
    drawVector3d(origin,yUnit*eyeRadius*1.5,'Color','g','LineWidth',3)
    drawVector3d(origin,zUnit*eyeRadius*1.5,'Color','b','LineWidth',3)
    
    drawSphere([eyeCardinalPointNorth eyeRadius*.025],'FaceColor','b'); %pupil
    drawSphere([eyeCardinalPointEast eyeRadius*.025],'FaceColor','y'); %pupil
    drawSphere([eyeCardinalPointSouth eyeRadius*.025],'FaceColor','c'); %pupil
    drawSphere([eyeCardinalPointWest eyeRadius*.025],'FaceColor','g'); %pupil
    
    plot3(saggitalEyeCirclePoints(:,1)-thisEyeXYZ(fr,1), saggitalEyeCirclePoints(:,2)-thisEyeXYZ(fr,2), saggitalEyeCirclePoints(:,3)-thisEyeXYZ(fr,3),'ko')
    plot3(coronalEyeCirclePoints(:,1)-thisEyeXYZ(fr,1), coronalEyeCirclePoints(:,2)-thisEyeXYZ(fr,2), coronalEyeCirclePoints(:,3)-thisEyeXYZ(fr,3),'ko')
    plot3(transverseEyeCirclePoints(:,1)-thisEyeXYZ(fr,1), transverseEyeCirclePoints(:,2)-thisEyeXYZ(fr,2), transverseEyeCirclePoints(:,3)-thisEyeXYZ(fr,3),'ko')
    
%         plot3(saggitalEyeCirclePoints(1:35,1)-thisEyeXYZ(fr,1), saggitalEyeCirclePoints(1:35,2)-thisEyeXYZ(fr,2), saggitalEyeCirclePoints(1:35,3)-thisEyeXYZ(fr,3),'ro','MarkerSize',30)

    view(45, 30)
    
    
    
    subplot(122)
    
    lLeg = [2 3 4 5 6 7 5];
    rLeg = [2 8 9 10 11 12 10];
    tors = [2 13 14 15 26 27 28];
    lArm = [15 16 17 26 17 18 19 20];
    rArm = [15 21 22 26 22 23 24 25];
    
    %%%Plotcherself up a nice little skeleetoon friend
    plot3(shadow_fr_mar_dim(fr,1:28,1),shadow_fr_mar_dim(fr,1:28,2),shadow_fr_mar_dim(fr,1:28,3),'ko','MarkerFaceColor','k','MarkerSize',4)
    hold on
    
    
    plot3(shadow_fr_mar_dim(fr,lLeg,1),shadow_fr_mar_dim(fr,lLeg,2),shadow_fr_mar_dim(fr,lLeg,3),'c','LineWidth',2)
    plot3(shadow_fr_mar_dim(fr,rLeg,1),shadow_fr_mar_dim(fr,rLeg,2),shadow_fr_mar_dim(fr,rLeg,3),'r','LineWidth',2)
    plot3(shadow_fr_mar_dim(fr,tors,1),shadow_fr_mar_dim(fr,tors,2),shadow_fr_mar_dim(fr,tors,3),'g','LineWidth',2)
    plot3(shadow_fr_mar_dim(fr,lArm,1),shadow_fr_mar_dim(fr,lArm,2),shadow_fr_mar_dim(fr,lArm,3),'c','LineWidth',2)
    plot3(shadow_fr_mar_dim(fr,rArm,1),shadow_fr_mar_dim(fr,rArm,2),shadow_fr_mar_dim(fr,rArm,3),'r','LineWidth',2)
    
    plot3(eyeXYZ(fr,1),eyeXYZ(fr,2), eyeXYZ(fr,3),'kp','MarkerFaceColor','r')
    plot3(pupilXYZ(1), pupilXYZ(2), pupilXYZ(3),'kp','MarkerFaceColor','b')
    plot3(thisFixPoint(1), thisFixPoint(2), thisFixPoint(3),'kp','MarkerFaceColor','m')
    %
    drawSphere(eyeSphere,'FaceColor',[.5 .5 .5], 'FaceAlpha',.8); %eyeball, at origin
    drawSphere([pupilXYZ eyeRadius*.1],'FaceColor','k'); %pupil
    
    
    
    drawVector3d(eyeXYZ(fr,:),eyeXvec*eyeRadius*5,'Color','m','LineWidth',3)
    drawCircle3d(coronalEyeCircle,'LineWidth', 5,'Color', 'm');%coronal eye circle, at origin
    
    drawVector3d(eyeXYZ(fr,:), eyeYvec*eyeRadius*5,'Color','y','LineWidth',3)
    drawCircle3d(saggitalEyeCircle,'LineWidth', 5,'Color', 'y');%Cyan eye circle, at origin
    
    drawVector3d(eyeXYZ(fr,:), eyeZvec*eyeRadius*5,'Color','c','LineWidth',3)
    drawCircle3d(transverseEyeCircle,'LineWidth', 5,'Color', 'c');%Cyan eye circle, at origin
    
    
    drawVector3d(eyeXYZ(fr,:),xUnit*eyeRadius*5,'Color','r','LineWidth',3)
    drawVector3d(eyeXYZ(fr,:),yUnit*eyeRadius*5,'Color','g','LineWidth',3)
    drawVector3d(eyeXYZ(fr,:),zUnit*eyeRadius*5,'Color','b','LineWidth',3)
    
    
    
    drawLine3d(gazeLine)
    axis equal
    view(45, 30)
    
    
    
    if spotcheck
        dbstack
        keyboard;
    end
    
end


