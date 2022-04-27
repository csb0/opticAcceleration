function  [eyeStruct] = defineHeadAndEyeTraj(eyeStruct, interocularDistance, debug)

%% blow apart eyeStruct into its component variables
varNames = fieldnames(eyeStruct);
for i=1:length(varNames)
    eval([varNames{i} '=eyeStruct.' varNames{i} ';']);
end

%%
%%takes in a 3D head trajectory and fixation point(s) and spits out the
%%paths of two eyes attached to said head trajectory, pointing at the
%%fixation points. Assumes the head and eyes are both pointing at the smae
%%place (for now).

rEye = [];
lEye = [];


for fr = 1:length(headXYZ)
    
    
    headX = headXYZ(fr,1);
    headY = headXYZ(fr,2);
    headZ = headXYZ(fr,3);
    
    headRadius = interocularDistance/2;
    
    headSphere = [headX headY headZ headRadius];
    
    
    thisFixPoint = fixXYZ(fr,:);
    
    %%%%%define line connecting eyeball center and fixated point
    gazeLine = createLine3d(headXYZ(fr,:), thisFixPoint);
    
    
    
    gazeHeadIntersections = intersectLineSphere(gazeLine, headSphere);
    
    
    gazeEyeIntersectionDists = sqrt( (thisFixPoint(1) - gazeHeadIntersections(:,1)).^2 +...
        (thisFixPoint(2) - gazeHeadIntersections(:,2)).^2 +...
        (thisFixPoint(3) - gazeHeadIntersections(:,3)).^2);
    
    assert(length(gazeEyeIntersectionDists) == 2,'Gaze intersects headball at more (or less) than 2 points! Your headball is screwy!!') %% if there are not 2 points, somethings screwy with your eyeballs!!
    
    frontHeadXYZ = gazeHeadIntersections(gazeEyeIntersectionDists < mean(gazeEyeIntersectionDists),:);%the pupil is closer to the fixation point than the fovea
    backHeadXYZ = gazeHeadIntersections(gazeEyeIntersectionDists > mean(gazeEyeIntersectionDists),:);%cringing as I name this point "foveaXYZ," but I mean... what else would I call it?!
    
    
    %%%
    %% define rEye and lEye positoins relative to this headXYZ location
    
    %begin with the eye at the origin, pointing towards the "fixation point"
    gaze_z = frontHeadXYZ - headXYZ(fr,:);
    
    %%%% define some world coordinates
    origin = [ 0 0 0];
    xUnit = [1 0 0]; %x unit vector
    yUnit = [0 1 0]; %y unit vector
    zUnit = [0 0 1]; %z unit vector
    
    %%%%define rEye and lEye initial positions
    rEyeInit = [0 interocularDistance/2 0];
    lEyeInit = [0 -interocularDistance/2 0];
    
    xyPlane = createPlane(origin, xUnit, yUnit); %the XY plane contains the origin, the xUnit vector and the yUnit vector (could also have defined it using zUnit as the normal vector to the origin)
    xzPlane = createPlane(origin, xUnit, zUnit); %etc
    yzPlane = createPlane(origin, yUnit, zUnit); %etc
    
    
    [xAz, xEl, xRho] = cart2sph(xUnit(1), xUnit(2), xUnit(3));
    [gAz, gEl, gRho] = cart2sph(gaze_z(1), gaze_z(2), gaze_z(3));
    
    headRotMat(fr,:,:) = eul2rotm([gAz, -gEl, 0]); %if you want to add torsion, replace that 0 with your desired torsion value (radians)
    thisHeadRotMat = squeeze(headRotMat(fr,:,:));
    %%%% rotate rEye intial position about the head center, then add the
    %%%% head position
    rEye_z = [thisHeadRotMat * rEyeInit']';
    rEyeXYZ(fr,:) = rEye_z+headXYZ(fr,:);
    
    
    %%%% ditto for the eye that's left
    lEye_z = [thisHeadRotMat * lEyeInit']';
    lEyeXYZ(fr,:) = lEye_z+headXYZ(fr,:);
    
    %%% rotate yr head vectors
    headXvec = [thisHeadRotMat * xUnit']';
    headYvec = [thisHeadRotMat * yUnit']';
    headZvec = [thisHeadRotMat * zUnit']';
    
    assert( ((pdist([lEyeXYZ(fr,:); rEyeXYZ(fr,:)]) - interocularDistance) + (pdist([rEye_z; lEye_z]) - interocularDistance)) <1e-10 );%make sure eyes remain 1 interocular distance apart
    
    if debug
        if mod(fr,10) == 0
            figure(743);clf
            hold on
            view(3)
            axis equal
            drawSphere([0 0 0 headRadius],'FaceColor',[.5 .5 .5], 'FaceAlpha',.8); %headball, at origin
            drawSphere([frontHeadXYZ - headXYZ(fr,:) 3],'FaceColor','w'); %nose, i guess?
            plot3(0,0,0,'kp','MarkerSize', 12,'MarkerFaceColor','y')
            
            drawSphere([rEyeInit 2],'FaceColor',[.5 0 0 ], 'FaceAlpha',.8); %right eyeInitial
            drawSphere([rEye_z 5],'FaceColor',[1 0 0 ], 'FaceAlpha',.8); %right eyeInitial
            
            drawSphere([lEyeInit 2],'FaceColor',[0 0 0.5 ], 'FaceAlpha',.8); %right eyeInitial
            drawSphere([lEye_z 5],'FaceColor',[0 0 1 ], 'FaceAlpha',.8); %right eyeInitial
            
            drawPlane3d(xyPlane,'FaceColor','w', 'FaceAlpha',.6)
            drawPlane3d(xzPlane,'FaceColor','w', 'FaceAlpha',.6)
            drawPlane3d(yzPlane,'FaceColor','w', 'FaceAlpha',.6)
            
            drawVector3d(origin, headXvec*headRadius*1.5,'Color','m','LineWidth',3)
            drawVector3d(origin, headYvec*headRadius*1.5,'Color','y','LineWidth',3)
            drawVector3d(origin, headZvec*headRadius*1.5,'Color','c','LineWidth',3)
            
            
            drawVector3d(origin,xUnit*headRadius*1.5,'Color','r','LineWidth',3)
            drawVector3d(origin,yUnit*headRadius*1.5,'Color','g','LineWidth',3)
            drawVector3d(origin,zUnit*headRadius*1.5,'Color','b','LineWidth',3)
            
            view(45, 30)
            drawnow
        end
        %         if spotcheck
        %             dbstack
        %             keyboard;
        %         end
    end
    
    %     [eyeBackMesh] = cutMeshByPlane(eyeMesh, midEyePlane,'part','above'); %asking for 'part','above' returns back of the eye?
end

eyeStruct.rEye.eyeXYZ = rEyeXYZ;
eyeStruct.lEye.eyeXYZ = lEyeXYZ;



end

