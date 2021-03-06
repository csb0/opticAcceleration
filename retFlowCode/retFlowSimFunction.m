% Title:            retFlowSimFunction.m
%
% Authors:          Oliver Xu, Charlie Burlingham
%
% Purpose:          This code generates a simulation of the optic flow
%                   experienced during locomotion of a cyclopean eye over a
%                   simplified, linear eye trajectory that maintains perfect
%                   fixation on a specified point in the direction of the plane
%
% Inputs:           heading_ang: heading of simulation subject, measured as
%                   degrees from the normal to the plane
%                   depth: depth of the plane (mm)
%                   fix_depth: depth of fixation (mm)
%                   save_dir: path for data
%
% Last updated:     April 27 2022
%
% License:          Copyright (C) 2022 Oliver Xu and Charlie Burlingham
%
%                   This program is free software: you can redistribute it and/or modify
%                   it under the terms of the GNU Affero General Public License as
%                   published by the Free Software Foundation, either version 3 of the
%                   License, or (at your option) any later version.
%
%                   This program is distributed in the hope that it will be useful,
%                   but WITHOUT ANY WARRANTY; without even the implied warranty of
%                   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%                   GNU Affero General Public License for more details.
%
%                   You should have received a copy of the GNU Affero General Public License
%                   along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function retFlowSimFunction(heading_ang, depth, fix_depth, save_dir)

    debug = false; %generate debug plots?
    spotcheck = false; %pause on each debug plot?

    % debug = false; %generate debug plots?
    % spotcheck = false; %pause on each debug plot?

    headFlowBool = false;

    plane_depth = depth; % vary this from -4.5 to 6 in 0.5 increments
    % 2/24/2022: added addition cases where the plane is behind fixation
    basePath = 'C:\Users\Oli\Documents\charlie_data\kate_data\simulation';
    dataPath = 'C:\Users\Oli\Documents\charlie_data\kate_data\simulation';
    toolboxPath = 'C:\Users\Oli\Documents\charlie_data\kate_data\simulation';
    savePathPath = [save_dir int2str(plane_depth)];
    addpath(genpath(dataPath))
    addpath(genpath(basePath))
    addpath(genpath(toolboxPath))

    sessionID = 'eyeSim'; terrainID = 'FlatGroundplane'; trialType = 'GazeToLeftOfPath'; walkNum = nan; %simulated data of an eyeball moving in a straight line, looking at a point on its path

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Step #1 - Initialize Stuff!!
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        if ~isnan(walkNum)   %if there is a number in the variable 'walkNum,' then we will base the simulation on real laser skeleton data.

            skelBool = true; %if true, we are use skeleton eye data, otherwise make your own!

            %load AllWalks data
            allWalksPath = [dataPath filesep terrainID '_allWalks.mat' ]; %absolute path to .mat containing allWalks
            load(allWalksPath)

            w = allWalks{walkNum}; %'w' is the struct with all the data for WalkNum included <3

            %savePath = [savePathPath w.subID '_' w.takeID  '_' w.trialType '_' num2str(walkNum)];
            savePath = savePathPath;
        else %purely simulated eyeball/optic flow data

            skelBool = false;

            savePath = [savePathPath sessionID '_' terrainID '_' trialType  '_' trialType '_' num2str(walkNum) ];
            savePath = savePathPath;
        end

        if headFlowBool
            savePath = [savePath '_headFlow']; %mostly vestigial. I think this bool paralyzes the eyeball
        end

        mkdir(savePath)

        %% PREPARE SKELETON DATA & DEFINE EYE PATH
        if skelBool
            %%% clean fixation data  - pin gaze to ground during fixation
            grHeight = fillmissing(w.rGazeGroundIntersection(:,2),'nearest'); %get

            %%%flip gaze data so Z is up (to better interface with Geom3D)

            fixXYZraw(:,1) = w.rGazeGroundIntersection(:,1);
            fixXYZraw(:,2) = -w.rGazeGroundIntersection(:,3);
            fixXYZraw(:,3) = fillmissing(w.rGazeGroundIntersection(:,2),'nearest');

            fixXYZraw(:,3) = fixXYZraw(:,3)-grHeight; %set ground Height to zero
            comXYZ = yUp2zUp(w.comXYZ);

            %%%remove distant fixations (they make the simulation blow up, as they
            %%%include too many dottos inthe FOV
            [gTh, gR] = cart2pol(fixXYZraw(:,1) - comXYZ(:,1), fixXYZraw(:,2)-comXYZ(:,2));
            distCutOff = 1e5; 
            fixXYZraw(gR>distCutOff,:) = nan;

            fixXYZ = pinGazeToGround(fixXYZraw, w, debug);

            %Flip skeleton data so Z is up (and smooth it)
            order = 4;
            cutoff = 5;
            samplingRate = mean(diff(w.syncedUnixTime))^-1; %120

            shadow_fr_mar_dim = w.shadow_fr_mar_dim;
            for mm = 1:length(shadow_fr_mar_dim(1,:,1))
                shadow_fr_mar_dim(:,mm,:) = butterLowZero(order,cutoff,samplingRate, yUp2zUp(shadow_fr_mar_dim(:,mm,:)));
                shadow_fr_mar_dim(:,mm,3) =  shadow_fr_mar_dim(:,mm,3)  - grHeight;
            end

            %load steps
            steps_HS_TO_StanceLeg_XYZ = w.steps_HS_TO_StanceLeg_XYZ;

            rHeel = squeeze(shadow_fr_mar_dim(:,strcmp('RightHeel',w.shadowMarkerNames),:));
            lHeel = squeeze(shadow_fr_mar_dim(:,strcmp('LeftHeel',w.shadowMarkerNames),:));
            for ss = 1:length(steps_HS_TO_StanceLeg_XYZ)

                switch steps_HS_TO_StanceLeg_XYZ(ss, 3)

                    case 1 %the RIGHT FOOT is on the ground during this step (...I think I should've named this variable "foothold...")
                        steps_HS_TO_StanceLeg_XYZ(ss, [4:6]) = rHeel(steps_HS_TO_StanceLeg_XYZ(ss, 1),:);


                    case 2 %the LEFT FOOT is on the ground during this step
                        steps_HS_TO_StanceLeg_XYZ(ss, [4:6]) = lHeel(steps_HS_TO_StanceLeg_XYZ(ss, 1),:);
                end
            end


            %%% Define eye path (either from skeleton data, or created)

            rEyeballCenterXYZ  = yUp2zUp(w.rEyeballCenterXYZ);
            lEyeballCenterXYZ  = yUp2zUp(w.lEyeballCenterXYZ);
            clear eyeXYZ
            eyeXYZ(:,1) = mean([rEyeballCenterXYZ(:,1) lEyeballCenterXYZ(:,1)],2);
            eyeXYZ(:,2) = mean([rEyeballCenterXYZ(:,2) lEyeballCenterXYZ(:,2)],2);
            eyeXYZ(:,3) = mean([rEyeballCenterXYZ(:,3) lEyeballCenterXYZ(:,3)],2) - grHeight;

            %     Smoove the carp out of that eye trajectory
            samplingRate = mean(diff(w.syncedUnixTime))^-1; %120
            eyeXYZ = butterLowZero(order,cutoff,samplingRate, eyeXYZ);

        elseif ~skelBool %simData
            % params: heading is at 30 deg, subject starts at -10e3 in the
            % x-axis, fixation is at 1e3 in the x-axis; plane is allowed to
            % vary in the x-axis
            eyeHeight = 1.5e3;
            samplingRate = 120;
            heading_angle = heading_ang;
            heading_ratio = tand(heading_angle)*0.2;
            theta = atand(heading_ratio);% initial heading
            switch trialType
                case 'GazeOnPath'
                    eyeStartPosition = [-5e3 0 0];
                    eyeEndPosition   = [-4.5e3 0 0]; % changed here

                case 'GazeToRightOfPath'
                    eyeStartPosition = [-5e3 0 0];
                    eyeEndPosition   = [-4e3 1e3 0];

                case 'GazeToLeftOfPath'
                    eyeStartPosition = [-5e3 0 0];
                    % eyeStartPosition = [-4e3 -1e3*heading_ratio 0];
                    eyeEndPosition   = [-4.8e3 -1e3*heading_ratio 0];
                    % eyeEndPosition = [-5e3 0 0];
            end
            if sum(strcmp(trialType, {'SinX', 'SinY', 'SinZ', 'SinYCosZ'}))>0
                eyeStartPosition = [-5e3 0 0];
                eyeEndPosition   = [-4e3 0 0];
            end

            %%%calc numFrames based on enforcing a 1000mm/s eye velocity(@120fps)
            targetSpeed = 1500; %mm/s <-doesn't really work, I think ? Who knows how fast that dang ol eyeball gonna move. wtfever
            % igdr5t
            eyeDistTraveled = pdist([eyeStartPosition; eyeEndPosition]);
            numFrames =round((eyeDistTraveled/targetSpeed) * samplingRate);

            % We can change the fixation here % added by Mengjian Hua
            fixXYZ = zeros(numFrames, 3);
            x_place = 7.5e3; % significantly behind fixation relative to observer
            % fixXYZ(:,1) = x_place;
            x_fix = 7.5e3;
            fixXYZ(:, 1) = fix_depth; % fixation in x-axis
            fixXYZ(:,2) = 0;%-(x_place-eyeStartPosition(1))*1;
            fixXYZ(:,3) = 1.5e3;
            % fixXYZ(:, 3) = 0;
            eyeX = linspace(eyeStartPosition(1),eyeEndPosition(1),numFrames)';
            eyeY = linspace(eyeStartPosition(2),eyeEndPosition(2),numFrames)';
            eyeZ = zeros(numFrames, 1) + eyeHeight;
            sec = numFrames/samplingRate;

            eyeXYZ = [eyeX eyeY eyeZ];
            % these variables sit on a throne of LIES!!!
            shadow_fr_mar_dim = nan(numFrames, 30, 3);
            steps_HS_TO_StanceLeg_XYZ = [nan nan nan nan nan nan]; %batman!

        end

        if headFlowBool
            fixXYZ_orig = fixXYZ;
            %         meanAllFixPoint_Xz =  mean(fixXYZ(:,1) - eyeXYZ(:,1));
            %         meanAllFixPoint_Yz =  mean(fixXYZ(:,2) - eyeXYZ(:,2));
            %         meanAllFixPoint_Zz =  mean(fixXYZ(:,3) - eyeXYZ(:,3));
            %         fixXYZ = [meanAllFixPoint_Xz+eyeXYZ(:,1) meanAllFixPoint_Yz+eyeXYZ(:,2) meanAllFixPoint_Zz+eyeXYZ(:,3)];

            fixXYZ = w.headVecX_fr_xyz+eyeXYZ;
        end
         %% DEFINE GROUND PLANE (flat for now, but functionality for bumpy terrain exists in \ResearchProjects\FlowSimulationProject\3dFlowSimulation_old)

        origin = [0 0 0];
        groundPlane = createPlane([x_place 0 0], [x_place 1 1], [x_place 1 0]); % changed here 2/24 by Mengjian Hua (yz plane or x=0)
        % groundPlane = createPlane([1.75e5 1 0], [1.75e5 0 1], [1.75e5 0 0]); % create plane at edge of point cloud 7.7.21
        gridRes = 200;
        gridRange = max(fixXYZ(:,1)+1e4)*10; %big ol bongo grid of dotto's
        % planeSize = 8000*abs(-5e3-plane_depth)/0.5e3;
        planeSize = 88000*2; % test (original 8000/scaling^)
        % gapSize = 1000;
%         planeRes = planeSize/200; % keep resolution constant
        planeRes = 100; % test
        % one plane
        yVec = (-1/2*planeSize):planeRes:(1/2*planeSize);
        zVec = (-1/2*planeSize)+eyeHeight:planeRes:(1/2*planeSize)+eyeHeight;
        [yy, zz] = meshgrid(yVec, zVec);
        % xx controls the depth of the plane in the x-axis; we vary it from +1
        % to -9 in intervals of 0.5 (subject to change)
        % temp_xx = linspace(1, -9, 21);
        xx = 1*plane_depth*ones(size(yy));
        distr = [xx(:), yy(:), zz(:), [1:length(xx(:))]'];

        groundPoints_xyz_id = distr;
        dotColors = lines(length(groundPoints_xyz_id(:,1)));


        %% Create eyeStruct - with all the basic information that will be shared across all frames

        baseEyeStruct = [];
        baseEyeStruct.eyeXYZ = eyeXYZ;
        baseEyeStruct.fixXYZ = fixXYZ;

        if skelBool
            baseEyeStruct.eyeRadius = 12; %anatomically realistic eyeball (12mm radius)
        else
            baseEyeStruct.eyeRadius = 200; %12; %Honkin huge eyeball for... visibility, lol
        end

        baseEyeStruct.samplingRate = samplingRate;
        baseEyeStruct.trialType = trialType;
        baseEyeStruct.terrainID = terrainID;
        baseEyeStruct.sessionID =sessionID;
        baseEyeStruct.origin = origin;
        baseEyeStruct.groundPlane = groundPlane;
        baseEyeStruct.groundPoints_xyz_id = groundPoints_xyz_id;
        baseEyeStruct.xx = xx;
        baseEyeStruct.yy = yy;
        baseEyeStruct.zz = zz;
        baseEyeStruct.dotColors = dotColors;
        baseEyeStruct.shadow_fr_mar_dim = shadow_fr_mar_dim;
        baseEyeStruct.steps_HS_TO_StanceLeg_XYZ = steps_HS_TO_StanceLeg_XYZ;
        baseEyeStruct.savePath = savePath;

        if skelBool
            baseEyeStruct.walkNum = w.ww;
        else
            baseEyeStruct.walkNum = 'nah';
        end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Step 2 - (LOOP #1)  - For each frame of the animation, project ground dottos onto the back of the simulated eyeball1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if skelBool
            startFrame = round(length(w.frames)/4);
            endFrame = startFrame+150;%  round(length(w.frames)/4)*3;

            if strcmp(savePath,     'D:\OpticFlowSimFiles\JSM_Woodchips_Ground_2')
                startFrame = 16192;
                endFrame = startFrame+1500;
            end
            numFrames = endFrame-startFrame;

        else
            startFrame = 1;
            endFrame = startFrame+numFrames-1;
        end

        baseEyeStruct.startFrame  = startFrame;
        baseEyeStruct.endFrame  = endFrame;
        baseEyeStruct.numFrames = numFrames;

        % Then construct a ParforProgressbar object:
    %     parforprogressbar1 = ParforProgressbar(numFrames, 'showWorkerProgress', true, 'title', 'First Parfor Progress') ;
         disp('Starting first for loop')

        tic
        %     parfor fr  = startFrame:endFrame %parallel processing version (debug     and spotcheck bools don't work in this mode
        for fr  = startFrame:endFrame %single threaded version
            if mod(fr,10)==0; disp(['Loop #1 - ', num2str(fr), ' of ', num2str(endFrame)]); end %turn this off when parallel processing
            
            %% new stuff (experimental?)
            if fr > startFrame % after fr loops past the first frame, store the previous frame's data as a separate var
                prevEyeStruct = thisEyeStruct;
            end
            newEyeStruct = [];
            thisEyeStruct = baseEyeStruct; %create a new eyeStruct for this frame
            thisEyeStruct.baseEyeStruct = thisEyeStruct;


            thisEyeStruct.fr = fr;
            %% define eyeball (including pupil, fovea, and  defining sagital, coronal, and transverse planes) and point it at the fixation point

            %         [eyeMesh, eyeSphere, thisFixPoint, gazeLine, pupilXYZ, foveaXYZ, rotXToGzEul, rotXToGz, eyeXvec, eyeYvec, eyeZvec, coronalEyePlane, saggitalEyePlane, transverseEyePlane , coronalEyeCircle, saggitalEyeCircle,transverseEyeCircle, eyeCardinalPointNorth, eyeCardinalPointSouth, eyeCardinalPointEast, eyeCardinalPointWest, saggitalEyeCirclePoints, coronalEyeCirclePoints, transverseEyeCirclePoints] ...
            %             = defineEyeball(thisEyeStruct, fr,  debug, spotcheck);
            [eyeMesh,...
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
                defineEyeball(thisEyeStruct,  fr, debug, spotcheck);

            thisEyeStruct.eyeMesh = eyeMesh;
            thisEyeStruct.eyeSphere = eyeSphere;
            thisEyeStruct.eyeFixPoint= thisFixPoint;
            thisEyeStruct.gazeLine = gazeLine;
            thisEyeStruct.pupilXYZ = pupilXYZ;
            thisEyeStruct.foveaXYZ = foveaXYZ;
            thisEyeStruct.rotXToGzEul = rotXToGzEul;
            thisEyeStruct.rotXToGz = rotXToGz;
            thisEyeStruct.eyeXvec = eyeXvec;
            thisEyeStruct.eyeYvec = eyeYvec;
            thisEyeStruct.eyeZvec = eyeZvec;
            thisEyeStruct.coronalEyePlane = coronalEyePlane;
            thisEyeStruct.saggitalEyePlane = saggitalEyePlane;
            thisEyeStruct.transverseEyePlane = transverseEyePlane;
            thisEyeStruct.coronalEyeCircle = coronalEyeCircle;
            thisEyeStruct.saggitalEyeCircle = saggitalEyeCircle;
            thisEyeStruct.transverseEyeCircle = transverseEyeCircle;

            thisEyeStruct.saggitalEyeCirclePoints = saggitalEyeCirclePoints;
            thisEyeStruct.coronalEyeCirclePoints = coronalEyeCirclePoints;
            thisEyeStruct.transverseEyeCirclePoints = transverseEyeCirclePoints;

            thisEyeStruct.verticalEyeAxisGroundIntersect =verticalEyeAxisGroundIntersect;
            thisEyeStruct.horizontalEyeAxisGroundIntersect =horizontalEyeAxisGroundIntersect;

            %eye cardinal points
            thisEyeStruct.eyeCardinalPointNorth = eyeCardinalPointNorth;
            thisEyeStruct.eyeCardinalPointSouth = eyeCardinalPointSouth;
            thisEyeStruct.eyeCardinalPointEast = eyeCardinalPointEast;
            thisEyeStruct.eyeCardinalPointWest = eyeCardinalPointWest;







            %% define eye field of view (fov) by projecting a cone out from the pupil and intersecting it with the groundplane
            distCutOff = 1e5; %Changed by Mengjian Hua Mar 26
            fovRadiusDeg = 60; %field of view radius in degrees
            [fovRadiusDeg,fovRadiusRad, fovGroundPoints] = defineFOV(fovRadiusDeg, distCutOff, thisEyeStruct, fr,debug, spotcheck);

            thisEyeStruct.fovRadDeg = fovRadiusDeg;
            thisEyeStruct.fovRadRad = fovRadiusRad;
            thisEyeStruct.fovGroundPoints = fovGroundPoints;
            % thisEyeStruct.horizontalEyeCrossHairGroundPoints = horizontalEyeCrossHairGroundPoints;
            % thisEyeStruct.verticalEyeCrossHairGroundPoints = verticalEyeCrossHairGroundPoints;
            %% make points for Fov-ground grid
            distCutOff = 1e5;

            fovGridRange = linspace(0,fovRadiusDeg,5);
            fovGridRange(1) = 2;%foveal ring
            for  ff = 1:length(fovGridRange)
                [thisFovRadiusDeg,thisFovRadiusRad, thisFovGroundPoints] = defineFOV(fovGridRange(ff), distCutOff, thisEyeStruct, fr,debug, spotcheck);

                thisEyeStruct.fovGridRadDeg{ff} = thisFovRadiusDeg;
                thisEyeStruct.fovGridRadRad{ff} = thisFovRadiusRad;
                thisEyeStruct.fovGridGroundPoints{ff} = thisFovGroundPoints;
            end
            %% define heading and T, dT/dt
            if(fr ==endFrame)
                heading_dir = (eyeXYZ(fr,:) - eyeXYZ(fr-1,:));
                heading_line = createLine3d(eyeXYZ(fr-1,1),eyeXYZ(fr-1,2),eyeXYZ(fr-1,3),heading_dir(1),heading_dir(2),heading_dir(3));
                heading = intersectLinePlane(heading_line, groundPlane);
            else
                heading_dir = (eyeXYZ(fr+1,:) - eyeXYZ(fr,:));
                heading_line = createLine3d(eyeXYZ(fr,1),eyeXYZ(fr,2),eyeXYZ(fr,3),heading_dir(1),heading_dir(2),heading_dir(3));
                heading = intersectLinePlane(heading_line, groundPlane);
            end
            T = heading_dir;
            thisEyeStruct.translation = T;
            %% Project ground dots onto pupil, and remove ground occlusions (projected points that intersect with the ground plane <--This part hella slow)
            switch trialType
                    case 'SinX'
                        translation_x = DsinPath(fr);
                        translation_y = DeyeY(fr);
                        translation_z = DsinPath(fr);
                    case 'SinY'
                        translation_y = DcosPath(fr);
                        translation_x = DeyeX(fr);
                        translation_z = DsinPath(fr);
            end

            thisEyeStruct.heading = heading; %+ [thisFovGroundPoints(1,1) thisFovGroundPoints(1,2) 0]; % added by Mengjian Hua
            [backEyeProjPoints_xyz_id,back_proj_heading] = projectDotsToPupil(thisEyeStruct, fr, debug, spotcheck);
            thisEyeStruct.backEyeProjPoints_xyz_id= backEyeProjPoints_xyz_id;


            %% resituate backProjPoints into polar coordinates centered on "fovea," using "great circle distance" for eccentricity (rather than just projecting to a plane)

    % 5/23       [backEyeProjPoints_xyz_id_rot, retinaProjPoints_fr_theta_rho_id, retinaProjPoints_fr_X_Y_id,headX,headY] = defineRetinalPolarCoordinates(thisEyeStruct,back_proj_heading, fr, debug, spotcheck);
            [~,retinaProjPoints_fr_X_Y_id,headX,headY] = defineRetinalPolarCoordinates(thisEyeStruct,back_proj_heading, fr, debug, spotcheck);
            %%%% points projected on sphere, rotated so fovea and pupil are on the X axis
    %         thisEyeStruct.backEyeProjPoints_xyz_id_rot = backEyeProjPoints_xyz_id_rot;
    %         thisEyeStruct.retinaProjPoints_fr_theta_rho_id = retinaProjPoints_fr_theta_rho_id;
            thisEyeStruct.retinaProjPoints_fr_X_Y_id = retinaProjPoints_fr_X_Y_id;
            % store fixation
            newEyeStruct.headX = headX;
            newEyeStruct.headY = headY;
            
            %% more new stuff
            if fr > startFrame
                [retinaProjVel_fr_xVel_yVel_id] = calcFlowPerFrame(thisEyeStruct, prevEyeStruct, fr, debug, spotcheck); %gets vel in X Y units...
            else
                [retinaProjVel_fr_xVel_yVel_id] = calcFlowPerFrame(thisEyeStruct, thisEyeStruct, fr, debug, spotcheck); %gets vel in X Y units...
            end
            
            thisEyeStruct.retinaProjVel_fr_xVel_yVel_id = retinaProjVel_fr_xVel_yVel_id;

            %% Calculate retinal vector field
            retGridRes = 200; %larger numbers = more grid points = higher resolution
            if mod(retGridRes,2)==0; retGridRes = retGridRes+1;end %make sure the res# odd to make sure the middle value passes through zero

%             gridRadius =  (2*pi*thisEyeStruct.eyeRadius)/4;
            % if things are in degrees
            gridRadius = fovRadiusDeg;
            newEyeStruct.samplingRate = samplingRate;
            newEyeStruct.gridRadius = gridRadius;

            [retGridxx, retGridyy, retGridVelxx, retGridVelyy, retGridIDs] = calcRetinalVectorField(thisEyeStruct, retGridRes, gridRadius, fr, debug, spotcheck);

            newEyeStruct.retGridxx = retGridxx;
            newEyeStruct.retGridyy = retGridyy;

            newEyeStruct.retGridVelxx = retGridVelxx;
            newEyeStruct.retGridVelyy = retGridVelyy;

            newEyeStruct.retGridIDs = retGridIDs;
            
            %% Save this eyestruct out to a file, which will be loaded later to determine flow per frame
            savePathFile = [savePath filesep 'eyeStructFrame' sprintf('%09d',fr) '.mat'];
            parsave(savePathFile, newEyeStruct)
        end

        t = toc;
        disp([ num2str(t) ' seconds elapsed. ' num2str(t/(endFrame-startFrame)) ' seconds per frame'])    

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% internal fucntion to allow saving of files within a parfor loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function parsave(fname, thisEyeStruct)
        save(fname, 'thisEyeStruct', '-v7.3')
    end
end