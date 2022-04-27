function [simDeets, eyeStruct] = loadSkelEyeData(eyeStruct, basePath, debug)


varNames = fieldnames(eyeStruct);
for i=1:length(varNames)
    eval([varNames{i} '=eyeStruct.' varNames{i} ';']);
end


cd([basePath, filesep, 'Data' ])

load('allWalks.mat')


switch pathType
    case 'ground'
        w = woodchips{5}; %ground
    case 'fix'
        w = woodchips{3}; %fix
    case 'free'
        w = woodchips{1}; %free
    case 'rocks'
        w = rocks{1}; %rocks
end


if isnan(numFrames)
    startFrame = 1;
    endFrame = length(w.comXYZ);
    numFrames = length(startFrame:endFrame);
else    
    startFrame = 1000;
    endFrame = startFrame+numFrames;    
end
    eyeStruct.startFrame = startFrame;
    eyeStruct.endFrame = endFrame;
    eyeStruct.numFrames = numFrames;


%% clean fixation data  - pin gaze to ground during fixation

grHeight = fillmissing(w.rGazeGroundIntersection(:,2),'nearest');


if strcmp(w.trialType,'Fix')
    len = length(w.rGazeGroundIntersection(:,3));
    
    fixXYZraw(:,1) =  ones(len,1)*w.rEyeballCenterXYZ(end,1);
    fixXYZraw(:,2) = -ones(len,1)*w.rEyeballCenterXYZ(end,3);
    fixXYZraw(:,3) =  ones(len,1)*w.rEyeballCenterXYZ(end,2);
else
    fixXYZraw(:,1) = w.rGazeGroundIntersection(:,1);
    fixXYZraw(:,2) = -w.rGazeGroundIntersection(:,3);
    fixXYZraw(:,3) = fillmissing(w.rGazeGroundIntersection(:,2),'nearest');
end

fixXYZraw(:,3) = fixXYZraw(:,3)-grHeight; %set ground Height to zero


comXYZ = yUp2zUp(w.comXYZ);

%%%remove distant fixations (they make the simulation blow up, as they
%%%include too many dottos inthe FOV
[gTh, gR] = cart2pol(fixXYZraw(:,1) - comXYZ(:,1), fixXYZraw(:,2)-comXYZ(:,2));
%     figure; subplot(1,2,1);polarscatter(gTh, gR);  subplot(1,2,2); histogram(gTh);

distCutOff = 5000;

fixXYZraw(gR>distCutOff,:) = nan;

fixXYZ = pinGazeToGround(fixXYZraw, w, debug);




shadow_fr_mar_dim = w.shadow_fr_mar_dim;
for mm = 1:length(shadow_fr_mar_dim(1,:,1))
    shadow_fr_mar_dim(:,mm,:) = yUp2zUp(shadow_fr_mar_dim(:,mm,:));
    shadow_fr_mar_dim(:,mm,3) =  shadow_fr_mar_dim(:,mm,3)  - grHeight;
end


%% load simulation details




simDeets.gridRange = max(fixXYZ(startFrame:endFrame,1)+1e4); %big ol grid

simDeets.gridRes = 150;
simDeets.jitterVal = 0;
simDeets.groundType = w.condID; %woodchips or rocks
simDeets.gazeType = w.trialType; %free, ground, fix, rocks
simDeets.pathType = 'Skel';
simDeets.sigma = 1000;
simDeets.terrainHeight = 0;
simDeets.fixatedPoint = fixXYZ;



eyeStruct.steps_HS_TO_StanceLeg_XYZ =     w.steps_HS_TO_StanceLeg_XYZ;
eyeStruct.shadow_fr_mar_dim = shadow_fr_mar_dim;

eyeStruct.simDeets = simDeets;
eyeStruct.skel = true;
eyeStruct.fixXYZ = fixXYZ;
eyeStruct.grHeight = grHeight;


eyeStruct.eyeXYZ(:,1) = w.rEyeballCenterXYZ(:,1);
eyeStruct.eyeXYZ(:,2) = -w.rEyeballCenterXYZ(:,3);
eyeStruct.eyeXYZ(:,3) = w.rEyeballCenterXYZ(:,2)-grHeight;


order = 4;
cutoff = 5;
samplingRate = mean(diff(w.syncedUnixTime))^-1; %120
eyeStruct.eyeXYZ = butterLowZero(order,cutoff,samplingRate,eyeStruct.eyeXYZ);

eyeStruct.skelData = w;