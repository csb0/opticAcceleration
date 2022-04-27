function [eyeStruct] = defineGround(eyeStruct)

varNames = fieldnames(eyeStruct);
for i=1:length(varNames)
    eval([varNames{i} '=eyeStruct.' varNames{i} ';']);
end

varNames = fieldnames(simDeets);
for i=1:length(varNames)
    eval([varNames{i} '=simDeets.' varNames{i} ';']);
end


eyeStruct.origin = [0 0 0];
eyeStruct.groundPlane = createPlane(eyeStruct.origin, [1 0 0], [0 1 0]);

% if skel
%     eyeStruct.xVec = [gridRange(1,1):gridRes:gridRange(2,1)];
%     eyeStruct.yVec = [gridRange(1,2):gridRes:gridRange(2,2)];
% else
    eyeStruct.xVec = [-gridRange:gridRes:gridRange];
    eyeStruct.yVec = [-gridRange:gridRes:gridRange];
% end

[eyeStruct.xx, eyeStruct.yy] = meshgrid(eyeStruct.xVec, eyeStruct.yVec);

% zz = zeros(size(xx));
% zVec = sin(linspace(0,20*pi,numel(xVec)))*5e2;
% [xz, yz] = meshgrid(linspace(0,50*pi,numel(xVec)),linspace(0,50*pi,numel(yVec)));

xc = 0;
yc = 0;


[eyeStruct.zz] = gaussianPlot(eyeStruct.xx, eyeStruct.yy, xc, yc, simDeets.sigma);
eyeStruct.zz = eyeStruct.zz * terrainHeight;


eyeStruct.groundPoints_xyz_id = [eyeStruct.xx(:), eyeStruct.yy(:), eyeStruct.zz(:), [1:length(eyeStruct.xx(:))]' ];
eyeStruct.groundPoints_xyz_id = [eyeStruct.xx(:)+rand(size(eyeStruct.xx(:)))*jitterVal, eyeStruct.yy(:)+rand(size(eyeStruct.xx(:)))*jitterVal, eyeStruct.zz(:) [1:length(eyeStruct.xx(:))]' ]; %jitter 'em dots

eyeStruct.dotColors = lines(length(eyeStruct.groundPoints_xyz_id(:,1)));
