function [eyeStruct] = defineEyePath(eyeStruct)

varNames = fieldnames(eyeStruct);
for i=1:length(varNames)
    eval([varNames{i} '=eyeStruct.' varNames{i} ';']);
end

varNames = fieldnames(simDeets);
for i=1:length(varNames)
    eval([varNames{i} '=simDeets.' varNames{i} ';']);
end

if skelBool
    eyeStruct.headXYZ = eyeStruct.eyeXYZ;
    eyeStruct.fixatedPointXYZ = eyeStruct.fixXYZ;
    eyeStruct.headRotMat  = [];
    eyeStruct.fixPoints_fr_XYZ =  eyeStruct.fixXYZ;
    
    
    if headFlowBool
        eyeStruct.head.eyeXYZ = eyeStruct.headXYZ;
    else
        if ~twoEyeBool
            eyeStruct.eye.eyeXYZ = eyeStruct.headXYZ ;
        end
    end
else
    
    eyePathType = simDeets.eyePathType;
    
                headX = linspace(-2000,0,numFrames);
            %         headX = [headX fliplr(headX)];% headX fliplr(headX)];
            
            headY = linspace(0,0,numFrames);
            
            headZ =  linspace(1.5e3,1.5e3,numFrames);
    switch eyePathType
        
        case 'straight'
            headX = linspace(-2000,0,numFrames);
            %         headX = [headX fliplr(headX)];% headX fliplr(headX)];
            
            headY = linspace(0,0,numFrames);
            
            headZ =  linspace(1.5e3,1.5e3,numFrames);
            
        case 'diagonal'
            headX = linspace(-2000,fixXYZ(1,1),numFrames);
            
            headY = linspace(0,1000,numFrames);
            
            headZ =  linspace(3e3,3e3,numFrames);
            
        case 'horizSin'
            headX = linspace(-4000,fixatedPointXYZ(1)-1000,numFrames);
            
            headY = sin(linspace(0,2*pi,numFrames))*simDeets.sigma;
            
            headZ =  linspace(8e3,8e3,numFrames);
            
        case'vertSin'
            headX = linspace(-4000,fixatedPointXYZ(1)-1000,numFrames);
            
            headY = linspace(0,2*pi,numFrames);
            
            headZ =   sin(linspace(0,2*pi,numFrames))*simDeets.sigma +8e3;
            
            
        case'corkscrew'
            headX = linspace(-4000,fixatedPointXYZ(1)-1000,numFrames);
            
            headY = sin(linspace(0,2*pi,numFrames))*simDeets.sigma;
            
            headZ = cos(linspace(0,2*pi,numFrames))*simDeets.sigma +8e3;
            
    end
    
    eyeStruct.headXYZ = [headX' headY' headZ'];
    
    
    
    eyeStruct.fixXYZ(:,3) = zLook(fixXYZ(:,1), fixXYZ(:,2));
    
    
    if headFlowBool
        eyeStruct.head.eyeXYZ = eyeStruct.headXYZ;
    else
        if ~twoEyeBool
            eyeStruct.eye.eyeXYZ = eyeStruct.headXYZ ;
        end
    end
end