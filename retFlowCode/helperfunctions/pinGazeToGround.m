function [fixPinnedXYZ] = pinGazeToGround(fixXYZraw,w, debug)

gFx = fixXYZraw(:,[1 2]);

%% stupid/simple velocity based saccade detector (saccades == Eye-In-Head velocity > 65 deg/sec)
pxPerDeg = 18.33;

height = 1920;
width = 1080;

porX = w.gaze_norm_pos_x * width;
porY = height - w.gaze_norm_pos_y * height; %gotta do this weirdness to porY b/c of image vs XY coordinates

porY(isnan(porY)) = 0;
porX(isnan(porX)) = 0;

porXdisp = porX - nanmean(porX);
porYdisp = ((-porY)-nanmean(-porY)); %do some flippidoo nonsense to porY to make the trace look right (e.g. "downward saccades" correspond to the pupil moving downward)
porXdisp = porXdisp./pxPerDeg;
porYdisp = porYdisp./pxPerDeg;

gVraw = zeros (length(gFx),1);
for fr = 2:length(gFx)
    
    gVraw(fr) = sqrt( (porXdisp(fr-1) - porXdisp(fr)).^2 + (porYdisp(fr-1) - porYdisp(fr)).^2);
end

gV = smooth(gVraw);
framerate = 1/mean(diff(w.syncedUnixTime)); %120
gV = gV * framerate;%convert from deg/frame to deg/sec

saccThresh = 65; %thresh chosen by histogram, fairly arbitrary

saccFrames = gV>saccThresh;
sacc = gV;
sacc(~saccFrames) = nan;

%% go through fixXYZraw and set fixation posiiton to the median gaze positoin during that fixatoon 
fixPinnedXYZ = nan(size(fixXYZraw));

thisFixXYZ_fr = [];

for fr = 2:length(fixPinnedXYZ) %I really hate writing if-net code like this, and yet... here we are    
     if saccFrames(fr-1) == 0 && saccFrames(fr) == 0 %in a fixation, not the first frame
        inSacc = false;
        inFix = true;
        firstInSacc = false;
        firstInFix = false;
        
     elseif (fr==2) || (saccFrames(fr-1) == 1 && saccFrames(fr) == 0) %the first frame of a FIXATION %also do this for the first frame
        inSacc = false;
        inFix = true;
        firstInSacc = false;
        firstInFix = true;
        thisFixXYZ_fr = [];
    
    elseif saccFrames(fr-1) == 1 && saccFrames(fr) == 1 %in a SACCADE, not the first frame
        inSacc = true;
        inFix = false;
        firstInSacc = false;
        firstInFix = false;
        
    elseif saccFrames(fr-1) == 0 && saccFrames(fr) == 1 %the first frame of a SACCAD
        inSacc = true;
        inFix = false;
        firstInSacc = true;
        firstInFix = false;
    end
        

    if inFix
        thisFixXYZ_fr(end+1,1:4) = [fixXYZraw(fr,:) fr];
    end
    
    if firstInSacc
        thisFixFrames = thisFixXYZ_fr(:,4);
        
        thisFixMedianX = nanmedian(thisFixXYZ_fr(:,1));
        thisFixMedianY = nanmedian(thisFixXYZ_fr(:,2));
        thisFixMedianZ = nanmedian(thisFixXYZ_fr(:,3));
        
        fixPinnedXYZ(thisFixFrames,1) = thisFixMedianX;
        fixPinnedXYZ(thisFixFrames,2) = thisFixMedianY;
        fixPinnedXYZ(thisFixFrames,3) = thisFixMedianZ;
        
    end
    
    

end

fixPinnedXYZgaps = fixPinnedXYZ;
fixPinnedXYZ = fillmissing(fixPinnedXYZ,'nearest'); %replace missing nan values with nearest non-nan value


%%

if debug
    figure(732)
    clf
    subplot(2,1,1)
    hold on
    plot(gV,'.-')
    refline(0,saccThresh)
    plot(sacc,'r.-')
    
    subplot(2,1,2)
    plot(fixPinnedXYZ,'.-')
    hold on
    plot(fixXYZraw)
end
