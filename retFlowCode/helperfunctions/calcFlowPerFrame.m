function [retinaProjVel_fr_xVel_yVel_id] = calcFlowPerFrame(thisEyeStruct, prevEyeStruct, fr, debug, spotcheck);


thisFramePoints = thisEyeStruct.retinaProjPoints_fr_X_Y_id;
prevFramePoints = prevEyeStruct.retinaProjPoints_fr_X_Y_id;

retinaProjVel_fr_xVel_yVel_id = nan(size(thisFramePoints));


    for pp = 1:length(thisFramePoints)
        
        thisPointCurr= thisFramePoints(pp,:);
        
        thisPointPrev = prevFramePoints(prevFramePoints(:,3) == thisPointCurr(3), :);%find this dotto's location on the previous frame
        
        if isempty(thisPointPrev)
            retinaProjVel_fr_xVel_yVel_id(pp,:) = [nan nan thisPointCurr(3)]; %if this guy wasn't around last frame, his speed is listed as "nan"
        else
            assert(length(thisPointPrev) == 3)
            assert(length(thisPointCurr) == 3)
            assert(thisPointCurr(3) == thisPointPrev(3))
            


    
    retinaProjVel_fr_xVel_yVel_id(pp,:) = ...
        [thisPointCurr(1) - thisPointPrev(1)... %x vel
        thisPointCurr(2) - thisPointPrev(2)... %y vel
        thisPointCurr(3)];

        end
    end
    
    assert(sum(thisEyeStruct.retinaProjPoints_fr_X_Y_id(:,3) - retinaProjVel_fr_xVel_yVel_id(:,3))==0); %make sure ID's are maintained in both position and velocity data
%     if(fr == 10)
%     keyboard
%     end

if debug
    %%
    figure(32987);clf
    
    dotX = thisFramePoints(:,1);
    dotY = thisFramePoints(:,2);
    
    dotVelX = retinaProjVel_fr_xVel_yVel_id(:,1);
    dotVelY = retinaProjVel_fr_xVel_yVel_id(:,2);
    
    scale = 5;
    
    plot(dotX, dotY,'k.')
    hold on
    quiver(dotX, dotY, dotVelX*scale, dotVelY*scale)
    
    axis equal
    
    m = max([dotX; dotY]);
    plot([-m m], [0 0],'r-');
    plot( [0 0], [-m m],'r-');
    viscircles([0,0],m);
    
    if spotcheck
        dbstack
        keyboard
    end
end

