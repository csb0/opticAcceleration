function [eyeStruct] = calcFlowPerFrameDVA(eyeStruct, thisEyeD, fr, debug)

varNames = fieldnames(eyeStruct);
for i=1:length(varNames)
    eval([varNames{i} '=eyeStruct.' varNames{i} ';']);
end

thisEyeStruct = eyeStruct.(thisEyeD);
varNames = fieldnames(thisEyeStruct);
for i=1:length(varNames)
    eval([varNames{i} '=thisEyeStruct.' varNames{i} ';']);
end


retinaProjVel_fr_azVel_elVel_id{fr} = nan(length(backEyeProjPoints_xyz_id_z{fr}),3);
retinaProjVel_fr_azVel_elVel_id{fr}(:,3) = backEyeProjPoints_xyz_id_z{fr}(:,4);%keep dot IDs

pupil = [eyeRadius, 0, 0];

if  fr > 1 && ~isempty(retinaProjVel_fr_azVel_elVel_id{fr-1}) % don't go through this loop the first time around
    prevFramePoints = backEyeProjPoints_xyz_id_z{fr-1};
    thisFramePoints = backEyeProjPoints_xyz_id_z{fr};
    
    for pp = 1:length(thisFramePoints)
        
        thisPointCurr= thisFramePoints(pp,:);
        
        thisPointPrev = prevFramePoints(prevFramePoints(:,4) == thisPointCurr(4), :);%find this dotto's location on the previous frame
        
        if isempty(thisPointPrev)
            retinaProjVel_fr_azVel_elVel_id{fr}(pp,:) = [nan nan thisPointCurr(4)]; %if this guy wasn't around last frame, his speed is listed as "nan"
        else
            
            assert(length(thisPointPrev) == 4)
            assert(length(thisPointCurr) == 4)
            assert(thisPointCurr(4) == thisPointPrev(4))
            
            %center on pupil
            thisPointCurr_z = thisPointCurr(1:3) - pupil;
            thisPointPrev_z = thisPointPrev(1:3) - pupil;
            
            
            %azimuth (XY) and Elevation (XZ) angle
            [thisPointPrevAz, ~] = cart2pol(thisPointPrev_z(1), thisPointPrev_z(2));
            [thisPointCurrAz, ~] = cart2pol(thisPointCurr_z(1), thisPointCurr_z(2));
            AzAngle = rad2deg(thisPointPrevAz - thisPointCurrAz);
            
            [thisPointPrevEl, ~] = cart2pol(thisPointPrev_z(1), thisPointPrev_z(3));
            [thisPointCurrEl, ~] = cart2pol(thisPointCurr_z(1), thisPointCurr_z(3));
            ElAngle = rad2deg(thisPointPrevEl - thisPointCurrEl);
            
            if abs(AzAngle) > 150
                if pp >1
                AzAngle = retinaProjVel_fr_azVel_elVel_id{fr}(pp-1,1)/framerate;%this is sloppy, but I don't feel like figuring out this wrapparound crap right now
                else
                    AzAngle = nan;
                end
            end
            
            if abs(ElAngle) > 150
                if pp > 1
                ElAngle = retinaProjVel_fr_azVel_elVel_id{fr}(pp-1,2)/framerate;%this is sloppy, but I don't feel like figuring out this wrapparound crap right now
                else
                    ElAngle = nan;
                end
            end
            
            retinaProjVel_fr_azVel_elVel_id{fr}(pp,:) = [AzAngle*framerate, ElAngle*framerate, thisPointCurr(4)]; %save as Degrees Per Second!
            
%                         %%%%%debug (shows point-by-point comparison)
%                         figure(547237);clf;
%                         subplot(1,3,1)
%                         hold on;view(3);axis equal
%                         %%%%% plot eye ball
%                         plot3(-eyeRadius,0,0, 'ko','MarkerFaceColor','g')% fovea
%                         plot3(eyeRadius,0,0, 'ko','MarkerFaceColor','r') % pupil
%                         s1 = drawSphere([0 0 0 eyeRadius]);
%                         s1.FaceColor = [1 1 1]*.75;
%             
%                         s1.FaceAlpha = .5;
%                         s1.LineStyle = '-';
%             
%                         %%% plot XYZ (RGB) axes
%                         plot3([0 eyeRadius*1.5], [0 0], [ 0 0], 'r:','LineWidth',3)
%                         plot3([0 0], [0 eyeRadius*1.5], [ 0 0], 'g:','LineWidth',3)
%                         plot3([0 0], [0 0], [ 0 eyeRadius*1.5], 'b:','LineWidth',3)
%                         plot3(-[0 eyeRadius*1.5], [0 0], [ 0 0], 'r:','LineWidth',3)
%                         plot3([0 0], -[0 eyeRadius*1.5], [ 0 0], 'g:','LineWidth',3)
%                         plot3([0 0], [0 0], -[ 0 eyeRadius*1.5], 'b:','LineWidth',3)
%                         plot3(-eyeRadius,0,0, 'ko','MarkerFaceColor','g')% fovea
%                         plot3(pupil(1), pupil(2), pupil(3), 'ko','MarkerFaceColor','r','MarkerSize',24)% pupil
%             
%                         %%%plots the dots
%                         plot3(backEyeProjPoints_xyz_id_z{fr}(:,1), backEyeProjPoints_xyz_id_z{fr}(:,2), backEyeProjPoints_xyz_id_z{fr}(:,3),'m.')
%                         plot3(backEyeProjPoints_xyz_id_z{fr-1}(:,1), backEyeProjPoints_xyz_id_z{fr-1}(:,2), backEyeProjPoints_xyz_id_z{fr-1}(:,3),'b.')
%                         plot3([eyeRadius thisPointPrev(:,1)], [0 thisPointPrev(:,2)], [0 thisPointPrev(:,3)],'r-p','LineWidth',2)
%                         plot3([eyeRadius thisPointCurr(:,1)], [0 thisPointCurr(:,2)], [0 thisPointCurr(:,3)],'c-p','LineWidth',2)
%             
%                         subplot(1,3,2)
%                         hold on; axis equal
%                         viscircles([0 0]-[eyeRadius 0], eyeRadius,'Color','k');
%                         plot([0 thisPointPrev_z(:,1)], [0 thisPointPrev_z(:,3)],'r-p','LineWidth',2)
%                         plot([0 thisPointCurr_z(:,1)], [0 thisPointCurr_z(:,3)],'c-p','LineWidth',2)
%                         title(['X-Z view - Elevation ' num2str(ElAngle)])
%             
%                         subplot(1,3,3)
%                         hold on; axis equal
%                         viscircles([0 0]-[eyeRadius 0], eyeRadius,'Color','k');
%                         plot([0 thisPointPrev_z(:,1)], [0 thisPointPrev_z(:,2)],'r-p','LineWidth',2)
%                         plot([0 thisPointCurr_z(:,1)], [0 thisPointCurr_z(:,2)],'c-p','LineWidth',2)
%             
%                         title(['X-Y view - Azimuth ' num2str(AzAngle)])
%             
%                         drawnow
%             
            
            
        end
    end
    
    assert(sum(backEyeProjPoints_xyz_id_z{fr}(:,4) - retinaProjVel_fr_azVel_elVel_id{fr}(:,3))==0); %make sure ID's are maintained in both position and velocity data
end

eyeStruct.(thisEyeD).retinaProjVel_fr_azVel_elVel_id{fr} = retinaProjVel_fr_azVel_elVel_id{fr};

if  fr > 1 && ~isempty(retinaProjVel_fr_azVel_elVel_id{fr-1}) % don't go through this loop the first time around
    if debug
        %%
        figure(32987);clf
        
        dotX = backEyeProjPoints_xyz_id_z{fr}(:,2);
        dotY = backEyeProjPoints_xyz_id_z{fr}(:,3);
        dotVelX  = retinaProjVel_fr_azVel_elVel_id{fr}(:,1);
        dotVelY  = retinaProjVel_fr_azVel_elVel_id{fr}(:,2);
        
        prevDotX = backEyeProjPoints_xyz_id_z{fr-1}(:,2);
        prevDotY = backEyeProjPoints_xyz_id_z{fr-1}(:,3);
        
        plot(dotX, dotY,'k.')
        hold on
        plot(prevDotX, prevDotY,'r.')
        
        q1 = quiver(dotX, dotY, dotVelX, dotVelY);
        
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
end
