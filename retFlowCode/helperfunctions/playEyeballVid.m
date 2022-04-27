%
close all
% clear all

% savePath = 'D:\OpticFlowSimFiles\JSM_Rocks_Rocks_1';
% savePath = 'D:\OpticFlowSimFiles\JSM_Woodchips_Ground_2';
eyeStructStore = fileDatastore(savePath,'ReadFcn',@load);

% addpath(genpath('C:\Users\jonma\Dropbox\ResearchProjects\FlowSimulationProject\3dFlowSimulation\HelperFunctions'))
%%
playVid = true;
recordVid = true;

% flowType = 'retGrid'; %displays retinal grid flow rather than ground dot flow
flowType = 'groundDots'; %displays projected ground dots, rather than retinal grid

if exist('vidObj')
    close(vidObj)
    clear vidObj
end

if recordVid
    vidName = split(savePath,filesep);
    fileMadeTime = char(datetime(datetime,'Format','yyyy-mm-dd_HH:mm'));
    vidName = vidName{end};
    vidPath = "outputVids";%[basePath filesep 'outputVids' filesep  'FlowSim_' vidName];
    save([vidName '.mat'])

    vidObj = VideoWriter(vidPath,'MPEG-4');
    
    vidObj.FrameRate = 30;
    open(vidObj);
end
%%
temp = load(eyeStructStore.Files{1});
thisEyeStruct = temp.thisEyeStruct;
shadow_fr_mar_dim = thisEyeStruct.shadow_fr_mar_dim;
comXYZ = squeeze(shadow_fr_mar_dim(:,1,:));
steps_HS_TO_StanceLeg_XYZ = thisEyeStruct.steps_HS_TO_StanceLeg_XYZ;



%% draw stuff!

startFrame = temp.thisEyeStruct.startFrame;
endFrame = temp.thisEyeStruct.endFrame;
numFrames = endFrame-startFrame;

%%find good color scale ranges for curl and div plots
curlClim = 0.002;%we'll need to set this later
divClim = [0 .01];

curlContNum = 8; %how many lines to draw on curl contour plots?
divContNum = 8; %how many lines to draw on div contour plots?

%%% find good limit for retinal projection plots (coloreddots)
retProjPlotRLim = max(thisEyeStruct.retGridxx(:));

%% Create a circluar mask
% Next create the circle in the image.
xmaxT =max(temp.thisEyeStruct.retGridxx(:));
ymaxT =max(temp.thisEyeStruct.retGridxx(:));
nD = numel(-ymaxT:ymaxT);

xmax = retProjPlotRLim;%max(temp.thisEyeStruct.retGridxx(:));
ymax = retProjPlotRLim;%max(temp.thisEyeStruct.retGridxx(:));
[maskxx, maskyy] = meshgrid(linspace(-xmax,xmax,nD),linspace(-ymax,ymax,nD));
center = 0;
radius = retProjPlotRLim; %mae radius slightly larger than it needs to be to give a nice black outline
circleMask = (maskxx - center).^2 + (maskyy - center).^2 <= (radius).^2;

% figure;
% zc1 = imagesc(min(thisEyeStruct.retGridxx(:)) ,min(thisEyeStruct.retGridyy(:)),reshape([circleMask  circleMask circleMask], [size(circleMask) 3])*circleMaskLum);
% zc1.AlphaData = circleMask*.75;
% hold on
%
% viscircles([0,0], retProjPlotRLim, 'Color','k','EnhanceVisibility',false,'LineWidth',3);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% BIG FRAME-BY-FRAME LOOP STARTS HERE!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

looptimer = nan(size(startFrame:endFrame));

fileNum = 0;
for fr = startFrame:endFrame
    tic
    
    
    disp(strcat({'Fr#'},num2str(fr),{'-PROGRESS: '}, num2str(fr-startFrame+1),'-of-',num2str(numFrames),...
        '- time remaining ~',num2str((nanmean(looptimer) * (numFrames-fileNum))/60 ),'mins - Mean Frame Dur ~',num2str(nanmean(looptimer)),...
        '- RecordVidset to:',num2str(recordVid)));
    
    fileNum = fileNum+1;
    temp = load(eyeStructStore.Files{fileNum});
    thisEyeStruct = temp.thisEyeStruct;
    
    %%% unpack judiciously
    xx = thisEyeStruct.xx;
    yy = thisEyeStruct.yy;
    zz = thisEyeStruct.zz;
    
    
    %%%%set up figure
    f = figure(1);
    %     if ispc
    %         f.Position = [1921 30 1920 1080];
    %     elseif ismac
    %         f.Position = [1 30 1440 775];
    %     end
    f.Units = 'pixels';
    
    if strcmp(getenv('computername'),     'DESKTOP-L29LOMC')
        f.Position = [-1919 0 1920 1080];
    else
        f.Position = [0 0 1920 1080];
    end
    hold on
    clf
    
    f.Color = 'w';
    
    
    
    for curlDivIter = 1:2
        
        if curlDivIter == 1
            sp1 = axes;cla
            sp1.Position = [.03 -.25 .45 1];
            hold on
            axis equal
            sp1.Title.String = strcat({'Retinal Curl','Projected onto Groundplane'});
            sp1.Title.FontSize = 16;
            sp1.Title.Units = 'normalized';
            sp1.Title.Position = [.8 .075];
            sp1.Color = 'none';
            sp1.XTick = [];
            sp1.YTick = [];
            sp1.ZTick = [];
            shading(sp1,'flat')
            
        elseif curlDivIter == 2
            sp2 = axes;cla
            sp2.Position =     sp1.Position + [.52 0 0 0 ];
            hold on
            axis equal
            sp2.Title.String = strcat({'Retinal Divergence','Projected onto Groundplane'});
            sp2.Title.FontSize = 16;
            sp2.Title.Units = 'normalized';
            sp2.Title.Position = [.8 .075];
            sp2.Color = 'none';
            sp2.XTick = [];
            sp2.YTick = [];
            sp2.ZTick = [];
            shading(sp2,'flat')
        end
        
        
        
        thisEyeX = thisEyeStruct.eyeXYZ(fr,1);
        thisEyeY = thisEyeStruct.eyeXYZ(fr,2);
        thisEyeZ = thisEyeStruct.eyeXYZ(fr,3);
        
        
        xlim([thisEyeX-1e3 thisEyeX+8e3])
        ylim([thisEyeY-4e3 thisEyeY+4e3])
        zlim([-1e2 max(thisEyeStruct.eyeXYZ(:,3))+.5e3])
        
        view(3)
        
        % get curl (div) colors for groundplane
        
        thisRetGridIDs = thisEyeStruct.retGridIDs;
        
        thisCurlFrame = thisEyeStruct.curlFr;
        thisDivFrame = thisEyeStruct.divFr;
        
        
        
        
        % figure out the mapping between retinal curl/divergence and ground locations
        % (using the retGridID's from the calculation of the ground-retinal projection)
        thisRetGridIDsList = thisRetGridIDs(:);
        thisCurlList = thisCurlFrame(:);
        thisDivList = thisDivFrame(:);
        
        delThese = isnan(thisRetGridIDsList) | isnan(thisCurlList) | isnan(thisDivList); %get rid of NaNs
        thisRetGridIDsList(delThese) = [];
        thisCurlList(delThese) = [];
        thisDivList(delThese) = [];
        
        assert(numel(thisCurlList) == numel(thisRetGridIDsList));
        assert(numel(thisCurlList) == numel(thisDivList));
        
        %build a scatteredInterpolant for ground curl colors
        theseGroundX = xx(thisRetGridIDsList);
        theseGroundY = yy(thisRetGridIDsList);
        
        [~, idx] = unique([ theseGroundX theseGroundY],'rows','stable'); %to get rid of annoying "duplicate datapoints detected" message
        
        F_curl = scatteredInterpolant(theseGroundX(idx), theseGroundY(idx), thisCurlList(idx), 'natural','none');
        F_div = scatteredInterpolant(theseGroundX(idx), theseGroundY(idx), thisDivList(idx), 'natural','none');
        
        
        %%%plot ground plane
        if curlDivIter == 1
            groundColors = F_curl(xx,yy);%color the groundpoints in xx,yy according to the curl scattered interpolant
        elseif curlDivIter == 2
            groundColors = F_div(xx,yy);  %ditto for divergence. ScatteredInterpolant is the coolest friggin function, man. dang!
        end
        
        gb1 = surface(xx, yy, zz-90);%dark backdrop ground plane to make them colors *pop*
        gb1.FaceColor = [0 0 0]+.5;
        gb1.EdgeColor =[0 0 0]+.45;
        
        showGroundColors = true;
        if ~isempty(groundColors) && showGroundColors
            g1 = surface(xx, yy, zz, groundColors);
        else
            g1 = surface(xx, yy, zz);
            
        end
        g1.EdgeColor =  'none';
        g1.EdgeAlpha = 1;
        g1.FaceAlpha = 1;%.9;
        
        
        
        
        
        
        if curlDivIter == 1
            curlCMap = (cbrewer('div', 'RdBu', 64));
            colormap(sp1,(curlCMap));
            caxis(sp1,[-curlClim curlClim]);
        elseif curlDivIter == 2
            divCMap = (plasma);
            colormap(sp2,divCMap);
            caxis(sp2, divClim);
        end
        
        if ~showGroundColors
            g1.EdgeColor = 'k';
        end
        
        %%%draw contour lines (thick black lines == Zero crossings, greenline = foveal isoline)
        midCurl = round(size(thisCurlFrame)/2);  %middle of the curl field, aka the 'fovea'
        curlFovea = thisCurlFrame(midCurl(1),midCurl(2));
        
        midDiv = round(size(thisDivFrame)/2); %middle of the divergence field, aka the 'fovea'
        divFovea = thisDivFrame(midDiv(1),midDiv(2));
        
        if isnan(curlFovea); curlFovea = 0; end
        if isnan(divFovea); divFovea = 0; end
        
        fovContColor = [0 0.6 0.1];
        
        if fr > startFrame && showGroundColors %skip the first frame
            
            contZero = contour(xx,yy,groundColors,[0,0],'k-','LineWidth',3)';
            %                 c1 =  plot3(contZero(:,1),contZero(:,2),zLook(contZero(:,1),contZero(:,2))+10,'ko','MarkerSize',2,'MarkerFaceColor','k');
            c1 =  plot3(contZero(:,1),contZero(:,2),contZero(:,2)*0+10,'ko','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','k');
            
            if curlDivIter == 1 %curl
                contFov = contour(xx,yy,groundColors,[curlFovea curlFovea],'Color',fovContColor,'LineWidth', 2)';
                %                 contFov = contour(xx,yy,groundColors,[curlFovea curlFovea],'Color','k','LineWidth',1)';
                
                %                     c1 =  plot3(contFov(:,1),contFov(:,2),zLook(contFov(:,1),contFov(:,2))+10,'o','Color',fovContColor,'MarkerSize',2,'MarkerFaceColor',fovContColor);
                %                 c1 =  plot3(contFov(:,1),contFov(:,2),contFov(:,2)*0+10,'o','Color',fovContColor,'MarkerSize',2,'MarkerFaceColor',fovContColor);
                
                contour(xx, yy, groundColors,linspace(-curlClim, curlClim,curlContNum),'k','LineWidth',.5);
                
            elseif curlDivIter == 2 %div
                
                contFov = contour(xx,yy,groundColors,[divFovea divFovea],'Color',fovContColor,'LineWidth', 2)';
                %                 contFov = contour(xx,yy,groundColors,[divFovea divFovea],'Color','k','LineWidth',1)';
                %                     c1 =  plot3(contFov(:,1),contFov(:,2),zLook(contFov(:,1),contFov(:,2))+10,'o','Color',fovContColor,'MarkerSize',2,'MarkerFaceColor',fovContColor);
                %                 c1 =  plot3(contFov(:,1),contFov(:,2),contFov(:,2)*0+10,'o','Color',fovContColor,'MarkerSize',2,'MarkerFaceColor',fovContColor);
                
                contour(xx, yy, groundColors,linspace(divClim(1), divClim(2),divContNum),'k','LineWidth',.5);
                
            end
            
        end
        
        %%%%%draw fixated point
        thisFixPoint = thisEyeStruct.fixXYZ(fr,:);
        %         fx1 = drawPoint3d(thisFixPoint);
        %         fx1.Marker = 'o';
        %         fx1.MarkerSize = 8;
        %         fx1.MarkerFaceColor = 'w';
        %         fx1.Color = 'k';
        
        %         %%%%%draw cardinals axes
        %         eyeNorth = thisEyeStruct.eyeCardinalPointNorthGroundIntersect;
        %         eyeSouth = thisEyeStruct.eyeCardinalPointSouthGroundIntersect;
        %         eyeEast = thisEyeStruct.eyeCardinalPointEastGroundIntersect;
        %         eyeWest = thisEyeStruct.eyeCardinalPointWestGroundIntersect;
        %         plot3([eyeNorth(1) thisFixPoint(1)],[eyeNorth(2) thisFixPoint(2)],[eyeNorth(3) thisFixPoint(3)],'k-')
        
        kittyRed = [217 81 87]/255;
        
        for ff = 1:length( thisEyeStruct.fovGridRadDeg)-1
            g1 = plot3(thisEyeStruct.fovGridGroundPoints{ff}(:,1), thisEyeStruct.fovGridGroundPoints{ff}(:,2), thisEyeStruct.fovGridGroundPoints{ff}(:,3));
            
            
            if thisEyeStruct.fovGridRadDeg{ff} == 45
                g1.LineWidth = 1;
                g1.Color = 'w';
            else
                g1.LineWidth = 1;
                g1.Color = 'w';
                
            end
            
            g1.LineStyle = '-';
        end
        
        
        
        plot3(thisEyeStruct.horizontalEyeCrossHairGroundPoints(:,1),thisEyeStruct.horizontalEyeCrossHairGroundPoints(:,2), [0; 0 ],'-','LineWidth',1,'Color','w')
        plot3(thisEyeStruct.verticalEyeCrossHairGroundPoints(:,1),thisEyeStruct.verticalEyeCrossHairGroundPoints(:,2),[0; 0 ],'-','LineWidth',1,'Color','w')
        
        
        
        %%%%%%Draw maximum divergence point
        if curlDivIter == 2
            %find max divergence
            [C,I] = max(groundColors(:));
            [maxGroundDivXind, maxGroundDivYind] = ind2sub(size(groundColors),I);
            
            maxGroundDivX = xx(maxGroundDivXind,maxGroundDivYind);
            maxGroundDivY = yy(maxGroundDivXind,maxGroundDivYind);
            
            plot3(maxGroundDivX, maxGroundDivY, 0, 'kp','MarkerFaceColor','y','MarkerSize',12)
        end
        
        
        %%%% draw head path
        %         plot3(thisEyeStruct.eyeXYZ(:,1), thisEyeStruct.eyeXYZ(:,2), thisEyeStruct.eyeXYZ(:,3),'k-','LineWidth',2)%head trajectory
        %         plot3(thisEyeStruct.eyeXYZ(1,1), thisEyeStruct.eyeXYZ(1,2), thisEyeStruct.eyeXYZ(1,3),'rp') %head start position
        %         plot3(thisEyeStruct.eyeXYZ(1,1), thisEyeStruct.eyeXYZ(1,2), 0,'rp') %head vertical projection
        %
        
        %         plot3([thisEyeStruct.eyeXYZ(fr,1)-1e4 thisEyeStruct.eyeXYZ(fr,1)+1e4],[0 0],[15 15],'k--','LineWidth',2)
        %         plot3([thisEyeStruct.eyeXYZ(1,1) thisEyeStruct.eyeXYZ(1,1)],[-1e4 1e4],[15 15],'k','LineWidth',2)
        
        
        
        
        lLeg = [2 3 4 5 6 7 5];
        rLeg = [2 8 9 10 11 12 10];
        tors = [2 13 14 15 26 27 28];
        lArm = [15 16 17 26 17 18 19 20];
        rArm = [15 21 22 26 22 23 24 25];
        
        %%%Plotcherself up a nice little skeleton friend
        plot3(shadow_fr_mar_dim(fr,1:28,1),shadow_fr_mar_dim(fr,1:28,2),shadow_fr_mar_dim(fr,1:28,3),'ko','MarkerFaceColor','k','MarkerSize',4)
        hold on
        
        
        plot3(shadow_fr_mar_dim(fr,lLeg,1),shadow_fr_mar_dim(fr,lLeg,2),shadow_fr_mar_dim(fr,lLeg,3),'c','LineWidth',2)
        plot3(shadow_fr_mar_dim(fr,rLeg,1),shadow_fr_mar_dim(fr,rLeg,2),shadow_fr_mar_dim(fr,rLeg,3),'r','LineWidth',2)
        plot3(shadow_fr_mar_dim(fr,tors,1),shadow_fr_mar_dim(fr,tors,2),shadow_fr_mar_dim(fr,tors,3),'g','LineWidth',2)
        plot3(shadow_fr_mar_dim(fr,lArm,1),shadow_fr_mar_dim(fr,lArm,2),shadow_fr_mar_dim(fr,lArm,3),'c','LineWidth',2)
        plot3(shadow_fr_mar_dim(fr,rArm,1),shadow_fr_mar_dim(fr,rArm,2),shadow_fr_mar_dim(fr,rArm,3),'r','LineWidth',2)
        
        
        %%% plot foothold locations
        rFootholds = steps_HS_TO_StanceLeg_XYZ(steps_HS_TO_StanceLeg_XYZ(:,3) == 1 ,:);
        lFootholds = steps_HS_TO_StanceLeg_XYZ(steps_HS_TO_StanceLeg_XYZ(:,3) == 2 ,:);
        
        rFootholds(rFootholds(:,1)<fr-2000 | rFootholds(:,1)>fr+2000,:) = [];
        lFootholds(lFootholds(:,1)<fr-2000 | lFootholds(:,1)>fr+2000,:) = [];
        
        %   plot vertical projection of foothold locations onto groundplane
        
        plot3(rFootholds(:,4),  rFootholds(:,5), rFootholds(:,6),'ko','MarkerSize', 6, 'MarkerFaceColor','r')
        plot3(lFootholds(:,4), lFootholds(:,5), lFootholds(:,6), 'ko','MarkerSize', 6, 'MarkerFaceColor','c')
        
        
        
        
        
        
        
        
        
        
        eyeCol = 'k';
        
        
        
        
        if fr == 1
            eyeVel = [ 0 0 0 ];
        else
            eyeVel = thisEyeStruct.eyeXYZ(fr,:) - thisEyeStruct.eyeXYZ(fr-1,:);
        end
        
        
        quiver3(thisEyeStruct.eyeXYZ(fr,1), thisEyeStruct.eyeXYZ(fr,2), 15, eyeVel(1)*500, eyeVel(2)*500, 0, 'w','LineWidth',3)%eye instantaneous velocity vector
        quiver3(thisEyeStruct.eyeXYZ(fr,1), thisEyeStruct.eyeXYZ(fr,2), 15, eyeVel(1)*500, eyeVel(2)*500, 0, eyeCol,'LineWidth',1.5)%eye instantaneous velocity vector
        
        %%%%draw region from which we're allowing projections
        theseFovGroundPoints =  thisEyeStruct.fovGroundPoints;
        
        %                 plot3(theseFovGroundPoints(:,1), theseFovGroundPoints(:,2), zLook(theseFovGroundPoints(:,1), theseFovGroundPoints(:,2)),'k:','LineWidth',2)
        plot3(theseFovGroundPoints(:,1), theseFovGroundPoints(:,2), theseFovGroundPoints(:,2)*0,'-','Color',kittyRed,'LineWidth',6)
        
        thisEyeXYZ = thisEyeStruct.eyeXYZ;
        
        %%% draw gaze line (between eye center and fixated point)
        gl1 = drawEdge3d([thisEyeStruct.pupilXYZ thisEyeStruct.fixXYZ(fr,:)]);
        gl1.Color = 'm';
        gl1.LineWidth = 2;
        gl2 = drawEdge3d([thisEyeXYZ(fr,:) thisEyeStruct.fixXYZ(fr,:)]); %draw gaze line again in blue for visibility
        
        
        
        
        
        %%%%draw projected points on eyeball
        %         theseBackEyeProjPoints_xyz_id = thisEyeStruct.backEyeProjPoints_xyz_id;
        %         ep1 = scatter3(theseBackEyeProjPoints_xyz_id(:,1), theseBackEyeProjPoints_xyz_id(:,2), theseBackEyeProjPoints_xyz_id(:,3), 'filled'); %all intersections
        %         ep1.SizeData = 1;
        %         ep1.MarkerEdgeColor = 'k';
        %         ep1.MarkerEdgeAlpha = .1;
        %                 ep1.CData = dotColors(theseBackEyeProjPoints_xyz_id(:,4),:);
        
        
        %%%%%draw eyeball
        %     e1 = drawSphere(eyeSphere);
        thisEyeMesh = thisEyeStruct.eyeMesh;
        e1 = drawMesh(thisEyeMesh);
        
        
        e1.FaceColor = 'c';%[.5 .5 .5];
        e1.EdgeColor = 'none';
        e1.FaceAlpha = .8;
        
        
        drawPoint3d(thisEyeXYZ(fr,:),'yp')
        
        
        %%% draw circle around midpoint of eye sphere
        ec1 = drawCircle3d(thisEyeStruct.coronalEyeCircle,   'LineWidth', 2,'Color', 'r');
        es1 = drawCircle3d(thisEyeStruct.saggitalEyeCircle,  'LineWidth', 2,'Color', 'g');
        et1 = drawCircle3d(thisEyeStruct.transverseEyeCircle,'LineWidth', 2,'Color', 'b');
        
        %%% draw eye axes
        drawVector3d(thisEyeXYZ(fr,:) ,thisEyeStruct.eyeXvec*thisEyeStruct.eyeRadius*1.5,'Color','r','LineWidth',3)
        drawVector3d(thisEyeXYZ(fr,:) ,thisEyeStruct.eyeYvec*thisEyeStruct.eyeRadius*1.5,'Color','g','LineWidth',3)
        drawVector3d(thisEyeXYZ(fr,:) ,thisEyeStruct.eyeZvec*thisEyeStruct.eyeRadius*1.5,'Color','b','LineWidth',3)
        
        if ~skelBool
            ex1 = plot3(thisEyeXYZ(:,1), thisEyeXYZ(:,2), thisEyeXYZ(:,3),'b-','LineWidth',1.5);
            exz = plot3(thisEyeXYZ(:,1), thisEyeXYZ(:,2), zeros(size(thisEyeXYZ(:,3))),'b-','LineWidth',1);
        end           
            
        
        
        %%% pupil
        pup1 = drawSphere([thisEyeStruct.pupilXYZ thisEyeStruct.eyeRadius*.1]);
        pup1.FaceColor = 'k';
        pup1.FaceAlpha = e1.FaceAlpha;
        
        %%% fovea
        fov1 = drawSphere([thisEyeStruct.foveaXYZ thisEyeStruct.eyeRadius*.05]);
        fov1.FaceColor = 'g';
        
        
        
        
        
        %         if fr == startFrame
        box off
        axis off
        %         end
        
        %     lighting gouraud
        %     l = light;
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Show Retinal Projection Polar plot (Plots the Dots)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    po1 = axes;cla
    po1.Position = [   0.3354    0.33   0.3292    0.5852];
    hold on
    %     po1.Title.String = strcat({'Retinal Projection Flow (FOV radius '},num2str(thisEyeStruct.fovRadDeg),...
    %         {' Degrees) - Frame # '},num2str(fr));%, {' - EyeXYZ = '},num2str(thisEyeXYZ(fr,:),'%.1f'));
    
    
    po1.XTickLabel = [];
    po1.YTickLabel = [];
    
    po1.XLim = [-retProjPlotRLim retProjPlotRLim];
    po1.YLim = [-retProjPlotRLim retProjPlotRLim];
    
    
    
    
    
    
    switch flowType
        case 'groundDots' %show dots on the ground projected onto retina
            
            dotX     = thisEyeStruct.retinaProjPoints_fr_X_Y_id(:,1);
            dotY     = thisEyeStruct.retinaProjPoints_fr_X_Y_id(:,2);
            dotVelX  = thisEyeStruct.retinaProjVel_fr_xVel_yVel_id(:,1);
            dotVelY  = thisEyeStruct.retinaProjVel_fr_xVel_yVel_id(:,2);
            
        case 'retGrid' %show retinal grid flow
            
            dotX     = thisEyeStruct.retGridxx;
            dotY     = thisEyeStruct.retGridyy;
            dotVelX  = thisEyeStruct.retGridVelxx;
            dotVelY  = thisEyeStruct.retGridVelyy;
            
    end
    
    %don't draw dottos that cross over retinal limit
    dotX(sqrt((dotX+dotVelX).^2 + (dotY+dotVelY).^2) > retProjPlotRLim*.96) = nan;
    dotY(sqrt((dotX+dotVelX).^2 + (dotY+dotVelY).^2) > retProjPlotRLim*.96) = nan;
    
    circleMaskLum = .75;
    circleMaskAlpha = .5;
    zdc1 = imagesc(min(thisEyeStruct.retGridxx(:)) ,min(thisEyeStruct.retGridyy(:)),reshape([circleMask  circleMask circleMask], [size(circleMask) 3])*circleMaskLum);
    zdc1.AlphaData = circleMask*circleMaskAlpha/5;
    hold on
    
    plot(dotX, dotY,'.','Color', eyeCol)
    flowScaleFactor = 15;
    
    q1 = quiver(dotX, dotY, dotVelX*flowScaleFactor, dotVelY*flowScaleFactor, eyeCol, 'AutoScale','off');
    q1.LineWidth = 1;
    q1.MaxHeadSize = .5;
    
    eyeCol = 'k';
    
    
    showStreamLines = true; %don't show streamlines if you've got two eyes (it makes things too cluttered)
    if showStreamLines
        theseStreamXY =  thisEyeStruct.streamXY;
        for ss = 1:length(theseStreamXY)
            if norm(theseStreamXY{ss}(1,:)) < retProjPlotRLim %don't draw streams that start outside of FOV limit
                streamColPurp = [    0.6    0.    0.9560];
                plot(theseStreamXY{ss}(:,1), theseStreamXY{ss}(:,2),'-','Color',streamColPurp, 'LineWidth',2)
            end
        end
        %             plot( thisEyeStruct.FOExy(fr,1), thisEyeStruct.FOExy(fr,2),'kp','MarkerSize',18,'MarkerFaceColor','y')
        
    end
    
    
    
    
    
    
    viscircles([0,0], retProjPlotRLim, 'Color',kittyRed,'EnhanceVisibility',false,'LineWidth',6);
    viscircles([0,0], retProjPlotRLim/4, 'Color','k','EnhanceVisibility',false,'LineWidth',1);
    viscircles([0,0], retProjPlotRLim/2, 'Color','k','EnhanceVisibility',false,'LineWidth',1);
    viscircles([0,0], 3*retProjPlotRLim/4, 'Color','k','EnhanceVisibility',false,'LineWidth',1);
    
    plot([0 0],[-retProjPlotRLim retProjPlotRLim],'k')
    plot([-retProjPlotRLim retProjPlotRLim],[0 0],'k')
    
    %     plot(po1, 0, 0, 'k^','MarkerFaceColor','w')
    
    %%make little cute little triangle guy like "yo, friend, this way up ^-^"
    triSize = diff(xlim)/50;
    triColor = 'w';
    tx = [1 0 -1 1]*triSize/2;
    ty = [1 -1 1 1]*triSize/2;
    tz = [0 0 0 0]*triSize/2;
    tri = patch(tx, ty, tz,triColor);
    
    
    axis off
    axis ij
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Show  Curl plot (just for the left eye, for now)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cu1 = axes; cla
    cu1.Position = [-0.02 .6 .33 .33];
    
    gridDist = sqrt(thisEyeStruct.retGridxx.^2 + thisEyeStruct.retGridyy.^2);
    outOfRangeDots = gridDist>retProjPlotRLim;
    thisEyeStruct.retGridxx( outOfRangeDots) = nan;
    
    zc1 = imagesc(min(thisEyeStruct.retGridxx(:)) ,min(thisEyeStruct.retGridyy(:)),reshape([circleMask  circleMask circleMask], [size(circleMask) 3])*circleMaskLum);
    zc1.AlphaData = circleMask*circleMaskAlpha;
    hold on
    
    c1 = surface(thisEyeStruct.retGridxx, thisEyeStruct.retGridyy, zeros(size(thisEyeStruct.retGridxx)),thisCurlFrame);
    c1.EdgeColor = 'none';
    
    if sum(isnan(thisCurlFrame(:))) ~= numel(thisCurlFrame)
        contour(thisEyeStruct.retGridxx, thisEyeStruct.retGridyy, thisCurlFrame,[0 0],'k','LineWidth',3);
        contour(thisEyeStruct.retGridxx, thisEyeStruct.retGridyy, thisCurlFrame,linspace(-curlClim, curlClim,curlContNum),'k','LineWidth',.5);
        
        contour(thisEyeStruct.retGridxx, thisEyeStruct.retGridyy, thisCurlFrame,[curlFovea curlFovea],'Color',fovContColor,'LineWidth', 2);
        %         contour(thisEyeStruct.retGridxx, thisEyeStruct.retGridyy, thisCurlFrame,[curlFovea curlFovea],'Color','k','LineWidth',.5);
        
    end
    
    curlZeroCrossings = contourc(thisCurlFrame,[0 0])';
    curlZeroCrossings(curlZeroCrossings(:,1) ==0 |curlZeroCrossings(:,2) ==0, :) = [];
    
    colormap(cu1, curlCMap);
    
    caxis([-curlClim curlClim]);
    
    ccr = colorbar;
    ccr.Position =  [0.2519 0.6 0.0139 0.3307];
    
    viscircles([0,0], retProjPlotRLim, 'Color',kittyRed,'EnhanceVisibility',false,'LineWidth',4);
    %     viscircles([0,0], retProjPlotRLim, 'Color','k','EnhanceVisibility',false,'LineWidth',1);
    viscircles([0,0], retProjPlotRLim/4, 'Color','w','EnhanceVisibility',false,'LineWidth',1);
    viscircles([0,0], retProjPlotRLim/2, 'Color','w','EnhanceVisibility',false,'LineWidth',1);
    viscircles([0,0], 3*retProjPlotRLim/4, 'Color','w','EnhanceVisibility',false,'LineWidth',1);
    
    plot([0 0],[-retProjPlotRLim retProjPlotRLim],'w')
    plot([-retProjPlotRLim retProjPlotRLim],[0 0],'w')
    %     plot(cu1, 0, 0, 'k^','MarkerFaceColor','w')
    %%make little cute little triangle guy like "yo, friend, this way up ^-^"
    tx = [1 0 -1 1]*triSize;
    ty = [1 -1 1 1]*triSize;
    tz = [0 0 0 0]*triSize;
    tri = patch(tx, ty, tz,triColor);
    
    
    xlim([-retProjPlotRLim retProjPlotRLim]);
    ylim([-retProjPlotRLim retProjPlotRLim]);
    
    axis equal
    axis off
    axis ij
    
    
    cu1.Title.String = strcat({'Retinal Curl'});%,'(subtractive normalization)'});
    
    cu1.Title.FontSize = 20;
    cu1.Title.FontWeight = 'normal';
    
    %%%% add lines denoting current value range
    maxCurl = min([curlClim max(thisCurlFrame(:))]);%either the largest Curl value, or the max value of the colormap, whichever is smallest
    minCurl = max([-curlClim min(thisCurlFrame(:))]); %simlarly for the bottom
    
    h_axes = axes('position', ccr.Position, 'ylim', ccr.Limits, 'color', 'none', 'visible','off');
    line(h_axes.XLim, maxCurl*[1 1], 'color', 'm', 'LineWidth', 2, 'parent', h_axes);
    line(h_axes.XLim, minCurl*[1 1], 'color', 'm', 'LineWidth', 2, 'parent', h_axes);
    line(h_axes.XLim, curlFovea*[1 1], 'color', fovContColor, 'LineWidth', 3, 'parent', h_axes);
    line(h_axes.XLim, [0 0], 'color', 'k', 'LineWidth', 2, 'parent', h_axes);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Show  Curl histogram plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ch1 = axes; cla
    ch1.Position = [0.195 0.44 0.125 0.15];
    nBins = 100;
    cHist1 = histogram(thisCurlFrame,linspace(-curlClim*2,curlClim*2,nBins), 'Normalization','probability');
    hold on
    
    plot([0 0], [0 1],'k-')
    
    curlYlim = .1;
    
    plot([0 0]+curlFovea, [0 curlYlim*.1],'-','Color',fovContColor,'LineWidth',3)
    plot([0 0]+maxCurl, [0 curlYlim*.1],'-','Color','m','LineWidth',2)
    plot([0 0]+minCurl, [0 curlYlim*.1],'-','Color','m','LineWidth',2)
    
    ch1.YLim = [0 curlYlim];
    ch1.XLim = [-curlClim curlClim]*2;
    
    
    ch1.YLabel.String = 'Probability';
    ch1.XLabel.String = 'Curl';
    box off
    hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Show  Curl surface plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cs1 = axes; cla
    cs1.Position = [0.005 0.25 0.1706 0.4];
    
    
    
    
    
    eyeCol = 'k';
    
    
    thisCurlFrame = thisEyeStruct.curlFr;
    
    %     thisCurlFrame(outOfRangeDots) = nan;
    
    zc1 = imagesc(min(thisEyeStruct.retGridxx(:)) ,min(thisEyeStruct.retGridyy(:)),reshape([circleMask  circleMask circleMask], [size(circleMask) 3])*circleMaskLum);
    zc1.AlphaData = circleMask*circleMaskAlpha;
    hold on
    contZero = contour(thisEyeStruct.retGridxx, thisEyeStruct.retGridyy, thisCurlFrame,[0,0],'k-','LineWidth',3)';
    contFov = contour(thisEyeStruct.retGridxx, thisEyeStruct.retGridyy, thisCurlFrame,[curlFovea curlFovea],'Color',fovContColor,'LineWidth', 2)';
    contour(thisEyeStruct.retGridxx, thisEyeStruct.retGridyy, thisCurlFrame,linspace(-curlClim, curlClim,curlContNum),'k','LineWidth',.5);
    
    
    viscircles([0,0], retProjPlotRLim, 'Color',kittyRed,'EnhanceVisibility',false,'LineWidth',4);
    viscircles([0,0], retProjPlotRLim/4, 'Color','w','EnhanceVisibility',false,'LineWidth',1);
    viscircles([0,0], retProjPlotRLim/2, 'Color','w','EnhanceVisibility',false,'LineWidth',1);
    viscircles([0,0], 3*retProjPlotRLim/4, 'Color','w','EnhanceVisibility',false,'LineWidth',1);
    
    plot([0 0],[-retProjPlotRLim retProjPlotRLim],'w')
    plot([-retProjPlotRLim retProjPlotRLim],[0 0],'w')
    %     plot(0,0,'k^','MarkerFaceColor','w')
    
    %%make little cute little triangle guy like "yo, friend, this way up ^-^"
    tx = [.5 0 -.5 .5]*5*triSize;
    ty = [.5 -.5 .5 .5]*5*triSize-.25;
    tz = [0 0 0 0]+curlClim/50;
    cs1tri = patch(tx, ty, tz,triColor);
    
    sc1 = surface(thisEyeStruct.retGridxx,thisEyeStruct.retGridyy,thisCurlFrame);
    sc1.EdgeColor = 'none';
    
    contZero = contour3(thisEyeStruct.retGridxx, thisEyeStruct.retGridyy, thisCurlFrame,[0,0],'k-','LineWidth',3)';
    contFov = contour3(thisEyeStruct.retGridxx, thisEyeStruct.retGridyy, thisCurlFrame,[curlFovea curlFovea],'Color',fovContColor,'LineWidth', 2)';
    contour3(thisEyeStruct.retGridxx, thisEyeStruct.retGridyy, thisCurlFrame,linspace(-curlClim, curlClim,curlContNum),'k','LineWidth',.5);
    
    
    az = 60;
    el = 10;
    cs1.View = [az el];
    %     curlLight = camlight;
    %     curlLight.Position = [ 0 0 2000];
    caxis(cs1, [-curlClim curlClim]);
    cs1.XLim= [-retProjPlotRLim retProjPlotRLim];
    cs1.YLim = [-retProjPlotRLim retProjPlotRLim];
    cs1.ZLim = [-curlClim curlClim]*3;
    axis off
    colormap(cs1,(curlCMap));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Show Divergence plot (just for the left eye, for now)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    div1 = axes; cla;
    div1.Position = [0.66 .6 .33 .33];
    
    
    
    %find max divergence
    [C,I] = max(thisDivFrame(:));
    [maxDivXind, maxDivYind] = ind2sub(size(thisDivFrame),I);
    
    maxDivX = squeeze(thisEyeStruct.retGridxx(maxDivXind,maxDivYind));
    maxDivY = squeeze(thisEyeStruct.retGridyy(maxDivXind,maxDivYind));
    
    zc1 = imagesc(min(thisEyeStruct.retGridxx(:)) ,min(thisEyeStruct.retGridyy(:)),reshape([circleMask  circleMask circleMask], [size(circleMask) 3])*circleMaskLum);
    zc1.AlphaData = circleMask*circleMaskAlpha;
    hold on
    
    d1 = surface(squeeze(thisEyeStruct.retGridxx), squeeze(thisEyeStruct.retGridyy), zeros(size(squeeze(thisEyeStruct.retGridyy))),thisDivFrame);
    d1.EdgeColor = 'none';
    d1.AlphaData = circleMask;
    
    
    
    if sum(isnan(thisDivFrame(:))) ~= numel(thisDivFrame)
        
        contour(thisEyeStruct.retGridxx, thisEyeStruct.retGridyy, thisDivFrame,[0 0],'k','LineWidth',3);
        
        contour(thisEyeStruct.retGridxx, thisEyeStruct.retGridyy, thisDivFrame,linspace(divClim(1), divClim(2),divContNum), 'k','LineWidth',.5);
        
        contour(thisEyeStruct.retGridxx, thisEyeStruct.retGridyy, thisDivFrame,[divFovea divFovea],'Color',fovContColor,'LineWidth', 2);
        %         contour(thisEyeStruct.retGridxx, thisEyeStruct.retGridyy, thisDivFrame,[divFovea divFovea],'Color','k','LineWidth',.5);
    end
    
    colormap(div1,divCMap);
    caxis(divClim);
    
    cdv = colorbar;
    cdv.Position = [0.9439 0.6 0.0139 0.3307];
    
    plot(maxDivX, maxDivY,'kp','MarkerFaceColor','y','MarkerSize',12)
    
    
    
    %     viscircles([0,0], retProjPlotRLim, 'Color','k','EnhanceVisibility',false,'LineWidth',1);
    viscircles([0,0], retProjPlotRLim, 'Color',kittyRed,'EnhanceVisibility',false,'LineWidth',4);
    viscircles([0,0], retProjPlotRLim/4, 'Color','w','EnhanceVisibility',false,'LineWidth',1);
    viscircles([0,0], retProjPlotRLim/2, 'Color','w','EnhanceVisibility',false,'LineWidth',1);
    viscircles([0,0], 3*retProjPlotRLim/4, 'Color','w','EnhanceVisibility',false,'LineWidth',1);
    
    plot([0 0],[-retProjPlotRLim retProjPlotRLim],'w')
    plot([-retProjPlotRLim retProjPlotRLim],[0 0],'w')
    
    
    %%make little cute little triangle guy like "yo, friend, this way up ^-^"
    tx = [1 0 -1 1]*triSize;
    ty = [1 -1 1 1]*triSize;
    tz = [0 0 0 0]*triSize;
    tri = patch(tx, ty, tz,triColor);
    
    
    
    
    xlim([-retProjPlotRLim retProjPlotRLim]);
    ylim([-retProjPlotRLim retProjPlotRLim]);
    axis equal
    axis off
    axis ij
    div1.Title.String = strcat({'Retinal Divergence' });
    
    div1.Title.FontSize = 20;
    div1.Title.FontWeight = 'normal';
    %%%% add lines denoting current value range
    
    maxDiv = min([divClim(2) max(thisDivFrame(:))]);%either the largest Div value, or the max value of the colormap, whichever is smallest
    minDiv = max([divClim(1) min(thisDivFrame(:))]); %simlarly for the bottom
    
    
    h_axes = axes('position', cdv.Position, 'ylim', cdv.Limits, 'color', 'none', 'visible','off');
    line(h_axes.XLim, maxDiv*[1 1], 'color', 'm', 'LineWidth', 2, 'parent', h_axes);
    line(h_axes.XLim, minDiv*[1 1], 'color', 'm', 'LineWidth', 2,'parent', h_axes);
    line(h_axes.XLim, divFovea*[1 1], 'color', fovContColor, 'LineWidth', 3,'parent', h_axes);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Show  Divergence histogram plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dh1 = axes; cla
    dh1.Position = [0.865 0.44 .125 .15];
    
    cDiv1 = histogram(thisDivFrame, linspace(divClim(1), divClim(2), 60), 'Normalization','probability');
    
    hold on
    
    divYlim = .07;
    
    hdg1 = plot([0 0]+divFovea, [0 divYlim*.1],'-','Color',fovContColor,'LineWidth',3);
    plot([0 0]+maxDiv, [0 curlYlim*.1],'-','Color','m','LineWidth',2)
    plot([0 0]+minDiv, [0 divYlim*.1],'-','Color','m','LineWidth',2)
    
    dh1.YLim = [0 divYlim];
    dh1.XLim = divClim*1.2;
    dh1.YLabel.String = 'Probability';
    dh1.XLabel.String = 'Divergence';
    box off
    hold on
    colormap(dh1,(divCMap));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Show  Divergence surface plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ds1 = axes; cla
    ds1.Position = [0.678 0.409 0.1706 0.4];
    
    eyeCol = 'k';
    
    %     thiDivFrame(outOfRangeDots) = nan;
    
    
    zc1 = imagesc(min(thisEyeStruct.retGridxx(:)) ,min(thisEyeStruct.retGridyy(:)),reshape([circleMask  circleMask circleMask], [size(circleMask) 3])*circleMaskLum);
    zc1.AlphaData = circleMask*circleMaskAlpha;
    hold on
    
    
    contZero = contour(thisEyeStruct.retGridxx, thisEyeStruct.retGridyy, thisDivFrame,[0,0],'k-','LineWidth',3)';
    contFov = contour(thisEyeStruct.retGridxx, thisEyeStruct.retGridyy, thisDivFrame,[divFovea divFovea],'Color',fovContColor,'LineWidth', 2)';
    contour(thisEyeStruct.retGridxx, thisEyeStruct.retGridyy, thisDivFrame,linspace(divClim(1), divClim(2), divContNum),'k','LineWidth',.5);
    
    
    viscircles([0,0], retProjPlotRLim, 'Color',kittyRed,'EnhanceVisibility',false,'LineWidth',4);
    viscircles([0,0], retProjPlotRLim/4, 'Color','w','EnhanceVisibility',false,'LineWidth',1);
    viscircles([0,0], retProjPlotRLim/2, 'Color','w','EnhanceVisibility',false,'LineWidth',1);
    viscircles([0,0], 3*retProjPlotRLim/4, 'Color','w','EnhanceVisibility',false,'LineWidth',1);
    
    plot([0 0],[-retProjPlotRLim retProjPlotRLim],'w')
    plot([-retProjPlotRLim retProjPlotRLim],[0 0],'w')
    
    %%make little cute little triangle guy like "yo, friend, this way up ^-^"
    tx = [.5 0 -.5 .5]*5*triSize;
    ty = [.5 -.5 .5 .5]*5*triSize-.25;
    tz = [0 0 0 0]+.0001;
    ds1tri = patch(tx, ty, tz,triColor);
    
    sd1 = surface(thisEyeStruct.retGridxx, thisEyeStruct.retGridyy, thisDivFrame);
    sd1.EdgeColor = 'none';
    
    
    
    
    
    contZero = contour3(thisEyeStruct.retGridxx, thisEyeStruct.retGridyy, thisDivFrame,[0,0],'k-','LineWidth',3)';
    contFov = contour3(thisEyeStruct.retGridxx, thisEyeStruct.retGridyy, thisDivFrame,[divFovea divFovea],'Color',fovContColor,'LineWidth', 2)';
    contour3(thisEyeStruct.retGridxx, thisEyeStruct.retGridyy, thisDivFrame,linspace(divClim(1), divClim(2), divContNum),'k','LineWidth',.5);
    
    ds1.View = [az, el];
    %     divLight = camlight;
    %     divLight.Position = [ 0 0 2000];
    caxis(ds1, divClim);
    ds1.XLim= [-retProjPlotRLim retProjPlotRLim];
    ds1.YLim = [-retProjPlotRLim retProjPlotRLim];
    ds1.ZLim = divClim*2.5;
    axis off
    
    
    colormap(ds1,(divCMap));
    caxis(ds1, divClim);
    
    
    
    %%        display frame info
    w = thisEyeStruct.baseEyeStruct;
    w.subID = w.sessionID(end-2:end);
    
    
    %%
    %%% infor text box
    ex(fr) = (thisEyeX-thisFixPoint(1))/1000;
    ey(fr) = (thisEyeY-thisFixPoint(2))/1000;
    ez(fr) = (thisEyeZ-thisFixPoint(3))/1000;
    
    if fr == startFrame
        dex(fr) = 0;         dey(fr) = 0;         dez(fr) = 0;
    else
        dex(fr) = ex(fr)-ex(fr-1)*w.samplingRate;  dey(fr) = ey(fr)-ey(fr-1)*w.samplingRate;   dez(fr) = ez(fr)-ez(fr-1)*w.samplingRate;
    end
    
    a1 = axes('Position',[0 0 1 1]);
    a1.Box = 'off';
    a1.Color = 'none';
    a1.XTick = [];
    a1.YTick = [];
    i1 = text(0,0,{...
        ['SessionID: ' w.sessionID]; ...
        ['TrialType: ' w.trialType];...
        ['Walk # ' num2str(w.walkNum)];...
        ['StartFrame: ' num2str(startFrame)];...
        ['EndFrame: ' num2str(endFrame) ];...
        ['DataFrame# ' num2str(fr)];...
        ['Frame# ' num2str(fr-startFrame) ' of ' num2str(numFrames)];...
        }, 'Interpreter','none');
    i1.Position = [0.005 .935];
    i1.FontSize = 10;
    %
    %             ['Eye-Fix(m): [' num2str(ex(fr),'%.1f') ' ' num2str(ey(fr),'%.1f') ' ' num2str(ex(fr),'%.1f') ']'];...
    %         ['Eye Vel(m/s): [' num2str(dex(fr),'%.1f') ' ' num2str(dey(fr),'%.1f') ' ' num2str(dex(fr),'%.1f') ']'
    
    tit = text; %main title
    tit2 = text;%subtitle
    tit3 = text;%subsubtitle
    tit4 = text;%subsubsubtitle
    
    %%
    tit.String = 'Simulated Retinal Flow';%, {' - EyeXYZ = '},num2str(thisEyeXYZ(fr,:),'%.1f'));
    tit.FontSize = 32;
    tit.Position =[.5 .978];
    tit.HorizontalAlignment = 'center';
    
    tit2.Position = tit.Position + [0 -.03 0 ];
    tit2.String = (['Field of View Radius: ' num2str(thisEyeStruct.fovRadDeg),' degrees']);
    tit2.FontSize = 12;
    tit2.HorizontalAlignment = 'center';
    
    
    
    tit3.Position = tit2.Position + [0 -.0182 0 ];
    tit3.String =  ['eyeXYZ - fixationXYZ: [' num2str(ex(fr),'%.2f') ', ' num2str(ey(fr),'%.2f') ', ' num2str(ez(fr),'%.2f') '] (m)'];
    %           ' | Eye Vel(m/s): [' num2str(dex(fr),'%.1f') ' ' num2str(dey(fr),'%.1f') ' ' num2str(dex(fr),'%.1f') ']'];
    tit3.HorizontalAlignment = 'center';
    %         tit3.FontSize = 12;
    
    %
    %             tit4.Position = tit3.Position + [0 -.02 0 ];
    %             tit4.String = ['Eye Velocity relative to Fixation [dx dy dz](mm/s): [' num2str([dex(fr) dey(fr) dez(fr)]),']'];
    %             tit4.HorizontalAlignment = 'center';
    %             tit4.FontSize = 12;
    %%
    drawnow
    
    
    %     if recordVid
    %         im_out = frame2im(getframe(gcf));
    %         writeVideo(vidObj,im_out);
    %     end

        f.Position = [0 0 1920 1080];
    if recordVid
        im_out = frame2im(getframe(gcf));
        
        writeVideo(vidObj,im_out);
    end
    
    
    % if recordVid
    %     %         im_out = frame2im(getframe(gcf));
    %     print('thisFrame', '-dpng', '-r300')
    %
    %     if fr == 1
    %         print('firstFrame', '-dpng', '-r300')
    %     end
    %
    %     im = imread('thisFrame.png');
    %
    %     writeVideo(vidObj,im);
    % end
    looptimer(end+1) = toc;
    
end


%%
if recordVid
    close(vidObj)
end
%     end


