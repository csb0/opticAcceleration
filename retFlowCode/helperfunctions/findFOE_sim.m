
function [FOExy, streamXY] = findFOE_sim(retGridxx, retGridyy,retGridVelxx, retGridVelyy, downResFactor, debug, spotcheck)




numRows = 1;
numCols = 2;
%%

minX = min(retGridxx(:));
maxX = max(retGridxx(:));
minY = min(retGridyy(:));
maxY = max(retGridyy(:));



%% Invert flow field
nvx = -retGridVelxx;
nvy = -retGridVelyy;

dwnRes = round(linspace(1,length(retGridxx(1,:)),round(length(retGridxx(1,:)))/downResFactor));


%in the simulation version, we just pick the starting points to be retinal grid
%remove start points where this is no retinal velocity data
startx = retGridxx(dwnRes,dwnRes);
starty = retGridyy(dwnRes,dwnRes);
nodata = isnan(retGridVelxx(dwnRes,dwnRes));
startx = startx(~nodata);
starty = starty(~nodata);
%%

stepsize = .01;

streamXY = stream2(retGridxx, retGridyy, nvx, nvy, startx(:), starty(:), stepsize);
%%
streamEndsXY = nan(length(streamXY),2);
delThese = nan(length(streamXY),1);
for ss = 1:length(streamXY)
    thisXY = streamXY{ss};
    
    
    
    if isnan(thisXY(end,1))
        delThese(ss) = true;
    else
        delThese(ss) = false;
        %         streamEndsXY(ss,:) = mean(thisXY(100:end,:));
        streamEndsXY(ss,:) = thisXY(end,:);
        
    end
    
end

%%
FOExy = nanmean(streamEndsXY);


if debug
    figure(99984)
    clf
    s1 = streamline(streamXY);
    
    for ll = 1:length(s1)
        s1(ll).LineWidth = .5;
        s1(ll).Color = [0 0 0]+.5;
    end
    hold on
    q1 = quiver(retGridxx, retGridyy, retGridVelxx, retGridVelyy);
    hold on
    q1.MaxHeadSize = .001;
    q1.Marker = '.';
    q1.MarkerEdgeColor = 'k';
    
    q1.Color = 'r';
    q1.LineWidth = 1;
    
    
    plot(FOExy(1), FOExy(2),'kp','MarkerFaceColor','y','MarkerSize',24)
    plot(streamEndsXY(:,1), streamEndsXY(:,2),'rp')
    
    plot([minX maxX],[0 0],'k-','LineWidth',1)
    plot([0 0],[minY maxY],'k-','LineWidth',1)
    
    axis equal
    axis([minX maxX minY maxY ])
    
    if spotcheck
        dbstack
        keyboard
    end
end

