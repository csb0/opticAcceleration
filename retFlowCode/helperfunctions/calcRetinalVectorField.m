function     [retGridxx, retGridyy, retGridVelxx, retGridVelyy, retGridIDs] = calcRetinalVectorField(thisEyeStruct, retGridRes, gridRadius, fr, debug, spotcheck);


varNames = fieldnames(thisEyeStruct);
for i=1:length(varNames)
    eval([varNames{i} '=thisEyeStruct.' varNames{i} ';']);
end



%% the plan here is to use ScatteredInterpolant to determin the flow values
%along a specified retinal grid based on the dottos that intersect the back of the eye

dotX = retinaProjPoints_fr_X_Y_id(:,1);
dotY = retinaProjPoints_fr_X_Y_id(:,2);

%%
dotID =retinaProjPoints_fr_X_Y_id(:,3);

dotVelX = retinaProjVel_fr_xVel_yVel_id(:,1);
dotVelY = retinaProjVel_fr_xVel_yVel_id(:,2);


%% get rid of NaNs, as those screw with the scattered interpolant stuff downstream

delThese =isnan(dotX);
dotX(delThese) = [];
dotY(delThese) = [];
dotVelX(delThese) = [];
dotVelY(delThese) = [];
dotID(delThese) = [];

delThese =isnan(dotVelX);
dotX(delThese) = [];
dotY(delThese) = [];
dotVelX(delThese) = [];
dotVelY(delThese) = [];
dotID(delThese) = [];

%% define retinal grid



[retGridxx, retGridyy] = meshgrid(linspace(-gridRadius,gridRadius,retGridRes));


%% define interpolant functions for X and Y velocities. the X and Y velocities at query point Xq, Yq is [F_x(Xq, Yq) F_y(Xq, Yq)]
if ~isempty(dotVelX)
    
    [dx_F] = scatteredInterpolant(dotX, dotY, dotVelX,'natural','none'); %linear interpolation, no extrapolation
    [dy_F] = scatteredInterpolant(dotX, dotY, dotVelY,'natural','none');
    [dotID_F] = scatteredInterpolant(dotX, dotY, dotID,'nearest','none'); %interpolant to determine dot ID, lol
    
else
    dx_F = @(x,y) nan(size(x)); %returns nans for any input :D
    dy_F = @(x,y) nan(size(x));
    dotID_F = @(x,y) nan(size(x));    
end



%% determine flow vel in X direction for each point on retgrid using interpolant functions!
retGridVelxx = dx_F(retGridxx, retGridyy);
retGridVelyy = dy_F(retGridxx, retGridyy);

retGridIDs = dotID_F(retGridxx, retGridyy);



if debug
    

    
    figure(49857);clf
    
    c1 = viscircles([0 0],max(gridRadius));
    hold on
    plot(retGridxx(1:10:retGridRes,1:10:retGridRes), retGridyy(1:10:retGridRes,1:10:retGridRes),'.','Color',[.5 .5 .5])
%     q1 = quiver(retGridxx(1:10:retGridRes,1:10:retGridRes), retGridyy(1:10:retGridRes,1:10:retGridRes), retGridVelxx, retGridVelyy);
    q1 = quiver(retGridxx, retGridyy, retGridVelxx, retGridVelyy);
    q1.Marker = '.';
    q1.MarkerEdgeColor = 'k';
    q1.Color = 'k';
    q1.ShowArrowHead = 'off';
    
    
    plot(dotX, dotY,'b.')
    q2 = quiver(dotX, dotY, dotVelX, dotVelY,'b','ShowArrowHead','off','Color','r','Marker','.','MarkerEdgeColor','b');
    q2.LineWidth = 2;
    
    m = gridRadius;
    plot([-m m], [0 0],'r-');
    plot( [0 0], [-m m],'r-');
    viscircles([0,0],m);
    
    xlim([-gridRadius gridRadius])
    ylim([-gridRadius gridRadius])
    
    axis equal
    
    if spotcheck
        dbstack
        keyboard
    end
    
end
