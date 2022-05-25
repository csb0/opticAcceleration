% Title:            plotFlowAccel.m
%
% Authors:          Oliver Xu and Charlie Burlingham
%
% Purpose:          Loads pre-computed flow data, calculates acceleration 
%                   data, plots flow and acceleration for middle frame
%
% Inputs:           foldername: name of directory containing flow data
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
function plotFlowAccel(foldername)
    cd(foldername);
    frame_ct = length(dir('eyeStructFrame*'));
    for i = 1:frame_ct

        filename = sprintf("eyeStructFrame%09d.mat",i);
        load(filename);

        retGridVelxx(:,:,i) = thisEyeStruct.retGridVelxx;
        retGridVelyy(:,:,i) = thisEyeStruct.retGridVelyy;
        dt = 1/thisEyeStruct.samplingRate;
        heading(1,i) = thisEyeStruct.headX;
        heading(2,i) = thisEyeStruct.headY;
    end
    retGridxx = thisEyeStruct.retGridxx;
    retGridyy = thisEyeStruct.retGridyy;
    % acceleration
    for frame = 3:(frame_ct-2)
        accel_x(:,:,frame) = retGridVelxx(:,:,frame) - retGridVelxx(:,:,frame-1);
        accel_y(:,:,frame) = retGridVelyy(:,:,frame) - retGridVelyy(:,:,frame-1);
    end
    middle_frame = idivide(frame_ct, int16(2), 'floor'); % middle frame of simulation, rounded down
    figure;
    quiver(retGridxx(1:5:end, 1:5:end), retGridyy(1:5:end, 1:5:end), ...
        retGridVelxx(1:5:end, 1:5:end, middle_frame), retGridVelyy(1:5:end, 1:5:end, middle_frame));
    xlim([-60 60]);
    ylim([-60 60]);
    xlabel('Retinal horizontal position ({\circ})');
    ylabel('Retinal vertical position ({\circ})');
    title(['Flow field at frame ' num2str(middle_frame) ', heading ' num2str(heading(1)) '{\circ}']);

    % plot accel and singularity/heading
    figure;
    quiver(retGridxx(1:5:end, 1:5:end), retGridyy(1:5:end, 1:5:end), ...
        accel_x(1:5:end, 1:5:end, middle_frame), accel_y(1:5:end, 1:5:end, middle_frame), 1);
    % hold on;
    % plot(median(estimator), 0, 'ro');
    xlim([-60 60]);
    ylim([-60 60]);
    xlabel('Retinal horizontal position ({\circ})');
    ylabel('Retinal vertical position ({\circ})');
    title(['Acceleration field at frame ' num2str(middle_frame) ', heading ' num2str(heading(1)) '{\circ}']);


