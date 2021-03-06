% Title:            plotHeadingBias.m
%
% Authors:          Oliver Xu and Charlie Burlingham
%
% Purpose:          Loads pre-computed flow/acceleration field singularity
%                   data for a series of headings, plots the bias over ground
%                   truth heading
%
% Inputs:           foldername: name of directory containing bias data
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
function plotHeadingBias(foldername)
    vel_pt = [];
    accel_pt = [];
    cd(foldername);
    datafolders = dir('*0*');
    for ii = 1:length(datafolders)
        folder = [foldername '/' datafolders(ii).name];
        velfile = [folder '/vel_bias.mat'];
        vel = load(velfile);
        accelfile = [folder '/accel_bias.mat'];
        accel = load(accelfile);
        heading = str2num(extractBefore(datafolders(ii).name, '-'));
        accel.sing_bias(:,2)
        % Set bias to be for the singularity further away from 0
        if median(abs((accel.sing_bias(:, 2)+heading))) >= median(abs((accel.sing_bias(:, 1)+heading)))
            bias = median(accel.sing_bias(:, 2));
        else
            bias = median(accel.sing_bias(:, 1));
        end
        % Check if the two singularities are the same, if so, remove from
        % plot
        if median(accel.sing_bias(:, 2)) == median(accel.sing_bias(:, 1))
            heading = NaN;
            bias = NaN;
        end
        vel_pt = [vel_pt; [heading, bias]];
        accel_pt = [accel_pt; [heading, bias]];
    end
    figure;
    plot(accel_pt(:,1), -accel_pt(:, 2), 'ro'); % invert the y-axis
    xlabel('Ground truth heading ({\circ})');
    ylabel('Bias ({\circ})');
    title(['Acceleration singularity bias over ground truth']);
    
%     figure;
%     plot(vel_pt(:,1), vel_pt(:, 2), 'bo');
%     xlabel('Ground truth heading ({\circ})');
%     ylabel('Bias ({\circ})');
%     title(['Velocity singularity bias over ground truth']);

    
    