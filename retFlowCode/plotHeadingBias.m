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
        vel_pt = [vel_pt; [heading, median(vel.sing_bias)]];
        accel_pt = [accel_pt; [heading, median(accel.sing_bias)]];
    end
    figure;
    plot(accel_pt(:,1), accel_pt(:, 2), 'ro');
    xlabel('Ground truth heading ({\circ})');
    ylabel('Bias ({\circ})');
    title(['Acceleration singularity bias over ground truth']);
    figure;
    plot(vel_pt(:,1), vel_pt(:, 2), 'bo');
    xlabel('Ground truth heading ({\circ})');
    ylabel('Bias ({\circ})');
    title(['Velocity singularity bias over ground truth']);

    
    