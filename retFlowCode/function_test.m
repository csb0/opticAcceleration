% Title:            function_test.m
%
% Authors:          Oliver Xu and Charlie Burlingham
%
% Purpose:          Uses various functions collected in this folder to
%                   generate data, calculate singularities, and generate
%                   figures. Functions used here are based on code provided
%                   by Matthis et. al, used in their paper "Retinal optic
%                   flow during natural locomotion"
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
%% set params
save_path = 'C:\Users\Oli\Documents\charlie_data\kate_data\function_test\';
save_path2 =  'C:\Users\Oli\Documents\charlie_data\kate_data\function_test2\';
heading = 30; % degrees
fix_depth = 6000; % mm

%% compute flows for various plane depths and save out data
% ii represents plane depth
for ii = -4500:500:6000
    retFlowSimFunction(heading, ii, fix_depth, save_path);
    clearvars -except save_path save_path2 ii heading fix_depth;
end
% add loop over various headings at a fixed depth (12.5m)(for the purpose of reproducing fig. 5 in
% the paper)
for h = 0:5:50 % heading varies from 0 to 50 deg in 5 deg intervals
    s_p = [save_path2 num2str(h) '-']; % heading (in degrees) dash fixation depth is the folder name
    retFlowSimFunction(h, 7.5e3, 7.5e3, s_p); % plane is 12.5m away from observer's initial position
    clearvars -except save_path save_path2 h heading fix_depth;
end
% plot heading bias over ground truth heading (second singularity) [part of
% fig. 5]

%% calculate singularity estimates from flow and acceleration fields
save_vel_and_accel_biases(save_path);
save_vel_and_accel_biases(save_path2);

%% display singularity, flow/accel data
plotFigs(save_path, fix_depth, heading); % plots fig 6
plotHeadingBias(save_path2); % plots fig 5
plotFlowAccel([save_path2 '/45-7500']); % example flow/accel field, can be adjusted