% Title:            plotFigs.m
%
% Authors:          Oliver Xu and Charlie Burlingham
%
% Purpose:          Loads pre-computed flow/acceleration field singularity
%                   data, plots the bias relative to ground truth heading
%                   over depth discrepancy (fixation depth-plane depth)
%
% Inputs:           foldername: name of directory containing bias data
%                   fix: fixation depth (mm)
%                   heading: heading angle (deg)
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
function plotFigs(foldername, fix, heading) % inputs: directory, fixation depth
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
        displacement = fix/1000-str2num(datafolders(ii).name)/1000; % fixation depth minus plane depth
        
        vel_pt = [vel_pt; [displacement, median(vel.sing_bias)]];
        accel_pt = [accel_pt; [displacement, median(accel.sing_bias)]];
    end

    vXs = vel_pt(:, 1);
    [sorted_vx, v_i] = sort(vXs);
    vYs = vel_pt(:, 2);
    for ii=1:length(v_i)
        sorted_vy(ii, :) = vYs(v_i(ii));
    end
    vYs2 = vel_pt(:, 3);
    for ii=1:length(v_i)
        sorted_vy2(ii, :) = vYs2(v_i(ii));
    end
    xs = accel_pt(:, 1);
    [sorted_x, a_i] = sort(xs);
    ys = accel_pt(:, 2);
    for ii=1:length(a_i)
        sorted_y(ii, :) = ys(a_i(ii));
    end
    ys2 = accel_pt(:, 3);
    for ii=1:length(a_i)
        sorted_y2(ii, :) = ys2(a_i(ii));
    end

    % %% ratio of accel over vel bias
    % %ratio = ys./vYs;
    % difference = ys-vYs;
    % difference2 = ys2-vYs2;
    % % avg_ratio = mean(ratio);
    % avg_difference = mean(difference);
    % avg_difference2 = mean(difference2);
    % plot(xs, ratio, 'ro');

    %% velocity and accel plot
    % need to update how the 'fixation singularity' is assigned, maybe use
    % minimal distance to zero?
    figure;
    plot(sorted_vx, -sorted_vy, 'b+-');
    hold on;
    plot(sorted_vx, -sorted_vy2, 'k+-');
    hold on;
    plot(sorted_x, -sorted_y, 'ro-');
    hold on;
    plot(sorted_x, -sorted_y2, 'go-');
    xlabel('Displacement (m)');
    ylabel('Bias ({\circ})');
    legend('Flow field fix', 'Flow field', 'Acceleration field fix', 'Acceleration field', 'Location', 'northeast');
    title([int2str(heading) ' deg']);

%     %% acceleration
%     % linear regression on log(x) against y
%     logx = [log(sorted_x), ones(size(sorted_x))];
%     sorted_y = sorted_y(2:end, :);
%     logx = logx(2:end, :);
%     b_accel = regress(sorted_y, logx);
%     [~,~,~,~,stats] = regress(sorted_y, logx);
%     r2 = stats(1);
% 
%     % remove the column of '1's from logx to facilitate plotting
%     logx = logx(:, 1);
% 
%     % graph data points with trendline
%     plot(logx, sorted_y, 'ro');
%     hold on;
%     a = xlim;
%     x = linspace(a(1), a(2));
%     y = b_accel(1)*x+b_accel(2);
%     plot(x, y);
%     hold on;
%     grid on;
%     % Place equation in upper left of graph.
%     b = ylim;
%     xt = 0.2 * (a(2)-a(1)) + a(1)
%     yt = 0.90 * (b(2)-b(1)) + b(1)
%     caption = sprintf('y = %f * x + %f \n r^2 = %f', b_accel(1), b_accel(2), r2);
%     text(xt, yt, caption, 'FontSize', 12, 'Color', 'black', 'FontWeight', 'bold');
%     xlabel('log(displacement)');
%     ylabel('Bias (degrees)');
