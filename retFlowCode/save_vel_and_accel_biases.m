% Title:            save_vel_and_accel_biases.m
%
% Authors:          Oliver Xu and Charlie Burlingham
%
% Purpose:          Loads pre-computed flow data, calculates acceleration
%                   data, and finds the first singularity, then a second
%                   singularity (if it exists)
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
function save_vel_and_accel_biases(foldername)
    cd(foldername);
    D = dir;
    D2 = D(arrayfun(@(x) ~strcmp(x.name(1),'.'),D));
    for ii = 1:length(D2)
        cd(foldername);
        currD = D2(ii).name;
        cd(currD);
        frame_ct = length(dir('eyeStructFrame*'));
        tic
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
        % for i = 1:frame_ct
        %     angle(i) = acosd(transaccel(i,1)/norm(transaccel(i,:)));
        % end
        % for ii = 10:150
        %     figure(4)
        %     quiver(retGridxx(:,:,ii-9),retGridyy(:,:,ii-9),retGridVelxx(:,:,ii-9),retGridVelyy(:,:,ii-9),1.5,'color','k');
        %     hold on
        %     plot(heading(1,ii-9),heading(2,ii-9),'*r');
        %     hold off
        %     pause(0.1)
        %     shg
        % end
        toc
        %% use lineXline
        % accel_x = dz3test(retGridVelxx);
        % accel_y = dz3test(retGridVelyy);
        spacing = linspace(-60,60,401);
        % acceleration
        for frame = 3:(frame_ct-2)

            % %% David's method
            accel_x(:,:,frame) = retGridVelxx(:,:,frame) - retGridVelxx(:,:,frame-1);
            accel_y(:,:,frame) = retGridVelyy(:,:,frame) - retGridVelyy(:,:,frame-1);
            [Gridres,~,~] = size(retGridyy);
            index = 1;
            for i = 1:Gridres
                for j = 1:Gridres
                    if(retGridyy(i,j)==0&& abs(retGridxx(i,j))>0)
                           list_pos(:,index) = [retGridxx(i,j);retGridyy(i,j)];
                           list_id(:,index) = [i;j];
                           list_v_x(index) = accel_x(i,j,frame);
                           index = index + 1;
                    end
                end
            end
            index_2 = 1;
            for i = 1:Gridres
                for j = 1:Gridres
                    if(retGridyy(i,j)==0)
                           list_v_x_2(index_2) = accel_x(i,j,frame);
                           index_2 = index_2 + 1;
                    end
                end
            end
            [min_value,min_index] = min(abs(list_v_x));
            estimator(frame-2) = list_pos(1,min_index);
            bias(frame-2) = -(list_pos(1,min_index) - heading(1,frame)) ;
            profile(frame-2,:) = list_v_x_2;
            mean_heading_over_time(frame-2,:) = sum(heading(1,3:frame))/(frame-3+1);
            mean_est_over_time(frame-2,:) = sum(estimator(1,1:frame-2))/(frame-3+1);
            %% find 2 singularities
            [out, idx] = sort(abs(list_v_x));
            sing_one(frame-2, :) = list_pos(1, idx(1));
            frame_two_idx = 2;
            while abs(list_pos(1, idx(1))-list_pos(1, idx(frame_two_idx))) < 10 % 10 degree tolerance for singularities
            % while abs(idx(1)-idx(frame_two_idx)) < 10
                frame_two_idx = frame_two_idx + 1;
            end
            [idx(1), idx(frame_two_idx)];
            sing_two(frame-2, :) = list_pos(1, idx(frame_two_idx));
            if out(frame_two_idx) - out(1) > 10*out(1) % was 0.5x
                sing_two(frame-2, :) = NaN;
            end
            sing(frame-2, 1) = max(sing_one(frame-2, :), sing_two(frame-2, :));
            sing(frame-2, 2) = min(sing_one(frame-2, :), sing_two(frame-2, :));
            sing_bias(frame-2, 1) = sing(frame-2, 1)-heading(1, frame);
            sing_bias(frame-2, 2) = sing(frame-2, 2)-heading(1, frame);
            end
        save(['accel_bias' '.mat'], 'sing_bias');
        save(['accel_estimatator' '.mat'], 'sing');
        clear list_pos;
        clear list_id;
        clear list_v_x
        clear index;
        clear index_2;
        clear list_v_x_2;
        clear sing;
        clear sing_one;
        clear sing_two;
        clear sing_bias;
        % velocity
        for frame = 3:(frame_ct-2)

            % %% David's method
            [Gridres,~,~] = size(retGridyy);
            index = 1;
            for i = 1:Gridres
                for j = 1:Gridres
                    if(retGridyy(i,j)==0&& abs(retGridxx(i,j))>0)
                           list_pos(:,index) = [retGridxx(i,j);retGridyy(i,j)];
                           list_id(:,index) = [i;j];
                           list_v_x(index) = retGridVelxx(i,j,frame);
                           index = index + 1;
                    end
                end
            end
            index_2 = 1;
            for i = 1:Gridres
                for j = 1:Gridres
                    if(retGridyy(i,j)==0)
                           list_v_x_2(index_2) = retGridVelxx(i,j,frame);
                           index_2 = index_2 + 1;
                    end
                end
            end
            [min_value,min_index] = min(abs(list_v_x));
            vel_estimator(frame-2) = list_pos(1,min_index);
            vel_bias(frame-2) = -(list_pos(1,min_index) - heading(1,frame)) ;
            vel_profile(frame-2,:) = list_v_x_2;
            vel_mean_heading_over_time(frame-2,:) = sum(heading(1,3:frame))/(frame-3+1);
            vel_mean_est_over_time(frame-2,:) = sum(vel_estimator(1,1:frame-2))/(frame-3+1);
            %% find 2 singularities
            [out, idx] = sort(abs(list_v_x));
            sing_one(frame-2, :) = list_pos(1, idx(1));
            frame_two_idx = 2;
            while abs(list_pos(1, idx(1))-list_pos(1, idx(frame_two_idx))) < 10 % 10 degree tolerance for singularities
            % while abs(idx(1)-idx(frame_two_idx)) < 10
                frame_two_idx = frame_two_idx + 1;
            end
            [idx(1), idx(frame_two_idx)];
            sing_two(frame-2, :) = list_pos(1, idx(frame_two_idx));
            if out(frame_two_idx) - out(1) > 10*out(1) % was 0.5x
                sing_two(frame-2, :) = NaN;
            end
            sing(frame-2, 1) = max(sing_one(frame-2, :), sing_two(frame-2, :));
            sing(frame-2, 2) = min(sing_one(frame-2, :), sing_two(frame-2, :));
            sing_bias(frame-2, 1) = sing(frame-2, 1)-heading(1, frame);
            sing_bias(frame-2, 2) = sing(frame-2, 2)-heading(1, frame);
        end

        save(['vel_bias' '.mat'], 'sing_bias');
        save(['vel_estimator' '.mat'], 'sing');
    end
end