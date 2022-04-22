% Title:            estAccel.m
%
% Authors:          Mengjian Hua and Charlie Burlingham
%
% Purpose:          Loads a video, computes optic flow using the method of
%                   Farneback (2000) and subsequently estimates optic acceleration using
%                   Farid & Simoncelli (FS) spatiotemporal derivative filters.
%
% Usage:            You must first compile MatlabPyrTools by running compilePyrTools.m
%                   in the subdirectory of MEX files called MatlabPyrTools.
%                   Then input the video name/path and number of frames for
%                   which you want to compute optic acceleration beyond
%                   frame 1. Use at least 30 frames or so, as the flow
%                   estimator takes 11 consecutive frames by default,
%                   and the derivative filters require 5 consecutive flow fields.
%
% Inputs:           NumberofFrames: the number of frames that we are going to use,
%                                   starting from frame 1 of video.
%                   vidName: video file name
%
% Last updated:     April 22 2022
%
% License:          Copyright (C) 2022 Mengjian Hua and Charlie Burlingham
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

function [accel_x,accel_y] = estAccel(vidName,NumberofFrames)

%% load videos from frame 1 to n (NumberofFrames)
vidObj = VideoReader(vidName);

count = 0;
while hasFrame(vidObj) && count<NumberofFrames
    frameRGB = readFrame(vidObj);
    frameGray = double(rgb2gray(frameRGB));
    count = count + 1;
    images(:,:,count) = frameGray;
end

%% Pre-process the video

% Blur and downsample (decimate) the video first to reduce temporal
% aliasing and enable estimation of larger optic flow

numPyrLevels = 1; % number of times to downsample

reduced_images{1} = images;
for kk = 2:1+numPyrLevels
    for frames = 1:NumberofFrames
        reduced_images{kk}(:,:,frames) = impyramid(reduced_images{kk-1}(:,:,frames),'reduce');
    end
end

imageSequence = reduced_images{end};

%% Estimate optical flow using Farneback (2000)

% default parameters
tic % time algorithm
opt.basisLen   = 11;
opt.basisSigma  = 1.6;
opt.gamma       = 1/256;
opt.tensorLen   = 41;
opt.tensorSigma = 6.5;

% need m video frames (m = opt.basisLen) to estimate one optic flow field

try % if you have parallel computing toolbox, use parfor to speed up compute time
    parfor i = 1:NumberofFrames/opt.basisLen
        disp(i);
        [flow2x(:,:,i), flow2y(:,:,i), D(:,:,i), cout(:,:,i)] = estimateOpticFlow2D(imageSequence(:,:,i:i+opt.basisLen-1));
    end
catch % if not use for loop
    for i = 1:NumberofFrames/opt.basisLen
        disp(i);
        [flow2x(:,:,i), flow2y(:,:,i), D(:,:,i), cout(:,:,i)] = estimateOpticFlow2D(imageSequence(:,:,i:i+opt.basisLen-1));
    end
end

%% Estimate optic acceleration using Farid-Simoncelli derivative filters
accel_x = dz3test(flow2x);
accel_y = dz3test(flow2y);

end
