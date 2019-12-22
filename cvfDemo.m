% cost volume filter demo
clear; close all; clc;

% parameters
r = 9;          % window size = 2 * r + 1
eps = 0.0001;
tao1 = 7 / 255;
tao2 = 2 / 255;
alpha = 0.89;
threshBorder = 3 / 255;

% read image
Il_ori = imread('data/left.ppm');
Ir_ori = imread('data/right.ppm');
maxDisp = 16;

% normalize
Il = double(Il_ori) / 255;
Ir = double(Ir_ori) / 255;

% conver to grayscale & compute gradient in X-direction
% gradient: 0.5 * (I(i+1) - I(i-1)), thus between -0.5~0.5
% we need the gradient image in range [0, 1]
Gl = gradient(rgb2gray(Il)) + 0.5;
Gr = gradient(rgb2gray(Ir)) + 0.5;

% construct disparity cost volume
% C_i_l = (1-alpha) * min{|Ir(i-l)-Il(i)|, tao1} 
%       + alpha * min{|Gr(i-l) - Gl(i)|, tao2}
[H, W, C] = size(Il);
dispVol = zeros(H, W, maxDisp);

for d = 1 : maxDisp
    
    % truncated SAD of color image for current disp
    tmp = zeros(H,W,C) + threshBorder;
    tmp(:, (d+1):W, :) = Ir(:, 1:(W-d), :);
    p_color = mean(abs(Il - tmp), 3);
    p_color = min(p_color, tao1);
    
    % truncated SAD of gradient image for current disp
    tmp = zeros(H,W) + threshBorder;
    tmp(:, (d+1):W) = Gr(:, 1:(W-d));
    p_grad = abs(Gl - tmp);
    p_grad = min(p_grad, tao2);
    
    dispVol(:,:,d) = (1-alpha) * p_color + alpha * p_grad;
    
end

% smooth using guided filter
for d = 1 : maxDisp
    q = guidedFilterColor(Il, dispVol(:,:,d), r, eps);
    dispVol(:,:,d) = q;
end

% winner takes all
[~, diss] = min(dispVol, [], 3);

% show
figure('Name', 'disparity for left image'), imagesc(diss);


%%%%%%%%%%%%%%%% post processing %%%%%%%%%%%%%%%%%%%%%%%%

%% right-left consistency check
dispVolR = zeros(H, W, maxDisp);

for d = 1 : maxDisp
    % truncated SAD of color image for current disp
    tmp = zeros(H, W, 3) + threshBorder;
    tmp(:, 1:(W-d), :) = Il(:, (1+d):W, :);
    p_color = mean(abs(tmp - Ir), 3);
    p_color = min(p_color, tao1);
    
    % truncated SAD of gradient image for current disp
    tmp = zeros(H,W) + threshBorder;
    tmp(:, 1:(W-d)) = Gl(:, (1+d):W);
    p_grad = abs(tmp - Gr);
    p_grad = min(p_grad, tao2);
    
    dispVolR(:,:,d) = (1-alpha) * p_color + alpha * p_grad;
end

% smooth using guided filter
for d = 1 : maxDisp
    q = guidedFilterColor(Ir, dispVolR(:,:,d), r, eps);
    dispVolR(:,:,d) = q;
end

% winner takes all
[~, dissR] = min(dispVolR, [], 3);

% show
figure('Name', 'disparity for right image'), imagesc(dissR)


% detect occlusion
for i = 1 : H
    for j = 1 : W
        jr = round(j - diss(i,j));
        if jr >= 1 && jr <= W && abs(dissR(i, jr) - diss(i,j)) > 1
            diss(i,j) = -1;
        end
    end
end
figure('Name', 'occlusion detected disparity'), imagesc(diss)
            
%% fix occlusions
occPix = diss < 0;

% fill from left
diss_filled = diss;
fillVals = ones(1,H)' * maxDisp;
for col = 1 : W
    cc = occPix(:, col);
    diss_filled(cc, col) = fillVals(cc);
    fillVals = diss_filled(:, col);
end

% fill from right
diss_filled2 = diss;
fillVals = ones(1, H)' * maxDisp;
for col = W : -1 : 1
    cc = occPix(:, col);
    diss_filled2(cc, col) = fillVals(cc);
    fillVals = diss_filled2(:, col);
end

% choose minimum disparity
diss = min(diss_filled, diss_filled2);
figure('Name', 'occlusion filled disparity'), imagesc(diss)

%% filter the occluded pixels
diss_filtered = guidedFilterColor(Il, diss, 15, eps);
diss(occPix) = diss_filtered(occPix);
figure('Name', 'filal disparity'), imagesc(diss)




