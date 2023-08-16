clear;close all;clc;

addpath(genpath('whyte_code'));
addpath(genpath('cho_code'));

%% parameters
opts.prescale = 1; 
opts.xk_iter = 5; 
opts.k_thresh = 20;
Blambda = 0.2;
lambda = 0.1;
filename = 'images\flower.jpg'; 
opts.kernel_size = 35; 
saturation = 0;
lambda_grad = 4e-3;      
opts.gamma_correct = 1.0; 

%% Note: lambda_tv, lambda_l0, weight_ring are non-necessary, they are not used in kernel estimation.
lambda_tv = 0.001; lambda_l0 = 1e-3; weight_ring = 1;
 
%% read image
y = imread(filename);
    
if size(y,3)==3
    yg = im2double(rgb2gray(y));
else
    yg = im2double(y);
end

tic;
[kernel, interim_latent] = blind_deconv(yg, lambda, Blambda, lambda_grad, opts);
toc

%% ============Non-blind deconvolution====================%%
y = im2double(y);
%% Final Deblur: 

if ~saturation
%% 1. TV-L2 denoising method
   Latent = ringing_artifacts_removal(y, kernel, lambda_tv, lambda_l0, weight_ring);
else
   %% 2. Whyte's deconvolution method (For saturated images)
   Latent = whyte_deconv(y, kernel);
end


figure; imshow(Latent)


k = kernel - min(kernel(:));
k = k./max(k(:));

   
imwrite(k,[filename(8:end-4), '_kernel.png']);
imwrite(Latent,[filename(8:end-4), '_deblurred.png']);

