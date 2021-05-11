% Author : Rishabh Srivastava

% Corner detection

%% Reading input
Im = imread('IITG.jpg');
Im1 = rgb2gray(Im); % Image is in RGB format, so it is converted to grayscale
[N,M] = size(Im1); % Image size

%% Padding images to handle calculation of final pixel intensities of boundary pixels 
Im1_pad = padarray(Im1,[1 1],'replicate'); % Replicates border pixel intensities 1 time
Im1_pad = double(Im1_pad); % Converting image from uint8 to double for required calculations

%% 3x3 Gaussian Window
mask1 = zeros(3,3);
for i = 1:3
    for j = 1:3
        mask1(i,j) = exp(-0.5 * ((i-2)^2+(j-2)^2)/2);
    end
end

%% Calculating gradients
sobelX = [-1 0 1;-2 0 2;-1 0 1]; % Sobel filter for calculating gradient in horizontal direction
sobelY = sobelX'; % Sobel filter for calculating gradient in vertical direction

Ix = zeros(N,M); % Will store gradients in horizontal direction
Iy = zeros(N,M); % Will store gradients in vertical direction

for i = 1:N
    for j = 1:M
        Ix(i,j) = calc_grad(i+1,j+1,sobelX,Im1_pad); % Horizontal gradient at (i,j)
        Iy(i,j) = calc_grad(i+1,j+1,sobelY,Im1_pad); % Vertical gradient at (i,j)
    end
end

Ix = padarray(Ix,[1 1],'replicate'); % Replicates border pixel intensities 1 time
Iy = padarray(Iy,[1 1],'replicate'); % Replicates border pixel intensities 1 time

%% Calculate corner response
corner_response = zeros(N,M); % Will store corner responses for all pixels
for i = 1:N
    for j = 1:M
        corner_response(i,j) = calcR(i+1,j+1,Ix,Iy,mask1); % Corner response at (i,j)
    end
end

%% Non-maximal suppression
corner_response_pad = padarray(corner_response, [1 1], 'replicate'); % Replicates border pixel intensities 1 time
modified_response = zeros(N,M); %% Will store modified pixel intensity after non-maximal suppression
for i = 1:N
    for j = 1:M
        modified_response(i,j) = non_max_sup(i+1,j+1,corner_response_pad); %% Modified pixel intensity after non-maximal suppression
    end 
end

%% Thresholding corner responses to detect corners
threshold = max(max(corner_response));
threshold = threshold/70; %% Chosen threshold
mask_threshold = (modified_response>=threshold); %% Will store 1 if pixel corner response greater than or equal to threshold, else 0

%% Plotting figures
figure();
imshow(uint8(Im)); %% Shows original image
hold on;
for i=1:N
    for j=1:M
        if(mask_threshold(i,j)==1)
            plot(j,i,'r+','Linewidth',0.9); %% Superimposes '+' on image to indicate corners
        end
    end
end
title('\fontsize{16}Image with corners marked');

%% Function to calculate gradients
function finalvalue = calc_grad(pi, pj, mask, Im1)
%CALC_GRAD calculates gradient value for pixel (pi,pj) by applying the given mask
A = Im1(pi-1:pi+1,pj-1:pj+1);
A = A.*mask;
finalvalue = sum(A,'all');
end

%% Function to calculate corner response
function finalvalue = calcR(pi, pj, Ix, Iy, mask)
 A = Ix(pi-1:pi+1,pj-1:pj+1);
 B = Iy(pi-1:pi+1,pj-1:pj+1);
 C = A.*B; % Ix * Iy
 A = A.^2; % Ix ^2
 B = B.^2; % Iy ^2
 w1 = sum(mask.*A,'all');
 w3 = sum(mask.*B,'all');
 w2 = sum(mask.*C,'all');
 M = [w1 w2;w2 w3];
 alpha = 0.05;
 finalvalue = det(M) - alpha*(w1+w3)^2; 
end
 
%% Function to perform non-maximal suppression in 3x3 neighbourhood
function finalvalue = non_max_sup(pi, pj, Im)
A = Im(pi-1:pi+1,pj-1:pj+1);
max1 = max(max(A));
finalvalue=0;
if(max1==A(2,2))
    finalvalue=max1;
end
end
