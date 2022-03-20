clear all; close all; clc;

% Load reference image
ref_img=imread('C:\Users\Suheda\Desktop\70.jpg');
ref_img_gray=rgb2gray(ref_img);

ref_pts=detectSURFFeatures(ref_img_gray);
[ref_features,ref_validPts]=extractFeatures(ref_img_gray,ref_pts);
figure,
imshow(ref_img);
hold on; plot(ref_pts.selectStrongest(50));

% Load comparison image 
image = imread('C:\Users\Suheda\Desktop\yol.png'); 
I = rgb2gray(image);

I_pts = detectSURFFeatures(I);
[I_features, I_validPts] = extractFeatures(I, I_pts);
figure,
imshow(image);
hold on; plot(I_pts.selectStrongest(50));

% Compare 
index_pairs = matchFeatures(ref_features, I_features);
%Görüntülerin karşılık geldiği noktalar bulunur
ref_matched_pts = ref_validPts(index_pairs(:,1)).Location;
I_matched_pts = I_validPts(index_pairs(:,2)).Location;

figure, showMatchedFeatures(image, ref_img, I_matched_pts, ref_matched_pts, 'montage');
title('Showing all matches');


