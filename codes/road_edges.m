clc;clear all;close all;
%Part1: CANNY EDGE DETECTION
img1 = imread ('C:\Users\Suheda\Desktop\1.jpg');
figure, imshow(img1); title('Original');
img= rgb2gray(img1);
img=im2double(img);
%Gauss Filtresi Katsayısı
B = [2, 4, 5, 4, 2; 4, 9, 12, 9, 4;5, 12, 15, 12, 5;4, 9, 12, 9, 4;2, 4, 5, 4, 2 ];
B = 1/159.* B;
%Görüntünün Gauss Katsayısı ile Konvolüsyonu
A=conv2(img, B, 'same');
figure, imshow(A); title('Gaussian Filter Result');
%SOBEL EDGE DETECTION
maskx = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
masky = [1, 2, 1; 0, 0, 0; -1, -2, -1];
%x in ve y nin gradyanı alındı.
Gx = conv2(A, maskx, 'same');
Gy = conv2(A, masky, 'same');
G =sqrt(Gx.^2+Gy.^2);
figure, imshow(G);
title('Sobel Edge Detection Result');
%Kenar inceltme,Non-maximal suppression
%Yönlendirmelerin hesaplabması, teta= gradyan türevinin yönü,Kenar yönü
teta = atan2(Gy,Gx);
teta = teta*180/pi;
row=size(A,1);
col=size(A,2); 
%Negatif yönler için ayar, bütün yönlerin pozitif yapılması
for i=1:row
    for j=1:col
        if (teta(i,j)<0) 
            teta(i,j)=360+teta(i,j);
        end
    end
end
 teta2=zeros(row,col);
 %Yönleri en yakın 0, 45, 90 veya 135 dereceye ayarlama
for i=1:row
    for j=1:col
        if (teta(i,j)>=0 && teta(i,j)<22.5) || (teta(i,j)>=157.5 && teta(i,j)<202.5) || (teta(i,j)>=337.5 && teta(i,j)<=360)
            teta2(i,j)=0;
        elseif (teta(i,j)>=22.5 && teta(i,j)<67.5) || (teta(i,j)>=202.5 && teta(i,j)<247.5)
            teta2(i,j)=45;
        elseif (teta(i,j)>=67.5 && teta(i,j)<112.5) || (teta(i,j)>=247.5 && teta(i,j)<292.5)
            teta2(i,j)=90;
        elseif (teta(i,j)>=112.5 && teta(i,j)<157.5) || (teta(i,j)>=292.5 && teta(i,j)<337.5)
            teta2(i,j)=135;
        end
    end
end
% figure, imagesc(teta2); colorbar;

M=zeros(row,col);
for i=2:row-1
    for j=2:col-1
        if (teta2(i,j)==0)
            M(i,j) = G(i,j) == max([G(i,j),G(i,j+1),G(i,j-1)]);
        elseif (teta2(i,j)==45)
            M(i,j) = G(i,j) == max([G(i,j),G(i+1,j-1),G(i-1,j+1)]);
        elseif (teta2(i,j)==90)
            M(i,j) = G(i,j) == max([G(i,j),G(i+1,j),G(i-1,j)]);
        elseif (teta2(i,j)==135)
            M(i,j) = G(i,j) == max([G(i,j),G(i+1,j+1),G(i-1,j-1)]);
        end
    end
end
M = M.*G;
figure, imshow(M);
title('After Non-maximal Suppression');
%İkili eşik yöntemi, Hysteresis
%Eşik Değeri belirledik
T_Low = 0.075;
T_High = 0.175;
T_Low = T_Low * max(max(M)); 
T_High = T_High * max(max(M));
T_res = zeros(row,col);
for i = 1  : row
    for j = 1 : col
        if (M(i,j)<T_Low) %Piksel değeri, düşük eşikten küçükse sıfıra çekilir
            T_res(i,j)=0;
        elseif (M(i,j)>T_High)%yüksek eşiğin üzerindeyse 1 
            T_res(i,j)=1;
%İki eşik arasındaysa, bu pikselden herhangi bir yol yoksa sıfır olarak ayarlanır.
        elseif ( M(i+1,j)>T_High || M(i-1,j)>T_High || M(i,j+1)>T_High || M(i,j-1)>T_High || M(i-1, j-1)>T_High || M(i-1, j+1)>T_High || M(i+1, j+1)>T_High || M(i+1, j-1)>T_High)
            T_res(i,j)=0;
        end
    end
end
figure, imshow(T_res);
title('Canny Edge Detection Result');
%------------------------------------------------------------------------------
%Part2: Hough Method
[H,Theta,Rho] = hough(T_res);      
P  = houghpeaks(H,20,'threshold',ceil(0.1*max(H(:))));         
lines = houghlines(T_res,Theta,Rho,P,'FillGap',12,'MinLength',40); 
 figure; imshow(img1); hold on; title('Hough Method Result');

 [w,h]=size(lines);
 for i=1:1:h
     xy=[lines(i).point1;lines(i).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',1.5,'Color','green');
 end

   
   
    


 