clear
clc
close all;
load('Series_Sig.mat');
%  deleteline=[1:2:79];      
%   Series1_Sig=Series1_Sig(:,deleteline)  ; 
  Imagecmplx = fftshift(fft2(Series1_Sig)); 
  figure
 
imshow(abs(Imagecmplx)',[]); 

 Series1_Sig(1:60,:)=Series1_Sig(1:60,:)*0;
%  Series1_Sig(:,[1:35 37:80])=Series1_Sig(:,[1:35 37:80])*0;  
  Imagecmplx = fftshift(fft2(Series1_Sig)); 
  figure
 
  imshow(abs(Imagecmplx)',[]); 