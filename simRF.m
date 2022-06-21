function [R,S,randomField,I1,I2,I,Inorm] = simRF(nbrOfSamples,part,sigma,imageSize, FWHM,t,randomField)
%Function to simulate two random samples on a surface.
%   randomField - a randomField to sample from, shiftet so min(randomField) = 0
%   nbrOfSamples - number of samples
%   part - part of sample 1 to be basis for sample 2
%   sigma - standarddiviation for normal shift in resample
%   imageSize - size of image [rows, cols]
%   FWHM - full width half maximum of PSF
%   t - threshold
%   randomField - Field to use as intensity, if empty a gaussian random
%   field is generated

%if no background is given, create random field
if nargin < 8 || isempty(randomField)
randomField = randomFieldSim(imageSize);
end

%sample from field
[s1,s2] = randomFieldSample(randomField,nbrOfSamples,part,sigma,imageSize);

sigma2 = FWHM/(2*(sqrt(2*log(2)))^2); %FWHM -> variance
%make into image by adding a psf of size sz
img1 = imgaussfilt(s1,sigma2);
img2 = imgaussfilt(s2,sigma2);

%thresholding
img1(img1 < t) = 0;
img2(img2 < t) = 0;

%save img1, img2...
img1c = img1;
img2c = img2;

%set zeros to NaN to leave out of correlation calculations
img1c(img1c == 0) = NaN;
img2c(img2c == 0) = NaN;

I = [img1c(:), img2c(:)]; %save for scatterplots

%compute correlations
R = zeros(1,6);
S = zeros(1,6);
[r1,~] = tiedrank(img1c(:));%compute ranks for spearman
[r2,~] = tiedrank(img2c(:));
[r3,~] = tiedrank(randomField(:));
[R1,~] = corrcoef(img1c(:),img2c(:),'Rows','complete'); %Pearson
[Rb11,~] = corrcoef(img1c(:),randomField(:),'Rows','complete'); %Pearson
[Rb12,~] = corrcoef(randomField(:),img2c(:),'Rows','complete'); %Pearson

[S1,~] = corrcoef(r1,r2,'Rows','complete'); %Spearman
[Sb11,~] = corrcoef(r1,r3,'Rows','complete'); %Spearman
[Sb12,~] = corrcoef(r3,r2,'Rows','complete'); %Spearman

%save computed correlations
R(1) = R1(1,2);
R(2) = Rb11(1,2);
R(3) = Rb12(1,2);
S(1) = S1(1,2);
S(2) = Sb11(1,2);
S(3) = Sb12(1,2);

%merge channels to one image
I1 = imfuse(img1,img2,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);

%Background normalization
img1Norm = bgNorm2(randomField,img1);
img2Norm = bgNorm2(randomField,img2);

%thresholding
img1Norm(img1Norm < t) = 0;
img2Norm(img2Norm < t) = 0;

%set zeros to NaN to leave out of correlation calculations
img1Norm(img1Norm == 0) = NaN;
img2Norm(img2Norm == 0) = NaN;

Inorm = [img1Norm(:),img2Norm(:)]; %save for scatterplots

%compute correlations
[r4,~] = tiedrank(img1Norm(:));%compute ranks for spearman
[r5,~] = tiedrank(img2Norm(:));
[R2,~] = corrcoef(img1Norm(:),img2Norm(:),'Rows','complete'); %Pearson
[Rb21,~] = corrcoef(img1Norm(:),randomField(:),'Rows','complete'); %Pearson
[Rb22,~] = corrcoef(randomField(:),img2Norm(:),'Rows','complete'); %Pearson
[S2,~] = corrcoef(r4,r5,'Rows','complete'); %Spearman
[Sb21,~] = corrcoef(r4,r3,'Rows','complete'); %Spearman
[Sb22,~] = corrcoef(r3,r5,'Rows','complete'); %Spearman

%save computed correlations
R(4) = R2(1,2);
R(5) = Rb21(1,2);
R(6) = Rb22(1,2);
S(4) = S2(1,2);
S(5) = Sb21(1,2);
S(6) = Sb22(1,2);

%Merge channels to one image
I2 = imfuse(img1Norm,img2Norm,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
end

