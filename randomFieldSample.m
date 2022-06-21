function [s1,s2] = randomFieldSample(randomField,nbrOfSamples,part,sigma,imageSize)
%Create two samples of size nbrOfSamples with intensity according to
%randomField
%   randomField - a randomField to sample from, shiftet so min(randomField) = 0
%   nbrOfSamples - number of samples
%   part - part of sample 1 to be basis for sample 2
%   sigma - standarddiviation for normal shift in resample
%   imageSize - size of image [rows, cols]

%Sample according to intensity exp(randomField)
randomField(mat2gray(randomField) == 0) = 0.001; %no zeros
lambda = mat2gray(exp(randomField)); %normalize to [0,1] to use as probability

%sample binomially with probability lambda
s1 = binornd(1,lambda); 
s2 = binornd(1,lambda);

%Englarge sample until correct number of samples
while nbrOfSamples > size(find(s1),1) || nbrOfSamples > size(find(s2),1)
    s11 = binornd(1,lambda);
    s1 = s1+s11;
    s1(s1 >1) = 1;
    
    s21 = binornd(1,lambda);
    s2 = s2+s21;
    s2(s2 >1) = 1;
end
%Remove samples from s1 for correct sample size
ind = find(s1); %locations of samples
msize = numel(find(s1)); %number of samples
idx = randperm(msize); %randomize samples
ind = ind(idx(1:msize-nbrOfSamples)); %only keep nbrOfSamples samples
s1(ind') = 0;

%Resample part number of samples from s1, shift by normal distribution with
%variance sigma^2
[i,j] = find(s1);
s2_coord = sampleFromSample(part,[i,j],sigma,imageSize);
%round sample to indexes
s2_round = zeros(size(s2_coord));
s2_round(:,1) = round(s2_coord(:,1));
s2_round(:,2) = round(s2_coord(:,2));

s22 = sub2ind(imageSize, s2_round(:,1), s2_round(:,2));

%create rest of sample2, size = (1-part)*nbrOfSamples
ind = find(s2);
msize = numel(find(s2));
idx = randperm(msize);
ind = ind(idx(1:msize-nbrOfSamples*(1-part)));
s2(ind') = 0;

s2(s22) = 1; %add resampled values


end

