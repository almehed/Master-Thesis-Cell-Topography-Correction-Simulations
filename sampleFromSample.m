function [sample2] = sampleFromSample(p,samp,sigma,imageSize)

msize = length(samp); %number of samples
idx = randperm(msize)'; %randomize indices
indi = samp(idx(1:msize*p),:); %choose correct amount of samples according to collocation prop.

sample2 = abs(normrnd(indi,sigma)); %shift by normal distribution, no negative coordinates
sample2(sample2 > max(imageSize)) = max(imageSize); %set max index to image size
sample2(sample2 < 1) = 1; %set lowest index to 1

end

