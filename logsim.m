%% Create background
%create 100x100 coordinates
x = linspace(0.01,1,100);
y = x;
[X,Y] = meshgrid(x,y);

Z = 1-(log(X+1)+log(Y+1)); %non-uniform background

Z2 = 0.5*ones(100); %uniform background

figure(2)
clf
imshow(Z)

%% For scattered images
%parameters
nbrOfSamples = 1000; %number of samples
sigma = 0.1; %std of normaldispacement when resampling
imageSize = [100,100]; %image size 
FWHM = 4; %fullwidth half maximum of PSF
t = 1e-4; %threshold
part = 0; %collocation proportion 

%simulate two random samples on a surface according to Z. If no Z, a
%Gaussian random field according to randomFieldSim(imageSize) is generated
[R,S,randomField,I1,I2,I,Inorm] = simRF(nbrOfSamples,part,sigma,imageSize,FWHM,t,Z);

%Display calculated correlation coefficients
disp('pearson(rg,bg,br)/spearman(rg,bg,br)')
disp([R(1:3),S(1:3)])
disp('pearson(rg,bg,br)/spearman(rg,bg,br) after norm')
disp([R(4:6),S(4:6)])

%Show simulated images and background
figure(3)
clf
imshow(I1)
figure(4)
clf
imshow(I2)
figure(5)
clf
imshow(mat2gray(randomField))

%scatterplots of 
figure(6)
clf
scatter(I(:,1),I(:,2),1,'k','MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3)
title('Orgiginal simulation intensities')
figure(7)
clf
scatter(Inorm(:,1),Inorm(:,2),1,'k','MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3)
title('Normalized simulation intensities')

%% Many times for all sigma and all parts
%parameters
nbrOfSamples = [100, 250, 500, 750, 1000, 1500, 2000,2500, 3000, 3500, 4000, 4500, 5000];
sigma = [0.01 0.025 0.05 0.075 0.1 0.25 0.5 0.75 1 2.5 5 7.5 10];
imageSize = [100,100];
FWHM = 4;
t = 1e-4;
part = linspace(0,1,11); 
n = 100; %numbers of simulations for each parameter combination

%initialize arrays for saving mean correlation coefficients
R1 = zeros(length(part),length(sigma), length(nbrOfSamples));
R2 = zeros(length(part),length(sigma), length(nbrOfSamples));
S1 = zeros(length(part),length(sigma), length(nbrOfSamples));
S2 = zeros(length(part),length(sigma), length(nbrOfSamples));

%initialize arrays for saving stds of correlation coefficients
stdR1 = zeros(length(part),length(sigma), length(nbrOfSamples));
stdR2 = zeros(length(part),length(sigma), length(nbrOfSamples));
stdS1 = zeros(length(part),length(sigma), length(nbrOfSamples));
stdS2 = zeros(length(part),length(sigma), length(nbrOfSamples));

%initialize arrays for saving correlation coefficients
r1 = zeros(1,n);
r2 = zeros(1,n);
s1 = zeros(1,n);
s2 = zeros(1,n);

flag = 0; 
tic %compute runtime
for l = 1:length(nbrOfSamples)
    for k = 1:length(sigma)
        for j = 1:length(part)
            for i = 1:n
                [R,S,randomField,I1,I2] = simRF(nbrOfSamples(l),part(j),sigma(k),imageSize, FWHM,t,randomField);
                if isnan(R(1))|| isnan(R(4)) || isnan(S(1)) || isnan(S(4)) %break if NaN
                    flag = 1;
                    break
                end
                %save correlation coefficient
                r1(i) = R(1);
                r2(i) = R(4);
                s1(i) = S(1);
                s2(i) = S(4);
            end
            if flag
                break
            end
                %calculate and save mean correlation coefficient for each parameter comb.
                R1(j,k,l) = mean(r1);
                R2(j,k,l) = mean(r2);
                S1(j,k,l) = mean(s1);
                S2(j,k,l) = mean(s2);
                
                %calculate and save std of correlation coefficients for
                %each parameter comb.
                stdR1(j,k,l) = std(r1);
                stdR2(j,k,l) = std(r2);
                stdS1(j,k,l) = std(s1);
                stdS2(j,k,l) = std(s2);
        end
        %break if nan
        if flag
                break
        end
        %display iteration
        disp('sigma iter')
        disp(k)
    end
    %break if nan
    if flag
        break
    end
    %display iteration
    disp('nbrOfSamples iter')
    disp(l)
end
toc

ts = abs((R1-R2)./(sqrt(stdR1.^2/n+stdR2.^2/n))); %t-test of difference for all parameter combs.
not_reject = find(ts <= 2.36); %find all combination not passing test
%% For T-tests
%parameters
nbrOfSamples = 1000;
sigma = 0.1;
imageSize = [100,100];
FWHM = 4;
t = 1e-4;
part = 0;
n=100;

%initialize arrays to save correlation coefficients
Rs = zeros(n,6);
Ss = zeros(n,6);

for i = 1:n
[R,S,randomField,I1,I2] = simRF(nbrOfSamples,part,sigma,imageSize, FWHM,t,randomField);
%save correlation coefficients
Rs(i,:) = R;
Ss(i,:) = S;
end

%t-tests
[hr,pr,cir] = ttest(Rs(:,1),Rs(:,4)) %pearson, between samples
[hs,ps,cis] = ttest(Ss(:,1),Ss(:,4)) %Spearman, between samples
[hrgb,prgb,cirgb] = ttest(Rs(:,2),Rs(:,5)) %Pearson, between sample1 and background
[hsgb,psgb,cisgb] = ttest(Ss(:,2),Ss(:,5)) %Spearman, between sample1 and background
[hrrb,prrb,cirrb] = ttest(Rs(:,3),Rs(:,6)) %Pearson, between sample2 and background
[hsrb,psrb,cisrb] = ttest(Ss(:,3),Ss(:,6)) %Spearman, between sample2 and background

%display mean correlation coeffs
mr = mean(Rs)
ms = mean(Ss)

%display stds of correlation coeffs
stdr = std(Rs)
stds = std(Ss) 

%% Some plots
R1s =squeeze(R1);
R2s = squeeze(R2);
stdR1 =squeeze(stdR1);
stdR2 = squeeze(stdR2);
figure
errorbar(part,R1(:,1),stdR1(:,1),'-ob')
hold on
errorbar(part,R1(:,2),stdR1(:,2),'-or')
errorbar(part,R1(:,3),stdR1(:,3),'-og')
errorbar(part,R1(:,4),stdR1(:,4),'-om')
errorbar(part,R1(:,5),stdR1(:,5),'-oy')
errorbar(part,R1(:,6),stdR1(:,6),'-ok')
errorbar(part,R1(:,7),stdR1(:,7),'-oc')
errorbar(part,R2(:,7), stdR2(:,7),'-.ob')
errorbar(part,R2(:,2),stdR2(:,2),'-.or')
errorbar(part,R2(:,3),stdR2(:,3),'-.og')
errorbar(part,R2(:,4),stdR2(:,4),'-.om')
errorbar(part,R2(:,5),stdR2(:,5),'-.oy')
errorbar(part,R2(:,6),stdR2(:,6),'-.ok')
errorbar(part,R2(:,7),stdR2(:,7),'-.oc')
grid on
xlabel('collocation proportion')
ylabel('correlation coefficient')
title('Pearson')
legend('Location','southeast')

%% More plots
R1s =squeeze(R1);
R2s = squeeze(R2);
S1s =squeeze(S1);
S2s = squeeze(S2);
stdR1s = squeeze(stdR1);
stdR2s = squeeze(stdR2);
figure
errorbar(part,R2s(:,1),stdR2s(:,1),'-o')
hold on
errorbar(part,R2s(:,2),stdR2s(:,2),'-o')
errorbar(part,R2s(:,3),stdR2s(:,3),'-o')
errorbar(part,R2s(:,4),stdR2s(:,4),'-o')
errorbar(part,R2s(:,5),stdR2s(:,5),'-o')
errorbar(part,R2s(:,6),stdR2s(:,6),'-o')
errorbar(part,R2s(:,7),stdR2s(:,7),'-o')
errorbar(part,R2s(:,8),stdR2s(:,8),'-o')
errorbar(part,R2s(:,9),stdR2s(:,9),'-o')
errorbar(part,R2s(:,10),stdR2s(:,10),'-o')
errorbar(part,R2s(:,11),stdR2s(:,11),'-o')
errorbar(part,R2s(:,12),stdR2s(:,12),'-o')
errorbar(part,R2s(:,13),stdR2s(:,13),'-o')

grid on
xlabel('collocation proportion')
ylabel('correlation coefficient')
title('Pearson correlation, 1000 points')
legend('Location','southeast')

