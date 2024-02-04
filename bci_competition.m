
clear;
clc;
close all;

load('BCICIV_calib_ds1b_100Hz.mat');

cnt = double(cnt) *0.1;
fs =nfo.fs;
str = nfo.clab;



%%
x = nfo.xpos;
y = nfo.ypos;
channelName = nfo.clab ; 
r = 1;
theta = 0:0.1:2*pi ;
plot(x , y ,"O" ,'LineWidth',2)
axis square
hold on 
text(x+0.05 ,y ,str)
%%
figure
plot(cnt(:,1))

hold on 
pwelch(cnt(:,21));figure(gcf)

figure 
hold on
pwelch(cnt(:,21),[],[],[],fs);figure(gcf);% plotting powerspectrum 
title("Welch Power Spectral Density Estimate for ch is 21");

cwt(cnt(:,1));
cwt(cnt(mrk.pos(1)+400 :mrk.pos(1)+800,20));

cwt(cnt(mrk.pos(1)+400 :mrk.pos(1)+800,20),[],[],[],nfo.fs);% showing just one trail in wavelet for channel is 15
cwt(cnt(mrk.pos(1) :mrk.pos(1)+800,45),[],[],[],nfo.fs);title("plot wavelet one trail for channel 45");


%%
% bands=[1,4 ;4,8; 8,12; 12,30];
% for band=1:4
%     for ch=1:59
% disp(['Working on Bands:',num2str(band),'channel:',num2str(ch)]);
%      
%             %[b,a] =buuter(order,wn,ftype)
%             if band ==1
%                 [b,a] = butter(4,bands(band,2)/(fs/2));%make IIR filter (band pass filter )
%                   else
%                   [b,a] = butter(4,bands(band,:)/(fs/2));%make IIR filter (band pass filter )
% 
%             end
%         temp_cnt = filtfilt(b,a,cnt(:,ch));%we use filtfilt for problem of filt that make phase shift
% 
%         temp_cnt = (temp_cnt-mean(temp_cnt))/std(temp_cnt);%Normalize
%         
%     end
% end
%% step1 : proprocessing
% CAR 
mean_cnt = mean(cnt,2);

for i=1:59 
 temp = cnt(:,i);
 cnt(:,i) = temp - mean_cnt;
end


%% Design filter 
% bandpass butterworth filter

fl= 0.1;
fh= 4;
wn=[fl fh] / (fs/2);
[b,a]= butter(3,wn,'bandpass'); % delta band
cnt_delta= filtfilt(b,a,cnt);
save('delta.mat' ,'cnt_delta');
plot(cnt_delta(:,3));title("plot channel 3 of delta rythm")

%%
[b2,a2]= butter(3,[4  8]/(fs/2),'bandpass');% theta band
cnt_theta= filtfilt(b2,a2,cnt);
save('theta.mat' ,'cnt_theta');
plot(cnt_theta(:,3));title("plot channel 3 of theta rythm")
%%
[b3,a3]= butter(3,[8  12]/(fs/2),'bandpass');% MU band
cnt_alpha= filtfilt(b3,a3,cnt);
plot(cnt_alpha(:,3));title("plot channel 3 of alpha rythm")
save('alpha.mat' ,'cnt_alpha');
%%
[b4,a4]= butter(3,[12 30]/(fs/2),'bandpass');% Beta band
cnt_betha= filtfilt(b4,a4,cnt);
plot(cnt_betha(:,3));title("plot channel 3 of betha rythm")
save('betha.mat' ,'cnt_betha');


%% Fourier transformer
N= length(cnt);
fx1= fft(cnt)';
fx1= abs(fx1(1:round(N/2)));
rf= linspace(0,fs/2,round(N/2)); 
figure
stem(rf,fx1,'b','linewidth',1,'marker','none')
grid on
grid minor
title(' showing FFT signal before filtering')

hold on

N= length(cnt_delta);
fx1= fft(cnt_delta)';
fx1= abs(fx1(1:round(N/2)));
rf= linspace(0,fs/2,round(N/2)); 
figure 
stem(rf,fx1,'b','linewidth',1,'marker','none')
grid on
grid minor
title('showing FFT signal of delta rythm after filtering')
%%
%removing 50 hz noise
%[b,a] = butter(4, [ 49.5 51.5 ]/((nfo.fs/2)) ,"stop");%design notch filter 

%cnt= filtfilt(b,a ,cnt);

%%
% select data in beta and mu band
[b,a] = butter(4, [8 30]/((nfo.fs/2)));

new_cnt= filtfilt(b,a ,cnt);

pwelch(new_cnt(:,1));figure(gcf);
pwelch(new_cnt(:,1),[],[],[],fs);figure(gcf);title("Welch Power Spectral Density Estimate for channel 1")

%%
% select data in  mu band
[b,a] = butter(4, [6 13]/((nfo.fs/2)));

mu_cnt= filtfilt(b,a ,cnt);
mu = mu_cnt;
save mu;
pwelch(mu_cnt(:,1));figure(gcf);
pwelch(mu_cnt(:,1),[],[],[],fs);figure(gcf);title("Welch Power Spectral Density Estimate of mu band for channel 1")

%%
% select data in  beta band
[b2,a2] = butter(4, [12 30]/((nfo.fs/2)));

beta_cnt= filtfilt(b2,a2 ,cnt);
beta = beta_cnt ;
save beta
pwelch(beta_cnt(:,1));figure(gcf);
pwelch(beta_cnt(:,1),[],[],[],fs);figure(gcf);title("Welch Power Spectral Density Estimate od beta band for channel 1")


%% feature selection
mean_cnt = mean(new_cnt);
ZeroMeanCnt = new_cnt - mean_cnt ;
CovMatrix = cov(ZeroMeanCnt);

Number_EigVec = 35; % we choose only 35 principal components
[EigVec , EigVal] = eig(CovMatrix);
[Sort_EigVec , Sort_Index] = sort(diag(EigVal),'descend');
chooseEigVec = Sort_Index(1:Number_EigVec);%only 35 principal components
cnt_pca = cnt * EigVec(:,chooseEigVec);
%% channel selection
channel_index = 26:32; % select channels :c5 c3 c1  cz c2 c4 c6
channel_index(4) = []; % ignore information of Cz
cnt_ch = new_cnt(:,channel_index);
%%
number_class1 = sum(mrk.y == 1); %right  hand
number_class2 = sum(mrk.y == -1); % left hand 

class1 = cell(1,number_class1);
class2 = cell(1 , number_class2);

%% sate 1 when using PCA
%length of mrk.pos = 200
for i =1:length(mrk.pos)
    if mrk.y(i) ==1
        class1{i}= cnt_pca(mrk.pos(i)+200 : mrk.pos(i)+600,:);
    elseif mrk.y(i) ==-1
        class2{i}= cnt_pca(mrk.pos(i)+200 : mrk.pos(i)+600,:);
    end
  
end
%% state 2 when using channel selection
for i =1:length(mrk.pos)
    if mrk.y(i) ==1
        class1{i}= cnt_ch(mrk.pos(i)+200 : mrk.pos(i)+600,:);
    elseif mrk.y(i) ==-1
        class2{i}= cnt_ch(mrk.pos(i)+200 : mrk.pos(i)+600,:);
    end
  
end

%%
class1 = class1(~cellfun ('isempty' , class1));
class2 = class2(~cellfun ('isempty' , class2));

%% halete 2 for feature extraction 
for i = 1:100
    var1_delta(i, :) = var(delta1{i});
    var1_theta(i, :) = var(theta1{i});
    var1_mu(i, :) = var(mu1{i});
    var1_betha(i, :) = var(beta1{i});
    
    kur1_delta(i, :) = kurtosis(delta1{i});
    kur1_theta(i, :) = kurtosis(theta1{i});
    kur1_mu(i, :) = kurtosis(mu1{i});
    kur1_beta(i, :) = kurtosis(beta1{i});
    
    power1_delta(i, :) = sum(abs(fft(delta1{i})) .^ 2);
    power1_theta(i, :) = sum(abs(fft(theta1{i})) .^ 2);
    power1_mu(i, :) = sum(abs(fft(mu1{i})) .^ 2);
    power1_beta(i, :) = sum(abs(fft(beta1{i})) .^ 2);
end

for i = 1:n2
    var2_delta(i, :) = var(delta2{i});
    var2_theta(i, :) = var(theta2{i});
    var2_mu(i, :) = var(mu2{i});
    var2_beta(i, :) = var(beta2{i});
    
    kur2_delta(i, :) = kurtosis(delta2{i});
    kur2_theta(i, :) = kurtosis(theta2{i});
    kur2_mu(i, :) = kurtosis(mu2{i});
    kur2_beta(i, :) = kurtosis(beta2{i});
    
    power2_delta(i, :) = sum(abs(fft(delta2{i})) .^ 2);
    power2_theta(i, :) = sum(abs(fft(theta2{i})) .^ 2);
    power2_mu(i, :) = sum(abs(fft(mu2{i})) .^ 2);
    power2_beta(i, :) = sum(abs(fft(beta2{i})) .^ 2);
    
end
%% halete 2 for feature extraction 
% create feature using variance
for i =1:length(class2)

 class1_FeatureVar (i,:) = var(class1{i});
 class1_FeatureKur(i,:)   = kurtosis(class1{i});

 class2_FeatureVar (i,:) = var(class2{i});
 class2_FeatureKur(i,:)   = kurtosis(class2{i});
end

%%
figure
scatter(var_delta, var_alpha)
xlabel('Delta Variance')
ylabel('Alpha Variance')
title('Scatter  Alpha Var vs Delta Var')

figure
scatter(kur_delta, kur_alpha)
xlabel('Delta Kur')
ylabel('Alpha Kur')
title('Scatter of Alpha Kur vs Delta Kur')

figure
scatter(power_delta, power_alpha)
xlabel('Delta Power')
ylabel('Alpha Power')
title('Scatter  of Alpha Power vs Delta Power')
%% creation and cancat features & labels
features = [class1_Feature ;class2_Feature];%concatition
label = [zeros(length(class1_Feature),1);
         ones(length(class2_Feature),1)];

%% state when using kurtosis
features_class1 = [class1_FeatureVar;class1_FeatureKur];
features_class2 = [class2_FeatureVar;class2_FeatureKur];%concatition

features = [features_class1 ; features_class2];

label = [zeros(length(features_class1),1);
         ones(length(features_class2),1)];
%% ploting & Scattering plot 

figure(1)
scatter(class1_FeatureVar, var_alpha)
xlabel('Delta Variance')
ylabel('Alpha Variance')
title('Scatter Plot of Alpha Variance vs Delta Variance')
%% %shuflling 

%randperm(200)
idx = randperm(200);
features = features(idx , :);
label = label(idx);
%% % Classification using SVM
SVMModel = fitcsvm(features , label , 'KernelFunction','gaussian', 'ClassNames',[0;1]);

% validation using k-fold cross validation
PartitionModel  = crossval(SVMModel , 'kfold' , 5);
ACC=1-kfoldLoss(PartitionModel , 'LossFun','classiferror');

disp(['ACC for using SVM model is :',num2str(ACC)])

%% classificationLearner
DataFrame = [features , label];

%% Classification using KNN
KnnModel = fitcknn(features , label , "NumNeighbors",5);
% validation
PartitionModel  = crossval(KnnModel , 'kfold' , 5);
ACC=1-kfoldLoss(PartitionModel , 'LossFun','classiferror');

disp(['ACC for using KNN model is :',num2str(ACC)])