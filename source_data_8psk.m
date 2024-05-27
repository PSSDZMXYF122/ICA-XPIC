% clear;clc;
% 
% n = 8;  % 设定调制数目
% 
% bitdata1=randi([0 n-1],5000,1);   %用于生成0,1原始符号，符号数为5000.
% bitdata2=randi([0 n-1],5000,1);
% 
% sn_unnoise = pskmod(bitdata1,n,pi/n);  % 生成有相位偏移的PSK信号
% gn_unnoise = pskmod(bitdata2,n,pi/n);
% 
% % sn_unnoise = pskmod(bitdata1,n);  % 生成无相位偏移的PSK信号
% % gn_unnoise = pskmod(bitdata2,n);
% 
% source_data_unnoised = [sn_unnoise.';gn_unnoise.'] ;  %非共轭转置
% 
% 
% alpha_mix = 10^(-10/20);    % 0.3162 10dB     可以用
% mix_matrix = [1,alpha_mix;alpha_mix,1];
% 
% mixed_data_unnoised = mix_matrix*source_data_unnoised;  %根据信号模型进行混合 X=A*S;
% 
% %% %计算分离前的evm
% evm = comm.EVM();
% dB_noise = 10:1:25;  %实验结果图像的横坐标，设定混合后的信号的信噪比
% 
% % mixed_data_noised = awgn(mixed_data_unnoised,dB_noise,'measured');   %给信号加噪声
% % rmsEVM_unseparated = evm(source_data_unnoised.',mixed_data_noised.');
% % rmsEVM_unseparated = 20*log10(0.01*rmsEVM_unseparated) ;
% %     % rmsEVM定义为平均误差矢量功率与平均基准功率的比值的平方根
% % %     disp(['混合信号加噪声后的evm分别为：  ' num2str(rmsEVM_unseparated(i,:)) ' dB']);
% 
% 
% psk_8_rmsEVM_unseparated = zeros(length(dB_noise),2);
% 
% 
% for i = 1:length(dB_noise)
%     mixed_data_noised(:,:,i) = awgn(mixed_data_unnoised,dB_noise(i),'measured');   %给信号加噪声
%     psk_8_rmsEVM_unseparated(i,:) = evm(source_data_unnoised.',mixed_data_noised(:,:,i).');
%     psk_8_rmsEVM_unseparated(i,:) = 20*log10(0.01*psk_8_rmsEVM_unseparated(i,:)) ;
%     % rmsEVM定义为平均误差矢量功率与平均基准功率的比值的平方根
% %     disp(['混合信号加噪声后的evm分别为：  ' num2str(rmsEVM_unseparated(i,:)) ' dB']);
% end
% 
% save source_data_8psk.mat psk_8_rmsEVM_unseparated mixed_data_noised sn_unnoise gn_unnoise  dB_noise source_data_unnoised alpha_mix mix_matrix







%% 改进后的模型

clear;clc;

%% n = 4;  % 设定调制数目，并生成 "原始信号序列"
n = 8;
N = 5000;       %符号数
bitdata1=randi([0 n-1],N,1);   %用于生成源调制信号，符号数为5000.
bitdata2=randi([0 n-1],N,1);

sn_unnoise = pskmod(bitdata1,n,pi/n);  % 生成有相位偏移的PSK信号
gn_unnoise = pskmod(bitdata2,n,pi/n);

% sn_unnoise = pskmod(bitdata1,n);  % 生成无相位偏移的PSK信号
% gn_unnoise = pskmod(bitdata2,n);
% 
% sn_unnoise = qammod(bitdata1,n);  % 生成无相位偏移的PSK信号
% gn_unnoise = qammod(bitdata2,n);

source_data_unnoised(1,:) = sn_unnoise.';
source_data_unnoised(2,:) = gn_unnoise.';
%% 生成特定信噪比的噪声。
P_sn = sum(sn_unnoise.*conj(sn_unnoise))/N;     %计算出信号功率
P_gn = sum(gn_unnoise.*conj(gn_unnoise))/N;

dB_noise = 10:1:25;  %实验结果图像的横坐标，设定混合后的信号的信噪比
num_noise = length(dB_noise);
% P_noise_sn = P_sn*10^(-15/10);      %设定需要信噪比
% P_noise_gn = P_gn*10^(-15/10);
% noise_real_img = normrnd(0,1,[5000,4]);   %用于生成均值为0，方差为1的高斯白噪声，符号数为5000.
% 
% noise_1 = sqrt(P_noise_sn)*(noise_real_img(:,1)+sqrt(-1)*noise_real_img(:,2))/sqrt(2);
% noise_2 = sqrt(P_noise_gn)*(noise_real_img(:,3)+sqrt(-1)*noise_real_img(:,4))/sqrt(2);
% mean(noise_1)
% disp(var(noise_1));

%% %计算分离前的evm
evm = comm.EVM();
mixed_data_noised = zeros(2,N,num_noise);
psk_8_rmsEVM_unseparated = zeros(length(dB_noise),2);


for i = 1:length(dB_noise)    
    P_noise_sn(i) = P_sn*10^(-dB_noise(i)/10);      %设定需要信噪比
    P_noise_gn(i) = P_gn*10^(-dB_noise(i)/10);
    noise_real_img = normrnd(0,1,[5000,4]);   %用于生成均值为0，方差为1的高斯白噪声，符号数为5000.
    
    noise_1(:,i) = sqrt(P_noise_sn(i))*(noise_real_img(:,1)+sqrt(-1)*noise_real_img(:,2))/sqrt(2);
    noise_2(:,i) = sqrt(P_noise_gn(i))*(noise_real_img(:,3)+sqrt(-1)*noise_real_img(:,4))/sqrt(2);
    
    %% 设定交叉极化混合系数
    % alpha_mix = 10^(-8/20);     % 0.3981 8dB    可以用
    % alpha_mix = 10^(-10/20)*exp(sqrt(-1)*pi/6);    % 0.3162 10dB     可以用
%     alpha_mix = 10^(-10/20)*exp(sqrt(-1)*pi/6);
%     alpha_mix = 10^(-10/20);
%     alpha_mix = 10^(-5/20);
    alpha_mix = 10^(-10/20)*exp(sqrt(-1)*pi/8);
    % alpha_mix = 10^(-15/20);    % 0.1778 15dB     可以用
    % alpha_mix = 10^(-20/20);    % 0.1000  20dB   可以用
    % alpha_mix = 10^(-25/20);    %  0.0562 25dB  不可能用于仿真
    % alpha_mix = 10^(-30/20);    %  0.0316  30dB  不可能用于仿真
    mix_matrix = [1,alpha_mix;alpha_mix,1];
    
    mixed_data_unnoised_sn = 1 * sn_unnoise + alpha_mix * gn_unnoise;  %根据信号模型进行混合 X=A*S;
    mixed_data_unnoised_gn = alpha_mix * sn_unnoise + 1 * gn_unnoise;
    
    mixed_data_noised_sn = mixed_data_unnoised_sn + noise_1(:,i);   %给信号加噪声
    mixed_data_noised_gn = mixed_data_unnoised_gn + noise_2(:,i);   %给信号加噪声

    mixed_data_noised(1,:,i) = mixed_data_noised_sn.';    %% 将两路信号用一个矩阵表示
    mixed_data_noised(2,:,i) = mixed_data_noised_gn.';

%     mixed_data_noised(:,:,i) = awgn(mixed_data_unnoised,dB_noise(i),'measured');   %给信号加噪声

    psk_8_rmsEVM_unseparated(i,:) = evm(source_data_unnoised.',mixed_data_noised(:,:,i).');
    psk_8_rmsEVM_unseparated(i,:) = 20*log10(0.01*psk_8_rmsEVM_unseparated(i,:)) ;
    % rmsEVM定义为平均误差矢量功率与平均基准功率的比值的平方根
%     disp(['混合信号加噪声后的evm分别为：  ' num2str(rmsEVM_unseparated(i,:)) ' dB']);


end




save source_data_8psk.mat psk_8_rmsEVM_unseparated mixed_data_noised sn_unnoise gn_unnoise  dB_noise source_data_unnoised alpha_mix mix_matrix













