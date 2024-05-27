
%% 改进后的模型  先混合再加入噪声。

clear;clc;
close all;

%% n = 4;  % 设定调制数目，并生成 "原始信号序列"
n = 4;
% N = 1e6;       %符号数

N = 5000;       %设定符号数

for m=1:100


bitdata1=randi([0 n-1],N,1);   %用于生成源调制信号，符号数为5000.
bitdata2=randi([0 n-1],N,1);
bitdata(1,:) = bitdata1';
bitdata(2,:) = bitdata2';

sn_unnoise = pskmod(bitdata1,n,pi/n);  % 生成有相位偏移的PSK信号
gn_unnoise = pskmod(bitdata2,n,pi/n);

% sn_unnoise = pskmod(bitdata1,n);  % 生成无相位偏移的PSK信号
% gn_unnoise = pskmod(bitdata2,n);

% sn_unnoise = qammod(bitdata1,n);  % 生成无相位偏移的PSK信号
% gn_unnoise = qammod(bitdata2,n);

source_data_unnoised(1,:) = sn_unnoise.';
source_data_unnoised(2,:) = gn_unnoise.';
%% 生成特定信噪比的噪声。
P_sn = sum(sn_unnoise.*conj(sn_unnoise))/N;     %计算出信号功率
P_gn = sum(gn_unnoise.*conj(gn_unnoise))/N;

% dB_noise = 10:1:25;  %实验结果图像的横坐标，设定混合后的信号的信噪比
dB_noise = 10:1:20;  %实验结果图像的横坐标，设定混合后的信号的信噪比
num_noise = length(dB_noise);

%% %计算分离前的evm
evm = comm.EVM();
mixed_data_noised = zeros(2,N,num_noise);
qpsk_rmsEVM_unseparated = zeros(length(dB_noise),2);


for i = 1:length(dB_noise)    
    P_noise_sn(i) = P_sn*10^(-dB_noise(i)/10);      %设定需要信噪比
    P_noise_gn(i) = P_gn*10^(-dB_noise(i)/10);
    noise_real_img = normrnd(0,1,[N,4]);   %用于生成均值为0，方差为1的高斯白噪声，符号数为5000.
    
    noise_1(:,i) = sqrt(P_noise_sn(i))*(noise_real_img(:,1)+sqrt(-1)*noise_real_img(:,2))/sqrt(2);
    noise_2(:,i) = sqrt(P_noise_gn(i))*(noise_real_img(:,3)+sqrt(-1)*noise_real_img(:,4))/sqrt(2);
    
    %% 设定交叉极化混合系数
    % alpha_mix = 10^(-8/20);     % 0.3981 8dB    可以用
    % alpha_mix = 10^(-10/20)*exp(sqrt(-1)*pi/6);    % 0.3162 10dB     可以用
%     alpha_mix = 10^(-10/20)*exp(sqrt(-1)*pi/6);
%     alpha_mix = 10^(-10/20);
%     alpha_mix = 10^(-5/20);
    alpha_mix = 10^(-10/20)*exp(sqrt(-1)*pi/8);       %2*PI/5
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

    qpsk_rmsEVM_unseparated(i,:) = evm(source_data_unnoised.',mixed_data_noised(:,:,i).');
    qpsk_rmsEVM_unseparated(i,:) = 20*log10(0.01*qpsk_rmsEVM_unseparated(i,:)) ;
    % rmsEVM定义为平均误差矢量功率与平均基准功率的比值的平方根
%     disp(['混合信号加噪声后的evm分别为：  ' num2str(rmsEVM_unseparated(i,:)) ' dB']);


end



for i = 1:length(dB_noise)
    %%%%%%%%%% 算法迭代部分 %%%%%%%%%%
    t1 = clock;%tic;
    [B,Q,iteration_num_f1] = fastica_achieve1(mixed_data_noised(:,:,i));
    t2 = clock;%toc;
    t_cfastica(i) = etime(t2,t1);
    B_FORLOOK(:,:,i) = B;
    Q_FORLOOK(:,:,i) = Q;
    aa(:,:,i) = B'*Q; 
end

%% %%%%%%%% 绘制第2张 ICA分离前后XPD变化图 %%%%%%%%%%
num = length(dB_noise);
qpsk_aa_2 = zeros(num,2);
qpsk_aa_3 = zeros(2,2,num);
for i = 1:num
%     aa_3(:,:,i) = real (aa(:,:,i)*mix_matrix); 
    qpsk_aa_3(:,:,i) = aa(:,:,i)*mix_matrix;  
    qpsk_aa_2(i,1) = qpsk_aa_3(2,2,i); 
    qpsk_aa_2(i,2) = qpsk_aa_3(2,1,i);

end

XPD = 20*log10(abs(qpsk_aa_2(:,1))./abs(qpsk_aa_2(:,2)));



% string_expression = 'A'+ string(m);
% filename = 'XPD_average.xlsx';
% xlswrite(filename,XPD',1,string_expression);

XPD_Array(m,:)=XPD;

end

XPD_Array_mean = mean(XPD_Array);




figure(1);
plot(dB_noise',XPD_Array_mean,'-og');
hold on;
alpha_db = ones(length(dB_noise),1)*alpha_mix;
plot(dB_noise',-20*log10(abs(alpha_db)),'-or');
ylim([10 45]);
grid on; grid minor;
xlabel('SNR (dB)');ylabel('separated XPD (dB)');
legend('CFastICA separated XPD','Unseparated XPD');
title('给定参数下,研究SNR对分离后XPD的影响');




