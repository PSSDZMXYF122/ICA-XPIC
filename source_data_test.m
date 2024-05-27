clear;clc;
close all;

%% n = 4;  % 设定调制数目，并生成 "原始信号序列"
% n = 2;
% N = 5000;       %符号数
% bitdata1=randi([0 n-1],N,1);   %用于生成源调制信号，符号数为5000.
% bitdata2=randi([0 n-1],N,1);
% 
% sn_unnoise = pskmod(bitdata1,n,pi/n);  % 生成有相位偏移的PSK信号
% gn_unnoise = pskmod(bitdata2,n,pi/n);

% sn_unnoise = pskmod(bitdata1,n);  % 生成无相位偏移的PSK信号
% gn_unnoise = pskmod(bitdata2,n);
% 
% sn_unnoise = qammod(bitdata1,n);  % 生成无相位偏移的PSK信号
% gn_unnoise = qammod(bitdata2,n);


% mixed_data_unnoised_sn = 1 * sn_unnoise + alpha_mix * gn_unnoise;  %根据信号模型进行混合 X=A*S;
% mixed_data_unnoised_gn = alpha_mix * sn_unnoise + 1 * gn_unnoise;
% 
ii = sqrt(-1);
mix_matrix = [0.5222 + 0.2605*ii,0.0872 + 0.4251*ii 0.7080+0.3895*ii;0.0696 + 0.8375*ii, 0.9118 + 0.0278*ii 0.6470+0.7243*ii;0.9201+0.3909*ii 0.9925+0.9317*ii 0.1794+0.5786*ii];

% mix_matrix = [1,0.2;0.2, 1];

t = linspace(0,0.02,1000);
signal_1_real = cos(100*pi*t);
signal_1_imag = sin(100*pi*t);
signal_1 = signal_1_real + sqrt(-1)*signal_1_imag;

n = 2;
N = 1000;       %符号数
bitdata1=randi([0 n-1],N,1);   %用于生成源调制信号，符号数为5000.
signal_2 = (pskmod(bitdata1,n)).';
n_2 = 4;
bitdata2=randi([0 n_2-1],N,1);   %用于生成源调制信号，符号数为5000.
signal_3 = qammod(bitdata2,4);

signal(1,:) = signal_1;
signal(3,:) = signal_2.';
signal(2,:) = sqrt(2)/2*signal_3.';

figure(1);
subplot(311);
plot(signal(1,:).','.');
subplot(312);
plot(signal(2,:).','.');
% axis([-1.5 1.5 -1.5 1.5]);

subplot(313);
plot(signal(3,:).','.');
axis([-1.5 1.5 -1.5 1.5]);
% xlim(-[2 2]);ylim([-1.5 1.5]);



mixed_signal = mix_matrix * signal;
figure(2);
subplot(311);
plot(mixed_signal(1,:).','.');
subplot(312);
plot(mixed_signal(2,:).','.');
% axis([-1.5 1.5 -1.5 1.5]);

subplot(313);
plot(mixed_signal(3,:).','.');




%%%%%%%%%% 算法迭代部分 %%%%%%%%%%
[B,Q,iteration_num_f1,initial_w] = fastica_achieve1(mixed_signal);
% B_FORLOOK = B;
aa = B'*Q; 
   

ICAedS = aa*mixed_signal;   %得到分离后的输出信号
IS_eye = aa*mix_matrix;

%     qpsk_rmsEVM_separated_f1(i,:) = evm(source_data_unnoised.',qpsk_ICAedS(:,:,i).');
%     qpsk_rmsEVM_separated_f1(i,:) = 20*log10(0.01*qpsk_rmsEVM_separated_f1(i,:));
%     rmsEVM_subtract_f1(i,:) = qpsk_rmsEVM_separated_f1(i,:) - qpsk_rmsEVM_unseparated(i,:);



figure(3);
% plot(ICAedS(1,:).','.r');hold on;
% plot(ICAedS(2,:).','.b');
% 
% plot(signal_1.','xg');
% plot(signal_2.','pk')
% axis([-2 2 -2 2]);
subplot(311);
plot(ICAedS(1,:).','.');
subplot(312);
plot(ICAedS(2,:).','.');
% axis([-1.5 1.5 -1.5 1.5]);

subplot(313);
plot(ICAedS(3,:).','.');


figure(4);
subplot(331);
plot(signal(1,:).','.');
subplot(332);
plot(signal(2,:).','.');
% axis([-1.5 1.5 -1.5 1.5]);
subplot(333);
plot(signal(3,:).','.');
axis([-1.5 1.5 -1.5 1.5]);
%% %性能指数（PI）
% PI_separation = 1/2*(abs(IS_eye(1,2))/abs(IS_eye(1,1)) +abs(IS_eye(2,1))/abs(IS_eye(2,2))+abs(IS_eye(2,1))/abs(IS_eye(1,1)) + abs(IS_eye(1,2))/abs(IS_eye(2,2)) )
% PI_separation = 


subplot(334);
plot(mixed_signal(1,:).','.');
subplot(335);
plot(mixed_signal(2,:).','.');
% axis([-1.5 1.5 -1.5 1.5]);

subplot(336);
plot(mixed_signal(3,:).','.');


subplot(337);
plot(ICAedS(1,:).','.');
subplot(338);
plot(ICAedS(2,:).','.');
% axis([-1.5 1.5 -1.5 1.5]);

subplot(339);
plot(ICAedS(3,:).','.');


% 
% 
% source_data_unnoised(1,:) = sn_unnoise.';
% source_data_unnoised(2,:) = gn_unnoise.';
% %% 生成特定信噪比的噪声。
% P_sn = sum(sn_unnoise.*conj(sn_unnoise))/N;     %计算出信号功率
% P_gn = sum(gn_unnoise.*conj(gn_unnoise))/N;
% 
% dB_noise = 10:1:25;  %实验结果图像的横坐标，设定混合后的信号的信噪比
% num_noise = length(dB_noise);
% % P_noise_sn = P_sn*10^(-15/10);      %设定需要信噪比
% % P_noise_gn = P_gn*10^(-15/10);
% % noise_real_img = normrnd(0,1,[5000,4]);   %用于生成均值为0，方差为1的高斯白噪声，符号数为5000.
% % 
% % noise_1 = sqrt(P_noise_sn)*(noise_real_img(:,1)+sqrt(-1)*noise_real_img(:,2))/sqrt(2);
% % noise_2 = sqrt(P_noise_gn)*(noise_real_img(:,3)+sqrt(-1)*noise_real_img(:,4))/sqrt(2);
% % mean(noise_1)
% % disp(var(noise_1));
% 
% %% %计算分离前的evm
% evm = comm.EVM();
% mixed_data_noised = zeros(2,N,num_noise);
% qpsk_rmsEVM_unseparated = zeros(length(dB_noise),2);
% 
% 
% for i = 1:length(dB_noise)    
%     P_noise_sn(i) = P_sn*10^(-dB_noise(i)/10);      %设定需要信噪比
%     P_noise_gn(i) = P_gn*10^(-dB_noise(i)/10);
%     noise_real_img = normrnd(0,1,[5000,4]);   %用于生成均值为0，方差为1的高斯白噪声，符号数为5000.
%     
%     noise_1(:,i) = sqrt(P_noise_sn(i))*(noise_real_img(:,1)+sqrt(-1)*noise_real_img(:,2))/sqrt(2);
%     noise_2(:,i) = sqrt(P_noise_gn(i))*(noise_real_img(:,3)+sqrt(-1)*noise_real_img(:,4))/sqrt(2);
%     
%     %% 设定交叉极化混合系数
%     % alpha_mix = 10^(-8/20);     % 0.3981 8dB    可以用
%     % alpha_mix = 10^(-10/20)*exp(sqrt(-1)*pi/6);    % 0.3162 10dB     可以用
% %     alpha_mix = 10^(-10/20)*exp(sqrt(-1)*pi/6);
%     alpha_mix = 10^(-10/20)*exp(sqrt(-1)*pi/30);
% %      alpha_mix = 10^(-10/20);
%     % alpha_mix = 10^(-15/20);    % 0.1778 15dB     可以用
%     % alpha_mix = 10^(-20/20);    % 0.1000  20dB   可以用
%     % alpha_mix = 10^(-25/20);    %  0.0562 25dB  不可能用于仿真
%     % alpha_mix = 10^(-30/20);    %  0.0316  30dB  不可能用于仿真
%     mix_matrix = [1,alpha_mix;alpha_mix,1];
%     
%     mixed_data_unnoised_sn = 1 * sn_unnoise + alpha_mix * gn_unnoise;  %根据信号模型进行混合 X=A*S;
%     mixed_data_unnoised_gn = alpha_mix * sn_unnoise + 1 * gn_unnoise;
%     
%     mixed_data_noised_sn = mixed_data_unnoised_sn + noise_1(:,i);   %给信号加噪声
%     mixed_data_noised_gn = mixed_data_unnoised_gn + noise_2(:,i);   %给信号加噪声
% 
%     mixed_data_noised(1,:,i) = mixed_data_noised_sn.';    %% 将两路信号用一个矩阵表示
%     mixed_data_noised(2,:,i) = mixed_data_noised_gn.';
% 
% %     mixed_data_noised(:,:,i) = awgn(mixed_data_unnoised,dB_noise(i),'measured');   %给信号加噪声
% 
%     qpsk_rmsEVM_unseparated(i,:) = evm(source_data_unnoised.',mixed_data_noised(:,:,i).');
%     qpsk_rmsEVM_unseparated(i,:) = 20*log10(0.01*qpsk_rmsEVM_unseparated(i,:)) ;
%     % rmsEVM定义为平均误差矢量功率与平均基准功率的比值的平方根
% %     disp(['混合信号加噪声后的evm分别为：  ' num2str(rmsEVM_unseparated(i,:)) ' dB']);
% 
% 
% end
% 
% 
% 
% 
% save source_data_qpsk.mat qpsk_rmsEVM_unseparated mixed_data_noised sn_unnoise gn_unnoise  dB_noise source_data_unnoised alpha_mix mix_matrix
% 
% 
% 
% % 
% % subplot(121);plot(mixed_data_noised(1,:,6),'.'); hold on;
% % plot(source_data_unnoised(:,1).','xr') 
% % % axis([-2 2 -1 1]); 
% % title('XPD=10dB,SNR=15dB,un-ICAed consequence');
% 
% % subplot(122);plot(qpsk_ICAedS(1,:,6).','.'); hold on;
% % plot(source_data_unnoised(1,:),'xr') % axis([-2 2 -1 1]); 
% % title('XPD=10dB,SNR=15dB,ICAed consequence');
% 


