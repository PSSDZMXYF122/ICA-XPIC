clear;clc;
close all;


load source_data_bpsk.mat ;

evm = comm.EVM();

%% 使用基于负熵的FastICA算法的分离  标记为 f1

for i = 1:length(dB_noise)
    %%%%%%%%%% 算法迭代部分 %%%%%%%%%%
    [B,Q,iteration_num_f1,initial_w,W_unchanged] = fastica_achieve1(mixed_data_noised(:,:,i));
    B_FORLOOK(:,:,i) = B;
    aa(:,:,i) = B'*Q; 
    bpsk_ICAedS(:,:,i) = aa(:,:,i)*mixed_data_noised(:,:,i);   %得到分离后的输出信号
    bpsk_rmsEVM_separated(i,:) = evm(source_data_unnoised.',bpsk_ICAedS(:,:,i).');
    bpsk_rmsEVM_separated(i,:) = 20*log10(0.01*bpsk_rmsEVM_separated(i,:));
    rmsEVM_subtract_f1(i,:) = bpsk_rmsEVM_separated(i,:) - bpsk_rmsEVM_unseparated(i,:);
end


%% %%%%%%%% 绘制第1张 alpha-EVM图 %%%%%%%%%%
figure(1);
plot(dB_noise',bpsk_rmsEVM_unseparated(:,2),'-or');hold on;
grid on;
grid minor;
plot(dB_noise',bpsk_rmsEVM_separated(:,2),'-ok');
% plot(dB_noise',rmsEVM_subtract_f1(:,1),'-*m',dB_noise',rmsEVM_subtract_f1(:,2),'-*k');

xlabel('SNR (dB)');ylabel('EVM (dB)');
legend('unseparation gn','f1 separated gn');
% title('xpd = 10dB, alpha = sqrt(0.1),different alpha && EVM');



% % figure(2); scatter(real(mixed_data_noised.'),imag(mixed_data_noised.'));
% figure(1); subplot(121);plot(mixed_data_noised(1,:,12),'.'); hold on;
% plot(source_data_unnoised(1,:),'xr') % axis([-2 2 -1 1]); 
% title('un-ICAed consequence');
% 
% % scatter(real(ICAedS.'),imag(ICAedS.'));
% subplot(122);plot(ICAedS(1,:,16),'.'); hold on;
% plot(source_data_unnoised(1,:),'xr') % axis([-2 2 -1 1]); 
% title('ICAed consequence');
% 
% scatter(real(sn_unnoise),imag(sn_unnoise));


%% %%%%%%%% 绘制第2张 ICA分离前后alpah变化图 %%%%%%%%%%
num = length(dB_noise);
bpsk_aa_2 = zeros(num,2);
bpsk_aa_3 = zeros(2,2,num);
for i = 1:num
%     aa_3(:,:,i) = real (aa(:,:,i)*mix_matrix); 
    bpsk_aa_3(:,:,i) = aa(:,:,i)*mix_matrix; 
    bpsk_aa_2(i,1) = bpsk_aa_3(2,2,i); 
    bpsk_aa_2(i,2) = bpsk_aa_3(2,1,i);

end



figure(2);

plot(dB_noise',20*log10(abs(bpsk_aa_2(:,1))./abs(bpsk_aa_2(:,2))),'-og');
hold on;
alpha_db = ones(length(dB_noise),1)*alpha_mix;
plot(dB_noise',-20*log10(abs(alpha_db)),'-or')
grid on; grid minor;
xlabel('SNR (dB)');ylabel('separated XPD (dB)');
legend('f1 separated','reference');

% title('xpd = 10dB alpha=sqrt(0.1), Input & Output alpha difference');



%% %%%%%%%% 绘制第3张 alpha-SINR图 %%%%%%%%%%
%计算分离的性能
P_signal = sum((source_data_unnoised.*conj(source_data_unnoised)).');%计算信号的功率P_signal的两个元素分别为sn和gn的功率
%计算SINR
% dB_SINR = zeros(length(dB_noise),2);

for i =1:length(dB_noise)
    % 计算未分离的SINR
    bpsk_interference_noise_unseparation(:,:,i) = mixed_data_noised(:,:,i) - source_data_unnoised;
    bpsk_P_ganrao_plus_noise_unseparation(i,:) = sum((bpsk_interference_noise_unseparation(:,:,i) .*conj(bpsk_interference_noise_unseparation(:,:,i) )).');    
    bpsk_dB_SINR_unseparation(i,:) = 10*log10(P_signal./bpsk_P_ganrao_plus_noise_unseparation(i,:));  

    % 计算 分离的SINR    f1
    interference_noise_separation_f1(:,:,i) = bpsk_ICAedS(:,:,i) - source_data_unnoised;
    bpsk_P_ganrao_plus_noise_separation_f1(i,:) = sum((interference_noise_separation_f1(:,:,i) .*conj(interference_noise_separation_f1(:,:,i) )).');    
    bpsk_dB_SINR_separation_f1(i,:) = 10*log10(P_signal./bpsk_P_ganrao_plus_noise_separation_f1(i,:)); 

end



figure(3);
plot(dB_noise',bpsk_dB_SINR_unseparation(:,2),'-.pr');
hold on;grid on;
grid minor;
plot(dB_noise',bpsk_dB_SINR_separation_f1(:,2),'-pg');

xlabel('SNR (dB)');ylabel('SINR (dB)');
legend('unseparated gn');
% title('SNR=20dB, different SNR && SINR');


% %% %%%%%%%% 绘制第4张 alpha-迭代次数图 %%%%%%%%%%
% %
% 
% figure(4);
% plot(dB_noise,iteration_num_f1,'-or');
% hold on;grid on;
% grid minor;
% 
% xlabel('SNR (dB)');ylabel('Time of iteration');
% legend('f1');
% % title('different alpha && 迭代次数');

%% %性能指数（PI）
for i = 1:length(dB_noise)
%     PI_unseparation(i) = 1/2*abs(mix_matrix(1,2,i))/abs(mix_matrix(1,1,i)) + ...
%         abs(mix_matrix(2,1,i))/abs(mix_matrix(2,2,i)) + ...
%         abs(mix_matrix(2,1,i))/abs(mix_matrix(1,1,i)) + ...
%         abs(mix_matrix(1,2,i))/abs(mix_matrix(2,2,i)) ;
%    PI_unseparation = 1/2*4*ones(1,length(dB_noise))*0.1^0.5;
     PI_unseparation = 1/2*4*ones(1,length(dB_noise))*alpha_mix;
    bpsk_PI_separation_1(i) = 1/2*(abs(bpsk_aa_3(1,2,i))/abs(bpsk_aa_3(1,1,i)) +abs(bpsk_aa_3(2,1,i))/abs(bpsk_aa_3(2,2,i)) + ...
        abs(bpsk_aa_3(2,1,i))/abs(bpsk_aa_3(1,1,i)) + abs(bpsk_aa_3(1,2,i))/abs(bpsk_aa_3(2,2,i)) );

end

%%%%%%%%%% 绘制第5张 SNR-PI(performance index)图 %%%%%%%%%%
figure(5);
plot(dB_noise,PI_unseparation,'-ko',dB_noise,bpsk_PI_separation_1,'-bp');
hold on;
plot(dB_noise,ones(1,length(dB_noise))*0.01,'-r');
grid on;grid minor;
% xlim([10 30]);ylim([0.05 0.4]);
% title('different SNR & PI');
xlabel('SNR (dB)');ylabel('PI');
legend('PI unseparation','f1 PI separation','performance boundry');

% save bpsk_0.1_20db_data.mat bpsk_ICAedS bpsk_rmsEVM_separated bpsk_rmsEVM_unseparated ...
%     bpsk_aa_2 ...
%     bpsk_dB_SINR_unseparation bpsk_dB_SINR_separation_f1 ...
%     bpsk_PI_separation_1



figure(6); 

subplot(121);plot(mixed_data_noised(1,:,6),'.'); hold on;
plot(source_data_unnoised(1,:),'xr') 
% axis([-2 2 -1 1]); 
% title('XPD=10dB,SNR=15dB,un-ICAed consequence');

subplot(122);plot(bpsk_ICAedS(1,:,6),'.'); hold on;
plot(source_data_unnoised(1,:),'xr') % axis([-2 2 -1 1]); 
% title('XPD=10dB,SNR=15dB,ICAed consequence');









