clear all; close all;

load qpsk_0.1_20db_data.mat
load qpsk_bit_data.mat
%误码率曲线的绘制

SNR=0:1:20;%信噪比变化范围
% snr=10.^(SNR/10);%将信噪比转化成直角坐标

[~,N] = size(bitdata);% N仿真符号数

M=4;%进制数
% num_bit = N*log2(M);    %比特数

% errorRate = comm.ErrorRate;



% x=randi([0,1],1,N);%产生随机信号
% y=pskmod(x,M);%调用matlab自带的psk调制函数

% y1=real(y);
% y2=imag(y);
x = bitdata(1,:);

for i=1:length(SNR)

    yAn_separated = qpsk_ICAedS(1,:,i);
    yA_separated=pskdemod(yAn_separated,M,pi/4);%调用matlab自带的psk解调函数
    bit_A_separated = symerr(bitdata(1,:),yA_separated);
%     bit_A_separated=length(find(x~=yA_separated));%统计错误比特数 
    QPSK_s_AWGN_separated(i)=bit_A_separated/N;%计算误码率
%     err_separated(i) = biterr(x,yA_separated);
%     dataOut(:,:,i) = de2bi(yA_separated,2);


    yAn_unseparated = mixed_data_noised(1,:,i);
    yA_unseparated=pskdemod(yAn_unseparated,M,pi/4);%调用matlab自带的psk解调函数
    bit_A_unseparated = symerr(bitdata(1,:),yA_unseparated);
%     bit_A_unseparated=length(find(x~=yA_unseparated));%统计错误比特数 
    QPSK_s_AWGN_unseparated(i)=bit_A_unseparated/N;%计算误码率
    err_unseparated(i) = biterr(x,yA_unseparated);
%     dataIn(:,:,i) = de2bi(yA_unseparated,2);
% 
%     nErrors(i) = biterr(dataIn(:,:,i),dataOut(:,:,i));
%     berrr(i) = nErrors(i)/(2*N);
%     berVec = step(2,Rx_bit,Tx_bit);
%     brate(i) = berVec(1);

end
% QPSK_t_AWGN=1/2*erfc(sqrt(snr));%AWGN信道下QPSK理论误码率

%绘制图形
figure(1);
semilogy(SNR,QPSK_s_AWGN_unseparated,'-bx');hold on;

semilogy(SNR,QPSK_s_AWGN_separated,'-g*');hold on;

% semilogy(SNR,QPSK_t_AWGN,'-go');hold on;

% semilogy(SNR,err_unseparated,'-go');hold on;
% semilogy(SNR,err_separated,'-ro');hold on;

grid on;grid minor;
axis([-1,20,10^-7,1]);
legend('分离前仿真','分离后仿真','理论','理论');
title('QPSK误码性能分析');
xlabel('信噪比（dB）');ylabel('BER');

% QPSK_s_AWGN_unseparated_xpd_15 = QPSK_s_AWGN_unseparated;
% QPSK_s_AWGN_separated_xpd_15 = QPSK_s_AWGN_separated;
% save qpsk_pe_xpd_5_pi_10.mat QPSK_s_AWGN_unseparated_xpd_5 QPSK_s_AWGN_separated_xpd_5 SNR;
% save qpsk_pe_xpd_10_pi_10.mat QPSK_s_AWGN_unseparated_xpd_10 QPSK_s_AWGN_separated_xpd_10;
% save qpsk_pe_xpd_15_pi_10.mat QPSK_s_AWGN_unseparated_xpd_15 QPSK_s_AWGN_separated_xpd_15;





































% 
% 
% 
% clear all;
% close all;
% %误码率曲线的绘制
% SNR=1:1:20;%信噪比变化范围
% snr=10.^(SNR/10);%将信噪比转化成直角坐标
% N=2000000;%仿真点数
% M=4;%进制数
% x=randi([0,1],1,N);%产生随机信号
% y=pskmod(x,M);%调用matlab自带的psk调制函数
% y1=real(y);
% y2=imag(y);
% for i=1:length(SNR)
%     N0=1/2/snr(i);%计算噪声功率
%     N0_dB=10*log10(N0);%将噪声功率转换为dBW
%     ni=wgn(1,N,N0_dB);%产生高斯噪声
%     h=raylrnd(1/sqrt(2),1,N);%产生瑞利信号
% 
%     yAn=y+ni;%通过高斯信道
%     yA=pskdemod(yAn,M);%调用matlab自带的psk解调函数
%     bit_A=length(find(x~=yA));%统计错误比特数 
%     QPSK_s_AWGN(i)=bit_A/N;%计算误码率
%     
%     yRn=y.*h+ni;%通过瑞利信道
%     yR=pskdemod(yRn,M);%调用matlab自带的psk解调函数
%     bit_R=length(find(x~=yR));%统计错误比特数
%     QPSK_s_Ray(i)=bit_R/N;%计算误码率 
% end
% QPSK_t_AWGN=1/2*erfc(sqrt(snr));%AWGN信道下QPSK理论误码率
% QPSK_t_Ray=1/2*(1-sqrt((snr)./(1+snr)));%Rayleigh信道下QPSK理论误码率
% 
% %绘制图形
% figure;
% semilogy(SNR,QPSK_s_AWGN,'-k*');hold on;
% semilogy(SNR,QPSK_t_AWGN,'-go');hold on;
% semilogy(SNR,QPSK_s_Ray,'-b*');hold on
% semilogy(SNR,QPSK_t_Ray,'-ro');grid on;
% axis([-1,20,10^-4,1]);
% legend('AWGN仿真','AWGN理论','瑞利仿真','瑞利理论');
% title('QPSK误码性能分析');
% xlabel('信噪比（dB）');ylabel('BER');
% 
