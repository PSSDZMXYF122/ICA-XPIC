clear all; close all;

load qpsk_pe_xpd_5_pi_10.mat;
load qpsk_pe_xpd_10_pi_10.mat;
load qpsk_pe_xpd_15_pi_10.mat;


SNR=0:1:20;%信噪比变化范围
snr=10.^(SNR/10);%将信噪比转化成直角坐标
QPSK_ideal=1/2*erfc(sqrt(snr));%AWGN信道下QPSK理论误码率



figure(1);
semilogy(SNR,QPSK_s_AWGN_unseparated_xpd_5,'-bx',LineWidth = 2);hold on;

semilogy(SNR,QPSK_s_AWGN_separated_xpd_5,'-b*',LineWidth = 2);hold on;

semilogy(SNR,QPSK_s_AWGN_unseparated_xpd_10,'-gx',LineWidth = 2);hold on;

semilogy(SNR,QPSK_s_AWGN_separated_xpd_10,'-g*',LineWidth = 2);hold on;

semilogy(SNR,QPSK_s_AWGN_unseparated_xpd_15,'-rx',LineWidth = 2);hold on;

semilogy(SNR,QPSK_s_AWGN_separated_xpd_15,'-r*',LineWidth = 2);hold on;

semilogy(SNR,QPSK_ideal,'-ko',LineWidth = 2);hold on;

grid on;grid minor;
axis([-1,20,10^-7,1]);
legend('XPD=5分离前仿真','XPD=5分离后仿真','XPD=10分离前仿真','XPD=10分离后仿真','XPD=15分离前仿真','XPD=15分离后仿真');
% 
legend('XPD=5分离前仿真','XPD=5分离后仿真','XPD=10分离前仿真','XPD=10分离后仿真','XPD=15分离前仿真','XPD=15分离后仿真','理想QPSK误码率');
title('CFastICA算法对QPSK误码信号性能分析');
xlabel('信噪比（dB）');ylabel('误码率 Pe');


