clear;close all;clc;

M = 64; % QAM调制阶数
Ns = 20000; % 符号数
sym = randi([0 M-1],1,Ns); % 500个符号
%生成64QAM基带信号，格雷码映射，画出星座图

Nb = Ns*log2(M); % 比特数
sbit = randi([0 1],Nb,1);
% figure;
sbas = qammod(sbit,M,'InputType','bit','PlotConstellation',true);
Fs = 1;
Rb = 20e6;
RB = Rb/log2(M);
EbN0=0:20;
ber_theo = berawgn(EbN0,'qam',M);
Eb=sum(abs(sbas(:)).^2)/Nb; % 比特能量
N0 = Eb./10.^(EbN0/10);
for i=1:length(EbN0)
    rng('shuffle');
    sigma_n2=N0(i); % 噪声功率
    noise = sqrt(sigma_n2/2)* (randn(size(sbas)) + 1i*randn(size(sbas))); % 生成复高斯噪声
    sbas_n = sbas+noise;
    sdbit = qamdemod(sbas_n,M,'OutputType','bit');
    bit_err_rate(i) = sum(sbit~=sdbit)/Nb;
end
figure;
semilogy(EbN0,ber_theo,'o-r');hold on;grid on;grid minor;
semilogy(EbN0,bit_err_rate,'s-b');
legend('64QAM理论误码率','64QAM实际误码率');
xlabel('Eb/N0(dB)');ylabel('BER');

EbN0=[0 5 15 25 35];
N0 = Eb./10.^(EbN0/10);
for i=1:length(EbN0)
    rng('shuffle');
    sigma_n2=N0(i);
    noise = sqrt(sigma_n2/2)* (randn(size(sbas)) + 1i*randn(size(sbas))); % 生成复高斯噪声
    sbas_n = sbas+noise;
    scatterplot(sbas_n);title(['Eb/N0 = ' num2str(EbN0(i)) ' dB']);
end
