function [W,Q,iteration_num,initial_w,W_unchanged] = fastica_achieve1(mixed_data_noised) 
% B 和 Q 分别是分离矩阵以及白化矩阵

% mixed_data_noised 为算法的输入

MixedS = mixed_data_noised;     % 得到混合信号的矩阵      

%%下面程序是使用FastICA方法对信号进行分离的过程
%%%%%%%%%%%%%%%%%%%%%%%%%%  标准化  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 计算MixedS的均值
MixedS_mean = mean(MixedS.');   % 计算MixedS的均值

MixedS = MixedS - MixedS_mean.';

%%%%%%%%%%%%%%%%%%%%%%%%%%%  白化  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MixedS_cov=cov(MixedS.');                    % cov为求协方差的函数
% MixedS_cov=cov(MixedS(1,:).',MixedS(2,:).'); 

[raw,column]=size(MixedS);
% E_x_real = (real(MixedS)*(real(MixedS))')/column;
% E_x_imag = (imag(MixedS)*(imag(MixedS))')/column;
% cov_x = E_x_real + E_x_imag;

[E1,D]=eig(MixedS_cov);                      % 对信号矩阵的协方差函数进行特征值分解
% Q=inv(sqrt(D))*(E1)';                        % Q为白化矩阵
Q=(D)^(-0.5)*(E1).'; 
% Q_1 = (E1)';
% Q1_x = Q_1*MixedS;
% 
% E_Q1_x_real = (real(Q1_x)*(real(Q1_x))')/column;
% E_Q1_x_imag = (imag(Q1_x)*(imag(Q1_x))')/column;
% 
% E_x_realimag = (real(Q1_x)*(imag(Q1_x))')/column;
% 
% pcov_Q1_x = E_Q1_x_real - E_Q1_x_imag + 2 * E_x_realimag;
% 
% %% 使用SVD奇异值分解；另一种白化方式
% 
% [U,D,V]=svd(pcov_Q1_x);                      % 对信号矩阵的协方差函数进行特征值分解
% % [U2,D2,V2]=eig(pcov_Q1_x)
% Q_2 = U;
% 
% 
% Q = Q_2'*Q_1;

% Q=inv(sqrt(D))*U';                        % Q为白化矩阵
MixedS_white=Q*MixedS;                      % MixedS_white为白化后的信号矩阵
IsI=cov(MixedS_white.');                     % IsI应为单位阵    



% MixedS_white= Q * MixedS;                      % MixedS_white为白化后的信号矩阵     这里的Q有两个注意点，一个是左右乘注意维度，另一个是复数时要乘以其共轭
% IsI=cov(MixedS_white.');                     % IsI应为单位阵            
% (MixedS_white') * MixedS_white/5000

%%%%%%%%%%%%%%%%%%%%%%%%　FASTICA算法  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = MixedS_white;                            % 以下算法将对X即白化后的数据进行操作
[VariableNum,~]=size(X);
numofIC = VariableNum;                       % 在此应用中，独立元个数等于变量个数
W=zeros(numofIC,VariableNum);              % 初始化列向量w的寄存矩阵,B=[b1  b2  ...   bd]
% 
% for r=1:numofIC
% %     maxIterationsNum=5000;               % 设置最大迭代次数（即对于每个独立分量而言迭代均不超过此次数）
% %     b=rand(numofIC,1);
%     b = rand(2,1) + sqrt(-1)*rand(2,1);
%     b=b/norm(b);                         % 对b标准化 norm(b):向量元素平方和开根号(向量范数)
%     iteration_num = 1;
%     while true  
%         b_Old = b;                          
%         w_x = b'*X;
% %         g= tanh(w_x*w_x');
% %         g= 0.5*(0.1+w_x*w_x')^(-0.5);        %对目标函数中的优化函数g()进行的设定。
%         g= (0.1+w_x*w_x')^(-1);
% %         dg = 1 - (tanh(w_x*w_x')).^2;
%         dg = -(0.1+w_x*w_x')^(-2);
% %         dg = -0.25*(0.1+w_x*w_x')^(-1.5);
%         b = (X*w_x')*g - (g+w_x*w_x'*dg)*b;          %迭代计算b
%         b=b-B*B'*b;                        % 对b正交化
%         b=b/norm(b);         
%         if abs(abs(b'*b_Old)-1)<1e-8        % 判断矩阵b收敛的依据，如果收敛，则
%              B(:,r)=b;      % 保存所得向量b
%              disp(iteration_num);               
%              break;
%         end
% %         if abs(abs(b ./ b_Old)-1)<1e-11        % 判断矩阵b收敛的依据，如果收敛，则
% %              B(:,r)=b;                     % 保存所得向量b
% %              break;
% %          end
%         iteration_num = iteration_num + 1;      
%     end
% % %      B11_real=abs(real(B(1,1)));
% % %      B11_imag=abs(imag(B(1,1)));
% % %      B(1,1)=B11_real+B11_imag*sqrt(-1);
% % % 
% % %      B12_real=-abs(real(B(1,2)));
% % %      B12_imag=abs(imag(B(1,2)));
% % %      B(1,2)=B12_real+B12_imag*sqrt(-1);
% % % 
% % %      B21_real=abs(real(B(2,1)));
% % %      B21_imag=abs(imag(B(2,1)));
% % %      B(2,1)=B21_real+B21_imag*sqrt(-1);
% % % 
% % %      B22_real=abs(real(B(2,2)));
% % %      B22_imag=abs(imag(B(2,2)));
% % %      B(2,2)=B22_real+B22_imag*sqrt(-1);
% % 
% end

% maxIterationsNum=5000;               % 设置最大迭代次数（即对于每个独立分量而言迭代均不超过此次数）


%% 估计一个独立成分
% w=rand(numofIC,1) + sqrt(-1)*rand(numofIC,1) ;
% b = rand(3,1) + sqrt(-1)*rand(3,1);     %随机生成一个列向量
% b = [0; 0];  %迭代不出来
% w = rand(3,1);    %测试用的初始值
% w = [0.1;0.1];
% w = [0.9;0.1];
% w=[1;10^(-10/20)*exp(sqrt(-1)*pi/4)];
% w=[0.3;10^(-10/20)];
w = [1;1];
% w = [0.1;0.2;0.3];

% w=rand(numofIC,1);
w=w/norm(w);        % 对b标准化 norm(b):向量元素平方和开根号(向量范数)
initial_w = w;

for r=1:numofIC

w = [1;1];

w=w/norm(w);        % 对b标准化 norm(b):向量元素平方和开根号(向量范数)
initial_w = w;


iteration_num=1;
   

while true  
    w_old = w;                          
    w_x = w'*X;
    w_x_square = w_x.*conj(w_x);

%     g= tanh(t);     %对目标函数中的优化函数G()进行的设定。

%% G2
%     g= (0.1+(w_x_square)).^(-1);
%     dg = -(0.1+(w_x_square)).^(-2);

%     g= 0.5*(0.1+w_x.*conj(w_x)).^(-0.5);
%     dg = -0.25*(0.1+w_x.*conj(w_x)).^(-1.5);
%% Tukey函数
%     for i=1:1000
%     if aaa(i)<4
      g = 0.5*(1- w_x_square/(2*2)).^2;
      dg = -1/(2*2)*(1-w_x_square/(2*2));
%     else    
%         g = 0;
%         dg = 0;
%     end
%     end
%     g= 0.5*(1-w_x*w_x'/4)^(2);
%     dg = -0.25*(1-(w_x*w_x')/4);

%       g = 0.5*(1- w_x.*conj(w_x)/(2*2)).^2;
%       dg = -1/(2*2)*(1 - w_x.*conj(w_x)/(2*2));
%       g = 0.5*(1- w_x*(w_x)'/(2*2))^2;
%       dg = -1/(2*2)*(1 - w_x*(w_x)'/(2*2));

%% G1        
%       g = 0.5*(0.1 +  w_x.*conj(w_x)).^(-1/2);
%       dg = -0.25*(0.1 +  w_x.*conj(w_x)).^(-3/2);
%       g = 0.5*(0.1 +  w_x_square).^(-1/2);
%       dg = -0.25*(0.1 +  w_x_square).^(-3/2);


% 
%     g = tanh(w_x_square);
%     dg = 1 - (tanh(w_x_square))^2;

%     g = (w_x.*conj(w_x)).*exp(-(w_x.*conj(w_x)).^2/2);
%     dg = (1 + w_x.*conj(w_x)).*exp(-(w_x.*conj(w_x)).^2/2);

%     dg = 1 - (tanh(t)).^2;      %优化函数G()的微分表达式
%     temp = sum((X.*conj(w_x)).');
%      b = (temp.'*g)/5000 - ((g+(w_x*w_x')*dg)/5000)*b;
    expect_1 = X.*conj(w_x).*g;
    expect_2 = mean(g) + mean(w_x_square.*dg);
    w = (mean((expect_1).')).' - expect_2*w;

% %     b = (X*w_x')/column*mean(g) - (mean(g) - w_x*w_x'/column*mean(dg))*b;
% %       b = (X*w_x'*mean(g.*conj(g))) - mean(g+w_x.*conj(w_x).*dg)*b;
% %       b = (X*w_x'*g) - mean( g + w_x*w_x'*dg)*b;
% %     b = (X*g) - mean(dg)*b;           %迭代计算b
% %     b=b-B*B'*b;                        % 对b正交化
    w=w-W*W'*w;                        % 对b正交化    
% 
%     b = b-W*W'*b;
    w=w/norm(w);    

    if abs(abs(w'*w_old)-1)<1e-9       % 判断矩阵b收敛的依据，如果收敛，则
%     if abs(abs(b-[0.7;-0.7])-1)<0.3       % 判断矩阵b收敛的依据，如果收敛，则
%     if abs(abs(b ./ b_Old)-1)<1e-8        % 判断矩阵b收敛的依据，如果收敛，则
         W(:,r)=w;      % 保存所得向量b
% 
% %          B(:,1)=b;                     % 保存所得向量b                 
%          disp(iteration_num);               
         break;
%     end
% %         if abs(abs(b ./ b_Old)-1)<1e-11        % 判断矩阵b收敛的依据，如果收敛，则
% %              B(:,r)=b;                     % 保存所得向量b
% %              break;
% %          end
     
    end
    iteration_num = iteration_num + 1;    

    end

end

W_unchanged = W;

% B_Primary = B;
% B11_real=abs(real(B_Primary(1,1)));
% B11_imag=0;%-abs(imag(B_Primary(1,1)));  虚部置零
% B(1,1)=B11_real+B11_imag*sqrt(-1);
% 
% B21_real=-abs(real(B_Primary(2,1)));
% B21_imag=0;%abs(imag(B_Primary(2,1)));
% B(2,1)=B21_real+B21_imag*sqrt(-1);
% 
% B12_real=abs(real(B_Primary(1,2)));
% B12_imag=0;%-abs(imag(B_Primary(1,2)));
% B(1,2)=B12_real+B12_imag*sqrt(-1);
% 
% B22_real=abs(real(B_Primary(2,2)));
% B22_imag=0;%abs(imag(B_Primary(2,2)));
% B(2,2)=B22_real+B22_imag*sqrt(-1);




% B_Primary = B;
% B11_real = abs(real(B_Primary(1,1)));
% B11_imag = abs(imag(B_Primary(1,1)));  
% % B11_imag = 0;
% B(1,1) = B11_real+B11_imag*sqrt(-1);
% 
% B21_real = -abs(real(B_Primary(2,1)));
% B21_imag = -abs(imag(B_Primary(2,1)));
% % B21_imag = -B11_imag;
% B(2,1) = B21_real+B21_imag*sqrt(-1);
% 
% B12_real = abs(real(B_Primary(1,2)));
% % B12_imag = abs(imag(B_Primary(1,2)));
% B12_imag = -B11_imag;
% B(1,2) = B12_real+B12_imag*sqrt(-1);
% 
% B22_real = abs(real(B_Primary(2,2)));
% B22_imag = -abs(imag(B_Primary(2,2)));
% % B22_imag = B11_imag;
% B(2,2) = B22_real+B22_imag*sqrt(-1);




%% 扭转相位 0-pi/2 左闭右开
W_Primary = W;
W11_real = abs(real(W_Primary(1,1)));
W11_imag = -abs(imag(W_Primary(1,1))); 
% W11_imag = 0;
W(1,1) = W11_real+W11_imag*sqrt(-1);

W21_real = abs(real(W_Primary(2,1)));
W21_imag = abs(imag(W_Primary(2,1)));
% W21_imag = 0;
W(2,1) = W21_real+W21_imag*sqrt(-1);

W12_real = -abs(real(W_Primary(2,1)));
W12_imag = abs(imag(W_Primary(2,1)));
% W12_imag = -W11_imag;
% W12_imag = 0;
W(1,2) = W12_real+W12_imag*sqrt(-1);

W22_real = abs(real(W_Primary(1,1)));
W22_imag = abs(imag(W_Primary(1,1)));
% W22_imag = 0;
W(2,2) = W22_real+W22_imag*sqrt(-1);



% % 
% % 
% %% 扭转相位 pi/2-pi 左闭右开
% W_Primary = W;
% W11_real = -abs(real(W_Primary(1,1)));
% W11_imag = -abs(imag(W_Primary(1,1))); 
% % W11_imag = 0;
% W(1,1) = W11_real+W11_imag*sqrt(-1);
% 
% W21_real = abs(real(W_Primary(2,1)));
% W21_imag = -abs(imag(W_Primary(2,1)));
% % W21_imag = 0;
% W(2,1) = W21_real+W21_imag*sqrt(-1);
% 
% W12_real = -abs(real(W_Primary(2,1)));
% W12_imag = -abs(imag(W_Primary(2,1)));
% % W12_imag = -W11_imag;
% % W12_imag = 0;
% W(1,2) = W12_real+W12_imag*sqrt(-1);
% 
% W22_real = -abs(real(W_Primary(1,1)));
% W22_imag = abs(imag(W_Primary(1,1)));
% % W22_imag = 0;
% W(2,2) = W22_real+W22_imag*sqrt(-1);
% % W=-W;
% 
% 



% %% 扭转相位 pi-pi*3/2 左闭右开  行不通
% W_Primary = W;
% W11_real = abs(real(W_Primary(1,1)));
% W11_imag = -abs(imag(W_Primary(1,1))); 
% % W11_imag = 0;
% W(1,1) = W11_real+W11_imag*sqrt(-1);
% 
% W21_real = -abs(real(W_Primary(2,1)));
% W21_imag = -abs(imag(W_Primary(2,1)));
% % W21_imag = 0;
% W(2,1) = W21_real+W21_imag*sqrt(-1);
% 
% W12_real = abs(real(W_Primary(2,1)));
% W12_imag = -abs(imag(W_Primary(2,1)));
% % W12_imag = -W11_imag;
% % W12_imag = 0;
% W(1,2) = W12_real+W12_imag*sqrt(-1);
% 
% W22_real = abs(real(W_Primary(1,1)));
% W22_imag = abs(imag(W_Primary(1,1)));
% % W22_imag = 0;
% W(2,2) = W22_real+W22_imag*sqrt(-1);
% % W=-W;

% %% 扭转相位 pi*3/2-pi*2 左闭右开
% W_Primary = W;
% W11_real = abs(real(W_Primary(1,1)));
% W11_imag = abs(imag(W_Primary(1,1))); 
% % W11_imag = 0;
% W(1,1) = W11_real+W11_imag*sqrt(-1);
% 
% W21_real = abs(real(W_Primary(2,1)));
% W21_imag = -abs(imag(W_Primary(2,1)));
% % W21_imag = 0;
% W(2,1) = W21_real+W21_imag*sqrt(-1);
% 
% W12_real = -abs(real(W_Primary(2,1)));
% W12_imag = -abs(imag(W_Primary(2,1)));
% % W12_imag = -W11_imag;
% % W12_imag = 0;
% W(1,2) = W12_real+W12_imag*sqrt(-1);
% 
% W22_real = abs(real(W_Primary(1,1)));
% W22_imag = -abs(imag(W_Primary(1,1)));
% % W22_imag = 0;
% W(2,2) = W22_real+W22_imag*sqrt(-1);
% % W=-W;

% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%　FASTICA算法  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X=MixedS_white;                            % 以下算法将对X进行操作
% [VariableNum,~]=size(X);
% numofIC=VariableNum;                       % 在此应用中，独立元个数等于变量个数
% B=zeros(numofIC,VariableNum);              % 初始化列向量w的寄存矩阵,B=[b1  b2  ...   bd]
% maxIterationsNum=500;               % 设置最大迭代次数（即对于每个独立分量而言迭代均不超过此次数）
% b=rand(numofIC,1);
% b=b/norm(b);   % 对b标准化 norm(b):向量元素平方和开根号(向量范数)
% i = 0;
% while i<maxIterationsNum+1  
%     i = i + 1;
%     b_Old = b;                          
%     t=X'*b;
% %         g=(exp(2*t)-1)./(exp(2*t)+1);     %函数g的选取这个函数挺重要
% %         dg=4*exp(2*t)./(exp(2*t)+1).^2;   %dg表示对函数g的微分后的表达式 
%     g= tanh(t);
%     dg = 1 - (tanh(t)).^2;
% %         g= t.*exp(-0.5*t.^2);
% %         dg = (1-t.^2).*exp(-0.5*t.^2);
%     b = mean(X*g) - mean(dg)*b;           %迭代计算b
%     b=b-B*B'*b;                        % 对b正交化
%     b=b/norm(b);  
% %         if(r==2)
% %             if(real(b(1))>0)
% %                 continue;
% %             end
% %         end
%     if abs(abs(b'*b_Old)-1)<1e-9        % 判断矩阵b收敛的依据，如果收敛，则
%          B(:,1)=b;                     % 保存所得向量b
%          break;
%     end
% 
% end
% 
% 
% B = [B(1,1),-B(1,1);B(2,1),B(2,1)];
