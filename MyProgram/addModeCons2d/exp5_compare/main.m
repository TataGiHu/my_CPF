%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Motion tracking framework based on fusion mode     %%%
%%%                                                        %%%
%%%     Version:2-dimensional                              %%%
%%%     Author:Tata                                        %%%
%%%     Last modified date:2021-10-31                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;

numOfExp = 100;%实验次数(Number of experiments)
EXP_delay = zeros(4,numOfExp);
EXP_RMSE = zeros(4,numOfExp);
EXP_SD = zeros(4,numOfExp);
for exp=1:numOfExp

close all;
%% Parameters Setting
N = 300;%粒子数（Number of particle）
Q = 5;%过程噪声指标（Process noise index）
R = 5;%观测噪声指标（Observation noise index）
Noi_enhance =2;%噪声增强系数（Noise Enhancement factor）
T = 60;%时间序列长度（Length of time series）
st = 4*80/T;%转移步长（Transfer step）
WorldSize = 100; %世界大小(Size of World)
numOfStage = 6;%阶段数（Number of stage）
%numOfMethod = 4;%用于实现的方法数（Number of methods）
%Method Gate
run_SIR = 1;
run_APF = 1;
run_MPF = 1;
run_CPF = 1;%Our Method
%%  Valuation Index
%Error evaluation
% RMSE_SIR = zeros(numOfExp,T);
% RMSE_APF = zeros(numOfExp,T);
% RMSE_MPF = zeros(numOfExp,T);
% RMSE_CPF = zeros(numOfExp,T);
%Standerd deviation dvaluation
% SD_SIR = zeros(numOfExp,T);
% SD_APF = zeros(numOfExp,T);
% SD_MPF = zeros(numOfExp,T);
% SD_CPF = zeros(numOfExp,T);
%Delay evaluation
delay_SIR = 0;
delay_APF = 0;
delay_MPF = 0;
delay_CPF = 0;

%% Real State of System
X = zeros(2, T);    %系统的真实状态（State of System）
Z = zeros(2, T);    %系统的观测状态（Observation of System）
X(:, 1) = [50; 20];     %初始系统状态(State Initialization)
Z(:, 1) = [50; 20] + wgn(2, 1, 10*log10(R));    %初始系统的观测状态(Observation Initialization)
for k=2:T
   if k<=round(T/numOfStage)
        X(1, k) = X(1, k-1) + st * 1.5 + Noi_enhance*wgn(1, 1, 10*log10(Q));     %状态方程
        X(2, k) = X(2, k-1) + st * 0 + Noi_enhance*wgn(1, 1, 10*log10(Q));
        Z(:, k) = X(:, k) + wgn(2, 1, 10*log10(R));     %观测方程    
    elseif (round(T/numOfStage)<k) && (k<=round(2*T/numOfStage))
        X(1, k) = X(1, k-1) + st *  1.5 + Noi_enhance*wgn(1, 1, 10*log10(Q));     %状态方程
        X(2, k) = X(2, k-1) + st * 1.5 + Noi_enhance*wgn(1, 1, 10*log10(Q));
        Z(:, k) = X(:, k) + wgn(2, 1, 10*log10(R));     %观测方程 
    elseif (round(2*T/numOfStage)<k) && (k<=round(3*T/numOfStage))
        X(1, k) = X(1, k-1) + st *  1.5 + Noi_enhance*wgn(1, 1, 10*log10(Q));     %状态方程
        X(2, k) = X(2, k-1) + st * -1.5 + Noi_enhance*wgn(1, 1, 10*log10(Q));
        Z(:, k) = X(:, k) + wgn(2, 1, 10*log10(R));     %观测方程 
    elseif (round(3*T/numOfStage)<k) && (k<=round(4*T/numOfStage))
        X(1, k) = X(1, k-1) + st * 1.5 + Noi_enhance*wgn(1, 1, 10*log10(Q));     %状态方程
        X(2, k) = X(2, k-1) + st * 0 + Noi_enhance*wgn(1, 1, 10*log10(Q));
        Z(:, k) = X(:, k) + wgn(2, 1, 10*log10(R));     %观测方程
    elseif (round(4*T/numOfStage)<k) && (k<=round(5*T/numOfStage))
        X(1, k) = X(1, k-1) + st *  1.5 + Noi_enhance*wgn(1, 1, 10*log10(Q));     %状态方程
        X(2, k) = X(2, k-1) + st * -1.5 + Noi_enhance*wgn(1, 1, 10*log10(Q));
        Z(:, k) = X(:, k) + wgn(2, 1, 10*log10(R));     %观测方程 
    else
        X(1, k) = X(1, k-1) + st *  1.5 + Noi_enhance*wgn(1, 1, 10*log10(Q));     %状态方程
        X(2, k) = X(2, k-1) + st * 1.5 + Noi_enhance*wgn(1, 1, 10*log10(Q));
        Z(:, k) = X(:, k) + wgn(2, 1, 10*log10(R));     %观测方程
    end 
end

%% SIR-PF Run
if(run_SIR ==1)
   disp('-- SIR Begin --');
   tic;
   [RMSE_SIR, SD_SIR, PCenter_SIR] = SIR_function(N,Q,R,T,st,WorldSize,X,Z);
   SIR_time = toc;
   delay_SIR = SIR_time/T;
   disp('-- SIR Already --');
end

%% APF Run
if(run_APF == 1)
   disp('-- APF Begin --');
   tic;
   [RMSE_APF, SD_APF, PCenter_APF] = APF_function(N,Q,R,T,st,WorldSize,X,Z);
   APF_time = toc;
   delay_APF = APF_time/T;
   disp('-- APF Already --');
end

%% 
if(run_MPF == 1)
   disp('-- MPF Begin --');
   tic;
   [RMSE_MPF, SD_MPF, PCenter_MPF] = MPF_function(N,Q,R,T,st,WorldSize,X,Z);
   MPF_time = toc;
   delay_MPF = MPF_time/T;
   disp('-- MPF Already --');
end

%%
if(run_CPF == 1)
   disp('-- CPF Begin --');
   tic;
   [RMSE_CPF, SD_CPF, PCenter_CPF] = CPF_function(N,Q,R,T,st,WorldSize,X,Z);
   CPF_time = toc;
   delay_CPF = CPF_time/T;
   disp('-- CPF Already --');
end

%% Visualization
figure();
%Trajectory
plot(X(1,:),X(2,:),'k-','LineWidth',1.5);hold on;
plot(PCenter_SIR(1,:),PCenter_SIR(2,:));hold on;
plot(PCenter_APF(1,:),PCenter_APF(2,:));hold on;
plot(PCenter_MPF(1,:),PCenter_MPF(2,:));hold on;
plot(PCenter_CPF(1,:),PCenter_CPF(2,:));hold on;
xlabel('X');
ylabel('Y');
title('Trajectory','FontWeight','bold');
legend('Real Transfer','SIR','APF','MPF','CPF');

figure();
%RMSE
plot(2:1:T,RMSE_SIR(2:T),'.-.');hold on;
plot(2:1:T,RMSE_APF(2:T),'.-.');hold on;
plot(2:1:T,RMSE_MPF(2:T),'.-.');hold on;
plot(2:1:T,RMSE_CPF(2:T),'.-.');hold on;
xlabel('Time Series');
ylabel('Value');
title('RMSE','FontWeight','bold');
legend('SIR','APF','MPF','CPF');

figure();
%Standard deviation
plot(2:1:T,SD_SIR(2:T),'.-.');hold on;
plot(2:1:T,SD_APF(2:T),'.-.');hold on;
plot(2:1:T,SD_MPF(2:T),'.-.');hold on;
plot(2:1:T,SD_CPF(2:T),'.-.');hold on;
xlabel('Time Series');
ylabel('Value');
title('Standard Deviation','FontWeight','bold');
legend('SIR','APF','MPF','CPF');

%Delay
fprintf('The average delay of SIR: %f\t\t',delay_SIR);
fprintf('The RMSE of SIR: %f\t\t',sum(RMSE_SIR(2:T))/(T-1));
fprintf('The standard deviation of SIR: %f\n',sum(SD_SIR(2:T))/(T-1));
fprintf('The average delay of APF: %f\t\t',delay_APF);
fprintf('The RMSE of APF: %f\t\t',sum(RMSE_APF(2:T))/(T-1));
fprintf('The standard deviation of AIR: %f\n',sum(SD_APF(2:T))/(T-1));
fprintf('The average delay of MPF: %f\t\t',delay_MPF);
fprintf('The RMSE of MPF: %f\t\t',sum(RMSE_MPF(2:T))/(T-1));
fprintf('The standard deviation of MPF: %f\n',sum(SD_MPF(2:T))/(T-1));
fprintf('The average delay of CPF: %f\t\t',delay_CPF);
fprintf('The RMSE of CPF: %f\t\t',sum(RMSE_CPF(2:T))/(T-1));
fprintf('The standard deviation of CPF: %f\n',sum(SD_CPF(2:T))/(T-1));


EXP_delay(1,exp)=delay_SIR;
EXP_delay(2,exp)=delay_APF;
EXP_delay(3,exp)=delay_MPF;
EXP_delay(4,exp)=delay_CPF;

EXP_RMSE(1,exp)=sum(RMSE_SIR(2:T))/(T-1);
EXP_RMSE(2,exp)=sum(RMSE_APF(2:T))/(T-1);
EXP_RMSE(3,exp)=sum(RMSE_MPF(2:T))/(T-1);
EXP_RMSE(4,exp)=sum(RMSE_CPF(2:T))/(T-1);

EXP_SD(1,exp)=sum(SD_SIR(2:T))/(T-1);
EXP_SD(2,exp)=sum(SD_APF(2:T))/(T-1);
EXP_SD(3,exp)=sum(SD_MPF(2:T))/(T-1);
EXP_SD(4,exp)=sum(SD_CPF(2:T))/(T-1);


exp
disp('-------------');
end
mean_EXP_delay = mean(EXP_delay,2);
mean_RMSE_delay = mean(EXP_RMSE,2);
mean_SD_delay = mean(EXP_SD,2);
disp('*********************************************************************************************************')
disp('--------------------------------Results of 100 experiments----------------------------------------------')
fprintf('The average delay of SIR: %f\t\t',mean_EXP_delay(1));
fprintf('The RMSE of SIR: %f\t\t',mean_RMSE_delay(1));
fprintf('The standard deviation of SIR: %f\n',mean_SD_delay(1));
fprintf('The average delay of APF: %f\t\t',mean_EXP_delay(2));
fprintf('The RMSE of APF: %f\t\t',mean_RMSE_delay(2));
fprintf('The standard deviation of APF: %f\n',mean_SD_delay(2));
fprintf('The average delay of MPF: %f\t\t',mean_EXP_delay(3));
fprintf('The RMSE of MPF: %f\t\t',mean_RMSE_delay(3));
fprintf('The standard deviation of MPF: %f\n',mean_SD_delay(3));
fprintf('The average delay of CPF: %f\t\t',mean_EXP_delay(4));
fprintf('The RMSE of CPF: %f\t\t',mean_RMSE_delay(4));
fprintf('The standard deviation of CPF: %f\n',mean_SD_delay(4));
disp('************************************************************************************************************')

EXP_delay=EXP_delay';
EXP_RMSE=EXP_RMSE';
EXP_SD=EXP_SD';

figure();
plot(EXP_delay(:,1),'--');hold on;
plot(EXP_delay(:,2),'--');hold on;
plot(EXP_delay(:,3),'--');hold on;
plot(EXP_delay(:,4),'-','LineWidth',1.2);hold on;
xlabel('Time Series');
ylabel('Value(s)');
title('Delay','FontWeight','bold');
legend('SIR','APF','MPF','CPF(Ours)');
figure();
plot(EXP_RMSE(:,1),'--');hold on;
plot(EXP_RMSE(:,2),'--');hold on;
plot(EXP_RMSE(:,3),'--');hold on;
plot(EXP_RMSE(:,4),'-','LineWidth',1.2);hold on;
xlabel('Time Series');
ylabel('Value');
title('RMSE','FontWeight','bold');
legend('SIR','APF','MPF','CPF(Ours)');
figure();
plot(EXP_SD(:,1),'--');hold on;
plot(EXP_SD(:,2),'--');hold on;
plot(EXP_SD(:,3),'--');hold on;
plot(EXP_SD(:,4),'-','LineWidth',1.2);hold on;
xlabel('Time Series');
ylabel('Value');
title('Standard Deviation','FontWeight','bold');
legend('SIR','APF','MPF','CPF(Ours)');





