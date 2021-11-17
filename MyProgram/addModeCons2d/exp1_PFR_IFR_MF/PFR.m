%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Motion tracking framework based on fusion mode     %%%
%%%                                                        %%%
%%%     Version:2-dimensional                              %%%
%%%     Author:Tata                                        %%%
%%%     Last modified date:2020-12-16                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
close all;
%clear all;

numberOfExp=10;
EXP_frameacc=zeros(1,numberOfExp);
EXP_acc=zeros(1,numberOfExp);
% EXP_framedelay=zeros(1,numberOfExp);
% EXP_recdelay=zeros(1,numberOfExp);
EXP_time=zeros(1,numberOfExp);
for e=1:numberOfExp
close all;
N = 300;   %粒子总数
Q = 5;      %过程噪声
R = 5;      %测量噪声
T = 60;     %测量时间
theta = pi/T;       %旋转角度
distance = 4*80/T;    %每次走的距离
WorldSize = 100;    %世界大小
X = zeros(2, T);    %存储系统状态
Z = zeros(2, T);    %存储系统的观测状态
P = zeros(2, N);    %建立粒子群
PCenter = zeros(2, T);  %所有粒子的中心位置
w = zeros(N, 1);         %每个粒子的权重
err = zeros(1,T);     %误差
X(:, 1) = [5; 20];     %初始系统状态
%wgn（m,n,p）产生一个m行n列的高斯白噪声的矩阵，p以dBW为单位指定输出噪声的强度。
Z(:, 1) = [5; 20] + wgn(2, 1, 10*log10(R));    %初始系统的观测状态
M=[1 0 0 0;%模式矩阵
   1 0 1 0;
   0 1 1 0;
   1 0 0 1;
   0 1 0 1];
Mq = zeros(N,1);%传递模式索引
Mqu = zeros(N,T);%存储模式索引
MPro(:,1) = 0.2*ones(5,1);%初始化模式权重
MNumber = length(MPro);%模式数量
MCum = cumsum(MPro);%模式权重累加
modeCorrect=0;
modeError=0;
timeOfAll=0;
stageNumber=6;%定义阶段数
recSpaceSize=3;%定义识别帧基数
recSpace=zeros(recSpaceSize,1);%定义当前识别基
timeEachFrame=zeros(T,1);%存放每一帧耗时
recChange=zeros(stageNumber,1);%记录识别变更帧
timeDelay=zeros(stageNumber,1);%记录识别变更时延
tempRec=0;%识别变更栈
u=1;%用于recChange的循环赋值

recRight=0;%识别准确数
recCount=0;%识别帧数（未丢失）
recMiss=0;%识别丢失帧
recACC=0;%识别准确率

outTime=0;
XX=zeros(2,8);%存储真实状态的不同转移(虚拟转移)
ZZ=zeros(2,8);%存储观测状态的不同转移（虚拟转移）
EE=zeros(1,8);%存储虚拟误差

recRecord=zeros(T+10,1);
% modeInStage1=zeros(T/6,1);
% modeInStage2=zeros(T/6,1);
% modeInStage3=zeros(T/6,1);
% modeInStage4=zeros(T/6,1);
% modeInStage5=zeros(T/6,1);
% modeInStage6=zeros(T/6,1);
%初始化粒子群
for i = 1 : N
    P(:, i) = [WorldSize*rand; WorldSize*rand];%在worldSize区域内随机生成N个粒子
    dist = norm(P(:, i)-Z(:, 1));     %与测量位置相差的距离，用于估算该粒子的权重
    %由于上面已经随机生成了N个粒子，现在将其与真实的测量值z进行比较，越接近则权重越大，或者说差值越小权重越大
    %这里的权重计算是关于p(z/x)的分布，即观测方程的分布，假设观测噪声满足高斯分布，那么w（i）=p(z/x)
    w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);   %求权重
end
PCenter(:, 1) = sum(P, 2) / N;      %所有粒子的几何中心位置

%%
err(1) = norm(X(:, 1) - PCenter(:, 1));     %粒子几何中心与系统真实状态的误差
% figure(1);
% %利用set(gca,'propertyname','propertyvalue'......)命令可以调整图形的坐标属性。
% %propertyname和property value分别为坐标属性名称（不区分字母大小写）及对应的属性值。
% set(gca,'FontSize',10);
% hold on%使当前及图形保持而不被刷新，准备接受此后将绘制的新曲线
% plot(X(1, 1), X(2, 1), 'r.', 'markersize',30)   %系统状态位置
% %axis[xmin,xmax,ymin,ymax]设定坐标范围，必须满足xmin<xmax，后面同理
% axis([0 100 0 100]);
% plot(P(1, :), P(2, :), 'kx', 'markersize',5);   %绘出各个粒子位置
% plot(PCenter(1, 1), PCenter(2, 1), 'b.', 'markersize',25); %所有粒子的几何中心位置
% %图形标注函数legend(string1,string2,string3...),在当前图中添加图例
%  
% legend('True State', 'Particles', 'The Center of Particles');
% %添加图形标题命令title('string'),在当前坐标系的顶部加一个文本串string,作为该图形的标题
% title('Initial State');
% %使当前轴及图形不再具备不被刷新的性质
% hold off

%%
%开始运动
for k = 2 : T

    %模拟一个弧线运动的状态
    %wgn（m,n,p）产生一个m行n列的高斯白噪声的矩阵，p以dBW为单位指定输出噪声的强度。
%     X(:, k) = X(:, k-1) + distance * [(-cos(k * theta)); sin(k * theta)] + wgn(2, 1, 10*log10(Q));     %状态方程
%     Z(:, k) = X(:, k) + wgn(2, 1, 10*log10(R));     %观测方程 
    if k<=round(T/stageNumber)
        X(1, k) = X(1, k-1) + distance * 1.5 + wgn(1, 1, 10*log10(Q));     %状态方程
        X(2, k) = X(2, k-1) + distance * 0 + wgn(1, 1, 10*log10(Q));
        Z(:, k) = X(:, k) + wgn(2, 1, 10*log10(R));     %观测方程    
    elseif (round(T/stageNumber)<k) && (k<=round(2*T/stageNumber))
        X(1, k) = X(1, k-1) + distance *  1.5 + wgn(1, 1, 10*log10(Q));     %状态方程
        X(2, k) = X(2, k-1) + distance * 1.5 + wgn(1, 1, 10*log10(Q));
        Z(:, k) = X(:, k) + wgn(2, 1, 10*log10(R));     %观测方程 
    elseif (round(2*T/stageNumber)<k) && (k<=round(3*T/stageNumber))
        X(1, k) = X(1, k-1) + distance *  1.5 + wgn(1, 1, 10*log10(Q));     %状态方程
        X(2, k) = X(2, k-1) + distance * -1.5 + wgn(1, 1, 10*log10(Q));
        Z(:, k) = X(:, k) + wgn(2, 1, 10*log10(R));     %观测方程 
    elseif (round(3*T/stageNumber)<k) && (k<=round(4*T/stageNumber))
        X(1, k) = X(1, k-1) + distance * 1.5 + wgn(1, 1, 10*log10(Q));     %状态方程
        X(2, k) = X(2, k-1) + distance * 0 + wgn(1, 1, 10*log10(Q));
        Z(:, k) = X(:, k) + wgn(2, 1, 10*log10(R));     %观测方程
    elseif (round(4*T/stageNumber)<k) && (k<=round(5*T/stageNumber))
        X(1, k) = X(1, k-1) + distance *  1.5 + wgn(1, 1, 10*log10(Q));     %状态方程
        X(2, k) = X(2, k-1) + distance * -1.5 + wgn(1, 1, 10*log10(Q));
        Z(:, k) = X(:, k) + wgn(2, 1, 10*log10(R));     %观测方程 
    else
        X(1, k) = X(1, k-1) + distance *  1.5 + wgn(1, 1, 10*log10(Q));     %状态方程
        X(2, k) = X(2, k-1) + distance * 1.5 + wgn(1, 1, 10*log10(Q));
        Z(:, k) = X(:, k) + wgn(2, 1, 10*log10(R));     %观测方程
    end

    %粒子滤波
    tic;
    %预测
    for i = 1 : N
        %将之前生成的粒子带入到状态方程，进行下一步状态的预测
%         P(:, i) = P(:, i) + distance * [-cos(k * theta); sin(k * theta)] + wgn(2, 1, 10*log10(Q));
       RanNumber = rand;
            for MIndex = 1:MNumber
               if RanNumber <= MCum(MIndex)
                   Mq(i) = MIndex;
                   break;
               end
            end
        thisMode = M(Mq(i),:);
        P(1, i) = P(1, i) + distance * (thisMode(1)*1.5+thisMode(2)*(-1.5)) + wgn(1, 1, 10*log10(Q));
        P(2, i) = P(2, i) + distance * (thisMode(3)*1.5+thisMode(4)*(-1.5)) + wgn(1, 1, 10*log10(Q));
        dist = norm(P(:, i)-Z(:, k));     %与测量位置相差的距离，用于估算该粒子的权重
        %由于上面已经随机生成了N个粒子，现在将其与真实的测量值z进行比较，越接近则权重越大，或者说差值越小权重越大
        %这里的权重计算是关于p(z/x)的分布，即观测方程的分布，假设观测噪声满足高斯分布，那么w（i）=p(z/x)
        w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);   %求权重
    end
%归一化权重
    wsum = sum(w);
    for i = 1 : N
        w(i) = w(i) / wsum;
    end

    %重采样（更新）
    for i = 1 : N
        wmax = 2 * max(w) * rand;  %另一种重采样规则,获取权重当中的最大值与均匀分布在0-1之间数的乘积再乘以2作为后面判断选出大权重粒子的依据
        %randi()函数生成均匀分布的伪随机整数，范围为imin--imax，如果没指定imin，则默认为1
        %r = randi([imin,imax],...)返回一个在[imin,imax]范围内的伪随机整数
        index = randi(N, 1);
        %在这里随机产生N个粒子当中的某个粒子，然后选择该粒子的权值，判断是否是大的粒子权重，是的话就把该粒子选出来
        while(wmax > w(index))
            %使vmax递减，不然的话单个粒子的权重是不可能大于vmax的值的
            wmax = wmax - w(index);
            index = index + 1;
            %如果从某个随机的粒子开始找，没找到从该粒子开始到最后的粒子权值总和大于vmax，那么就重新从第一个粒子开始查找
            if index > N
                index = 1;
            end          
        end
        P(:, i) = P(:, index);     %从上面的index中得到新粒子
%         Mq(i) = Mq(index); %选出粒子对应的模式索引
    end

%     table = tabulate(Mq);
%     [maxCount,idx] = max(table(:,2));
%     bestMode=table(idx);
%     Mqu(:,k)=Mq;
%     MPro(1,1)=length(find(Mq==1));
%     MPro(2,1)=length(find(Mq==2));
%     MPro(3,1)=length(find(Mq==3));
%     MPro(4,1)=length(find(Mq==4));
%     MPro(5,1)=length(find(Mq==5));
%     MPro = MPro/sum(MPro);
%     for i=1:MNumber
%        if MPro(i) <= 0.05
%           MPro= 0.2*ones(5,1);
%        end
%     end
%     MCum = cumsum(MPro);
    

    %sum(X,2)表示把X按行求和;如果是sum(X),那就是按列求和
    PCenter(:, k) = sum(P, 2) / N;      %所有粒子的几何中心位置

    %计算误差
    err(k) = norm(X(:, k) - PCenter(:, k));     %粒子几何中心与系统真实状态的误差
    
    %计算不同转移的误差
    XX(1,1)= X(1, k-1) + distance * 1.5 + wgn(1, 1, 10*log10(Q));     %状态方程
    XX(2,1) = X(2, k-1) + distance * 0 + wgn(1, 1, 10*log10(Q));   
    ZZ(:,1) = XX(:, 1) + wgn(2, 1, 10*log10(R));     %观测方程
    EE(1)=norm(Z(:,k)-ZZ(:,1));
    
    XX(1,2)= X(1, k-1) + distance * 1.5 + wgn(1, 1, 10*log10(Q));     %状态方程
    XX(2,2) = X(2, k-1) + distance * 1.5 + wgn(1, 1, 10*log10(Q));
    ZZ(:,2) = XX(:, 2) + wgn(2, 1, 10*log10(R));     %观测方程
    EE(2)=norm(Z(:,k)-ZZ(:,2));
    
    XX(1,3)= X(1, k-1) + distance * -1.5 + wgn(1, 1, 10*log10(Q));     %状态方程
    XX(2,3) = X(2, k-1) + distance * 1.5 + wgn(1, 1, 10*log10(Q));
    ZZ(:,3) = XX(:, 3) + wgn(2, 1, 10*log10(R));     %观测方程
    EE(3)=norm(Z(:,k)-ZZ(:,3));
    
    XX(1,4)= X(1, k-1) + distance * 1.5 + wgn(1, 1, 10*log10(Q));     %状态方程
    XX(2,4) = X(2, k-1) + distance * -1.5 + wgn(1, 1, 10*log10(Q));
    ZZ(:,4) = XX(:, 4) + wgn(2, 1, 10*log10(R));     %观测方程
    EE(4)=norm(Z(:,k)-ZZ(:,4));
    
    XX(1,5)= X(1, k-1) + distance * -1.5 + wgn(1, 1, 10*log10(Q));     %状态方程
    XX(2,5) = X(2, k-1) + distance * -1.5 + wgn(1, 1, 10*log10(Q));
    ZZ(:,5) = XX(:, 5) + wgn(2, 1, 10*log10(R));     %观测方程 
    EE(5)=norm(Z(:,k)-ZZ(:,5));
    
    

%     
%     
%      
%     
    %获取最佳模式
    bestMode = 1;
%     temp = err(k)-EE(1);
%     for i=2:MNumber
%         if err(k)-EE(i) < temp
%             temp = err(k)-EE(i);
%             bestMode = i;
%         end
%     end
    temp=EE(1);
    for i=2:MNumber
        if EE(i) < temp
            temp = EE(i);
            bestMode = i;
        end
    end
    outTime=outTime+toc;
%     timeOfAll=toc+timeOfAll;
%     timeEachFrame(k)=toc;
    %统计模式准确帧
    if (k<=round(T/stageNumber)) && (bestMode==1)
        modeCorrect=modeCorrect+1;
    elseif (round(T/stageNumber)<k) && (k<=round(2*T/stageNumber)) && (bestMode==2)
        modeCorrect=modeCorrect+1;
    elseif (round(2*T/stageNumber)<k) && (k<=round(3*T/stageNumber)) && (bestMode==4)
        modeCorrect=modeCorrect+1;
    elseif (round(3*T/stageNumber)<k) && (k<=round(4*T/stageNumber)) && (bestMode==1)
        modeCorrect=modeCorrect+1;
    elseif (round(4*T/stageNumber)<k) && (k<=round(5*T/stageNumber)) && (bestMode==4)
        modeCorrect=modeCorrect+1;
    elseif (round(5*T/stageNumber)<k) && (k<=round(6*T/stageNumber)) && (bestMode==2)
        modeCorrect=modeCorrect+1;
    else
        modeError=modeError+1;
    end

    %每T/3帧效果
%     if k<=round(T/6)
%         modeInStage1(k)=bestMode;
%     elseif (round(T/6)<k) && (k<=round(2*T/6))
%         modeInStage2(k-round(T/6))=bestMode;
%     elseif (round(2*T/6)<k) && (k<=round(3*T/6))
%         modeInStage3(k-round(2*T/6))=bestMode;
%     elseif (round(3*T/6)<k) && (k<=round(4*T/6))
%         modeInStage4(k-round(3*T/6))=bestMode;
%     elseif (round(4*T/6)<k) && (k<=round(5*T/6))
%         modeInStage5(k-round(4*T/6))=bestMode;
%     else
%         modeInStage6(k-round(5*T/6))=bestMode;
%     end
    if k<(recSpaceSize+1)
        recSpace(k)=bestMode;
    else
        for p=1:recSpaceSize-1
            recSpace(p)=recSpace(p+1);
        end
        recSpace(recSpaceSize)=bestMode;
    end
    tableNow=tabulate(recSpace);
    [maxCountNow,idxNow] = max(tableNow(:,2));
    recNow=tableNow(idxNow);
    
    %%
    %记录模式变更帧
    if recNow~=tempRec
        recChange(u)=k;
        u=u+1;
        tempRec=recNow;
    end
    %%
    
    %%%统计识别帧数（未丢失帧）
    if recNow~=0
       recCount=recCount+1;
    else
        recMiss=recMiss+1;
    end
    %%%统计识别准确帧
%     if (k<=round(T/6)) && (recNow==1) 
%         recRight=recRight+1;
%     elseif (round(T/stageNumber)<k) && (k<=round(2*T/stageNumber)) && (recNow==2)
%         recRight=recRight+1;
%     elseif (round(2*T/stageNumber)<k) && (k<=round(3*T/stageNumber)) && (recNow==4)
%         recRight=recRight+1;
%     elseif (round(3*T/stageNumber)<k) && (k<=round(4*T/stageNumber)) && (recNow==1)
%         recRight=recRight+1;
%     elseif (round(4*T/stageNumber)<k) && (k<=round(5*T/stageNumber)) && (recNow==4)
%         recRight=recRight+1;
%     elseif (round(4*T/stageNumber)<k) && (k<=round(5*T/stageNumber)) && (recNow==2)
%         recRight=recRight+1;
%     end
    %记录识别结果
    recRecord(k)=recNow;
    %%
    figure(2)
    set(gca,'FontSize',12);
    set(gcf,'position',[50 100 450 400]);
    %clf; 用来清除图形的命令。一般在画图之前用
    clf;
    title('Trace Simulation');
    hold on
    plot(X(1, k), X(2, k), 'r.', 'markersize',20);  %系统状态位置
    grid on;
    axis([0 600 -50 300]);
    plot(P(1, :), P(2, :), 'k.', 'markersize',5);   %各个粒子位置
    plot(PCenter(1, k), PCenter(2, k), 'b.', 'markersize',3); %所有粒子的中心位置
    legend('True State', 'Particle', 'The Center of Particles');
    hold off
    %pause(0.2);
    figure(3);
    set(gcf,'position',[515 100 450 400]);
    %hold on;
    
    subplot(211);plot(k,bestMode,'rx');hold on;grid on;title('Mode Selection');
    axis([0 T+5 0 MNumber+1]);
    subplot(212);plot(k,recNow,'go');hold on;grid on;title('Mode Recgnition');
    axis([0 T+5 0 MNumber+1]);
    
    drawnow;
end
%%
%统计识别准确帧
for n=1:T
    if (n<=round(T/stageNumber)+recMiss) && (recRecord(n)==1) 
        recRight=recRight+1;
    elseif (round(T/stageNumber)+recMiss<n) && (n<=round(2*T/stageNumber)+recMiss) && (recRecord(n)==2)
        recRight=recRight+1;
    elseif (round(2*T/stageNumber)+recMiss<n) && (n<=round(3*T/stageNumber)+recMiss) && (recRecord(n)==4)
        recRight=recRight+1;
    elseif (round(3*T/stageNumber)+recMiss<n) && (n<=round(4*T/stageNumber)+recMiss) && (recRecord(n)==1)
        recRight=recRight+1;
    elseif (round(4*T/stageNumber)+recMiss<n) && (n<=round(5*T/stageNumber)+recMiss) && (recRecord(n)==4)
        recRight=recRight+1;
    elseif (round(5*T/stageNumber)+recMiss<n) && (n<=round(6*T/stageNumber)+recMiss) && (recRecord(n)==2)
        recRight=recRight+1;
    end
end
%%
%识别准确率计算
recACC=recRight/recCount;



%%
%%统计几次识别延时
for i=1:recChange(1)
    timeDelay(1)=timeDelay(1)+timeEachFrame(i);
end
for i=10:recChange(2)
    timeDelay(2)=timeDelay(2)+timeEachFrame(i);
end
for i=20:recChange(3)
    timeDelay(3)=timeDelay(3)+timeEachFrame(i);
end
for i=30:recChange(4)
    timeDelay(4)=timeDelay(4)+timeEachFrame(i);
end
for i=40:recChange(5)
    timeDelay(5)=timeDelay(5)+timeEachFrame(i);
end
for i=50:recChange(6)
    timeDelay(6)=timeDelay(6)+timeEachFrame(i);
end

%%
%统计各阶段最多的模式
% table1 = tabulate(modeInStage1);
% [maxCount1,idx1] = max(table1(:,2));
% firRec=table1(idx1);
% 
% table2 = tabulate(modeInStage2);
% [maxCount2,idx2] = max(table2(:,2));
% secRec=table2(idx2);
% 
% table3 = tabulate(modeInStage3);
% [maxCount3,idx3] = max(table3(:,2));
% thiRec=table3(idx3);
% 
% table4 = tabulate(modeInStage4);
% [maxCount4,idx4] = max(table4(:,2));
% fouRec=table4(idx4);
% 
% table5 = tabulate(modeInStage5);
% [maxCount5,idx5] = max(table5(:,2));
% fifRec=table5(idx5);
% 
% table6 = tabulate(modeInStage6);
% [maxCount6,idx6] = max(table6(:,2));
% sixRec=table6(idx6);
%%
figure(4);
set(gca,'FontSize',12);
plot(X(1,:), X(2,:), 'r', Z(1,:), Z(2,:), 'g', PCenter(1,:), PCenter(2,:), 'b-');
grid on;
%axis([0 250 -50 200]);
legend('True State', 'Measurement', 'Particle Filter');
xlabel('x', 'FontSize', 20); ylabel('y', 'FontSize', 20);

%%
figure(5);
set(gca,'FontSize',12);
plot(err,'.-');
xlabel('t', 'FontSize', 20);
title('The err');

%
fprintf('The accuracy of mode recognition in each frame is:%.2f\n',modeCorrect/T);
% fprintf('The recognition delay of each frame is:%.4f\n',timeOfAll/T);
% fprintf('The mean time delay of recognition is:%.4f\n',sum(timeDelay)/stageNumber);
fprintf('The accuracy of recognition is:%.4f\n',recACC);
fprintf('The average output time after frame is:%.4f\n',outTime/T);
% 
% 
% 
EXP_frameacc(e)=modeCorrect/T;
EXP_acc(e)=recACC;
EXP_time(e)= outTime/T;
% EXP_framedelay(e)=timeOfAll/T;
% EXP_recdelay(e)=sum(timeDelay)/stageNumber;
e
end

meanEXP_frameacc=mean(EXP_frameacc);
meanEXP_acc=mean(EXP_acc);
meanEXP_outputTimeAfterEND=mean(EXP_time);
% meanEXP_framedelay=mean(EXP_framedelay);
% meanEXP_recdelay=mean(EXP_recdelay);
fprintf('The mean frame accuracy in 100 times experiment is:%.4f\n',meanEXP_frameacc);
fprintf('The mean accuracy of recognition in 100 times experiment is:%.4f\n',meanEXP_acc);
fprintf('The mean output time after end in 100 times experiment is:%.4f\n',meanEXP_outputTimeAfterEND);
% fprintf('The mean frame delay in 50 times experiment is:%.4f\n',meanEXP_framedelay);
% fprintf('The mean recognition delay in 50 times experiment is:%.4f\n',meanEXP_recdelay);
figure(6);
plot(1:numberOfExp,EXP_frameacc,'bo');
xlabel('Experiment');ylabel('Mean frame accuracy');
title('Mean Frame accuracy in 100 Experiments');
figure(7);
plot(1:numberOfExp,EXP_acc,'b+');
xlabel('Experiment');ylabel('Mean accuracy of recognition');
title('Mean accuracy of Recognition in 100 Experiments');
figure(8);
plot(1:numberOfExp,EXP_time,'bd');
xlabel('Experiment');ylabel('Mean frame delay');