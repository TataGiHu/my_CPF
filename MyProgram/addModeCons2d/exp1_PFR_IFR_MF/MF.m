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
EXP_delay=zeros(1,numberOfExp);
err_Mean=zeros(1,numberOfExp);
err_BAR=zeros(60,numberOfExp);
for e=1:numberOfExp
close all;
N = 500;   %粒子总数
Q = 5;      %过程噪声
R = 5;      %测量噪声
T = 60;     %测量时间
stageNumber=6;%定义阶段数

featureExtractureCoefficient = 10;%设置特征提取系数（针对多少帧进行特征提取）
container_OF_featureFrame = zeros(2, featureExtractureCoefficient);%设置特征帧的容器
referenceMode = [0 0 0 0];%初始化参考模式

% theta = pi/T;       %旋转角度
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
MPro(:,1) = 0.2*ones(5,1);%模式权重
MNumber = length(MPro);%模式数量
recSequence = zeros(1,T/featureExtractureCoefficient);%识别结果存放
modeMark = 1;%模式标记
Delay = zeros(1,featureExtractureCoefficient);
% MCum = cumsum(MPro);%模式权重累加



% modeCorrect=0;
% modeError=0;
% timeOfAll=0;
% outTime=0;%帧内输出时间
% 
% recSpaceSize=3;%定义识别帧基数
% recSpace=zeros(recSpaceSize,1);%定义当前识别基
% timeEachFrame=zeros(T,1);%存放每一帧耗时
% recChange=zeros(stageNumber,1);%记录识别变更帧
% timeDelay=zeros(stageNumber,1);%记录识别变更时延
% tempRec=0;%识别变更栈
% u=1;%用于recChange的循环赋值
% recRight=0;%识别准确数
% recCount=0;%识别帧数（未丢失）
% recMiss=0;%识别丢失帧
% recACC=0;%识别准确率



%初始化粒子群
for i = 1 : N
    P(:, i) = [WorldSize*rand; WorldSize*rand];%在worldSize区域内随机生成N个粒子
    dist = norm(P(:, i)-Z(:, 1));     %与测量位置相差的距离，用于估算该粒子的权重
    %由于上面已经随机生成了N个粒子，现在将其与真实的测量值z进行比较，越接近则权重越大，或者说差值越小权重越大
    %这里的权重计算是关于p(z/x)的分布，即观测方程的分布，假设观测噪声满足高斯分布，那么w（i）=p(z/x)
    w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);   %求权重
end

% wsum = sum(w);
% for i = 1 : N
%    w(i) = w(i) / wsum;
% end
% 
% for i = 1 : N
%        wmax = 2 * max(w) * rand;  %另一种重采样规则,获取权重当中的最大值与均匀分布在0-1之间数的乘积再乘以2作为后面判断选出大权重粒子的依据
%        %randi()函数生成均匀分布的伪随机整数，范围为imin--imax，如果没指定imin，则默认为1
%        %r = randi([imin,imax],...)返回一个在[imin,imax]范围内的伪随机整数
%        index = randi(N, 1);
%        %在这里随机产生N个粒子当中的某个粒子，然后选择该粒子的权值，判断是否是大的粒子权重，是的话就把该粒子选出来
%        while(wmax > w(index))
%            %使vmax递减，不然的话单个粒子的权重是不可能大于vmax的值的
%            wmax = wmax - w(index);
%            index = index + 1;
%            %如果从某个随机的粒子开始找，没找到从该粒子开始到最后的粒子权值总和大于vmax，那么就重新从第一个粒子开始查找
%            if index > N
%                index = 1;
%            end          
%        end
%        P(:, i) = P(:, index);     %从上面的index中得到新粒子  
% end
% PCenter(:, 1) = sum(P, 2) / N;      %所有粒子的几何中心位置
container_OF_featureFrame(:,1) = Z(:,1);
%%
% err(1) = norm(X(:, 1) - PCenter(:, 1));     %粒子几何中心与系统真实状态的误差
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

    tic;
    %模拟一个运动的状态
    %wgn（m,n,p）产生一个m行n列的高斯白噪声的矩阵，p以dBW为单位指定输出噪声的强度。

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

    
    %%
    %观测特征提取
    if(mod(k,featureExtractureCoefficient)~=0)
        
        container_OF_featureFrame(:,mod(k,featureExtractureCoefficient)) = Z(:,k);
        Delay(mod(k,featureExtractureCoefficient))=toc;
    else
        
        container_OF_featureFrame(:,featureExtractureCoefficient) = Z(:,k);
        deltaX = container_OF_featureFrame(1,featureExtractureCoefficient) -  container_OF_featureFrame(1,1);
        deltaY = container_OF_featureFrame(2,featureExtractureCoefficient) -  container_OF_featureFrame(2,1);
        %参考模式设置*（参考模式可能是虚模式）
        if(deltaX>5*10 && abs(deltaY)<5*10)
           referenceMode(1)=1;
           referenceMode(2)=0;
           referenceMode(3)=0;
           referenceMode(4)=0;
           modeMark = 1;
       elseif(deltaX>5*10 && deltaY>5*10)
           referenceMode(1)=1;
           referenceMode(2)=0;
           referenceMode(3)=1;
           referenceMode(4)=0;
           modeMark = 2;
       elseif(deltaX<-5*10 && deltaY>5*10)
           referenceMode(1)=0;
           referenceMode(2)=1;
           referenceMode(3)=1;
           referenceMode(4)=0;
           modeMark = 3;
       elseif(deltaX>5*10 && deltaY<-5*10)
           referenceMode(1)=1;
           referenceMode(2)=0;
           referenceMode(3)=0;
           referenceMode(4)=1;
           modeMark = 4;
       elseif(deltaX<-5*10 && deltaY<-5*10)
           referenceMode(1)=0;
           referenceMode(2)=1;
           referenceMode(3)=0;
           referenceMode(4)=1;
           modeMark = 5;
       elseif(abs(deltaX)<-5*10 && deltaY>5*10)
           referenceMode(1)=0;
           referenceMode(2)=0;
           referenceMode(3)=1;
           referenceMode(4)=0;
           modeMark = 6;
       elseif(deltaX<-5*10 && abs(deltaY)<5*10)
           referenceMode(1)=0;
           referenceMode(2)=1;
           referenceMode(3)=0;
           referenceMode(4)=0;
           modeMark = 7;
       elseif(abs(deltaX)<5*10 && deltaY<-5*10)
           referenceMode(1)=0;
           referenceMode(2)=0;
           referenceMode(3)=0;
           referenceMode(4)=1;
           modeMark = 8;
       else
           referenceMode(1)=0;
           referenceMode(2)=0;
           referenceMode(3)=0;
           referenceMode(4)=0;
           modeMark = 9;
        end
        %识别
        H1=[referenceMode;M(1,:)];
        MPro(1)=1-pdist(H1,'hamming');
        H2=[referenceMode;M(2,:)];
        MPro(2)=1-pdist(H2,'hamming');
        H3=[referenceMode;M(3,:)];
        MPro(3)=1-pdist(H3,'hamming');
        H4=[referenceMode;M(4,:)];
        MPro(4)=1-pdist(H4,'hamming');
        H5=[referenceMode;M(5,:)];
        MPro(5)=1-pdist(H5,'hamming');
        
        [m,pos]=max(MPro);
        recNow=pos;
        recSequence(k/featureExtractureCoefficient)=recNow;
        for i=1:MNumber
            if(i==recNow)
                MPro(i)=0.8;
            else
                MPro(i)=0.05;
            end
        end
        MPro=MPro./sum(MPro);%归一化模式权重
        %预测
        MCum = cumsum(MPro);%模式权重累加
        for n = 1:featureExtractureCoefficient
            
          for i = 1 : N
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
            dist = norm(P(:, i)-container_OF_featureFrame(:,n));    
            w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);  
          end
          %归一化权重
          wsum = sum(w);
          for i = 1 : N
             w(i) = w(i) / wsum;
          end
          %重采样（更新）
          for i = 1 : N
              wmax = 2 * max(w) * rand; 
              index = randi(N, 1);  
              while(wmax > w(index)) 
                  wmax = wmax - w(index);
                  index = index + 1;
                  if index > N
                      index = 1;
                  end          
              end
              P(:, i) = P(:, index);     %从上面的index中得到新粒子  
          end
          %所有粒子的几何中心位置(estimation)
          PCenter(:, k-10+n) = sum(P, 2) / N;    
          
          %计算误差
          err(k-10+n) = norm(X(:, k-10+n) - PCenter(:, k-10+n));     %粒子几何中心与系统真实状态的误差
          
        end
        Delay(featureExtractureCoefficient)=toc;
 
        
        
       

        
        
        

        
    end
    
    
    
         
    
     

     
  
    %%

    
%     
%     table = tabulate(Mq);
%     [maxCount,idx] = max(table(:,2));
%     bestMode=table(idx);
    

    

    %sum(X,2)表示把X按行求和;如果是sum(X),那就是按列求和
    
    
    



    %%
    figure(2)
    set(gca,'FontSize',12);
    set(gcf,'position',[50 100 450 400]);
    %clf; 用来清除图形的命令。一般在画图之前用
    clf;
    title('Trace Simulation');
    hold on
    plot(X(1, k), X(2, k), 'r.', 'markersize',20);  %系统状态位置
    %grid on;
    %axis([0 500 -60 120]);
    plot(P(1, :), P(2, :), 'k.', 'markersize',5);   %各个粒子位置
    %plot(PCenter(1, k), PCenter(2, k), 'b.', 'markersize',20); %所有粒子的中心位置
    legend('True State', 'Particle');

    hold off
    pause(0.2);
%     figure(3);
%     set(gcf,'position',[515 100 450 400]);
%     subplot(211);plot(k,bestMode,'rx');hold on;grid on;title('Mode Selection');
%     axis([0 T+5 0 MNumber+1]);
%     subplot(212);plot(k,recNow,'go');hold on;grid on;title('Mode Recgnition');
%     axis([0 T+5 0 MNumber+1]);
    
    drawnow;
end



%%
figure(4);
set(gca,'FontSize',12);
plot(X(1,:), X(2,:), 'r', Z(1,:), Z(2,:), 'g', PCenter(1,:), PCenter(2,:), 'b');
%grid on;
%axis([0 250 -50 200]);
legend('True State', 'Measurement', 'Particle Filter');
xlabel('x', 'FontSize', 20); ylabel('y', 'FontSize', 20);

%%
figure(5);
set(gca,'FontSize',12);
plot(err,'r+-');
xlabel('t', 'FontSize', 20);
title('The err');

recSequence

err_Mean(e)=mean(err);
EXP_delay(e) = sum(Delay);
mean(err)
mean(Delay)
e
end

meanEXP_err=mean(err_Mean)
meanEXP_delay = mean(EXP_delay)
% meanEXP_frameacc=mean(EXP_frameacc);
% meanEXP_acc=mean(EXP_acc);
% meanEXP_outputTimeAfterSTART=mean(EXP_time);
% 
% fprintf('The mean frame accuracy in 100 times experiment is:%.4f\n',meanEXP_frameacc);
% fprintf('The mean accuracy of recognition in 100 times experiment is:%.4f\n',meanEXP_acc);
% fprintf('The mean output time after start in 100 times experiment is:%.4f\n',meanEXP_outputTimeAfterSTART);
% 
% figure(6);
% plot(1:numberOfExp,EXP_frameacc,'bo');
% xlabel('Experiment');ylabel('Mean frame accuracy');
% title('Mean Frame accuracy in 100 Experiments');
% figure(7);
% plot(1:numberOfExp,EXP_acc,'b+');
% xlabel('Experiment');ylabel('Mean accuracy of recognition');
% title('Mean accuracy of Recognition in 100 Experiments');
% figure(8);
% plot(1:numberOfExp,EXP_time,'bd');
% xlabel('Experiment');ylabel('Mean frame delay');
