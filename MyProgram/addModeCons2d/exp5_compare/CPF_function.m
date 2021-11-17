function [RMSE_CPF, SD_CPF, PCenter_CPF] = CPF_function(N,Q,R,T,st,WorldSize,X,Z)
RMSE_CPF = zeros(1,T);
SD_CPF = zeros(1,T);
P = zeros(2, N);    %建立粒子群
PCenter = zeros(2, T);  %所有粒子的中心位置
w = zeros(N, 1);         %每个粒子的权重
Dis = zeros(N,1);     %每个粒子到真值的距离

M=[1 0 0 0;%模式矩阵
   1 0 1 0;
   0 1 1 0;
   1 0 0 1;
   0 1 0 1];
Mq = zeros(N,1);%传递模式索引
Mqu = zeros(N,T);%存储模式索引
MPro(:,1) = 0.2*ones(5,1);%模式权重
MNumber = length(MPro);%模式数量
MCum = cumsum(MPro);%模式权重累加

%初始化粒子群
for i = 1 : N
    P(:, i) = [WorldSize*rand; WorldSize*rand];
    dist = norm(P(:, i)-Z(:, 1));    
    Dis(i) = dist;
    w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);   %求权重
end
PCenter(:, 1) = sum(P, 2) / N;      %初始化所有粒子的几何中心位置
RMSE_CPF(1) = norm(PCenter(:,1)-Z(:,1));  %初始化RMSE
S = 0;
for j = 1:N
    S = S+Dis(j)^2;  
end
SD_CPF(1) = sqrt(S/N);   %初始化SD
%时间序列开始
for k=2:T
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
        P(1, i) = P(1, i) + st * (thisMode(1)*1.5+thisMode(2)*(-1.5)) + wgn(1, 1, 10*log10(Q));
        P(2, i) = P(2, i) + st * (thisMode(3)*1.5+thisMode(4)*(-1.5)) + wgn(1, 1, 10*log10(Q));
        dist = norm(P(:, i)-Z(:, k));     %与测量位置相差的距离，用于估算该粒子的权重
        %由于上面已经随机生成了N个粒子，现在将其与真实的测量值z进行比较，越接近则权重越大，或者说差值越小权重越大
        %这里的权重计算是关于p(z/x)的分布，即观测方程的分布，假设观测噪声满足高斯分布，那么w（i）=p(z/x)
        w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);   %求权重
    end
    %Normalization
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
        P(:, i) = P(:, index); 
        Dis(i) = norm(P(:,i)-Z(:,k));%统计重采样后的粒子群和真值的误差
        Mq(i) = Mq(index); %选出粒子对应的模式索引
    end
    
%     table = tabulate(Mq);
%     [maxCount,idx] = max(table(:,2));
%     bestMode=table(idx);
    
    Mqu(:,k)=Mq;
    MPro(1,1)=length(find(Mq==1));
    MPro(2,1)=length(find(Mq==2));
    MPro(3,1)=length(find(Mq==3));
    MPro(4,1)=length(find(Mq==4));
    MPro(5,1)=length(find(Mq==5));
    MPro = MPro/sum(MPro);
    for i=1:MNumber
       if MPro(i) < 0.05
          MPro= 0.2*ones(5,1);
       end
    end
    MCum = cumsum(MPro);
    
    
    PCenter(:, k) = sum(P, 2) / N;
    RMSE_CPF(k) = norm(PCenter(:,k)-X(:,k));
    S = 0;
    for j = 1:N
        S = S+Dis(j)^2;  
    end
    SD_CPF(k) = sqrt(S/N);
end

PCenter_CPF = PCenter;
