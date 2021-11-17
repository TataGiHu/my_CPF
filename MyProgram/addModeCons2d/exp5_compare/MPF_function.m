function [RMSE_MPF, SD_MPF, PCenter_MPF] = MPF_function(N,Q,R,T,st,WorldSize,X,Z)
RMSE_MPF = zeros(1,T);
SD_MPF = zeros(1,T);
P = zeros(2, N);    %建立粒子群
PCenter = zeros(2, T);  %所有粒子的中心位置
w = zeros(N, 1);         %每个粒子的权重
Dis = zeros(N,1);     %每个粒子到真值的距离
key = 5;   %聚类类数
%初始化粒子群
for i = 1 : N
    P(:, i) = [WorldSize*rand; WorldSize*rand];
    dist = norm(P(:, i)-Z(:, 1));    
    Dis(i) = dist;
    w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);   %求权重
end
PCenter(:, 1) = sum(P, 2) / N;      %初始化所有粒子的几何中心位置
RMSE_MPF(1) = norm(PCenter(:,1)-Z(:,1));  %初始化RMSE
S = 0;
for j = 1:N
    S = S+Dis(j)^2;  
end
SD_MPF(1) = sqrt(S/N);   %初始化SD
%时间序列开始
for k=2:T
    %预测
    for i = 1 : N
        P(1, i) = P(1, i)  + st + 2.5*wgn(1, 1, 10*log10(Q));
        P(2, i) = P(2, i)  + st + 2.5*wgn(1, 1, 10*log10(Q));
        dist = norm(P(:, i)-Z(:, k));    
        w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);
    end
    
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
        P(:, i) = P(:, index)+ 0.01*wgn(1, 1, 10*log10(Q)); 
        Dis(i) = norm(P(:,i)-Z(:,k));%统计重采样后的粒子群和真值的误差
    end
     %%%%%%%%%%%%%%%%%%%%%%   聚类为k类 并分配类权重   %%%%%%%%%%%%
    [IDX,C] = kmeans(P',key);
    C=C';
    modeWeight = zeros(1,key);
    for j=1:key
        modeWeight(j) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(norm(C(:,j)-Z(:,k)))^2 / 2 / R);
    end
    for j=1:key
       modeWeight(j) = modeWeight(j)/sum(modeWeight); 
    end
    cumModeW = cumsum(modeWeight);
    %%%%%%%%%%%%%%%%%%%%    基于各类权重作模拟重要性采样（仅用于输出，放大大权重类中的粒子影响力）    %%%%%%
    P_temp = zeros(2,N);
    for i=1:N
        id =0;
        while(id==0)
            temp_id = round(rand*N);
            if(temp_id>0 && rand<cumModeW(IDX(temp_id)))
               id = temp_id; 
            end
        end
        P_temp(:,i) = P(:,id);
    end
    PCenter(:, k) = sum(P_temp, 2) / N;
    RMSE_MPF(k) = norm(PCenter(:,k)-X(:,k));
    S = 0;
    for j = 1:N
        S = S+Dis(j)^2;  
    end
    SD_MPF(k) = sqrt(S/N);
end

PCenter_MPF = PCenter;
