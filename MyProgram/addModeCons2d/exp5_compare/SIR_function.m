function [RMSE_SIR, SD_SIR, PCenter_SIR] = SIR_function(N,Q,R,T,st,WorldSize,X,Z)
RMSE_SIR = zeros(1,T);
SD_SIR = zeros(1,T);
P = zeros(2, N);    %建立粒子群
PCenter = zeros(2, T);  %所有粒子的中心位置
w = zeros(N, 1);         %每个粒子的权重
Dis = zeros(N,1);     %每个粒子到真值的距离
%初始化粒子群
for i = 1 : N
    P(:, i) = [WorldSize*rand; WorldSize*rand];
    dist = norm(P(:, i)-Z(:, 1));    
    Dis(i) = dist;
    w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);   %求权重
end
PCenter(:, 1) = sum(P, 2) / N;      %初始化所有粒子的几何中心位置
RMSE_SIR(1) = norm(PCenter(:,1)-Z(:,1));  %初始化RMSE
S = 0;
for j = 1:N
    S = S+Dis(j)^2;  
end
SD_SIR(1) = sqrt(S/N);   %初始化SD
%时间序列开始
for k=2:T
   %预测
    for i = 1 : N
        P(1, i) = P(1, i)  + st + 2.25*wgn(1, 1, 10*log10(Q));
        P(2, i) = P(2, i)  + st + 2.25*wgn(1, 1, 10*log10(Q));
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
        P(:, i) = P(:, index); 
        Dis(i) = norm(P(:,i)-Z(:,k));%统计重采样后的粒子群和真值的误差
    end
    PCenter(:, k) = sum(P, 2) / N;
    RMSE_SIR(k) = norm(PCenter(:,k)-X(:,k));
    S = 0;
    for j = 1:N
        S = S+Dis(j)^2;  
    end
    SD_SIR(k) = sqrt(S/N);
end

PCenter_SIR = PCenter;
