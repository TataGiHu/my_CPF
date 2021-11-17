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
N = 300;   %��������
Q = 5;      %��������
R = 5;      %��������
T = 60;     %����ʱ��
theta = pi/T;       %��ת�Ƕ�
distance = 4*80/T;    %ÿ���ߵľ���
WorldSize = 100;    %�����С
X = zeros(2, T);    %�洢ϵͳ״̬
Z = zeros(2, T);    %�洢ϵͳ�Ĺ۲�״̬
P = zeros(2, N);    %��������Ⱥ
PCenter = zeros(2, T);  %�������ӵ�����λ��
w = zeros(N, 1);         %ÿ�����ӵ�Ȩ��
err = zeros(1,T);     %���
X(:, 1) = [5; 20];     %��ʼϵͳ״̬
%wgn��m,n,p������һ��m��n�еĸ�˹�������ľ���p��dBWΪ��λָ�����������ǿ�ȡ�
Z(:, 1) = [5; 20] + wgn(2, 1, 10*log10(R));    %��ʼϵͳ�Ĺ۲�״̬
M=[1 0 0 0;%ģʽ����
   1 0 1 0;
   0 1 1 0;
   1 0 0 1;
   0 1 0 1];
Mq = zeros(N,1);%����ģʽ����
Mqu = zeros(N,T);%�洢ģʽ����
MPro(:,1) = 0.2*ones(5,1);%��ʼ��ģʽȨ��
MNumber = length(MPro);%ģʽ����
MCum = cumsum(MPro);%ģʽȨ���ۼ�
modeCorrect=0;
modeError=0;
timeOfAll=0;
stageNumber=6;%����׶���
recSpaceSize=3;%����ʶ��֡����
recSpace=zeros(recSpaceSize,1);%���嵱ǰʶ���
timeEachFrame=zeros(T,1);%���ÿһ֡��ʱ
recChange=zeros(stageNumber,1);%��¼ʶ����֡
timeDelay=zeros(stageNumber,1);%��¼ʶ����ʱ��
tempRec=0;%ʶ����ջ
u=1;%����recChange��ѭ����ֵ

recRight=0;%ʶ��׼ȷ��
recCount=0;%ʶ��֡����δ��ʧ��
recMiss=0;%ʶ��ʧ֡
recACC=0;%ʶ��׼ȷ��

outTime=0;
XX=zeros(2,8);%�洢��ʵ״̬�Ĳ�ͬת��(����ת��)
ZZ=zeros(2,8);%�洢�۲�״̬�Ĳ�ͬת�ƣ�����ת�ƣ�
EE=zeros(1,8);%�洢�������

recRecord=zeros(T+10,1);
% modeInStage1=zeros(T/6,1);
% modeInStage2=zeros(T/6,1);
% modeInStage3=zeros(T/6,1);
% modeInStage4=zeros(T/6,1);
% modeInStage5=zeros(T/6,1);
% modeInStage6=zeros(T/6,1);
%��ʼ������Ⱥ
for i = 1 : N
    P(:, i) = [WorldSize*rand; WorldSize*rand];%��worldSize�������������N������
    dist = norm(P(:, i)-Z(:, 1));     %�����λ�����ľ��룬���ڹ�������ӵ�Ȩ��
    %���������Ѿ����������N�����ӣ����ڽ�������ʵ�Ĳ���ֵz���бȽϣ�Խ�ӽ���Ȩ��Խ�󣬻���˵��ֵԽСȨ��Խ��
    %�����Ȩ�ؼ����ǹ���p(z/x)�ķֲ������۲ⷽ�̵ķֲ�������۲����������˹�ֲ�����ôw��i��=p(z/x)
    w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);   %��Ȩ��
end
PCenter(:, 1) = sum(P, 2) / N;      %�������ӵļ�������λ��

%%
err(1) = norm(X(:, 1) - PCenter(:, 1));     %���Ӽ���������ϵͳ��ʵ״̬�����
% figure(1);
% %����set(gca,'propertyname','propertyvalue'......)������Ե���ͼ�ε��������ԡ�
% %propertyname��property value�ֱ�Ϊ�����������ƣ���������ĸ��Сд������Ӧ������ֵ��
% set(gca,'FontSize',10);
% hold on%ʹ��ǰ��ͼ�α��ֶ�����ˢ�£�׼�����ܴ˺󽫻��Ƶ�������
% plot(X(1, 1), X(2, 1), 'r.', 'markersize',30)   %ϵͳ״̬λ��
% %axis[xmin,xmax,ymin,ymax]�趨���귶Χ����������xmin<xmax������ͬ��
% axis([0 100 0 100]);
% plot(P(1, :), P(2, :), 'kx', 'markersize',5);   %�����������λ��
% plot(PCenter(1, 1), PCenter(2, 1), 'b.', 'markersize',25); %�������ӵļ�������λ��
% %ͼ�α�ע����legend(string1,string2,string3...),�ڵ�ǰͼ�����ͼ��
%  
% legend('True State', 'Particles', 'The Center of Particles');
% %���ͼ�α�������title('string'),�ڵ�ǰ����ϵ�Ķ�����һ���ı���string,��Ϊ��ͼ�εı���
% title('Initial State');
% %ʹ��ǰ�ἰͼ�β��پ߱�����ˢ�µ�����
% hold off

%%
%��ʼ�˶�
for k = 2 : T

    %ģ��һ�������˶���״̬
    %wgn��m,n,p������һ��m��n�еĸ�˹�������ľ���p��dBWΪ��λָ�����������ǿ�ȡ�
%     X(:, k) = X(:, k-1) + distance * [(-cos(k * theta)); sin(k * theta)] + wgn(2, 1, 10*log10(Q));     %״̬����
%     Z(:, k) = X(:, k) + wgn(2, 1, 10*log10(R));     %�۲ⷽ�� 
    if k<=round(T/stageNumber)
        X(1, k) = X(1, k-1) + distance * 1.5 + wgn(1, 1, 10*log10(Q));     %״̬����
        X(2, k) = X(2, k-1) + distance * 0 + wgn(1, 1, 10*log10(Q));
        Z(:, k) = X(:, k) + wgn(2, 1, 10*log10(R));     %�۲ⷽ��    
    elseif (round(T/stageNumber)<k) && (k<=round(2*T/stageNumber))
        X(1, k) = X(1, k-1) + distance *  1.5 + wgn(1, 1, 10*log10(Q));     %״̬����
        X(2, k) = X(2, k-1) + distance * 1.5 + wgn(1, 1, 10*log10(Q));
        Z(:, k) = X(:, k) + wgn(2, 1, 10*log10(R));     %�۲ⷽ�� 
    elseif (round(2*T/stageNumber)<k) && (k<=round(3*T/stageNumber))
        X(1, k) = X(1, k-1) + distance *  1.5 + wgn(1, 1, 10*log10(Q));     %״̬����
        X(2, k) = X(2, k-1) + distance * -1.5 + wgn(1, 1, 10*log10(Q));
        Z(:, k) = X(:, k) + wgn(2, 1, 10*log10(R));     %�۲ⷽ�� 
    elseif (round(3*T/stageNumber)<k) && (k<=round(4*T/stageNumber))
        X(1, k) = X(1, k-1) + distance * 1.5 + wgn(1, 1, 10*log10(Q));     %״̬����
        X(2, k) = X(2, k-1) + distance * 0 + wgn(1, 1, 10*log10(Q));
        Z(:, k) = X(:, k) + wgn(2, 1, 10*log10(R));     %�۲ⷽ��
    elseif (round(4*T/stageNumber)<k) && (k<=round(5*T/stageNumber))
        X(1, k) = X(1, k-1) + distance *  1.5 + wgn(1, 1, 10*log10(Q));     %״̬����
        X(2, k) = X(2, k-1) + distance * -1.5 + wgn(1, 1, 10*log10(Q));
        Z(:, k) = X(:, k) + wgn(2, 1, 10*log10(R));     %�۲ⷽ�� 
    else
        X(1, k) = X(1, k-1) + distance *  1.5 + wgn(1, 1, 10*log10(Q));     %״̬����
        X(2, k) = X(2, k-1) + distance * 1.5 + wgn(1, 1, 10*log10(Q));
        Z(:, k) = X(:, k) + wgn(2, 1, 10*log10(R));     %�۲ⷽ��
    end

    %�����˲�
    tic;
    %Ԥ��
    for i = 1 : N
        %��֮ǰ���ɵ����Ӵ��뵽״̬���̣�������һ��״̬��Ԥ��
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
        dist = norm(P(:, i)-Z(:, k));     %�����λ�����ľ��룬���ڹ�������ӵ�Ȩ��
        %���������Ѿ����������N�����ӣ����ڽ�������ʵ�Ĳ���ֵz���бȽϣ�Խ�ӽ���Ȩ��Խ�󣬻���˵��ֵԽСȨ��Խ��
        %�����Ȩ�ؼ����ǹ���p(z/x)�ķֲ������۲ⷽ�̵ķֲ�������۲����������˹�ֲ�����ôw��i��=p(z/x)
        w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);   %��Ȩ��
    end
%��һ��Ȩ��
    wsum = sum(w);
    for i = 1 : N
        w(i) = w(i) / wsum;
    end

    %�ز��������£�
    for i = 1 : N
        wmax = 2 * max(w) * rand;  %��һ���ز�������,��ȡȨ�ص��е����ֵ����ȷֲ���0-1֮�����ĳ˻��ٳ���2��Ϊ�����ж�ѡ����Ȩ�����ӵ�����
        %randi()�������ɾ��ȷֲ���α�����������ΧΪimin--imax�����ûָ��imin����Ĭ��Ϊ1
        %r = randi([imin,imax],...)����һ����[imin,imax]��Χ�ڵ�α�������
        index = randi(N, 1);
        %�������������N�����ӵ��е�ĳ�����ӣ�Ȼ��ѡ������ӵ�Ȩֵ���ж��Ƿ��Ǵ������Ȩ�أ��ǵĻ��ͰѸ�����ѡ����
        while(wmax > w(index))
            %ʹvmax�ݼ�����Ȼ�Ļ��������ӵ�Ȩ���ǲ����ܴ���vmax��ֵ��
            wmax = wmax - w(index);
            index = index + 1;
            %�����ĳ����������ӿ�ʼ�ң�û�ҵ��Ӹ����ӿ�ʼ����������Ȩֵ�ܺʹ���vmax����ô�����´ӵ�һ�����ӿ�ʼ����
            if index > N
                index = 1;
            end          
        end
        P(:, i) = P(:, index);     %�������index�еõ�������
%         Mq(i) = Mq(index); %ѡ�����Ӷ�Ӧ��ģʽ����
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
    

    %sum(X,2)��ʾ��X�������;�����sum(X),�Ǿ��ǰ������
    PCenter(:, k) = sum(P, 2) / N;      %�������ӵļ�������λ��

    %�������
    err(k) = norm(X(:, k) - PCenter(:, k));     %���Ӽ���������ϵͳ��ʵ״̬�����
    
    %���㲻ͬת�Ƶ����
    XX(1,1)= X(1, k-1) + distance * 1.5 + wgn(1, 1, 10*log10(Q));     %״̬����
    XX(2,1) = X(2, k-1) + distance * 0 + wgn(1, 1, 10*log10(Q));   
    ZZ(:,1) = XX(:, 1) + wgn(2, 1, 10*log10(R));     %�۲ⷽ��
    EE(1)=norm(Z(:,k)-ZZ(:,1));
    
    XX(1,2)= X(1, k-1) + distance * 1.5 + wgn(1, 1, 10*log10(Q));     %״̬����
    XX(2,2) = X(2, k-1) + distance * 1.5 + wgn(1, 1, 10*log10(Q));
    ZZ(:,2) = XX(:, 2) + wgn(2, 1, 10*log10(R));     %�۲ⷽ��
    EE(2)=norm(Z(:,k)-ZZ(:,2));
    
    XX(1,3)= X(1, k-1) + distance * -1.5 + wgn(1, 1, 10*log10(Q));     %״̬����
    XX(2,3) = X(2, k-1) + distance * 1.5 + wgn(1, 1, 10*log10(Q));
    ZZ(:,3) = XX(:, 3) + wgn(2, 1, 10*log10(R));     %�۲ⷽ��
    EE(3)=norm(Z(:,k)-ZZ(:,3));
    
    XX(1,4)= X(1, k-1) + distance * 1.5 + wgn(1, 1, 10*log10(Q));     %״̬����
    XX(2,4) = X(2, k-1) + distance * -1.5 + wgn(1, 1, 10*log10(Q));
    ZZ(:,4) = XX(:, 4) + wgn(2, 1, 10*log10(R));     %�۲ⷽ��
    EE(4)=norm(Z(:,k)-ZZ(:,4));
    
    XX(1,5)= X(1, k-1) + distance * -1.5 + wgn(1, 1, 10*log10(Q));     %״̬����
    XX(2,5) = X(2, k-1) + distance * -1.5 + wgn(1, 1, 10*log10(Q));
    ZZ(:,5) = XX(:, 5) + wgn(2, 1, 10*log10(R));     %�۲ⷽ�� 
    EE(5)=norm(Z(:,k)-ZZ(:,5));
    
    

%     
%     
%      
%     
    %��ȡ���ģʽ
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
    %ͳ��ģʽ׼ȷ֡
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

    %ÿT/3֡Ч��
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
    %��¼ģʽ���֡
    if recNow~=tempRec
        recChange(u)=k;
        u=u+1;
        tempRec=recNow;
    end
    %%
    
    %%%ͳ��ʶ��֡����δ��ʧ֡��
    if recNow~=0
       recCount=recCount+1;
    else
        recMiss=recMiss+1;
    end
    %%%ͳ��ʶ��׼ȷ֡
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
    %��¼ʶ����
    recRecord(k)=recNow;
    %%
    figure(2)
    set(gca,'FontSize',12);
    set(gcf,'position',[50 100 450 400]);
    %clf; �������ͼ�ε����һ���ڻ�ͼ֮ǰ��
    clf;
    title('Trace Simulation');
    hold on
    plot(X(1, k), X(2, k), 'r.', 'markersize',20);  %ϵͳ״̬λ��
    grid on;
    axis([0 600 -50 300]);
    plot(P(1, :), P(2, :), 'k.', 'markersize',5);   %��������λ��
    plot(PCenter(1, k), PCenter(2, k), 'b.', 'markersize',3); %�������ӵ�����λ��
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
%ͳ��ʶ��׼ȷ֡
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
%ʶ��׼ȷ�ʼ���
recACC=recRight/recCount;



%%
%%ͳ�Ƽ���ʶ����ʱ
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
%ͳ�Ƹ��׶�����ģʽ
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