clc;
close all;


%%
x_err=[50 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500];
y_err=[3.2666 2.8616 2.6886 2.6366 2.5873 2.5552 2.5478 2.5197 2.5419 2.5294 2.5109 2.5152 2.5176 2.5077 2.5080 2.5124];

x_FACC=[50 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500];
y_FACC=[74.20 75.50 77.60 77.50 78.15 78.30 79.25 79.35 79.40 79.47 79.40 78.45 79.42 79.53 79.17 78.70];

x_GACC=[50 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500];
y_GACC=[83.86 84.27 86.46 86.59 87.20 87.18 87.44 87.89 87.60 87.24 87.59 87.53 88.06 87.29 87.50 87.56];

x_IMP=[50 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500];
y_IMP=[9.66 8.77 8.86 9.09 9.05 8.88 8.19 8.54 8.20 7.77 8.19 8.58 8.64 7.76 8.34 8.86];

x_RD=[50 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500];
y_RD=[0.0011 0.0021 0.0049 0.0057 0.0097 0.0098 0.0132 0.0173 0.0174 0.0219 0.0258 0.0255 0.0318 0.0307 0.0328 0.0369];


%%
f_err=polyfit(x_err,y_err,7);
xi_err=0:0.05:1500;
yi_err=polyval(f_err,xi_err);

f_FACC=polyfit(x_FACC,y_FACC,3);
xi_FACC=0:0.05:1500;
yi_FACC=polyval(f_FACC,xi_FACC);

f_GACC=polyfit(x_GACC,y_GACC,3);
xi_GACC=0:0.05:1500;
yi_GACC=polyval(f_GACC,xi_GACC);

f_IMP=polyfit(x_IMP,y_IMP,12);
xi_IMP=0:0.05:1500;
yi_IMP=polyval(f_IMP,xi_IMP);

f_RD=polyfit(x_RD,y_RD,3);
xi_RD=0:0.05:1500;
yi_RD=polyval(f_RD,xi_RD);

%%
fig1 = figure(1);
% left_color=[0 0 0];
% right_color=[0 0 0];
% set(fig1,'defaultAxesColorOrder',[left_color; right_color]);


yyaxis left
plot(xi_FACC,yi_FACC,'-.',x_FACC,y_FACC,'r+');hold on;set(gca,'XTick',0:200:1500);
plot(xi_GACC,yi_GACC,'--',x_GACC,y_GACC,'ro');hold on;set(gca,'XTick',0:200:1500);
xlabel('Swarm size');
ylabel('Accuracy(%)');
axis([50 1500 72 97]);
% set(gca,'YTick',[98.8 99 99.2 99.4]);

yyaxis right
plot(xi_IMP,yi_IMP,'-',x_IMP,y_IMP,'kx');
ylabel('Improvement(%)');

legend('FACC','','GACC','','IMP','','Location','Northeast','Orientation','Vertical');

%%
fig2 = figure(2);
% left_color=[0 0 0];
% right_color=[0 0 0];
% set(fig1,'defaultAxesColorOrder',[left_color; right_color]);


yyaxis left
plot(xi_err,yi_err,'-.',x_err,y_err,'r+');hold on;

xlabel('Swarm size');
ylabel('Error');



yyaxis right
plot(xi_RD,yi_RD,'-',x_RD,y_RD,'kx');
ylabel('Recognition delay(s)');

legend('Err','','RD','','Location','Northeast','Orientation','Vertical');
