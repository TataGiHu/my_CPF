clc;
close all;

fittingNumber=4;
x=[3 4 5 6 7 8 9 10 11 ...
   3 4 5 6 7 8 9 10 11 ...
   3 4 5 6 7 8 9 10 11 ...
   3 4 5 6 7 8 9 10 11 ...
   3 4 5 6 7 8 9 10 11 ...
   3 4 5 6 7 8 9 10 11];
y=[9.30 7.19 8.89 7.33 6.85 6.06 4.17 2.42 1.49 ...
   9.22 6.76 8.64 6.29 4.58 5.68 3.83 1.42 2.51 ...
   9.30 5.94 8.04 6.72 6.59 5.89 4.20 2.18 2.19 ...
   8.25 6.70 7.71 8.45 7.87 4.18 5.56 3.30 1.82 ...
   8.71 5.96 8.62 7.51 7.47 4.30 5.16 3.24 -0.48 ...
   8.46 6.52 9.42 7.81 6.84 4.99 5.09 3.22 1.26];


f=polyfit(x,y,fittingNumber);
xi=3:0.02:11;
yi=polyval(f,xi);

figure();
plot(xi,yi,x,y,'r+');hold on;xlabel('Recognition base size');ylabel('Improvement degree');

x1=[3 5 7 9 11 ...
   3 5 7 9 11 ...
   3 5 7 9 11 ...
   3 5 7 9 11 ...
   3 5 7 9 11 ...
   3 5 7 9 11];
y1=[9.30 8.89 6.85 4.17 1.49 ...
   9.22 8.64 4.58 3.83 2.51 ...
   9.30 8.04 6.59 4.20 2.19 ...
   8.25 7.71 7.87 5.56 1.82 ...
   8.71 8.62 7.47 5.16 -0.48 ...
   8.46 9.42 6.84 5.09 1.26];
x2=[4 6 8 10 ...
    4 6 8 10 ...
    4 6 8 10 ...
    4 6 8 10 ...
    4 6 8 10 ...
    4 6 8 10];
y2=[7.19 7.33 6.06 2.42 ...
   6.76 6.29 5.68 1.42 ...
   5.94 6.72 5.89 2.18 ...
   6.70 8.45 4.18 3.30 ...
   5.96 7.51 4.30 3.24 ...
   6.52 7.81 4.99 3.22];

fittingNumber2=2;
f1=polyfit(x1,y1,fittingNumber2);
xi1=3:0.02:11;
yi1=polyval(f1,xi1);

f2=polyfit(x2,y2,fittingNumber2);
xi2=3:0.02:11;
yi2=polyval(f2,xi2);

figure();

plot(xi1,yi1,x1,y1,'^');hold on;xlabel('Recognition base size');ylabel('Improvement degree','linewidth',2);
plot(xi2,yi2,x2,y2,'v');hold on;xlabel('Recognition base size');ylabel('Improvement degree','linewidth',2);
legend('Odd fitting','','Even fitting','');
