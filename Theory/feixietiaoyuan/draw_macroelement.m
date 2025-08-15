%%画宏单元
clear; clc;
r = 2.0;
N = 6;

theta = 2*pi*(0:N-1)/N;
theta = theta';

a = zeros(N,2);
a(:,1) = r*cos(theta);
a(:,2) = r*sin(theta);


for i = 1:N
    
   fprintf('\\coordinate (c%d) at (%.4f, %.4f);\n',i,a(i,1),a(i,2)) 
    
    
end