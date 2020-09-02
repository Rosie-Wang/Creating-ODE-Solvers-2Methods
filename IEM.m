function [T,Y]= IEM(t0,tN,y0,h,f) 
n=(tN-t0)/h; %h is the uniform step size; n is the number of steps 
T=t0;
Y=y0;
for i=1:n
   k1= f(t0,y0);
   k2= f(t0+h,y0+h*k1);
   k = (k1+k2)/2;
   t0=t0+h;
   y0=y0+h*k;
   T=[T;t0]; %horizontal vector 
   Y=[Y;y0];
end

