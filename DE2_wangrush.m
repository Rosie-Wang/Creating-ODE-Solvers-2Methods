function [t,y] = DE2_wangrush(t0,tN,y0,y1,h,p,q,g)

   N = (tN-t0)/h; % number of steps, floor will round to next smallest integer 

   % preallocating arrays
   t = t0:h:tN;
   
   y = zeros(1, length(t));
   y_p = zeros(1, length(t));
   y_pp = zeros(1, length(t));

   % ICs
   y(1) = y0;
   y(2)  = y0 + (h*y1);
   y_p(1) = y1;
   y_p(2) = (y(2)-y(1))/h;

   for (i = 2:N) %loop through    
      y_pp = g(t(i)) - q(t(i)) * y(i) - p(t(i)) * y_p(i);
      y(i+1) =  (y_pp*h^2) + 2*y(i) - y(i-1);
      y_p(i+1) = (y(i+1) - y(i))/h; 
   end

end