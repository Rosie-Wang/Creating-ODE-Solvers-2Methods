function [T,Y] = Aeuler(t0,tN,y0,h,f)
    i = 1;
    T = zeros(1,i); % preallocating t vector
    Y = zeros(1,i); % preallocating y vector
    Y(1) = y0; % initial conditions
    T(1) = t0;
    H = h; % allowing h to reset
    tol = 1e-8; % tolerance
    while (T(i) < tN) % until vectors are full
        h = H; % reset h
        y = Y(i) + h*f(T(i), Y(i)); % creating y and z using euler
        z1 = Y(i) + (0.5*h)*f(T(i), Y(i));
        z = z1 + (h*0.5)*(f(T(i) + (h/2),z1));
        D = abs(z-y); % comparing 
        while D>= tol % if diff too large, try new h until it isn't
            h = 0.9*h*min(max(tol/abs(D),0.3),2); % new h
            y = Y(i) + h*f(T(i), Y(i)); % recalculating y and z
            z1 = Y(i) + (0.5*h)*f(T(i), Y(i));
            z = z1 + (h*0.5)*(f(T(i) + (h/2),z1));
            D = abs(z-y); % compare again
        end
        Y(i+1) = z; % add the correct approx. to y vec
        T(i+1) = T(i) + h; % move onto next t value
        i = i+1; % up the index
    end
end