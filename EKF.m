function [w, s, lambda] = EKF(x,y,w,s,q,r,b)
    
    % Extended Kalman filter.
    %
    % USAGE: [w, s, lambda] = EKF(x,y,w,s,q,r,[b])
    
    if nargin < 7; b = 10; end
    
    for t = 1:length(x)
        decay = 1./(1+exp(-b*(w.^2)));
        w = w.*decay;
        z = 2*w.*decay.*(1-decay) + decay;
        s2 = (z.^2).*s + q;
        lambda = s2.*(x(t)^2) + r;
        alpha = s2./lambda;
        w = w + alpha.*x(t).*(y(t) - w*x(t));
        s = (1-alpha*(x(t)^2)).*s2;
    end