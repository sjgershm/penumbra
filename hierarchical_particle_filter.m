function [W, q, Phi, ww, S] = hierarchical_particle_filter(x,y,Q,r,DA,PRP,b,theta)
    
    % Hierarchical particle filter.
    %
    % USAGE: [W q Phi ww S] = hierarchical_particle_filter(x,y,Q,r,DA,PRP,b,theta)
    
    % initialize
    [N, D] = size(x);
    P = length(Q);
    phi = theta(1)./((1:P).^theta(2)); phi=phi';
    w = zeros(P,D);               % synaptic weights
    s = ones(P,D);                      % synaptic variances
    Phi = zeros(N,P);
    W = zeros(N,D);
    q = zeros(N,1);
    if nargin < 7; b = zeros(N,1)+10; end
    
    for t = 1:N
        w0 = w;
        for d = 1:D
            [w(:,d), s(:,d), lambda] = EKF(x(t,d),y(t,d),w(:,d),s(:,d),Q,r,b(t));
            v = y(t,d)-w0(:,d).*x(t,d);
            v = v.*DA(t);
            phi = phi - (v.^2)./(2*lambda) - 0.5*log(lambda);
        end
        Phi(t,:) = exp(phi-logsumexp(phi));
        Phi(t,:) = PRP(t).*Phi(t,:) + (1-PRP(t)).*[1 zeros(1,P-1)];
        W(t,:) = Phi(t,:)*w;
        q(t,1) = Phi(t,:)*Q;
        S(t,:) = Phi(t,:)*s;
        ww(t,:,:) = w;
    end