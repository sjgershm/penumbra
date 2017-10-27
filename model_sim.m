function [r,lik,data] = model_sim(conds,model)
    
    % Simulate hierarchical Kalman filter model
    
    data = []; C = [];
    for cond = 1:length(conds)
        d = load_data(conds(cond));
        C = [C zeros(1,length(d))+cond];
        data = [data d];
    end
    
    R = [25 25 100];
    Q = [25 1 25];
    
    for s = 1:length(data)
        p = data(s).p;
        ix = ~isnan(data(s).p);
        i = C(s);
        results = run_model(data(s),model,R(i),Q(i));
        r(s,1) = corr(results.p(ix),p(ix));
        sd = sqrt(mean((results.p(ix)-p(ix)).^2));
        lik(s,1) = sum(log(normpdf(results.p(ix),p(ix),sd)));
        data(s).p = results.p;
    end

end

function results = run_model(data,model,r,Q)
    N = size(data.x,1);
    x = zeros(N,2);
    x(data.x==1,1) = 1;
    x(data.x==2,2) = 1;
    data.y = data.y-50;
    y = zeros(N,2);
    y(data.x==1,1) = data.y(data.x==1);
    y(data.x==2,2) = data.y(data.x==2);
    
    if model == 1
        Q = linspace(1,30,50)';
    end
    if nargin < 3; r = 25; end
    DA = ones(N,1);
    PRP = ones(N,1);
    b = zeros(N,1)+10;
    theta = [8 5];
    
    [results.w, results.q, results.Phi, results.ww] = hierarchical_particle_filter(x,y,Q,r,DA,PRP,b,theta);
    
    w = [0 0; results.w(1:end-1,:)];
    p = x.*w;
    results.p = zeros(N,1);
    results.p(data.x==1) = p(data.x==1,1);
    results.p(data.x==2) = p(data.x==2,2);
    results.p = results.p+50;
end