function lik = fit_model(model)
    
    q = 0:30;
    data = load_data(1:3);
    
    for i = 1:length(q)
        disp(num2str(i));
        [~,L] = model_sim(data,model,q(i));
        lik(i) = sum(L);
    end
    [~,ix] = max(lik)