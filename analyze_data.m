function varargout = analyze_data(mode,cond)
    
    if nargin > 1
        if isstruct(cond)
            data = cond;
        else
            data = load_data(cond);
        end
    end
    
    switch mode
        
        case 'results'
            figure;
            subplot(2,2,1);
            analyze_data('update',cond);
            mytitle('A','Left','FontSize',25,'FontWeight','Bold');
            subplot(2,2,2);
            analyze_data('update_back',cond);
            mytitle('B','Left','FontSize',25,'FontWeight','Bold');
            subplot(2,2,3);
            analyze_data('cross_update',cond);
            mytitle('C','Left','FontSize',25,'FontWeight','Bold');
            subplot(2,2,4);
            analyze_data('cross_update_back',cond);
            mytitle('D','Left','FontSize',25,'FontWeight','Bold');
            set(gcf,'Position',[200 200 750 800]);
            
        case 'results_glm'
            figure;
            subplot(1,2,1);
            analyze_data('update_glm',cond);
            mytitle('A','Left','FontSize',25,'FontWeight','Bold');
            subplot(1,2,2);
            analyze_data('update_glm_back',cond);
            mytitle('B','Left','FontSize',25,'FontWeight','Bold');
            set(gcf,'Position',[200 200 900 400]);
            
        case 'update_glm'
            for s = 1:length(data)
                d = abs(diff(data(s).p(1:2:end)));
                lr = abs(diff(data(s).p(2:2:end)));
                ix = ~isnan(lr) & ~isnan(d);
                d = d(ix); lr = lr(ix);
                K = 5;
                x = ones(length(d)-K,1);
                for k=1:K; x = [x d((K-k)+1:end-k)]; end
                b(s,:) = regress(lr(K+1:end),x);
            end
            b = b(:,2:end);
            m = mean(b);
            se = std(b)./sqrt(length(data));
            [~,p] = ttest(b)
            xlim = [0 K+1];
            errorbar(m,se,'-ok','MarkerFaceColor','k','MarkerSize',12,'LineWidth',4);
            hold on;
            plot(xlim,[0 0],'--k','LineWidth',3);
            set(gca,'FontSize',25,'XTick',1:5,'XLim',xlim,'YLim',[-0.05 0.15])
            xlabel('Trials back','FontSize',25);
            ylabel('Regression coefficient (a.u.)','FontSize',25);
            
        case 'update_glm_back'
            for s = 1:length(data)
                d = abs(diff(data(s).p(2:2:end)));
                lr = abs(diff(data(s).p(1:2:end)));
                ix = ~isnan(lr) & ~isnan(d);
                d = d(ix); lr = lr(ix);
                K = 5;
                x = ones(length(d)-K,1);
                for k=1:K; x = [x d(k:end-(K-k)-1)]; end
                b(s,:) = regress(lr(K+1:end),x);
            end
            b = fliplr(b(:,2:end));
            m = mean(b);
            se = std(b)./sqrt(length(data));
            [~,p] = ttest(b)
            xlim = [0 K+1];
            errorbar(m,se,'-ok','MarkerFaceColor','k','MarkerSize',12,'LineWidth',4);
            hold on;
            plot(xlim,[0 0],'--k','LineWidth',3);
            set(gca,'FontSize',25,'XTick',1:5,'XLim',xlim,'YLim',[-0.05 0.15])
            xlabel('Trials forward','FontSize',25);
            ylabel('Regression coefficient (a.u.)','FontSize',25);
            
        case 'volatility'
            M = [1 4];
            for m = 1:2
                data = load_data(M(m));
                LRm = zeros(length(data),1);
                for s = 1:length(data)
                    d = abs(diff(data(s).p(1:2:end)));
                    ix = ~isnan(d);
                    d = d(ix);
                    LRm(s,1) = nanmean(d);
                end
                x(m) = mean(LRm);
                se(m) = std(LRm)./sqrt(length(data));
                X{m} = LRm;
            end
            [~,p,~,stat] = ttest2(X{1},X{2})
            barerrorbar(x',se'); colormap bone
            set(gca,'XLim',[0.5 2.5],'XTickLabel',{'Low' 'High'},'FontSize',25);
            ylabel('|\Delta_2|','FontSize',25);
            xlabel('Volatility (q_1)','FontSize',25);
            
        case 'update'
            for s = 1:length(data)
                d = abs(diff(data(s).p(1:2:end)));
                lr = abs(diff(data(s).p(2:2:end)));
                ix = ~isnan(lr) & ~isnan(d);
                %d = d(ix); lr = lr(ix);
                d = detrend(d(ix)); lr = detrend(lr(ix));
                m = nanmean(d);
                LR(s,:) = [nanmean(lr(d>m)) nanmean(lr(d<m))];
                LRm(s,1) = nanmean(d);
            end
            [se,m] = wse(LR);
            myerrorbar(m,se,'-k','LineWidth',3);
            set(gca,'FontSize',25,'YLim',[min(m)-2*max(se) max(m)+2*max(se)],'XLim',[0.5 2.5],'XTickLabel',{'Large' 'Small'})
            xlabel('|\Delta_1|','FontSize',25);
            ylabel('|\Delta_2|','FontSize',25);
            [~,p,~,stat] = ttest(LR(:,1),LR(:,2))
            d = LR(:,1)-LR(:,2);
            disp(['Cohen d = ',num2str(abs(mean(d))./std(d))]);
            
        case 'update_back'
            for s = 1:length(data)
                d = abs(diff(data(s).p(2:2:end)));
                lr = abs(diff(data(s).p(1:2:end)));
                ix = ~isnan(lr) & ~isnan(d);
                %d = d(ix); lr = lr(ix);
                d = detrend(d(ix)); lr = detrend(lr(ix));
                m = nanmean(d);
                LR(s,:) = [nanmean(lr(d>m)) nanmean(lr(d<m))];
            end
            [se,m] = wse(LR);
            myerrorbar(m,se,'-k','LineWidth',3);
            set(gca,'FontSize',25,'YLim',[min(m)-2*max(se) max(m)+2*max(se)],'XLim',[0.5 2.5],'XTickLabel',{'Large' 'Small'})
            xlabel('|\Delta_2|','FontSize',25);
            ylabel('|\Delta_1|','FontSize',25);
            [~,p,~,stat] = ttest(LR(:,1),LR(:,2))
            d = LR(:,1)-LR(:,2);
            disp(['Cohen d = ',num2str(abs(mean(d))./std(d))]);
            
        case 'cross_update'
            for s = 1:length(data)
                d = diff(data(s).p(1:2:end));
                lr = diff(data(s).p(2:2:end));
                ix = ~isnan(lr) & ~isnan(d);
                %d = d(ix); lr = lr(ix);
                d = detrend(d(ix)); lr = detrend(lr(ix));
                m1 = nanmean(d(d<0)); m2 = nanmean(d(d>0));
                LR2(s,:,:) = [nanmean(lr(d<m1&lr>0)) nanmean(lr(d>m1&d<m2&lr>0)); nanmean(lr(d>m2&lr<0)) nanmean(lr(d<m2&d>m1&lr<0))];
            end
            c = [1 1 -1 -1];
            C = abs(LR2(:,:))*c';
            [~,p,~,stat] = ttest(C)
            [se,m] = wse(LR2);
            myerrorbar(m',se','-','LineWidth',3);
            set(gca,'FontSize',25,'YLim',[min(m(:))-2*max(se(:)) max(m(:))+2*max(se(:))],'XLim',[0.5 2.5],'XTickLabel',{'Large' 'Small'})
            xlabel('|\Delta_1|','FontSize',25);
            ylabel('\Delta_2','FontSize',25);
            legend({'\Delta_1<0, \Delta_2>0' '\Delta_1>0, \Delta_2<0'},'FontSize',20,'Location','East');
            
        case 'cross_update_back'
            for s = 1:length(data)
                d = diff(data(s).p(2:2:end));
                lr = diff(data(s).p(1:2:end));
                ix = ~isnan(lr) & ~isnan(d);
                %d = d(ix); lr = lr(ix);
                d = detrend(d(ix)); lr = detrend(lr(ix));
                m1 = nanmean(d(d<0)); m2 = nanmean(d(d>0));
                LR2(s,:,:) = [nanmean(lr(d<m1&lr>0)) nanmean(lr(d>m1&d<m2&lr>0)); nanmean(lr(d>m2&lr<0)) nanmean(lr(d<m2&d>m1&lr<0))];
            end
            c = [1 1 -1 -1];
            C = abs(LR2(:,:))*c';
            [~,p,~,stat] = ttest(C)
            [se,m] = wse(LR2);
            myerrorbar(m',se','-','LineWidth',3);
            set(gca,'FontSize',25,'YLim',[min(m(:))-2*max(se(:)) max(m(:))+2*max(se(:))],'XLim',[0.5 2.5],'XTickLabel',{'Large' 'Small'})
            xlabel('|\Delta_2|','FontSize',25);
            ylabel('\Delta_1','FontSize',25);
            legend({'\Delta_2<0, \Delta_1>0' '\Delta_2>0, \Delta_1<0'},'FontSize',20,'Location','East');
            
        case 'corr'
            for s = 1:length(data)
                p1 = data(s).p(1:2:end);
                p2 = data(s).p(2:2:end);
                ix = ~isnan(p1) & ~isnan(p2);
                r(s,1) = corr(p1(ix),p2(ix));
            end
            hist(r); colormap linspecer
            set(gca,'FontSize',25);
            ylabel('Frequency','FontSize',25);
            xlabel('Correlation','FontSize',25);
            nanmedian(r)
            p = signrank(r)
            
        case 'bic'
            [~,~,model1] = model_sim(1:3,1);
            [~,~,model2] = model_sim(1:3,2);
            for i = 1:length(data)
                N = sum(~isnan(data(i).p));
                bic(i,1) = N*log(nanmean((model1(i).p-data(i).p).^2)) + log(N);
                bic(i,2) = N*log(nanmean((model2(i).p-data(i).p).^2)) + log(N);
            end
            mean(bic(:,2)-bic(:,1))
            [~,p] = ttest(bic(:,1),bic(:,2))
            [~,~,xp] = bms(-bic)
            
        case 'model_scatter'
            [~,~,model] = model_sim(1:3,1);
            x = linspace(0,100,20);
            for i = 1:length(x)-1
                for j = 1:length(data)
                    ix = model(j).p>=x(i) & model(j).p<x(i+1);
                    p(j,i) = nanmean(data(j).p(ix));
                end
            end
            [se,m] = wse(p);
            x = x(1:end-1)+diff(x)/2;
            errorbar(x,m,se,'ok','LineWidth',4,'MarkerFaceColor','w','MarkerSize',10);
            set(gca,'FontSize',25,'XLim',[-1 101],'YLim',[-1 101],'XTick',0:20:100,'YTick',0:20:100);
            ylabel('Data','FontSize',25);
            xlabel('Model','FontSize',25);
            axis square
    end