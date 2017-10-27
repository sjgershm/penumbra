function data = load_data(cond)
    
    if length(cond) > 1
        data = [];
        for i = 1:length(cond)
            data = [data load_data(cond(i))];
        end
    else
        
        switch cond
            case 1
                fdir = fullfile('..','data');
            case 2
                fdir = fullfile('..','data2');
            case 3
                fdir = fullfile('..','data3');
            case 4
                fdir = fullfile('..','data4');
        end
        files = dir(fullfile(fdir,'*.csv'));
        
        n = 0;
        for i = 1:length(files)
            try
                d = csvread(fullfile(fdir,files(i).name));
            catch
                disp(files(i).name);
            end
            if size(d,1)==100
                n=n+1;
                data(n).C = zeros(size(d,1))+cond;
                data(n).x = d(:,1);
                data(n).t = d(:,2);
                data(n).y = d(:,3);
                data(n).m = d(:,4);
                data(n).p = d(:,5);
                data(n).p(data(n).p<=0|data(n).p>=100) = nan;
                err(n) = nanmean(abs(data(n).m-data(n).p));
                na(n) = sum(isnan(data(n).p));
                u(n) = length(unique(data(n).p));
                if u(n)<5
                    disp(files(i).name);
                end
            end
        end
        
        %data(na>25) = [];
        %data(err>20) = [];
    end