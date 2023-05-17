%%

addpath(genpath('../PolytopeSamplerMatlab-develop'));
addpath(genpath('../CHRR-test'));
maxNumCompThreads(1);
assert(maxNumCompThreads == 1);

curFolder = fileparts(mfilename('fullpath'));
instances = {'box', 'simplex', 'birkhoff'};


dimBoxSimplex = [10];
dimBirkhoff = [9 16];
numSamples = 1e6;
maxIter = 1e6;
nb_repet = 3;
thin = 1;
thinOutput=false;
idx = 2; %% change idx to 1 for box, 2 for simplex, 3 for birkhoff



for c = [idx]
    inst = instances{c};
    if strcmp('box', inst) || strcmp('simplex', inst)
        dimInput = dimBoxSimplex;
    elseif strcmp('birkhoff', inst)
        dimInput = dimBirkhoff;
    else
        error('Not supported')
    end

    for dim = dimInput
        rng(1); 
        stat_mid = [];
        stat_GL = [];
        stat_mala= [];
        fprintf('%s: %d\n', inst, dim);
        disp(datetime('now'));
        unit_vec = 1:dim;
        unit_vec(1) = 0;
        unit_vec = unit_vec'/norm(unit_vec);
        unit_vec(2) = 1;
        unit_vec = 20*unit_vec;
        f = @(x) sum((x-unit_vec).^2)/2;
        df = @(x) (x-unit_vec); 
        for i=1:nb_repet
            filename = fullfile(curFolder, strcat('bias_analysis_local/_', inst, int2str(dim),'_rep_',num2str(i),'.mat'));
            if isfile(filename)
                o= load(filename);
                samples_used = o.exps.total_samples;
                stat_final =sum((unit_vec')*samples_used)/length(samples_used(1,:));
                stat_mid = [stat_mid stat_final];
            end
            filename = fullfile(curFolder, strcat('bias_analysis_local/GL_', inst, int2str(dim),'_rep_',num2str(i),'.mat'));
            if isfile(filename)
                o= load(filename);
                samples_used = o.exps.total_samples;
                stat_final =sum((unit_vec')*samples_used)/length(samples_used(1,:));
                stat_GL = [stat_GL stat_final];
            end

            if c == 1 
                n_iter_mala = 1e6;
                x_0 = zeros(dim,1);
                h=0.05;
                [samples_mala,accept_ratio] = mala_box(f,df,x_0,h,n_iter_mala,dim);
                accept_ratio(1,end)
                trunc = 1e5;
                u = samples_mala(:,trunc+1:end);
                unit_vec_mat = unit_vec'*ones(dim,length(samples_mala(1,trunc+1:end)));
                cumU = sum(((unit_vec')*u))/length(u);
                stat_mala = [stat_mala cumU];
            end
            if (c == 2) && (dim <100)
                n_iter_imh = 1e6;
                x_0 = abs(unit_vec)/sum(abs(unit_vec));
                [samples_imh,accept_ratio] =  imh_simplex(f,x_0,n_iter_imh,dim);
                accept_ratio(1,end)
                trunc = 1;
                u = samples_imh(:,trunc+1:end);
                unit_vec_mat = unit_vec'*ones(dim,length(samples_imh(1,trunc+1:end)));
                cumU = sum(((unit_vec')*u))/length(u);
                stat_mala = [stat_mala cumU];
            end
        end
        figure(dim);
        g1 = repmat({'CRHMC'},length(stat_mid),1);
        g2 = repmat({'n-BHMC'},length(stat_GL),1);
        g3 = repmat({'IMH'},length(stat_mala),1);
        g = [g1; g2; g3];
        boxplot([stat_mid'; stat_GL'; stat_mala'],g);
        set(findobj(gca,'Type','text'),'FontSize',32,'fontWeight','bold') % to set Xaxis
        set(gca,'FontSize',18,'fontWeight','bold') % to set Yaxis
        labelx = strcat('$$d=',  int2str(dim) ,'$$');
        title(labelx,'Interpreter','latex', 'FontSize',32);
        saveas(gcf,strcat('boxplot_',inst,'_dim_',int2str(dim)),'png')
        saveas(gcf,strcat('boxplot_',inst,'_dim_',int2str(dim)),'fig')
    end

end

