% Visualize results from the lowest temperature PT chain





close all; clear variables;
%%
set(0, 'DefaultLineLineWidth', 2);
set(0, 'defaultAxesFontSize', 14);
addpath ./plot_functions

%%
nchains = 1;        %indicate the number of chains you want to plot

runID = 'J10';      %indicate the runID of the desired result

save = false;       %(not reccomended, i.e. leave false) save each plotted figures as a .png

burn = 500;         %for easier interpretation, exclude some burn-in from plots

for i = 1:nchains   %load results
    try
  load(['./results_cont/',runID,'/','PTChain',num2str(i),'_cont','.mat']); 
    catch
  load(['../results/',runID,'/','PTChain',num2str(i),'.mat']); 
    end
  chains(i) = chain;
  clear chain
end

chains(1) = hungarian(chains(1),burn);
%% Save Files/Names--------------------------------------------------------
if save
    figpath = ['./IMAGES/5_27/',runID];
%     figpath = ['./IMAGES/8_20/',num2str(ground.rates(1:3),'%.0e_'),'/',num2str(ground.rates(7:9),'%.0e_'),'/'];  % G3 save
    if ~exist(figpath,'dir')
        mkdir(figpath)
    end 
end

%%
for h = 1:nchains
    rates     = 10.^chains(h).rates;
    loads     = chains(h).loads;
    sigma_ast = chains(h).sigma_ast;
    logpost   = chains(h).logpost;
    loglike   = chains(h).loglike;
    q         = chains(h).q;


    params = chains(h).params;

    burn     = floor(size(rates,2)*.3)+1;

    burn_rates = rates(:,burn:end);


    % Plotting----------------------------------------------------------------
    col = [0 1 0;%green
           0 0 1;%blue
           1 0 0];%red
    n_bins = 25;

    indices   = [nchoosek(1:params.L,2); flip(nchoosek(1:params.L,2),2)];

    figure
    sp1 = subplot(1,10,1:5);
    for i = 1:2*nchoosek(params.L,2)
        ind     = indices(i,:);
        bb      = loads(ind,:);
        load    = prod(bb);
        ON      = find(load);
        burn_ON = ON(ON > burn);

        rate_tmp = nan(size(rates(i,:)));
        rate_tmp(burn_ON) = rates(i,burn_ON);

        if ~isempty(burn_ON)
            kpp = plot(1:size(rate_tmp,2),rate_tmp,'LineWidth',2);
        end
        hold on
    end
    xlim([burn inf])
    title('Transition Rates')
    xlabel('Iteration')
    ylabel('(s^{-1})')

    % lg = legend(kp(lg_log),params.names(lg_log),'Orientation','Horizontal');
    % lg.NumColumns = 2;

    sp2 = subplot(1,10,6:10);
    ll = 1;
    for i = 1:2*nchoosek(params.L,2)
        ind  = indices(i,:);
        bb   = loads(ind,:);
        load = prod(bb);

        ON    = find(load);
        burn_ON = ON(find(ON > burn));

        if ~isempty(burn_ON)
            bp(ll) = histogram(rates(i,burn_ON),linspace(min(rates(i,burn_ON)),max(rates(i,burn_ON)),n_bins),...
                                            'Orientation','Horizontal',...
                                            'Normalization','pdf');
            lg_log(ll) = i;
            ll = ll+1;
        end
        hold on
    end
    yy = ylim;
    y  = linspace(yy(1),yy(2),100);
    prior = gampdf(y,params.phi(1),...
                    params.psi(1)/params.phi(1));
    pp = plot(prior,y,'Color','m');
    set(gca,'Ytick',[])

    linkaxes([sp1 sp2],'y')
    if exist('kp')
        legend([bp pp],[params.names(lg_log), {'Prior dist.'}])
    end

    figname = [ runID,'_BBrandinit_','trans_',...
                 '.png'];
    if save
        saveas(gcf,fullfile(figpath,figname))
    end


    %
    figure

    clear lg_log
    sp1 = subplot(1,2,1);
    k = 1;
    ll = 1;
    for i = 2*nchoosek(params.L,2)+1:size(rates,1)-1

        ON      = find(loads(k,:));
        burn_ON = ON(find(ON > burn));

        rate_tmp = nan(size(rates(i,:)));
        rate_tmp(burn_ON) = rates(i,burn_ON);

        if ~isempty(burn_ON)
            hh = plot(1:size(rate_tmp,2),rate_tmp,'LineWidth',2);
        end
        hold on
        k = k+1;
    end
    % legend(h([lg_log params.L+1]),params.names([lg_log+2*nchoosek(params.L,2) params.n_params]))
    xlim([burn inf])
    % ylim([0-.1 inf])
    title('Production Rates')
    xlabel('Iteration')

    sp2 = subplot(1,2,2);
    k = 1;
    for i = 2*nchoosek(params.L,2)+1:size(rates,1)-1

        ON      = find(loads(k,:));
        burn_ON = ON(find(ON > burn));

        if ~isempty(burn_ON)
            hb(ll) = histogram(rates(i,burn_ON),linspace(min(rates(i,burn_ON)),max(rates(i,burn_ON)),n_bins),'Orientation','Horizontal',...
                                          'Normalization','pdf');
            lgndB(ll) = i; 
            ll = ll+1;
        end
        hold on
        k = k+1;
    end
    yy = ylim;
    y  = linspace(yy(1),yy(2),100);
    prior = gampdf(y,params.phi(1),...
                    params.psi(1)/params.phi(1));
    pr = plot(prior,y,'Color','m');
    % ylim([0-.1 inf])
    set(gca,'Ytick',[])
    linkaxes([sp1 sp2],'y')
    legend([hb pr pp],[params.names(lgndB) {'Prior dist.'} {'GT'}])

    
    % Filename
    figname = [ runID,'_BBrandinit_','beta_',...
                 '_chain',num2str(h),...
                 '.png'];
    if save
        saveas(gcf,fullfile(figpath,figname))
    end

    % Degradation

    figure
    sp1 = subplot(1,10,1:5);
    h(k) = plot(burn:size(rates,2),rates(end,burn:end));
    xlim([burn inf])
    title('Degradation Rate')
    xlabel('Iteration')
    ylabel('\gamma (1/s)')

    sp2 = subplot(1,10,6:10);
    h(k) = histogram(rates(end,burn:end),linspace(min(rates(end,burn:end)),max(rates(end,burn:end)),n_bins),'Orientation','Horizontal',...
                                          'Normalization','pdf');
    hold on
    yy = ylim;
    y  = linspace(yy(1),yy(2),100);
    prior = gampdf(y,params.phi(end),...
                    params.psi(end)/params.phi(end));
    plot(prior,y,'Color','m')
    xlabel('Post. prob. dist.')
    set(gca,'Ytick',[])

    linkaxes([sp1 sp2],'y')

    % Filename
    figname = [ runID,'_BBrandinit_','gamma_',...
                 '.png'];
    if save
        saveas(gcf,fullfile(figpath,figname))
    end
    % Plot Loads
    figure('Units','inch','Position',[0 0 9 10]) % 8 width 7 height
    set(gcf,'color','w')
    movegui(gcf,'center')

    miniburn = burn;
    % 
    ax10 = subplot(6,1,1);
    plot(miniburn:size(logpost,2),logpost(miniburn:end))
    xlabel('Iteration')
    ylabel('Log Posterior')

    ax11 = subplot(6,1,2);
    plot(miniburn:size(loglike,2),loglike(miniburn:end))
    xlabel('Iteration')
    ylabel('Log Likelihood')


    ax2 = subplot(6,1,3);
    stairs(1:size(loads,2),sum(loads),'-o')
    % ylim([min(sum(loads))-.1 max(sum(loads))+.1])
    hold on
    % patch([2 2 41 41],[0 params.L params.L 0],'c','FaceAlpha',0.2)
    ylim([1 params.L])
    xlabel('Iteration')
    ylabel('Total # Active States')
    % legend({'Total # Active States','Adaptation Period'},'Location','NorthEast')

    ax3 = subplot(6,1,4);
    bar(sum(loads,2)/size(loads,2),1)
    ylabel('Active / Iter.')
    xlabel('Gene State')

    ax4 = subplot(6,1,5);
    histogram(sum(loads(:,burn:end)),'Normalization','pdf')
    xticks(1:params.L)
    xlabel('Total number gene states')
    ylabel('post. prob. dist.')

    ax5 = subplot(6,1,6);
    stairs(1:size(sigma_ast,2),sigma_ast,'-o')
    ylim([1 params.L+.5])
    xlabel('Iteration')
    ylabel('\sigma_{\ast}')

    linkaxes([ax3,ax4],'x');
    linkaxes([ax11,ax10,ax2,ax5],'x')

    % Filename
    figname = [ runID,'_BBrandinit_','loads_',...
                 '.png'];
    if save
        saveas(gcf,fullfile(figpath,figname))
        file = fullfile(figpath,figname);
    end

    figure('Units','inch','Position',[0 0 8 6]) % 8 width 7 height
    set(gcf,'color','w')
    movegui(gcf,'center')

    n_bins = 25;

    for i = 1:params.L
        histogram(q(i,:),linspace(min(q(i,burn:end)),max(q(i,burn:end)),n_bins),'Normalization','pdf')
        hold on
        lginf{i} = ['q_{',num2str(i),'}'];
    end
    lg = legend(lginf,'Orientation','horizontal','Position',[0.12760415608712,0.935185185185184,0.756944444444444,0.056712962962963]);

    box off
    xlabel('Post. prob. dist.')

    % linkaxes([sp1 sp2],'y')

    % Filename
    figname = [ runID,'_BB_','q_',...
                 '.png'];
    if save
        saveas(gcf,fullfile(figpath,figname))
    end

end

%% Logpost/loglike swapping

%%

figure('Units','inch','Position',[0 0 8 6]) % 8 width 7 height
set(gcf,'color','w')
movegui(gcf,'center')

for i = 1:nchains

ax10 = subplot(2,1,1);
plot(burn:size(chains(i).logpost,2),chains(i).logpost(burn:end))
hold on
xlabel('Iteration')
ylabel('Log Posterior')


ax11 = subplot(2,1,2);
plot(burn:size(chains(i).loglike,2),chains(i).loglike(burn:end))
hold on
xlabel('Iteration')
ylabel('Log Likelihood')

linkaxes([ax10 ax11],'x')
end

%%
figure('Units','inch','Position',[0 0 8 6]) % 8 width 7 height
set(gcf,'color','w')
movegui(gcf,'center')

for i = 1:nchains

ax1 = subplot(2,1,1);
stairs(1:size(chains(i).loads,2),sum(chains(i).loads),'-o')
hold on
ylim([1 chains(i).params.L+.5])
xlabel('Iteration')
ylabel('Total # Active States')

ax2 = subplot(2,1,2);
stairs(1:size(chains(i).sigma_ast,2),chains(i).sigma_ast,'-o')
hold on
ylim([1 chains(i).params.L+.5])
xlabel('Iteration')
ylabel('\sigma_{\ast}')

linkaxes([ax1 ax2],'x')

end

