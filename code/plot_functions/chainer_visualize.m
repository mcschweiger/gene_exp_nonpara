function Gim = chainer_visualize(Gim,chain)

%% -------------------------------------------------------------------------
if isempty(Gim)
   
figure(10);
set(gcf,'windowstyle','docked');
subplot(2,2,[1 2 3 4]);

% LABELS AND AXES

GTrates     = chain.params.ground.rates;
                 
p1 = plot(1,chain.sample.rates',...
           '-o',...
          'MarkerEdgeColor','b',...
          'MarkerFaceColor',[0 .7 .0],...
          'MarkerSize',5);

line(xlim,GTrates.*[1 1],'linestyle','-','color','g','linewidth',2);

Gim.p1 = p1;


xlim_0   =  0;
xlim_end =  size(chain.rates,2);


xlim([xlim_0 xlim_end]);
ylim([0,max(max([chain.rates]))+1]);

line(xlim,GTrates.*[1 1],'linestyle','-','color','g','linewidth',2);

xlabel('Iteration');
ylabel('Rate (1/s)');

set(gca);

box on

else

for j = 1:chain.params.nparams
    set(Gim.p1(j),'XData',1:max(chain.i),'YData',chain.rates(j,1:max(chain.i)));
end
xlim_0   =  0;
xlim_end =  size(chain.rates,2);

xlim([xlim_0 xlim_end]);
ylim([0,max(max([chain.rates]))+1]);
        
end

drawnow
