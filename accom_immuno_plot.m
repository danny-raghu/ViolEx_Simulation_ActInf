function accom_immuno_plot(MDP)

Ng = numel(MDP(1).A);
Ne    = size(MDP(1).V,1) + 1;
Nt    = length(MDP);

for i=1:length(MDP)
    
    act_prob(:,i) = MDP(i).P(:,:,2)';
    act(:,i) = MDP(i).u(2,2);
    out(:,i) = MDP(i).o(2,:)';
    
    p(i)  = 0;
    for g = 1:Ng
        try
            U = spm_softmax(MDP(i).C{g});
        catch
            U = spm_softmax(MDP(i).C);
        end
        for t = 1:Ne
            p(i) = p(i) + log(U(MDP(i).o(g,t),t))/Ne;
        end
    end
    q(i)   = sum(MDP(i).rt(2:end));
end

%% 
%Actions and Action probs

col   = {'r.','g.','b.','c.','m.','k.'};
subplot(3,1,1)
if Nt < 64
    MarkerSize = 16;
else
    MarkerSize = 12;
end

colormap gray

image(256*(1 - act_prob)),  hold on;

plot(act,col{4},'MarkerSize',MarkerSize)

title('Action selection and action probabilities')
xlabel('Trial')
ylabel('Action') 
hold off
set(gca, 'YTick', [1:1:4], 'YTickLabel', {'Start','Hint','Choose Left','Choose Right'})
set(gca,'FontSize',9)
%axis off;

%% 

%Wins and Losses

subplot(3,1,2), bar(p,'k'),   hold on

for i = 1:size(out,2)
    if MDP(i).o(3,2) == 2
        j(i,1) = MDP(i).o(2,3)-1;
    else
        j(i,1) = MDP(i).o(2,2)-1;
    end
    if j(i,1) == 1
        jj(i,1) = 1;
    else
        jj(i,1) = -2;
    end
end

plot((j),col{2},'MarkerSize',MarkerSize);
plot((jj),col{6},'MarkerSize',MarkerSize);


title('Win/Loss and Free energies')
ylabel('Value and Win/Loss') 
spm_axis tight, hold off, box off
%axis off;
set(gca,'FontSize',9)
set(gca,'YTick',[-4:1:3])
yticklabels({'','','','Free Energy','','Loss','Win'})

%% 

%Priors and posteriors

subplot(3,1,3), hold on

%fig3 = figure('color','w');
%hold on

%xlim([0.5 (num_trials+0.5)]);
%ylim([0 1])
%ax2 = gca;


[feedback_traj,card_traj,advised_card,context] = extract_behaviour(MDP);

num_trials = length(MDP);

prior_card      =  zeros(1,num_trials);
posterior_card  =  zeros(1,num_trials);
switching = zeros(1,num_trials);
correct_card_inf=  zeros(1,num_trials);

max_pos = length(MDP);
max_total = round(max_pos * 1.5);
sep = 10;

card_1_col        = '#104c91';
card_2_col        = '#efc9af';
posterior_col     = 'k';
poor_inf_col      = '#bdfff6';

for trial = 1:num_trials
    cur_trial = MDP(trial);
    outcomes   = cur_trial.o;
    policy_prior_left = cur_trial.un(1,32);
    policy_prior_right = cur_trial.un(2,32);
    inf_over_card = cur_trial.xn{1,1};
    
    if(context(1,trial) == 1) && (outcomes(3,3) == 3)
        if(feedback_traj(1,trial) == 4)
            switching(1,trial) = 1;
        else
            switching(1,trial) = 0;
        end
    elseif(context(1,trial) == 2) && (outcomes(3,3) == 4)
        if(feedback_traj(1,trial) == 4)
            switching(1,trial) = 1;
        else
            switching(1,trial) = 0;
        end
    end
    
    prior_card(1,trial) = (inf_over_card(end,1,3,2)*policy_prior_left) / ...
                          (inf_over_card(end,2,3,2)*policy_prior_right + ...
                          inf_over_card(end,1,3,2)*policy_prior_left) ;
    
    posterior_card(1,trial) = (inf_over_card(end,1,3,3)*prior_card(1,trial)/...
                                                inf_over_card(end,1,3,2)) / ...
                                 ((inf_over_card(end,1,3,3)*prior_card(1,trial)/...
                                                inf_over_card(end,1,3,2)) + ...
                                   (inf_over_card(end,2,3,3)*(1-prior_card(1,trial))/...
                                                inf_over_card(end,2,3,2)));               
    
                                            
    %posterior_card(1,trial) = inf_over_card(end,1,3,3);
    
     if (context(1,trial) == 1) && (posterior_card(1,trial) > 0.5)
         correct_card_inf(1,trial) = 1;
     elseif (context(1,trial) == 2) && (posterior_card(1,trial) < 0.5)
         correct_card_inf(1,trial) = 1;
     elseif (switching(1,trial) == 1)   
         correct_card_inf(1,trial) = 1;
     else
         correct_card_inf(1,trial) = 0;
     end
end

%xlim([0.5 (num_trials+0.5)]);
%ylim([0 1])
%ax2 = gca;

for trial = 1:num_trials
    if context(1,trial) == 1
        trustworthy_colour = card_1_col;
    else
        trustworthy_colour = card_2_col;
    end
    rectangle('Position',[trial-0.5,0,1,1],...
        'FaceColor',trustworthy_colour,'EdgeColor','none');
    if correct_card_inf(1,trial) == 0
        rectangle('Position',[trial-0.5,0.5,1,0.5],...
            'FaceColor',poor_inf_col,'EdgeColor','none');
    end
end

plot(1:num_trials,prior_card,'color','k','LineWidth',2);
plot(1:num_trials,posterior_card,'.','color',posterior_col,'MarkerSize',10);

title('Tracking Beliefs')
ylabel('Probability') 
spm_axis tight;
set(gca,'YTick',[0:0.25:1],'YTickLabel', [0:0.25:1])
set(gca, 'Layer', 'Top')
set(gca,'FontSize',9)
set(gca,'TickLength',[0.005, 0.01])


end
%% 


