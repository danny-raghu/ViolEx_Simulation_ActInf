
function [feedback_traj,card_traj,advised_card,context] = behaviour(mdp)

% Takes in and MDP output tile, returns the feedback that was recieved and which card was
% selected

num_trials = length(mdp);
card_traj     = zeros(1,num_trials);
feedback_traj    = zeros(1,num_trials);
advised_card = zeros(1,num_trials);
context  = zeros(1,num_trials);

for trial = 1:num_trials
    cur_trial = mdp(trial);
    cur_outcomes=cur_trial.o;
    cur_states   = cur_trial.s;
    feedback_traj(1,trial) = cur_outcomes(2,3);
    card_traj(1,trial)     = cur_outcomes(3,3);
    advised_card(1,trial)= cur_outcomes(1,2);
    context(1,trial) = cur_states(1,3);
    
end

