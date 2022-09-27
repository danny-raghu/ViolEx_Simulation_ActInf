clear all
close all 

%For exact replication of figures in our paper:

rs1 = 4.5; % Risk-seeking parameter: 4 for RAverse, 4.5 for RSeeking

% Simulation options after model building below:


% If Sim = 1, Simulate multiple trials where the left context is active.
             
% If Sim = 2, Simulate reversal learning, where the left context is active 
            % in early trials and then reverses in later 
            % trials.
            
            
Sim = 1;

%Change this to 1 if hint learning is allowed. Exhibited in the paper in
%3rd plot.
hint_learning = 0;

if rs1==4
    rand_seed=1;
elseif rs1==4.5
    rand_seed=2;
end

rng default



%% 1. Set up model structure

% Number of time points or 'epochs' within a trial: T
% =========================================================================

% Here, we specify 3 time points (T), in which the agent 1) starts in a 'Start'
% state, 2) first moves to either a 'Hint' state or a 'Choose Left' or 'Choose
% Right' slot machine state, and 3) either moves from the Hint state to one
% of the choice states or moves from one of the choice states back to the
% Start state.

T = 3;

% Priors about initial states: D and d
% =========================================================================

%--------------------------------------------------------------------------
% Specify prior probabilities about initial states in the generative 
% process (D)
% Note: By default, these will also be the priors for the generative model
%--------------------------------------------------------------------------

% For the 'context' state factor, we can specify that the 'left better' context 
% (i.e., where the left slot machine is more likely to win) is the true context:

D{1} = [1 0]';  % {'left better','right better'}

% For the 'behavior' state factor, we can specify that the agent always
% begins a trial in the 'start' state (i.e., before choosing to either pick
% a slot machine or first ask for a hint:

D{2} = [1 0 0 0]'; % {'start','hint','choose-left','choose-right'}

%--------------------------------------------------------------------------
% Specify prior beliefs about initial states in the generative model (d)
% Note: This is optional, and will simulate learning priors over states 
% if specified.
%--------------------------------------------------------------------------

% For context beliefs, we can specify that the agent starts out believing 
% that both contexts are equally likely, but with somewhat low confidence in 
% these beliefs:

d{1} = [.25 .25]';  % {'left better','right better'}

% For behavior beliefs, we can specify that the agent expects with 
% certainty that it will begin a trial in the 'start' state:

d{2} = [1 0 0 0]'; % {'start','hint','choose-left','choose-right'}


% State-outcome mappings and beliefs: A and a
% =========================================================================

%--------------------------------------------------------------------------
% Specify the probabilities of outcomes given each state in the generative 
% process (A)
% This includes one matrix per outcome modality
% Note: By default, these will also be the beliefs in the generative model
%--------------------------------------------------------------------------

% First we specify the mapping from states to observed hints (outcome
% modality 1). Here, the rows correspond to observations, the columns
% correspond to the first state factor (context), and the third dimension
% corresponds to behavior. Each column is a probability distribution
% that must sum to 1.

% We start by specifying that both contexts generate the 'No Hint'
% observation across all behavior states:

Ns = [length(D{1}) length(D{2})]; % number of states in each state factor (2 and 4)

for i = 1:Ns(2) 

    A{1}(:,:,i) = [1 1; % No Hint
                   0 0; % Machine-Left Hint
                   0 0];% Machine-Right Hint
end

% Then we specify that the 'Get Hint' behavior state generates a hint that
% either the left or right slot machine is better, depending on the context
% state. In this case, the hints are accurate with a probability of pHA. 

pHA = 0.9; %The Hint specifier is 100% accurate. Not probabilistic.

A{1}(:,:,2) = [0     0;      % No Hint
               pHA 1-pHA;    % Machine-Left Hint
               1-pHA pHA];   % Machine-Right Hint

% Next we specify the mapping between states and wins/losses. The first two
% behavior states ('Start' and 'Get Hint') do not generate either win or
% loss observations in either context:

for i = 1:2

    A{2}(:,:,i) = [1 1;  % Null
                   0 0;  % Loss
                   0 0]; % Win
end
           
% Choosing the left machine (behavior state 3) generates wins with
% probability pWin, which differs depending on the context state (columns):

pWin = 0.9; %Only probabilistic winning even if in the right context
           
A{2}(:,:,3) = [0      0;     % Null        
               1-pWin pWin;  % Loss
               pWin 1-pWin]; % Win

% Choosing the right machine (behavior state 4) generates wins with
% probability pWin, with the reverse mapping to context states from 
% choosing the left machine:
           
A{2}(:,:,4) = [0      0;     % Null
               pWin 1-pWin;  % Loss
               1-pWin pWin]; % Win
           
% Finally, we specify an identity mapping between behavior states and
% observed behaviors, to ensure the agent knows that behaviors were carried
% out as planned. Here, each row corresponds to each behavior state.
           
for i = 1:Ns(2) 

    A{3}(i,:,i) = [1 1];

end

%--------------------------------------------------------------------------
% Specify prior beliefs about state-outcome mappings in the generative model 
% (a)
% Note: This is optional, and will simulate learning state-outcome mappings 
% if specified.
%--------------------------------------------------------------------------
           
% One simple way to set up a matrix is by:
 
% 1. initially identifying it with the generative process 
% 2. multiplying the values by a large number to prevent learning all
%    aspects of the matrix (so the shape of the distribution changes very slowly)
% 3. adjusting the elements you want to differ from the generative process.

% To simulate learning the reward probabilities, we could specify:
    
    % a{1} = A{1}*200;
    % a{2} = A{2}*200;
    % a{3} = A{3}*200;
    % 
    % a{2}(:,:,3) =  [0  0;  % Null        
    %                .5 .5;  % Loss
    %                .5 .5]; % Win
    % 
    % 
    % a{2}(:,:,4) =  [0  0;  % Null        
    %                .5 .5;  % Loss
    %                .5 .5]; % Win

% We simulate learning the hint accuracy:
%Comment out the next few lines for Sim 1&2 of both conditions

  if (hint_learning == 1)
      a{1} = A{1}*200;
    a{2} = A{2}*200;
    a{3} = A{3}*200;
     
    a{1}(:,:,2) =  [0     0;     % No Hint
                   .25   .25;    % Machine-Left Hint
                   .25   .25];   % Machine-Right Hint
  else
      
  end
    
% Controlled transitions and transition beliefs : B{:,:,u} and b(:,:,u)
%==========================================================================

%--------------------------------------------------------------------------
% Next, we have to specify the probabilistic transitions between hidden states
% under each action (sometimes called 'control states'). 
% Note: By default, these will also be the transitions beliefs 
% for the generative model
%--------------------------------------------------------------------------

% Columns are states at time t. Rows are states at t+1.

% The agent cannot control the context state, so there is only 1 'action',
% indicating that contexts remain stable within a trial:

B{1}(:,:,1) = [1 0;  % 'Left Better' Context
               0 1]; % 'Right Better' Context
           
% The agent can control the behavior state, and we include 4 possible 
% actions:

% Move to the Start state from any other state
B{2}(:,:,1) = [1 1 1 1;  % Start State
               0 0 0 0;  % Hint
               0 0 0 0;  % Choose Left Machine
               0 0 0 0]; % Choose Right Machine
           
% Move to the Hint state from any other state
B{2}(:,:,2) = [0 0 0 0;  % Start State
               1 1 1 1;  % Hint
               0 0 0 0;  % Choose Left Machine
               0 0 0 0]; % Choose Right Machine

% Move to the Choose Left state from any other state
B{2}(:,:,3) = [0 0 0 0;  % Start State
               0 0 0 0;  % Hint
               1 1 1 1;  % Choose Left Machine
               0 0 0 0]; % Choose Right Machine

% Move to the Choose Right state from any other state
B{2}(:,:,4) = [0 0 0 0;  % Start State
               0 0 0 0;  % Hint
               0 0 0 0;  % Choose Left Machine
               1 1 1 1]; % Choose Right Machine        
           
%--------------------------------------------------------------------------
% Specify prior beliefs about state transitions in the generative model
% (b). This is a set of matrices with the same structure as B.
% Note: This is optional, and will simulate learning state transitions if 
% specified.
%--------------------------------------------------------------------------
          
% For this example, we will not simulate learning transition beliefs. 
% But, similar to learning d and a, this just involves accumulating
% Dirichlet concentration parameters. Here, transition beliefs are updated
% after each trial when the agent believes it was in a given state at time
% t and and another state at t+1.

% Preferred outcomes: C and c
%==========================================================================

%--------------------------------------------------------------------------
% Next, we have to specify the 'prior preferences', encoded here as log
% probabilities. 
%--------------------------------------------------------------------------

% Negative values indicate lower preference,
% positive values indicate a high preference. Stronger preferences promote
% risky choices and reduced information-seeking.

% We can start by setting a 0 preference for all outcomes:

No = [size(A{1},1) size(A{2},1) size(A{3},1)]; % number of outcomes in 
                                               % each outcome modality

C{1}      = zeros(No(1),T); % Hints
C{2}      = zeros(No(2),T); % Wins/Losses
C{3}      = zeros(No(3),T); % Observed Behaviors

% Then we can specify a 'loss aversion' magnitude (la) at time points 2 
% and 3, and a 'reward seeking' (or 'risk-seeking') magnitude (rs). Here,
% rs is divided by 2 at the third time point to encode a smaller win ($2
% instead of $4) if taking the hint before choosing a slot machine.

la = 1; % By default we set this to 1.

rs = rs1; % We set this value at the top of the script.

C{2}(:,:) =    [0  0   0   ;  % Null
                0 -la -la  ;  % Loss
                0  rs  rs/2]; % win
            
% Expanded out, this means that the other C-matrices will be:

% C{1} =      [0 0 0;     % No Hint
%              0 0 0;    % Machine-Left Hint
%              0 0 0];   % Machine-Right Hint
% 
% C{3} =      [0 0 0;  % Start State
%              0 0 0;  % Hint
%              0 0 0;  % Choose Left Machine
%              0 0 0]; % Choose Right Machine

            
%--------------------------------------------------------------------------
% One can also optionally choose to simulate preference learning by
% specifying a Dirichlet distribution over preferences (c). 
%--------------------------------------------------------------------------

% This will not be simulated here. However, this works by increasing the
% preference magnitude for an outcome each time that outcome is observed.
% The assumption here is that preferences naturally increase for entering
% situations that are more familiar.

% Allowable policies: U or V. 
%==========================================================================

%--------------------------------------------------------------------------
% Each policy is a sequence of actions over time that the agent can 
% consider. 
%--------------------------------------------------------------------------

% Policies can be specified as 'shallow' (looking only one step
% ahead), as specified by U. Or policies can be specified as 'deep' 
% (planning actions all the way to the end of the trial), as specified by
% V. Both U and V must be specified for each state factor as the third
% matrix dimension. This will simply be all 1s if that state is not
% controllable.

% For example, specifying U could simply be:

    % Np = 4; % Number of policies
    % Nf = 2; % Number of state factors
    % 
    % U         = ones(1,Np,Nf);
    % 
    % U(:,:,1) = [1 1 1 1]; % Context state is not controllable
    % U(:,:,2) = [1 2 3 4]; % All four actions in B{2} are allowed

% For our simulations, we will specify V, where rows correspond to time 
% points and should be length T-1 (here, 2 transitions, from time point 1
% to time point 2, and time point 2 to time point 3):

Np = 5; % Number of policies
Nf = 2; % Number of state factors

V         = ones(T-1,Np,Nf);

V(:,:,1) = [1 1 1 1 1;
            1 1 1 1 1]; % Context state is not controllable

V(:,:,2) = [1 2 2 3 4;
            1 3 4 1 1];
        
% For V(:,:,2), columns left to right indicate policies allowing: 
% 1. staying in the start state 
% 2. taking the hint then choosing the left machine
% 3. taking the hint then choosing the right machine
% 4. choosing the left machine right away (then returning to start state)
% 5. choosing the right machine right away (then returning to start state)


% Habits: E and e. 
%==========================================================================

%--------------------------------------------------------------------------
% Optional: a column vector with one entry per policy, indicating the 
% prior probability of choosing that policy (i.e., independent of other 
% beliefs). 
%--------------------------------------------------------------------------

% We will not equip our agent with habits in our example simulations, 
% but this could be specified as follows if one wanted to include a
% strong habit to choose the 4th policy:

% E = [.1 .1 .1 .6 .1]';

% To incorporate habit learning, where policies become more likely after 
% each time they are chosen, one can also specify concentration parameters
% by specifying e. For example:

% e = [1 1 1 1 1]';

% Yet to explore in depth - Shwartenbeck et al. code 
%==========================================================================

% Eta: learning rate (0-1) controlling the magnitude of concentration parameter
% updates after each trial (if learning is enabled).

    eta = 0.5; % By default we here set this to 0.5, but try changing its value  
               % to see how it affects model behavior

% Omega: forgetting rate (0-1) controlling the reduction in concentration parameter
% magnitudes after each trial (if learning is enabled). This controls the
% degree to which newer experience can 'over-write' what has been learned
% from older experiences. It is adaptive in environments where the true
% parameters in the generative process (priors, likelihoods, etc.) can
% change over time. A low value for omega can be seen as a prior that the
% world is volatile and that contingencies change over time.

    omega = 1; % By default we here set this to 1 (indicating no forgetting. 
               % Values below 1 indicate greater rates of forgetting.
               
% Beta: Expected precision of expected free energy (G) over policies (a 
% positive value, with higher values indicating lower expected precision).
% Lower values increase the influence of habits (E) and otherwise make
% policy selection less deteriministic. For our example simulations we will
% simply set this to its default value of 1:

     beta = 1; % By default this is set to 1, but try increasing its value 
               % to lower precision and see how it affects model behavior

% Alpha: An 'inverse temperature' or 'action precision' parameter that 
% controls how much randomness there is when selecting actions (e.g., how 
% often the agent might choose not to take the hint, even if the model 
% assigned the highest probability to that action. This is a positive 
% number, where higher values indicate less randomness. Here we set this to 
% a high value:

    alpha = 32;  % Any positive number. 1 is very low, 32 is fairly high; 
                 % an extremely high value can be used to specify
                 % deterministic action (e.g., 512)

% ERP: This parameter controls the degree of belief resetting at each 
% time point in a trial when simulating neural responses. A value of 1
% indicates no resetting, in which priors smoothly carry over. Higher
% values indicate degree of loss in prior confidence at each time step.

    erp = 1; % By default we here set this to 1, but try increasing its value  
             % to see how it affects simulated neural (and behavioral) responses
                          
% tau: Time constant for evidence accumulation. This parameter controls the
% magnitude of updates at each iteration of gradient descent. Larger values 
% of tau will lead to smaller updates and slower convergence time, 
% but will also promote greater stability in posterior beliefs. 

    tau = 12; % Here we set this to 12 to simulate smooth physiological responses,   
              % but try adjusting its value to see how it affects simulated
              % neural (and behavioral) responses
              
% Note: If these values are left unspecified, they are assigned default
% values when running simulations. These default values can be found within
% the spm_MDP_VB_X script (and in the spm_MDP_VB_X_tutorial script we
% provide in this tutorial).

% Other optional constants. 
%==========================================================================

% Chi: Occam's window parameter for the update threshold in deep temporal 
% models. In hierarchical models, this parameter controls how quickly
% convergence is 'cut off' during lower-level evidence accumulation. 
% specifically, it sets an uncertainty threshold, below which no additional 
% trial epochs are simulated. By default, this is set to 1/64. Smaller 
% numbers (e.g., 1/128) indicate lower uncertainty (greater confidence) is
% required before which the number of trial epochs are shortened.

% zeta: Occam's window for policies. This parameter controls the threshold
% at which a policy ceases to be considered if its free energy
% becomes too high (i.e., when it becomes too implausible to consider
% further relative to other policies). It is set to default at a value of 
% 3. Higher values indicate a higher threshold. For example, a value of 6
% would indicate that a greater difference between a given policy and the
% best policy before that policy was 'pruned' (i.e., ceased to be
% considered). Policies will therefore be removed more quickly with smaller
% zeta values.
         
% Note: The spm_MDP_VB_X function is also equipped with broader functionality
% allowing incorporation of mixed (discrete and continuous) models,
% plotting, simulating Bayesian model reduction during simulated
% rest/sleep, among others. We do not describe these in detail here, but
% are described in the documentation at the top of the function.

% True states and outcomes: s and o. 
%==========================================================================

%--------------------------------------------------------------------------
% Optionally, one can also specify true states and outcomes for some or all
% time points with s and o. If not specified, these will be 
% generated by the generative process. 
%--------------------------------------------------------------------------

% For example, this means the true states at time point 1 are left context 
% and start state:

    %      s = [1;
    %           1]; % the later time points (rows for each state factor) are 0s,
    %               % indicating not specified.
      

% And this means the observations at time point 1 are the No Hint, Null,
% and Start behavior observations.

    %      o = [1;
    %           1;
    %           1]; % the later time points (rows for each outcome modality) are 
    %               % 0s, indicating not specified
 
%% 2. Define MDP Structure
%==========================================================================
%==========================================================================

mdp.T = T;                    % Number of time steps
mdp.V = V;                    % allowable (deep) policies

mdp.A = A;                    % state-outcome mapping
mdp.B = B;                    % transition probabilities
mdp.C = C;                    % preferred states
mdp.D = D;                    % priors over initial states

mdp.d = d;                    % enable learning priors over initial states
    
mdp.eta = eta;                % learning rate
mdp.omega = omega;            % forgetting rate
mdp.alpha = alpha;            % action precision
mdp.beta = beta;              % expected precision of expected free energy over policies
mdp.erp = erp;                % degree of belief resetting at each timestep
mdp.tau = tau;                % time constant for evidence accumulation


if hint_learning == 1
    mdp.a = a;
else
    
end


% Note, here we are not including habits:

    % mdp.E = E;

% or learning other parameters:
    % mdp.a = a;                    
    % mdp.b = b;
    % mdp.c = c;
    % mdp.e = e;         

% or specifying true states or outcomes:

    % mdp.s = s;
    % mdp.o = o;
    
% or specifying other optional parameters (described above):

    % mdp.chi = chi;    % confidence threshold for ceasing evidence
                        % accumulation in lower levels of hierarchical models
    % mdp.zeta = zeta;  % occams window for ceasing to consider implausible
                        % policies
      
% We can add labels to states, outcomes, and actions for subsequent plotting:

label.factor{1}   = 'contexts';   label.name{1}    = {'left-better','right-better'};
label.factor{2}   = 'choice states';     label.name{2}    = {'start','hint','choose left','choose right'};
label.modality{1} = 'hint';    label.outcome{1} = {'null','left hint','right hint'};
label.modality{2} = 'win/lose';  label.outcome{2} = {'null','lose','win'};
label.modality{3} = 'observed action';  label.outcome{3} = {'start','hint','choose left','choose right'};
label.action{2} = {'start','hint','left','right'};
mdp.label = label;

clear beta
clear alpha
clear eta
clear omega
clear la
clear rs % We clear these so we can re-specify them in later simulations,
         % but is already stored in mdp for usage

%--------------------------------------------------------------------------
% Use a script to check if all matrix-dimensions are correct:
%--------------------------------------------------------------------------
mdp = spm_MDP_check(mdp);

%%
if Sim == 1
%% 4. Stationary Environment

% Next, we can expand the mdp structure to include multiple trials

N = 30; % number of trials

MDP = mdp;

[MDP(1:N)] = deal(MDP);

MDP = spm_MDP_VB_X_rand(MDP,rand_seed);

plot_experiment_assimilate(MDP); 

%accom_immuno_plot(MDP)
%subtitle('Risk Averse - Stationary Environment')


% We can again visualize simulated neural responses

%spm_figure('GetWin','Figure 3'); clf    % display behavior

elseif Sim == 3
%% 5. Non-stationary environment

N = 32; % number of trials (must be multiple of 8)

MDP = mdp;

[MDP(1:N)] = deal(MDP);

    for i = 1:N/8
        MDP(i).D{1}   = [1 0]'; % Start in the 'left-better' context for 
                                % early trials
    end

    for i = (N/8)+1:N
        MDP(i).D{1}   = [0 1]'; % Switch to 'right-better' context for 
                                % the remainder of the trials
    end
    
MDP = spm_MDP_VB_X_rand(MDP,rand_seed);

plot_experiment_assimilate(MDP);
subtitle('Risk-Seeking - Non stationary, Hint Uncertainty')

end
%% 
