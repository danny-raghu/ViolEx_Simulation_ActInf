clear all
close all   
dbstop if error


% Simulation options after model building below:

% If Sim = 1, Simulate multiple trials in stationary setting (i.e. only
% left is correct slot machine) with high reliability.

%e_dir = 600; pwin = 0.9;

% If Sim = 2, Simulate multiple trials in non-stationary setting (correct
% slot switches from left to right) w/o any habit learning and low reliability.

%e_dir = 600; pwin = 0.6;

% If Sim = 3, Simulate multiple trials in non-stationary setting (correct
% slot switches from left to right) with habit learning.

%e_dir = 2; pwin = 0.6;

% If Sim = 4, Simulate multiple trials in non-stationary setting (correct
% slot switches from left to right) with habit learning and high reliability.

%e_dir = 2; pwin = 0.9

pWin = 0.9; %This formalizes the reliability of evidence, ranging from 0.6 to 0.9 (low to high reliability).
e_dir = 2; %The constrains whether habits are learned or not. Higher values inhibit learning. 

Sim = 4;

 if (pWin == 0.9) && (e_dir == 2)
     rand_seed=1;
 else
     rand_seed=2;
 end

rng default


load('trust_sequence_250.mat') %For Sim 1
load('switch_sequence_250.mat') % For Sim 2 and 3.


%% 1. Set up model structure

% Number of time points or 'epochs' within a trial: T
% =========================================================================

% Here, we specify 3 time points (T), in which the agent 1) starts in a 'Start'
% state, 2) moves to a 'Hint' state, and 3) moves from the Hint state to one
% of the choice states.

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

D{1} = [1 1]';  % {'left better','right better'}

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

%d{1} = [1 1]';  % {'left better','right better'}
%d{1} = d{1}*200;

%The agent doesn't know what is the right context, but cannot learn it
%either. 

% For behavior beliefs, we can specify that the agent expects with 
% certainty that it will begin a trial in the 'start' state:

%d{2} = [1 0 0 0]'; % {'start','hint','choose-left','choose-right'}


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

pHA = 0.5; %The Hint specifier is 50% accurate. He is either trustworthy, or not in a trial. 

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

A{2}(:,:,3) = [0    0;    % Null        
               0    1;  % Loss
               1    0]; % Win

% Choosing the right machine (behavior state 4) generates wins with
% probability pWin, with the reverse mapping to context states from 
% choosing the left machine:
           
A{2}(:,:,4) = [0    0;     % Null
               1    0;  % Loss
               0    1]; % Win
          

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
    
     a{1} = A{1}*1200;
     a{3} = A{3}*1200;
     
     a{2}(:,:,3) =  [0      0;     % Null        
                     1-pWin pWin;  % Loss
                     pWin 1-pWin]; % Win

     
     
     a{2}(:,:,4) =  [0      0;     % Null
                     pWin 1-pWin;  % Loss
                     1-pWin pWin]; % Win
           
     a{2} = a{2}*1200;
                 
  %In our case, we are going to preclude learning of any sort. We are
  %changing the likelihood matrix
  
% Controlled transitions and transition beliefs : B{:,:,u} and b(:,:,u)
%==========================================================================

%--------------------------------------------------------------------------
% Next, we have to specify the probabilistic transitions between hidden states
% under each action (sometimes called 'control states'). 
% Note: By default, these will also be the transitions beliefs 
% for the generative model
%--------------------------------------------------------------------------

% Columns are states at time t. Rows are states at t+1.

% We introduce some randomness in the trials. Even if it at the start of
% the trial, the left slot was 'correct', it will change to right being correct 10% of the time.
% We can make sure that the agent doesn't know this by specifying an
% invariant 'b' (see below).

B{1}(:,:,1) = [0.9 0.1;  % 'Left Better' Context
               0.1 0.9]; % 'Right Better' Context
           
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
          
%Making sure that the agent doesn't know this randomness and precluding learning.

b{1}(:,:,1) = [1 0;  % 'Left Better' Context
               0 1]; % 'Right Better' Context
   
b{1} = b{1}*1200;

b{2} = B{2}*1200;

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


C{2}(:,:) =    [0  0   0;  % Null
                0  0  -2;  % Loss
                0  0   2]; % win
            
            
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

% For our simulations, we will specify V, where rows correspond to time 
% points and should be length T-1 (here, 2 transitions, from time point 1
% to time point 2, and time point 2 to time point 3):

Np = 2; % Number of policies
Nf = 2; % Number of state factors

V         = ones(T-1,Np,Nf);

V(:,:,1) = [1 1;
            1 1]; % Context state is not controllable

V(:,:,2) = [2 2;
            3 4];
        
% For V(:,:,2), columns left to right indicate policies allowing: 
% 1. staying in the start state 
% 2. taking the hint then choosing the left machine
% 3. taking the hint then choosing the right machine


% Habits: E and e. 
%==========================================================================

%--------------------------------------------------------------------------
% Optional: a column vector with one entry per policy, indicating the 
% prior probability of choosing that policy (i.e., independent of other 
% beliefs). 
%--------------------------------------------------------------------------

E = ones(2,1) .* e_dir;

%Single-valued parameters
%==========================================================================

% Beta: Expected precision of expected free energy (G) over policies (a 
% positive value, with higher values indicating lower expected precision).
% Lower values increase the influence of habits (E) and otherwise make
% policy selection less deteriministic. For our example simulations we will
% simply set this to its default value of 1:

     beta = 32; % By default this is set to 1, but try increasing its value 
               % to lower precision and see how it affects model behavior

% Alpha: An 'inverse temperature' or 'action precision' parameter that 
% controls how much randomness there is when selecting actions (e.g., how 
% often the agent might choose not to take the hint, even if the model 
% assigned the highest probability to that action. This is a positive 
% number, where higher values indicate less randomness. Here we set this to 
% a low value:

    alpha = 1.5;  % Any positive number. 1 is very low, 32 is fairly high; 
                 % an extremely high value can be used to specify
                 % deterministic action (e.g., 512)


% Other optional constants. 
%==========================================================================

% zeta: Occam's window for policies. This parameter controls the threshold
% at which a policy ceases to be considered if its free energy
% becomes too high (i.e., when it becomes too implausible to consider
% further relative to other policies). It is set to default at a value of 
% 3. Higher values indicate a higher threshold. For example, a value of 6
% would indicate that a greater difference between a given policy and the
% best policy before that policy was 'pruned' (i.e., ceased to be
% considered). Policies will therefore be removed more quickly with smaller
% zeta values.
         
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
      

%% 2. Define MDP Structure
%==========================================================================
%==========================================================================

mdp.T = T;                    % Number of time steps
mdp.V = V;                    % allowable (deep) policies  
mdp.A = A;                    % state-outcome mapping
mdp.a = a;
mdp.B = B;                    % transition probabilities
mdp.b = b;                    %   
mdp.C = C;                    % preferred states
mdp.D = D;                    % priors over initial states
% mdp.d = d;                    % 
mdp.e = E;                    %Habit learning

%Single-valued parameters
mdp.alpha = alpha;            % action precision
mdp.beta = beta;              % expected precision of expected free energy over policies
mdp.zeta = inf;

stationary = trust_sequence([1,5],1:50);
non_stationary = switch_sequence([1,5],[1:32,125:142]);

% or specifying true outcomes:

    % mdp.o = o;
    

% We can add labels to states, outcomes, and actions for subsequent plotting:

label.factor{1}   = 'contexts';   label.name{1}    = {'left-better','right-better'};
label.factor{2}   = 'choice states';     label.name{2}    = {'start','hint','choose left','choose right'};
label.modality{1} = 'hint';    label.outcome{1} = {'null','left hint','right hint'};
label.modality{2} = 'win/lose';  label.outcome{2} = {'null','lose','win'};
label.modality{3} = 'observed action';  label.outcome{3} = {'start','hint','choose left','choose right'};
label.action{2} = {'start','hint','left','right'};
mdp.label = label;

%clear pwin;
%clear e_dir; 
% We clear these so we can re-specify them in later simulations,
% but is already stored in mdp for usage

%--------------------------------------------------------------------------
% Use a script to check if all matrix-dimensions are correct:
%--------------------------------------------------------------------------
mdp = spm_MDP_check(mdp);

%%
if Sim ==1
%% 1. Multi trial simulations

N = size(stationary,2); % number of trials

MDP = mdp;

[MDP(1:N)] = deal(MDP);

for decided = 1:size(stationary,2)
    cur_state = stationary(:,decided);
    MDP(decided).s = cur_state;
end



MDP = spm_MDP_VB_X_rand(MDP, rand_seed);

%save('MDP','MDP_tutorial')

accom_immuno_plot(MDP)
subtitle('Accommodation (X = 90%, No Habits)')




elseif Sim == 2
%% 5. Simulate reversal learning

N = size(non_stationary,2); % number of trials

MDP = mdp;

[MDP(1:N)] = deal(MDP);

for decided = 1:size(non_stationary,2)
    cur_state = non_stationary(:,decided);
    MDP(decided).s = cur_state;
end



MDP = spm_MDP_VB_X_rand(MDP, rand_seed);

%save('MDP','MDP_tutorial')

accom_immuno_plot(MDP)
subtitle('No Immunization(X = 90%, Habit Learning)')



end


