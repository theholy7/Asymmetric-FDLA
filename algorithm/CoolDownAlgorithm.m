function [ weightMatrix ] = CoolDownAlgorithm( adjacencyMatrix)
%CoolDownAlgorithm Solves the Asymmetric FDLA optimization problem
%   Input adjacency matrix like:
% w =
%
%      1     0     0     1     1     1     1     0     0     0
%      0     1     0     0     1     0     0     0     0     1
%      0     1     1     0     0     0     0     0     0     1
%      1     0     0     1     0     0     0     0     1     0
%      1     1     1     0     1     0     0     1     0     1
%      1     0     0     0     0     1     1     0     0     0
%      1     0     0     0     1     1     1     0     1     0
%      0     0     0     0     1     0     0     1     0     0
%      1     0     0     0     0     0     1     0     1     0
%      0     1     0     0     0     0     0     0     0     1
%
%   For the graph represented by this adjacency matrix, the algorithm finds
%   the weight matrix that optimizes the communication weights between each
%   node, thus solving the FDLA problem for asymmetric graphs.

addpath('../requiredObjects/')
addpath('../requiredFunctions/')
global isAborted

%Inform user that previous files will be overritten
% Construct a questdlg with three options
choice = questdlg('Previous files will be overwritten. Proceed?', ...
	'Warning!');
% Handle response
switch choice
    case 'Yes'
        ;
    case 'No'
        error('Stopped by user');
    case 'Cancel'
        error('Stopped by user');
end

% Is the input value a matrix? DIM > 1?
if size(adjacencyMatrix) == 1
    error('Not a matrix')
end



% Is the matrix square? Needs to be.
% Dimension is nodes x nodes.
% Entries are: "1" if node i connected to j, "0" if not.
% This algorithm only makes sense with asymmetric matrices.

if ~( size(adjacencyMatrix) == size(adjacencyMatrix') )
    error('Adjacency Matrix is not square.')
end

%Is the matrix representing a connected graph? It needs to be, it is the
%only way that the info in any node can reach any other.

if ~isStrongly( adjacencyMatrix )
    error('Adjacency Matrix and graph are not strongly connected.')
end


% Now we can start the algorithm.
disp('Adjacency matrix is square and graph is strongly connected.')

% Step 1 of the algorithm is to minimize the Maximum Singular Value of W
% (the weight matrix).

% From the adjacency matrix get the number of nodes.
n = size(adjacencyMatrix, 1);



% I introduced an auxiliary variable, eigenOne, to be the vector of 1's
% whose dimension is (nodes x 1). It makes it easy to impose the doubly
% stochastic constraints of W.
eigenOne = ones(n, 1);



% Here we begin the optimization of W - find W with the smallest maximum
% singular value. We are now using the CVX package.
tic
% begin CVX environment
cvx_begin quiet

% define matrix of variables
variable W(n,n)

% define the problem: minimize maximum singluar value of W
minimize norm(W- (1/n * ones(n,n)))

% define the constraints on W as given by Boyd'd FDLA
subject to

%W equals its transpose
eigenOne'*W == eigenOne';
W*eigenOne == eigenOne;

%W entries that are 0 remain 0 - no physical link between nodes
W(~adjacencyMatrix)==0;

% end CVX environment and calculate W
cvx_end

% Round any entry (in abs value) smaller than 10^-6 to 0. They make no
% difference when calculating the FDLA
closeToZero = abs(W) <= 10^-6;
W(closeToZero) = 0;

% Save the resulting matrix W as an weightMatrixObj
wm = WeightMatrix(W);

% Print the status of the initial problem
fprintf('Initial problem status: %s\n', cvx_status)

% Print the iteration (W0), how long it took, and the Spectral Radius
% obtained via the initial minimization problem.
fprintf('W0 \t Time: %d \t Spectral Radius: %f\n', toc, wm.spectralRadius)


% Now we solve the second part of the algorithm, the iterative part, as
% described in my MSc Thesis or in the Manual of this program.


% Variable initiation for the algorithm
% We have an initial weight matrix, it comes from the 1st problem
Wmat = wm.weightMatrix; %ideal weight matrix

%this is the one shot matrix - if all nodes were connected, they would just
%need to calculate the average one time to reach the goal of the FDLA
Jmat = 1/n * ones(n,n); %one shot matrix

%This is the objective we are setting. We want to find a matrix with a
%spectral radius lower than eta
eta = wm.spectralRadius;%initial eta value

% this is the step we use to reduce eta at each iteration of the problem
epsilon = 10^-6;

% here we save the initial result in a results vector, for later plotting
results(1) = wm.spectralRadius;

% in the first step of this second part of the algorithm we need to
% initialize P and Q - matrices we are using in the convex approximation we
% have defined

% We are going to find the matrix P that solves que equality defined in the
% thesis with the smallest condition number. We want it to be stable for our
% computations.

% Solve for P, k

tic
% again we initialize the CVX environment in SDP mode
cvx_begin quiet sdp

%define variables to minimize
variable k
variable Pmat(n,n) symmetric

% define problem - minimize K
minimize k

subject to
% these are the constraints on P - small condition number
eye(n,n) <= Pmat <= k*eye(n,n)

% P needs to solve this inequality - It is explained in the MSc Thesis and
% manual
(Wmat - Jmat)'*Pmat*(Wmat - Jmat) <= (eta^2+epsilon)*Pmat

%close environment and calculate everything
cvx_end

% print CVX status of the problem of finding P - to ensure it ran properly
fprintf('Pzero status: %s\n', cvx_status)

%print time it took to calculate it - curiosity
fprintf('Pzero \t Time: %d\n', toc)

% Now we start the part of the problem where we are going to perform
% iterations and successive approximations to to the matrix with the best
% possible (lowest) spectral radius

% Set iterable variables
% P is also here - the variable Pmat will be overwritten at each iteration
% eta is the new objective for the spectral radius we want to hit or be
% under
eta = 0.95*eta;
% Qmat is another auxiliary variable - also overwritten at each iteration
Qmat = inv(Pmat);

%this flag determines when to stop the program
exitFlag = 0;


% Start iterations for the algorithm
for m = 1:20
    % Start iterations for the P-Q optimization
    for l = 1:20
        % check if stop button has been pressed
        % isAborted = true?
        pause(0.0001)
        %disp(isAborted)

        if isAborted
            error('Aborted by User');
        end

        % solve the problem of finding small variations to P and Q that
        % still respect the constraints on W
        tic
        %start the CVX environment
        cvx_begin quiet sdp
        % define the variables to work with
        variable W(n,n)
        variable deltaP(n,n) symmetric
        variable deltaQ(n,n) symmetric

        % minimize the trace of P*Q - the objective is a linearized version
        % of this
        minimize trace( (Pmat + deltaP)*Qmat + Pmat*(Qmat + deltaQ) )

        subject to

        %define the constraints on W and deltaP and deltaQ
        [ eta^2*(Pmat+deltaP) , (W-Jmat)';
            (W-Jmat), (Qmat+deltaQ) ] >= 0;

        [ (Pmat + deltaP) , eye(n,n);
            eye(n,n), (Qmat + deltaQ) ] >= 0;

        %define the constraints on W
        eigenOne'*W == eigenOne';
        W*eigenOne == eigenOne;
        W(~adjacencyMatrix)==0;

        %end cvx and calculate
        cvx_end

        % save the results as before, rounding to zero entreis smaller than
        % abs(10^-6)
        closeToZero = abs(W) <= 10^-6;
        W(closeToZero) = 0;

        %save the weight matrix as a weightMatrixObj
        wm = WeightMatrix(W);

        % again, print some information about the cvx status and program
        fprintf('P-Q Update status: %s \t Time: %d\n', cvx_status, toc);

        Pmat = Pmat + deltaP;
        Qmat = Qmat + deltaQ;

        % check if the problem was properly solved by CVX
        if strcmp(cvx_status, 'Solved') == 1 ...
                || strcmp(cvx_status, 'Inaccurate/Solved') == 1

            % set exit flag to zero - all is ok :)
            exitFlag = 0;

            % Check if P-Q optimization has stabilized - meaning that the
            % frobenius norm of the matrices is varying less than 1pc per
            % iteration - if it has, break the P-Q optimization cycle
            if (norm(deltaP,'fro')/norm(Pmat, 'fro') <= 0.01 && ...
                    norm(deltaQ,'fro')/norm(Qmat, 'fro') <= 0.01)

                %break the P-Q optimization cycle
                break
            end

        % if the problem wasn't properly solved
        else
            % set exitFlag to 1 and break P-Q cycle
            exitFlag = 1;
            break;
        end
    end

    % Since exitFlag is now = 1 (if P-Q wasn't properly solved) exit the
    % program completely
    if exitFlag == 1
        % break outer for-loop
        break;
    end

    % Everything was OK up to now - We have a matrix W that solves our
    % problem, it might not be the best yet, but it works
    % check if P and Q are OK and follow the constraints that we set :)
    % also check if the new spectral radius is lower than the previous
    if trace(Pmat*Qmat) <= n+1 && wm.spectralRadius <= results(m)

        %save matrix of m iteration
        filename = sprintf('weightMatrixAtIteration_%d.mat',m);

        %save(filename, 'wm');

        %save results to results vector for easy plotting
        results(m+1) = wm.spectralRadius;
        %write that vector as a CSV file
        csvwrite('weightMatrixSpectralRadius.csv', results);

        %set output of fucntion
        weightMatrix = wm;


        %traceResult = trace(Pmat*Qmat);
        %set a new goal for the spectral radius to be under
        eta = 0.95*eta;
        %print some more information
        fprintf('Result achieved: %f \t New goal: %d\n',results(m+1), eta)
    else
        break;
    end

end


end
