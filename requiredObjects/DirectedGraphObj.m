classdef DirectedGraphObj
    %GraphObj Creates a Graph Object
    %   A GraphObj recieves the number of Nodes, maximum number of Edges
    %   and Node distance, and creates a graph respecting those
    %   conditions. It outputs a Laplacian matrix of the graph - L - and a
    %   vector with the initial values of each node - initX.
    
    properties (SetAccess = private)
        %inputs
        numNodes; %number of nodes in Graph
        numMaxEdges; %maximum number of node connections
        nodeDistance; %maximum distance of connection (increases until all nodes are connected)
    end
    properties (Hidden = true)
        %outputs
        adjacencyMatrix;
        graphConnection;
        initX;
        nodePositions;
        nodeNumbers;
    end
    
    methods
        %Constructor
        function obj = DirectedGraphObj(varargin)
            if nargin == 1 && isa(varargin{1}, 'DirectedGraphObj')
                %create a copy of a graph
                obj.numNodes = varargin{1}.numNodes; %n;
                obj.numMaxEdges = varargin{1}.numMaxEdges; %maxE;
                obj.nodeDistance = varargin{1}.nodeDistance; %nodeD;
                
                obj.initX = varargin{1}.getInitX;
                obj.nodePositions = varargin{1}.getNodePositions;
                obj.nodeNumbers = 1:obj.numNodes;
                obj.adjacencyMatrix = varargin{1}.adjacencyMatrix;
                
                obj.graphConnection = obj.adjacencyMatrix;
            end
            
            if nargin == 2 && isa(varargin{1}, 'GraphObj') && isa(varargin{2}, 'logical')
                %create a graph based on a new connection matrix
                assert(issymmetric(double(varargin{2})), 'Not a symmetric matrix')
                
                obj.numNodes = varargin{1}.numNodes; %n;
                obj.numMaxEdges = varargin{1}.numMaxEdges; %maxE;
                obj.nodeDistance = varargin{1}.nodeDistance; %nodeD;
                obj.nodeNumbers = 1:obj.numNodes;
                obj.initX = varargin{1}.getInitX;
                obj.nodePositions = varargin{1}.getNodePositions;
                
                obj.adjacencyMatrix = varargin{2};
                obj.graphConnection = obj.adjacencyMatrix;
                
            end
            
            
            
            if nargin == 3
                % create a new graph
                tStart = tic;
                
                obj.numNodes = varargin{1}; %n;
                obj.numMaxEdges = varargin{2}; %maxE;
                obj.nodeDistance = varargin{3}; %nodeD;
                obj.nodeNumbers = 1:obj.numNodes;
                obj.initX = 10*rand(obj.numNodes,1);
                obj.nodePositions = 10*rand(2,obj.numNodes);
                
                obj = obj.createGraph();
                obj.graphConnection = obj.adjacencyMatrix;
                
                tElapsed = toc(tStart);
                s = sprintf('Graph in: %.1fs', tElapsed);
                disp(s)
            end
        end
    end
    
    methods (Access = private)
        %Create Graph Connections
        function obj = createGraph(obj)
            %init Laplacian Matrix
            graphLaplacian = zeros(obj.numNodes, obj.numNodes);
            
            
            %init Identity Matrix
            I = eye(obj.numNodes);
            
            %flag to iterate while and iteration counter
            to_iterate = 1;
            iterationCounter = 0;
            
            while to_iterate == 1
                
                %iterate over all lines and collumns of the matrix L (all nodes)
                for k = 1:obj.numNodes
                    for l = 1:obj.numNodes
                        
                        %number of edges of each node
                        nodeKEdges = graphLaplacian(k,k);
                        nodeLEdges = graphLaplacian(l,l);
                        
                        %If node distance < r AND
                        %If node K edges less or equal to 5
                        %If node L edges less or equal to 5
                        %Ensure K-L == L-K node connection (upper part of L matrix)
                        %If node K and L are not connected
                        if (norm(obj.nodePositions(:,k)-obj.nodePositions(:,l)) < obj.nodeDistance) && ...
                                (nodeKEdges <= obj.numMaxEdges) && ...
                                (nodeLEdges <= obj.numMaxEdges) && ...
                                (abs(k-l) > 0) && ...
                                (graphLaplacian(k,l) == 0)
                            
                            %create Edge K-L
                            e =  I(:,k)-I(:,l);
                            %create matrix with K-L and L-K and add to Edge Matrix L
                            graphLaplacian = graphLaplacian + e*e';
                            pause(0.01);
                        end;
                    end;
                end;
                
                s = svd(graphLaplacian);
                
                if s(obj.numNodes-1) > 1e-3
                    to_iterate = 0;
                    break;
                end;
                
                %increase distance that nodes can connect
                obj.nodeDistance = obj.nodeDistance+0.1;
                
                iterationCounter = iterationCounter +1;
                if(iterationCounter > 1000)
                    break;
                end;
            end;
            
            % Now we have a Laplacian of an undirected graph. We want to
            % make it a directed Graph.
            initAdjacencyMat = graphLaplacian ~= 0;
            finalAdjacencyMat = eye(obj.numNodes, obj.numNodes);
            
            probBothDir = 0.1;
            while(1)
                finalAdjacencyMat = eye(obj.numNodes, obj.numNodes);
                for i = 1:obj.numNodes
                    for j = i+1:obj.numNodes
                        if initAdjacencyMat(i,j) == 1
                            bothDirections = rand;
                            
                            if bothDirections <= probBothDir
                                %do nothing, keep (i,j) = (j,i) = 1
                                finalAdjacencyMat(i,j) = 1;
                                finalAdjacencyMat(j,i) = 1;
                            else
                                %check which direction to keep
                                keepDir_ij = rand;
                                if keepDir_ij >= .5
                                    finalAdjacencyMat(j,i) = 1;
                                    finalAdjacencyMat(i,j) = 0;
                                else
                                    finalAdjacencyMat(i,j) = 1;
                                    finalAdjacencyMat(j,i) = 0;
                                end
                            end
                        end
                    end
                end
                %check if strongly connected
                if isStrongly(finalAdjacencyMat)
                    break;
                end
                probBothDir = probBothDir+0.1;
            end
            obj.adjacencyMatrix = finalAdjacencyMat;
            
            
        end
        
    end
    
    methods (Access = public)
        %Get Node Positions
        function nodepos = getNodePositions(obj)
            nodepos = [obj.nodeNumbers; obj.nodePositions];
        end
        
        %Get Graph Laplacian Matrix
        function adjMat = getAdjacencyMatrix(obj)
            adjMat = obj.adjacencyMatrix;
        end
        
        %Get Graph Connections Matrix
        function graphCM = getGraphConnections(obj)
            graphCM = obj.graphConnection;
        end
        
        %Get initial node values
        function initx = getInitX(obj)
            initx = obj.initX;
        end
        
        %Plot Graph
        function plotGraph(obj)
            wb = waitbar(0, 'Please wait... Drawing');
            h = figure;
            set(h,'Visible','off')
            hold on;
            for l = 1:obj.numNodes
                for k = l:obj.numNodes
                    if obj.graphConnection(l, k) == 1 && obj.graphConnection(k, l) == 1
                        a = obj.nodePositions(:,k);
                        b = obj.nodePositions(:,l);
                        plot([a(1) b(1)],...
                            [a(2) b(2)],...
                            'LineStyle', '-', 'Color', 'black',...
                            'LineWidth', 1,...
                            'Marker','o','MarkerFaceColor','red',...
                            'MarkerEdgeColor', 'red',...
                            'MarkerSize', 8 ...
                            );
                    else if obj.graphConnection(l, k) == 1 && obj.graphConnection(k, l) == 0 ||...
                                obj.graphConnection(l, k) == 0 && obj.graphConnection(k, l) == 1
                            
                            a = obj.nodePositions(:,k);
                            b = obj.nodePositions(:,l);
                            
                            if b(1) > a(1)
                                %From left to right
                                plot([a(1) b(1)],...
                                    [a(2) b(2)],...
                                    'LineStyle', '--', 'Color', 'blue',...
                                    'LineWidth', 1,...
                                    'Marker','o','MarkerFaceColor','red',...
                                    'MarkerEdgeColor', 'red',...
                                    'MarkerSize', 8 ...
                                    );
                            else
                                %from right to left
                                plot([a(1) b(1)],...
                                    [a(2) b(2)],...
                                    'LineStyle', '-.', 'Color', 'green',...
                                    'LineWidth', 1,...
                                    'Marker','o','MarkerFaceColor','red',...
                                    'MarkerEdgeColor', 'red',...
                                    'MarkerSize', 8 ...
                                    );
                            end
                        end
                    end
                end
                waitbar(l/obj.numNodes)
            end
            title('-- LR -. RL',...
                'FontUnits', 'points',...
                'FontWeight', 'normal',...
                'FontSize', 12,...
                'FontName', 'Times');
            hold off;
            delete(wb)
            set(h,'Visible','on')
            pause(0.01);
        end
        
        %Redefined plot Graph for new connect Matrix
        function plotCustomGraph(obj, connectMatrix)
            
            %assert that impossivle channels are kept at zero
            zeroentry = obj.graphConnection == 0;
            if ~isequal(obj.graphConnection(zeroentry), connectMatrix(zeroentry))
                disp('Connect Matrix uses unpredicted connections')
            end
            
            %plot graph nodes and custom edges based on logical matrix
            wb = waitbar(0, 'Please wait... Drawing');
            h = figure;
            set(h,'Visible','off')
            hold on;
            for l = 1:obj.numNodes
                for k = l:obj.numNodes
                    if connectMatrix(l, k) == 1
                        a = obj.nodePositions(:,k);
                        b = obj.nodePositions(:,l);
                        plot([a(1) b(1)],...
                            [a(2) b(2)],...
                            'LineStyle', '-', 'Color', 'blue',...
                            'LineWidth', 1,...
                            'Marker','o','MarkerFaceColor','red',...
                            'MarkerEdgeColor', 'red',...
                            'MarkerSize', 8 ...
                            );
                    end
                end
                waitbar(l/obj.numNodes)
            end
            hold off;
            delete(wb)
            set(h,'Visible','on')
            pause(0.01);
        end
        
        %Plot convergence of Nodes
        function plotConvergence(obj, W)
            
            finalX = zeros(obj.numNodes, 3*obj.numNodes);
            finalX(:,1) = obj.initX;
            
            for i = 2:3*obj.numNodes
                finalX(:,i)=W*finalX(:,i-1);
            end
            
            prompt = 'Plot how many nodes?\n';
            plots = input(prompt)
            if plots > 0 && plots <= obj.numNodes
                
                figure;
                hold on;
                for i = 1:plots
                    plot(finalX(i,:))
                end
                hold off;
            end
            
        end
        
        %Mean value of Initial Nodes
        function avgNodes = initAverage(obj)
            avgNodes = mean(obj.initX);
        end
        
        function pUsedEdges = percenteUsedEdges(obj, connectMatrix)
            pUsedEdges = (nnz(triu(connectMatrix-diag(diag(connectMatrix)))))...
                /((obj.numNodes*(obj.numNodes-1))/2)...
                *100;
        end
    end
end

