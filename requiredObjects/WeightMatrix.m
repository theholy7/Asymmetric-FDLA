classdef WeightMatrix
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private, GetAccess = public)
        weightMatrix
        connectionMatrix
        spectralNorm
        l1norm
        maxWeight
        maxAbsWeight
        minWeight
        minAbsWeight
        cardinality
        edges
        sizeN
    end
    
    properties (Hidden = true)
        eigenOne
    end
    
    methods
        %Constructor
        function obj = WeightMatrix(weightMatrix)
            if nargin > 0
                weightMatrix = full(weightMatrix);
                %Assert if matrix is square
                assert(size(weightMatrix, 1) == size(weightMatrix, 2), 'Matrix is not square')
                
                %get size
                obj.sizeN = size(weightMatrix, 1);
                
                %create vector of ones with lenght n
                obj.eigenOne = ones(obj.sizeN, 1);
                
                obj.weightMatrix = weightMatrix;
                
                %calculate spectral norm
                obj.spectralNorm = norm(obj.weightMatrix - (1/obj.sizeN...
                    * obj.eigenOne*transpose(obj.eigenOne)));
                
                %calculate cardinality
                obj.cardinality = nnz(weightMatrix);
                
                %edges
                halfMatrix = triu(obj.weightMatrix,1) ~=0;
                obj.edges = nnz(halfMatrix);
                
                %connectionMatrix
                obj.connectionMatrix = weightMatrix ~=0;
                
                %l1 norm
                obj.l1norm = norm(obj.weightMatrix(:),1);
                
                %max entry
                obj.maxWeight = max(obj.weightMatrix(:));
                
                %min entry
                notZero = obj.weightMatrix ~= 0;
                tempWeightMatrix = obj.weightMatrix(notZero);
                obj.minWeight = min(tempWeightMatrix);
                
                %max abs entry
                obj.maxAbsWeight = max(abs(obj.weightMatrix(:)));
                
                %min abs entry
                obj.minAbsWeight = min(abs(tempWeightMatrix));
                
            end
        end
        
        %Conservation of Mass
        function errCons = getConsMassError(obj)
            leftEigVector = (obj.eigenOne)'*obj.weightMatrix;
            
            errCons = norm(leftEigVector-obj.eigenOne')/norm(obj.eigenOne');
        end
        
        %Fixed Point Iter
        function errCons = getFixPointErr(obj)
            RightEigVector = obj.weightMatrix*obj.eigenOne;
            
            errCons = norm(RightEigVector-obj.eigenOne)/norm(obj.eigenOne);
        end
        
        function saveWm(obj, filename)
            results = [obj.cardinality obj.spectralNorm]
            csvwrite(filename, results)
        end
        
        function specRad = spectralRadius(obj)
            e = eig(obj.weightMatrix - (1/obj.sizeN...
                    * obj.eigenOne*transpose(obj.eigenOne)));
            abs_e = abs(e);
            specRad = max(abs_e);
        end
        
        
    end
    
end

