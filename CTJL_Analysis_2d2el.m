classdef CTJL_Analysis_2d2el < RC_Analysis_2d1el
% Replace XYZ by your initials and rename the file accordingly before proceeding
% This is the child class of RC_Analysis_2d1el.m.
% All the protected & public data properties of parent class will be inherited to this class 


% Analysis class for 2nd order analysis of a 2-dimensional structure
    
    % Protected properties go here
    properties (Access = protected)
        % Matrices received from Mastan2. Refer to comments in ud_2d1el.m for details.
        nnodes
        nele
        ends
        truss
        
        % Transposes of the corresponding matrices received from Mastan2. Refer to comments in ud_2d1el.m for 
        % details.
        coord_t
        fixity_t
        concen_t
        
        % Total number of degrees of freedom in the structure
        num_dof_total
        
        % Total number of degrees of freedom that are free, have specified displacements, and are supported
        num_dof_free
        num_dof_disp
        num_dof_supp
        
        % Vectors of the free, displaced, and support degree of freedom numbers
        dof_free
        dof_disp
        dof_supp
        
        % nnodes x 1 vector of node objects representing all the nodes in the structure
        nodes
        
        % nele x 1 vector of element objects representing all the elements in the structure
        elements
        
        % Global stiffness matrix for the structure (sparse)
        K
        
        % Sub-matrices of K (sparse)
        Kff
        Kfn
        Knf
        Knn
        Ksf
        Ksn
        
        % Sub-vectors of the force vector
        Pf
        Pn
        Ps
        
        % Vector of forces applied directly on the supports
        Psupp
        
        % Sub-vectors of the displacement vector
        delf
        deln
        
        % Matrices to be returned to Mastan2. Refer to comments in ud_2d1el.m for details.
        DEFL
        REACT
        ELE_FOR
        AFLAG

        % requested maximum number of load steps or increments
        numsteps

        % requested load step or increment size
        ratio_req

        % requested maximum applied load ratio
        stop_ratio

    end
    
    % Public methods go here
    methods (Access = public)
        %% Constructor
        %    Arguments are all matrices received from Mastan2. Refer to comments in ud_2d1el.m for details.
        function self = CTJL_Analysis_2d2el(nnodes, coord, fixity, concen, nele, ends, A, Ayy, Izz, E, v, truss)
            self.num_dof_total = nnodes*self.num_dof_node;
            
            self.nnodes = nnodes;
            self.coord_t = coord';
            self.fixity_t = fixity';
            self.concen_t = concen';
            self.nele = nele;
            self.ends = ends;
            self.truss = truss;
            
            self.CreateNodes();
            self.CreateElements(A, Ayy, Izz, E, v);
            self.ClassifyDOF();
        end

        %% Run Analysis
        %  Run 2nd order analysis
        function RunAnalysis(self)
            self.InitializeOutputVariables();
            
            % Conduct 2nd order analysis
            self.SecondOrderAnalysis;
            
            % TODO: Run the analysis only if the structure is stable, i.e. the Kff matrix is well conditioned
            if self.AFLAG
                self.CreateLoadVectors();
                self.ComputeDisplacementsReactions();
                self.RecoverElementForces();
            end
        end
    end
    
    % Protected methods go here
    methods (Access = protected)
        %% Main 2nd Order Analysis Calculation
        function SecondOrderAnalysis(self)
            for i = 1:self.numsteps
                % Construct Kt_(i-1) and check Kff
                self.CreateStiffnessMatrix();
                self.delf = self.Kff \ self.ratio_req;

                % Update geometry
                self.coord_t(:,self.dof_free) = self.delf;
                self.CreateNodes();
                self.CreateElements();


            end


        end

        %% Create Stiffness Matrix
        %  Create the global stiffness matrix for the structure and store it in sparse format
        function CreateStiffnessMatrix(self)
            
            % Initialize the vectors that will be store the coordinates and values of the non-zero elements
            % in K
            K_row = [];
            K_col = [];
            K_data = [];
            
            % Loop over all elements and append the contribution of each element's global stiffness matrix to 
            % the "K_row", "K_col", and "K_data" vectors
            for i = 1:self.nele
                self.elements(i).ComputeGlobalStiffnessMatrix();
                [row, col, data] = find(self.elements(i).GetKGlobal());
                
                element_dof = self.elements(i).GetElementDOF();
                K_row = [K_row; element_dof(row)];
                K_col = [K_col; element_dof(col)];
                K_data = [K_data; data];
            end
            
            % Convert the "K_row", "K_col", and "K_data" vectors to a sparse matrix
            self.K = sparse(K_row, K_col, K_data, self.num_dof_total, self.num_dof_total);
            
            self.ComputeStiffnessSubMatrices();
            self.CheckKffMatrix();
        end

    end
end
