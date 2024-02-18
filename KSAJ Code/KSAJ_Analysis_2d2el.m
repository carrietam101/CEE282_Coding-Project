classdef KSAJ_Analysis_2d2el < RC_Analysis_2d1el
%% Class Definition
% Class KSAJ_Analysis_2d2el inherits from RC_Analysis_2d1el
% Analysis class for 2nd order analysis of a 2-dimensional structure
    
%% Protected Properties
    properties (Access = protected)
        % properties from ud2d file
        APRATIOS double = []           % applied load ratio at end of load step i
        LIMIT_STATE double = 0         % flag indicating post-limit state analysis

        % properties created in the Analysis file
        
        DEFL_step double = []          % deflection at each step
        REACT_step double = []         % reactions at each step
        load_norm_errors double = []   % normalized error for applied loads
        energy_norm_errors double = [] % nomalized error for energy 
           
    end
    
%% Public Methods
    methods (Access = public)
        %% Constructor
        function self = KSAJ_Analysis_2d2el(nnodes, coord, fixity, concen, nele, ends, A, Ayy, Izz, E, v, truss)
             % calls the parent class constructor
            self = self@RC_Analysis_2d1el(nnodes, coord, fixity, concen, nele, ends, A, Ayy, Izz, E, v, truss);
        end

%% Run Analysis 
        % Define number of steps for while loop and run analysis for
        % specified number of steps
        function [DEFL, REACT, ELE_FOR, AFLAG, APRATIOS, LIMIT_STATE] = RunAnalysis(self, numsteps, ratio_req, stop_ratio)
            % Define number of steps for while loop
            maxStep = self.DefineStep(numsteps, ratio_req, stop_ratio);

            % Initialize variables before starting the loop
            InitializeOutputVariables(self, maxStep)

            CreateStiffnessMatrix(self)
            CreateLoadVectors(self, ratio_req)

            step = 1;
                % while loop to run through all the steps
                while step <= maxStep && self.LIMIT_STATE == 0 && self.AFLAG == 1
                    self.ComputeStep(ratio_req, step)
                    % ComputeStep will contain all functions needed in the
                    % while loop so it performs actions each step
                    step = step + 1;
                end
    
            RC_Plot_Errors(self.load_norm_errors, self.energy_norm_errors, self.APRATIOS) 
            % plot error vs applied load ratio
            [DEFL, REACT, ELE_FOR, AFLAG, APRATIOS, LIMIT_STATE] = GetMastan2Returns(self, step);
            % get returns to MASTAN2
        end

%% Compute Step (computes all functions for the While Loop)
        % Comprised of all the functions that will be used in the iteration
        % loop in RunAnalysis

        function ComputeStep(self, ratio_req, step)
                % calls each of the functions below for each step
                self.APRATIOS = [self.APRATIOS; (step * ratio_req)]; 
                % make a column and continuously add to vector
                ComputeDisplacementsReactions(self, step)
                RecoverElementForces(self, step) % from element class
                UpdateNodeCoord(self) %from node class
                UpdateStiffnessMatrix(self)
                CreateStiffnessMatrix(self)
                CheckLimitState(self)  
                ComputeError(self, step)
        end

%% Get Mastan2 Returns
        % Return DEFL, REACT, ELE_FOR, AFLAG, APRATIOS, and LIMIT_STATE to
        % MASTAN2
        function [DEFL, REACT, ELE_FOR, AFLAG, APRATIOS, LIMIT_STATE] = GetMastan2Returns(self, step)
            DEFL = self.DEFL(:, :, 1:step - 1);
            REACT = self.REACT(:, :, 1:step - 1);
            ELE_FOR = self.ELE_FOR(:, :, 1:step - 1);
            AFLAG = self.AFLAG;
            APRATIOS = self.APRATIOS;
            LIMIT_STATE = self.LIMIT_STATE;
        end
    end

%% Protected Methods
    methods (Access = protected)

%% Define Steps
        % Define the number of steps needed for the iteration loop in
        % RunAnalysis. Use the minimum of numsteps and
        % stop_ratio/ratio_req
        function maxStep = DefineStep(self, numsteps, ratio_req, stop_ratio)
            maxStep = int64(min(numsteps, (stop_ratio/ratio_req)));
        end

%% Initialize Output Variables (updated from RC file)
        % Initializes matrices to be returned to MASTAN2 of size
        % maxStepsize
        % if this max is not reached, the rest of it is deleted

        function InitializeOutputVariables(self, maxStep)
            self.DEFL = zeros(self.nnodes, self.num_dof_node, maxStep); % 3D
            self.DEFL_step = zeros(self.num_dof_node, self.nnodes); % 2D
            self.REACT = zeros(self.nnodes, self.num_dof_node, maxStep); % 3D
            self.REACT_step = zeros(self.num_dof_node, self.nnodes); % 2D
            self.ELE_FOR = zeros(self.nele, self.num_dof_node*2, maxStep); % 3D
        end
       
%% Create Nodes
        % Create the node vector with all nodes in the structure
        function CreateNodes(self)
            % makes empty node vector
            self.nodes = KSAJ_Node_2d2el.empty;
            for i = 1:self.nnodes
                % Makes a node and adds it to the node vector
                self.nodes(i) = KSAJ_Node_2d2el(i, self.coord_t(:, i));
            end
        end
        
%% Create Elements
        % Creates an element vector with all the elements in the structure

        function CreateElements(self, A, Ayy, Izz, E, v)
            % makes empty element vector
            self.elements = KSAJ_Element_2d2el.empty;
                 for i = 1:self.nele  
                % Create an Element object and adds it to the element
                % vector
                self.elements(i) = KSAJ_Element_2d2el(self.nodes(self.ends(i, 1:2)), A(i), ...
                                    Ayy(i), Izz(i), E(i), v(i), self.truss);
            end
        end

%% Create Load Vectors (updates from the RC file to become incremental)
         %  Update the function from RC and create the applied load vectors
        function CreateLoadVectors(self, ratio_req)
            % Compute vector of concentrated loads applied at the free and support degrees of freedom using
            % linear indexing of the "concen_t" matrix
            self.Pf = self.concen_t(self.dof_free) * ratio_req; 
            self.Psupp = self.concen_t(self.dof_supp) * ratio_req;
            
            % Compute vector of specified displacements using linear indexing of the "fixity_t" matrix
            self.deln = self.fixity_t(self.dof_disp);
        end

%% Compute Displacements Reactions
        function ComputeDisplacementsReactions(self, step)
            
            % Compute the displacements
            self.delf = self.Kff \ (self.Pf - self.Kfn*self.deln); % delf = deltaDelta
            
            % Compute the reactions, accounting for loads applied directly on the supports
            self.Ps = self.Ksf*self.delf + self.Ksn*self.deln - self.Psupp;
            self.Pn = self.Knf*self.delf + self.Knn*self.deln;
            
            % Format the computed displacements using linear indexing of the "DEFL" matrix
            self.DEFL_step(self.dof_free) = self.delf;
            self.DEFL_step(self.dof_disp) = self.deln;

            if step == 1
                self.DEFL(:, :, step) = self.DEFL_step.';
            else
                self.DEFL(:, :, step) = self.DEFL(:, :, step - 1) + self.DEFL_step.';
            end
            
            % Format the computed reactions using linear indexing of the "REACT" matrix
            self.REACT_step(self.dof_supp) = self.Ps;
            self.REACT_step(self.dof_disp) = self.Pn;
            
            if step == 1
                self.REACT(:, :, step) = self.REACT_step.';
            else
                self.REACT(:, :, step) = self.REACT(:, :, step - 1) + self.REACT_step.';
            end
        end

%% Recover Element Forces
        %  Recover the local element forces and format them to return to Mastan2
        function RecoverElementForces(self, step)
            for i = 1:self.nele 
                % Obtain the displacements at the degrees of freedom corresponding to element i using linear
                % indexing of the "DEFL_t" matrix
                self.elements(i).ComputeForces(self.DEFL_step(self.elements(i).GetElementDOF()));
                self.ELE_FOR(i, :, step) = self.elements(i).GetFLocal();
            end
        end

%% Update Node Coordinates
        function UpdateNodeCoord(self)
            DEFL_step_t = self.DEFL_step.';
            for i = 1:self.nnodes
                self.nodes(i).UpdateNodes(DEFL_step_t(i, [1;2])); % 1, 2 = x, y
            end
        end

%% Update Stiffness Matrix
        % Update the transformation matrix and the geometric stiffness
        % matrix
        function UpdateStiffnessMatrix(self)
            for i = 1:self.nele
                self.elements(i).UpdateTransformationMatrix();
                self.elements(i).UpdateGeometricStiffnessMatrix();
            end
        end

%% Check Limit State
        % Set the limit state using cholesky
        function CheckLimitState(self)
            [~, p] = chol(self.Kff);
            if p == 0
                self.LIMIT_STATE = 0; 
            else
                self.LIMIT_STATE = 1;
            end
        end

%% Compute Error (new function)
        % Compute errors between applied loads and energy
        function ComputeError(self, step)

            % initialize R matrix
            R = zeros(self.num_dof_total, 1);

            % iterate and set R equal to internal forces
            for i = 1:self.nele
                f_global = self.elements(i).GetGlobalElementForces();
                element_dof = self.elements(i).GetElementDOF();
                R(element_dof) = R(element_dof) + f_global;
            end

            % find the error for each step
            appliedLoads = self.APRATIOS(step) * self.concen_t(:);
            E_total = appliedLoads - R;

            % separate the error at free DOF
            % track the error
            Error_step = E_total(self.dof_free);

            % calculate load normal error at each step
            appliedLoadsFreeDOF = appliedLoads(self.dof_free);
            loadNormErrorIndex = norm(Error_step) / norm(appliedLoadsFreeDOF);
            self.load_norm_errors = [self.load_norm_errors; loadNormErrorIndex];

            % find energy norm error index
            energyNormErrorIndex = (abs(Error_step).' * abs(self.delf)) / (abs(appliedLoadsFreeDOF).' * abs(self.delf));

            % find load normal error again after recalculating error index
            self.energy_norm_errors = [self.energy_norm_errors; energyNormErrorIndex];
        end
    end
end
