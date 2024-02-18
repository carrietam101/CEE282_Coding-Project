classdef CTJL_Analysis_2d2el < RC_Analysis_2d1el
% Replace XYZ by your initials and rename the file accordingly before proceeding
% This is the child class of RC_Analysis_2d1el.m.
% All the protected & public data properties of parent class will be inherited to this class 


% Analysis class for 2nd order analysis of a 2-dimensional structure
    
    % Protected properties go here
    properties (Access = protected)

        % requested maximum number of load steps or increments
        numsteps

        % max step number
        maxsteps

        % requested load step or increment size
        ratio_req

        % requested maximum applied load ratio
        stop_ratio

        % applied load ratio at previous step i
        APRATIOS

        % flag indicating post-limit state analysis
        limit_state

        % deflection at each step i
        DEFL_i

        % reactions at each step i
        REACT_i

        % load norm errors
        load_norm_errors

        % energy norm errors
        energy_norm_errors

    end
    
    % Public methods go here
    methods (Access = public)
        %% Constructor
        %    Arguments are all matrices received from Mastan2. Refer to comments in ud_2d1el.m for details.
        function self = CTJL_Analysis_2d2el(nnodes, coord, fixity, concen, nele, ends, A, Ayy, Izz, E, v, truss)
            % Inherit from the parent class
            self = self@RC_Analysis_2d1el(nnodes, coord, fixity, concen, nele, ends, A, Ayy, Izz, E, v, truss)
        end

        %% Run Analysis
        %  Run 2nd order analysis
        function [DEFL, REACT, ELE_FOR, AFLAG, APRATIOS, LIMIT_STATE] = RunAnalysis(self, numsteps, ratio_req, stop_ratio)
            clc;

            % Initialize output variables
            self.InitializeOutputVariables();
            
            % Record variables as global variables
            self.numsteps = numsteps;
            self.ratio_req = ratio_req;
            self.stop_ratio = stop_ratio;
            
            % Calculate how many steps needed
            if stop_ratio / ratio_req > numsteps
                self.maxsteps = numsteps;
            else
                self.maxsteps = stop_ratio / ratio_req;
            end

            % Create stiffness matrix
            self.CreateStiffnessMatrix();

            % Create load vectors
            self.CreateLoadVectors();

            % Initialize the first step
            i = 1;
            self.limit_state = 0;
            self.AFLAG = 1;

            % Conduct 2nd order analysis
            while i <= self.maxsteps && self.limit_state == 0 && self.AFLAG == 1
                self.SecondOrderAnalysis(i);
                % Update steps
                i = i + 1;
            end

            % Plot errors
            RC_Plot_Errors(self.load_norm_errors, self.energy_norm_errors, self.APRATIOS); 
            
            % Return to MASTAN2
            [DEFL, REACT, ELE_FOR, AFLAG, APRATIOS, LIMIT_STATE] = GetMastan2Returns(self, i);

        end
        %% Get Mastan2 Returns
        % Return DEFL, REACT, ELE_FOR, AFLAG, APRATIOS, and LIMIT_STATE to
        % MASTAN2
        function [DEFL, REACT, ELE_FOR, AFLAG, APRATIOS, LIMIT_STATE] = GetMastan2Returns(self, i)
            DEFL = self.DEFL(:, :, 1:i - 1);
            REACT = self.REACT(:, :, 1:i - 1);
            ELE_FOR = self.ELE_FOR(:, :, 1:i - 1);
            AFLAG = self.AFLAG;
            APRATIOS = self.APRATIOS;
            LIMIT_STATE = self.limit_state;
        end
     
    end
    
    % Protected methods go here
    methods (Access = protected)
        %% Initialize Output Variables
        %  Initialize the matrices to be returned to Mastan2 with zeros
        function InitializeOutputVariables(self)
            % DEFL and REACT are nodes x 6 dofs x maxsteps
            self.DEFL = zeros(self.nnodes, self.num_dof_node, self.maxsteps);
            self.REACT = zeros(self.nnodes, self.num_dof_node, self.maxsteps);

            % DEFL_i and REACT_i are 6 dofs x nodes
            self.DEFL_i = zeros(self.num_dof_node, self.nnodes);
            self.REACT_i = zeros(self.num_dof_node, self.nnodes);

            % ELE_FOR is nele x 6 dofs x maxsteps
            self.ELE_FOR = zeros(self.nele, self.num_dof_node*2, self.maxsteps);
        end

        %% Main 2nd Order Analysis Calculation
        function SecondOrderAnalysis(self, step)
            % Update applied load ratio
            self.APRATIOS = [self.APRATIOS; (step * self.ratio_req)];

            % Construct Kt_(i-1) and check Kff
            self.CreateStiffnessMatrix();

            % Update displacement
            self.ComputeDisplacementsReactions(step)

            % Recover forces
            DEFL_t = self.DEFL_i'; % nodes x 3 dof vector
            for j = 1:self.nele
                self.elements(j).ComputeForces(self.DEFL_i(self.elements(j).GetElementDOF()));
                self.ELE_FOR(j, :, step) = self.elements(j).GetFLocal();
            end

            % Update geometry
            for j = 1:self.nnodes
                self.nodes(j).UpdateNodes(DEFL_t(j, [1;2]));
            end
            
            % Update gamma and Kg
            for k = 1:self.nele
                self.elements(k).UpdateTransformationMatrix();
                self.elements(k).UpdateGeometricStiffnessMatrix();
            end

            % Create stiffness matrix
            self.CreateStiffnessMatrix();

            % Check limit state
            self.CheckLimitState(); 

            % Compute Error
            ComputeError(self, step);
        end
        
        %% Create Nodes
        %  Create the nnodes x 1 vector of node objects representing all the nodes in the structure
        function CreateNodes(self)
            self.nodes = CTJL_Node_2d2el.empty; % Ensure it's empty before starting
            for j = 1:self.nnodes
                % Directly create and store the CTJL_Node_2d2el object
                self.nodes(j) = CTJL_Node_2d2el(j, self.coord_t(:,j));
            end
        end

        %% Create Elements
        %  Create the nele x 1 vector of element objects representing all the elements in the structure
        %    Arguments are all matrices received from Mastan2. Refer to comments in ud_2d1el.m for details.
        function CreateElements(self, A, Ayy, Izz, E, v)
            self.elements = CTJL_Element_2d2el.empty;
            for j = 1:self.nele
                % Create an Element object and append it to the "elements" vector
                self.elements = [self.elements; CTJL_Element_2d2el(self.nodes(self.ends(j, 1:2)), A(j), ...
                                    Ayy(j), Izz(j), E(j), v(j), self.truss)];
            end
        end

        %% Create Load Vectors
        %  Create the applied load vectors for each steps
        function CreateLoadVectors(self)
            % Compute vector of concentrated loads applied at the free and support degrees of freedom using
            % linear indexing of the "concen_t" matrix
            self.Pf = self.concen_t(self.dof_free) * self.ratio_req;
            self.Psupp = self.concen_t(self.dof_supp) * self.ratio_req;
            
            % Compute vector of specified displacements using linear indexing of the "fixity_t" matrix
            self.deln = self.fixity_t(self.dof_disp);
        end

        %% Compute Displacements Reactions
        %  Compute the displacements and reactions and format them to return to Mastan2
        function ComputeDisplacementsReactions(self, step)
            
            % Compute the displacements
            self.delf = self.Kff \ (self.Pf - self.Kfn*self.deln);
            
            % Compute the reactions, accounting for loads applied directly on the supports
            self.Ps = self.Ksf*self.delf + self.Ksn*self.deln - self.Psupp;
            self.Pn = self.Knf*self.delf + self.Knn*self.deln;
            
            % Format the computed displacements using linear indexing of the "DEFL" matrix
            self.DEFL_i(self.dof_free) = self.delf;
            self.DEFL_i(self.dof_disp) = self.deln;

            if step == 1
                self.DEFL(:, :, step) = self.DEFL_i.';
            else
                self.DEFL(:, :, step) = self.DEFL(:, :, step - 1) + self.DEFL_i.';
            end

            % Format the computed reactions using linear indexing of the "REACT" matrix
            self.REACT_i(self.dof_supp) = self.Ps;
            self.REACT_i(self.dof_disp) = self.Pn;
            if step == 1
                self.REACT(:, :, step) = self.REACT_i.';
            else
                self.REACT(:, :, step) = self.REACT(:, :, step - 1) + self.REACT_i.';
            end
        end

        %% Check Limit State
        % Set the limit state using cholesky
        function CheckLimitState(self)
            % Cholesky factorization
            [~, p] = chol(self.Kff); 

            if p == 0
                % limit state not yet reached, system is loading
                self.limit_state = 0; 
            else
                % limit state reached or exceeded, system is unloading
                self.limit_state = 1;
            end
        end

        %% Check Errors
        % Compute the errors of the analysis
        function ComputeError(self, step)

            % Initialize resultant matrix
            R = zeros(self.num_dof_total, 1);

            % For each element
            for i = 1:self.nele
                % Get element forces & degrees of freedom
                f_global = self.elements(i).GetFGlobal();
                element_dof = self.elements(i).GetElementDOF();

                % Index resultant matrix to equal to associated internal forces
                R(element_dof) = R(element_dof) + f_global;
            end

            % Calculate total error
            P = self.APRATIOS(step) * self.concen_t(:);
            E_total = P - R;

            % Extract errors & loads at free DOFs
            E_free = E_total(self.dof_free);
            P_free = P(self.dof_free);

            % Calculate load norm error index at each step
            loadNormErrorIndex = norm(E_free) / norm(P_free);
            self.load_norm_errors = [self.load_norm_errors; loadNormErrorIndex];

            % Calculate energy norm error index
            energyNormErrorIndex = (abs(E_free).' * abs(self.delf)) / (abs(P_free).' * abs(self.delf));

            % Find load normal error again after recalculating error index
            self.energy_norm_errors = [self.energy_norm_errors; energyNormErrorIndex];
        end
    end
end
