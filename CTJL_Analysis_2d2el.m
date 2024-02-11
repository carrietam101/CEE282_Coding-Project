classdef CTJL_Analysis_2d2el < RC_Analysis_2d1el
% Replace XYZ by your initials and rename the file accordingly before proceeding
% This is the child class of RC_Analysis_2d1el.m.
% All the protected & public data properties of parent class will be inherited to this class 


% Analysis class for 2nd order analysis of a 2-dimensional structure
    
    % Protected properties go here
    properties (Access = protected)

        % requested maximum number of load steps or increments
        numsteps

        % total step number
        step_num

        % requested load step or increment size
        ratio_req

        % requested maximum applied load ratio
        stop_ratio

        % Force resultant
        R

        % Error
        E

        % Update Nodes
        nodes2nd

        % Update Elements
        elements2nd

        restart
        defl
        react
        ele_for
        apratios

        parameters


    end
    
    % Public methods go here
    methods (Access = public)
        %% Constructor
        %    Arguments are all matrices received from Mastan2. Refer to comments in ud_2d1el.m for details.
        function self = CTJL_Analysis_2d2el(nnodes, coord, fixity, concen, nele, ends, A, Ayy, Izz, E, v, truss)
            % Inherit from the parent class
            self = self@RC_Analysis_2d1el(nnodes, coord, fixity, concen, nele, ends, A, Ayy, Izz, E, v, truss)
            self.parameters = [A, Ayy, Izz, E, v];
        end

        %% Run Analysis
        %  Run 2nd order analysis
        function RunAnalysis(self, numsteps, ratio_req, stop_ratio)
            self.InitializeOutputVariables();
            
          
            self.numsteps = numsteps;
            self.ratio_req = ratio_req;
            self.stop_ratio = stop_ratio;
            
            % Calculate how many steps needed
            if stop_ratio / ratio_req > numsteps
                self.step_num = numsteps;
            else
                self.step_num = stop_ratio / ratio_req;
            end

            % Create load vectors
            self.CreateLoadVectors;

            % Conduct 2nd order analysis
            self.SecondOrderAnalysis;
            
            % TODO: Run the analysis only if the structure is stable, i.e. the Kff matrix is well conditioned
            if self.AFLAG
                self.CreateLoadVectors();
                self.ComputeDisplacementsReactions();
                self.RecoverElementForces(); % Is covered in SecondOrderAnalysis
            end
        end
    end
    
    % Protected methods go here
    methods (Access = protected)
        %% Initialize Output Variables
        %  Initialize the matrices to be returned to Mastan2 with zeros
        function InitializeOutputVariables(self)
            self.DEFL = zeros(self.num_dof_node, self.nnodes);
            self.REACT = zeros(self.num_dof_node, self.nnodes);
            self.ELE_FOR = zeros(self.nele, self.num_dof_node*2);
        end

        %% Main 2nd Order Analysis Calculation
        function SecondOrderAnalysis(self)
            for i = 1:self.step_num
                % Construct Kt_(i-1) and check Kff
                self.CreateStiffnessMatrix();

                % Load increment
                self.Pf(1) = self.ratio_req * i;
                self.delf = self.Kff \ self.Pf;

                % Update geometry
                % TODO: this is just for HW3P2 case, I heard from Adam that
                % we need to have UpdateNodes() at node class.
                self.coord_t(:,2) = self.nodes(2).GetNodeCoord + self.delf(1:2);
                self.CreateNodes2nd;
                p = self.parameters;
                self.CreateElements2nd(p(1),(2),p(3),p(4),p(5));

                % Recover forces and compile R
                DEFL_t = self.DEFL';

                for j = 1:self.nele
                    % Obtain the displacements at the degrees of freedom corresponding to element i using linear
                    % indexing of the "DEFL_t" matrix
                    self.elements2nd(j).ComputeForces(DEFL_t(self.elements2nd(j).GetElementDOF()));
                    self.ELE_FOR(j,:) = self.elements2nd(j).GetFLocalNew();
                    self.ELE_FOR
                end
                % write R and E later
%                 self.R(:,i) = self.ELE_FOR(:,self.dof_free);
%                 self.E

            end
        end
        
        %% Create Nodes
        %  Create the nnodes x 1 vector of node objects representing all the nodes in the structure
        function CreateNodes2nd(self)
            self.nodes2nd = []; % Ensure it's empty before starting
            for i = 1:self.nnodes
                % Directly create and store the CTJL_Node_2d2el object
                self.nodes2nd = [self.nodes2nd; CTJL_Node_2d2el(i, self.coord_t(:,i))];
            end
        end

        %% Create Elements
        %  Create the nele x 1 vector of element objects representing all the elements in the structure
        %    Arguments are all matrices received from Mastan2. Refer to comments in ud_2d1el.m for details.
        function CreateElements2nd(self, A, Ayy, Izz, E, v)
            self.elements2nd = [];
            for i = 1:self.nele
                % Create an Element object and append it to the "elements" vector
                self.elements2nd = [self.elements2nd; CTJL_Element_2d2el(self.nodes2nd(self.ends(i, 1:2)), A, ...
                                    Ayy, Izz, E, v, self.truss)];
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
