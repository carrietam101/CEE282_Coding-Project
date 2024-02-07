classdef CTJL_Analysis_2d2el < RC_Analysis_2d1el
% Replace XYZ by your initials and rename the file accordingly before proceeding
% This is the child class of RC_Analysis_2d1el.m.
% All the protected & public data properties of parent class will be inherited to this class 


% Analysis class for 2nd order analysis of a 2-dimensional structure
    
    % Protected properties go here
    properties (Access = protected)

        % requested maximum number of load steps or increments
        numsteps

        % requested load step or increment size
        ratio_req

        % requested maximum applied load ratio
        stop_ratio

        % Force resultant
        R

        % Error
        E

        restart
        defl
        react
        ele_for
        apratios


    end
    
    % Public methods go here
    methods (Access = public)
        %% Constructor
        %    Arguments are all matrices received from Mastan2. Refer to comments in ud_2d1el.m for details.
        function self = CTJL_Analysis_2d2el(nnodes, coord, fixity, concen, nele, ends, A, Ayy, Izz, E, v, truss, ...
                                numsteps, ratio_req, stop_ratio, restart, defl, react, ele_for, apratios, ...
                                limit_state, h_stat_mes)
            
            self.num_dof_total = nnodes*self.num_dof_node;
            self.nnodes = nnodes;
            self.coord_t = coord';
            self.fixity_t = fixity';
            self.concen_t = concen';
            self.nele = nele;
            self.ends = ends;
            self.truss = truss;

            self.numsteps = numsteps;
            self.ratio_req = ratio_req;
            self.stop_ratio = stop_ratio;

            self.restart = restart;
            self.defl = defl;
            self.react = react;
            self.ele_for = ele_for;
            self.apratios = apratios;

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
%             self.R = zeros(self.dof_free,self.numsteps);
        end

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

                % Recover forces and compile R
                DEFL_t = self.DEFL';
                for j = 1:self.nele
                    % Obtain the displacements at the degrees of freedom corresponding to element i using linear
                    % indexing of the "DEFL_t" matrix
                    self.elements(j).ComputeForces(DEFL_t(self.elements(j).GetElementDOF()));
                    self.ELE_FOR(j,:) = self.elements(j).GetFLocalNew();
                end
                % write R and E later
%                 self.R(:,i) = self.ELE_FOR(:,self.dof_free);
%                 self.E



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
