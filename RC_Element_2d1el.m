classdef RC_Element_2d1el < handle

% Element class for 1st order analysis of a 2-dimensional structure
    
    properties (Access = protected)
        % 2x1 vector containing node objects representing the two nodes at either end of the element
        element_nodes
        
        % 6x1 vector containing the numbers assigned to the six degrees of freedom corresponding to the
        % element (first three corresponding to the first node and the next three corresponding to the second
        % node)
        element_dof
        
        % Flag that indicates whether the element comes from a truss or not
        truss
        
        % Length of the element
        L
        
        % 6x6 coordinate transformation matrix (sparse)
        gamma
        
        % 6x6 elastic stiffness matrix of the element in local coordinates (sparse)
        ke_local
        
        % 6x6 stiffness matrix of the element in global coordinates (sparse)
        k_global
        
        % 6x1 element force vector in local coordinates
        f_local
        
        % 6x1 element displacement vector in global coordinates
        del_global
        
        % 6x1 element displacement vector in local coordinates
        del_local
    end
    
    methods (Access = public)
        %% Constructor
        %    Arguments
        %      element_nodes: 2x1 vector containing node objects representing the two nodes at either end of
        %                     the element
        %      A:             Area of the element's cross section
        %      Ayy:           Shear area of the element's cross section along its local y-y axis
        %      Izz:           Moment of inertia of the element's cross section about its local z-z axis
        %      E:             Young's modulus of the element
        %      v:             Poisson's ratio of the element
        %      truss:         Flag that indicates whether the element comes from a truss or not
        function self = RC_Element_2d1el(element_nodes, A, Ayy, Izz, E, v, truss)
            self.element_nodes = element_nodes;
            self.truss = truss;
            
            self.RetrieveDOF();
            self.ComputeLength();
            self.ComputeLocalElasticStiffnessMatrix(A, Ayy, Izz, E, v);
            self.ComputeTransformationMatrix();
        end
        
        %% Get Element DOF
        %  Return "element_dof"
        function element_dof = GetElementDOF(self)
            element_dof = self.element_dof;
        end
        
        %% Compute Global Stiffness Matrix
        %  Compute the stiffness matrix of the element in global coordinates and store it in sparse format
        function ComputeGlobalStiffnessMatrix(self)
            self.k_global = self.gamma' * self.ke_local * self.gamma;
        end
        
        %% Get K Global
        %  Return "k_global"
        function k_global = GetKGlobal(self)
            k_global = self.k_global;
        end
        
        %% Compute Forces
        %  Compute the element force vector in local coordinates
        %    Arguments
        %      del_global: 6x1 element displacement vector in global coordinates
        function ComputeForces(self, del_global)
            self.del_global = del_global;
            
            % Compute the element displacement vector in local coordinates
            self.del_local = self.gamma * self.del_global;
            
            % Compute the element force vector in local coordinates
            self.f_local = self.ke_local * self.del_local;
        end
        
        %% Get F Local
        %  Return "f_local"
        function f_local = GetFLocal(self)
            f_local = self.f_local;
        end
    end
    
    methods (Access = protected)
        %% Retrieve DOF
        %  Compile the "element_dof" vector containing the numbers assigned to the six degrees of freedom
        %  corresponding to the element from the numbers assigned to the degrees of freedom of its two nodes
        function RetrieveDOF(self)
            for i = 1:2
                self.element_dof = [self.element_dof; self.element_nodes(i).GetNodeDOF()];
            end
        end
        
        %% Compute Length
        %  Compute the length of the element from the cartesian coordinates of its nodes
        function ComputeLength(self)
            axis = self.element_nodes(2).GetNodeCoord() - self.element_nodes(1).GetNodeCoord();
            self.L = norm(axis);
        end
        
        %% Compute Transformation Matrix
        %  Compute the coordinate transformation matrix of the element and store it in sparse format
        function ComputeTransformationMatrix(self)
            axis = self.element_nodes(2).GetNodeCoord() - self.element_nodes(1).GetNodeCoord();
            
            theta = cart2pol(axis(1), axis(2));
            c = cos(theta);
            s = sin(theta);
            
            gamma3 = [ c,  s,  0;
                      -s,  c,  0;
                       0,  0,  1];
            zeros3 = zeros(3);
            
            self.gamma = sparse([gamma3, zeros3;
                                 zeros3, gamma3]);
        end
        
        %% Compute Local Elastic Stiffness Matrix
        %  Check whether the element is part of a truss or not, and compute its local elastic stiffness matrix
        %  accordingly. Store the computed matrix in sparse format.
        function ComputeLocalElasticStiffnessMatrix(self, A, Ayy, Izz, E, v)
            if self.truss
                self.ke_local = E*A/self.L*sparse([1; 1; 4; 4], [1; 4; 1; 4], [1; -1; -1; 1], 6, 6);
            else
                
                % Factor that the Young's modulus needs to be multiplied by to give the Shear modulus
                Gf = 1/(2*(1 + v));

                d = [self.L/A, 0                             , 0             ;
                     0       , self.L^3/3/Izz + self.L/Gf/Ayy, self.L^2/2/Izz;
                     0       , self.L^2/2/Izz                , self.L/Izz] /E;

                T = [-1,  0,  0     , 1, 0, 0;
                      0, -1, -self.L, 0, 1, 0;
                      0,  0, -1     , 0, 0, 1];

                self.ke_local = sparse(T' / d * T);
            end
        end
    end
end