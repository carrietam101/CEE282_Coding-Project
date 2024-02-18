classdef KSAJ_Element_2d2el < RC_Element_2d1el
%% Class Definition
% Class KSAJ_Element_2d2el inherits from RC_Element_2d1el
% Element class for 2nd order analysis of a 2-dimensional structure
    
%% Protected Properties
    properties (Access = protected)
        kg_local % geometric stiffness matrix in local coordinates
    end
    
%% Public Methods
    methods (Access = public)
        %% Constructor
        % Defines the constructor of the 2nd Order element class
        % and calls the constructor of the parent element class
        function self = KSAJ_Element_2d2el(element_nodes, A, Ayy, Izz, E, v, truss)
            self = self@RC_Element_2d1el(element_nodes, A, Ayy, Izz, E, v, truss);
            self.f_local = zeros(6,1);
            self.ComputeLocalGeometricStiffnessMatrix();
        end
        
        %% Compute Global Stiffness Matrix
        %  Compute the stiffness matrix of the element in global 
        %  coordinates and store it in sparse format
        %  overwrites the parent class ComputeGlobalStiffnessMatrix

        function ComputeGlobalStiffnessMatrix(self)
            % Compute matrix in global coordinates  
            self.k_global = self.gamma.' * (self.ke_local + self.kg_local) * self.gamma;
        end

        %% Update Transformation Matrix
        function UpdateTransformationMatrix(self)
            self.ComputeTransformationMatrix();
        end
        
        %% Compute Forces for Force Recovery
        % compute the local element incremental vector for displacement
        function ComputeForces(self, del_global)
            self.del_global = del_global;
            % del_global 6x1 element displacement vector in global
            % coordinates from RC file
            
            % finds the element displacement vector
            self.del_local = self.gamma * self.del_global;

            % Using the method of natural deformation
            theta_R = atan((self.del_local(5) - self.del_local(2)) / (self.L + self.del_local(4) - self.del_local(1)));
            theta_an = self.del_local(3) - theta_R;
            theta_bn = self.del_local(6) - theta_R;
            u_n = (self.del_local(4) - self.del_local(1)) + ((self.del_local(4) - self.del_local(1))^2 + (self.del_local(5) - self.del_local(2))^2)/2 / self.L;
            dDelta_n = [0; 0; theta_an; u_n; 0; theta_bn];

            self.f_local = self.f_local + (self.ke_local + self.kg_local) * dDelta_n;
        end

        %% Update Local Geometric Stiffness Matrix
        % updates the local geometric stiffness matrix using the 
        % ComputeLocalGeometricStiffnessMatrix from the RC file
        function UpdateGeometricStiffnessMatrix(self)
            self.ComputeLocalGeometricStiffnessMatrix();
        end

        %% Get Global Element Forces
        % gets and returns transformed element forces for use in other
        % classes
        function f_global = GetGlobalElementForces(self)
            f_global = self.gamma.' * self.f_local;
        end
    end
    
%% Protected Methods
    methods (Access = protected)

        %% Compute Local Geometric Stiffness Matrix
        % computes the local geometric stiffness matrix of the element in
        % local coordinates, checks whether the element is part of a truss
        % or not and computes its local geometric stiffness matrix
        
        function ComputeLocalGeometricStiffnessMatrix(self)
            P = self.f_local(4); % defines load P using f_local from Rc files
            % f_local is a 6x1 element force vector in local coordinates
            % from the RC files
            self.kg_local = P/self.L * [1, 0, 0, -1, 0, 0;
                                        0, 6/5, self.L/10, 0, -6/5, self.L/10;
                                        0, self.L/10, (2*self.L^2)/15, 0, -self.L/10, -(self.L^2)/30;
                                        -1, 0, 0, 1, 0, 0;
                                        0, -6/5, -self.L/10, 0, 6/5, -self.L/10;
                                        0, self.L/10, -(self.L^2)/30, 0, -self.L/10, (2*self.L^2)/15];
        end
        
    end
end
