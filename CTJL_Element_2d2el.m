classdef CTJL_Element_2d2el < RC_Element_2d1el
% Replace XYZ by your initials and rename the file accordingly before proceeding
% This is the child class of RC_Element_2d1el.m.
% All the protected & public data properties of parent class will be inherited to this class 


% Element class for 2nd order analysis of a 2-dimensional structure
    
    % Protected properties go here
    properties (Access = protected)

        % 6x6 geometric stiffness matrix of the element in local coordinates (sparse)
        kg_local
        
        % 6x6 stiffness matrix (ke_local + kg_local) of the element in global coordinates (sparse)
        k_global
        
        % dDeltap = gamma * dDelta_global
        dDeltap

        % 6x6 coordinate transformation matrix (sparse)
        gamma

        % 6x1 element (i-1) force vector in local coordinates
        f_local

        % 6x1 element (i) force vector in local coordinates
        f_local_new

        % 6x1 natural deformation 
        natDef
        
    end
    
    % Public methods go here
    methods (Access = public)
        
        % Define the constructor here and use it to call the parent class.
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
        function self = CTJL_Element_2d2el(element_nodes, A, Ayy, Izz, E, v, truss)
            self.element_nodes = element_nodes;
            
            self.RetrieveDOF();
            self.ComputeLength();
            self.ComputeLocalElasticStiffnessMatrix(A, Ayy, Izz, E, v);
            self.ComputeLocalGeometricStiffnessMatrix(A, Ayy, Izz, E, v);
            self.ComputeTransformationMatrix();
        end
        
        %% Compute Global Stiffness Matrix
        %  Compute the stiffness matrix (ke_local + kg_local) of the element in global coordinates and store it in sparse format
        function ComputeGlobalStiffnessMatrix(self)
            self.k_global = self.gamma' * (self.ke_local + self.kg_local) * self.gamma;
        end

        
        %% Get K Global
        %  Return "k_global"
        function k_global = GetKGlobal(self)
            k_global = self.k_global;
        end

        %% Recover Forces (Force Recovery)
        %  Recover the element force vector in local coordinates
        %    Arguments
        %      del_global: 6x1 element displacement vector in global coordinates
        function ComputeForces(self, del_global)
            self.del_global = del_global;
            
            % Compute the element displacement vector in local coordinates
            self.del_local = self.gamma * self.del_global;

            % Compute natural deformation
            d = del_local;
            theta_R = atan((d(5) - d(2))/(self.L + d(4) - d(1)));
            theta_an = d(3) - theta_R;
            theta_bn = d(6) - theta_R;
            u_n = (d(4) - d(1)) + ((d(4) - d(1))^2 + (d(5) - d(2))^2) / (2* self.L);
            self.natDef = [0;0;theta_an;u_n;0;theta_bn];
            
            % Update gamma
            self.ComputeTransformationMatrix;
            
            % Compute the element force vector in local coordinates
            self.f_local_new = self.gamma' * (self.f_local + self.k_global * self.natDef);

        end

        %% Get F Local New
        %  Return "f_local_new"
        function f_local_new = GetFLocalNew(self)
            f_local_new = self.f_local_new;
        end
        
    end
    
    % Protected methods go here
    methods (Access = protected)
        
        % Use this space to create any protected functions (if required)
        %% Update Transformation Matrix
        %  Update the coordinate transformation matrix of the element and store it in sparse format
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


        %% Compute Local Geometric Stiffness Matrix
        %  Check whether the element is part of a truss or not, and compute its local geometric stiffness matrix
        %  accordingly. Store the computed matrix in sparse format.
        function ComputeLocalGeometricStiffnessMatrix(self)
            if self.truss
                self.kg_local = self.f_local(4)/self.L*[1, 0, -1, 0; 0, 1, 0, -1; -1, 0, 1, 0; 0, -1, 0, 1];
            else
                self.kg_local = P/self.L * [1 , 0        , 0            , -1, 0         , 0;
                                            0 , 6/5      , self.L/10    , 0 , -6/5      , self.L/10;...
                                            0 , self.L/10, 2*self.L^2/15, 0 , -self.L/10, -self.L^2/30;...
                                            -1, 0        , 0            , 1 , 0         , 0;...
                                            0 , -6/5     , -self.L/10   , 0 , 6/5       , -self.L/10;...
                                            0 , self.L/10, -self.L^2/30 , 0 , -self.L/10, 2*self.L^2/15;];
            end
        end  
    end
end
