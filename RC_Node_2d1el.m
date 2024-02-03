classdef RC_Node_2d1el < handle

% Node class for 1st order analysis of a 2-dimensional structure
    
    properties (Access = protected)
        % 2x1 vector containing the x and y coordinates of the node
        node_coord
        
        % 3x1 vector containing the numbers assigned to the three degrees of freedom corresponding to the node
        node_dof
    end
    
    methods (Access = public)
        %% Constructor
        %   Arguments
        %     node_number: Node number in the structure
        %     node_coord:  2x1 vector containing the x and y coordinates of the node
        function self = RC_Node_2d1el(node_number, node_coord)
            self.node_coord = node_coord;
            self.AssignDOF(node_number);
        end
        
        %% Get Node Coordinates
        %  Return "node_coord"
        function node_coord = GetNodeCoord(self)
            node_coord = self.node_coord;
        end
        
        %% Get Node DOF
        %  Return "node_dof"
        function node_dof = GetNodeDOF(self)
            node_dof = self.node_dof;
        end
    end
    
    methods (Access = protected)
        %% Assign DOF
        %  Assign numbers to the three degrees of freedom corresponding to the node
        %    Arguments
        %      node_number: Node number in the structure
        function AssignDOF(self, node_number)
            start_dof = 3*node_number - 2;
            self.node_dof = start_dof + [0; 1; 2];
        end
    end
end