classdef CTJL_Node_2d2el < RC_Node_2d1el
% Replace XYZ by your initials and rename the file accordingly before proceeding
% This is the child class of RC_Node_2d1el.m.
% All the protected & public data properties of parent class will be inherited to this class 

% Node class for 2nd order analysis of a 2-dimensional structure
    
    % Protected properties go here
    properties (Access = protected)
    end
    
    % Public methods go here
    methods (Access = public)
        %% Constructor
        %   Arguments
        %     node_number: Node number in the structure
        %     node_coord:  2x1 vector containing the x and y coordinates of the node
        function self = CTJL_Node_2d2el(node_number, node_coord)
            self = self@RC_Node_2d1el(node_number, node_coord);
        end
        
        %% Update Nodes
        function UpdateNodes(self, Disp)
            self.node_coord = self.GetNodeCoord + Disp;
        end

%         %% Get Node Coordinates
%         %  Return "node_coord"
%         function node_coord = GetNodeCoord(self)
%             node_coord = self.node_coord;
%         end
%         
%         %% Get Node DOF
%         %  Return "node_dof"
%         function node_dof = GetNodeDOF(self)
%             node_dof = self.node_dof;
%         end

    end
    
    % Protected methods go here
    methods (Access = protected)
%         %% Assign DOF
%         %  Assign numbers to the three degrees of freedom corresponding to the node
%         %    Arguments
%         %      node_number: Node number in the structure
%         function AssignDOF(self, node_number)
%             start_dof = 3*node_number - 2;
%             self.node_dof = start_dof + [0; 1; 2];
%         end
    end
end
