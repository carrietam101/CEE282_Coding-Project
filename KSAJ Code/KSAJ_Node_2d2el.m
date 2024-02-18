classdef KSAJ_Node_2d2el < RC_Node_2d1el
%% Class Definition
% Class KSAJ_Node_2d2el inherits from RC_Node_2d1el.m. 
% Node class for 2nd order analysis of a 2-dimensional structure

%% Protected properties
    properties (Access = protected)        
    end
    
%% Public Methods
    methods (Access = public)
        %% Constructor
        function self = KSAJ_Node_2d2el(node_number, node_coord)
            % calls parent class constructor and uses properties from ud2d
            % file to create Node objects

            self = self@RC_Node_2d1el(node_number, node_coord);
        end

        %% Update Node Coordinates
        % update node coordinates for each iteration using ddelta property 
        % from Analysis class
        function UpdateNodes(self, delta)
            self.node_coord = self.node_coord + delta.';
        end 
    end
    
%% Protected Methods
    methods (Access = protected)
    end
end
