function [DEFL,REACT,ELE_FOR,AFLAG] = ud_2d1el(...
		nnodes,coord,concen,fixity,nele,ends,A,Izz,Zzz,Ayy,...
		E,v,Fy,YldSurf,Wt,webdir,w,thermal,truss,anatype)

% UD_2D1EL performs a user defined two-dimensional
% first-order elastic analysis of a structural system.
% It is assumed that all information is defined in the
% X-Y plane.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Dictionary of Variables
%     Input Information:
%       nnodes         ==  total number of nodes
%       coord(i,1:2)   ==  node i's coordinates
%                            coord(i,1) = X coordinate
%                            coord(i,2) = Y coordinate
%       concen(i,1:3)  ==  concentrated loads for node i's 3 d.o.f.
%                            concen(i,1) = force in global X direction
%                            concen(i,2) = force in global Y direction
%                            concen(i,3) = moment about global Z axis
%       fixity(i,1:3)  ==  prescribed displacements for node i's 3 d.o.f.
%                          Note: A free d.o.f. will have a value of NaN
%                          and hence, you will find the Matlab function
%                          isnan very useful.
%                          Examples: If fixity(15,1) is set to NaN, then node 15's
%                                      X-disp component is free;
%                                    If fixity(2,3) is set to 0.0, then node 2's
%                                      Z-rotation component is supported;
%                                    If fixity(5,2) is set to -2.1, then node 5's
%                                      Y-disp component is supported and defined
%                                      with a settlement of -2.1 units.
%                            fixity(i,1) = prescribed disp. in global X direction
%                            fixity(i,2) = prescribed disp. in global Y direction
%                            fixity(i,3) = prescribed rotation about global Z axis
%       nele           ==  total number of elements
%       ends(i,1:8)   ==  element i's nodal information
%                            ends(i,1) = start node #
%                            ends(i,2) = finish node #
%                            ends(i,3) = flag to indicate whether or not flexural
%                            moment is released at start node.  ends(i,3)=0 not
%                            released (rigid connection); ends(i,3)=1 flexural
%                            moment is released (pinned connection); ends(i,3)=3
%                            at least one of the flexural moments are partially or fully
%                            released (see below for connection stiffness attributes)
%                            ends(i,4) = flag to indicate whether or not flexural
%                            moment is released at finish node.  ends(i,4)=0 not
%                            released (rigid connection); ends(i,4)=1 flexural
%                            moment is released (pinned connection); ends(i,4)=2
%                            at least one of the flexural moments are partially or fully
%                            released (see below for connection stiffness attributes)
%                            ends(i,5) = rotational spring stiffness at the start
%                            node and about element i's local z-z axis.
%                            ends(i,6) = rotational spring stiffness at the finish
%                            node and about element i's local z-z axis.
%                            ends(i,7) = connection moment capacity Mpz at the start
%                            node and about element i's local z-z axis.
%                            ends(i,8) = connection moment capacity Mpz at the finish
%                            node and about element i's local z-z axis.
%       A(i)           ==  element i's cross sectional area
%       Izz(i)         ==  element i's moment of inertia about its local z-z axis
%       Zzz(i)         ==  element i's plastic section modulus about its local z-z axis
%       Ayy(i)         ==  element i's effective shear area along its local y-y axis
%       E(i)           ==  element i's material elastic modulus, Young's Modulus
%       v(i)           ==  element i's material Poisson's ratio
%       Fy(i)          ==  element i's material yield strength
%       YldSurf(i)     ==  element i's yield surface maximum values
%                              YldSurf(i,1) = maximum P/Py value
%                              YldSurf(i,2) = maximum Mz/Mpz value
%                              YldSurf(i,3) = maximum My/Mpy value
%       Wt(i)          ==  element i's material weight density
%       webdir(i,1:3)  ==  element i's unit web vector.  This is a unit vector
%                          that defines the element's local y-y axis with respect
%                          to the global coordinate system.  It is based only on the
%                          structures undeformed geometry.
%                            webdir(i,1) = x component of element's unit web vector
%                            webdir(i,2) = y component of element's unit web vector
%                            webdir(i,3) = z component of element's unit web vector
%                                          (which will always be zero for 2D analysis)
%                          NOTE: An element's 3x3 rotation matrix, [g], is constructed
%                          as follows: First, calculate a unit vector, x_vect, that
%                          describes the element's local x-axis (be sure to include
%                          all three components with the z component always being zero).
%                          Second, take the cross product of x_vect and webdir(i,:) to
%                          obtain z_vect, i.e. z_vect = cross(x_vect,webdir(i,:)). Third,
%                          set z_vect to a unit vector, i.e. z_vect = z_vect/norm(z_vect).
%                          For 2D analysis, z_vect will be either [0 0 1] or [0 0 -1].
%                          Finally, the first row of [g] is x_vect, its second row is
%                          webdir(i,:), and its third row is z_vect.
%       w(i,1:2)         ==  element i's uniform load which references its
%                            local coordinate system
%                              w(i,1) = x component of uniform load
%                              w(i,2) = y component of uniform load
%       thermal(i,1:4)   ==  element i's thermal strain effects which reference its
%                            local coordinate system
%                              thermal(i,1) = coefficient of thermal expansion
%                              thermal(i,2) = change in temperature at centroid
%                              thermal(i,3) = linear temperature gradient in local y-dir
%                                           = (T_up_y - T_btm_y) / depth_y
%                              thermal(i,4) = linear temperature gradient in local z-dir
%                                           = (T_up_z - T_btm_z) / width_z
%       truss            ==  flag to indicate if structure is a truss or not
%                              truss = 0   System is not a truss
%                              truss = 1   System is a truss
%       anatype          ==  flag to indicate which type of analysis is requested
%                              anatype = 1  First-Order Elastic
%                              anatype = 2  Second-Order Elastic
%                              anatype = 3  First-Order Inelastic
%                              anatype = 4  Second-Order Inelastic
%                              anatype = 5  Elastic Buckling (Eigenvalue)
%                              anatype = 6  Inelastic Buckling (Eigenvalue)
%
%     Output Information:
%       DEFL(i,1:3)      ==  node i's calculated 3 d.o.f. deflections
%                              DEFL(i,1) = displacement in X direction
%                              DEFL(i,2) = displacement in Y direction
%                              DEFL(i,3) = rotation about Z direction
%       REACT(i,1:3)     ==  reactions for supported node i's 3 d.o.f.
%                              REACT(i,1) = force in X direction
%                              REACT(i,2) = force in Y direction
%                              REACT(i,3) = moment about Z direction
%       ELE_FOR(i,1:6)  ==  element i's internal forces and moments
%                            Note: All values reference the element's local
%                                  coordinate system.
%                              ELE_FOR(i,1)  = x-force at start node
%                              ELE_FOR(i,2)  = y-force at start node
%                              ELE_FOR(i,3)  = z-moment at start node
%                              ELE_FOR(i,4)  = x-force at end node
%                              ELE_FOR(i,5)  = y-force at end node
%                              ELE_FOR(i,6)  = z-moment at end node
%       AFLAG            ==  logical flag to indicate if a successful
%                            analysis has been completed
%                              AFLAG = 1     Successful
%                              AFLAG = 0     Unstable Structure
%                              AFLAG = inf   No analysis code available
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   For this code to work, ensure that the following files are in Matlab's path:
%       RC_Analysis_2d1el.m
%       RC_Element_2d1el.m
%       RC_Node_2d1el.m

% Instantiate an object of the Analysis class
analysis = RC_Analysis_2d1el(nnodes, coord, fixity, concen, nele, ends, A, Ayy, Izz, E, v, truss);

% Run the 1st order analysis
analysis.RunAnalysis();

% Extract the matrices to be returned to Mastan2
[DEFL, REACT, ELE_FOR, AFLAG] = analysis.GetMastan2Returns();

end
