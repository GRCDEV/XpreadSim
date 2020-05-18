function [ DN, NDM, DMT, WTT ] = Epidemic_Simulation_ONETrace(ONETrace,c_radio, Q, NF, X0, T_setup, T_items, BO)
%  Simulation of Epidemic diffusion using the ONE contact trace considering also
%  the tranmission time.
%  Parameters:
%   Trace: trace with the position of the nodes in ONE format [Time, node, posx, posy]
%   c_radio: contact radio (meters)
%   Q: Number of items to transmit.
%   NF: Number of fixed nodes 
%   X0: Initial state : Matrix of NxQ indicating the nodes that have the
%   message.
%   T_setup: connection setup time (the time required to start a
%       transmission). If T_setup < 0, then no time is considered in the interchange.
%   T_trans: a vector of size N with the time required to transmit each item
%   BO = Buffer order --> order in which the items are transmitted.
%  Return:
%   DN: Matrix of NxQ with the delivery times of each message.
%   NDM: Not delivered -> number of messages initiated but not delivered
%   DMT: Vector with the number of messages interchanged depending on time 
%   WTT: The whole transmission time - excluding the setup time (succesfully or not succesfully message transmissions) 


function [t_Trans, nm] = GetMinInterchangeItemsTime(node_i,node_j)    
% Calculates the total time to interchange the messages
% Returns 
%    t_Trans = time of transmission (or -1 if no transmission)
%    nm = number of messages interchanged.
% 
%  ***** Global variables used ***********************************
%
%      X: a matrix of NxQ size indicating if a node have an item.
%      node_trans: a vector of size N indicating if a node is transmitting
%      N, Q: number of nodes and items.
%      T_setup: connection setup time
%      T_items: a vector of size N with the time required to transmit each item
%
%  ***************************************************************
    nm = 0;
    t_Trans = T_setup;
    if node_i <= NF 
        % Only send the message from the fixed nodes.
         for q = BO(node_i,:)
            if X(node_i,q) == 1 && X(node_j,q) == 0 
                nm = nm + 1;
                t_Trans = t_Trans + T_items(q);
            end
         end
    else
        for q = BO(node_i,:)
            % If one of the nodes has the item.
            if (X(node_i,q) + X(node_j,q)) == 1
                nm = nm + 1;
                t_Trans = t_Trans + T_items(q);
            end
        end
    end
end


function [t_Trans, nm] = InterchangeItems(node_i,node_j,t_Contact,t)    
% This local function interchange items between node i and j with a
% duration of t_contact at time t.
% Interchange messages depending on time
% Returns 
%    t_Trans = time of transmission (or -1 if no transmission)
%    nm = number of messages interchanged.
% 
%  ***** Global variables used ***********************************
%
%      X: a matrix of NxQ size indicating if a node have an item.
%      node_trans: a vector of size N indicating if a node is transmitting
%      N, Q: number of nodes and items.
%      T_setup: connection setup time
%      T_items: a vector of size N with the time required to transmit each item
%
%    ---> Perfomance issues: to avoid passing it each time the function is called 
%  ***************************************************************

    nm = 0;
    if T_setup < t_Contact 
        t_Trans = T_setup;
        if node_i <= NF 
            % Only send the message from the fixed nodes
            for q = BO(node_i,:)
                if (X(node_i,q)) == 1 && X(node_j,q) == 0 
                    if t_Trans + T_items(q) <= t_Contact
                        nm = nm + 1;
                        t_Trans = t_Trans + T_items(q);
                        X(node_j,q) = 1;
                        DN(node_j,q) = t+t_Trans;
                        DMT(ceil(t+t_Trans)) = DMT(ceil(t+t_Trans)) + 1;  
                    else
                        NDM = NDM + 1;
                        t_Trans = t_Contact;
                        break;  % No more items can be interchanged!!!
                    end
                end
            end     
        else 
            for q = BO(node_i,:)
                % If one of the nodes has the item.
                if (X(node_i,q) + X(node_j,q)) == 1
                    if t_Trans + T_items(q) <= t_Contact
                        nm = nm + 1;
                        t_Trans = t_Trans + T_items(q);
                        if X(node_i,q) == 0
                            X(node_i,q) = 1;
                            DN(node_i,q) = t+t_Trans;
                        else 
                            X(node_j,q) = 1;
                            DN(node_j,q) = t+t_Trans;
                        end
                        DMT(ceil(t+t_Trans)) = DMT(ceil(t+t_Trans)) + 1;  
                    else
                        NDM = NDM + 1;
                        t_Trans = t_Contact;
                        break;  % No more items can be interchanged!!!
                    end
                end
            end     
        end
        WTT = WTT + t_Trans - T_setup;
    else
        t_Trans = t_Contact; % Not enough time for the setup, but the time is consumed...
    end  
end


N = max(ONETrace(:,2));
t_max = max(ONETrace(:,1));

NDM = 0;

% Set nodes initial information (no node has items): 0: no item; 1: item
X = X0;
WTT = 0;

ONETrace = sortrows(ONETrace,2);   % Sort by node number

%
% PART ONE: we interpolete the position in order to all the points.
%
nodeIndex = 0;
i_start = 1;
nodeID = ONETrace(1,2);
num_points = length(ONETrace); 
time_step = 1;
Sim_time = max(ONETrace(:,1));

v_t = 0:time_step:Sim_time;

for i = 1:num_points
    if nodeID ~= ONETrace(i,2) || i == num_points
        nodeIndex = nodeIndex + 1;
        V_TIME = ONETrace(i_start:i-1,1);
        V_POSITION_Y = ONETrace(i_start:i-1,3);
        V_POSITION_X = ONETrace(i_start:i-1,4);
        % Assure that first and last time has position
        if V_TIME(1) > 0
            V_TIME = [0; V_TIME; ]; 
            V_POSITION_Y = [V_POSITION_Y(1); V_POSITION_Y];
            V_POSITION_X = [V_POSITION_X(1); V_POSITION_X];
        end

        if V_TIME(end) < Sim_time
            V_TIME = [V_TIME; Sim_time]; 
            V_POSITION_Y = [V_POSITION_Y; V_POSITION_Y(end)];
            V_POSITION_X = [V_POSITION_X; V_POSITION_X(end)];
        end


        %Simple interpolation (linear) to get the position, anytime.
        %Remember that "interp1" is the matlab function to use in order to
        %get nodes' position at any continuous time.
        vs_node(nodeIndex).v_x = interp1(V_TIME,V_POSITION_X,v_t);
        vs_node(nodeIndex).v_y = interp1(V_TIME,V_POSITION_Y,v_t);   

        nodeID = ONETrace(i,2);
        i_start = i;
    end

end

%
% PART TWO: detect de contacts
%
% Return value
DN = zeros(nodeIndex,Q);
DMT = zeros(1,t_max+1);
PosXYNodos = zeros(nodeIndex,2)*inf;

TransInfo = zeros(nodeIndex,3);   % Contains the information about the started transmission
                          % [ start_time, pair node, end_of_transmision ]
TransInfo(:,1) = ones(nodeIndex,1).*-1;

iContact = 1;
contacts = [];
fprintf('       ');
Tam_v_t = size(v_t,2);
for timeIndex = 1:Tam_v_t
    if mod(timeIndex,100) == 0
        fprintf('\b\b\b\b\b\b\b%6.2f%%',100*timeIndex/Tam_v_t);
    end
    t = v_t(timeIndex);
    
    %  First update positions of nodes for time t
    for i = 1:nodeIndex
        PosXYNodos(i,1) = vs_node(i).v_x(timeIndex);
        PosXYNodos(i,2) = vs_node(i).v_y(timeIndex);  
    end
    
    % Second, test nodes that end the transmission or the contact
    for n = 1:nodeIndex
        Node1 =n;
        if TransInfo(Node1,1) ~= -1                
            X1 = PosXYNodos(Node1,1);
            Y1 = PosXYNodos(Node1,2);
            Node2 = TransInfo(Node1,2);
            X2 = PosXYNodos(Node2,1);
            Y2 = PosXYNodos(Node2,2); 
            dist = sqrt((X2-X1)^2 + (Y2-Y1)^2);
            if t >= TransInfo(Node1,3) || dist > c_radio 
                % Transmission ends if the transmission time finish or the contact finish
                [t_Trans, nm] = InterchangeItems(Node1,Node2,t-TransInfo(Node1,1),TransInfo(Node1,1));
                % fprintf('%5.2f ITEMS Contact (%d, %5.2f, %5.2f) - (%d, %5.2f, %5.2f)  tTrans= %5.2f\n', t, Node1, X1, Y1, Node2, X2, Y2,t-TransInfo(Node1,1));               
                TransInfo(Node1,1) = -1;
                TransInfo(Node2,1) = -1;
            end
        end
    end
    
    % Third: detect contacts
    for n = 1:nodeIndex
        Node1 = n;
        X1 = PosXYNodos(Node1,1);
        Y1 = PosXYNodos(Node1,2);
        if TransInfo(Node1,1) == -1 
            % Find a node to transmit a message
            n_Min = max(NF+1,n+1);
            for n2 = n_Min:nodeIndex  
                if TransInfo(n2,1) == -1
                    Node2 = n2;
                    X2 = PosXYNodos(Node2,1);
                    Y2 = PosXYNodos(Node2,2);                
                    dist = sqrt((X2-X1)^2 + (Y2-Y1)^2);
                    if dist <= c_radio 
                        t_Min = GetMinInterchangeItemsTime(Node1,Node2);
                        if t_Min > 0
                            TransInfo(Node1,:) = [t, Node2, t+t_Min];
                            TransInfo(Node2,:) = [t, Node1, t+t_Min];
                            break; % node found -> end of for loop
                        end
                    end
                end
            end
        end
    end
end

fprintf('\b\b\b\b\b\b\b');

end
