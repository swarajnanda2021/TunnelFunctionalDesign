function [value,isterminal,direction] = eventfcn(t,y)
% Locate the time when height passes through zero in a decreasing
% direction
% and stop integration.

value = y(1) - 1e-8; % detect y-1/2 = 0
isterminal = 1; % stop the integration
direction = -1; % negative direction


end