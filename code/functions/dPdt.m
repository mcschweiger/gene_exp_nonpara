function dPdt = dPdt(t,P)
% Chemical Master Equation Setup
% keyboard
global A
% Chemical Master Equations-Proabability state vector
dPdt = (sparse(A)*sparse(P));
% Below is where we check for stopping conditions
MAXTIME = 1e4;   % Max time in seconds
MINSTEP = 1e-17; % Minimum step
persistent tprev
if isempty(tprev)
    tprev = -inf;
end
% keyboard
% if isempty(elapsedtime)
%     elapsedtime = tic;
% end
timestep = t - tprev;
tprev = t;
if (timestep > 0) && (timestep < MINSTEP)
    error(['Stopped. Time step is too small: ' num2str(timestep)])
% elseif toc(elapsedtime) > MAXTIME
%     keyboard
%     error('Stopped. Taking too long.')
end
end