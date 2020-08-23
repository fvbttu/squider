function dy = squid2dt2d(t, parms, y)
% dy = squid2dt2d(t, params, y)
%
% For general help with the SQUIDER model see help for the original
% squid2dt function (which has not been superceded). This function
% implements the bifurcated death rate version of the model (ie. the death
% rate changes at some point during the period of evaluation). t and y are
% as for squid2dt, difference is that params now have the form:
%  [beta,a,epsilon,delta,alpha,gamma_1,t_g,gamma_2,rho,q1,t1,...]
% where gamma_1 and gamma_2 are the two death rates, and t_g the time of
% switchover. Note: if gamma_1 = gamma_2 or t_g is out of time range then
% this will work exactly as original squid2dt.
%   Compartments [S,U,E,I,R,G,Q] remain unchanged, as does implementation of
% Q sequestration and output format.

% Evaluate one step of the expanded SQID+ model, which simulates an epidemic
% with a) general quarantine measures, b) unknown vs. known infecteds, and
% c) resusceptibility (loss of immunity for those once infected). This
% model divides the population into the following categories (normalized to
% sum to 1):
%     S (susceptible)
%     U (currently infected but not diagnosed; infectuous)
%     E (recovered from undiagnosed/uncounted infection)
%     I (detected infected; under quarantine)
%     R (recovered from diagnosed infection)
%     G (deaths)
%     Q (quarantined from general population)
% Parameters
%     beta : infection rate (of susceptible population)
%     a : superspreader effect (increases effective prevalence of infected population)
%     epsilon : recovery rate for undetected infecteds
%     delta : detection rate; detecteds are assumed to be no longer circulating
%        among the general population and hence no longer infective
%     alpha : recovery rate for dectected cases
%     gamma1,t_g,gamma2 : death rates for detected cases
%     rho : rate at which recovered etc become susceptible again
%     q_1 : 1st quarantine param (max quarantine rate)
%     t_1 : 2nd q param (time t at which to apply pulse)
%     q_2 : [TEMP: rate at which to release from quarantine ]
% ...

dy = zeros(length(y),1);
[S, UI, E, DI, R, X, Q] = deal(1,2,3,4,5,6,7);
b_inf = parms(1); b_2 = parms(2); e_recov = parms(3); d_tect = parms(4); a_recov = parms(5);
if t < parms(7)
    g_death = parms(6);
else
    g_death = parms(8);
end
r_susc = parms(9); q_curr = 0; r_q = 0;
if length(parms) > 9    % apply sequestration to susceptibles and undetected infecteds
    if length(parms) == 10   % constant steady rate of sequestration
        q_curr = parms(10);
    else                    % pulse sequestration (about a day depending on height of pulse)
        if length(parms) == 11
            i = 10; j = 11;
        else    % multiple pulses specified, get one closest to time t.
            [~, k] = min(abs(t - parms(11:2:end)));
            j = 9 + 2*k; i = j - 1;
        end
        if parms(i) > 0
            q_curr = pqcalc(parms(i), parms(j), t);
        elseif parms(i) < 0
            r_q = pqcalc(abs(parms(i)), parms(j), t);
        end
    end
end

% y = min(max(y, 0), 1);
newinf = b_inf * y(S) * real(y(UI)^b_2);
dy(S) = -newinf - q_curr*y(S) + r_susc*(y(E) + y(R)) + r_q*y(Q);   % change in Susceptible population
% dy(UI) = newinf - (e_recov + d_tect) * y(UI);    % new infecteds
dy(UI) = newinf - (e_recov + d_tect + q_curr) * y(UI);    % new infecteds
dy(E) = e_recov * y(UI) - r_susc * y(E);                 % unknown recoveries
dy(DI) = d_tect * y(UI) - (a_recov + g_death)*y(DI);       % dectections (confirmations)
dy(R) = a_recov * y(DI) - r_susc * y(R);                 % confirmed recovered
dy(X) = g_death * y(DI);                                 % confirmed deaths (cumulative)
% dy(Q) = q_curr * y(S) - r_susc*y(Q);                  % new quarantines
dy(Q) = q_curr * (y(S) + y(UI)) - r_q*y(Q);           % new quarantines
% dy = max(dy, -y(:));
if length(y) > Q
    dy(Q+1) = r_susc * y(R);  % log recycled detected infecteds
end

end

function qval = pqcalc(pulseval, t_pulse, t_curr)
% Calculate the current value of the q_curr or r_q. Function we use is
% designed to remove or restore 100*pulseval % of the relevant population
% by (mostly) adjusting the width of a normalized (peak val = 1) gaussian curve.

if pulseval < 0.625  % set pulse width to get full requested reduction in population
    sig1 = 0.0433*(2*pulseval-1)^3 + 0.484*pulseval^2 + 0.434;
    qval = pulseval*normpdf(t_curr,t_pulse,sig1)/normpdf(0,0,sig1);
else    % at approx .625 above equals pulseval, we need to give extra juice to get full effect
    sig1 = 2*pulseval*(1.042*pulseval - 1) + 1.068;
    if sig1 > 0.95   % kluge to handle difficulty scraping away last of population
        if sig1 <= 1.05
            sig1 = 0.95 + 1.5*(sig1 - 0.95);
        elseif sig1 < 1.1
            sig1 = 1.1 + 2*(sig1 - 1.05);
        elseif pulseval < 1
            sig1 = 1.2 + 7.5*(sig1 - 1.1);
        else
            sig1 = 1.6;
        end
    end
    qval = sig1*normpdf(t_curr,t_pulse,sig1)/normpdf(0,0,sig1);
end

end

