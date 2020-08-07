function dy = squid2dt(t, parms, y)
% dy = squid2dt(t, [iudrg], y)
%
% Evaluate one step of the expanded SQID+ model, which simulates an epidemic
% with a) general quarantine measures, b) unknown vs. known infecteds, and
% c) resusceptibility (loss of immunity for those once infected). This
% model divides the population into the following categories (normalized to
% sum to 1):
%     S (susceptible)
%     I (currently infected but not diagnosed; infectuous)
%     E (recovered from undiagnosed/uncounted infection)
%     D (detected infected; under quarantine)
%     R (recovered from diagnosed infection)
%     X (deaths)
%     Q (quarantined from general population)
% Our current implementation uses 10? parameters governing transfer between
% populations [b1,b2,e,d,a,g,r,q1,q2,q...]:
%     b_1 : infection rate (of susceptible population)
%     b_2 : superspreader effect (increases effective prevalence of infected population)
%     e : recovery rate for undetected infecteds
%     d : detection rate; detecteds are assumed to be no longer circulating
%        among the general population and hence no longer infective
%     a : recovery rate for dectected cases
%     g : death rate for detected cases
%     r : rate at which recovered etc become susceptible again
%     q_1 : 1st quarantine param (max quarantine rate)
%     q_2 : 2nd q param (time t at which to apply pulse)
%     q_3 : [TEMP: rate at which to release from quarantine ]

%     q_3 : 3rd q param (time t at which to release quarantine)

dy = zeros(7,1);
[S, I, E, D, R, X, Q] = deal(1,2,3,4,5,6,7);
b_inf = parms(1); b_2 = parms(2); e_recov = parms(3);
d_tect = parms(4); a_recov = parms(5); g_death = parms(6);
r_susc = parms(7); q_curr = 0; r_q = 0;
if length(parms) > 7    % apply sequestration to susceptibles and undetected infecteds
    if length(parms) == 8   % constant steady rate of sequestration
        q_curr = parms(8);
    else                    % pulse sequestration (about a day depending on height of pulse)
        if length(parms) == 9
            i = 8; j = 9;
        else    % multiple pulses specified, get one closest to time t.
            ptmp = parms(9:2:end); ptmp = ptmp(1:end-1) + diff(ptmp)/2;
            j = 7 + 2*find(t < [ptmp(:)',Inf], 1); i = j - 1;
        end
        if parms(i) > 0
            q_curr = pqcalc(parms(i), parms(j), t);
        elseif parms(i) < 0
            r_q = pqcalc(abs(parms(i)), parms(j), t);
        end
    end
end

% y = min(max(y, 0), 1);
newinf = b_inf * y(S) * real(y(I)^b_2);
dy(S) = -newinf - q_curr*y(S) + r_susc*(y(E) + y(R)) + r_q*y(Q);   % change in Susceptible population
% dy(I) = newinf - (e_recov + d_tect) * y(I);    % new infecteds
dy(I) = newinf - (e_recov + d_tect + q_curr) * y(I);    % new infecteds
dy(E) = e_recov * y(I) - r_susc * y(E);                 % unknown recoveries
dy(D) = d_tect * y(I) - (a_recov + g_death)*y(D);       % dectections (confirmations)
dy(R) = a_recov * y(D) - r_susc * y(R);                 % confirmed recovered
dy(X) = g_death * y(D);                                 % confirmed deaths (cumulative)
% dy(Q) = q_curr * y(S) - r_susc*y(Q);                  % new quarantines
dy(Q) = q_curr * (y(S) + y(I)) - r_q*y(Q);           % new quarantines
% dy = max(dy, -y(:));

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

