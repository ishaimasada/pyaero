function [Vexit, Texit, Pexit, Mexit] = converging_nozzle(NPR, NPRcrit, Tt, Pt, gamma, R, Cp)
    if NPR >= NPRcrit
        fprintf("Converging Nozzle is choked\n")
        Mexit = 1;
        Texit = (2 / (gamma + 1)) * Tt;
        Pexit = Pt * (1 / NPR_crit);
        Vexit = sqrt(gamma * R * Texit);
    else
        fprintf("Converging Nozzle is unchoked\n")
        Pexit = Pt * (1 / NPR);
        Texit = Tt * (1 / NPR) ^ ((gamma - 1)/ gamma);
        Vexit = sqrt(2 * Cp * (Tt - Texit));
        Mexit = Vexit / sqrt(gamma * R * Texit);
    end
end