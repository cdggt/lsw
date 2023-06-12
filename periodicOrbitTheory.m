function [weights,getAverage] = periodicOrbitTheory(orbitSymbolLengths,orbitPeriods,orbitMonodromyMatrices)
%% [weights,getAverage] = periodicOrbitTheory(orbitSymbolLengths,orbitPeriods,orbitMonodromyMatrices)
%
%  Inputs:
%  ======
%
%  orbitSymbolLengths: P-many positive integers giving the symbolic
%                      dynamics of each orbit
%
%  orbitPeriods: P-many positive real numbers giving the temporal period of
%                each orbit
%
%  orbitMonodromyMatrices: Cell array of P-many matrices with no
%                          eigenvalues near positive one. Scalar values are
%                          accepted in place of 1x1 arrays.
%
%  Outputs:
%  ========
%
%  weights: The coefficients w_p in the expansion <a>_t = sum_p w_p <a>_p
%
%  getAverage: A function getAverage(orbitAverages) which computes <a>_t
%              using the vector of P-many orbit averages <a>_p. This
%              function uses Periodic Orbit Theoretical techniques WITHOUT
%              weights of any kind

    N = max(orbitSymbolLengths);
    P = numel(orbitSymbolLengths);
    
    C = zeros(1,N);
    dCds = zeros(1,N);
    dCdb = zeros(1,N);
    dC2dBdap = zeros(N,P);
    
    Q = zeros(1,N);
    dQds = zeros(1,N);
    dQdb = zeros(1,N);
    dQ2dBdap = zeros(N,P);
    
    if numel(orbitPeriods) ~= P
        error("Number of orbit periods provided (%g) does not match number of orbits (%g)",numel(orbitPeriods),P)
    end
    if numel(orbitMonodromyMatrices) ~= P
        error("Number of orbit monodromy matrices provided (%g) does not match number of orbits (%g)",numel(orbitMonodromyMatrices),P)
    end
    
    function [F,dFds,weights] = spectralDeterminantAtBetaEqualsZero(s)
        C = zeros(1,N);
        dCds = zeros(1,N);
        dC2dBdap = zeros(N,P);
        % These four lines implement the sum over (j = np * r)
        for p = 1:P
            np = orbitSymbolLengths(p);
            for r = 1:floor(N/np)
                j = np * r;

                Tp = orbitPeriods(p);
                numerator = exp(-r * Tp * s);
                denominator = abs(det(1 - mpower(orbitMonodromyMatrices{p},r)));
                C(j) = C(j) - (1/r) * (numerator/denominator); 
                dCds(j) = dCds(j) + (numerator/denominator) * Tp;

                if nargout > 2
                    % Second derivative of C with respect to beta and <a>_p
                    dC2dBdap(j,p) = dC2dBdap(j,p) - (numerator/denominator) * Tp;
                end
            end
        end
        
        Q = zeros(1,N);
        dQds = zeros(1,N);
        dQ2dBdap = zeros(N,P);
        % Compute [Q = exp(C)] treating C and Q as polynomials
        for j = 1:N
            i = 1:(j-1);
            jMinusI = j - i;
            prefactor = jMinusI ./ j;
            
            % "Cumulants" formula for polynomial exponential
            Q(j) = C(j) + sum(prefactor .* C(jMinusI) .* Q(i));

            % Using product rule on the above sum(...):
            dQds(j) = dCds(j) + sum(prefactor .* (dCds(jMinusI) .* Q(i) + C(jMinusI) .* dQds(i)));

            if nargout > 2
                dQ2dBdap(j,:) = dC2dBdap(j,:) + sum(prefactor' .* (dC2dBdap(jMinusI,:) .* Q(i)' + C(jMinusI)' .* dQ2dBdap(i,:)));
            end
        end
        
        F = 1 + sum(Q);
        dFds = sum(dQds);
        
        if nargout > 2
            dF2dBdap = sum(dQ2dBdap,1);
            weights = (- dF2dBdap ./ dFds)';
        end
    end

    % Newton's method to find F(s) = 0
    maxM = 1e5;
    sCur = -0.1; % Initial guess
    for m = 1:maxM
        [F,dFds] = spectralDeterminantAtBetaEqualsZero(sCur);
        update = -F/dFds;
        if abs(update) < 1e2 * eps
            break
        end
        sCur = sCur + update;
    end
    if m == maxM
        error("Failed to converge during newton iteration!")
    end
    
    s0 = sCur;
    [~,~,weights] = spectralDeterminantAtBetaEqualsZero(s0);
    
    
    function [temporalAverage,dFdb] = computeTemporalAverageWithoutWeights(orbitAverages)
        C = zeros(1,N);
        dCds = zeros(1,N);
        dCdb = zeros(1,N);
        
        % These four lines implement the sum over (j = np * r)
        for p = 1:P
            np = orbitSymbolLengths(p);
            for r = 1:floor(N/np)
                j = np * r;

                Tp = orbitPeriods(p);
                numerator = exp(-r * Tp * s0);
                denominator = abs(det(1 - mpower(orbitMonodromyMatrices{p},r)));
                C(j) = C(j) - (1/r) * (numerator/denominator); 
                dCds(j) = dCds(j) + (numerator/denominator) * Tp;
                dCdb(j) = dCdb(j) - (numerator/denominator) * Tp * orbitAverages(p);
            end
        end
        
        Q = zeros(1,N);
        dQds = zeros(1,N);
        dQdb = zeros(1,N);
        
        % Compute [Q = exp(C)] treating C and Q as polynomials
        for j = 1:N
            i = 1:(j-1);
            jMinusI = j - i;
            prefactor = jMinusI ./ j;
            
            % "Cumulants" formula for polynomial exponential
            Q(j) = C(j) + sum(prefactor .* C(jMinusI) .* Q(i));

            % Using product rule on the above sum(...):
            dQds(j) = dCds(j) + sum(prefactor .* (dCds(jMinusI) .* Q(i) + C(jMinusI) .* dQds(i)));
            dQdb(j) = dCdb(j) + sum(prefactor .* (dCdb(jMinusI) .* Q(i) + C(jMinusI) .* dQdb(i)));
        end
        
        F = 1 + sum(Q);
        dFds = sum(dQds);
        dFdb = sum(dQdb);
        
        temporalAverage = - dFdb / dFds;
        
        fprintf("\nDifference from average computed with weights: %.3e\n",temporalAverage - orbitAverages * weights)
    end
    getAverage = @computeTemporalAverageWithoutWeights;
end