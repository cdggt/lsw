% Download Viswanath's data from his website.
% Hopefully it's still online!
fn = './viswanath.zip';
if ~isfile(fn)
    websave(fn,'https://dept.math.lsa.umich.edu/~divakar/lorenz/data.zip')
    unzip(fn)
end

% Load orbit metadata from DATA2
dataPath = './data/DATA2/WinnerInfo.dat';
fid = fopen(dataPath);
tline = fgetl(fid);
i = 0;
orbitPeriods = nan(1,5e4);
orbitSymbolLength = nan(1,5e4);
orbitSymbolNames = cell(1,5e4);
orbitLyapExponent = nan(1,5e4);
while ischar(tline)
    ix = mod(i,5);
    if ix == 0
        orbitSymbolNames{1 + floor(i/5)} = tline;
        orbitSymbolLength(1 + floor(i/5)) = length(tline);
    end
    if ix == 1
        orbitPeriods(1 + floor(i/5)) = sscanf(tline,'%e');
    end
    if ix == 3
        orbitLyapExponent(1 + floor(i/5)) = sscanf(tline,'%e');
    end
    tline = fgetl(fid);
    i = i+1;
end
fclose(fid);

% Approximation for Lorenz:
% The orbit Monodromy matrix is approximately a 1x1 matrix containing the
% lyapunov multiplier
multiplier = exp(orbitLyapExponent .* orbitPeriods);
orbitMonodromyMatrices = num2cell(multiplier);

% We only need orbit metadata (not trajectories) to compute POT weights!
% [w,getAvg] = periodicOrbitTheory(orbitSymbolLength,orbitPeriods,orbitMonodromyMatrices);

% Load orbit trajectories from DATA1
dataPath = './data/DATA1/WINNERS.txt';
fid = fopen(dataPath);
tline = fgetl(fid);
orbitTraj = {};
ixInDATA2WithOrbitTraj = [];
while ischar(tline)
    if(startsWith(tline,'orbit'))
        x = sscanf(tline,'orbit %d: %s');
        orbitNumber = x(1);
        orbitName = num2str(x(2:end),'%s')';
        
        for k = 1:length(orbitName)
            rotName = circshift(orbitName,k-1);
            ixInDATA2 = find(strcmp(rotName,orbitSymbolNames),1);
            if ~isempty(ixInDATA2)
                break
            end
        end
        if ~isempty(ixInDATA2)
            ixInDATA2WithOrbitTraj(end+1) = ixInDATA2;
            orbitTraj{end+1} = readmatrix(sprintf('./data/DATA1/orbit%d.dat',orbitNumber));
        end
    end
    tline = fgetl(fid);
end
fclose(fid);

% Compute POT weights for the orbits for which we have trajectory data
[w,getAvg] = periodicOrbitTheory(orbitSymbolLength(ixInDATA2WithOrbitTraj),orbitPeriods(ixInDATA2WithOrbitTraj),orbitMonodromyMatrices(ixInDATA2WithOrbitTraj));

% Define an observable (e.g. state space speed)
aObs = @(x) vecnorm(v(x));


% Compute orbit averages of the observable
orbitAverages = zeros(numel(ixInDATA2WithOrbitTraj),1);
aMin = Inf;
aMax = -Inf;
for k = 1:numel(orbitTraj)
    aAvg = 0;
    for t = 1:size(orbitTraj{k},2)
        x = orbitTraj{k}(:,t);
        a = aObs(x);
        aMin = min(aMin,a);
        aMax = max(aMax,a);
        aAvg = aAvg + (a - aAvg)/t;
    end
    orbitAverages(k) = aAvg;
end

%% Plotting
potAverage = w' * orbitAverages;
trueAverage = 95.69514; % Computed by J. Pughe-Sanford
histogram(orbitAverages)
xline(potAverage,'r','LineWidth',1)
xline(trueAverage,'g','LineWidth',1)
xline(aMin,'k--')
xline(aMax,'k--')
legend({'Orbit Averages','POT average','True average','Minimum/Maximum observed'},'location','bestoutside')
xlabel("Observable value")

figure
cutLog = -4;
transform = @(x) sign(x - trueAverage) .* abs(log10(abs(x - trueAverage)) - cutLog);
histogram(transform(orbitAverages))
xline(transform(potAverage),'r','LineWidth',1)
xline(transform(aMin),'k--')
xline(transform(aMax),'k--')
xt = xticks;
xl = {};
for i = 1:numel(xt)
    if xt(i) > 0
        xl{i} = sprintf('10^{%g}',cutLog + abs(xt(i)));
    elseif xt(i) == 0
        xl{i} = "0";
    else
        xl{i} = sprintf('-10^{%g}',cutLog + abs(xt(i)));
    end
end
xticklabels(xl);
legend({'Orbit Averages','POT average','Minimum/Maximum observed'},'location','bestoutside')
xlabel("Difference from True Average")

% Helper function to compute Lorenz state space speed
function s = v(s)
a = 10;
b = 28;
c = 8/3;

s = [
    a.*(s(2,:)-s(1,:));
    s(1,:).*(b-s(3,:))-s(2,:); 
    s(1,:).*s(2,:)-c.*s(3,:);
];
    
end

