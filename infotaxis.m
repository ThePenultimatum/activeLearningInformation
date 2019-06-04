clear;
global mapMat measurementLocation doorLocation numUniqueVisited totalSpaces entropy epsilonEntropy T_t prior posterior iters valsAtPos visited;
mapMat = [];
totalSpaces = 25*25;
prior = ones(25,25) / totalSpaces;
posterior = ones(25,25) / totalSpaces;
measurementLocation = [round(rand() * 25) round(rand() * 25)];
doorLocation = [round(rand() * 25) round(rand() * 25)];
numUniqueVisited = 0;
epsilonEntropy = 0.1;
%entropy = 0;

visited = zeros(25, 25);

T_t = measurementLocation;

%%%%
entropies = [];
entropy = entropyBoard();
entropies = [entropies, entropy];
valsAtPos = [];
iters = 0;











while entropy > 9 %< 5 %entropy > epsilonEntropy
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% take a measurement x
    
    x_i = takeMeasurement(); % 1 or 0 represents the value at the measurementLocation
    visited(measurementLocation(1), measurementLocation(2)) = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% update posterior p_x(theta)
    
    p_r_i = getPriorDoor(measurementLocation(1), measurementLocation(2)); % this is the prob. of the door being at measurementLocation according to the prior
    
    % this is the likelihood that the measurement x_i would turn out as the value we measured
    if (x_i == 1)%0)
        likelihood_i = getLikelihoodNotDoor(measurementLocation(1), measurementLocation(2));
    else
        likelihood_i = getLikelihoodDoor(measurementLocation(1), measurementLocation(2));
    end
    
    posterior_i = (prior .* likelihood_i); %(likelihood_i * p_r_i) / p_x_i;                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    posterior_i = posterior_i ./ (sum(sum(posterior_i)));
    posterior = posterior_i;
    prior = posterior_i;
    
    entropy_i = entropyBoard();
    
    entropy = entropy_i;
    entropies = [entropies, entropy];
    
    % [0, 0], [0, 1], [1, 0], [-1, 0], [0, -1] are the 5 possible controls u_i
    
    u_i = getBestControls();%[0, 0];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% apply input u_i
    
    valsAtPos = [valsAtPos, x_i];
    
    prior = posterior;
    
    measurementLocation = measurementLocation + u_i;
    T_t = [T_t; measurementLocation];
    
    
    
    iters = iters + 1;
end

prior
T_t
entropies
plot(T_t(:,1), T_t(:,2));
title("Trajectory for Infotaxis");
xlim([-1 26]); 
ylim([-1 26]);
xlabel("x");
ylabel("y");

%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS

function n = sumUniqueVisited()
    global prior visited;
    sumsofar = 0;
    for i=1:25
        for j=1:25
            if (visited(i,j) == 1)
                addon = prior(i,j);
                sumsofar = sumsofar + addon;
            end
        end
    end
    n = sumsofar;
end

function n = uniqueVisited()
    global visited;
    n = sum(sum(visited));
end

function us = getBestControls()
    global entropy measurementLocation prior totalSpaces numUniqueVisited;
    controlsPotential = [1, 0; 0, 1; 0, -1; -1, 0; 0, 0];
    tmp = [];
    for row=1:5
        controls_itry = controlsPotential(row,:);
        newPostry = measurementLocation + controls_itry;
        if ((newPostry(1) > 0) && (newPostry(1) < 26) && (newPostry(2) > 0) && (newPostry(2) < 26))
            tmp = [tmp; controls_itry(1) controls_itry(2)];
        end
    end
    controlsPotential = tmp;
    bestInd = 1;
    bestReduction = 0;
    sameReduction = [];
    sameInds = [];
    indices = [];
    entropyReductions = [];
    for i=1:length(controlsPotential(:,1))%5
        controls_i = controlsPotential(i,:);
        newPos = measurementLocation + controls_i;
        if ((newPos(1) > 0) && (newPos(1) < 26) && (newPos(2) > 0) && (newPos(2) < 26))
            p_t_rj = prior(newPos(1), newPos(2)); %1/(totalSpaces - numUniqueVisited);
            eds = entropy * p_t_rj + (1-p_t_rj) * 1;
            indices = [indices, i];
            entropyReductions = [entropyReductions, eds];
            if eds < bestReduction
                bestInd = i;
                bestReduction = eds;
                sameReduction = [];
            else
                if eds == bestReduction
                    sameReduction(i) = i;
                    sameInds = [sameInds, i-1, i];
                end
            end
        end
    end
    if (length(sameReduction) > 0)
        sameReduction = unique(sameReduction);
        lengthmat = length(sameInds);
        randval = sameInds(randperm(lengthmat));
        bestInd = randval(1);
    end
    bestInd;
    bestReduction
    entropyReductions
    us = controlsPotential(bestInd, :);
end

function s = entropyBoard()
  global posterior;
  total = 0;
  for i=1:25
      for j=1:25
          ijval = posterior(i, j);
          total = total + ijval * log2(ijval);
      end
  end
  s = -1*total;
end

function x = takeMeasurement()
  global measurementLocation doorLocation;
  if ((measurementLocation(1) == doorLocation(1)) && (measurementLocation(2) == doorLocation(2)))
      x = 1;
  else
      r = rand();
      l = getLikelihoodDoor(measurementLocation(1), measurementLocation(2));
      if (r < l)
          x = 1;
      else
          x = 0;
      end
  end
end

function pn = getPriorNotDoor(i, j) %%%%%%%%%%% MAKE SURE THAT THIS IS ONLY USED WHEN CHECKING FOR NON_VISITED SPACE
  pn = 1-(getPriorDoor(i, j));
end

function pr = getPriorDoor(i, j) %%%%%%%%%%% MAKE SURE THAT THIS IS ONLY USED WHEN CHECKING FOR NON_VISITED SPACE
  global prior;
  pr = prior(i, j);%1/(totalSpaces - numUniqueVisited);
end

% function likelihood = getLikelihoodDoor()
%   global doorLocation measurementLocation;
%   %global mapMat measurementLocation doorLocation;
%   rowdiff = abs(doorLocation(1) - measurementLocation(1));
%   coldiff = abs(doorLocation(2) - measurementLocation(2));
%   if ((rowdiff >= coldiff) && (rowdiff < 4))
%       likelihood = 1/(rowdiff + 1);
%   else
%       likelihood = 0.01;
%   end
% end

% function likelihood = getLikelihoodNotDoor()
%   global doorLocation measurementLocation;
%   %global mapMat measurementLocation doorLocation;
%   rowdiff = abs(doorLocation(1) - measurementLocation(1));
%   coldiff = abs(doorLocation(2) - measurementLocation(2));
%   if ((rowdiff >= coldiff) && (rowdiff < 4))
%       likelihood = 1-(1/(rowdiff + 1));
%   else
%       likelihood = 0.99;
%   end
% end

function l = getLikelihoodDoor(r, c)
    tmp = ones(25, 25);
    tmp(r,c) = 1;
    for i=0:3
        if (((r-i) > 0) && ((r+i)<26))
            for j=0:3
                if (((j-i) > 0) && ((j+i)<26))
                    if (i == 3)
                        tmp(r-i, j) = 1/4;
                        tmp(r+i, j) = 1/4;
                    end
                    if ((i == 2) && ((j > 1) && (j < 7)))
                        tmp(r-i, j) = 1/3;
                        tmp(r+i, j) = 1/3;
                    end
                    if ((i == 1) && ((j > 2) && (j < 6)))
                        tmp(r-i, j) = 1/2;
                        tmp(r+i, j) = 1/2;
                    end
                    if ((i == 0) && (j == 0))
                        tmp(r,c) = 1;
                    end
                end
            end
        end
    end
    l = tmp;
end

function l = getLikelihoodNotDoor(r, c)
    tmp = ones(25, 25);
    tmp(r,c) = 1;
    for i=0:3
        if (((r-i) > 0) && ((r+i)<26))
            for j=0:3
                if (((j-i) > 0) && ((j+i)<26))
                    if (i == 3)
                        tmp(r-i, j) = 3/4;
                        tmp(r+i, j) = 3/4;
                    end
                    if ((i == 2) && ((j > 1) && (j < 7)))
                        tmp(r-i, j) = 2/3;
                        tmp(r+i, j) = 2/3;
                    end
                    if ((i == 1) && ((j > 2) && (j < 6)))
                        tmp(r-i, j) = 1/2;
                        tmp(r+i, j) = 1/2;
                    end
                    if ((i == 0) && (j == 0))
                        tmp(r,c) = 0;
                    end
                end
            end
        end
    end
    l = tmp;
end