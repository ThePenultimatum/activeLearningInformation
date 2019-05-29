global mapMat measurementLocation doorLocation numUniqueVisited totalSpaces visitedMat entropy epsilonEntropy T_t prior posterior iters valsAtPos;
mapMat = [];
totalSpaces = 25*25;
visitedMat = ones(25,25) / totalSpaces;
prior = ones(25,25) / totalSpaces;
posterior = ones(25,25) / totalSpaces;
measurementLocation = [1, 1];
doorLocation = [15, 15];
numUniqueVisited = 0;
epsilonEntropy = 0.1;
%entropy = 0;

T_t = measurementLocation;

%%%%
entropy = entropyBoard();
valsAtPos = [];
iters = 0;











while iters < 3 %entropy > epsilonEntropy
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% take a measurement x
    
    x_i = takeMeasurement(); % 1 or 0 represents the value at the measurementLocation
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% update posterior p_x(theta)
    
    p_r_i = getPriorDoor(measurementLocation(1), measurementLocation(2)); % this is the prob. of the door being at measurementLocation according to the prior
    
    % this is the likelihood that the measurement x_i would turn out as the value we measured
    if (x_i == 1)%0)
        likelihood_i = getLikelihoodNotDoor();
        p_x_i = 1/totalSpaces; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ??????????????????????????
    else
        likelihood_i = getLikelihoodDoor();
        p_x_i = 1-(1/totalSpaces);
    end
    
    posterior_i = (likelihood_i * p_r_i) / p_x_i;                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% ????? I need to update the posterior entire matrix
    %%%%%%%%%%%%%%%% ????? I need to update the posterior entire matrix
    %%%%%%%%%%%%%%%% ????? I need to update the posterior entire matrix%%%%%%%%%%%%%%%% ????? I need to update the posterior entire matrix
    %%%%%%%%%%%%%%%% ????? I need to update the posterior entire matrix%%%%%%%%%%%%%%%% ????? I need to update the posterior entire matrix
    %%%%%%%%%%%%%%%% ????? I need to update the posterior entire matrix%%%%%%%%%%%%%%%% ????? I need to update the posterior entire matrix
    %%%%%%%%%%%%%%%% ????? I need to update the posterior entire matrix
    
    
    
    posterior(measurementLocation(1), measurementLocation(2)) = posterior_i;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% recalculate S from posterior
    
    entropy_i = entropyBoard();
    
    entropy = entropy_i;
    
    
    
    
    %%%%%%%%% EVERYTHING ABOVE SEEMS TO CHECK OUT........... CHECK THAT THE
    %%%%%%%%% BOARDS ARE BEING UPDATED PROPERLY AND THE WAY THAT THEY
    %%%%%%%%% SHOULD BE, PRIOR VS POSTERIOR UPDATES
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calculate control input u_i = argmin over u of E(delta_S(u))
    
    % [0, 0], [0, 1], [1, 0], [-1, 0], [0, -1] are the 5 possible controls u_i
    
    u_i = getBestControls();%[0, 0];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% apply input u_i
    
    valsAtPos = [valsAtPos, x_i];
    
    prior = posterior %updatePrior();
    
    measurementLocation = measurementLocation + u_i;
    T_t = [T_t, measurementLocation];
    
    
    
    %entropy = 0;
    iters = iters + 1;
    prior(1:5,1:5);
end

%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS

function p = updatePrior()
    global numUniqueVisited totalSpaces iters T_t valsAtPos;
    numUniqueVisited = numUniqueVisited + 1;
    tmp = ones(25, 25) / (totalSpaces - numUniqueVisited);
    for i=1:iters
        pos = T_t(i, :);
        tmp(pos(1), pos(2)) = valsAtPos(i);
    end
    p = tmp;
end

function us = getBestControls()
    global entropy measurementLocation prior totalSpaces numUniqueVisited;
    controlsPotential = [1, 0; 0, 1; 0, -1; -1, 0; 0, 0];
    bestInd = 1;
    bestReduction = 0;
    sameReduction = [];
    sameInds = [];
    indices = [];
    entropyReductions = [];
    for i=1:5
        controls_i = controlsPotential(i,:);
        newPos = measurementLocation + controls_i;
        if ((newPos(1) > 0) && (newPos(1) < 26) && (newPos(2) > 0) && (newPos(2) < 26))
            p_t_rj = prior(newPos(1), newPos(2)); %1/(totalSpaces - numUniqueVisited);
            eds = -1 * entropy * p_t_rj + (1-p_t_rj) * 0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5 this  is multiplied by -1....... should this change min vs max?
            %%%%%%%%%%%%%%%%%%%%%%%%% THIS ADDED TERM ABOVE THAT I AM
            %%%%%%%%%%%%%%%%%%%%%%%%% MAKING ZERO SHOULD INSTEAD ACCOUNT
            %%%%%%%%%%%%%%%%%%%%%%%%% FOR THE POTENTIAL NEXT STEPS, SO SEE
            %%%%%%%%%%%%%%%%%%%%%%%%% WHICH ONES BEYOND THE CURRENT ONE CAN
            %%%%%%%%%%%%%%%%%%%%%%%%% ACTUALLY HAVE A MOVE OR GET A
            %%%%%%%%%%%%%%%%%%%%%%%%% REDUCTION
            indices = [indices, i];
            entropyReductions(i) = eds;
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
    bestReduction;
    entropyReductions;
    us = controlsPotential(bestInd, :);
end

function s = entropyBoard()
  global posterior;
  %priorDoor = getPriorDoor();
  %s = (totalSpaces - numUniqueVisited) * priorDoor * log2(priorDoor);
  total = 0;
  for i=1:25
      for j=1:25
          ijval = posterior(i, j);
          total = total + ijval * log2(ijval);
      end
  end
  s = -1*total;
end

function expDSu = expDSU(u)
  u;
  expDSu = 0;
end

function x = takeMeasurement()
  global measurementLocation doorLocation;
  if ((measurementLocation(1) == doorLocation(1)) && (measurementLocation(2) == doorLocation(2)))
      x = 1;
  else
      r = rand();
      l = getLikelihoodDoor();
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

function likelihood = getLikelihoodDoor()
  global doorLocation measurementLocation;
  %global mapMat measurementLocation doorLocation;
  rowdiff = abs(doorLocation(1) - measurementLocation(1));
  coldiff = abs(doorLocation(2) - measurementLocation(2));
  if ((rowdiff >= coldiff) && (rowdiff < 4))
      likelihood = 1/(rowdiff + 1);
  else
      likelihood = 0.01;
  end
end

function likelihood = getLikelihoodNotDoor()
  global doorLocation measurementLocation;
  %global mapMat measurementLocation doorLocation;
  rowdiff = abs(doorLocation(1) - measurementLocation(1));
  coldiff = abs(doorLocation(2) - measurementLocation(2));
  if ((rowdiff >= coldiff) && (rowdiff < 4))
      likelihood = 1-(1/(rowdiff + 1));
  else
      likelihood = 0.99;
  end
end