global mapMat measurementLocation doorLocation numUniqueVisited totalSpaces visitedMat entropy epsilonEntropy T_t prior posterior;
mapMat = [];
totalSpaces = 25*25;
visitedMat = ones(25,25) / totalSpaces;
prior = ones(25,25) / totalSpaces;
posterior = ones(25,25) / totalSpaces;
measurementLocation = [0, 0];
doorLocation = [15, 15];
numUniqueVisited = 0;
epsilonEntropy = 0.1;
entropy = 0;

T_t = measurementLocation;

%%%%
%entropy = 1000000000;
while entropy > epsilonEntropy
    
    %%%% take a measurement x
    
    x_i = takeMeasurement();
    
    %%%% update posterior p_x(theta)
    
    p_r_i = getPriorDoor();   
    if (x == 0)
        likelihood_i = getLikelihoodNotDoor();
        p_x_i = 1/totalSpaces;
    else
        likelihood_i = getLikelihoodDoor();
        p_x_i = 1-(1/totalSpaces);
    end
    
    posterior_i = (likelihood_i * p_r_i) / p_x_i;
    posterior(measurementLocation(1), measurementLocation(2)) = posterior_i;
    
    %%%% recalculate S from posterior
    
    entropy = entropyBoard();
    
    %%%% calculate control input u_i = argmin over u of E(delta_S(u))
    
    % [0, 0], [0, 1], [1, 0], [-1, 0], [0, -1] are the 5 possible controls
    % u_i
    
    u_i = [0, 0];
    
    %%%% apply input u_i
    
    prior(measurementLocation(1), measurementLocation(2)) = x;
    
    measurementLocation = measurementLocation + u_i;
    T_t = [T_t, measurementLocation];
end

%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS

function s = entropyBoard()
  global totalSpaces numUniqueVisited prior posterior
  %priorDoor = getPriorDoor();
  %s = (totalSpaces - numUniqueVisited) * priorDoor * log2(priorDoor);
  total = 0;
  for i=1:25
      for j=1:25
          ijval = posterior(i, j);
          total = total + ijval * log2(ijval);
      end
  end
  s = total;
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
  global prior totalSpaces numUniqueVisited;
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