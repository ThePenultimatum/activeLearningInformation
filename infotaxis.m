global mapMat measurementLocation doorLocation numUniqueVisited totalSpaces visitedMat entropy epsilonEntropy T_t;
mapMat = [];
totalSpaces = 25*25;
visitedMat = ones(25,25) / totalSpaces;
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
    
    x = takeMeasurement();
    
    %%%% update posterior p_x(theta)
    
    %l_r_x = ;
    %p_r = getPriorDoor();   
    %%%%   = Likelihood_theta0(x) * p(theta0) / integral
    %%%%   (likelihood_theta(x) * p(theta) dtheta
    
    %%%% recalculate S from posterior
    
    entropy = entropyBoard();
    
    %%%% calculate control input u_i = argmin over u of E(delta_S(u))
    
    % [0, 0], [0, 1], [1, 0], [-1, 0], [0, -1] are the 5 possible controls
    % u_i
    
    u_i = [0, 0];
    
    %%%% apply input u_i
    
    measurementLocation = measurementLocation + u_i;
    T_t = [T_t, measurementLocation];
end

%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS

function s = entropyBoard()
  global totalSpaces numUniqueVisited
  priorDoor = getPriorDoor();
  s = (totalSpaces - numUniqueVisited) * priorDoor * log2(priorDoor);
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

function pn = getPriorNotDoor()
  pn = 1-(getPriorDoor());
end

function pr = getPriorDoor()
  global totalSpaces numUniqueVisited;
  pr = 1/(totalSpaces - numUniqueVisited);
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