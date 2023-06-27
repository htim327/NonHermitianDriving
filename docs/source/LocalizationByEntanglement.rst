=================================================================================================================================
Presentation of the Code that Demonstrates Localization when the Time Evolution is Frequently Interrupted by Quantum Entanglement
=================================================================================================================================

As stated by the title, this document presents the code that shows how systems become localized when their time evolution is frequently interrupted by quantum entanglement. The explanation for how this process works is given both by the comments within the code as well as in the corresponding article itself. First, we start with the main file that runs the file that obtains the data and then plots the information from that data.

.. code-block:: matlab

   clear;
   clc;
   % Determine the number of noise realizations you want to use
   NTot = 100;
   % Determine the amount of time you want the system evolved for
   NTime = 100;
   save('NTime.mat','NTime')
   % Generate the matrices that store the information about the participation
   % ratio
   PR1 = [];
   PR2 = [];
   % Generate the matrices that store the information about the probability of
   % particles occupying each of the j-indices
   jprobs1 = zeros(NTot,NTime,4);
   jprobs2 = zeros(NTot,NTime,4);
   % Iterate over the total number of noise realizations you want to use
   for zamp = 1:NTot
       tic
       % Run TwoDimxyQ
       TwoDimxyQ
       % Load the matrices that contain the relevant information
       load('PRveca.mat')
       load('PRvecb.mat')
       load('jprobsa.mat')
       load('jprobsb.mat')
       % Store this information in the appropriate locations
       PR1 = [PR1; PRveca];
       PR2 = [PR2; PRvecb];
       for i = 1:NTime
           for j = 1:4
               jprobs1(zamp,i,j) = jprobsa(1,j,i);
               jprobs2(zamp,i,j) = jprobsb(1,j,i);
           end
       end
       % Save this information
       save('PR1.mat','PR1')
       save('PR2.mat','PR2')
       save('jprobs1.mat','jprobs1')
       save('jprobs2.mat','jprobs2')
       clearvars -except zamp PR1 PR2 jprobs1 jprobs2 NTot NTime
       % Be continuously updated as the the progress the algorithm has made
       clc
       zamp
       toc
   end
   % Calculate the averages and standard errors associated with the
   % probabilities to have the particle occupy each of the j-indices
   jprobs1ave = zeros(NTime,4);
   jprobs2ave = zeros(NTime,4);
   jprobs1std = zeros(NTime,4);
   jprobs2std = zeros(NTime,4);
   for i = 1:NTime
       for j = 1:4
           jprobs1ave(i,j) = sum(jprobs1(:,i,j))/NTot;
           jprobs2ave(i,j) = sum(jprobs2(:,i,j))/NTot;
           jprobs1std(i,j) = std(jprobs1(:,i,j))/sqrt(NTot-1);
           jprobs2std(i,j) = std(jprobs2(:,i,j))/sqrt(NTot-1);
       end
   end
   % Plot this information
   figure('units','normalized','outerposition',[0 0 1 1]);
   errorbar(1:NTime,sum(PR1,1)/NTot,std(PR1)/sqrt(NTot),'Color','b')
   hold on
   errorbar(1:NTime,sum(PR2,1)/NTot,std(PR2)/sqrt(NTot),'Color','g')
   hold off
   title(['Participation Ratio for a System that is Measured 1000 Times for Each Time Step'],'FontSize',40,'Interpreter','latex')
   for i = 1:4
       figure('units','normalized','outerposition',[0 0 1 1]);
       errorbar(1:NTime,jprobs1ave(:,i),jprobs1std(:,i),'Color','b')
       hold on
       errorbar(1:NTime,jprobs2ave(:,i),jprobs2std(:,i),'Color','g')
       hold off
       title(['Probability of Occupying J-index ' num2str(i-1)],'FontSize',40,'Interpreter','latex')
   end

Now to present the file :math:`$\mathrm{TwoDimxyQ.m}$`, which is run in the file above and is the main file that obtains all of the data for this particular analysis.

.. code-block:: matlab

   % Determine the system size
   Li = 2;
   Lj = 4;
   LSquared = 2*Li*Lj;
   % Determine how many qubits are needed to define this system
   nqubits = log2(LSquared);
   % Determine how often the state of the AFAI system is quantum entangled
   % with the external qubit. measint = 1000 would mean that the AFAI is
   % quantum entangled with the external qubit 1000 times per driving step
   % whereas if measint = 1/100, then it would be entangled after 100 driving
   % steps.
   measint = 1000;
   % Determine how the code is going to run depending on whether the AFAI
   % system is going to be entangled multiple times per driving step or if it
   % is going to be entangled after a multiple of a single driving step.
   if (measint<1)
       timeinterupt = '0';
   else
       timeinterupt = '1';
   end
   % Set the number of particles you want in your AFAI. This algorithm was
   % really only intended to use ntimes = 1. This bug is irrelevant to this
   % proof of concept.
   ntimes = 1;
   % Determine the strength of the chemical potential
   del = 0.4;
   % Set the strength of the chemical potential disorder
   Noise = 1.5;
   % Set the strength of the temporal disorder
   tchaos = 0.2;
   % Set the size of the energy needed for hopping between sites
   J = 1.25;
   % Determine how many Floquet cycles the AFAI is going to be evolved for
   load('NTime.mat')
   NVec = 1:NTime;
   N = max(NVec);
   % Use a random seed for the random number generator
   rng('shuffle');
   % Generate the matrices for the chemical potential disorder
   magi = TwoDxyNoiseHamiltonians(Li,Lj,Noise);
   % Generate the Hamiltonians and the velocity matrices
   [H1, H2, H3, H4, H5, V1, V3] = FastTwoDxyHamiltonians(Li,Lj,J,del);
   % Generate the wave functions
   W = eye(LSquared);
   wave = W(:,1:ntimes);
   rng('shuffle');
   % Generate the variables that implement the temporal disorder
   TimeDisorder1 = -tchaos + 2*tchaos*rand(1,N);
   TimeDisorder2 = -tchaos + 2*tchaos*rand(1,N);
   TimeDisorder3 = -tchaos + 2*tchaos*rand(1,N);
   TimeDisorder4 = -tchaos + 2*tchaos*rand(1,N);
   TimeDisorder5 = -tchaos + 2*tchaos*rand(1,N);
   TimeDisorder1 = [-1 TimeDisorder1];
   TimeDisorder2 = [-1 TimeDisorder2];
   TimeDisorder3 = [-1 TimeDisorder3];
   TimeDisorder4 = [-1 TimeDisorder4];
   TimeDisorder5 = [-1 TimeDisorder5];
   wave2 = wave;
   % Generate the vectors that store information about the participation ratio
   PRveca = [];
   PRvecb = [];
   % Generate the vectors that store information about the probability the the
   % particles occupying each j-index
   jprobsa = zeros(1,Lj,N);
   jprobsb = zeros(1,Lj,N);
   aph = 0;
   % Generate the matrices that store the operator that entangle the AFAI with
   % an external particle depending on whether a particle occupies a specific
   % site or not
   measmats = zeros(2^(ntimes*nqubits+1),2^(ntimes*nqubits+1),2*Li*Lj);
   % Iterate over all of the j and i indices
   for j = (Lj-1):(-1):0
       for i = 0:(Li-1)
           aph = aph + 1;
           % locmat is the matrix that is populated corresponding to the site
           % of interest
           locmat = zeros(2^(ntimes*nqubits),2^(ntimes*nqubits));
           % notlocmat has every site populate except for the site of interest
           notlocmat = eye(2^(ntimes*nqubits),2^(ntimes*nqubits));
           % Make the entry of locmat corresponding to the site of interest 1
           locmat(1+2*i+2*Li*j,1+2*i+2*Li*j) = 1;
           % Make the entry of notlocmat corresponding to the site of interest
           % 0
           notlocmat(1+2*i+2*Li*j,1+2*i+2*Li*j) = 0;
           % Have the external particle flip its spin if the site of interest
           % is populate, otherwise do not flip the spin
           measmats(:,:,aph) = measmats(:,:,aph) + kron(locmat,[0 1; 1 0]) + kron(notlocmat,[1 0; 0 1]);
           aph = aph + 1;
           % locmat is the matrix that is populated corresponding to the site
           % of interest
           locmat = zeros(2^(ntimes*nqubits),2^(ntimes*nqubits));
           % notlocmat has every site populate except for the site of interest
           notlocmat = eye(2^(ntimes*nqubits),2^(ntimes*nqubits));
           % Make the entry of locmat corresponding to the site of interest 1
           locmat(2+2*i+2*Li*j,2+2*i+2*Li*j) = 1;
           % Make the entry of notlocmat corresponding to the site of interest
           % 0
           notlocmat(2+2*i+2*Li*j,2+2*i+2*Li*j) = 0;
           % Have the external particle flip its spin if the site of interest
           % is populate, otherwise do not flip the spin
           measmats(:,:,aph) = measmats(:,:,aph) + kron(locmat,[0 1; 1 0]) + kron(notlocmat,[1 0; 0 1]);
       end
   end
   % Count how many sites there are
   num = aph;
   % Iterate over all of the Floquet cycles
   for z = 1:N
       % Calculate the wave function up until the current Floquet cycle
       wave2 = expm(-1i*(H5+diag(magi))*(1+TimeDisorder5(z))*2*pi/5)*expm(-1i*(H4+diag(magi))*(1+TimeDisorder4(z))*2*pi/5)*expm(-1i*(H3+diag(magi))*(1+TimeDisorder3(z))*2*pi/5)*expm(-1i*(H2+diag(magi))*(1+TimeDisorder2(z))*2*pi/5)*expm(-1i*(H1+diag(magi))*(1+TimeDisorder1(z))*2*pi/5)*wave2;
       % Calculate the participation ratio for this particular Floquet cycle
       PRveca = [PRveca sum(abs(wave2).^4)];
       % Iterate over all of the sites and determine the probabilities of the
       % particles occupying each of the j-indices
       for j = 0:(Lj-1)
           probnow = 0;
           for i = 0:(Li-1)
               for k = 1:2
                   probnow = probnow + abs(wave2(k+2*i+2*Li*j))^2;
               end
           end
           jprobsa(1,j+1,z) = probnow;
       end
   end
   % Calculate the initial density matrix
   if (ntimes==1)
       density = wave(:,1)*ctranspose(wave(:,1));
   else
       density = kron(wave(:,1)*ctranspose(wave(:,1)),wave(:,2)*ctranspose(wave(:,2)));
       for i = 3:ntimes
           density = kron(density,wave(:,i)*ctranspose(wave(:,i)));
       end
   end
   % If the AFAI is entangled multiple times per driving step
   if (timeinterupt=='1')
       % Iterate over all of the Floquet cycles
       for z = 1:N
           % Generate the matrix that time evolves the system for a fraction
           % of the first driving step
           unitnow = expm(-1i*(H1+diag(magi))*(1+TimeDisorder1(z))*2*pi/(5*measint));
           for t = 2:ntimes
               unitnow = kron(unitnow,expm(-1i*(H1+diag(magi))*(1+TimeDisorder1(z))*2*pi/(5*measint)));
           end
           % Iterate over all of the times that we entangle the AFAI system
           % with external particles for this driving step
           for t = 1:measint
               % Time evolve this system for a fraction of the driving
               % step
               density = unitnow*density*ctranspose(unitnow);
               % Iterate over all of the sites
               for t2 = 1:num
                   % Add an external particle to the system
                   density = kron(density,[1 0; 0 0]);
                   % Entangle the AFAI with the external particle such that if
                   % a particle is present at the site of interest flip the
                   % external qubit, otherwise leave the qubit alone
                   density = measmats(:,:,t2)*density*ctranspose(measmats(:,:,t2));
                   % Remove the external particle
                   [rdensity] = ReducedDensity(density,ntimes*nqubits+1,1:(ntimes*nqubits));
                   density = rdensity;
               end
           end
           %%%
           % Generate the matrix that time evolves the system for a fraction
           % of the second driving step
           unitnow = expm(-1i*(H2+diag(magi))*(1+TimeDisorder2(z))*2*pi/(5*measint));
           for t = 2:ntimes
               unitnow = kron(unitnow,expm(-1i*(H2+diag(magi))*(1+TimeDisorder2(z))*2*pi/(5*measint)));
           end
           % Iterate over all of the times that we entangle the AFAI system
           % with external particles for this driving step
           for t = 1:measint
               % Time evolve this system for a fraction of the driving
               % step
               density = unitnow*density*ctranspose(unitnow);
               % Iterate over all of the sites
               for t2 = 1:num
                   % Add an external particle to the system
                   density = kron(density,[1 0; 0 0]);
                   % Entangle the AFAI with the external particle such that if
                   % a particle is present at the site of interest flip the
                   % external qubit, otherwise leave the qubit alone
                   density = measmats(:,:,t2)*density*ctranspose(measmats(:,:,t2));
                   % Remove the external particle
                   [rdensity] = ReducedDensity(density,ntimes*nqubits+1,1:(ntimes*nqubits));
                   density = rdensity;
               end
           end
           %%%
           % Generate the matrix that time evolves the system for a fraction
           % of the third driving step
           unitnow = expm(-1i*(H3+diag(magi))*(1+TimeDisorder3(z))*2*pi/(5*measint));
           for t = 2:ntimes
               unitnow = kron(unitnow,expm(-1i*(H3+diag(magi))*(1+TimeDisorder3(z))*2*pi/(5*measint)));
           end
           % Iterate over all of the times that we entangle the AFAI system
           % with external particles for this driving step
           for t = 1:measint
               % Time evolve this system for a fraction of the driving
               % step
               density = unitnow*density*ctranspose(unitnow);
               % Iterate over all of the sites
               for t2 = 1:num
                   % Add an external particle to the system
                   density = kron(density,[1 0; 0 0]);
                   % Entangle the AFAI with the external particle such that if
                   % a particle is present at the site of interest flip the
                   % external qubit, otherwise leave the qubit alone
                   density = measmats(:,:,t2)*density*ctranspose(measmats(:,:,t2));
                   % Remove the external particle
                   [rdensity] = ReducedDensity(density,ntimes*nqubits+1,1:(ntimes*nqubits));
                   density = rdensity;
               end
           end
           %%%
           % Generate the matrix that time evolves the system for a fraction
           % of the fourth driving step
           unitnow = expm(-1i*(H4+diag(magi))*(1+TimeDisorder4(z))*2*pi/(5*measint));
           for t = 2:ntimes
               unitnow = kron(unitnow,expm(-1i*(H4+diag(magi))*(1+TimeDisorder4(z))*2*pi/(5*measint)));
           end
           % Iterate over all of the times that we entangle the AFAI system
           % with external particles for this driving step
           for t = 1:measint
               % Time evolve this system for a fraction of the driving
               % step
               density = unitnow*density*ctranspose(unitnow);
               % Iterate over all of the sites
               for t2 = 1:num
                   % Add an external particle to the system
                   density = kron(density,[1 0; 0 0]);
                   % Entangle the AFAI with the external particle such that if
                   % a particle is present at the site of interest flip the
                   % external qubit, otherwise leave the qubit alone
                   density = measmats(:,:,t2)*density*ctranspose(measmats(:,:,t2));
                   % Remove the external particle
                   [rdensity] = ReducedDensity(density,ntimes*nqubits+1,1:(ntimes*nqubits));
                   density = rdensity;
               end
           end
           %%%
           % Generate the matrix that time evolves the system for a fraction
           % of the fifth driving step
           unitnow = expm(-1i*(H5+diag(magi))*(1+TimeDisorder5(z))*2*pi/(5*measint));
           for t = 2:ntimes
               unitnow = kron(unitnow,expm(-1i*(H5+diag(magi))*(1+TimeDisorder5(z))*2*pi/(5*measint)));
           end
           % Iterate over all of the times that we entangle the AFAI system
           % with external particles for this driving step
           for t = 1:measint
               % Time evolve this system for a fraction of the driving
               % step
               density = unitnow*density*ctranspose(unitnow);
               % Iterate over all of the sites
               for t2 = 1:num
                   % Add an external particle to the system
                   density = kron(density,[1 0; 0 0]);
                   % Entangle the AFAI with the external particle such that if
                   % a particle is present at the site of interest flip the
                   % external qubit, otherwise leave the qubit alone
                   density = measmats(:,:,t2)*density*ctranspose(measmats(:,:,t2));
                   % Remove the external particle
                   [rdensity] = ReducedDensity(density,ntimes*nqubits+1,1:(ntimes*nqubits));
                   density = rdensity;
               end
           end
           % After each Floquet cycle calculate the participation ratio of the
           % system
           PRnow = 0;
           for t = 1:(2^(ntimes*nqubits))
               proj = zeros(2^(ntimes*nqubits),2^(ntimes*nqubits));
               proj(t,t) = 1;
               PRnow = PRnow + abs(trace(density*proj))^2;
           end
           PRvecb = [PRvecb PRnow];
           % After each Floquet cycle calculate the probability for the
           % particle to occupy each j-index
           for j = 0:(Lj-1)
               probnow = 0;
               for i = 0:(Li-1)
                   for k = 1:2
                       probnow = probnow + abs(density(k+2*i+2*Li*j,k+2*i+2*Li*j));
                   end
               end
               jprobsb(1,j+1,z) = probnow;
           end
       end
   % If the AFAI is entangled at times that are multiple of a single driving step    
   else
       % Determine the multiple that we are supposed to entangle the AFAI
       measint2 = round(1/measint);
       aph = 0;
       % Iterate over all of the Floquet cycles
       for z = 1:N
           % Iterate over all of the driving steps
           for z2 = 1:5
               aph = aph + 1;
               if (z2==1)
                   % Time evolve the system for the first driving step
                   unitnow = expm(-1i*(H1+diag(magi))*(1+TimeDisorder1(z))*2*pi/5);
                   for z3 = 2:ntimes
                       unitnow = kron(unitnow,expm(-1i*(H1+diag(magi))*(1+TimeDisorder1(z))*2*pi/5));
                   end
                   density = unitnow*density*ctranspose(unitnow);
               elseif (z2==2)
                   % Time evolve the system for the second driving step
                   unitnow = expm(-1i*(H2+diag(magi))*(1+TimeDisorder2(z))*2*pi/5);
                   for z3 = 2:ntimes
                       unitnow = kron(unitnow,expm(-1i*(H2+diag(magi))*(1+TimeDisorder2(z))*2*pi/5));
                   end
                   density = unitnow*density*ctranspose(unitnow);
               elseif (z2==3)
                   % Time evolve the system for the third driving step
                   unitnow = expm(-1i*(H3+diag(magi))*(1+TimeDisorder3(z))*2*pi/5);
                   for z3 = 2:ntimes
                       unitnow = kron(unitnow,expm(-1i*(H3+diag(magi))*(1+TimeDisorder3(z))*2*pi/5));
                   end
                   density = unitnow*density*ctranspose(unitnow);
               elseif (z2==4)
                   % Time evolve the system for the fourth driving step
                   unitnow = expm(-1i*(H4+diag(magi))*(1+TimeDisorder4(z))*2*pi/5);
                   for z3 = 2:ntimes
                       unitnow = kron(unitnow,expm(-1i*(H4+diag(magi))*(1+TimeDisorder4(z))*2*pi/5));
                   end
                   density = unitnow*density*ctranspose(unitnow);
               elseif (z2==5)
                   % Time evolve the system for the fifth driving step
                   unitnow = expm(-1i*(H5+diag(magi))*(1+TimeDisorder5(z))*2*pi/5);
                   for z3 = 2:ntimes
                       unitnow = kron(unitnow,expm(-1i*(H5+diag(magi))*(1+TimeDisorder5(z))*2*pi/5));
                   end
                   density = unitnow*density*ctranspose(unitnow);
               end
               % If the current driving step (taking into account the total
               % number of driving steps since the very beginning of the
               % evolution) is of the correct multiple, entangle the AFAI with
               % the external particles.
               if (mod(aph,measint2)==0)
                   % Iterate over all of the sites
                   for t = 1:num
                       % Add an external particle to the system
                       density = kron(density,[1 0; 0 0]);
                       % Entangle the AFAI with the external particle such that if
                       % a particle is present at the site of interest flip the
                       % external qubit, otherwise leave the qubit alone
                       density = measmats(:,:,t)*density*ctranspose(measmats(:,:,t));
                       % Remove the external particle
                       [rdensity] = ReducedDensity(density,ntimes*nqubits+1,1:(ntimes*nqubits));
                       density = rdensity;
                   end
               end
               % If we have reached the fifth driving step, calculate the
               % participation ratio as well as the probability for the
               % particle to occupy a given j-index
               if (z2==5)
                   PRnow = 0;
                   for t = 1:(2^(ntimes*nqubits))
                       proj = zeros(2^(ntimes*nqubits),2^(ntimes*nqubits));
                       proj(t,t) = 1;
                       PRnow = PRnow + abs(trace(density*proj))^2;
                   end
                   PRvecb = [PRvecb PRnow];
                   for j = 0:(Lj-1)
                       probnow = 0;
                       for i = 0:(Li-1)
                           for k = 1:2
                               probnow = probnow + abs(density(k+2*i+2*Li*j,k+2*i+2*Li*j));
                           end
                       end
                       jprobsb(1,j+1,z) = probnow;
                   end
               end
           end
       end
   end
   % Save the information with respect to the participation ratios and the
   % probabilities for the particles to occupy a certain j-index
   save('PRveca.mat','PRveca')
   save('PRvecb.mat','PRvecb')
   save('jprobsa.mat','jprobsa')
   save('jprobsb.mat','jprobsb')

Here is the helper function that generates the matrices that implement the temporal disorder when added to the Hamiltonians.

.. code-block:: matlab

   function magi = TwoDxyNoiseHamiltonians(Li,Lj,chaos)
   % This function calculates a matrix that implements a particular
   % configuration of chemical potential disorder. This matrix is given by
   % magi and is added to each of the Hamiltonians during the time evolution.
   % Li defines the number of sites in the x-direction, Lj defines the number
   % of sites in the y-direction, and chaos defines the strength of the
   % chemical potential disorder.
   %%%
   % Calculate the total number of sites in the system and store the value in
   % LSquared.
   LSquared = 2*Li*Lj;
   % Iterate over all of the sites of the system.
   for i = 1:(LSquared)
       % Generate a random number that is drawn between -W and W.
       candy = -chaos + 2*chaos*rand;
       % Use the random number to apply a random on-site potential
       ioph(i) = candy;
       clear candy
   end
   % Return the resulting matrix as output.
   magi = ioph;
   end

Then there is the helper function that generates the Hamiltonians for the five driving steps as well as the velocity matrices for the first and third driving steps.

.. code-block:: matlab

   function [Ham1, Ham2, Ham3, Ham4, Ham5, Vel1, Vel3] = FastTwoDxyHamiltonians(Li,Lj,J,del)
   % This function generates the Hamiltonians that implement the five step
   % Floquet drive as well as the velocity matrices that are used to measure
   % the topological current during the first and third driving steps. The
   % system is defined by Li sites in the x-direction and Lj sites in the
   % y-direction, the hopping strength is given by J, and the strength of the
   % on-site potential implemented during step 5 is given by del.
   %%%
   % Define the total number of sites that defines the system with LSquared
   LSquared = 2*Li*Lj;
   % Initialize all of the Hamiltonians and the velocity matrices as matrices
   % of zeros
   Muy = zeros(LSquared);
   H1 = Muy;
   H2 = Muy;
   H3 = Muy;
   H4 = Muy;
   H5 = Muy;
   V1 = Muy;
   V3 = Muy;
   % Populate all of the Hamiltonians and the velocity matrices in the
   % appropriate locations such that they perform that actions they were
   % intended to.
   for i = 2:2:LSquared
       H1(i,(i-1)) = -J;
       H1((i-1),i) = -J;
       V1((i-1),i) = -1i*J;
       V1(i,(i-1)) = 1i*J;
   end
   clear i
   for i = 0:(Li-1)
       for j = 0:(Lj-2)
           H2((2+2*i+2*Li*(j+1)),(1+2*rem((i+1),Li)+2*Li*j)) = -J;
           H2((1+2*rem((i+1),Li)+2*Li*j),(2+2*i+2*Li*(j+1))) = -J;
           H4((2+2*i+2*Li*j),(1+2*i+2*Li*(j+1))) = -J;
           H4((1+2*i+2*Li*(j+1)),(2+2*i+2*Li*j)) = -J;
       end
       clear j
       for j = 0:(Lj-1)
           H3((1+2*rem((i+1),Li)+2*Li*j),(2+2*i+2*Li*j)) = -J;
           H3((2+2*i+2*Li*j),(1+2*rem((i+1),Li)+2*Li*j)) = -J;
           V3((1+2*rem((i+1),Li)+2*Li*j),(2+2*i+2*Li*j)) = -1i*J;
           V3((2+2*i+2*Li*j),(1+2*rem((i+1),Li)+2*Li*j)) = 1i*J;
       end
   end
   for k = 1:LSquared
       H5(k,k) = ((-1)^(k-1))*del;
   end
   % Give the results as output.
   Ham1 = H1;
   Ham2 = H2;
   Ham3 = H3;
   Ham4 = H4;
   Ham5 = H5;
   Vel1 = V1;
   Vel3 = V3;
   end

Finally, we have the helper function that calculates the reduced density matrices given an input matrix as well as a list of target qubits to construct the reduced density matrix from.

.. code-block:: matlab

   function [rdensity] = ReducedDensity(densityi,size,targets)
   % This function takes the density matrix densityi composed of size qubits
   % and calculates the reduced density matrix for the qubits given by targets
   % and returns this reduced density matrix as rdensity
   %%%
   % Determine the number of qubits that compose targets
   nq = length(targets);
   % Determine the number of qubits in densityi that are not going to compose
   % the outputted reduced density matrix
   nq2 = size - nq;
   % Initialize the matrix that will store the reduced density matrix
   redden = zeros(2^nq);
   % Iterate over all possible configurations of the qubits that will not
   % compose the reduced density matrix
   for i = 0:(2^nq2-1)
       % Express the number for the current iteration as a bitstring of length
       % nq2
       const = dec2bin(i);
       const2 = nq2 - length(const);
       for j = 1:const2
           const = ['0' const];
       end
       % count is used to determine how far across the bitstring we have gone
       % when using the information in the bitstring to generate the matrix
       % opmat that will be used to create the reduced density matrix.
       count = 0;
       % If 1 is an entry of targets, then make the first matrix that composes
       % the set of Kronecker products that generates opmat be the 2 by 2
       % identity matrix
       if sum(1==targets)
           opmat = eye(2);
       % Otherwise make the first matrix that composes this set of Kronecker
       % products be the appropriate single qubit spin vector
       else
           count = count+1;
           if (const(count)=='1')
               opmat = [0; 1];
           else
               opmat = [1; 0];
           end
       end
       % Iterate through all of the rest of the qubits (both the target qubits
       % for the reduced density matrix as well as all of the other qubits)
       % and determine whether the next matrix in the set of Kronecker
       % products should be an identity matrix or the spin up or down state
       % vector. If the qubit of interest is a target qubit for the reduced
       % density matrix then use the identity matrix otherwise use the
       % appropriate state vector.
       for j = 2:size
           if sum(j==targets)
               opmat = kron(opmat,eye(2));
           else
               count = count + 1;
               if (const(count)=='1')
                   opmat = kron(opmat,[0; 1]);
               else
                   opmat = kron(opmat,[1; 0]);
               end
           end
       end
       % Use opmat to perform operations on densityi in order to obtain the
       % appropriate information about the reduced density matrix and add this
       % information to redden.
       redden = redden + ctranspose(opmat)*densityi*opmat;
   end
   % Normalize redden
   redden = redden/trace(abs(redden));
   % Return the reduced density matrix as rdensity
   rdensity = redden;
   end
