=========================================================================
Presenting code for the AFAI under the influence of non-Hermitian driving
=========================================================================

This document is simply going to give you the code that used to implement the non-Hermitian driving protocol along with the standard driving steps for the Anomalous Floquet-Anderson Insulator (AFAI). The explanation for how the code works is given by the comments within the code. First, we start with the main file that runs the appropriate sub-file :math:`$\mathrm{TwoDimxyQ.m}$` and then plots the corresponding information obtained.

.. code-block:: matlab

   clear;
   clc;
   % Generate the matrix that will store the information about the charge
   % pumped per Floquet cycle that is obtained without the application of
   % the non-Hermitian drive.
   QTot1 = [];
   % Generate the matrix that will store the information about the charge
   % pumped per Floquet cycle that is obtain using the non-Hermitian driving
   % methods
   QTot2 = [];
   % Determine the number of noise realizations you want to use
   NTot = 100;
   % Iterate over all of the noise realizations
   for zamp = 1:NTot
       tic
       % Run TwoDimxyQ
       TwoDimxyQ
       % Load the information about the charge pumped per Floquet cycle that
       % is obtained without the application of the non-Hermitian drive.
       load('Q.mat')
       % Store this information in QTot1
       QTot1 = [QTot1; Q];
       % Load the information about the charge pumped per Floquet cycle that
       % is obtained using the non-Hermitian driving methods.
       load('Qb.mat')
       % Store this information in QTot2
       QTot2 = [QTot2; Qb];
       % Save these two matrices
       save('QTot1.mat','QTot1')
       save('QTot2.mat','QTot2')
       clearvars -except zamp QTot1 QTot2 NTot
       % The following lines will help you see how much progress your
       % algorithm has made
       clc
       zamp
       toc
   end
   % Plot all of this information with errorbars
   figure('units','normalized','outerposition',[0 0 1 1]);
   errorbar(1:1000,sum(QTot1,1)/NTot,std(QTot1)/sqrt(NTot-1),'Color','b')
   hold on
   errorbar(1:1000,sum(QTot2,1)/NTot,std(QTot2)/sqrt(NTot-1),'Color','g')
   hold off
   title(['NTot = ' num2str(NTot) ', Gamma = 0.01 and Gamma2 = 0.01' ...
       ''],'FontSize',40,'Interpreter','latex')

Now the file :math:`$\mathrm{TwoDimxyQ.m}$` is presented.

.. code-block:: matlab

   % Determine the system size
   Li = 2;
   Lj = 4;
   % Determine the total number of sites of the system
   LSquared = 2*Li*Lj;
   % Determine how many qubits are needed to define this system
   nqubits = log2(LSquared);
   % Set the number of particles you want in your AFAI
   ntimes = 2;
   % Determine the strength of the chemical potential
   del = 0.4;
   % Set the strength of the chemical potential disorder
   Noise = 1.5;
   % Set the strength of the temporal disorder
   tchaos = 0.2;
   % Set the size of the energy needed for hopping between sites
   J = 1.25;
   % Set the strength of the noise determined by the two noise models
   gamma = 0.01;
   gamma2 = 0.01;
   % Determine how many Floquet cycles the AFAI is going to be evolved for
   NVec = 1:1000;
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
   % Generate the variables that store the information about the charge pumped
   % for every Floquet cycle
   P1 = 0;
   P3 = 0;
   QVec = [0];
   Q = [];
   P1b = 0;
   P3b = 0;
   QVecb = [0];
   Qb = [];
   % Iterate over all of the Floquet cycles
   for z = 1:N
       % Calculate the wave functions up until the first and third driving
       % steps of the current Floquet cycle
       wave2 = expm(-1i*(H5+diag(magi))*(1+TimeDisorder5(z))*2*pi/5)*expm(-1i*(H4+diag(magi))*(1+TimeDisorder4(z))*2*pi/5)*expm(-1i*(H3+diag(magi))*(1+TimeDisorder3(z))*2*pi/5)*expm(-1i*(H2+diag(magi))*(1+TimeDisorder2(z))*2*pi/5)*expm(-1i*(H1+diag(magi))*(1+TimeDisorder1(z))*2*pi/5)*wave2;
       wave3 = expm(-1i*(H2+diag(magi))*(1+TimeDisorder2(z+1))*2*pi/5)*expm(-1i*(H1+diag(magi))*(1+TimeDisorder1(z+1))*2*pi/5)*wave2;
       % Generate the matrices that express the velocity matrices integrated
       % over the appropriate time intervals for the first and third driving
       % steps
       PMatrix3 = zeros(LSquared);
       for t = 1:2:length(V1(1,:))
           u1 = magi(t);
           u2 = magi((t+1));
           PMatrix1(t:(t+1),t:(t+1)) = TwoDPxyMatrix(V1(t:(t+1),t:(t+1)),u1,u2,(1+TimeDisorder1(z+1))*2*pi/5,0,J,1,Li,Lj,magi);
       end
       clear t
       PMatrix3 = TwoDPxyMatrix([0 -1i*J; 1i*J 0],0,0,(1+TimeDisorder3(z+1))*2*pi/5,0,J,3,Li,Lj,magi);
       % Calculate the charge pumped for the first and third Floquet cycles
       for s = 1:length(wave(1,:))
           P1 = P1 + ctranspose(wave2(:,s))*PMatrix1*wave2(:,s);
           P3 = P3 + ctranspose(wave3(:,s))*PMatrix3*wave3(:,s);
       end
       if sum(z==NVec)
           QVec = [QVec real(P1 - P3)/(2*Li)];
           Q = [Q (QVec(end)-QVec(end-1))/Step];
           save('Q.mat','Q')
       end
   end
   % Generate the density matrix that defines the AFAI with all of the
   % particles in it
   if (ntimes==1)
       density = wave(:,1)*ctranspose(wave(:,1));
   else
       density = kron(wave(:,1)*ctranspose(wave(:,1)),wave(:,2)*ctranspose(wave(:,2)));
       for i = 3:ntimes
           density = kron(density,wave(:,i)*ctranspose(wave(:,i)));
       end
   end
   aph = 0;
   % Generate the matrix measmats which will determine how many particles are
   % located at the site that we want to push particles away from vs how many
   % particles are located at the site that we want to push particles to
   measmats = zeros(2^(ntimes*nqubits+1),2^(ntimes*nqubits+1),Li*Lj,1);
   % Generate rotmats which will rotate particles between the site that we
   % want to push particles away from and the site that we want to push
   % particles to
   rotmats = zeros(2^(nqubits),2^(nqubits),Li*Lj,8);
   meascheck = [];
   rotcheck = [];
   % Iterate over all j indices corresponding to the bottom half of the
   % cylinder
   for j = (Lj-1):(-1):round(Lj/2)
       % Iterate over all of the i indices
       for i = 0:(Li-1)
           aph = aph + 1;
           % If ind1 describes the site where we want to push particles away
           % from with alpha=1 describing the A site, i describing the index in the
           % i-direction, and j describing the index in the j-direction, then
           % calculate measmat which determines how many particles are located
           % at this site vs. how many particles are located at the site above
           % that one. Also calculate the rotation matrix that transfers
           % particles between these to sites.
           ind1 = [1 i j];
           ind2 = [2 i (j-1)];
           [rotmat,measmat] = PresenceRevealed2(Li,Lj,ntimes,ind1,ind2);
           rotmats(:,:,aph,1) = rotmat;
           measmats(:,:,aph,1) = measmat;
           % Calculate the rotation matrix that transfers particles directly
           % to the right of ind1 = [1 i j];
           ind1 = [1 i j];
           ind2 = [2 i j];
           [rotmat,measmat] = PresenceRevealed2(Li,Lj,ntimes,ind1,ind2);
           rotmats(:,:,aph,2) = rotmat;
           % Calculate the rotation matrix that transfers particles directly
           % to the left of ind1 = [1 i j];
           i2 = mod(i-1,Li);
           ind1 = [1 i j];
           ind2 = [2 i2 j];
           [rotmat,measmat] = PresenceRevealed2(Li,Lj,ntimes,ind1,ind2);
           rotmats(:,:,aph,3) = rotmat;
           % Calculate the rotation matrix that transfers particles directly
           % above ind1 = [1 i j];
           i2 = mod(i-1,Li);
           j2 = j+1;
           if j2>(Lj-1)
               rotmats(:,:,aph,4) = eye(2^(nqubits));
           else
               ind1 = [1 i j];
               ind2 = [2 i2 j2];
               [rotmat,measmat] = PresenceRevealed2(Li,Lj,ntimes,ind1,ind2);
               rotmats(:,:,aph,4) = rotmat;
           end
           % Now we are going to calculate the rotation matrices corresponding
           % to the site where we want to push particles to so that we can
           % construct the appropriate Kraus operators
           j3 = j - 1;
           i3 = i;
           % Generate the appropriate rotation matrix that transfers particles
           % between this B site to the site directly below it
           ind1 = [2 i3 j3];
           ind2 = [1 i j];
           [rotmat,measmat] = PresenceRevealed2(Li,Lj,ntimes,ind1,ind2);
           rotmats(:,:,aph,5) = rotmat;
           % Transfer particles between this site and the site directly to the
           % left of it
           ind1 = [2 i3 j3];
           ind2 = [1 i3 j3];
           [rotmat,measmat] = PresenceRevealed2(Li,Lj,ntimes,ind1,ind2);
           rotmats(:,:,aph,6) = rotmat;
           % Transfer particles between this site and the site directly to the
           % right of it
           ind1 = [2 i3 j3];
           i2 = mod(i3+1,Li);
           ind2 = [1 i2 j3];
           [rotmat,measmat] = PresenceRevealed2(Li,Lj,ntimes,ind1,ind2);
           rotmats(:,:,aph,7) = rotmat;
           % Transfer particles between this site and the site directly above it
           j2 = j3-1;
           i2 = mod(i3+1,Li);
           ind1 = [2 i3 j3];
           ind2 = [1 i2 j2];
           if j2<(0)
               rotmats(:,:,aph,8) = eye(2^(nqubits));
           else
               [rotmat,measmat] = PresenceRevealed2(Li,Lj,ntimes,ind1,ind2);
               rotmats(:,:,aph,8) = rotmat;
           end
           aph = aph + 1;
           % If ind1 describes the site where we want to push particles away
           % from with alpha=2 describing the B site, i describing the index in the
           % i-direction, and j describing the index in the j-direction, then
           % calculate measmat which determines how many particles are located
           % at this site vs. how many particles are located at the site above
           % that one. Also calculate the rotation matrix that transfers
           % particles between these to sites.
           i2 = mod(i+1,Li);
           ind1 = [2 i j];
           ind2 = [1 i2 (j-1)];
           [rotmat,measmat] = PresenceRevealed2(Li,Lj,ntimes,ind1,ind2);
           rotmats(:,:,aph,1) = rotmat;
           measmats(:,:,aph,1) = measmat;
           % Calculate the rotation matrix that transfers particles directly
           % to the left of ind1 = [2 i j];
           ind1 = [2 i j];
           ind2 = [1 i j];
           [rotmat,measmat] = PresenceRevealed2(Li,Lj,ntimes,ind1,ind2);
           rotmats(:,:,aph,2) = rotmat;
           % Calculate the rotation matrix that transfers particles directly
           % to the right of ind1 = [2 i j];
           i2 = mod(i+1,Li);
           ind1 = [2 i j];
           ind2 = [1 i2 j];
           [rotmat,measmat] = PresenceRevealed2(Li,Lj,ntimes,ind1,ind2);
           rotmats(:,:,aph,3) = rotmat;
           % Calculate the rotation matrix that transfers particles directly
           % below ind1 = [2 i j];
           j2 = j + 1;
           if j2>(Lj-1)
               rotmats(:,:,aph,4) = eye(2^(nqubits));
           else
               ind1 = [2 i j];
               ind2 = [1 i j2];
               [rotmat,measmat] = PresenceRevealed2(Li,Lj,ntimes,ind1,ind2);
               rotmats(:,:,aph,4) = rotmat;
           end
           % Now we are going to calculate the rotation matrices corresponding
           % to the site where we want to push particles to from this B site to so that we can
           % construct the appropriate Kraus operators
           j3 = j - 1;
           i3 = mod(i+1,Li);
           % Generate the rotation matrix for pushing particles from the site
           % where we want to push particles so that they are pushed to the
           % site below
           ind1 = [1 i3 j3];
           ind2 = [2 i j];
           [rotmat,measmat] = PresenceRevealed2(Li,Lj,ntimes,ind1,ind2);
           rotmats(:,:,aph,5) = rotmat;
           % Rotation matrix for pushing particles to the right
           ind1 = [1 i3 j3];
           ind2 = [2 i3 j3];
           [rotmat,measmat] = PresenceRevealed2(Li,Lj,ntimes,ind1,ind2);
           rotmats(:,:,aph,6) = rotmat;
           % Rotation matrix for pushing particles to the left
           ind1 = [1 i3 j3];
           i2 = mod(i3-1,Li);
           ind2 = [2 i2 j3];
           [rotmat,measmat] = PresenceRevealed2(Li,Lj,ntimes,ind1,ind2);
           rotmats(:,:,aph,7) = rotmat;
           % Produce the rotation matrix for pushing particles up above this
           % site
           ind1 = [1 i3 j3];
           j2 = j3-1;
           ind2 = [2 i3 j2];
           if j2<(0)
               rotmats(:,:,aph,8) = eye(2^(nqubits));
           else
               [rotmat,measmat] = PresenceRevealed2(Li,Lj,ntimes,ind1,ind2);
               rotmats(:,:,aph,8) = rotmat;
           end
       end
   end
   % Count the total number of sites on the bottom half of the cylinder
   num = aph;
   % Iterate over all of the Floquet cycles
   for z = 1:N
       % Evolve the density matrix for the first driving step
       unitnow = expm(-1i*(H1+diag(magi))*(1+TimeDisorder1(z))*2*pi/5);
       for t = 2:ntimes
           unitnow = kron(unitnow,expm(-1i*(H1+diag(magi))*(1+TimeDisorder1(z))*2*pi/5));
       end
       density = unitnow*density*ctranspose(unitnow);
       % Evolve the density matrix for the secon driving step
       unitnow = expm(-1i*(H2+diag(magi))*(1+TimeDisorder2(z))*2*pi/5);
       for t = 2:ntimes
           unitnow = kron(unitnow,expm(-1i*(H2+diag(magi))*(1+TimeDisorder2(z))*2*pi/5));
       end
       density = unitnow*density*ctranspose(unitnow);
       % Evolve the density matrix for the third driving step
       unitnow = expm(-1i*(H3+diag(magi))*(1+TimeDisorder3(z))*2*pi/5);
       for t = 2:ntimes
           unitnow = kron(unitnow,expm(-1i*(H3+diag(magi))*(1+TimeDisorder3(z))*2*pi/5));
       end
       density = unitnow*density*ctranspose(unitnow);
       % Evolve the density matrix for the fourth driving step
       unitnow = expm(-1i*(H4+diag(magi))*(1+TimeDisorder4(z))*2*pi/5);
       for t = 2:ntimes
           unitnow = kron(unitnow,expm(-1i*(H4+diag(magi))*(1+TimeDisorder4(z))*2*pi/5));
       end
       density = unitnow*density*ctranspose(unitnow);
       % Evolve the density matrix for the fifth driving step
       unitnow = expm(-1i*(H5+diag(magi))*(1+TimeDisorder5(z))*2*pi/5);
       for t = 2:ntimes
           unitnow = kron(unitnow,expm(-1i*(H5+diag(magi))*(1+TimeDisorder5(z))*2*pi/5));
       end
       density = unitnow*density*ctranspose(unitnow);
       % Iterate over all of the sites on the bottom half of the cylinder
       for s = 1:num
           % Generate the Kraus operators for moving particles from the two sites
           % of interest (the site where we are trying to move particles to and
           % the site where we are trying to move particles away from) to the four
           % sites that surround each of these two sites.
           correctnows = zeros(2^(ntimes*nqubits),2^(ntimes*nqubits),4);
           aph = 0;
           for s1 = 1:8
               for s2 = 1:ntimes
                   aph = aph + 1;
                   if (s2==1)
                       correctnow = rotmats(:,:,s,s1);
                   else
                       correctnow = eye(2^nqubits);
                   end
                   for s3 = 2:ntimes
                       if (s3==s2)
                           correctnow = kron(correctnow,rotmats(:,:,s,s1));
                       else
                           correctnow = kron(correctnow,eye(2^nqubits));
                       end
                   end
                   correctnows(:,:,aph) = correctnow;
               end
           end
           num2 = aph;
           % Implement all of these Kraus operators
           densityb = (1-num2*gamma2)*eye(2^(ntimes*nqubits))*density*ctranspose(eye(2^(ntimes*nqubits)));
           for s1 = 1:num2
               densityb = densityb + gamma2*correctnows(:,:,s1)*density*ctranspose(correctnows(:,:,s1));
           end
           density = densityb;
           % Add an external qubit
           density = kron(density,[1 0; 0 0]);
           % Determine whether we want to move particles between the two sites
           % of interest by transfering this information through entanglement
           % with the external particle
           density = measmats(:,:,s,1)*density*ctranspose(measmats(:,:,s,1));
           % Transform the particles in the AFAI depending on the state of the
           % external particle
           correctnow = rotmats(:,:,s,1);
           for t = 2:ntimes
               correctnow = kron(correctnow,rotmats(:,:,s,1));
           end
           correction = zeros(2^(ntimes*nqubits+1),2^(ntimes*nqubits+1));
           correction = correction + kron(eye(2^(ntimes*nqubits)),[1 0; 0 0]);
           correction = correction + kron(correctnow,[0 0; 0 1]);
           density = correction*density*ctranspose(correction);
           % Separate the external qubit from the system
           [rdensity] = ReducedDensity(density,ntimes*nqubits+1,1:(ntimes*nqubits));
           density = rdensity;
           % Apply noise to the operation that is supposed to transfer
           % particles from the site that particles are supposed to be
           % transferred away from to the site that they are supposed to be
           % transferred to
           for t = 1:ntimes
               if (t==1)
                   matnow = rotmats(:,:,s,1);
               else
                   matnow = eye(2^(nqubits));
               end
               for t2 = 2:ntimes
                   if (t==t2)
                       matnow = kron(matnow,rotmats(:,:,s,1));
                   else
                       matnow = kron(matnow,eye(2^(nqubits)));
                   end
               end
               density = (1-gamma)*eye(2^(ntimes*nqubits))*density*ctranspose(eye(2^(ntimes*nqubits))) + gamma*matnow*density*ctranspose(matnow);
           end
           % Apply noise to the two sites of interest that transfers particles
           % from these two sites to the four sites that surround each othe
           % these two sites
           densityb = (1-num2*gamma2)*eye(2^(ntimes*nqubits))*density*ctranspose(eye(2^(ntimes*nqubits)));
           for s1 = 1:num2
               densityb = densityb + gamma2*correctnows(:,:,s1)*density*ctranspose(correctnows(:,:,s1));
           end
           density = densityb;
       end
       % Create a new density matrix from the original density matrix and
       % evolve this density matrix through the first driving step
       unitnow = expm(-1i*(H1+diag(magi))*(1+TimeDisorder1(z+1))*2*pi/5);
       for t = 2:ntimes
           unitnow = kron(unitnow,expm(-1i*(H1+diag(magi))*(1+TimeDisorder1(z+1))*2*pi/5));
       end
       density2 = unitnow*density*ctranspose(unitnow);
       % Evolve the newly created density matrix through the second driving
       % step
       unitnow = expm(-1i*(H2+diag(magi))*(1+TimeDisorder2(z+1))*2*pi/5);
       for t = 2:ntimes
           unitnow = kron(unitnow,expm(-1i*(H2+diag(magi))*(1+TimeDisorder2(z+1))*2*pi/5));
       end
       density2 = unitnow*density2*ctranspose(unitnow);
       % Generate the velocity matrices integrated over time for the first and
       % third driving steps
       PMatrix3 = zeros(LSquared);
       for t = 1:2:length(V1(1,:))
           u1 = magi(t);
           u2 = magi((t+1));
           PMatrix1(t:(t+1),t:(t+1)) = TwoDPxyMatrix(V1(t:(t+1),t:(t+1)),u1,u2,(1+TimeDisorder1(z+1))*2*pi/5,0,J,1,Li,Lj,magi);
       end
       clear t
       PMatrix3 = TwoDPxyMatrix([0 -1i*J; 1i*J 0],0,0,(1+TimeDisorder3(z+1))*2*pi/5,0,J,3,Li,Lj,magi);
       % Calculate the charge pumped for the first and third driving steps
       for s = 1:ntimes
           % The charge pumped for these driving steps is calculated using the
           % eigenvectors of the appropriate reduced density matrices that are
           % scaled using the appropriate eigenvalues
           [rdensity] = ReducedDensity(density,ntimes*nqubits,(((s-1)*nqubits+1):(s*nqubits)));
           [V,D] = eig(rdensity);
           lennow = length(V(1,:));
           for t = 1:lennow
               P1b = P1b + abs(D(t,t))*ctranspose(V(:,t))*PMatrix1*V(:,t);
           end
           [rdensity] = ReducedDensity(density2,ntimes*nqubits,(((s-1)*nqubits+1):(s*nqubits)));
           [V,D] = eig(rdensity);
           lennow = length(V(1,:));
           for t = 1:lennow
               P3b = P3b + abs(D(t,t))*ctranspose(V(:,t))*PMatrix3*V(:,t);
           end
       end
       % Determine the charge pumped for the current Floquet cycle
       if sum(z==NVec)
           QVecb = [QVecb real(P1b - P3b)/(2*Li)];
           Qb = [Qb (QVecb(end)-QVecb(end-1))/Step];
           save('Qb.mat','Qb')
       end
   end

The following helper function generate the matrices that implements the chemical potential disorder when added to the Hamiltonians.

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
        % Generate a random number that is drawn between -W and W where W is 
        % given by the value of chaos.
        candy = -chaos + 2*chaos*rand;
        % Use the random number to apply a random on-site potential
        ioph(i) = candy;
        clear candy
    end
    % Return the resulting matrix as output.
    magi = ioph;
    end

Here is the helper function that generates the Hamiltonians for the five driving steps as well as the velocity matrices for the first and third driving steps.

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

This is the helper function that generates the velocity matrices (for both the first and third driving steps) integrated over the appropriate time intervals.

.. code-block:: matlab

    function y = TwoDPxyMatrix(VelMat,u1,u2,tf,ti,J,step,Li,Lj,magi)
    % This function is used to calcuate the charge pumped integrated over a certain time
    % period so that the computationally expensive method of numerical
    % integration is unneeded. Here, VelMat is the velocity matrix of interest,
    % u1 and u2 are the terms that define the chemical potential disorder, ti
    % is the starting time, tf is when the driving step of interest is over, J
    % is the hopping strength, step defines whether the driving step of
    % interest is the first or third driving step, and magi is the full matrix
    % that implements the chemical potential disorder. Unfortunately, I have
    % lost the notes that derive the math to form this algorithm, which is why
    % I am showing plots that justify that this is in fact the correct way
    % integrating the expectation value of the velocity matrix over time. If it
    % is deemed necessary for me to find the notes or rederive them, then I
    % will do that.
    LSquared = 2*Li*Lj;
    if step == 1
        SigY = [0 -1i; 1i 0];
        theta = atan(-2*J/(u1-u2));
        Ry = expm(-1i*SigY*theta/2);
        B = ctranspose(Ry)*VelMat*Ry;
        Sec = (u1-u2)*cos(theta)/2-J*sin(theta);
        y = Ry*([B(1,1)*tf (B(1,2)/(2*1i*Sec))*exp(2*1i*Sec*tf); (-B(2,1)/(2*1i*Sec))*exp(-2*1i*Sec*tf) B(2,2)*tf] - [B(1,1)*ti (B(1,2)/(2*1i*Sec))*exp(2*1i*Sec*ti); (-B(2,1)/(2*1i*Sec))*exp(-2*1i*Sec*ti) B(2,2)*ti])*ctranspose(Ry);
    elseif step == 3
        Muy = zeros(LSquared);
        for i = 0:(Li-1)
            for j = 0:(Lj-1)
                u3 = magi((1+2*rem((i+1),Li)+2*Li*j));
                u4 = magi((2+2*i+2*Li*j));
                SigY = [0 -1i; 1i 0];
                theta = atan(-2*J/(u3-u4));
                Ry = expm(-1i*SigY*theta/2);
                B = ctranspose(Ry)*VelMat*Ry;
                Sec = (u3-u4)*cos(theta)/2-J*sin(theta);
                Result = Ry*([B(1,1)*tf B(1,2)*exp(2*1i*Sec*tf)/(2*1i*Sec); -B(2,1)*exp(-2*1i*Sec*tf)/(2*1i*Sec) B(2,2)*tf] - [B(1,1)*ti B(1,2)*exp(2*1i*Sec*ti)/(2*1i*Sec); -B(2,1)*exp(-2*1i*Sec*ti)/(2*1i*Sec) B(2,2)*ti])*ctranspose(Ry);
                Muy((1+2*rem((i+1),Li)+2*Li*j),(1+2*rem((i+1),Li)+2*Li*j)) = Result(1,1);
                Muy((2+2*i+2*Li*j),(2+2*i+2*Li*j)) = Result(2,2);
                Muy((1+2*rem((i+1),Li)+2*Li*j),(2+2*i+2*Li*j)) = Result(1,2);
                Muy((2+2*i+2*Li*j),(1+2*rem((i+1),Li)+2*Li*j)) = Result(2,1);
            end
        end
        y = Muy;
    end
    end

Finally, we have the helper function that calculates the reduced density matrix given an input matrix as well as some target qubits that will define the reduced density matrix.

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
        else
        % Otherwise make the first matrix that composes this set of Kronecker
        % products be the appropriate single qubit spin vector
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
