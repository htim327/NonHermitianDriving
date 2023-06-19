===============================================================
Presenting the code for the Anomalous Floquet-Anderson Insulator
===============================================================

This file presents how the MATLAB code is set up for simulating the Anomalous Floquet-Anderson Insulator and shows that the method used for calculating the pumped charge integrated over the appropriate time intervals using analytical integration approaches produces the correct results when compared to using numerical integration approaches. We start with the main file that runs the appropriate files and then plots the data.

.. code-block:: matlab

   clear;
   clc;
   % Store the error between the numerical and analytical integration
   % approaches used for calculating the charge pumped in errormat
   errormat = [];
   % Define how many noise realizations we use
   NTot = 100;
   % Define how accurate we want our numerical integration approach to be. The
   % higher the number, the more accurate.
   iter = 100;
   save('iter.mat','iter')
   for zamp = 1:NTot
       tic
       % Run TwoDimxyQ
       TwoDimxyQ
       % Load the information related to the charge pump calculated using the
       % numerical and analytical integration approaches
       load(['QTot_' num2str(iter) '.mat'])
       % Calculate the error between the analytical and numerical integration
       % approaches
       errormat = [errormat; (QTot(1,:)-QTot(2,:))];
       clearvars -except zamp errormat NTot iter
       % If you're impatient, the following lines let you become more patient.
       clc
       zamp
       toc
   end
   % Plot the data
   figure('units','normalized','outerposition',[0 0 1 1]);
   errorbar(1:length(errormat(1,:)),sum(errormat)/NTot,std(errormat)/sqrt(NTot),'LineWidth',2,'Color','m')
   title(['Error for iter=' num2str(iter)],'FontSize',40,'Interpreter','latex')

We now present the code that makes up the TwoDimxyQ.m file is run during each iteration in the file above.

.. code-block:: matlab

    % This file calculates the charge pumped using both the numerical
    % integration and analytical integration approaches
    %%%
    % Define the number of unit cells in both the i (x) and j (y) directions
    Li = 2;
    Lj = 4;
    % Define the number of sites of the system.
    LSquared = 2*Li*Lj;
    % Define the strength of the chemical potential implemented during the
    % fifth driving step
    del = 0.4;
    % Define the strength of the chemical potential disorder
    Noise = 1.5;
    % Define the strength of the temporal disorder
    tchaos = 0.2;
    % Define the strength of the hopping term
    J = 1.25;
    % Determine the number of Floquet cycles implemented
    NVec = 1:100;
    N = max(NVec);
    % Use a random seed for the random number generator
    rng('shuffle');
    % Generate the matrix that implements the chemical potential disorder
    magi = TwoDxyNoiseHamiltonians(Li,Lj,Noise);
    % Set up the Hamiltonians for your system
    [H1, H2, H3, H4, H5, V1, V3] = FastTwoDxyHamiltonians(Li,Lj,J,del);
    % Set up the wave functions that are used for our systems
    W = eye(LSquared);
    wave = W(:,1:round(LSquared/2));
    rng('shuffle');
    % Set up the variables used to implement the temporal disorder
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
    % Store the information related to how much charge is pumped in the first
    % and third driving steps as well as the total charge pumped during each
    % Floquet cycle
    P1 = 0;
    P3 = 0;
    QVec = [0];
    Q = [];
    P1a = 0;
    P3a = 0;
    QVeca = [0];
    Qa = [];
    % Load the information related to how accurate we want our numerical
    % integration to be.
    load('iter.mat')
    % Iterate over all of the Floquet cycles
    for z = 1:N
        % Time evolve the system to the current Floquet cycle
        wave2 = expm(-1i*(H5+diag(magi))*(1+TimeDisorder5(z))*2*pi/5)*expm(-1i*(H4+diag(magi))*(1+TimeDisorder4(z))*2*pi/5)*expm(-1i*(H3+diag(magi))*(1+TimeDisorder3(z))*2*pi/5)*expm(-1i*(H2+diag(magi))*(1+TimeDisorder2(z))*2*pi/5)*expm(-1i*(H1+diag(magi))*(1+TimeDisorder1(z))*2*pi/5)*wave2;
        % Time evolve the system to the beginning of the third driving step of
        % the current Floquet cycle
        wave3 = expm(-1i*(H2+diag(magi))*(1+TimeDisorder2(z+1))*2*pi/5)*expm(-1i*(H1+diag(magi))*(1+TimeDisorder1(z+1))*2*pi/5)*wave2;
        % Generate the time evolution matrices used for the numerical
        % integration
        Unit1 = expm(-1i*(H1+diag(magi))*(1+TimeDisorder1(z+1))*2*pi/(5*iter));
        Unit3 = expm(-1i*(H3+diag(magi))*(1+TimeDisorder3(z+1))*2*pi/(5*iter));
        % Generate the matrices PMatrix1 and PMatrix3 that are used for the
        % analytical integration of the pumped charge for the first and third
        % driving steps
        PMatrix3 = zeros(LSquared);
        for t = 1:2:length(V1(1,:))
            u1 = magi(t);
            u2 = magi((t+1));
            PMatrix1(t:(t+1),t:(t+1)) = TwoDPxyMatrix(V1(t:(t+1),t:(t+1)),u1,u2,(1+TimeDisorder1(z+1))*2*pi/5,0,J,1,Li,Lj,magi);
        end
        PMatrix3 = TwoDPxyMatrix([0 -1i*J; 1i*J 0],0,0,(1+TimeDisorder3(z+1))*2*pi/5,0,J,3,Li,Lj,magi);
        % Iterate through all of the wave functions
        for s = 1:length(wave(1,:))
            % Calculate the information related to the analytically derived
            % calculation of the integrated current for the first and third
            % driving steps
            P1 = P1 + ctranspose(wave2(:,s))*PMatrix1*wave2(:,s);
            P3 = P3 + ctranspose(wave3(:,s))*PMatrix3*wave3(:,s);
            % Set up the matrices that are used for the time evolution used for
            % the numerical integration
            Unit1a = eye(length(W(1,:)));
            Unit3a = eye(length(W(1,:)));
            % Iterate over the number of steps used for the numerical
            % integration
            for t = 1:iter
                % Set up the time evolution matrices corresponding to the
                % current iteration
                Unit1a = Unit1*Unit1a;
                Unit3a = Unit3*Unit3a;
                % Calculate the information related to the charge pumped using
                % the numerical integration approach for both the first and
                % third driving steps.
                P1a = P1a + ctranspose(wave2(:,s))*ctranspose(Unit1a)*V1*Unit1a*wave2(:,s)*2*pi*(1+TimeDisorder1(z+1))/(5*iter);
                P3a = P3a + ctranspose(wave3(:,s))*ctranspose(Unit3a)*V3*Unit3a*wave3(:,s)*2*pi*(1+TimeDisorder3(z+1))/(5*iter);
            end
        end
        if sum(z==NVec)
            % Calculate the total charge pumped using the numerical integration
            % approach
            QVeca = [QVeca real(P1a - P3a)/(2*Li)];
            % Calculate the total charge pumped using the analytical integration
            % approach
            QVec = [QVec real(P1 - P3)/(2*Li)];
            % Calculate the charge pumped during the current Floquet cycle
            % using the analytical integration approach
            Q = [Q (QVec(end)-QVec(end-1))];
            % Calculate the charge pumped during the current Floquet cycle
            % using the numerical integration approach
            Qa = [Qa (QVeca(end)-QVeca(end-1))];
            save('Q.mat','Q')
            save('Qa.mat','Qa')
        end
    end
    % Store the information related to the charge pumped for each Floquet cycle
    % in QTot
    QTot = [Q; Qa];
    save(['QTot_' num2str(iter) '.mat'],'QTot')

Then we have the helper function that creates the matrices that when added to the Hamiltonians, implements the chemical potential disorder.

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

For the Hamiltonians, what we do is divide the cylindrical lattice into unit cells such that each unit cell has two sites in the x (i) direction and one site in the y (j) direction. The leftmost site of each unit cell can be the A site and the rightmost site of the unit cell can be the B site. This division is important for the implementation of the chemical potential, which is done in step five, where site A is evolved to have the opposite phase added as that of site B. It is important to remember that the A sites are always surrounded by B sites in all four directions and the B sites are always surrounded by A sites in all four directions. The wave function is set up to have indices :math:`$|\Psi(\alpha+2\times i+ 2\times \mathrm{L}_i \times j)\rangle$`, where :math:`$\alpha$` can be 1 or 2 depending on whether we are referring to an A or B site, respectively, :math:`$i$` defines the unit cell of interest in the i direction, :math:`$\mathrm{L}_i$` defines how many unit cells there are in the i direction, and :math:`$j$` defines the unit cell of interest in the j direction. 

If this is the case and :math:`$\mathrm{L}_\mathrm{squared}$` defines the number of sites of the system, then the Hamiltonian for the Floquet driving step five is given by:

:math:`$$H_5 = \sum_{k=1}^{\mathrm{L}_\mathrm{squared}} (-1)^{k-1} \times \Delta$$`

where :math:`$\Delta$` defines the strength of the chemical potential. Meanwhile, if we rewrite our wave function expressed above as :math:`$\Psi(i,j,\alpha)$`, then the Hamiltonian for the first four driving steps can be expressed as:

:math:`$$H_{1-4} = -J \sum_{i,j} (|i,j,1\rangle\langle i+a, j+b, 2| + h.c.)$$`

where for :math:`$H_1$`, :math:`$a=b=0$`, for :math:`$H_2$`, :math:`$a = -1$` and :math:`$b = 1$`, for :math:`$H_3$` :math:`$a = -1$` and :math:`$b = 0$`, and for :math:`$H_4$`, :math:`$a = 0$` and :math:`$b = -1$`. The function that implements this is given by:

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

Finally, we have the function that gives us the analytical calculation of the charge pumped integrated over a certain time interval:

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

If we go back to the file that we presented in the beginning, we are going to see what happens when we set :math:`$iter=100$`, :math:`$iter=1000$`, and then :math:`$iter=10000$`. For :math:`$iter=100$`, we have:

