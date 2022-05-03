%% Implementation of the R = diag(1,1,-1) gate unitarily in qutrit Clifford+T
% accompanying the paper "Qutrit metaplectic gates are a subset of Clifford+T"
% Andrew Glaudell, Neil J. Ross, John van de Wetering, Lia Yeh

%% identity and zero matrices
I = eye(3);     % single-qutrit identity gate, i.e. 3 x 3 identity matrix
II = eye(3^2);  % two-qutrit identity gate
zero = zeros([3]);  % 3 x 3 zero matrix

%% Z and X basis states
k0 = I(:,1);    % Z-basis states
k1 = I(:,2);    % "
k2 = I(:,3);    % "
w = exp(2*pi*i/3);      % 3rd root of unity
kpl = i * (k0 + k1 + k2)/sqrt(3);               % X-basis states // Eqs.3-5
kw = i * (k0 + w * k1 + w^2 * k2)/sqrt(3);      % "
kwsq = i * (k0 + w^2 * k1 + w * k2)/sqrt(3);    % "

%% qutrit generalization of Pauli gates
x = [0 0 1; 1 0 0; 0 1 0];    % X, i.e. tau(0 1 2) // Def.1
z = [1 0 0; 0 w 0; 0 0 w^2];  % Z = H * X * H' // Def.1

%% qutrit S, H, and CX, which generate the Clifford group
s = [1 0 0; 0 1 0; 0 0 w];  % S // Def.7
h = (kron(kpl, k0') + kron(kw, k1') + kron(kwsq, k2'));   % H // Def.8
cx = [I zero zero; zero x zero; zero zero x'];   % CX, i.e. tau(10 11 12)(20 22 21) // Def.10

%% various other Clifford gates
tau0_1 = [0 1 0; 1 0 0; 0 0 1]; % tau(0 1)  // Sec.2.1
tau0_2 = [0 0 1; 0 1 0; 1 0 0]; % tau(0 2)  // Sec.2.1
tau1_2 = [1 0 0; 0 0 1; 0 1 0]; % tau(1 2) = H * H up to a global phase of -1 // Sec.2.1; Def.8
z01 = diag([1,1,w]);            % Z(0,1)    // Def.6
z22 = diag([1,w^2,w^2]);        % Z(2,2)    // Def.6
swap = II(:, [1 4 7 2 5 8 3 6 9]);  % two-qutrit SWAP gate

%% the T gate
t = diag([1,w^(1/3),w^(-1/3)]); % T // Def.13

%%  Controlled gates from the paper "Factoring with Qutrits: Shor’s
%   Algorithm on Ternary and Metaplectic Quantum Architectures"
zcz = cx * kron(I,t) * cx * kron(I,t) * cx * kron(I,t); % |0>-controlled Z // Eq.10
ocx = kron(x,h') * zcz * kron(x',h); % tau(10 11 12)
tcx = kron(x,I) * ocx * kron(x',I); % tau(20 21 22)
socxs = swap * ocx * swap; % tau(01 11 21)
tau02_20 = swap * socxs * ocx * socxs * ocx * socxs; % tau(02 20)
map21_22to02_20 = swap * cx' * swap * kron(I,tau0_1);
tctau1_2 = map21_22to02_20' * tau02_20 * map21_22to02_20;   % |2>-controlled tau(1 2), i.e. tau(21 22) // Lem.17

%% |2>-controlled w^(-1/3) * S gate where S = Z(0,1)
tcsphase = tcx' * kron(I,tau0_1*t*tau0_1) * tcx * kron(I,tau0_1*t'*tau0_1); % // Lem.18
%% |2>-controlled w^(2/3) * Z(2,2) gate
tcz22phase = kron(I,tau0_2) * tcsphase * kron(I,tau0_2); % Cor.19
%% |2>-controlled -H gate
tcmh = tcz22phase * kron(I,h')*tcz22phase*kron(I,h) * tcz22phase * kron(s,I);   % // Lem.20
%% |2>-controlled -tau(1 2) gate, i.e. tctau1_2 but when the control is |2>, the target gains a -1 phase
tcmtau1_2 = kron(z01,z22*h') * tcz22phase * kron(I,h*z22) * tcz22phase * kron(I,h') * tcz22phase * kron(I,h*z22); % // Thm.21

%% the below two-qutrit Clifford+T unitary is the R gate on one qutrit and the identity gate on the other
rtensorid = tcmtau1_2 * tctau1_2;
%% to verify: subtracting the R gate tensor the single-qutrit identity from it equals zero
r = diag([1,1,-1]); % R gate, i.e. the metaplectic gate
assert(all(all(round(rtensorid,10) - kron(r,I) == 0))) % rounding to the 10^(-10) decimal place due to floating point error
