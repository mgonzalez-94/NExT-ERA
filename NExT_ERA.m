function [Wn,zn,Ph,ModelOrder] = NExT_ERA(y,fs,RefCh,MxLgs,Windw,Nvrlp,N,NModes,FnApr)
%
% ID = NExT_ERA(y,fs,RefCh,MxLgs,Windw,Nvrlp,N,NModes,FnApr,ZMax)
%
% Función para identificar las propiedades dinámicas de un sistema
% empleando NExT y ERA.
%
% INPUTS:
%   y: señales en tiempo con tantas filas como canales y columnas como
%       muestras.
%   fs: frecuencia de muestreo (Hz).
%   RefCh: canal de referencia indicado como escalar o como vector.
%   MxLgs: cantidad de muestras para NExT en dominio de tiempo.
%   Windw: número de puntos del ventaneo a implementar en NExT en domionio
%       de frecuencia.
%   Nvrlp: traslapo indicado de 0 a 1 a implementar en NExT en domionio de
%       frecuencia.
%   N: cantidad de ventanas a implementar en NExT en domionio de
%       frecuencia.
%   NModes: cantidad de modos a identificar especificados como un escalar o
%       un vector.
%   FnApr: matriz con rangos de frecuencia a proximada, sino [].
%
% OUTPUTS:
%   Wn,z,Ph,ModelOrder.
%
% Mateo González H.
% 2020/05/02.
% Adapatado de:
% https://www.mathworks.com/matlabcentral/fileexchange/72170-next-era.
%
tic_NExT_ERA = tic;
%%% -----------------------------------------------------------------------
MxLgs = min([MxLgs,length(y(1,:))-1]);
NCols = floor((2*Windw*Nvrlp)/3);
NCols = min([NCols,ceil(Windw/2+1)-1]);
Shift = 10;
EMACO = 1;
Wn = [];
zn = [];
Ph = [];
ModelOrder = [];
%%% -----------------------------------------------------------------------
for RefCh_i = 1:length(RefCh)
    %%% -------------------------------------------------------------------
    for NModes_i = 1:length(NModes)
        %%% ---------------------------------------------------------------
        NRows = 20*NModes(NModes_i);
        Cut   = 2*NModes(NModes_i);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% NExT-ERA en Frecuencia %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        IRF = NExTF(y,RefCh(RefCh_i),Windw,N,Nvrlp);
        Inp = length(RefCh(RefCh_i));
        ID1 = ERA(IRF,fs,NCols,NRows,Inp,Cut,Shift,EMACO);
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% NExT-ERA en Tiempo %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        IRF = NExTT(y,RefCh(RefCh_i),MxLgs);
        Inp = length(RefCh(RefCh_i));
        ID2 = ERA(IRF,fs,NCols,NRows,Inp,Cut,Shift,EMACO);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Seleccionar Resultados %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Wn_i = [ID1.Wn(:);ID2.Wn(:)];
        zn_i = [ID1.z(:);ID2.z(:)];
        Ph_i = [ID1.Ph,ID2.Ph];
        Wn = [Wn;Wn_i];
        zn = [zn;zn_i];
        Ph = [Ph,Ph_i]; %#ok<*AGROW>
        ModelOrder = [ModelOrder(:);NModes(NModes_i)*ones(length(Wn_i),1)];
        %%% ---------------------------------------------------------------
    end
    %%% -------------------------------------------------------------------
end
%%% -----------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%
%%% DismissThis %%%
%%%%%%%%%%%%%%%%%%%
%%% zn>0.2
DismissThis = false(length(Wn),1);
DismissThis(zn>0.2) = true;
%
Wn(DismissThis)   = [];
zn(DismissThis)   = [];
Ph(:,DismissThis) = [];
%%% Wn./(2*pi)>0.4*fs
DismissThis = false(length(Wn),1);
DismissThis(Wn./(2*pi)>0.4*fs) = true;
%
Wn(DismissThis)   = [];
zn(DismissThis)   = [];
Ph(:,DismissThis) = [];
%%% Wn<=1e-4
DismissThis = false(length(Wn),1);
DismissThis(Wn<=1e-4) = true;
%
Wn(DismissThis)   = [];
zn(DismissThis)   = [];
Ph(:,DismissThis) = [];
%%% ~( Fn>(FnApr(:,1)') & Fn<(FnApr(:,2)') )
if ~isempty(FnApr)
    Fn = Wn./(2*pi);
    DismissThis = ~any(Fn>(FnApr(:,1)') & Fn<(FnApr(:,2)'),2);
    Wn(DismissThis)   = [];
    zn(DismissThis)   = [];
    Ph(:,DismissThis) = [];
end
%%% -----------------------------------------------------------------------
disp(['NExT_ERA: ',num2str(toc(tic_NExT_ERA),'%.2f')])
end
function ID = ERA(Y,fs,NCols,NRows,Inp,Cut,Shift,EMACO)
%
% ID = ERA(Y,fs,NCols,NRows,Inp,Cut,Shift,EMACO)
%
% INPUTS:
%   Y: Free vibration output data in a form of Y=[Y1 Y2 ... Y_Ndata] Yi is
%       Markov Parameter of size (outputs,inputs) and the total size is
%       (outputs,inputs*Ndata). Were outputs is the number of output
%       channels, inputs is the number of inputs which equals to 1 unless
%       free vibration data comes from Multi-reference channels NExT.
%       Ndata is the length of the data samples
%   fs: Sampling frequency.
%   NCols: The number of columns in hankel matrix (more than 2/3 of No. of
%       data).
%   NRows: The number of rows in hankel matrix (more than 20 * number of
%       modes).
%   Inp: The number of inputs which equals to 1 unless free vibration data
%       comes from Multi-reference channels NExT.
%   Cut: cutoff value=2*no of modes.
%   Shift: Shift value in the final row and column blocks (Increase EMAC
%       sensitivity) usually =10.
%   EMACO: if this value equals to 1, EMAC will be independent of the
%       number of columns (calculated only from observability matrix not from
%       controllability).
%
% OUTPUTS:
%   ID: A structure consist of the below components.
%   Parameters: NaFreq : Natural frequencies vector.
%   DampRatio: Damping ratios vector.
%   ModeShape: Mode shape matrix.
%   Indicators: MAmC : Modal Amplitude Coherence.
%   EMAC: Extended Modal Amplitude Coherence.
%   MPC: Modal Phase Collinearity.
%   CMI: Consistent Mode Indicator.
%   partfac: Participation factor.
%   Matrices A,B,C: Discrete A,B and C matrices.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Verificar Dimensiones Y %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nc,nm] = size(Y);       % (número de canales y número de muestras)
if nc > nm             % Check if Y matrix size is proper or should be transposed
    Y = Y';
    [nc,~] = size(Y);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Determinar Parámetros %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BRows = fix(NRows/nc);  % BRows = how many output blocks. ¿max([1, floor(NRows/nc)])?
NRows = nc*BRows;       % Redefine the row numbers.
BCols = fix(NCols/Inp); % BCols = how many time steps.
NCols = Inp*BCols;      % Redefine the column numbers.
m = Inp;                % Inputs number.
q = nc;                 % Outputs number.
%%%%%%%%%%
%%% H0 %%% ¿?
%%%%%%%%%%
H0=zeros(NRows,NCols);
for ii = 1:BRows
    for jj = 1:BCols
        if ii == BRows || jj == BCols
            sh = Shift;
        else
            sh = 1;
        end
        if ii == BRows && jj == BCols
            sh = 2*Shift-1;
        end
        H0((ii-1)*q+1:ii*q,(jj-1)*m+1:jj*m) = Y(:,(sh-1+(jj-1)+ii-1)*m+1:(sh-1+(jj)+ii-1)*m);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Descomposición de H0 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[P,D,Q] = svd(H0,0);
Dn = diag(sqrt(diag(D(1:Cut,1:Cut))));
Pn = P(:,1:Cut);
Qn = Q(:,1:Cut);
%%%%%%%%%%
%%% H1 %%%
%%%%%%%%%%
H1=zeros(NRows,NCols);
for ii=1:BRows
    for jj=1:BCols
        if (ii==BRows || jj==BCols)
            sh=Shift;
        else
            sh=1;
        end
        if ii==BRows && jj==BCols
            sh=2*Shift-1;
        end
        H1((ii-1)*q+1:ii*q,(jj-1)*m+1:jj*m)=Y(:,(sh-1+(jj-1)+ii)*m+1:(sh-1+(jj)+ii)*m);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Espacio de Estados %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
A = Dn\Pn'*H1*Qn/Dn;
B = Dn*Qn';
B = B(:,1:m);
C = Pn*Dn;
C = C(1:q,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Clean Workspace so far %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear H0 H1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Propiedades de la Matriz de Estados %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Ph,W2] = eig(A);                   % Valores y vectores propios
Lambda  = diag(W2);                 % Raíces de Laplace en el plano Z.
s       = log(Lambda).*fs;        	% Raíces de Laplace.
Zn      = (-1)*real(s)./abs(s)*100;	% Razones de amortiguamiento.
fd      = (imag(s)./2./pi);         % Frecuencias naturales a mortiguadas.
Shapes  = C*Ph;                     % Mode shapes.
PartFac = Ph\Dn*Qn(1:m,:)';         % Facotores de participación.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate MAmC Values %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bhat  = Ph\B;
qhati = zeros(1,NCols);
qhat  = zeros(Cut,NCols);
for ii = 1:Cut
    for jj = 1:BCols
        qhati(1,(jj-1)*m+1:jj*m) = Bhat(ii,:)*Lambda(ii)^(jj-1);
    end
    qhat(ii,:) = qhati;
end
Qbar = Ph\Dn*Qn';
qbar = zeros(Cut,length(Qbar(1,:)));
for ii = 1:Cut
    qbar(ii,:) = Qbar(ii,:);
end
EMAC = zeros(Cut);
for ii = 1:Cut
    qhati = qhat(ii,:);
    qbari = qbar(ii,:);
    EMAC(ii,1) = (abs(qbari*qhati')/sqrt(abs(qbari*qbari')*abs(qhati*qhati')))*100;
end
%--------------------------------------------------------------------------
% Calculate EMAC values
% Con=Dn*Qn';   Con=Con(:,1:m);    %Con=Con(:,ncols-2*m+1:ncols-m);
partinput22  = Ph\A^(BCols+Shift-2)*B;
partinput2   = Ph\Dn*Qn'; partinput2=partinput2(:,NCols-m+1:NCols);
partoutput22 = (C*A^(BRows+Shift-2)*Ph)';   %partoutput22=(Rn(nrows-2*q+1:nrows-q,:)*Dn*A^shift*Ph)';
partoutput2  = (Pn(NRows-q+1:NRows,:)*Dn*Ph)';
EMACC = zeros(Cut,1);
for i=1:1:Cut
    
    sum1=0;
    sum2=0;
    
    for j=1:1:m
        for k=1:1:q
            
            Rij=abs(partinput22(i,j))/abs(partinput2(i,j));
            if Rij>1
                Rij=1/Rij;
            end
            Wij=1-(4/pi)*abs(angle(partinput22(i,j)/partinput2(i,j)));
            
            if Wij<0
                Wij=0;
            end
            
            EMACin=Rij*Wij;
            
            if EMACO==1
                EMACin=1; % we assume that in order to make value independent of number of columns
            end
            
            Rik=abs(partoutput22(i,k))/abs(partoutput2(i,k));
            if Rik>1
                Rik=1/Rik;
            end
            Wik=1-(4/pi)*abs(angle(partoutput22(i,k)/partoutput2(i,k)));
            
            if Wik<0
                Wik=0;
            end
            
            EMACout=Rik*Wik;
            
            EMACijk=EMACin*EMACout;
            sum1=sum1+EMACijk*(abs(partinput2(i,j)))^2*(abs(partoutput2(i,k)))^2*100;
            sum2=sum2+(abs(partinput2(i,j)))^2*(abs(partoutput2(i,k)))^2;
        end
    end
    EMACC(i)=sum1/sum2;
end
%--------------------------------------------------------------------------
% Calculate MPC values
MPC = zeros(ii,1);
for ii=1:1:Cut
    cd=sum(Shapes(:,ii))/q;
    cj=Shapes(:,ii);
    Sxx=(real(cj-cd))'*(real(cj-cd));
    Sxy=(real(cj-cd))'*(imag(cj-cd));
    Syy=(imag(cj-cd))'*(imag(cj-cd));
    mu=(Syy-Sxx)/(2*Sxy);
    Beta=mu+sign(Sxy)*sqrt(mu^2+1);
    Tau=atan(Beta);
    l1=Sxx+(Sxy*(2*(mu^2+1)*(sin(Tau))^2-1))/mu;
    l2=Syy-(Sxy*(2*(mu^2+1)*(sin(Tau))^2-1))/mu;
    MPC(ii)=100*(2*(l1/(l1+l2)-0.5))^2;
end
%--------------------------------------------------------------------------
% Sort into ascending order
[fd,I]  = sort(fd);
Zn      = Zn(I);
Shapes  = Shapes(:,I);
PartFac = PartFac(I,:);
EMAC    = EMAC(I);
EMACC   = EMACC(I);
MPC     = MPC(I);
%--------------------------------------------------------------------------
% Remove the negative frequencies and frequencies>fs/2
lower=1;
upper=Cut;
for ii=1:Cut
    if fd(ii) <= 0
        lower=ii+1;
    end
    if fd(Cut-ii+1) >= 0.499*fs
        upper=Cut-ii;
    end
end
fd1=fd(lower:upper);
zeta1=Zn(lower:upper);
fd1=fd1./sqrt(1-(zeta1/100).^2);
Shapes=Shapes(:,lower:upper);
PartFac=PartFac(lower:upper,:);
EMAC1=EMAC(lower:upper);
EMACC1=EMACC(lower:upper);
MPC1=MPC(lower:upper);
%--------------------------------------------------------------------------
[fd1,ii]=sort(fd1);
zeta1=zeta1(ii);
EMAC1=EMAC1(ii);
EMACC1=EMACC1(ii);
MPC1=MPC1(ii);
Shapes=Shapes(:,ii);
PartFac=PartFac(ii,:);
%--------------------------------------------------------------------------
Phi = Shapes;
HH  = size(Phi);
[~,II]     = max(abs(Phi),[],1);
ModeShapeS = zeros(size(Phi));
for jj = 1:HH(2)
    b = -angle(Phi(II(jj),jj));
    ModeShapeS(:,jj) = real(Phi(:,jj)*exp(1i*b));
    ModeShapeS(:,jj) = ModeShapeS(:,jj)/norm(ModeShapeS(:,jj));
end
Shapes=ModeShapeS;
%--------------------------------------------------------------------------
MAmC=EMAC1;
EMAC=EMACC1';
MPC=MPC1';
CMI=EMAC.*MPC/100;
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%
%%% Output format %%%
%%%%%%%%%%%%%%%%%%%%%
NaFreq    = fd1;          % Natural frequencies vector
DampRatio = zeta1/100;    % Damping ratios vector
ModeShape = Shapes;       % Mode shape matrix
%%%
ID.Fn  = NaFreq;
ID.Wn  = NaFreq*2*pi;
ID.z   = DampRatio;
ID.Ph  = ModeShape;
ID.Ass = A;
ID.Bss = B;
ID.Css = C;
ID.Indicators.MAmC  = MAmC;
ID.Indicators.EMAC  = EMAC;
ID.Indicators.MPC   = MPC;
ID.Indicators.CMI   = CMI;
ID.Indicators.partfac = PartFac;
end
function IRF = NExTF(y,RefCh,Windw,N,Nvrlp)
%
% IRF = NExTF(y,RefCh,Window,N,p)
%
% INPUTS:
%	y: matriz de señales con tantas filas como canales y columnas como
%       muestras.
%  	RefCh: canal de referencia.
%  	Window: tamaño de la ventana para realizar la densidad espectral.
%  	N: número de ventanas.
%  	p: translapo entre ventanas indicado ente 0 (0%) y 1 (100%).
%
% OUPUTS:
%  	IRF: respuesta impulsiva con tantas filas como canales y columnas como MxLgs+1.
%
y(:,Windw*Nvrlp*(N+1)) = 0;         % En caso que el valor exceda la dimensión de y.
y  = y(:,1:Windw*Nvrlp*(N+1));      % Tomar solo los valores de interés.
nc = size(y,1);                 	% Número de canales.
Y  = zeros(nc,ceil(Windw/2+1)-1);   % Impulse Response Function matrix.
for jj=1:nc	% Loop for all channels.
    %%% Generation of vectors before applying cross-correlation to them
    yr = y(RefCh,:);	% Señal de referencia.
    yc = y(jj,:);       %
    wn = ceil((1-Nvrlp)*Windw);
    for i=1:N
        ind_i = 1+wn*(i-1);
        ind_f = Windw+wn*(i-1);
        xx = yr(ind_i:ind_f);
        yy = yc(ind_i:ind_f);
        M1 = fft(xx);
        M2 = fft(yy);
        MX1(i,:) = M1;
        MX2(i,:) = M2;
    end
    MX  = (ones(1,N)*(MX1.*conj(MX2)))/N;   % Espectro de densidad cruzada.
    XCF = ifft(MX);                         % Respuesta en vibración libre.
    XCR = XCF(1:ceil(Windw/2+1)-1);
    MaxLags = length(XCR)-1;
    Y(jj,1:MaxLags+1) = XCR;                % Array for cross-correlation.
end
IRF = Y;
end
function IRF = NExTT(y,RefCh,MxLgs)
%
% IRF = NExTT(y,RefCh,MxLgs)
%
% INPUTS:
% 	y: matriz de señales con tantas filas como canales y columnas como
%       muestras.
%  	RefCh: canal de referencia.
%  	MxLgs: número de muestras en la función de correlación cruzada.
%
% OUPUTS:
% 	IRF: respuesta impulsiva con tantas filas como canales y columnas como
% 	ceil(Windw/2+1)-1.
%
nc = size(y,1);
Y  = zeros(nc,MxLgs+1);
for jj = 1:nc
    yr  = y(RefCh,:);
    yc  = y(jj,:);
    XCR = xcorr(yr,yc,MxLgs,'unbiased');
    XCF = XCR(MxLgs+1:2*MxLgs+1);
    Y(jj,1:1+MxLgs) = XCF;
end
IRF = Y;
end