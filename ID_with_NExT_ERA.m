function [P50,ID,GlobalChannels] = ID_with_NExT_ERA(varargin)
%
% [P50,ID] = ID_with_NExT_ERA(varin)
%
% Función para realizar el procesamiento de señales.
%
% INPUTS:
%   data: matriz de datos organizados en columnas con aceleraciones y
%       (opcionalmente) la primera columna con un vector de tiempo.
%   IsTimeVectorIncluded: valor lógico que indica si la primera columna de
%       data es el vector de tiempo.
%   fs: frecuencia de muestreo de data.
%   fr: frecuencia de *remuestreo* o [].
%   ffi: frecuencia de corte del filtro pasa alto o [].
%   fff: frecuencia de corte del filtro pasa bajo o [].
%   Wndw: puede ser [] o vector con dos valores en *número de puntos* que
%       limitan la ventana de datos seleccionada. Ejemplo: Wndw = [1, 578];
%       (seleccionar los datos entre el segundo 1/fs y el 578/fs).
%   Trend: valor lógico (1:true, 0:false) para realizar detrend.
%   pAcel: valor lógico (1:true, 0:false) indicando si graficar la señal en
%       tiempo o no.
%   RefChannel: escalar o vector indicando los canales de referencia en
%       data. Es decir, si la primera columna de aceleración es de
%       referencia, entonces RefChannel=1.
%   ModelOrderRange: escalar o vector indicando el orden del modelo que se
%       empleará en NExT-ERA para la identificación.
%   NumClusters: número de grupos (modos de vibración) que el Cluster
%       intentará identificar.
%   MAC_Treshold: valor mínimo aceptable del MAC.
%   PLOT_MAC: valor lógico (1:true, 0:false) para graficar modos de vibración.
%   PLOT_FFT: valor lógico (1:true, 0:false) para graficar FFT.
%   FreqAprox: en caso que no se conozcan indicar [], en caso que se deseen
%       incluir indicar rango de cada frecuencia como vectores fila. Es
%       decir, FreqAprox(2,1)=frecuencia mínima asociada al segundo modo de
%       vibración y FreqAprox(2,2)= frecuencia máxima asociada al segundo
%       modo de vibración.
%   GlobalChannels: vector indicando las columnas de aceleración en data a
%       qué orden global corresponden. La finalidad de este vector es poder
%       comparar los resultados obtenidos a partir de diferentes registros.
%       Por ejemplo, en caso que se tengan los sensores _023 (cada sensor
%       con las direcciones), entonces GlobalChannels = [1:3,7:12];
%       En caso que no se cuente con esta información, indicar [].
%
% OUTPUTS:
%   SignOut: estructura con las propiedades de la señal procesada.
%
% %%%%%%%%%%%%%%%%%%
% %%%%% M.G.H. %%%%%
% %%% 2020/06/01 %%%
% %%%%%%%%%%%%%%%%%%
tic_ID_with_NExT_ERA = tic;
%%% Assign input variables ------------------------------------------------
data = varargin{1};
IsTimeVectorIncluded = varargin{2};
fs = varargin{3};
fr = varargin{4};
ffi = varargin{5};
fff = varargin{6};
Wndw = varargin{7};
Trend = varargin{8};
pAcel = varargin{9};
RefChannel = varargin{10};
ModelOrderRange = varargin{11};
NumClusters = varargin{12};
MAC_Treshold = varargin{13};
PLOT_MAC = varargin{14};
PLOT_FFT = varargin{15};
FreqAprox = varargin{16};
GlobalChannels = varargin{17};
%%% Prepare input signals -------------------------------------------------
if IsTimeVectorIncluded == 0
    t = (0:size(data,1)-1)./fs; % (s)
    data = [t(:),data];
end
%%%
if isempty(fr)
    fr = fs;
elseif fr == 0
    fr = fs;
end
%%%
if ~isempty(fff)
    if fff >= 0.4*fr
        warning('Cutoff low-pass filter will be redefined as fff = 0.4*fr.')
        fff = 0.4*fr; % (Hz)
    end
end
%%%
if ~isempty(pAcel)
    if pAcel==1
        pAcel = 'Acel. (m/{s^2})';
    else
        pAcel = ' ';
    end
end
%%%
SignOut = SigPro(...
    data(:,1),... % ti
    data(:,2:end),... % yi
    fr,... % fs
    Wndw,... % Wndw
    Trend,... % Trend
    1,... % ts
    ffi,... % ffi
    fff,... % fff
    pAcel,... % plotts
    PLOT_FFT,... % plotfft
    0,... % plotpsd
    0); % plotspt
%%% Check zero signals ----------------------------------------------------
if ~isempty(GlobalChannels)
    ZeroSignals = sum(data(:,2:end),1)==0;
    if any(ZeroSignals)
        DissmisThisChannel = find(ZeroSignals);
        warning(['Los siguientes canales serán excluidos del análisis por no tener información: ',num2str(DissmisThisChannel)])
        SignOut.y(:,DissmisThisChannel) = [];
        GlobalChannels(DissmisThisChannel) = [];
        [ZeroRef,DissmisThisRef] = ismember(DissmisThisChannel,RefChannel);
        %%% RefChannel ----------------------------------------------------
        if any(ZeroRef)
            for ii = 1:length(ZeroRef)
                if ZeroRef(ii)
                    RefChannel(DissmisThisRef(ii)) = [];
                end
            end
        end
        %%% ---------------------------------------------------------------
    end
end
%%% -----------------------------------------------------------------------
%%% ID with NExT_ERA
PropNERA.L     = size(SignOut.y,1);
PropNERA.MxLgs = floor(PropNERA.L/2);
PropNERA.N     = 10; % Número de ventanas
PropNERA.nfft  = pow2(nextpow2(PropNERA.L/PropNERA.N))*PropNERA.N;
PropNERA.Nvrlp = 0.50; % Razón de traslapo (0.5=50%)
PropNERA.Windw = floor(PropNERA.L/((PropNERA.N+1)*PropNERA.Nvrlp)/2)*2;
%%%
[Wn,Zn,Ph,ModelOrder] = NExT_ERA(...
    SignOut.y',...
    SignOut.fs,...
    RefChannel,...
    PropNERA.MxLgs,...
    PropNERA.Windw,...
    PropNERA.Nvrlp,...
    PropNERA.N,...
    ModelOrderRange,...
    FreqAprox);
%%% Sort
[Wn,Ind] = sort(Wn);
Ph = Ph(:,Ind);
Zn = Zn(Ind);
Fn = Wn./(2*pi);
ModelOrder = ModelOrder(Ind);
%%% -----------------------------------------------------------------------
%%% Clustering
ClusterIndex = getCluster(Fn,Ph,NumClusters,MAC_Treshold,PLOT_MAC,FreqAprox);
plotCluster(ModelOrder,Fn,ClusterIndex,'Freq. (Hz)','Model Order',[SignOut.f,abs(SignOut.Y)])
%%% -----------------------------------------------------------------------
%%% Organize data
jj = 0;
for ii = 1:size(ClusterIndex,2)
    if sum(ClusterIndex(:,ii))>0
        jj = jj+1;
        %%%
        Location = ClusterIndex(:,ii);
        %%%
        Ph_ii = Ph(:,Location);
        Ph_ii_max = ones(size(Ph_ii,1),1)*(max(abs(Ph_ii)).*sign(Ph_ii(1,:)));
        Ph_ii = Ph_ii./Ph_ii_max;
        %%%
        ID.Fn{jj,1}  = Fn(Location);
        ID.Wn{jj,1}  = Wn(Location);
        ID.Zn{jj,1}  = Zn(Location);
        ID.Ph{1,jj}  = Ph_ii;
        if sum(Location)==1
            P50_Fn(jj,1) = Fn(Location); %#ok<AGROW>
            P50_Zn(jj,1) = Zn(Location); %#ok<AGROW>
            P50_Ph(:,jj) = Ph_ii; %#ok<AGROW>
        else
            P50_Fn(jj,1) = prctile(Fn(Location),50); %#ok<AGROW>
            P50_Zn(jj,1) = prctile(Zn(Location),50); %#ok<AGROW>
            P50_Ph(:,jj) = prctile(Ph_ii',50)'; %#ok<AGROW>
        end
    end
end
%%%
if PLOT_MAC==1
    Ph = cell2mat(ID.Ph);
    MAC = getMAC(Ph,Ph);
    plotMAC(MAC,'Después del Clúster');
end
%%%
[~,Ind] = sort(P50_Fn);
P50_Fn = P50_Fn(Ind);
P50_Zn = P50_Zn(Ind);
P50_Ph = P50_Ph(:,Ind);
ID.Fn = ID.Fn(Ind,1);
ID.Wn = ID.Wn(Ind,1);
ID.Zn = ID.Zn(Ind,1);
ID.Ph = ID.Ph(:,Ind);
%%%
disp(' '); disp(table(P50_Fn,P50_Zn)); disp(' ');
%%%
P50.Fn = P50_Fn;
P50.Zn = P50_Zn;
P50.Ph = P50_Ph;
%%% -----------------------------------------------------------------------
disp(['ID_with_NExT_ERA: ',num2str(toc(tic_ID_with_NExT_ERA),'%.3f')])
end
