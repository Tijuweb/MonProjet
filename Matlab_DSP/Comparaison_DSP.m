% Julien WEBER
%
% Projet - MEDISOV
%
% Dernière modification :
% 06/10/16
% 17h

clear variables;
close all;
clc;

[Entree, Fs]=audioread('10h.wav');  % Acquisition d'un signal audio avec Fs=16kHz


%% Déclaration variables
Hop_Size=160;       %160 pour Fs=16000    & 147 pour Fs=44100;
Wnd_Len=(2 * Hop_Size);
Nfft=1024;          % Taille de la FFT

Complex_data=zeros(1,Wnd_Len);
OverlapIn=zeros(1,Wnd_Len);
BufferZer=zeros(1,Wnd_Len);

Mean_FFT=zeros(1,Nfft/2);
Half_FFT=zeros(1,Nfft/2);
BufferFFT=zeros(1,Nfft);


%% Détermination des fréquences pour les bandes d'octavess
% freq_centre = [   31.5    40    50    63    80    100   125  160   200   250   315 400  500  630  800  1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 ];
% freq_min = freq_centre/(2^(1/6));
% freq_max = freq_centre*(2^(1/6));
% freq_min(1)=0;
% 
% axe_fft = 0:Fs/Nfft:Fs/2;
% axefOctave=zeros(2,length(freq_centre));   
% 
% for k=1:length(freq_centre)
% ind = find(axe_fft<=freq_max(k) & axe_fft>=freq_min(k));
%     if(~isempty(ind))   
%         axefOctave(1,k)=max(ind);
%         axefOctave(2,k)=length(ind);
%     end
% end

%%
Window = hanning(Wnd_Len,'periodic');       % Création de la fenetre de Hanning

load A_Lin.mat;             % Chargement de la pondération en A

 N=floor(length(Entree)/Wnd_Len);    % Nombre de boucle de calculs pour chaque frame
 Nframes = 2*N-1;                    % Nombre de frames totale
 Nframe_1s=Fs/Hop_Size;            % Nombre de frame pour 1 seconde


if(round(Nframe_1s)-Nframe_1s~=0)
    disp('WARNING, Nframe_1s n''est pas entier (Fs/Hop_size)');
end

meanEnergy=zeros(1,Nframe_1s);
Leq=zeros(1,floor(Nframes/(Nframe_1s)));
LeqMinute=zeros(1,60);
LA10=zeros(1,2);
LA50=zeros(1,2);
LA90=zeros(1,2);

C=96.5;     % Constante de calibration

k=1;
sec=1;
min=1;

%% Boucle de process
for iFrame=1:Nframes    % Boucle sur le nombre de frame
       
  OverlapIn=Entree((iFrame-1)*Hop_Size+1:(iFrame+1)*Hop_Size);  % Acquisition du signal

  BufferZer = OverlapIn.*Window;        % Application de la fenetre de Hanning
  BufferFFT= [ BufferZer' , zeros(1,Nfft-length(BufferZer)) ];  % Zero padding
  % Magnitude = BufferFFT*2^15;
    
    Complex_data=fft(BufferFFT,Nfft)/Nfft; % <=> rfft(BufferFFT,Nfft*2,SCALE)
   % Real=real(Complex_data(1:Nfft/2));     % Partie Réelle de la FFT
   % Imag=imag(Complex_data(1:Nfft/2));     % Partie Imaginaire de la FFT
    Half_FFT=abs(Complex_data(1:Nfft/2));   % Valeur absolue de la FFT
    HalfA_FFT=Half_FFT.*A;                  % Application de la pondération en A
    TotalEnergy= sum(HalfA_FFT.^2)/length(Half_FFT);    % Calcul de l'energie totale de la frame
    meanEnergy(k) = TotalEnergy/((1/Fs)*length(OverlapIn)); %Calcul de l'energie moyenne de la frame

     if (k==Nframe_1s)      % Si on a traité 1 seconde de frame
     k=1;        
     Leq(sec)=10*log10(sum(meanEnergy))+C;      % Calcul du log energy 10*log10
     sec=sec+1;                                 % Nouvelle seconde
     meanEnergy(:)=0;                           % Remise à zéro du tableau meanEnergy
     
     if(mod(sec,60)==0 && min~=60)                  % Si on a traité 1 minute de frame
       LeqMinute=sort(Leq(60*(min-1)+1:60*min));    % Rangement des Leq sur 1 minute
       LA10(min)=LeqMinute(6);                      % Récupération de LA10, LA50, LA90,..  
       LA50(min)=LeqMinute(30);
       LA90(min)=LeqMinute(54);
       min=min+1;         
     end
     
   
    else 
        k=k+1;
    end           
end

figure(1);
plot(Leq);title('Valeur du LAeq');ylabel('LAeq (dBA)');xlabel('temps (s)');


dBA=analyzeSignal(Entree(:,1),Fs);      % Fonction Matlab trouvée sur FileExchange

figure(40);
plot(Leq(1:length(dBA))'-dBA);
title('Difference des deux methodes de calcul');ylabel('Amplitude (dBA)');xlabel('Temps (s)');
