function [Data, SampFreq, metaData] = ExtractMERdat(fileName, directory)
% ************************************************************************* 
% Extract file from .dat files
% Modified from function provided by INOMED - 'MER_bigFile'   
% ************************************************************************* 
% function [Data, SampFreq, metaData] = MER_dat(fileName, directory)
% Written on 10-04-2014     

% direc_init=cd;
fid = fopen(fileName);

% '-----------------------------------------------------------------------%
% ' Header-Structur of MER Software                                     %
%            
% header of version 3.1
% struct PlotWinInfo
% {
%   char PatName [40];  // Patientenkürzel
%   char Position [40];  // Position (für EMG nicht relevant)
%   int SamplingRate; 
%   int SpikeCount;      // für EMG nicht relevant
%   int TriggerLinie;    // für EMG nicht relevant
%   int SiteNr;          // für EMG nicht relevant
%   int KanalNr;          // Nummer des Messkanals
%   int op_id;           // ID der Operation
%   int SiteID;          // ID der Site
%   int MaxYValue;        // Maximale Auflösung (Amplitude) der Y - Achse (ISIS1: 4095, ISIS2: 0xFFFF)
%   int typ;             // Typ der Datei (0 = Normale Messdatendatei, 255 = EMG - Datei
%   int Multiplexing;    // Zur Zeit nicht benutzt
%   int EMGChannels;     // Zugeordnete EMG - Kanäle
%   char EMGChannelsDesc [3][40];  // Bezeichner der EMG - Kanäle
%   int BitsPerValue;            // Compression, wenn Wert <> 12  keine Kompression
%   int fcHPInHz;                     // Grenzfrequenz des Hochpass Hardware Filters
%   int fcLPInHz;                      // Grenzfrequenz des Tiefpass Hardware Filters
%   int vu;                                                  // Verstaerkung Hardware
%   int fcHPInHzSW;                              // Grenzfrequenz des Hochpass Software Filters
%   int fcLPInHzSW;                               // Grenzfrequenz des Tiefpass Software Filters
%   double bitResolution;// Bit-Aufloesung in nV / Bit
%   char versionNrMER[ 32 ];        // Versionsnummer der MER Software
%   char filler [294];               // Reserve für zukünftige Nutzung
%};

Header=fread(fid,80,'*char')';
SampFreq=double(fread(fid,1,'*int')');             % 80..81

Header=fread(fid,2,'*int')';                      
SiteNr=double(fread(fid,1,'*int')');            
KanalNr=double(fread(fid,1,'*int')');            
op_id=double(fread(fid,1,'*int')');  
SiteID=double(fread(fid,1,'*int')');
MaxYValue=double(fread(fid,1,'*int')');
typ=double(fread(fid,1,'*int')');
Multiplexing=double(fread(fid,1,'*int')');
EMGChannels=double(fread(fid,1,'*int')');

Header=fread(fid,120,'*char')';
 
BitsPerValue=double(fread(fid,1,'*int')');
fcHPInHz=double(fread(fid,1,'*int')');
fcLPInHz=double(fread(fid,1,'*int')');
vu=double(fread(fid,1,'*int')');
fcHPInHzSW=double(fread(fid,1,'*int')');
fcLPInHzSW=double(fread(fid,1,'*int')');
bitResolution=double(fread(fid,1,'*double')');
versionNrMER=char(fread(fid,32,'*char')');

fclose(fid); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read the data values and plot them:
fid = fopen(fileName);

Data=(double(fread(fid,20000000,'*ushort'))-(MaxYValue/2))*(0.4/MaxYValue);
% t=0:(1/SampFreq):((length(Data)-1)/SampFreq);

Header=fread(fid,604,'*char')';
metaData.SiteNr = SiteNr;
metaData.KanalNr = KanalNr;
metaData.op_id = op_id;
metaData.SampFreq = SampFreq;

% cd(direc_init);
fclose(fid);         