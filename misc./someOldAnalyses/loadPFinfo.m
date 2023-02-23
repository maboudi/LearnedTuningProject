function  [activeUnits, tuningRL, tuningLR, RLtemplate, LRtemplate, positionBins]= loadPFinfo(fileinfo, longORshort)

currentdir = pwd;
FileBase = [currentdir '/' fileinfo.name '/data part' num2str(longORshort) '/PlaceFields'];

activeUnits = load([FileBase '/leftward/activeUnits.mat']); %% all the active units with qClu in a range 
activeUnits = activeUnits.allActUnits;

tuningRL = load([FileBase '/leftward/AllCellsTuning_leftward.mat']);
tuningRL = tuningRL.MeanoverLaps;

tuningLR = load([FileBase '/rightward/AllCellsTuning_rightward.mat']);
tuningLR = tuningLR.MeanoverLaps;

RLtemplate = load([FileBase '/leftward/sortedUnits_leftward.mat']);
RLtemplate = RLtemplate.trajTemplate;

LRtemplate = load([FileBase '/rightward/sortedUnits_rightward.mat']);
LRtemplate = LRtemplate.trajTemplate;

load([FileBase '/leftward/PositionBinCenters.mat']);

end
